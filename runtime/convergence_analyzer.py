"""Sprint 2 — Convergence and stability analysis for GROMACS XVG data.

Provides:
    analyze_rmsd_convergence(data, threshold_nm) → ConvergenceResult
    analyze_energy_stability(data) → StabilityResult
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from runtime.xvg_parser import XVGData


# ─────────────────────────────────────────────────────────────────────────────
# Result dataclasses
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class ConvergenceResult:
    converged:       bool
    plateau_ns:      Optional[float]   # time (ns) when plateau was reached, or None
    mean_last20pct:  float             # mean of last 20% of data
    std_last20pct:   float             # std of last 20% of data
    drift:           float             # slope of linear fit of last 20% (≈ 0 = converged)
    verdict:         str               # human-readable summary


@dataclass
class StabilityResult:
    stable:        bool
    mean:          float
    std:           float
    drift_per_ns:  float    # slope of linear fit (units per ns)
    verdict:       str


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

def _mean(vals: list[float]) -> float:
    if not vals:
        return 0.0
    return sum(vals) / len(vals)


def _std(vals: list[float]) -> float:
    if len(vals) < 2:
        return 0.0
    m = _mean(vals)
    variance = sum((v - m) ** 2 for v in vals) / (len(vals) - 1)
    return variance ** 0.5


def _linear_slope(xs: list[float], ys: list[float]) -> float:
    """Least-squares slope (dy/dx). Returns 0.0 if degenerate."""
    n = len(xs)
    if n < 2:
        return 0.0
    sx  = sum(xs)
    sy  = sum(ys)
    sxy = sum(x * y for x, y in zip(xs, ys))
    sxx = sum(x * x for x in xs)
    denom = n * sxx - sx * sx
    if abs(denom) < 1e-12:
        return 0.0
    return (n * sxy - sx * sy) / denom


# ─────────────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────────────

def analyze_rmsd_convergence(
    data:           XVGData,
    threshold_nm:   float = 0.3,
) -> ConvergenceResult:
    """Analyze RMSD convergence from XVG data.

    RMSD is considered converged if std of the last 20% of values < threshold_nm.
    The plateau_ns is estimated as the time at which the RMSD first enters the
    final-20% mean ± threshold band and stays there.
    """
    if not data.series or not data.time_ps:
        return ConvergenceResult(
            converged      = False,
            plateau_ns     = None,
            mean_last20pct = 0.0,
            std_last20pct  = 0.0,
            drift          = 0.0,
            verdict        = "No data available for RMSD convergence analysis.",
        )

    values = data.series[0].values   # first data column = RMSD
    times  = data.time_ps

    # Align lengths (should be equal but be defensive)
    n = min(len(values), len(times))
    if n == 0:
        return ConvergenceResult(
            converged=False, plateau_ns=None,
            mean_last20pct=0.0, std_last20pct=0.0, drift=0.0,
            verdict="Empty data.",
        )

    values = values[:n]
    times  = times[:n]

    # Last 20% window
    window_start = max(1, int(n * 0.8))
    last_vals    = values[window_start:]
    last_times   = times[window_start:]

    mean_last = _mean(last_vals)
    std_last  = _std(last_vals)

    # Convert times to ns for the slope
    last_times_ns = [t / 1000.0 for t in last_times]
    drift = _linear_slope(last_times_ns, last_vals)   # nm/ns

    converged = std_last < threshold_nm

    # Estimate plateau_ns: first time RMSD enters [mean±threshold] and stays
    plateau_ns: Optional[float] = None
    if converged:
        lo = mean_last - threshold_nm
        hi = mean_last + threshold_nm
        for i in range(n - 1, -1, -1):
            if not (lo <= values[i] <= hi):
                # Last point outside band → plateau started just after i
                if i + 1 < n:
                    plateau_ns = times[i + 1] / 1000.0
                break
        else:
            # All points in band
            plateau_ns = times[0] / 1000.0

    # Build verdict
    if converged:
        parts = [
            f"RMSD converged: std={std_last:.4f} nm < {threshold_nm} nm threshold.",
            f"Mean (last 20%): {mean_last:.4f} nm.",
        ]
        if plateau_ns is not None:
            parts.append(f"Plateau reached at ~{plateau_ns:.2f} ns.")
        if abs(drift) > 0.01:
            parts.append(f"Residual drift: {drift:+.4f} nm/ns.")
    else:
        parts = [
            f"RMSD NOT converged: std={std_last:.4f} nm >= {threshold_nm} nm.",
            f"Mean (last 20%): {mean_last:.4f} nm.",
            f"Drift: {drift:+.4f} nm/ns.",
            "Consider longer simulation or verify equilibration.",
        ]

    return ConvergenceResult(
        converged      = converged,
        plateau_ns     = plateau_ns,
        mean_last20pct = mean_last,
        std_last20pct  = std_last,
        drift          = drift,
        verdict        = " ".join(parts),
    )


def analyze_energy_stability(data: XVGData) -> StabilityResult:
    """Analyze energy stability from XVG data.

    Energy is stable if drift < 1% of |mean| per ns.
    Uses the first data series (typically potential energy).
    """
    if not data.series or not data.time_ps:
        return StabilityResult(
            stable=False, mean=0.0, std=0.0, drift_per_ns=0.0,
            verdict="No data available for energy stability analysis.",
        )

    values = data.series[0].values
    times  = data.time_ps

    n = min(len(values), len(times))
    if n == 0:
        return StabilityResult(
            stable=False, mean=0.0, std=0.0, drift_per_ns=0.0,
            verdict="Empty data.",
        )

    values = values[:n]
    times  = times[:n]

    mean_val = _mean(values)
    std_val  = _std(values)

    times_ns   = [t / 1000.0 for t in times]
    drift_per_ns = _linear_slope(times_ns, values)

    # Stability criterion: |drift| < 1% of |mean| per ns
    threshold = abs(mean_val) * 0.01 if abs(mean_val) > 1e-6 else 1.0
    stable = abs(drift_per_ns) < threshold

    if stable:
        verdict = (
            f"Energy stable: drift={drift_per_ns:+.2f} kJ/mol/ns "
            f"({abs(drift_per_ns)/max(abs(mean_val),1e-6)*100:.2f}% of mean). "
            f"Mean={mean_val:.1f} kJ/mol, std={std_val:.1f} kJ/mol."
        )
    else:
        verdict = (
            f"Energy NOT stable: drift={drift_per_ns:+.2f} kJ/mol/ns "
            f"({abs(drift_per_ns)/max(abs(mean_val),1e-6)*100:.2f}% of mean). "
            f"Mean={mean_val:.1f} kJ/mol, std={std_val:.1f} kJ/mol. "
            "Verify equilibration and simulation parameters."
        )

    return StabilityResult(
        stable       = stable,
        mean         = mean_val,
        std          = std_val,
        drift_per_ns = drift_per_ns,
        verdict      = verdict,
    )
