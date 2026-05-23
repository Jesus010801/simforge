"""Phase A — Scientific Interpretation Layer: quality classifier.

Classifies a molecular dynamics simulation run from XVG data into one of
five quality tiers, with confidence, evidence, warnings, and recommendations.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from runtime.xvg_parser import XVGData


class RunQuality(str, Enum):
    CONVERGED           = "converged"
    PARTIALLY_CONVERGED = "partially_converged"
    NOT_CONVERGED       = "not_converged"
    PROBLEMATIC         = "problematic"
    INSUFFICIENT_DATA   = "insufficient_data"


# ─────────────────────────────────────────────────────────────────────────────
# Internal statistics helpers (pure Python, no numpy)
# ─────────────────────────────────────────────────────────────────────────────

def _mean(vals: list[float]) -> float:
    if not vals:
        return 0.0
    return sum(vals) / len(vals)


def _std(vals: list[float]) -> float:
    if len(vals) < 2:
        return 0.0
    m = _mean(vals)
    return (sum((v - m) ** 2 for v in vals) / (len(vals) - 1)) ** 0.5


def _linear_slope(xs: list[float], ys: list[float]) -> float:
    """Least-squares slope dy/dx. Returns 0.0 if degenerate."""
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


def _clamp(val: float, lo: float = 0.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, val))


# ─────────────────────────────────────────────────────────────────────────────
# QualityReport
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class QualityReport:
    quality:         RunQuality
    confidence:      float           # 0.0–1.0
    evidence:        list[str]       # human-readable evidence
    warnings:        list[str]       # anomalies detected
    recommendations: list[str]       # what to do next
    metrics:         dict            # raw numeric values

    def as_dict(self) -> dict:
        return {
            "quality":         self.quality.value,
            "confidence":      round(self.confidence, 4),
            "evidence":        self.evidence,
            "warnings":        self.warnings,
            "recommendations": self.recommendations,
            "metrics":         self.metrics,
        }

    def as_markdown(self) -> str:
        badge = {
            RunQuality.CONVERGED:           "CONVERGED ✓",
            RunQuality.PARTIALLY_CONVERGED: "PARTIALLY CONVERGED ~",
            RunQuality.NOT_CONVERGED:       "NOT CONVERGED ✗",
            RunQuality.PROBLEMATIC:         "PROBLEMATIC ⚠",
            RunQuality.INSUFFICIENT_DATA:   "INSUFFICIENT DATA ?",
        }.get(self.quality, self.quality.value)

        lines = [
            f"# Quality Classification: {badge}",
            f"",
            f"**Confidence:** {self.confidence:.0%}",
            "",
        ]

        if self.evidence:
            lines += ["## Evidence", ""]
            for e in self.evidence:
                lines.append(f"- {e}")
            lines.append("")

        if self.warnings:
            lines += ["## Warnings", ""]
            for w in self.warnings:
                lines.append(f"- {w}")
            lines.append("")

        if self.recommendations:
            lines += ["## Recommendations", ""]
            for r in self.recommendations:
                lines.append(f"- {r}")
            lines.append("")

        if self.metrics:
            lines += ["## Key Metrics", ""]
            for k, v in self.metrics.items():
                lines.append(f"- **{k}:** {v}")
            lines.append("")

        return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# Internal analysis helpers
# ─────────────────────────────────────────────────────────────────────────────

def _min_data_ns() -> float:
    return 1.0   # minimum simulation time (ns) for meaningful analysis


def _duration_ns(data: "XVGData") -> float:
    """Return total simulation time in ns from an XVGData object."""
    if not data.time_ps:
        return 0.0
    return data.time_ps[-1] / 1000.0


def _has_enough_data(data: Optional["XVGData"]) -> bool:
    """True when data has >= 1 ns of simulation time and at least one series."""
    if data is None or not data.series or not data.time_ps:
        return False
    return _duration_ns(data) >= _min_data_ns()


def _last20(vals: list[float]) -> list[float]:
    n = len(vals)
    start = max(1, int(n * 0.8))
    return vals[start:]


def _analyze_rmsd(data: "XVGData") -> dict:
    """Return dict of RMSD statistics for the last 20% window."""
    vals  = data.series[0].values
    times = data.time_ps
    n     = min(len(vals), len(times))
    vals  = vals[:n]
    times = times[:n]

    last_vals  = _last20(vals)
    last_times = _last20(times)

    mean_last = _mean(last_vals)
    std_last  = _std(last_vals)
    last_ns   = [t / 1000.0 for t in last_times]
    drift     = _linear_slope(last_ns, last_vals)   # nm/ns
    max_rmsd  = max(vals) if vals else 0.0
    total_ns  = times[-1] / 1000.0 if times else 0.0

    return {
        "mean_last20": mean_last,
        "std_last20":  std_last,
        "drift":       drift,
        "max_rmsd":    max_rmsd,
        "total_ns":    total_ns,
    }


def _analyze_energy(data: "XVGData") -> dict:
    """Return dict of energy statistics for the full trajectory."""
    vals  = data.series[0].values
    times = data.time_ps
    n     = min(len(vals), len(times))
    vals  = vals[:n]
    times = times[:n]

    mean_val   = _mean(vals)
    std_val    = _std(vals)
    times_ns   = [t / 1000.0 for t in times]
    drift      = _linear_slope(times_ns, vals)

    pct_std  = abs(std_val  / mean_val * 100) if abs(mean_val) > 1e-6 else 0.0
    pct_drift = abs(drift   / mean_val * 100) if abs(mean_val) > 1e-6 else 0.0
    total_ns  = times[-1] / 1000.0 if times else 0.0

    return {
        "mean":       mean_val,
        "std":        std_val,
        "drift":      drift,
        "pct_std":    pct_std,
        "pct_drift":  pct_drift,
        "total_ns":   total_ns,
    }


# ─────────────────────────────────────────────────────────────────────────────
# classify_run — public API
# ─────────────────────────────────────────────────────────────────────────────

def classify_run(
    rmsd_data:        Optional["XVGData"] = None,
    energy_data:      Optional["XVGData"] = None,
    pressure_data:    Optional["XVGData"] = None,
    temperature_data: Optional["XVGData"] = None,
    context:          Optional[str] = None,
) -> QualityReport:
    """Classify simulation quality from XVG data.

    Works with partial data (only RMSD, only energy, or both).
    Confidence is higher when both sources agree.

    When ``context`` is provided (e.g. "globular_protein", "membrane_protein"),
    uses the md_knowledge base for context-aware soft ranges and evidence.
    """
    if context is not None:
        return _classify_with_context(
            context=context,
            rmsd_data=rmsd_data,
            energy_data=energy_data,
            pressure_data=pressure_data,
            temperature_data=temperature_data,
        )
    evidence:        list[str] = []
    warnings:        list[str] = []
    recommendations: list[str] = []
    metrics:         dict      = {}

    rmsd_ok   = _has_enough_data(rmsd_data)
    energy_ok = _has_enough_data(energy_data)

    # ── Check for insufficient data ──────────────────────────────────────────
    if not rmsd_ok and not energy_ok:
        return QualityReport(
            quality         = RunQuality.INSUFFICIENT_DATA,
            confidence      = 1.0,
            evidence        = ["No RMSD or energy data with at least 1 ns of simulation time."],
            warnings        = [],
            recommendations = [
                "Run the simulation analysis steps to generate XVG files.",
                "Ensure the simulation has completed at least 1 ns.",
            ],
            metrics         = {},
        )

    # ── Collect statistics ───────────────────────────────────────────────────
    rmsd_stats:   Optional[dict] = None
    energy_stats: Optional[dict] = None

    if rmsd_ok:
        rmsd_stats = _analyze_rmsd(rmsd_data)   # type: ignore[arg-type]
        metrics["rmsd_mean_last20_nm"]  = round(rmsd_stats["mean_last20"], 4)
        metrics["rmsd_std_last20_nm"]   = round(rmsd_stats["std_last20"],  4)
        metrics["rmsd_drift_nm_per_ns"] = round(rmsd_stats["drift"],       4)
        metrics["rmsd_max_nm"]          = round(rmsd_stats["max_rmsd"],    4)
        metrics["rmsd_total_ns"]        = round(rmsd_stats["total_ns"],    2)

    if energy_ok:
        energy_stats = _analyze_energy(energy_data)   # type: ignore[arg-type]
        metrics["energy_mean_kJ_mol"]        = round(energy_stats["mean"],      1)
        metrics["energy_std_kJ_mol"]         = round(energy_stats["std"],       1)
        metrics["energy_drift_kJ_mol_per_ns"]= round(energy_stats["drift"],     2)
        metrics["energy_pct_std"]            = round(energy_stats["pct_std"],   2)
        metrics["energy_pct_drift_per_ns"]   = round(energy_stats["pct_drift"], 2)
        metrics["energy_total_ns"]           = round(energy_stats["total_ns"],  2)

    # ── Check for PROBLEMATIC conditions (highest priority) ──────────────────
    # These override everything else.
    problematic = False

    if rmsd_stats is not None:
        if rmsd_stats["max_rmsd"] > 1.0:
            problematic = True
            warnings.append(
                f"RMSD exceeded 1.0 nm (max={rmsd_stats['max_rmsd']:.3f} nm) — "
                "system may have exploded or undergone a major conformational change."
            )
            recommendations.append("Inspect trajectory visually for aggregation or explosion.")
            recommendations.append("Check periodic boundary conditions and forcefield parameters.")

    if energy_stats is not None:
        if energy_stats["pct_drift"] > 5.0:
            problematic = True
            warnings.append(
                f"Energy drift exceeds 5% per ns ({energy_stats['pct_drift']:.1f}%/ns) — "
                "system is likely not equilibrated."
            )
            recommendations.append("Extend equilibration phase significantly.")
            recommendations.append("Check thermostat and barostat coupling constants.")

        if energy_stats["pct_std"] > 10.0:
            problematic = True
            warnings.append(
                f"Energy fluctuations are very large (std={energy_stats['pct_std']:.1f}% of mean) — "
                "possible numerical instability."
            )
            recommendations.append("Check timestep (reduce if > 2 fs) and constraints.")
            recommendations.append("Verify forcefield parameter compatibility.")

    if problematic:
        # Add any positive evidence we do have
        if rmsd_stats is not None:
            evidence.append(
                f"RMSD analysis available ({rmsd_stats['total_ns']:.1f} ns); "
                f"mean={rmsd_stats['mean_last20']:.3f} nm in last 20%."
            )
        if energy_stats is not None:
            evidence.append(
                f"Energy analysis available ({energy_stats['total_ns']:.1f} ns); "
                f"mean={energy_stats['mean']:.1f} kJ/mol."
            )

        # Confidence based on how clearly problematic
        n_problems = len(warnings)
        confidence = _clamp(0.8 + 0.1 * (n_problems - 1))

        if not recommendations:
            recommendations.append("Review system setup and equilibration protocol.")

        return QualityReport(
            quality         = RunQuality.PROBLEMATIC,
            confidence      = confidence,
            evidence        = evidence,
            warnings        = warnings,
            recommendations = recommendations,
            metrics         = metrics,
        )

    # ── Evaluate convergence signals ─────────────────────────────────────────
    rmsd_converged:   Optional[bool] = None
    rmsd_drifting:    bool           = False
    energy_stable:    Optional[bool] = None

    if rmsd_stats is not None:
        std = rmsd_stats["std_last20"]
        drift = rmsd_stats["drift"]

        if std < 0.15 and abs(drift) < 0.01:
            rmsd_converged = True
            evidence.append(
                f"RMSD plateau detected: std={std:.4f} nm (< 0.15 nm threshold), "
                f"drift={drift:+.4f} nm/ns (≈ 0)."
            )
        elif std < 0.30:
            rmsd_converged = False   # partially
            if drift > 0.01:
                rmsd_drifting = True
                evidence.append(
                    f"RMSD partially stabilized (std={std:.4f} nm < 0.30 nm) "
                    f"but still drifting ({drift:+.4f} nm/ns)."
                )
            else:
                evidence.append(
                    f"RMSD partially stabilized: std={std:.4f} nm, drift={drift:+.4f} nm/ns."
                )
        else:
            # std >= 0.30 — still growing
            rmsd_converged = False
            rmsd_drifting  = True
            evidence.append(
                f"RMSD not converged: std={std:.4f} nm (>= 0.30 nm threshold), "
                f"drift={drift:+.4f} nm/ns."
            )
            recommendations.append("Extend simulation to allow RMSD to plateau.")

    if energy_stats is not None:
        pct_std = energy_stats["pct_std"]
        if pct_std < 1.0:
            energy_stable = True
            evidence.append(
                f"Energy stable: fluctuations={pct_std:.2f}% of mean "
                f"(< 1% threshold), mean={energy_stats['mean']:.1f} kJ/mol."
            )
        else:
            energy_stable = False
            evidence.append(
                f"Energy not yet fully stable: fluctuations={pct_std:.2f}% of mean "
                f"(>= 1% threshold), mean={energy_stats['mean']:.1f} kJ/mol."
            )
            recommendations.append("Consider extending equilibration phase.")

    # ── NOT_CONVERGED: RMSD growing ──────────────────────────────────────────
    if rmsd_stats is not None and rmsd_drifting and rmsd_stats["std_last20"] >= 0.30:
        # std >= 0.30 nm is clear evidence of non-convergence.
        # Confidence grows with std magnitude; base = 0.7 for std=0.30.
        std_factor = _clamp((rmsd_stats["std_last20"] - 0.30) / 0.50)  # 0→1 for std 0.30→0.80
        confidence = _clamp(0.70 + std_factor * 0.25)
        if not recommendations:
            recommendations.append("Extend simulation length to reach RMSD plateau.")
        return QualityReport(
            quality         = RunQuality.NOT_CONVERGED,
            confidence      = confidence,
            evidence        = evidence,
            warnings        = warnings,
            recommendations = recommendations,
            metrics         = metrics,
        )

    # ── CONVERGED: RMSD converged AND energy stable (or only one available) ──
    if rmsd_converged is True and (energy_stable is True or energy_stable is None):
        # Confidence scales with data length and agreement
        total_ns = rmsd_stats["total_ns"] if rmsd_stats else 0.0
        if energy_stats:
            total_ns = max(total_ns, energy_stats["total_ns"])
        length_bonus = _clamp((total_ns - 1.0) / 99.0) * 0.15   # up to +0.15 for 100 ns
        both_agree   = 0.05 if (rmsd_stats and energy_stats and energy_stable) else 0.0
        confidence   = _clamp(0.85 + length_bonus + both_agree)

        if not recommendations:
            recommendations.append("Simulation appears well-converged. Proceed with analysis.")

        return QualityReport(
            quality         = RunQuality.CONVERGED,
            confidence      = confidence,
            evidence        = evidence,
            warnings        = warnings,
            recommendations = recommendations,
            metrics         = metrics,
        )

    # ── PARTIALLY_CONVERGED ───────────────────────────────────────────────────
    # Covers:
    #   - rmsd partially OK but energy not stable
    #   - rmsd converged but energy not stable
    #   - only energy data that is not fully stable
    #   - only rmsd data that is partially converged

    total_ns = 0.0
    if rmsd_stats:
        total_ns = max(total_ns, rmsd_stats["total_ns"])
    if energy_stats:
        total_ns = max(total_ns, energy_stats["total_ns"])

    length_factor = _clamp(total_ns / 50.0)   # approaches 1 at 50 ns
    confidence    = _clamp(0.5 + length_factor * 0.3)

    if not recommendations:
        recommendations.append("Run longer to confirm convergence.")
        recommendations.append("Check RMSD trend over additional time windows.")

    return QualityReport(
        quality         = RunQuality.PARTIALLY_CONVERGED,
        confidence      = confidence,
        evidence        = evidence,
        warnings        = warnings,
        recommendations = recommendations,
        metrics         = metrics,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Context-aware classify (uses core/md_knowledge)
# ─────────────────────────────────────────────────────────────────────────────

def _classify_with_context(
    context: str,
    rmsd_data:        Optional["XVGData"] = None,
    energy_data:      Optional["XVGData"] = None,
    pressure_data:    Optional["XVGData"] = None,
    temperature_data: Optional["XVGData"] = None,
) -> QualityReport:
    """Use the md_knowledge scientific knowledge base for context-aware classification."""
    from core.md_knowledge.contexts import SystemContext
    from core.md_knowledge.interpreter import interpret_simulation, ObservableResult
    from core.md_knowledge.states import SimulationState

    # Resolve context string → enum, falling back to UNKNOWN
    try:
        sys_ctx = SystemContext(context.lower().replace("-", "_").replace(" ", "_"))
    except ValueError:
        sys_ctx = SystemContext.UNKNOWN

    observables: list[ObservableResult] = []

    def _xvg_to_obs(name: str, data: Optional["XVGData"]) -> None:
        if data is None or not data.series or not data.time_ps:
            return
        vals  = data.series[0].values
        times = [t / 1000.0 for t in data.time_ps]
        n = min(len(vals), len(times))
        if n < 5:
            return
        observables.append(ObservableResult(name=name, values=vals[:n], times_ns=times[:n]))

    _xvg_to_obs("rmsd",             rmsd_data)
    _xvg_to_obs("potential_energy", energy_data)
    _xvg_to_obs("pressure",         pressure_data)
    _xvg_to_obs("temperature",      temperature_data)

    result = interpret_simulation(observables, sys_ctx)

    # Map InterpretationResult → QualityReport
    quality_map = {
        "EXCELLENT":  RunQuality.CONVERGED,
        "GOOD":       RunQuality.CONVERGED,
        "ACCEPTABLE": RunQuality.PARTIALLY_CONVERGED,
        "POOR":       RunQuality.NOT_CONVERGED,
        "FAILED":     RunQuality.PROBLEMATIC,
    }
    if not observables:
        quality = RunQuality.INSUFFICIENT_DATA
    else:
        quality = quality_map.get(result.quality_tier, RunQuality.PARTIALLY_CONVERGED)

    # Collapse evidence items → strings
    evidence_strs = [ev.message for ev in result.evidence.items]

    # Build metrics from pattern results
    metrics: dict = {"context": sys_ctx.value, "state": result.state.value}
    for obs_name, pr in result.pattern_results.items():
        metrics[f"{obs_name}_pattern"] = pr.pattern.value
        metrics[f"{obs_name}_mean"] = round(pr.mean_last, 4)
        metrics[f"{obs_name}_std"]  = round(pr.plateau_std, 4)

    return QualityReport(
        quality         = quality,
        confidence      = result.confidence,
        evidence        = evidence_strs,
        warnings        = result.warnings,
        recommendations = result.recommendations,
        metrics         = metrics,
    )
