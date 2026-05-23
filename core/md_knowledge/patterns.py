"""core/md_knowledge/patterns.py — Temporal pattern detection for MD observables."""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from typing import Sequence


class TemporalPattern(str, Enum):
    PLATEAU        = "plateau"          # stable convergence
    DRIFT          = "drift"            # monotonic trend, no plateau
    OSCILLATING    = "oscillating"      # periodic variation, no net drift
    JUMP_PLATEAU   = "jump_plateau"     # step change then new plateau
    NOISY          = "noisy"            # high variance, no clear pattern
    INSUFFICIENT   = "insufficient"     # too few points to classify


@dataclass
class PatternResult:
    pattern:    TemporalPattern
    confidence: float               # [0, 1]
    slope:      float               # units/ns (linear regression slope)
    plateau_std: float              # std of last 20% of values
    mean_last:  float               # mean of last 20%
    mean_first: float               # mean of first 20%
    n_points:   int
    details:    dict = field(default_factory=dict)


# ── Internal helpers ──────────────────────────────────────────────────────────

def _linreg(x: list[float], y: list[float]) -> tuple[float, float]:
    """Return (slope, intercept) via least squares."""
    n = len(x)
    if n < 2:
        return 0.0, y[0] if y else 0.0
    sx  = sum(x);  sy  = sum(y)
    sxx = sum(xi * xi for xi in x)
    sxy = sum(xi * yi for xi, yi in zip(x, y))
    denom = n * sxx - sx * sx
    if abs(denom) < 1e-12:
        return 0.0, sy / n
    slope = (n * sxy - sx * sy) / denom
    intercept = (sy - slope * sx) / n
    return slope, intercept


def _std(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    m = sum(values) / len(values)
    return math.sqrt(sum((v - m) ** 2 for v in values) / (len(values) - 1))


def _mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def _r_squared(x: list[float], y: list[float], slope: float, intercept: float) -> float:
    """R² of a linear fit."""
    if len(y) < 2:
        return 0.0
    y_mean = sum(y) / len(y)
    ss_tot = sum((yi - y_mean) ** 2 for yi in y)
    ss_res = sum((yi - (slope * xi + intercept)) ** 2 for xi, yi in zip(x, y))
    return 1.0 - ss_res / ss_tot if abs(ss_tot) > 1e-12 else 0.0


def _autocorr_lag1(values: list[float]) -> float:
    """Lag-1 autocorrelation — positive near 1 → trend/plateau; near 0 → noise."""
    if len(values) < 3:
        return 0.0
    m = _mean(values)
    numer = sum((values[i] - m) * (values[i - 1] - m) for i in range(1, len(values)))
    denom = sum((v - m) ** 2 for v in values)
    return numer / denom if abs(denom) > 1e-12 else 0.0


# ── Public API ────────────────────────────────────────────────────────────────

MIN_POINTS = 10  # fewer → INSUFFICIENT


def detect_temporal_pattern(
    values: Sequence[float],
    times_ns: Sequence[float] | None = None,
    *,
    plateau_slope_threshold: float = 0.05,   # units/ns; below → plateau candidate
    drift_slope_threshold:   float = 0.01,   # units/ns; above → drift
    noise_cv_threshold:      float = 0.30,   # CV of last 20%; above → noisy
    jump_delta_sigma:        float = 3.0,    # jump ≥ N×std of first 80%
) -> PatternResult:
    """
    Classify the temporal behaviour of a scalar observable.

    Parameters
    ----------
    values      : observable values (same order as simulation frames)
    times_ns    : corresponding times in nanoseconds (optional; indices used if absent)
    plateau_slope_threshold : |slope| below this → plateau candidate
    drift_slope_threshold   : |slope| above this → drift
    noise_cv_threshold      : coefficient of variation threshold for noisy classification
    jump_delta_sigma        : jump detection sensitivity (multiples of std)
    """
    vals = list(values)
    n = len(vals)

    if n < MIN_POINTS:
        return PatternResult(
            pattern=TemporalPattern.INSUFFICIENT,
            confidence=1.0,
            slope=0.0,
            plateau_std=0.0,
            mean_last=_mean(vals),
            mean_first=_mean(vals),
            n_points=n,
        )

    if times_ns is not None:
        ts = list(times_ns)
    else:
        ts = list(range(n))

    # Segment indices
    tail_start = max(1, int(n * 0.80))
    head_end   = max(1, int(n * 0.20))

    last_vals  = vals[tail_start:]
    first_vals = vals[:head_end]

    mean_last  = _mean(last_vals)
    mean_first = _mean(first_vals)
    std_last   = _std(last_vals)
    std_first  = _std(first_vals)
    std_all    = _std(vals)

    slope, intercept = _linreg(ts, vals)
    r2 = _r_squared(ts, vals, slope, intercept)
    ac = _autocorr_lag1(vals)

    # ── NONPHYSICAL guard (delegated to caller via extreme values) ────────────

    # ── DRIFT: linear signal check (before jump, to avoid misclassifying ramps) ─
    if r2 >= 0.90 and abs(slope) >= drift_slope_threshold:
        confidence = min(1.0, r2 * (abs(slope) / (drift_slope_threshold * 5 + 1e-9)))
        return PatternResult(
            pattern=TemporalPattern.DRIFT,
            confidence=min(1.0, confidence),
            slope=slope,
            plateau_std=std_last,
            mean_last=mean_last,
            mean_first=mean_first,
            n_points=n,
            details={"r2": r2},
        )

    # ── JUMP detection ────────────────────────────────────────────────────────
    jump_threshold = jump_delta_sigma * (std_first + 1e-9)
    delta = abs(mean_last - mean_first)

    if delta > jump_threshold and std_last < delta * 0.4:
        # Large offset between first/last segments AND last segment is tight
        confidence = min(1.0, delta / jump_threshold / 2)
        return PatternResult(
            pattern=TemporalPattern.JUMP_PLATEAU,
            confidence=confidence,
            slope=slope,
            plateau_std=std_last,
            mean_last=mean_last,
            mean_first=mean_first,
            n_points=n,
            details={"delta": delta, "jump_threshold": jump_threshold},
        )

    # ── PLATEAU check ─────────────────────────────────────────────────────────
    cv_last = std_last / (abs(mean_last) + 1e-9)
    rel_slope = abs(slope) / (abs(mean_last) + 1e-9)  # relative slope

    if abs(slope) <= plateau_slope_threshold and cv_last < noise_cv_threshold:
        confidence = 1.0 - min(1.0, rel_slope / (plateau_slope_threshold + 1e-9) * 0.5)
        return PatternResult(
            pattern=TemporalPattern.PLATEAU,
            confidence=min(1.0, confidence),
            slope=slope,
            plateau_std=std_last,
            mean_last=mean_last,
            mean_first=mean_first,
            n_points=n,
        )

    # ── DRIFT check (lower R², but still trending) ───────────────────────────
    if abs(slope) >= drift_slope_threshold and ac > 0.6:
        confidence = min(1.0, abs(slope) / (drift_slope_threshold * 5))
        return PatternResult(
            pattern=TemporalPattern.DRIFT,
            confidence=confidence,
            slope=slope,
            plateau_std=std_last,
            mean_last=mean_last,
            mean_first=mean_first,
            n_points=n,
        )

    # ── OSCILLATING check ─────────────────────────────────────────────────────
    # Low autocorrelation AND slope ~0 AND moderate variance
    if ac < 0.3 and abs(slope) < drift_slope_threshold * 2:
        zero_crossings = sum(
            1 for i in range(1, len(vals))
            if (vals[i] - mean_last) * (vals[i - 1] - mean_last) < 0
        )
        if zero_crossings > n * 0.1:
            return PatternResult(
                pattern=TemporalPattern.OSCILLATING,
                confidence=0.6,
                slope=slope,
                plateau_std=std_last,
                mean_last=mean_last,
                mean_first=mean_first,
                n_points=n,
                details={"zero_crossings": zero_crossings},
            )

    # ── NOISY fallback ────────────────────────────────────────────────────────
    if cv_last >= noise_cv_threshold or (std_all / (abs(mean_last) + 1e-9)) > noise_cv_threshold:
        return PatternResult(
            pattern=TemporalPattern.NOISY,
            confidence=0.5,
            slope=slope,
            plateau_std=std_last,
            mean_last=mean_last,
            mean_first=mean_first,
            n_points=n,
        )

    # Default: weak plateau
    return PatternResult(
        pattern=TemporalPattern.PLATEAU,
        confidence=0.4,
        slope=slope,
        plateau_std=std_last,
        mean_last=mean_last,
        mean_first=mean_first,
        n_points=n,
    )
