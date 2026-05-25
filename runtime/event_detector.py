"""Time-resolved event detection for MD trajectories.

Scans individual ObservableSeries for patterns that indicate dynamic events:
  late_destabilization  — observable drifts significantly in the second half
  abrupt_transition     — sudden large jump relative to overall variance
  contact_loss          — sustained reduction in contacts during the simulation
  ligand_drift          — monotonic increase in ligand RMSD suggesting departure
"""
from __future__ import annotations

from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from runtime.study_models import Study, ObservableSeries

from runtime.synthesis_models import TemporalEvent


# ─── Helpers ─────────────────────────────────────────────────────────────────

def _mean(vals: list[float]) -> float:
    return sum(vals) / len(vals) if vals else 0.0


def _std(vals: list[float]) -> float:
    if len(vals) < 2:
        return 0.0
    m = _mean(vals)
    return (sum((v - m) ** 2 for v in vals) / (len(vals) - 1)) ** 0.5


def _linear_slope(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    if n < 2:
        return 0.0
    sx = sum(xs); sy = sum(ys)
    sxy = sum(x * y for x, y in zip(xs, ys))
    sxx = sum(x * x for x in xs)
    denom = n * sxx - sx * sx
    return (n * sxy - sx * sy) / denom if abs(denom) > 1e-12 else 0.0


# ─── Lower-is-better observables ─────────────────────────────────────────────

_LOWER_IS_BETTER = frozenset({
    "protein_rmsd", "ligand_rmsd", "rmsf", "mindist", "catalytic_distance", "distance",
})
_HIGHER_IS_BETTER = frozenset({"contacts", "hydrogen_bonds"})


# ─── Event detectors ─────────────────────────────────────────────────────────

_MIN_POINTS = 20   # minimum data points needed for any event detection


def _detect_late_destabilization(series: "ObservableSeries") -> Optional[TemporalEvent]:
    """Detect significant drift in the second half of the trajectory."""
    vals  = series.values
    times = series.time_ns
    n = len(vals)
    if n < _MIN_POINTS:
        return None

    mid = n // 2
    first_half  = vals[:mid]
    second_half = vals[mid:]
    t_mid       = times[mid] if mid < len(times) else 0.0

    mean_first  = _mean(first_half)
    mean_second = _mean(second_half)

    if abs(mean_first) < 1e-9:
        return None

    rel_change = (mean_second - mean_first) / abs(mean_first)

    obs = series.observable
    # For lower-is-better: bad if second half increases significantly
    if obs in _LOWER_IS_BETTER and rel_change > 0.25:
        disp = obs.replace("_", " ")
        return TemporalEvent(
            system     = series.system,
            replica    = series.replica,
            observable = obs,
            event_type = "late_destabilization",
            time_ns    = t_mid,
            description = (
                f"{disp} increases by {rel_change:.0%} in second half "
                f"(from {mean_first:.3g} to {mean_second:.3g}) — possible late destabilization"
            ),
        )
    # For higher-is-better: bad if second half decreases significantly
    if obs in _HIGHER_IS_BETTER and rel_change < -0.25:
        disp = obs.replace("_", " ")
        return TemporalEvent(
            system     = series.system,
            replica    = series.replica,
            observable = obs,
            event_type = "contact_loss",
            time_ns    = t_mid,
            description = (
                f"{disp} drops by {abs(rel_change):.0%} in second half "
                f"(from {mean_first:.3g} to {mean_second:.3g}) — possible interaction loss"
            ),
        )
    return None


def _detect_abrupt_transition(series: "ObservableSeries") -> Optional[TemporalEvent]:
    """Detect a sudden large change (> 3σ) using sliding window means."""
    vals  = series.values
    times = series.time_ns
    n = len(vals)
    if n < _MIN_POINTS * 2:
        return None

    overall_std = _std(vals)
    if overall_std < 1e-10:
        return None

    win = max(5, n // 15)
    window_means: list[float] = []
    window_times: list[float] = []

    for i in range(0, n - win, max(1, win // 2)):
        window_means.append(_mean(vals[i: i + win]))
        mid_t = times[i + win // 2] if (i + win // 2) < len(times) else 0.0
        window_times.append(mid_t)

    if len(window_means) < 3:
        return None

    max_delta = 0.0
    max_t     = 0.0
    for i in range(len(window_means) - 1):
        delta = abs(window_means[i + 1] - window_means[i])
        if delta > max_delta:
            max_delta = delta
            max_t     = window_times[i]

    if max_delta > 3.0 * overall_std:
        disp = series.observable.replace("_", " ")
        return TemporalEvent(
            system      = series.system,
            replica     = series.replica,
            observable  = series.observable,
            event_type  = "abrupt_transition",
            time_ns     = max_t,
            description = (
                f"Abrupt {disp} change (Δ={max_delta:.3g}, "
                f"{max_delta / overall_std:.1f}σ) detected at ~{max_t:.1f} ns"
            ),
        )
    return None


def _detect_ligand_drift(series: "ObservableSeries") -> Optional[TemporalEvent]:
    """Detect monotonically increasing ligand RMSD (ligand departure signal)."""
    if series.observable != "ligand_rmsd":
        return None
    vals  = series.values
    times = series.time_ns
    n = len(vals)
    if n < _MIN_POINTS:
        return None

    # Split into thirds; monotonically increasing slope is the signal
    third = n // 3
    mean_early  = _mean(vals[:third])
    mean_middle = _mean(vals[third: 2 * third])
    mean_late   = _mean(vals[2 * third:])

    if mean_early < 1e-9:
        return None

    increasing = mean_early < mean_middle < mean_late
    total_rise  = (mean_late - mean_early) / mean_early

    if increasing and total_rise > 0.30:
        return TemporalEvent(
            system      = series.system,
            replica     = series.replica,
            observable  = "ligand_rmsd",
            event_type  = "ligand_drift",
            time_ns     = times[2 * third] if 2 * third < len(times) else 0.0,
            description = (
                f"Ligand RMSD shows monotonic increase "
                f"({mean_early:.3g}→{mean_middle:.3g}→{mean_late:.3g}) — "
                f"possible ligand departure"
            ),
        )
    return None


# ─── Main entry point ────────────────────────────────────────────────────────

def detect_all_events(study: "Study") -> list[TemporalEvent]:
    """Scan all series in a study for time-resolved events."""
    events: list[TemporalEvent] = []
    seen: set[tuple[str, str, str, str]] = set()  # dedup by (system, replica, obs, type)

    for sg in study.systems.values():
        for replica in sg.replicas.values():
            for obs_name, series in replica.observables.items():
                for detector in (
                    _detect_late_destabilization,
                    _detect_abrupt_transition,
                    _detect_ligand_drift,
                ):
                    event = detector(series)
                    if event is None:
                        continue
                    key = (event.system, event.replica, event.observable, event.event_type)
                    if key not in seen:
                        seen.add(key)
                        events.append(event)

    # Sort by severity: abrupt_transition first, then late_destabilization, others
    _severity = {
        "abrupt_transition":    0,
        "late_destabilization": 1,
        "contact_loss":         2,
        "ligand_drift":         3,
    }
    events.sort(key=lambda e: (_severity.get(e.event_type, 9), e.system, e.replica))
    return events
