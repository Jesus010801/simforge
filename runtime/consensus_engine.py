"""Replica consensus engine.

For each (system, observable), computes how consistently the replicas agree
on the direction of the observable relative to the study-wide median.

Consensus labels:
  full        — all replicas on same side of study median
  strong      — ≥ 80 % agree
  moderate    — ≥ 60 % agree
  weak        — ≥ 50 % agree
  conflicting — < 50 % agree (random behavior)
"""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from runtime.study_models import Study

from runtime.synthesis_models import ConsensusResult


def _mean(vals: list[float]) -> float:
    return sum(vals) / len(vals) if vals else 0.0


def _median(vals: list[float]) -> float:
    if not vals:
        return 0.0
    s = sorted(vals)
    n = len(s)
    mid = n // 2
    return s[mid] if n % 2 else (s[mid - 1] + s[mid]) / 2.0


def _label_and_multiplier(frac: float, n: int) -> tuple[str, float]:
    if n < 2:
        return "insufficient_data", 0.5
    if frac >= 1.0:
        return "full",       1.00
    if frac >= 0.80:
        return "strong",     0.90
    if frac >= 0.60:
        return "moderate",   0.75
    if frac >= 0.50:
        return "weak",       0.60
    return "conflicting",    0.45


def evaluate_consensus(study: "Study") -> list[ConsensusResult]:
    """Return notable consensus results across all (system, observable) pairs."""
    results: list[ConsensusResult] = []

    for obs in study.observables_detected:
        # Study-wide median for this observable (all replica means pooled)
        all_replica_means: list[float] = []
        for sg in study.systems.values():
            if obs not in sg.aggregate:
                continue
            for replica in sg.replicas.values():
                if obs in replica.observables:
                    all_replica_means.append(replica.observables[obs].mean)

        if not all_replica_means:
            continue
        study_median = _median(all_replica_means)

        for sys_name, sg in study.systems.items():
            if obs not in sg.aggregate:
                continue

            replica_means = [
                r.observables[obs].mean
                for r in sg.replicas.values()
                if obs in r.observables
            ]
            n = len(replica_means)
            if n < 2:
                continue

            n_above = sum(1 for v in replica_means if v > study_median)
            n_below = n - n_above
            n_agreeing = max(n_above, n_below)
            frac = n_agreeing / n

            label, mult = _label_and_multiplier(frac, n)

            # Only include non-trivial consensus results
            if label not in ("full",):
                results.append(ConsensusResult(
                    system     = sys_name,
                    observable = obs,
                    n_total    = n,
                    n_agreeing = n_agreeing,
                    label      = label,
                    multiplier = mult,
                ))

    # Sort: most problematic (conflicting/weak) first
    _order = {"conflicting": 0, "weak": 1, "moderate": 2, "strong": 3, "full": 4}
    results.sort(key=lambda r: _order.get(r.label, 5))
    return results


def consensus_multiplier_map(study: "Study") -> dict[str, float]:
    """Return per-system average consensus multiplier (for confidence scaling)."""
    raw = evaluate_consensus(study)
    sys_mults: dict[str, list[float]] = {}
    for r in raw:
        sys_mults.setdefault(r.system, []).append(r.multiplier)
    return {
        sys: sum(ms) / len(ms)
        for sys, ms in sys_mults.items()
    }
