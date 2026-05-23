"""core/md_knowledge/evidence.py — Multi-observable evidence accumulation."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from core.md_knowledge.states import SimulationState
from core.md_knowledge.patterns import TemporalPattern, PatternResult


@dataclass
class Evidence:
    """A single piece of interpretive evidence from one observable."""
    observable:    str
    value:         float | None          # representative scalar (plateau mean)
    pattern:       TemporalPattern
    state_vote:    SimulationState       # which state this evidence supports
    confidence:    float                 # [0, 1] — how certain this evidence is
    message:       str                   # human-readable interpretation
    severity:      str = "info"          # "info" | "warning" | "error"
    recommendation: str = ""


@dataclass
class EvidenceBundle:
    """Collection of evidence pieces for one simulation run."""
    items: list[Evidence] = field(default_factory=list)

    # Tallies updated by accumulate_evidence()
    state_votes:   dict[str, float] = field(default_factory=dict)   # state → weighted vote sum
    overall_confidence: float = 0.0

    def add(self, evidence: Evidence) -> None:
        self.items.append(evidence)
        key = evidence.state_vote.value
        self.state_votes[key] = self.state_votes.get(key, 0.0) + evidence.confidence

    @property
    def warnings(self) -> list[str]:
        return [e.message for e in self.items if e.severity in ("warning", "error")]

    @property
    def errors(self) -> list[str]:
        return [e.message for e in self.items if e.severity == "error"]

    @property
    def recommendations(self) -> list[str]:
        return [e.recommendation for e in self.items if e.recommendation]

    def dominant_state(self) -> SimulationState:
        if not self.state_votes:
            return SimulationState.INSUFFICIENT_SAMPLING
        top = max(self.state_votes, key=lambda k: self.state_votes[k])
        return SimulationState(top)

    def compute_confidence(self) -> float:
        """Aggregate confidence: weighted average across items."""
        if not self.items:
            return 0.0
        total_weight = sum(e.confidence for e in self.items)
        if total_weight < 1e-9:
            return 0.0
        # Errors pull confidence down
        error_penalty = sum(0.15 for e in self.items if e.severity == "error")
        raw = total_weight / len(self.items)
        return max(0.0, min(1.0, raw - error_penalty))


# ── Pattern → state mapping ───────────────────────────────────────────────────

_PATTERN_STATE: dict[TemporalPattern, SimulationState] = {
    TemporalPattern.PLATEAU:      SimulationState.STABLE_EQUILIBRATED,
    TemporalPattern.DRIFT:        SimulationState.DRIFTING,
    TemporalPattern.OSCILLATING:  SimulationState.METASTABLE,
    TemporalPattern.JUMP_PLATEAU: SimulationState.CONFORMATIONAL_TRANSITION,
    TemporalPattern.NOISY:        SimulationState.UNSTABLE,
    TemporalPattern.INSUFFICIENT: SimulationState.INSUFFICIENT_SAMPLING,
}

# Observables that are always expected to plateau (pattern matters most)
_CONVERGENCE_SENSITIVE = {"rmsd", "rg", "potential_energy", "temperature", "pressure"}
# Observables where absolute value matters more
_VALUE_SENSITIVE = {"ligand_pocket_distance", "secondary_structure"}


def accumulate_evidence(
    observable_patterns: dict[str, tuple[PatternResult, float | None]],
    *,
    heuristics: dict | None = None,
) -> EvidenceBundle:
    """
    Build an EvidenceBundle from per-observable PatternResults.

    Parameters
    ----------
    observable_patterns : mapping from observable name →
                          (PatternResult, representative_value | None)
    heuristics          : optional pre-fetched heuristic dict (for testability)
    """
    bundle = EvidenceBundle()

    for obs_name, (pattern_result, rep_value) in observable_patterns.items():
        # Determine base state vote from pattern
        state_vote = _PATTERN_STATE.get(pattern_result.pattern, SimulationState.STABLE_EQUILIBRATED)
        confidence = pattern_result.confidence

        # Override based on absolute value for critical observables
        severity = "info"
        message_parts = [
            f"{obs_name}: {pattern_result.pattern.value} "
            f"(mean={pattern_result.mean_last:.3g}, std={pattern_result.plateau_std:.3g})"
        ]
        recommendation = ""

        if obs_name == "ligand_pocket_distance" and rep_value is not None:
            if rep_value > 0.80:
                state_vote  = SimulationState.LIGAND_DISSOCIATED
                severity    = "error"
                message_parts.append(f"Ligand has dissociated (dist={rep_value:.2f} nm)")
                recommendation = "Check force-field parameters and binding pose preparation."
            elif rep_value > 0.50:
                state_vote  = SimulationState.PARTIALLY_CONVERGED
                severity    = "warning"
                message_parts.append(f"Ligand at pocket periphery (dist={rep_value:.2f} nm)")

        if obs_name == "secondary_structure" and rep_value is not None:
            if rep_value < 0.50:
                state_vote  = SimulationState.UNSTABLE
                severity    = "error"
                message_parts.append(f"Severe secondary structure loss ({rep_value:.0%} retained)")
                recommendation = "Verify force field, pH, and disulfide bonds."
            elif rep_value < 0.75:
                state_vote  = SimulationState.DRIFTING
                severity    = "warning"
                message_parts.append(f"Significant SS loss ({rep_value:.0%} retained)")

        # Temperature / pressure — catch physically implausible values
        if obs_name == "temperature" and rep_value is not None:
            if rep_value > 400.0:
                state_vote  = SimulationState.NONPHYSICAL_BEHAVIOR
                severity    = "error"
                message_parts.append(f"Physically implausible temperature: {rep_value:.0f} K")
                recommendation = "Simulation is likely exploding. Restart with better minimization."
                confidence = 1.0

        if obs_name == "potential_energy" and pattern_result.pattern == TemporalPattern.NOISY:
            if abs(pattern_result.slope) > 10.0:   # kJ/mol/atom per ns
                state_vote  = SimulationState.NONPHYSICAL_BEHAVIOR
                severity    = "error"
                message_parts.append("Energy spike detected — numerical instability")
                recommendation = "Reduce timestep. Check for clashes."

        ev = Evidence(
            observable=obs_name,
            value=rep_value,
            pattern=pattern_result.pattern,
            state_vote=state_vote,
            confidence=confidence,
            message=" | ".join(message_parts),
            severity=severity,
            recommendation=recommendation,
        )
        bundle.add(ev)

    bundle.overall_confidence = bundle.compute_confidence()
    return bundle
