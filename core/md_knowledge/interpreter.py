"""core/md_knowledge/interpreter.py — Context-aware MD simulation interpreter."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Sequence

from core.md_knowledge.states import SimulationState, STATE_DESCRIPTIONS, STATE_SEVERITY
from core.md_knowledge.patterns import detect_temporal_pattern, PatternResult, TemporalPattern
from core.md_knowledge.contexts import SystemContext, SYSTEM_CONTEXTS
from core.md_knowledge.evidence import Evidence, EvidenceBundle, accumulate_evidence
from core.md_knowledge.heuristics import get_heuristic


@dataclass
class ObservableResult:
    """Parsed data for a single observable to be interpreted."""
    name:    str
    values:  list[float]
    times_ns: list[float] | None = None


@dataclass
class InterpretationResult:
    """Full scientific interpretation of a simulation run."""
    state:          SimulationState
    confidence:     float                    # [0, 1]
    context:        SystemContext
    evidence:       EvidenceBundle
    summary:        str                      # 1-2 sentence human-readable verdict
    warnings:       list[str] = field(default_factory=list)
    recommendations: list[str] = field(default_factory=list)
    pattern_results: dict[str, PatternResult] = field(default_factory=dict)

    @property
    def state_description(self) -> str:
        return STATE_DESCRIPTIONS.get(self.state, "")

    @property
    def quality_tier(self) -> str:
        """Map state+confidence to a 5-tier quality label."""
        sev = STATE_SEVERITY.get(self.state, 0)
        if sev >= 4:
            return "FAILED"
        if sev >= 3:
            return "POOR"
        if self.confidence >= 0.75 and sev == 0:
            return "EXCELLENT"
        if self.confidence >= 0.55 and sev <= 1:
            return "GOOD"
        return "ACCEPTABLE"


# ── Default pattern detection thresholds (context-aware overrides possible) ───

_PLATEAU_SLOPE_BY_OBSERVABLE: dict[str, float] = {
    "rmsd":                  0.005,   # nm/ns
    "rg":                    0.010,   # nm/ns
    "temperature":           0.100,   # K/ns
    "pressure":              1.000,   # bar/ns
    "potential_energy":      0.500,   # kJ/mol/atom per ns
    "rmsf":                  0.010,
    "sasa":                  0.500,   # nm²/ns
    "hbonds":                0.500,   # count/ns
    "ligand_pocket_distance": 0.005,  # nm/ns
    "secondary_structure":   0.005,   # fraction/ns
}

_DRIFT_SLOPE_BY_OBSERVABLE: dict[str, float] = {
    "rmsd":                  0.002,
    "rg":                    0.005,
    "temperature":           0.050,
    "pressure":              0.500,
    "potential_energy":      0.100,
    "rmsf":                  0.005,
    "sasa":                  0.200,
    "hbonds":                0.200,
    "ligand_pocket_distance": 0.002,
    "secondary_structure":   0.002,
}


def _build_summary(
    state: SimulationState,
    confidence: float,
    context: SystemContext,
    bundle: EvidenceBundle,
) -> str:
    """Generate a concise human-readable verdict."""
    ctx_name = SYSTEM_CONTEXTS[context].name if context in SYSTEM_CONTEXTS else context.value
    state_label = state.value.replace("_", " ").title()

    conf_label = (
        "high confidence" if confidence >= 0.75 else
        "moderate confidence" if confidence >= 0.50 else
        "low confidence"
    )

    errors   = bundle.errors
    warnings = bundle.warnings

    if errors:
        return (
            f"[{ctx_name}] {state_label} ({conf_label}). "
            f"Critical issues detected: {errors[0][:120]}."
        )
    if warnings:
        return (
            f"[{ctx_name}] {state_label} ({conf_label}). "
            f"Notable: {warnings[0][:120]}."
        )
    return (
        f"[{ctx_name}] {state_label} ({conf_label}). "
        f"No critical issues detected across {len(bundle.items)} observables."
    )


def interpret_simulation(
    observables: list[ObservableResult],
    context: SystemContext = SystemContext.UNKNOWN,
    *,
    custom_thresholds: dict[str, dict[str, float]] | None = None,
) -> InterpretationResult:
    """
    Interpret a MD simulation run from per-observable time series.

    Parameters
    ----------
    observables         : list of ObservableResult (name + values + optional times)
    context             : system context for context-aware soft ranges
    custom_thresholds   : override plateau/drift slopes per observable
                          e.g. {"rmsd": {"plateau": 0.01, "drift": 0.003}}

    Returns
    -------
    InterpretationResult with state, confidence, evidence, and recommendations.
    """
    custom_thresholds = custom_thresholds or {}

    # ── 1. Detect temporal patterns ───────────────────────────────────────────
    pattern_results: dict[str, PatternResult] = {}
    observable_patterns: dict[str, tuple[PatternResult, float | None]] = {}

    for obs in observables:
        if not obs.values:
            continue

        plateau_slope = custom_thresholds.get(obs.name, {}).get(
            "plateau", _PLATEAU_SLOPE_BY_OBSERVABLE.get(obs.name, 0.05)
        )
        drift_slope = custom_thresholds.get(obs.name, {}).get(
            "drift", _DRIFT_SLOPE_BY_OBSERVABLE.get(obs.name, 0.01)
        )

        pr = detect_temporal_pattern(
            obs.values,
            obs.times_ns,
            plateau_slope_threshold=plateau_slope,
            drift_slope_threshold=drift_slope,
        )
        pattern_results[obs.name] = pr

        # Representative value: mean of last 20% (plateau mean)
        rep_value: float | None = pr.mean_last if pr.n_points >= 5 else None

        observable_patterns[obs.name] = (pr, rep_value)

    # ── 2. Accumulate evidence ────────────────────────────────────────────────
    bundle = accumulate_evidence(observable_patterns)

    if not bundle.items:
        return InterpretationResult(
            state=SimulationState.INSUFFICIENT_SAMPLING,
            confidence=0.0,
            context=context,
            evidence=bundle,
            summary="No observables provided — cannot assess simulation quality.",
            pattern_results=pattern_results,
        )

    # ── 3. Apply context-aware soft range scoring ─────────────────────────────
    # Adjust confidence based on whether plateau values fall in expected ranges
    range_adjustments: list[float] = []
    for obs in observables:
        heuristic = get_heuristic(obs.name, context)
        if heuristic is None:
            continue
        pr = pattern_results.get(obs.name)
        if pr is None or pr.pattern == TemporalPattern.INSUFFICIENT:
            continue
        rep = pr.mean_last
        score = heuristic.soft_range.score(rep)
        range_adjustments.append(score)

    if range_adjustments:
        avg_range_score = sum(range_adjustments) / len(range_adjustments)
        # Blend pattern confidence with range score
        bundle.overall_confidence = 0.60 * bundle.overall_confidence + 0.40 * avg_range_score
        bundle.overall_confidence = max(0.0, min(1.0, bundle.overall_confidence))

    # ── 4. Determine dominant state ───────────────────────────────────────────
    state = bundle.dominant_state()
    confidence = bundle.overall_confidence

    # Promote to NONPHYSICAL if any error-severity evidence votes for it
    for ev in bundle.items:
        if ev.state_vote == SimulationState.NONPHYSICAL_BEHAVIOR and ev.severity == "error":
            state = SimulationState.NONPHYSICAL_BEHAVIOR
            break

    # ── 5. Collect warnings / recommendations ─────────────────────────────────
    warnings = list(dict.fromkeys(bundle.warnings))          # deduplicate, preserve order
    recommendations = list(dict.fromkeys(bundle.recommendations))

    # Context-specific recommendations for drifting RMSD
    rmsd_pr = pattern_results.get("rmsd")
    if rmsd_pr and rmsd_pr.pattern == TemporalPattern.DRIFT:
        if context in (SystemContext.IDR, SystemContext.PEPTIDE):
            recommendations.append(
                "Consider enhanced sampling (REST2, metadynamics) for disordered systems."
            )
        else:
            recommendations.append(
                "Extend equilibration or production. Check thermostat/barostat parameters."
            )

    # ── 6. Build summary ──────────────────────────────────────────────────────
    summary = _build_summary(state, confidence, context, bundle)

    return InterpretationResult(
        state=state,
        confidence=confidence,
        context=context,
        evidence=bundle,
        summary=summary,
        warnings=warnings,
        recommendations=recommendations,
        pattern_results=pattern_results,
    )
