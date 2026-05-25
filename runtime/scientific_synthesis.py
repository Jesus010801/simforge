"""Scientific synthesis engine — multi-observable reasoning for comparative MD studies.

Orchestrates: normalization → rule evaluation → consensus → events → narrative.

Entry point: synthesize_study(study) → SynthesisResult
"""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from runtime.study_models import Study

from runtime.synthesis_models import (
    SynthesisResult, SystemSynthesis, ConsensusResult, TemporalEvent,
)


# ─── Narrative templates ──────────────────────────────────────────────────────

_STATE_NARRATIVE = {
    "stable_binding":            "stable, well-maintained binding",
    "interaction_persistent":    "consistent interaction persistence",
    "weak_binding":              "elevated ligand mobility with maintained contacts",
    "transient_binding":         "transient, inconsistently maintained interaction",
    "ligand_destabilization":    "progressive loss of ligand stability",
    "possible_dissociation":     "evidence consistent with ligand dissociation",
    "structurally_stable":       "well-equilibrated structural stability",
    "flexible_but_stable":       "elevated flexibility without unfolding",
    "conformational_rearrangement": "significant structural reorganization",
    "uncertain_behavior":        "mixed or insufficient signals for confident classification",
}


def _state_label(state: str) -> str:
    return _STATE_NARRATIVE.get(state, state.replace("_", " "))


def _format_val(obs: str, mean: float, study: "Study") -> str:
    units = study.observable_units.get(obs, "")
    return f"{mean:.3g}{' ' + units if units else ''}"


def _build_explanation(syn: SystemSynthesis, study: "Study") -> str:
    """Generate a per-system explanation sentence."""
    sg = study.systems.get(syn.system)
    if not sg:
        return syn.explanation

    # Primary state label
    state_label = _state_label(syn.primary_state)

    # Build supporting evidence phrase
    top_evidence = syn.evidence[:3]
    if top_evidence:
        ev_str = "; ".join(top_evidence)
    else:
        ev_str = "limited observable coverage"

    return (
        f"{syn.system} classified as {state_label} "
        f"(confidence: {syn.primary_confidence:.2f}). "
        f"Key evidence: {ev_str}."
    )


def _generate_narrative(
    result: SynthesisResult,
    study:  "Study",
) -> str:
    """Produce a multi-paragraph scientific narrative from synthesis results."""
    if not result.systems:
        return ""

    ranked = result.ranking  # already (system, score) sorted descending

    paragraphs: list[str] = []

    # ── Paragraph 1: top system summary ──────────────────────────────────────
    if ranked:
        top_sys, top_score = ranked[0]
        top_syn = result.systems.get(top_sys)
        if top_syn:
            state_label = _state_label(top_syn.primary_state)
            ev_strs = top_syn.evidence[:3]
            ev_part = (
                ", ".join(ev_strs) if ev_strs else "limited observable data"
            )
            paragraphs.append(
                f"{top_sys} systems showed the most favorable interaction profile "
                f"(composite score: {top_score:.2f}), characterized by {state_label}. "
                f"Supporting observations: {ev_part}."
            )

    # ── Paragraph 2: comparative description ─────────────────────────────────
    comparison_parts: list[str] = []
    for sys_name, score in ranked[1:]:
        syn = result.systems.get(sys_name)
        if not syn:
            continue
        state_label = _state_label(syn.primary_state)
        diff = ranked[0][1] - score if ranked else 0.0
        gap_desc = (
            "slightly lower" if diff < 0.10
            else "notably lower" if diff < 0.25
            else "substantially lower"
        )
        ev_strs = syn.evidence[:2]
        ev_part = (", ".join(ev_strs) if ev_strs else "limited observable data")
        comparison_parts.append(
            f"{sys_name} exhibited {state_label} ({gap_desc} score: {score:.2f}), "
            f"with {ev_part}"
        )

    if comparison_parts:
        paragraphs.append(". ".join(comparison_parts) + ".")

    # ── Paragraph 3: consensus and uncertainty ────────────────────────────────
    conflict_sys = [
        r for r in result.consensus
        if r.label in ("conflicting", "weak")
    ]
    if conflict_sys:
        conflict_strs = [
            f"{r.system} ({r.observable.replace('_', ' ')}: {r.n_agreeing}/{r.n_total} replicas agree)"
            for r in conflict_sys[:3]
        ]
        paragraphs.append(
            f"Replica consensus was inconsistent for: {'; '.join(conflict_strs)}. "
            "These systems warrant additional replicas or extended simulations."
        )

    # ── Paragraph 4: time-resolved events ────────────────────────────────────
    if result.events:
        event_strs = [e.description for e in result.events[:3]]
        paragraphs.append(
            "Time-resolved analysis detected: "
            + "; ".join(event_strs) + "."
        )

    # ── Paragraph 5: outlier note ─────────────────────────────────────────────
    if study.summary and study.summary.outlier_replicas:
        out_strs = [
            f"{s} {r}"
            for s, r, _ in study.summary.outlier_replicas[:3]
        ]
        paragraphs.append(
            f"Potential outlier replicas detected: {', '.join(out_strs)}. "
            "These should be inspected visually before inclusion in downstream analyses."
        )

    return "\n\n".join(paragraphs)


# ─── Main entry point ────────────────────────────────────────────────────────

def synthesize_study(study: "Study") -> SynthesisResult:
    """Run the full scientific synthesis pipeline on a Study object.

    Returns a SynthesisResult with system interpretations, consensus,
    time-resolved events, a composite ranking, and a narrative summary.
    """
    from runtime.interaction_interpreter import interpret_all
    from runtime.consensus_engine import evaluate_consensus, consensus_multiplier_map
    from runtime.event_detector import detect_all_events

    if not study.systems or not study.observables_detected:
        return SynthesisResult()

    # Step 1 — Consensus (used to modulate rule confidence)
    consensus      = evaluate_consensus(study)
    consensus_mmap = consensus_multiplier_map(study)

    # Step 2 — System interpretation (rules + composite scores)
    system_synths  = interpret_all(study, consensus_mmap)

    # Update explanations with the finalized synthesis
    for sys_name, syn in system_synths.items():
        syn.explanation = _build_explanation(syn, study)

    # Step 3 — Time-resolved events
    events = detect_all_events(study)

    # Step 4 — Ranking by composite score
    ranking = sorted(
        [(sn, syn.composite_score) for sn, syn in system_synths.items()],
        key=lambda x: -x[1],
    )

    # Step 5 — Narrative
    result = SynthesisResult(
        systems   = system_synths,
        consensus = consensus,
        events    = events,
        ranking   = ranking,
    )
    result.narrative = _generate_narrative(result, study)

    return result
