"""
core/semantic_inference.py — Semantic normalization and scientific intent inference.

Stage 0.5 of the parse pipeline — runs after YAML load, before run_inference().

Responsibilities:
  1. Normalize simulation_objectives (aliases → canonical)
  2. Apply simulation_profile presets
  3. Detect membrane protein signals and auto-inject objectives / hints
  4. Emit structured NormalizationWarning into global_reasoning
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from core.models import SystemState


# ── Normalization warning (lightweight, no Pydantic) ──────────────────────────

@dataclass
class NormalizationNote:
    kind:    str   # "alias" | "fuzzy" | "unknown" | "preset" | "membrane_inferred"
    original: str
    resolved: list[str]
    message:  str


# ── Membrane protein signal detection ─────────────────────────────────────────

def detect_membrane_protein_signals(state: "SystemState") -> list[str]:
    """
    Inspect SystemState for evidence that this is a membrane protein system.

    Returns a list of human-readable signal strings (non-empty = membrane protein
    context detected).
    """
    signals: list[str] = []

    for comp in state.components:
        ctx = getattr(comp, "biological_context", None) or []
        if "transmembrane" in ctx:
            signals.append(f"{comp.id}: biological_context=transmembrane")
        if "membrane_associated" in ctx:
            signals.append(f"{comp.id}: biological_context=membrane_associated")

    if getattr(state.environment, "membrane", None) and state.environment.membrane.enabled:
        signals.append("environment.membrane.enabled=true")

    # Lipid components explicitly declared
    lipid_comps = [c for c in state.components if c.role == "lipid"]
    if lipid_comps:
        signals.append(f"lipid component(s): {', '.join(c.id for c in lipid_comps)}")

    # Membrane-related objectives already present
    mem_objectives = {"membrane_perturbation", "membrane_insertion", "permeability"}
    current = set(state.simulation_objectives)
    if current & mem_objectives:
        signals.append(f"membrane objective(s): {current & mem_objectives}")

    return signals


# ── Objective normalization ────────────────────────────────────────────────────

def _normalize_objectives_list(
    raw_objectives: list[str],
) -> tuple[list[str], list[NormalizationNote]]:
    """
    Normalize a list of raw objective strings to canonical names.

    Returns (normalized_list, notes).
    Unknowns are retained as-is so the caller can surface suggestions.
    """
    from core.semantic_objectives import normalize_objective

    normalized: list[str] = []
    notes: list[NormalizationNote] = []

    for obj in raw_objectives:
        canonicals, note_text = normalize_objective(obj)

        if canonicals:
            for c in canonicals:
                if c not in normalized:
                    normalized.append(c)
            if note_text:
                notes.append(NormalizationNote(
                    kind="fuzzy" if "fuzzy" in (note_text or "") else "alias",
                    original=obj,
                    resolved=canonicals,
                    message=note_text,
                ))
        else:
            # Unknown — keep original so the caller can surface suggestions
            if obj not in normalized:
                normalized.append(obj)
            notes.append(NormalizationNote(
                kind="unknown",
                original=obj,
                resolved=[],
                message=f"Unknown objective: '{obj}'",
            ))

    return normalized, notes


# ── Preset application ─────────────────────────────────────────────────────────

def apply_preset(state: "SystemState", profile_name: str) -> tuple["SystemState", NormalizationNote | None]:
    """
    Expand a simulation_profile name into objectives + workflow hints.

    Returns (modified_state, note_or_None).
    """
    from core.semantic_objectives import SIMULATION_PRESETS, suggest_objectives

    key = profile_name.lower().strip().replace(" ", "_").replace("-", "_")
    preset = SIMULATION_PRESETS.get(key)

    if preset is None:
        # Suggest closest preset names
        import difflib
        matches = difflib.get_close_matches(key, SIMULATION_PRESETS.keys(), n=3, cutoff=0.4)
        note = NormalizationNote(
            kind="unknown",
            original=profile_name,
            resolved=[],
            message=(
                f"Unknown simulation_profile: '{profile_name}'."
                + (f"  Closest: {matches}" if matches else "")
            ),
        )
        return state, note

    # Inject preset objectives (deduplicate)
    for obj in preset["objectives"]:
        if obj not in state.simulation_objectives:
            state.simulation_objectives.append(obj)

    # Populate structured WorkflowHints — replaces text storage in global_reasoning.notes.
    hints = preset.get("hints", {})
    if hints:
        for field_name, value in hints.items():
            if hasattr(state.workflow_hints, field_name):
                setattr(state.workflow_hints, field_name, value)
        # Keep a minimal audit trace (key names only, not values)
        state.global_reasoning.notes.append(
            f"[semantic] preset:{key} → hints={sorted(hints.keys())}"
        )

    note = NormalizationNote(
        kind="preset",
        original=profile_name,
        resolved=preset["objectives"],
        message=f"simulation_profile '{profile_name}' → {preset['objectives']}",
    )
    return state, note


# ── Membrane protein auto-inference ───────────────────────────────────────────

def _inject_membrane_objectives(state: "SystemState") -> list[NormalizationNote]:
    """
    If membrane protein signals are detected AND no membrane objectives are present,
    auto-inject membrane_perturbation + stability.
    """
    signals = detect_membrane_protein_signals(state)
    if not signals:
        return []

    mem_objectives = {"membrane_perturbation", "stability"}
    current = set(state.simulation_objectives)
    missing = mem_objectives - current

    notes: list[NormalizationNote] = []
    if missing:
        for obj in sorted(missing):
            state.simulation_objectives.append(obj)

        # Set bilayer-specific hints only when a membrane is actually present in
        # the simulation box (membrane.enabled=true). A protein with
        # biological_context=membrane_associated but no bilayer in the box
        # should use standard protein parameters — not semiisotropic coupling.
        bilayer_in_box = (
            getattr(state.environment, "membrane", None) is not None
            and state.environment.membrane.enabled
        )
        if bilayer_in_box:
            state.workflow_hints.membrane_required      = True
            state.workflow_hints.semiisotropic_coupling = True
            state.workflow_hints.conservative_timestep  = True
            state.workflow_hints.membrane_equilibration = True
        else:
            # Membrane context detected but no bilayer simulated — mark awareness only.
            state.workflow_hints.membrane_required = True

        notes.append(NormalizationNote(
            kind="membrane_inferred",
            original="(auto-detected)",
            resolved=sorted(missing),
            message=(
                f"Membrane protein context detected — auto-injecting: "
                f"{sorted(missing)}\n"
                f"  Signals: {'; '.join(signals)}"
                + ("" if bilayer_in_box else
                   "\n  Note: membrane.enabled=false — bilayer-specific params NOT applied")
            ),
        ))

    return notes


# ── Main pipeline stage ────────────────────────────────────────────────────────

def run_semantic_normalization(state: "SystemState") -> "SystemState":
    """
    Stage 0.5: semantic normalization and intent inference.

    Runs after YAML load, before run_inference().

    Steps:
      1. Apply simulation_profile preset (if specified)
      2. Normalize simulation_objectives (aliases → canonical)
      3. Detect membrane protein signals → auto-inject objectives
      4. Record all normalization notes in global_reasoning
    """
    from core.models import Warning, Severity

    all_notes: list[NormalizationNote] = []

    # ── 1. Preset ─────────────────────────────────────────────────────────────
    profile = getattr(state, "simulation_profile", None)
    if profile:
        state, preset_note = apply_preset(state, profile)
        if preset_note:
            all_notes.append(preset_note)

    # ── 2. Normalize objectives ───────────────────────────────────────────────
    state.simulation_objectives, obj_notes = _normalize_objectives_list(
        state.simulation_objectives
    )
    all_notes.extend(obj_notes)

    # Emit warnings for unknowns with suggestions
    from core.semantic_objectives import suggest_objectives
    for note in obj_notes:
        if note.kind == "unknown":
            suggestions = suggest_objectives(note.original)
            hint = (
                f"Closest matches: {suggestions}" if suggestions
                else "Check core.semantic_objectives.CANONICAL_OBJECTIVES for valid names."
            )
            state.warnings.append(Warning(
                message=(
                    f"Unknown simulation objective: '{note.original}'\n"
                    f"  {hint}"
                ),
                target="simulation_objectives",
                severity=Severity.MEDIUM,
            ))

    # ── 3. Membrane auto-inference ────────────────────────────────────────────
    mem_notes = _inject_membrane_objectives(state)
    all_notes.extend(mem_notes)

    # ── 4. Record normalization trace in global_reasoning ────────────────────
    for note in all_notes:
        if note.kind != "unknown":
            state.global_reasoning.notes.append(
                f"[semantic] {note.message}"
            )

    return state
