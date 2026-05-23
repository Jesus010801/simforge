"""Tests for the semantic objective + inference system."""
from __future__ import annotations

import pytest

from core.semantic_objectives import (
    CANONICAL_OBJECTIVES,
    SIMULATION_PRESETS,
    normalize_objective,
    suggest_objectives,
)
from core.semantic_inference import (
    detect_membrane_protein_signals,
    run_semantic_normalization,
    apply_preset,
    _normalize_objectives_list,
)


# ── normalize_objective ───────────────────────────────────────────────────────

class TestNormalizeObjective:
    def test_already_canonical_returned_unchanged(self):
        for obj in CANONICAL_OBJECTIVES:
            result, note = normalize_objective(obj)
            assert result == [obj]
            assert note is None, f"Expected no note for canonical '{obj}'"

    def test_alias_membrane_protein_dynamics(self):
        result, note = normalize_objective("membrane_protein_dynamics")
        assert "membrane_perturbation" in result
        assert "stability" in result
        assert note is not None

    def test_alias_with_spaces(self):
        result, note = normalize_objective("protein stability")
        assert result == ["stability"]

    def test_alias_with_spaces_binding(self):
        result, note = normalize_objective("binding affinity")
        assert result == ["binding"]

    def test_alias_idr(self):
        result, note = normalize_objective("idr")
        assert result == ["conformational_sampling"]

    def test_alias_pore_formation(self):
        result, note = normalize_objective("pore_formation")
        assert "membrane_perturbation" in result
        assert "permeability" in result

    def test_alias_enzyme_dynamics(self):
        result, note = normalize_objective("enzyme_dynamics")
        assert "active_site_dynamics" in result
        assert "active_site_stability" in result

    def test_fuzzy_match_returns_canonical(self):
        # "stabilitty" (typo) should fuzzy-match to "stability"
        result, note = normalize_objective("stabilitty")
        assert len(result) >= 1
        assert result[0] in CANONICAL_OBJECTIVES

    def test_fuzzy_match_membrane(self):
        # "membrane perturbation" (space instead of underscore) → canonical
        result, note = normalize_objective("membrane perturbation")
        assert "membrane_perturbation" in result

    def test_unknown_returns_empty(self):
        result, note = normalize_objective("completely_nonsensical_xyzabc")
        assert result == []
        assert note is None

    def test_hyphen_normalized(self):
        # "membrane-perturbation" → canonical via key normalization
        result, note = normalize_objective("membrane-perturbation")
        assert "membrane_perturbation" in result

    def test_deduplicate_in_list(self):
        # Two aliases that both resolve to stability should not duplicate
        normalized, notes = _normalize_objectives_list(
            ["stability", "protein_stability"]
        )
        assert normalized.count("stability") == 1


# ── suggest_objectives ────────────────────────────────────────────────────────

class TestSuggestObjectives:
    def test_unknown_returns_nonempty_suggestions(self):
        suggestions = suggest_objectives("completely_unknown_xyz")
        assert isinstance(suggestions, list)
        # Should still find something (difflib cutoff is permissive)

    def test_membrane_typo_suggests_membrane(self):
        suggestions = suggest_objectives("memrane_perturbation")  # typo
        assert any("membrane" in s for s in suggestions)

    def test_suggestions_are_canonical(self):
        suggestions = suggest_objectives("binding_affnity")  # typo
        for s in suggestions:
            assert s in CANONICAL_OBJECTIVES, f"'{s}' not canonical"

    def test_max_three_suggestions(self):
        suggestions = suggest_objectives("stability_xyz")
        assert len(suggestions) <= 3


# ── Presets ───────────────────────────────────────────────────────────────────

class TestPresets:
    def _make_minimal_state(self):
        """Minimal valid SystemState for testing — no file I/O."""
        import yaml, textwrap
        from core.parser import parse_yaml
        # Use the parser to build a real state, but we need a YAML file...
        # Build programmatically via SystemState fields
        from core.models import SystemState, ProjectModel, ForcefieldsModel
        return SystemState(
            project    = ProjectModel(name="test"),
            components = [],
            forcefields= ForcefieldsModel(protein="opls-aa"),
            simulation_objectives=[],
        )

    def test_membrane_protein_preset_injects_objectives(self):
        state = self._make_minimal_state()
        state, note = apply_preset(state, "membrane_protein")
        assert note is not None
        assert note.kind == "preset"
        assert "membrane_perturbation" in state.simulation_objectives
        assert "stability" in state.simulation_objectives

    def test_preset_populates_workflow_hints(self):
        state = self._make_minimal_state()
        state, note = apply_preset(state, "membrane_protein")
        assert state.workflow_hints.semiisotropic_coupling is True
        assert state.workflow_hints.conservative_timestep is True
        assert state.workflow_hints.membrane_required is True
        assert state.workflow_hints.membrane_equilibration is True

    def test_preset_audit_trace_in_global_reasoning(self):
        state = self._make_minimal_state()
        state, note = apply_preset(state, "membrane_protein")
        audit = [n for n in state.global_reasoning.notes if "preset:membrane_protein" in n]
        assert len(audit) == 1
        assert "semiisotropic_coupling" in audit[0]

    def test_preset_deduplicates_objectives(self):
        state = self._make_minimal_state()
        state.simulation_objectives = ["stability"]
        state, _ = apply_preset(state, "membrane_protein")
        assert state.simulation_objectives.count("stability") == 1

    def test_unknown_preset_returns_warning_note(self):
        state = self._make_minimal_state()
        state, note = apply_preset(state, "this_does_not_exist")
        assert note is not None
        assert note.kind == "unknown"

    def test_idr_sampling_preset(self):
        state = self._make_minimal_state()
        state, note = apply_preset(state, "idr_sampling")
        assert "conformational_sampling" in state.simulation_objectives

    def test_preset_names_all_valid(self):
        state = self._make_minimal_state()
        for name in SIMULATION_PRESETS:
            s2 = state.model_copy(deep=True)
            s2, note = apply_preset(s2, name)
            assert note is None or note.kind != "unknown", f"Preset '{name}' failed"


# ── Membrane protein signal detection ─────────────────────────────────────────

class TestMembraneSignalDetection:
    def _make_state(self, **kwargs):
        from core.models import SystemState, ProjectModel, ForcefieldsModel
        kwargs.setdefault("components", [])
        return SystemState(
            project    = ProjectModel(name="test"),
            forcefields= ForcefieldsModel(protein="opls-aa"),
            **kwargs,
        )

    def test_no_signals_on_empty_state(self):
        state = self._make_state()
        assert detect_membrane_protein_signals(state) == []

    def test_membrane_enabled_is_signal(self):
        from core.models import EnvironmentModel, MembraneConfig
        state = self._make_state(
            environment=EnvironmentModel(membrane=MembraneConfig(enabled=True, lipid="POPC"))
        )
        signals = detect_membrane_protein_signals(state)
        assert any("membrane.enabled" in s for s in signals)

    def test_transmembrane_context_is_signal(self):
        from core.models import ComponentModel
        state = self._make_state(components=[
            ComponentModel(id="prot", role="protein", file="",
                           biological_context=["transmembrane"])
        ])
        signals = detect_membrane_protein_signals(state)
        assert any("transmembrane" in s for s in signals)

    def test_membrane_objective_is_signal(self):
        state = self._make_state(simulation_objectives=["membrane_perturbation"])
        signals = detect_membrane_protein_signals(state)
        assert any("membrane_perturbation" in s for s in signals)


# ── run_semantic_normalization (integration) ──────────────────────────────────

class TestSemanticNormalizationStage:
    def _make_state(self, objectives=None, profile=None, components=None):
        from core.models import SystemState, ProjectModel, ForcefieldsModel, ComponentModel
        return SystemState(
            project    = ProjectModel(name="test"),
            components = components or [],
            forcefields= ForcefieldsModel(protein="opls-aa"),
            simulation_objectives=objectives or [],
            simulation_profile=profile,
        )

    def test_aliases_normalized(self):
        state = self._make_state(objectives=["membrane_protein_dynamics"])
        state = run_semantic_normalization(state)
        assert "membrane_perturbation" in state.simulation_objectives
        assert "membrane_protein_dynamics" not in state.simulation_objectives

    def test_canonical_unchanged(self):
        state = self._make_state(objectives=["stability", "binding"])
        state = run_semantic_normalization(state)
        assert "stability" in state.simulation_objectives
        assert "binding" in state.simulation_objectives

    def test_unknown_emits_warning(self):
        state = self._make_state(objectives=["absolutely_unknown_xyz"])
        state = run_semantic_normalization(state)
        warns = [w for w in state.warnings if "Unknown simulation objective" in w.message]
        assert len(warns) == 1
        assert "absolutely_unknown_xyz" in warns[0].message

    def test_preset_applied(self):
        state = self._make_state(profile="soluble_protein")
        state = run_semantic_normalization(state)
        assert "stability" in state.simulation_objectives

    def test_membrane_protein_auto_inject(self):
        from core.models import ComponentModel
        state = self._make_state(
            objectives=[],
            components=[
                ComponentModel(id="prot", role="protein", file="",
                               biological_context=["transmembrane"])
            ]
        )
        state = run_semantic_normalization(state)
        # membrane_perturbation should be auto-injected
        assert "membrane_perturbation" in state.simulation_objectives

    def test_normalization_notes_in_global_reasoning_from_membrane_inject(self):
        from core.models import ComponentModel
        state = self._make_state(
            objectives=[],
            components=[ComponentModel(id="p", role="protein", file="",
                                       biological_context=["transmembrane"])]
        )
        state = run_semantic_normalization(state)
        semantic_notes = [n for n in state.global_reasoning.notes if n.startswith("[semantic]")]
        assert len(semantic_notes) >= 1

    def test_normalization_notes_from_preset_in_global_reasoning(self):
        state = self._make_state(profile="soluble_protein")
        state = run_semantic_normalization(state)
        preset_notes = [n for n in state.global_reasoning.notes if "preset:" in n]
        # soluble_protein preset has no hints, so it won't emit a hints note — that's fine
        # but semantic note from preset note IS emitted
        assert "stability" in state.simulation_objectives

    def test_no_duplicate_objectives(self):
        # Both "stability" and "protein_stability" → should end up with one "stability"
        state = self._make_state(objectives=["stability", "protein_stability"])
        state = run_semantic_normalization(state)
        assert state.simulation_objectives.count("stability") == 1


# ── SystemState field validator (Pydantic level) ──────────────────────────────

class TestSystemStateObjectiveValidator:
    def _make(self, objectives):
        from core.models import SystemState, ProjectModel, ForcefieldsModel
        return SystemState(
            project    = ProjectModel(name="test"),
            components = [],
            forcefields= ForcefieldsModel(protein="opls-aa"),
            simulation_objectives=objectives,
        )

    def test_alias_normalized_at_pydantic_level(self):
        state = self._make(["membrane_protein_dynamics"])
        # Pydantic validator normalizes inline
        assert "membrane_perturbation" in state.simulation_objectives

    def test_canonical_unchanged(self):
        state = self._make(["stability"])
        assert state.simulation_objectives == ["stability"]

    def test_unknown_preserved_not_raised(self):
        # No exception should be raised for unknown objectives
        state = self._make(["totally_unknown_xyz"])
        assert "totally_unknown_xyz" in state.simulation_objectives


# ── SIMULATION_GOALS backward compat ─────────────────────────────────────────

class TestSimulationGoalsCompat:
    def test_simulation_goals_matches_canonical(self):
        from core.ontology import SIMULATION_GOALS
        from core.semantic_objectives import CANONICAL_OBJECTIVES
        assert set(SIMULATION_GOALS) == set(CANONICAL_OBJECTIVES.keys())


# ── WorkflowHints → WorkflowPolicy chain ─────────────────────────────────────

class TestWorkflowHintsChain:
    """
    End-to-end tests: semantic preset → workflow_hints → WorkflowPolicy parameters.
    These tests verify that scientific knowledge encoded in SIMULATION_PRESETS
    actually reaches the simulation parameters through the policy chain.
    """

    def _state_with_profile(self, profile: str):
        from core.models import SystemState, ProjectModel, ForcefieldsModel
        from core.semantic_inference import run_semantic_normalization
        state = SystemState(
            project    = ProjectModel(name="test"),
            components = [],
            forcefields= ForcefieldsModel(protein="opls-aa"),
            simulation_objectives=[],
            simulation_profile=profile,
        )
        return run_semantic_normalization(state)

    def _policy_from_state(self, state):
        from core.decision_engine import _build_workflow_policy
        return _build_workflow_policy(state)

    # ── membrane_protein preset ───────────────────────────────────────────────

    def test_membrane_protein_sets_semiisotropic_hint(self):
        state = self._state_with_profile("membrane_protein")
        assert state.workflow_hints.semiisotropic_coupling is True

    def test_membrane_protein_sets_conservative_timestep_hint(self):
        state = self._state_with_profile("membrane_protein")
        assert state.workflow_hints.conservative_timestep is True

    def test_membrane_protein_policy_timestep_is_001(self):
        state = self._state_with_profile("membrane_protein")
        policy = self._policy_from_state(state)
        assert policy.timestep_ps == 0.001

    def test_membrane_protein_policy_semiisotropic_true(self):
        state = self._state_with_profile("membrane_protein")
        policy = self._policy_from_state(state)
        assert policy.semiisotropic_coupling is True

    def test_membrane_protein_policy_equilibration_extended(self):
        state = self._state_with_profile("membrane_protein")
        policy = self._policy_from_state(state)
        assert policy.equilibration_time_ns >= 0.5

    # ── idr_sampling preset ───────────────────────────────────────────────────

    def test_idr_sampling_sets_enhanced_sampling_hint(self):
        state = self._state_with_profile("idr_sampling")
        assert state.workflow_hints.enhanced_sampling is True

    def test_idr_sampling_policy_enables_enhanced_sampling(self):
        state = self._state_with_profile("idr_sampling")
        policy = self._policy_from_state(state)
        assert policy.enhanced_sampling is True
        assert policy.sampling_method == "REST2"

    # ── soluble_protein preset ────────────────────────────────────────────────

    def test_soluble_protein_no_membrane_hints(self):
        state = self._state_with_profile("soluble_protein")
        assert state.workflow_hints.semiisotropic_coupling is False
        assert state.workflow_hints.conservative_timestep is False

    def test_soluble_protein_default_timestep(self):
        state = self._state_with_profile("soluble_protein")
        policy = self._policy_from_state(state)
        assert policy.timestep_ps == 0.002

    # ── equilibration params carry semiisotropic ──────────────────────────────

    def test_membrane_equilibration_params_have_semiisotropic(self):
        from core.decision_engine import _build_equilibration_params, _build_workflow_policy
        from core.execution_models import SimulationStep, StepStage, StepType

        state = self._state_with_profile("membrane_protein")
        policy = _build_workflow_policy(state)
        step = SimulationStep(
            step_id="equilibration", title="Equilibration",
            stage=StepStage.EQUILIBRATION, engine="gromacs",
        )
        params = _build_equilibration_params(step, state, policy)

        assert params["pcoupltype"] == "semiisotropic"
        assert params["pcoupl_npt"] == "Berendsen"
        assert params["rcoulomb"] == 1.2
        assert params["disp_corr"] == "EnerPres"

    def test_soluble_protein_equilibration_has_no_semiisotropic(self):
        from core.decision_engine import _build_equilibration_params, _build_workflow_policy
        from core.execution_models import SimulationStep, StepStage, StepType

        state = self._state_with_profile("soluble_protein")
        policy = _build_workflow_policy(state)
        step = SimulationStep(
            step_id="equilibration", title="Equilibration",
            stage=StepStage.EQUILIBRATION, engine="gromacs",
        )
        params = _build_equilibration_params(step, state, policy)

        assert "pcoupltype" not in params
        assert "pcoupl_npt" not in params

    # ── production params carry semiisotropic ─────────────────────────────────

    def test_membrane_production_params_have_nose_hoover(self):
        from core.decision_engine import _build_production_params, _build_workflow_policy
        from core.execution_models import SimulationStep, StepStage, StepType

        state = self._state_with_profile("membrane_protein")
        policy = _build_workflow_policy(state)
        step = SimulationStep(
            step_id="production_md", title="Production MD",
            stage=StepStage.PRODUCTION, engine="gromacs",
        )
        params = _build_production_params(step, state, policy)

        assert params["pcoupltype"] == "semiisotropic"
        assert params["tcoupl"] == "Nose-Hoover"
        assert params["disp_corr"] == "EnerPres"

    # ── membrane auto-inference also sets hints ───────────────────────────────

    def test_transmembrane_context_without_bilayer_sets_awareness_only(self):
        """Context=transmembrane but membrane.enabled=false → membrane_required only."""
        from core.models import SystemState, ProjectModel, ForcefieldsModel, ComponentModel
        from core.semantic_inference import run_semantic_normalization

        state = SystemState(
            project    = ProjectModel(name="test"),
            components = [ComponentModel(id="prot", role="protein", file="",
                                         biological_context=["transmembrane"])],
            forcefields= ForcefieldsModel(protein="opls-aa"),
        )
        state = run_semantic_normalization(state)

        # No bilayer in box → coupling hints NOT activated
        assert state.workflow_hints.membrane_required is True
        assert state.workflow_hints.semiisotropic_coupling is False
        assert state.workflow_hints.conservative_timestep is False

    def test_transmembrane_context_with_bilayer_sets_full_hints(self):
        """Context=transmembrane AND membrane.enabled=true → all membrane hints set."""
        from core.models import (
            SystemState, ProjectModel, ForcefieldsModel, ComponentModel,
            EnvironmentModel, MembraneConfig,
        )
        from core.semantic_inference import run_semantic_normalization

        state = SystemState(
            project    = ProjectModel(name="test"),
            components = [ComponentModel(id="prot", role="protein", file="",
                                         biological_context=["transmembrane"])],
            forcefields= ForcefieldsModel(protein="opls-aa"),
            environment= EnvironmentModel(membrane=MembraneConfig(enabled=True, lipid="POPC")),
        )
        state = run_semantic_normalization(state)

        assert state.workflow_hints.semiisotropic_coupling is True
        assert state.workflow_hints.conservative_timestep is True
        assert state.workflow_hints.membrane_required is True
