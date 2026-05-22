# core/test_decision_engine.py
from __future__ import annotations

from core.execution_models import StepStage


# ── WorkflowPolicy ────────────────────────────────────────────────────────────

def test_competitive_binding_sets_50ns(plan):
    assert plan.workflow_policy.production_time_ns == 50.0


def test_competitive_binding_enables_enhanced_sampling(plan):
    assert plan.workflow_policy.enhanced_sampling is True


def test_sampling_method_is_rest2(plan):
    assert plan.workflow_policy.sampling_method == "REST2"


def test_temperature_set_from_config(plan, state):
    assert plan.workflow_policy.temperature_K == state.environment.temperature_K


# ── Step count & structure ────────────────────────────────────────────────────

def test_plan_has_steps(plan):
    assert len(plan.steps) > 0


def test_plan_has_all_expected_stages(plan):
    stages = {s.stage for s in plan.steps}
    expected = {
        StepStage.PREPARATION,
        StepStage.PARAMETRIZATION,
        StepStage.ASSEMBLY,
        StepStage.MINIMIZATION,
        StepStage.EQUILIBRATION,
        StepStage.PRODUCTION,
        StepStage.ANALYSIS,
    }
    assert expected.issubset(stages)


# ── step.params propagation ───────────────────────────────────────────────────

def test_minimization_params_emtol_set(plan):
    step = next(s for s in plan.steps if s.stage == StepStage.MINIMIZATION)
    assert "emtol" in step.params
    assert step.params["emtol"] > 0


def test_minimization_flexible_ligand_raises_emtol(plan, state):
    # hmg_competition has a flexible ligand → emtol should be 100 (lenient)
    step = next(s for s in plan.steps if s.stage == StepStage.MINIMIZATION)
    assert step.params["emtol"] == 100.0


def test_equilibration_temperature_matches_policy(plan):
    step = next(s for s in plan.steps if s.stage == StepStage.EQUILIBRATION)
    assert step.params["temperature"] == plan.workflow_policy.temperature_K


def test_equilibration_has_nvt_nsteps(plan):
    step = next(s for s in plan.steps if s.stage == StepStage.EQUILIBRATION)
    assert step.params["nvt_nsteps"] > 0


def test_production_has_nsteps(plan):
    step = next(s for s in plan.steps if s.stage == StepStage.PRODUCTION)
    assert step.params["nsteps"] > 0


def test_production_nsteps_matches_50ns_policy(plan):
    # 50ns / 0.002ps timestep = 25,000,000 steps
    step = next(s for s in plan.steps if s.stage == StepStage.PRODUCTION)
    assert step.params["nsteps"] == 25_000_000


def test_analysis_steps_have_analysis_type(plan):
    analysis_steps = [s for s in plan.steps if s.stage == StepStage.ANALYSIS]
    assert len(analysis_steps) > 0
    for step in analysis_steps:
        assert "analysis_type" in step.params


def test_add_ions_concentration_is_physiological(plan):
    step = next((s for s in plan.steps if s.step_id == "add_ions"), None)
    assert step is not None
    assert step.params["concentration"] == 0.15


# ── DAG dependencies ──────────────────────────────────────────────────────────

def test_production_depends_on_equilibration(plan):
    step = next(s for s in plan.steps if s.stage == StepStage.PRODUCTION)
    assert any("equilibration" in d for d in step.depends_on)


def test_analysis_steps_depend_on_production(plan):
    analysis_steps = [s for s in plan.steps if s.stage == StepStage.ANALYSIS]
    for step in analysis_steps:
        assert any("production" in d for d in step.depends_on)
