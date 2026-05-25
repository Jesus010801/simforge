# executors/test_executor.py
from __future__ import annotations

from executors.execution_state import StepStatus


# ── Dry-run completion ────────────────────────────────────────────────────────

def test_dry_run_not_complete_with_manual_steps(execution_state):
    # The hmg_competition workflow has a MANUAL review step that blocks
    # downstream automatic steps. is_complete is False until the user
    # completes the manual step and re-runs.
    assert not execution_state.is_complete


def test_dry_run_has_no_failed_steps(execution_state):
    assert execution_state.n_failed() == 0


def test_dry_run_done_count(execution_state):
    # Steps before the MANUAL review run: prepare × 3 + parametrize × 2 = 5
    assert execution_state.n_done() == 5


def test_dry_run_skipped_count(execution_state):
    # review_parametrization_substrate_1 is SKIPPED (manual, no script)
    skipped = sum(
        1 for s in execution_state.steps
        if s.status == StepStatus.SKIPPED
    )
    assert skipped == 1


# ── Step record correctness ───────────────────────────────────────────────────

def test_all_steps_have_step_id(execution_state):
    for record in execution_state.steps:
        assert record.step_id


def test_done_steps_have_elapsed_time(execution_state):
    for record in execution_state.steps:
        if record.status == StepStatus.DONE:
            assert record.elapsed_s is not None
            assert record.elapsed_s >= 0


def test_steps_have_depends_on_populated(execution_state):
    # manifest-driven execution must populate depends_on for non-root steps
    non_root = [
        r for r in execution_state.steps
        if r.step_id not in (
            "prepare_protein_1", "prepare_ligand_1", "prepare_substrate_1"
        )
    ]
    with_deps = [r for r in non_root if r.depends_on]
    assert len(with_deps) > 0


def test_analysis_steps_depend_on_production(execution_state):
    analysis_records = [
        r for r in execution_state.steps
        if r.step_id.startswith("analysis_")
    ]
    assert len(analysis_records) > 0
    for record in analysis_records:
        assert any("production" in d for d in record.depends_on)


# ── Blocked steps downstream of manual step ──────────────────────────────────

def test_assemble_and_downstream_are_blocked(execution_state):
    # assemble_system depends on the MANUAL review step → it must be BLOCKED.
    # All steps that follow it transitively are also BLOCKED.
    blocked = {r.step_id for r in execution_state.steps if r.status == StepStatus.BLOCKED}
    assert "assemble_system" in blocked, "assemble_system must be blocked (depends on manual review)"
    # production and analysis steps must also be blocked transitively
    assert "production_md" in blocked
    assert any(s.startswith("analysis_") for s in blocked)
