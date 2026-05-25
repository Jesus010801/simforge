# executors/test_dag_blocking.py
"""
Tests for DAG blocking propagation in BaseExecutor._is_blocked().

Critical invariant: a BLOCKED step must itself block all downstream steps
that depend on it (transitive blocking). If minimization fails, equilibration
must be BLOCKED, and production must also be BLOCKED — not run.
"""
from __future__ import annotations

import pytest

from executors.base_executor import BaseExecutor
from executors.execution_state import (
    StepExecutionRecord,
    StepStatus,
    WorkspaceExecutionState,
)


# ── Minimal concrete subclass for unit testing ────────────────────────────────

class _NullExecutor(BaseExecutor):
    """BaseExecutor subclass that does nothing — only exposes _is_blocked()."""

    def _run_step(self, record):
        pass

    def _fake_state(self, steps: list[tuple[str, list[str], StepStatus]]) -> None:
        """Build a synthetic execution state from (step_id, depends_on, status) tuples."""
        records = [
            StepExecutionRecord(step_id=sid, step_dir="/fake/" + sid, depends_on=deps, status=status)
            for sid, deps, status in steps
        ]
        self.state = WorkspaceExecutionState(
            workspace_path="/fake",
            steps=records,
        )


@pytest.fixture
def executor(tmp_path):
    ex = _NullExecutor(workspace_path=tmp_path)
    return ex


# ── Direct blocking: FAILED dep → downstream BLOCKED ─────────────────────────

def test_failed_dep_blocks_direct_child(executor):
    executor._fake_state([
        ("minimization", [],               StepStatus.FAILED),
        ("equilibration", ["minimization"], StepStatus.PENDING),
    ])
    equil = executor.state.steps[1]
    assert executor._is_blocked(equil)


def test_done_dep_does_not_block(executor):
    executor._fake_state([
        ("minimization", [],               StepStatus.DONE),
        ("equilibration", ["minimization"], StepStatus.PENDING),
    ])
    equil = executor.state.steps[1]
    assert not executor._is_blocked(equil)


def test_skipped_dep_blocks_downstream(executor):
    """SKIPPED steps (manual/external not yet done) must block downstream steps.

    A manual step that was skipped means its outputs don't exist yet.
    Running downstream steps would fail with missing-input errors from GROMACS.
    """
    executor._fake_state([
        ("orient_protein", [],                    StepStatus.SKIPPED),
        ("embed_in_bilayer", ["orient_protein"],   StepStatus.PENDING),
    ])
    embed = executor.state.steps[1]
    assert executor._is_blocked(embed), (
        "embed_in_bilayer must be blocked when orient_protein (manual) is SKIPPED"
    )


# ── Transitive blocking: FAILED → BLOCKED → BLOCKED ──────────────────────────

def test_blocked_dep_blocks_grandchild(executor):
    """
    Core regression: equilibration=BLOCKED must block production.

    Before the fix, _is_blocked() only checked == FAILED,
    so a BLOCKED equilibration would NOT block production.
    """
    executor._fake_state([
        ("minimization", [],                StepStatus.FAILED),
        ("equilibration", ["minimization"],  StepStatus.BLOCKED),
        ("production",   ["equilibration"],  StepStatus.PENDING),
    ])
    production = executor.state.steps[2]
    assert executor._is_blocked(production), (
        "production must be blocked when equilibration is BLOCKED (transitive)"
    )


def test_transitive_chain_three_levels(executor):
    """minimization FAILED → equilibration BLOCKED → production BLOCKED → analysis BLOCKED."""
    executor._fake_state([
        ("minimization", [],                StepStatus.FAILED),
        ("equilibration", ["minimization"],  StepStatus.BLOCKED),
        ("production",   ["equilibration"],  StepStatus.BLOCKED),
        ("analysis",     ["production"],     StepStatus.PENDING),
    ])
    analysis = executor.state.steps[3]
    assert executor._is_blocked(analysis)


def test_partial_failure_blocks_when_one_dep_failed(executor):
    """production depends on [equilibration, rest2]. If rest2 FAILED → block production."""
    executor._fake_state([
        ("equilibration", [],                         StepStatus.DONE),
        ("rest2",         [],                         StepStatus.FAILED),
        ("production",    ["equilibration", "rest2"],  StepStatus.PENDING),
    ])
    production = executor.state.steps[2]
    assert executor._is_blocked(production)


def test_all_deps_done_does_not_block(executor):
    """production depends on [equilibration, rest2]. Both DONE → not blocked."""
    executor._fake_state([
        ("equilibration", [],                         StepStatus.DONE),
        ("rest2",         [],                         StepStatus.DONE),
        ("production",    ["equilibration", "rest2"],  StepStatus.PENDING),
    ])
    production = executor.state.steps[2]
    assert not executor._is_blocked(production)


# ── Fallback path: empty depends_on (filesystem-scan workspace) ───────────────

def test_fallback_failed_blocks_sequential(executor):
    """Without depends_on, any prior FAILED blocks all subsequent steps."""
    executor._fake_state([
        ("minimization", [], StepStatus.FAILED),
        ("equilibration", [], StepStatus.PENDING),
    ])
    equil = executor.state.steps[1]
    assert executor._is_blocked(equil)


def test_fallback_blocked_also_blocks_sequential(executor):
    """Without depends_on, a BLOCKED step should block subsequent steps too."""
    executor._fake_state([
        ("minimization", [], StepStatus.FAILED),
        ("equilibration", [], StepStatus.BLOCKED),
        ("production",    [], StepStatus.PENDING),
    ])
    production = executor.state.steps[2]
    assert executor._is_blocked(production)


def test_fallback_no_failure_does_not_block(executor):
    executor._fake_state([
        ("minimization", [], StepStatus.DONE),
        ("equilibration", [], StepStatus.PENDING),
    ])
    equil = executor.state.steps[1]
    assert not executor._is_blocked(equil)
