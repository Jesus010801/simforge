# core/test_compiler.py
from __future__ import annotations


# ── CompilationResult structure ───────────────────────────────────────────────

def test_compilation_result_has_execution_order(compilation_result):
    assert len(compilation_result.execution_order) > 0


def test_execution_order_starts_with_preparation(compilation_result):
    first = compilation_result.execution_order[0]
    assert first.stage.value == "preparation"


def test_execution_order_ends_with_analysis(compilation_result):
    last = compilation_result.execution_order[-1]
    assert last.stage.value == "analysis"


def test_mermaid_graph_is_non_empty(compilation_result):
    assert compilation_result.mermaid_graph
    assert "graph" in compilation_result.mermaid_graph


def test_user_view_is_non_empty(compilation_result):
    assert len(compilation_result.user_view) > 0


def test_compilation_result_has_plan(compilation_result):
    assert compilation_result.plan is not None


def test_compilation_result_has_state(compilation_result):
    assert compilation_result.state is not None


def test_execution_order_has_no_duplicate_step_ids(compilation_result):
    ids = [s.step_id for s in compilation_result.execution_order]
    assert len(ids) == len(set(ids))


def test_execution_order_count_matches_plan(compilation_result):
    plan_count = len(compilation_result.plan.steps)
    order_count = len(compilation_result.execution_order)
    assert plan_count == order_count
