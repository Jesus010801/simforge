# workflows/test_workflow_graph.py
from __future__ import annotations

import pytest

from workflows.workflow_graph import WorkflowGraph


@pytest.fixture(scope="module")
def graph(plan):
    return WorkflowGraph(plan)


# ── DAG integrity ─────────────────────────────────────────────────────────────

def test_dag_validates_without_error(graph):
    graph.validate()  # must not raise


def test_dag_has_root_steps(graph):
    assert len(graph.root_steps()) > 0


def test_dag_has_leaf_steps(graph):
    assert len(graph.leaf_steps()) > 0


def test_root_steps_have_no_predecessors(graph):
    for step in graph.root_steps():
        assert graph.predecessors(step.step_id) == []


def test_leaf_steps_have_no_successors(graph):
    for step in graph.leaf_steps():
        assert graph.successors(step.step_id) == []


# ── Critical path ─────────────────────────────────────────────────────────────

def test_critical_path_is_non_empty(graph):
    assert len(graph.critical_path()) > 0


def test_critical_path_starts_at_root(graph):
    cp = graph.critical_path()
    roots = {s.step_id for s in graph.root_steps()}
    assert cp[0].step_id in roots


def test_critical_path_ends_at_leaf(graph):
    cp = graph.critical_path()
    leaves = {s.step_id for s in graph.leaf_steps()}
    assert cp[-1].step_id in leaves


# ── Parallel waves ────────────────────────────────────────────────────────────

def test_parallel_waves_cover_all_steps(graph):
    waves = graph.to_parallel_waves()
    all_ids = {s.step_id for wave in waves for s in wave}
    plan_ids = set(graph.nodes.keys())
    assert all_ids == plan_ids


def test_analysis_steps_are_in_same_wave(graph):
    waves = graph.to_parallel_waves()
    analysis_ids = {
        s.step_id
        for wave in waves
        for s in wave
        if s.stage.value == "analysis"
    }
    # all analysis steps should appear in the same wave (they share the same dep)
    for wave in waves:
        wave_ids = {s.step_id for s in wave}
        if analysis_ids & wave_ids:
            assert analysis_ids.issubset(wave_ids)
            break


# ── Graph query API ───────────────────────────────────────────────────────────

def test_contains_known_step(graph):
    assert "assemble_system" in graph


def test_predecessors_of_assemble_are_preparation_steps(graph):
    preds = graph.predecessors("assemble_system")
    pred_stages = {p.stage.value for p in preds}
    assert pred_stages.issubset({"preparation", "parametrization", "validation"})


def test_successors_of_assemble_include_solvation(graph):
    succs = graph.successors("assemble_system")
    succ_ids = {s.step_id for s in succs}
    assert "solvate_system" in succ_ids


def test_ancestors_of_production_include_minimization(graph):
    ancs = graph.ancestors("production_md")
    anc_ids = {a.step_id for a in ancs}
    assert "energy_minimization" in anc_ids


def test_descendants_of_preparation_include_analysis(graph):
    descs = graph.descendants("prepare_protein_1")
    desc_stages = {d.stage.value for d in descs}
    assert "analysis" in desc_stages


# ── Stats ─────────────────────────────────────────────────────────────────────

def test_stats_returns_expected_keys(graph):
    s = graph.stats()
    for key in ("n_steps", "n_edges", "n_root_steps", "n_leaf_steps",
                "n_parallel_waves", "critical_path_len"):
        assert key in s


def test_stats_n_steps_matches_plan(graph, plan):
    assert graph.stats()["n_steps"] == len(plan.steps)


# ── Mermaid export ────────────────────────────────────────────────────────────

def test_mermaid_grouped_is_non_empty(graph):
    out = graph.render_mermaid(group_by_stage=True)
    assert "graph" in out
    assert len(out) > 50


def test_mermaid_flat_is_non_empty(graph):
    out = graph.render_mermaid(group_by_stage=False)
    assert "graph" in out
