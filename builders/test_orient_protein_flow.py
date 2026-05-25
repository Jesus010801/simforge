# builders/test_orient_protein_flow.py
"""
Integration test: compile membrane_orient_test.yaml → workspace → orient_protein.

Verifies the full chain:
    YAML with membrane.orientation
    → pipeline sets automation_level=AUTOMATED
    → MembraneOrientBuilder generates run.sh + orient_helper.py
    → metadata.json has automation_level=automated
    → dry-run executor does NOT skip orient_protein
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from core.parser import parse_yaml
from core.execution_models import AutomationLevel
from core.decision_engine import build_simulation_plan
from core.compiler import SimulationCompiler
from builders.workspace_builder import WorkspaceBuilder
from executors.shell_executor import ShellExecutor
from executors.execution_state import StepStatus


ORIENT_YAML = "configs/membrane_orient_test.yaml"
NO_ORIENT_YAML = "configs/membrane_test.yaml"


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def orient_workspace(tmp_path_factory):
    result = SimulationCompiler().compile(ORIENT_YAML)
    tmp = tmp_path_factory.mktemp("orient_ws")
    return WorkspaceBuilder().build(result, output_dir=str(tmp))


@pytest.fixture(scope="module")
def no_orient_workspace(tmp_path_factory):
    result = SimulationCompiler().compile(NO_ORIENT_YAML)
    tmp = tmp_path_factory.mktemp("no_orient_ws")
    return WorkspaceBuilder().build(result, output_dir=str(tmp))


@pytest.fixture(scope="module")
def orient_execution(orient_workspace):
    return ShellExecutor(orient_workspace, dry_run=True).run()


@pytest.fixture(scope="module")
def no_orient_execution(no_orient_workspace):
    return ShellExecutor(no_orient_workspace, dry_run=True).run()


# ── Pipeline: automation_level in SimulationStep ─────────────────────────────

def test_orient_protein_step_has_automated_level_when_orientation_present():
    result = SimulationCompiler().compile(ORIENT_YAML)
    orient_step = next(s for s in result.plan.steps if s.step_id == "orient_protein")
    assert orient_step.automation_level == AutomationLevel.AUTOMATED


def test_orient_protein_step_has_guided_level_when_orientation_absent():
    result = SimulationCompiler().compile(NO_ORIENT_YAML)
    orient_step = next(s for s in result.plan.steps if s.step_id == "orient_protein")
    assert orient_step.automation_level == AutomationLevel.GUIDED


def test_effective_automation_level_automated(orient_workspace):
    result = SimulationCompiler().compile(ORIENT_YAML)
    orient_step = next(s for s in result.plan.steps if s.step_id == "orient_protein")
    assert orient_step.effective_automation_level() == AutomationLevel.AUTOMATED


# ── Workspace: files generated ────────────────────────────────────────────────

def test_orient_protein_generates_run_sh(orient_workspace):
    step_dir = _orient_dir(orient_workspace)
    assert (step_dir / "run.sh").exists(), "run.sh should be generated when orientation is present"


def test_orient_protein_generates_orient_helper(orient_workspace):
    step_dir = _orient_dir(orient_workspace)
    assert (step_dir / "orient_helper.py").exists()


def test_orient_protein_no_run_sh_when_no_orientation(no_orient_workspace):
    step_dir = _orient_dir(no_orient_workspace)
    assert not (step_dir / "run.sh").exists(), "No run.sh without orientation — only README.md"
    assert (step_dir / "README.md").exists()


# ── Workspace: metadata.json automation_level ────────────────────────────────

def test_orient_protein_metadata_automation_level_automated(orient_workspace):
    meta = _orient_meta(orient_workspace)
    assert meta.get("automation_level") == "automated"


def test_orient_protein_metadata_automation_level_guided(no_orient_workspace):
    meta = _orient_meta(no_orient_workspace)
    assert meta.get("automation_level") == "guided"


def test_orient_protein_metadata_step_type_automatic(orient_workspace):
    meta = _orient_meta(orient_workspace)
    assert meta.get("step_type") == "automatic"


# ── Execution manifest ────────────────────────────────────────────────────────

def test_manifest_records_automation_level(orient_workspace):
    manifest = json.loads(
        (orient_workspace / "metadata" / "execution_manifest.json").read_text()
    )
    entry = next(e for e in manifest["steps"] if e["step_id"] == "orient_protein")
    assert entry.get("automation_level") == "automated"


# ── Dry-run execution ─────────────────────────────────────────────────────────

def test_orient_protein_is_not_skipped_with_orientation(orient_execution):
    record = _orient_record(orient_execution)
    assert record.status != StepStatus.SKIPPED, (
        f"orient_protein should NOT be skipped when automation_level=automated; "
        f"got status={record.status}"
    )
    assert record.status == StepStatus.DONE


def test_orient_protein_is_skipped_without_orientation(no_orient_execution):
    record = _orient_record(no_orient_execution)
    assert record.status == StepStatus.SKIPPED, (
        f"orient_protein should be SKIPPED when automation_level=guided (no orientation data)"
    )


def test_orient_helper_contains_ec_residues(orient_workspace):
    helper = ((_orient_dir(orient_workspace)) / "orient_helper.py").read_text()
    assert "EC_RESIDUES" in helper
    assert "1" in helper   # residue 1 from "1-30" range


def test_orient_helper_is_valid_python(orient_workspace):
    helper = ((_orient_dir(orient_workspace)) / "orient_helper.py").read_text()
    compile(helper, "orient_helper.py", "exec")  # raises SyntaxError if broken


def test_orient_helper_has_ec_target_side(orient_workspace):
    helper = ((_orient_dir(orient_workspace)) / "orient_helper.py").read_text()
    assert "EC_TARGET_SIDE" in helper
    assert '"+z"' in helper or "'+z'" in helper


def test_orient_metadata_has_extracellular_side(orient_workspace):
    meta = _orient_meta(orient_workspace)
    assert "extracellular_side" in meta
    assert meta["extracellular_side"] == "+z"


def test_orient_helper_inverted_produces_minus_z(tmp_path_factory):
    """When extracellular_side=-z, the helper should write -Z orientation logic."""
    import yaml as pyyaml
    from core.compiler import SimulationCompiler
    from builders.workspace_builder import WorkspaceBuilder

    content = """
project:
  name: inverted_orient_test

components:
  - id: protein_1
    role: protein
    file: protein.pdb

structural_annotation:
  membrane_topology:
    extracellular_regions: ["1-30"]
    intracellular_regions: ["70-100"]
    transmembrane_segments: ["31-69"]
  orientation:
    extracellular_side: "-z"
    intracellular_side: "+z"
    source: "user_annotation"
    confidence: 0.9

environment:
  membrane:
    enabled: true
    type: DPPC
  solvent:
    water_model: spce
  ions:
    concentration: 0.154
  temperature_K: 300.0
  duration_ns: 10.0

forcefields:
  protein: opls-aa

simulation_objectives:
  - membrane_protein_dynamics
"""
    yaml_dir = tmp_path_factory.mktemp("inv_yaml")
    p = yaml_dir / "inverted.yaml"
    p.write_text(content)
    (yaml_dir / "protein.pdb").write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00\n")

    result = SimulationCompiler().compile(str(p))
    tmp = tmp_path_factory.mktemp("inv_ws")
    ws = WorkspaceBuilder().build(result, output_dir=str(tmp))

    step_dir = next(d for d in (ws / "steps").iterdir() if "orient_protein" in d.name)
    helper = (step_dir / "orient_helper.py").read_text()

    assert '"-z"' in helper or "'-z'" in helper, (
        "orient_helper.py should embed EC_TARGET_SIDE = '-z' for inverted orientation"
    )

    meta = json.loads((step_dir / "metadata.json").read_text())
    assert meta.get("extracellular_side") == "-z"


# ── orient_protein does NOT block downstream when automated ──────────────────

def test_downstream_not_blocked_when_orient_automated(orient_execution):
    # match_box_to_bilayer depends on orient_protein; if orient_protein is DONE,
    # match_box_to_bilayer should not be BLOCKED (it may be SKIPPED if it's GUIDED,
    # but not BLOCKED due to orient_protein)
    status_map = {r.step_id: r.status for r in orient_execution.steps}
    orient_status = status_map.get("orient_protein")
    assert orient_status == StepStatus.DONE, "orient_protein must be DONE to unblock DAG"

    match_status = status_map.get("match_box_to_bilayer")
    if match_status is not None:
        assert match_status != StepStatus.BLOCKED, (
            "match_box_to_bilayer must not be BLOCKED when orient_protein completed"
        )


# ── Helpers ───────────────────────────────────────────────────────────────────

def _orient_dir(workspace: Path) -> Path:
    return next(
        d for d in (workspace / "steps").iterdir()
        if "orient_protein" in d.name
    )


def _orient_meta(workspace: Path) -> dict:
    return json.loads((_orient_dir(workspace) / "metadata.json").read_text())


def _orient_record(execution_state):
    return next(r for r in execution_state.steps if r.step_id == "orient_protein")
