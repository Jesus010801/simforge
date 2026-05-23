"""
Tests for shared topology lifecycle across the GROMACS workflow.

Critical invariants:
  1. solvate_system produces a LOCAL topol.top (copy-on-modify pattern)
  2. add_ions produces a LOCAL topol.top (copy-on-modify pattern)
  3. minimization/equilibration/production resolve topology from add_ions,
     not from assemble_system (which lacks ion/water counts)
  4. Membrane workflows follow the same chain via solvate_membrane and clean_water

The error that motivated these tests:
  '../03_solvate_system/topol.top' does not exist
  → solvate_system was modifying assemble_system/topol.top in-place
    without keeping a local copy, so add_ions had no topology to read.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from core.compiler import SimulationCompiler
from builders.workspace_builder import WorkspaceBuilder


# ── Fixtures ──────────────────────────────────────────────────────────────────

PROTEIN_YAML   = "configs/lysozyme_test.yaml"
MEMBRANE_YAML  = "configs/membrane_test.yaml"


@pytest.fixture(scope="module")
def protein_workspace(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("protein")
    result = SimulationCompiler().compile(PROTEIN_YAML)
    return WorkspaceBuilder().build(result, output_dir=str(tmp))


@pytest.fixture(scope="module")
def protein_steps(protein_workspace):
    """Return {step_id: metadata_dict} for the protein workspace."""
    steps_dir = protein_workspace / "steps"
    return {
        p.name.split("_", 1)[1]: json.loads((p / "metadata.json").read_text())
        for p in sorted(steps_dir.iterdir())
        if p.is_dir() and (p / "metadata.json").exists()
    }


def _step(protein_steps, step_id):
    """Find a step by partial name match (ignores numeric prefix)."""
    for name, meta in protein_steps.items():
        if meta.get("step_id") == step_id:
            return meta
    return None


# ── solvate_system ─────────────────────────────────────────────────────────────

class TestSolvateTopology:
    def test_solvate_declares_topol_in_expected_outputs(self, protein_steps):
        meta = _step(protein_steps, "solvate_system")
        assert meta is not None, "solvate_system step not found in workspace"
        assert "topol.top" in meta.get("expected_outputs", []), (
            "solvate_system must declare topol.top as an expected output "
            "(it copies topology locally before gmx solvate modifies it)"
        )

    def test_solvate_run_script_copies_topology(self, protein_workspace):
        steps_dir = protein_workspace / "steps"
        solvate_dir = next(
            (p for p in steps_dir.iterdir() if "solvate_system" in p.name),
            None,
        )
        assert solvate_dir is not None
        script = (solvate_dir / "run.sh").read_text()
        assert "cp " in script and "topol.top" in script, (
            "solvate_system/run.sh must copy topol.top before calling gmx solvate"
        )

    def test_solvate_uses_local_topology_in_gromacs_call(self, protein_workspace):
        steps_dir = protein_workspace / "steps"
        solvate_dir = next(
            (p for p in steps_dir.iterdir() if "solvate_system" in p.name),
            None,
        )
        script = (solvate_dir / "run.sh").read_text()
        # gmx solvate must use '-p topol.top' (local), not '-p $ASSEMBLE_DIR/topol.top'
        assert "-p topol.top" in script, (
            "gmx solvate must use the local topol.top copy, "
            "not '$ASSEMBLE_DIR/topol.top'"
        )


# ── add_ions ──────────────────────────────────────────────────────────────────

class TestAddIonsTopology:
    def test_add_ions_declares_topol_in_expected_outputs(self, protein_steps):
        meta = _step(protein_steps, "add_ions")
        assert meta is not None, "add_ions step not found in workspace"
        assert "topol.top" in meta.get("expected_outputs", []), (
            "add_ions must declare topol.top as an expected output "
            "(gmx genion modifies it in-place)"
        )

    def test_add_ions_copies_topology_before_genion(self, protein_workspace):
        steps_dir = protein_workspace / "steps"
        ions_dir = next(
            (p for p in steps_dir.iterdir() if "add_ions" in p.name),
            None,
        )
        assert ions_dir is not None
        script = (ions_dir / "run.sh").read_text()
        assert "cp " in script and "topol.top" in script, (
            "add_ions/run.sh must copy topol.top locally before genion modifies it"
        )

    def test_add_ions_reads_from_solvate_not_assemble(self, protein_workspace):
        """Topology source must be solvate_system (post-SOL), not assemble_system."""
        steps_dir = protein_workspace / "steps"
        ions_dir = next(
            (p for p in steps_dir.iterdir() if "add_ions" in p.name),
            None,
        )
        script = (ions_dir / "run.sh").read_text()
        # The TOPOL_SRC variable must point toward solvate, not assemble
        assert "solvate" in script.lower(), (
            "add_ions must read topology from solvate_system "
            "(which has SOL entries), not directly from assemble_system"
        )

    def test_add_ions_genion_uses_local_topology(self, protein_workspace):
        steps_dir = protein_workspace / "steps"
        ions_dir = next(
            (p for p in steps_dir.iterdir() if "add_ions" in p.name),
            None,
        )
        script = (ions_dir / "run.sh").read_text()
        # Both grompp and genion must use the local copy
        assert script.count("-p topol.top") >= 2, (
            "Both gmx grompp and gmx genion in add_ions must use '-p topol.top' "
            "(the local copy), not a path into another step's directory"
        )


# ── minimization topology chain ───────────────────────────────────────────────

class TestMinimizationTopologyChain:
    def test_minimization_required_inputs_includes_topol(self, protein_steps):
        meta = _step(protein_steps, "energy_minimization")
        assert meta is not None, "energy_minimization step not found"
        required = meta.get("required_inputs", [])
        has_topol = any("topol.top" in r for r in required)
        assert has_topol, (
            f"energy_minimization required_inputs must reference topol.top; got {required}"
        )

    def test_minimization_topol_path_points_to_add_ions(self, protein_steps):
        meta = _step(protein_steps, "energy_minimization")
        required = meta.get("required_inputs", [])
        topol_path = next((r for r in required if "topol.top" in r), None)
        assert topol_path is not None
        assert "add_ions" in topol_path, (
            f"energy_minimization must read topol.top from add_ions "
            f"(topology with both SOL and ions), got: {topol_path}"
        )

    def test_minimization_run_script_topol_dir(self, protein_workspace):
        steps_dir = protein_workspace / "steps"
        em_dir = next(
            (p for p in steps_dir.iterdir() if "energy_minimization" in p.name),
            None,
        )
        assert em_dir is not None
        script = (em_dir / "run.sh").read_text()
        assert "add_ions" in script, (
            "energy_minimization/run.sh TOPOL_DIR must point to add_ions"
        )


# ── equilibration topology chain ──────────────────────────────────────────────

class TestEquilibrationTopologyChain:
    def test_equilibration_run_script_topol_dir(self, protein_workspace):
        steps_dir = protein_workspace / "steps"
        eq_dir = next(
            (p for p in steps_dir.iterdir() if "equilibration" in p.name),
            None,
        )
        assert eq_dir is not None
        # Check NVT script (run_nvt.sh)
        nvt = (eq_dir / "run_nvt.sh").read_text()
        assert "add_ions" in nvt, (
            "equilibration/run_nvt.sh TOPOL_DIR must point to add_ions"
        )

    def test_npt_run_script_topol_dir(self, protein_workspace):
        steps_dir = protein_workspace / "steps"
        eq_dir = next(
            (p for p in steps_dir.iterdir() if "equilibration" in p.name),
            None,
        )
        npt = (eq_dir / "run_npt.sh").read_text()
        assert "add_ions" in npt, (
            "equilibration/run_npt.sh TOPOL_DIR must point to add_ions"
        )


# ── production topology chain ─────────────────────────────────────────────────

class TestProductionTopologyChain:
    def test_production_run_script_topol_dir(self, protein_workspace):
        steps_dir = protein_workspace / "steps"
        prod_dir = next(
            (p for p in steps_dir.iterdir() if "production_md" in p.name),
            None,
        )
        assert prod_dir is not None
        # production builder writes run_md.sh (not run.sh)
        script_file = prod_dir / "run_md.sh"
        if not script_file.exists():
            script_file = prod_dir / "run.sh"
        script = script_file.read_text()
        assert "add_ions" in script, (
            f"{script_file.name} TOPOL_DIR must point to add_ions"
        )


# ── topology chain integrity (all steps read from the same canonical source) ──

class TestTopologyChainIntegrity:
    def test_assemble_does_not_appear_as_topol_source_in_late_steps(
        self, protein_workspace
    ):
        """
        After solvate+ions, no late step should pull topol.top directly from
        assemble_system. That topology lacks SOL and ion entries.
        """
        steps_dir = protein_workspace / "steps"
        late_steps = ["energy_minimization", "equilibration", "production_md"]

        for step_name in late_steps:
            step_dir = next(
                (p for p in steps_dir.iterdir() if step_name in p.name),
                None,
            )
            if step_dir is None:
                continue

            for script_file in step_dir.glob("*.sh"):
                text = script_file.read_text()
                # Detect lines setting TOPOL_DIR to assemble_system
                for line in text.splitlines():
                    if "TOPOL_DIR" in line and "assemble_system" in line:
                        pytest.fail(
                            f"{step_dir.name}/{script_file.name} sets TOPOL_DIR "
                            f"to assemble_system — should use add_ions: {line.strip()}"
                        )
