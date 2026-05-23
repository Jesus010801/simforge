"""
Tests for stage-aware mdrun flag generation.

GROMACS >=2024 compatibility matrix:
  minimization (steep/cg/l-bfgs — non-dynamical):
    - -nb gpu    → allowed
    - -pme gpu   → NOT ALLOWED  (causes: PME GPU does not support Non-dynamical integrator)
    - -bonded gpu→ NOT ALLOWED
    - -pmefft gpu→ NOT ALLOWED

  equilibration / production (md — dynamical):
    - all GPU flags allowed

These tests ensure the generated mdrun commands respect these constraints,
regardless of the hardware setting (gpu / cpu / auto).
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from builders.step_builders._utils import mdrun_block
from core.compiler import SimulationCompiler
from builders.workspace_builder import WorkspaceBuilder


PROTEIN_YAML = "configs/lysozyme_test.yaml"


# ── Unit tests for mdrun_block ────────────────────────────────────────────────

class TestMdrunBlockMinimization:
    """Stage='minimization' must never emit -pme gpu, -bonded gpu, or -pmefft gpu."""

    @pytest.mark.parametrize("hardware", ["gpu", "auto"])
    def test_no_pme_gpu(self, hardware):
        block = mdrun_block("em", hardware, stage="minimization")
        assert "-pme gpu" not in block, (
            f"hardware={hardware}: minimization must not use -pme gpu\n{block}"
        )

    @pytest.mark.parametrize("hardware", ["gpu", "auto"])
    def test_no_bonded_gpu(self, hardware):
        block = mdrun_block("em", hardware, stage="minimization")
        assert "-bonded gpu" not in block, (
            f"hardware={hardware}: minimization must not use -bonded gpu\n{block}"
        )

    @pytest.mark.parametrize("hardware", ["gpu", "auto"])
    def test_no_pmefft_gpu(self, hardware):
        block = mdrun_block("em", hardware, stage="minimization")
        assert "-pmefft gpu" not in block, (
            f"hardware={hardware}: minimization must not use -pmefft gpu\n{block}"
        )

    @pytest.mark.parametrize("hardware", ["gpu", "auto"])
    def test_pme_cpu_explicit(self, hardware):
        """PME must run on CPU during minimization."""
        block = mdrun_block("em", hardware, stage="minimization")
        assert "-pme cpu" in block, (
            f"hardware={hardware}: minimization must set -pme cpu explicitly\n{block}"
        )

    @pytest.mark.parametrize("hardware", ["gpu", "auto"])
    def test_nb_gpu_allowed(self, hardware):
        """Neighbor-list on GPU is valid even for non-dynamical integrators."""
        block = mdrun_block("em", hardware, stage="minimization")
        assert "-nb gpu" in block, (
            f"hardware={hardware}: minimization should offload -nb to GPU\n{block}"
        )

    def test_cpu_hardware_no_gpu_flags(self):
        block = mdrun_block("em", "cpu", stage="minimization")
        assert "-nb gpu" not in block
        assert "-pme gpu" not in block
        assert "$(nproc)" in block


class TestMdrunBlockMD:
    """Stage='md' (equilibration / production) gets full GPU offloading."""

    @pytest.mark.parametrize("hardware", ["gpu", "auto"])
    def test_pme_gpu_present(self, hardware):
        block = mdrun_block("nvt", hardware, stage="md")
        assert "-pme gpu" in block, (
            f"hardware={hardware}: md stage should use -pme gpu\n{block}"
        )

    @pytest.mark.parametrize("hardware", ["gpu", "auto"])
    def test_bonded_gpu_present(self, hardware):
        block = mdrun_block("nvt", hardware, stage="md")
        assert "-bonded gpu" in block

    @pytest.mark.parametrize("hardware", ["gpu", "auto"])
    def test_nb_gpu_present(self, hardware):
        block = mdrun_block("nvt", hardware, stage="md")
        assert "-nb gpu" in block

    def test_cpu_hardware_no_gpu_flags(self):
        block = mdrun_block("nvt", "cpu", stage="md")
        assert "-pme gpu" not in block
        assert "-nb gpu" not in block


class TestMdrunBlockAutoDetect:
    """auto hardware must produce a bash if/else that branches to GPU or CPU."""

    def test_minimization_auto_has_nvidia_check(self):
        block = mdrun_block("em", "auto", stage="minimization")
        assert "nvidia-smi" in block, "auto must check nvidia-smi"
        assert "if " in block and "else" in block, "auto must have if/else branch"

    def test_md_auto_has_nvidia_check(self):
        block = mdrun_block("md", "auto", stage="md")
        assert "nvidia-smi" in block

    def test_minimization_auto_gpu_branch_no_pme_gpu(self):
        block = mdrun_block("em", "auto", stage="minimization")
        # Parse the if-branch (before 'else')
        if_branch = block.split("else")[0]
        assert "-pme gpu" not in if_branch, (
            "GPU branch in auto/minimization must not use -pme gpu"
        )

    def test_md_auto_gpu_branch_has_pme_gpu(self):
        block = mdrun_block("md", "auto", stage="md")
        if_branch = block.split("else")[0]
        assert "-pme gpu" in if_branch


# ── Integration: workspace scripts match compatibility matrix ─────────────────

@pytest.fixture(scope="module")
def protein_workspace(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("stage_compat")
    result = SimulationCompiler().compile(PROTEIN_YAML)
    return WorkspaceBuilder().build(result, output_dir=str(tmp))


class TestWorkspaceScriptsCompat:
    def _step(self, workspace, name_fragment):
        return next(
            (p for p in (workspace / "steps").iterdir() if name_fragment in p.name),
            None,
        )

    # ── energy_minimization ───────────────────────────────────────────────────

    def test_em_script_no_pme_gpu(self, protein_workspace):
        em = self._step(protein_workspace, "energy_minimization")
        assert em is not None
        script = (em / "run.sh").read_text()
        assert "-pme gpu" not in script, (
            "energy_minimization/run.sh must not contain -pme gpu"
        )

    def test_em_script_no_bonded_gpu(self, protein_workspace):
        em = self._step(protein_workspace, "energy_minimization")
        script = (em / "run.sh").read_text()
        assert "-bonded gpu" not in script

    def test_em_script_pme_cpu_explicit(self, protein_workspace):
        em = self._step(protein_workspace, "energy_minimization")
        script = (em / "run.sh").read_text()
        assert "-pme cpu" in script, (
            "energy_minimization/run.sh must use -pme cpu in all branches"
        )

    # ── equilibration ─────────────────────────────────────────────────────────

    def test_equil_nvt_has_pme_gpu(self, protein_workspace):
        eq = self._step(protein_workspace, "equilibration")
        assert eq is not None
        script = (eq / "run_nvt.sh").read_text()
        # In auto-mode the GPU branch must have -pme gpu
        assert "-pme gpu" in script, (
            "equilibration/run_nvt.sh GPU branch must have -pme gpu"
        )

    def test_equil_npt_has_pme_gpu(self, protein_workspace):
        eq = self._step(protein_workspace, "equilibration")
        script = (eq / "run_npt.sh").read_text()
        assert "-pme gpu" in script

    # ── production ────────────────────────────────────────────────────────────

    def test_production_has_pme_gpu(self, protein_workspace):
        prod = self._step(protein_workspace, "production_md")
        assert prod is not None
        script = (prod / "run_md.sh").read_text()
        assert "-pme gpu" in script, (
            "production_md/run_md.sh GPU branch must have -pme gpu"
        )

    # ── cross-check: minimization script differs from production ─────────────

    def test_em_and_production_gpu_flags_differ(self, protein_workspace):
        em   = self._step(protein_workspace, "energy_minimization")
        prod = self._step(protein_workspace, "production_md")
        em_script   = (em   / "run.sh").read_text()
        prod_script = (prod / "run_md.sh").read_text()
        assert "-pme gpu" not in em_script
        assert "-pme gpu" in prod_script
