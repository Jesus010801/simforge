"""Tests for checkpoint-aware recovery in RuntimeExecutor."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from executors.execution_state import StepStatus


# ── _detect_checkpoint ────────────────────────────────────────────────────────

class TestDetectCheckpoint:
    """Unit tests for the static _detect_checkpoint helper."""

    def test_no_files_returns_none(self, tmp_path):
        from runtime.executor import RuntimeExecutor
        assert RuntimeExecutor._detect_checkpoint(tmp_path) is None

    def test_md_cpt_without_gro_returns_md_resume(self, tmp_path):
        from runtime.executor import RuntimeExecutor
        (tmp_path / "md.cpt").touch()
        result = RuntimeExecutor._detect_checkpoint(tmp_path)
        assert result is not None
        cpt_path, script = result
        assert cpt_path.name == "md.cpt"
        assert script == "run_md_resume.sh"

    def test_md_cpt_with_gro_returns_none(self, tmp_path):
        from runtime.executor import RuntimeExecutor
        (tmp_path / "md.cpt").touch()
        (tmp_path / "md.gro").touch()
        assert RuntimeExecutor._detect_checkpoint(tmp_path) is None

    def test_nvt_cpt_without_gro_returns_nvt_resume(self, tmp_path):
        from runtime.executor import RuntimeExecutor
        (tmp_path / "nvt.cpt").touch()
        result = RuntimeExecutor._detect_checkpoint(tmp_path)
        assert result is not None
        _, script = result
        assert script == "run_nvt_resume.sh"

    def test_npt_cpt_without_gro_returns_npt_resume(self, tmp_path):
        from runtime.executor import RuntimeExecutor
        (tmp_path / "npt.cpt").touch()
        result = RuntimeExecutor._detect_checkpoint(tmp_path)
        assert result is not None
        _, script = result
        assert script == "run_npt_resume.sh"

    def test_npt_takes_priority_over_nvt(self, tmp_path):
        """Both nvt.cpt and npt.cpt present → prefer npt (more advanced)."""
        from runtime.executor import RuntimeExecutor
        (tmp_path / "nvt.cpt").touch()
        (tmp_path / "nvt.gro").touch()   # NVT completed
        (tmp_path / "npt.cpt").touch()   # NPT interrupted
        result = RuntimeExecutor._detect_checkpoint(tmp_path)
        assert result is not None
        _, script = result
        assert script == "run_npt_resume.sh"

    def test_nvt_cpt_with_gro_skipped_falls_to_next(self, tmp_path):
        """nvt.cpt + nvt.gro (done) + no npt.cpt → no recovery needed."""
        from runtime.executor import RuntimeExecutor
        (tmp_path / "nvt.cpt").touch()
        (tmp_path / "nvt.gro").touch()
        assert RuntimeExecutor._detect_checkpoint(tmp_path) is None


# ── mdrun_resume_block ────────────────────────────────────────────────────────

class TestMdrunResumeBlock:
    def test_cpu_contains_cpi_append(self):
        from builders.step_builders._utils import mdrun_resume_block
        block = mdrun_resume_block("md", hardware="cpu", stage="md")
        assert "-cpi md.cpt" in block
        assert "-append" in block
        assert "gmx mdrun" in block

    def test_gpu_contains_cpi_append(self):
        from builders.step_builders._utils import mdrun_resume_block
        block = mdrun_resume_block("md", hardware="gpu", stage="md")
        assert "-cpi md.cpt" in block
        assert "-append" in block

    def test_auto_mode_contains_cpi_in_both_branches(self):
        from builders.step_builders._utils import mdrun_resume_block
        block = mdrun_resume_block("nvt", hardware="auto", stage="md")
        assert block.count("-cpi nvt.cpt") == 2  # once in GPU, once in CPU branch
        assert block.count("-append") == 2

    def test_deffnm_propagated(self):
        from builders.step_builders._utils import mdrun_resume_block
        block = mdrun_resume_block("npt", hardware="cpu", stage="md")
        assert "-deffnm npt" in block
        assert "-cpi npt.cpt" in block


# ── ProductionBuilder generates run_md_resume.sh ─────────────────────────────

class TestProductionBuilderResumeScript:
    def _build(self, tmp_path):
        from core.execution_models import SimulationStep, StepStage, StepType
        from builders.step_builders.production_builder import ProductionBuilder

        step = SimulationStep(
            step_id="production_md",
            title="Production MD",
            stage=StepStage.PRODUCTION,
            step_type=StepType.AUTOMATIC,
            engine="gromacs",
            depends_on=[],
            params={"hardware": "cpu"},
        )
        ProductionBuilder().build(step, tmp_path, step_dir_map={})

    def test_run_md_resume_sh_generated(self, tmp_path):
        self._build(tmp_path)
        assert (tmp_path / "run_md_resume.sh").exists()

    def test_resume_script_has_no_grompp(self, tmp_path):
        self._build(tmp_path)
        content = (tmp_path / "run_md_resume.sh").read_text()
        assert "gmx grompp" not in content

    def test_resume_script_has_cpi_append(self, tmp_path):
        self._build(tmp_path)
        content = (tmp_path / "run_md_resume.sh").read_text()
        assert "-cpi md.cpt" in content
        assert "-append" in content

    def test_main_script_still_has_grompp(self, tmp_path):
        self._build(tmp_path)
        content = (tmp_path / "run_md.sh").read_text()
        assert "grompp" in content


# ── EquilibrationBuilder generates resume scripts ─────────────────────────────

class TestEquilibrationBuilderResumeScripts:
    def _build(self, tmp_path):
        from core.execution_models import SimulationStep, StepStage, StepType
        from builders.step_builders.equilibration_builder import EquilibrationBuilder

        step = SimulationStep(
            step_id="equilibration",
            title="Equilibration",
            stage=StepStage.EQUILIBRATION,
            step_type=StepType.AUTOMATIC,
            engine="gromacs",
            depends_on=[],
            params={"hardware": "cpu"},
        )
        EquilibrationBuilder().build(step, tmp_path, step_dir_map={})

    def test_nvt_resume_sh_generated(self, tmp_path):
        self._build(tmp_path)
        assert (tmp_path / "run_nvt_resume.sh").exists()

    def test_npt_resume_sh_generated(self, tmp_path):
        self._build(tmp_path)
        assert (tmp_path / "run_npt_resume.sh").exists()

    def test_nvt_resume_has_cpi_nvt(self, tmp_path):
        self._build(tmp_path)
        content = (tmp_path / "run_nvt_resume.sh").read_text()
        assert "-cpi nvt.cpt" in content
        assert "-append" in content
        assert "gmx grompp" not in content

    def test_npt_resume_has_cpi_npt(self, tmp_path):
        self._build(tmp_path)
        content = (tmp_path / "run_npt_resume.sh").read_text()
        assert "-cpi npt.cpt" in content
        assert "-append" in content
        assert "gmx grompp" not in content


# ── _plan_resume checkpoint detection (integration) ──────────────────────────

def _make_workspace(tmp_path: Path, step_id: str = "production_md") -> Path:
    """Create a minimal workspace with one step for testing _plan_resume."""
    steps_dir = tmp_path / "steps" / f"01_{step_id}"
    steps_dir.mkdir(parents=True)
    metadata_dir = tmp_path / "metadata"
    metadata_dir.mkdir()

    # Minimal metadata.json for the step
    (steps_dir / "metadata.json").write_text(json.dumps({
        "step_id": step_id,
        "step_type": "automatic",
        "expected_outputs": ["md.xtc", "md.edr", "md.log", "md.gro", "md.cpt"],
        "required_inputs": [],
        "params": {},
    }))

    # Minimal execution_state.json
    state = {
        "workspace_path": str(tmp_path),
        "is_complete": False,
        "dry_run": False,
        "steps": [{
            "step_id": step_id,
            "step_dir": str(steps_dir),
            "depends_on": [],
            "status": "running",
        }],
    }
    (tmp_path / "execution_state.json").write_text(json.dumps(state))

    # Minimal manifest so BaseExecutor loads
    manifest = {
        "workspace_name": "test",
        "steps": [{
            "step_id": step_id,
            "step_dir": str(steps_dir),
        }],
        "build_signature": "abc",
        "builder_signature": "abc",
    }
    (metadata_dir / "execution_manifest.json").write_text(json.dumps(manifest))

    return tmp_path


class TestPlanResumeCheckpoint:
    def test_step_with_cpt_no_gro_marked_recoverable(self, tmp_path):
        from runtime.executor import RuntimeExecutor

        ws = _make_workspace(tmp_path)
        step_dir = ws / "steps" / "01_production_md"
        (step_dir / "md.cpt").touch()          # checkpoint present
        (step_dir / "run_md_resume.sh").touch()  # resume script exists

        executor = RuntimeExecutor(ws, dry_run=False)
        executor.state = executor._initialize_state()
        plan = executor._plan_resume()

        assert "production_md" in plan.recoverable
        assert "production_md" in executor._resumable_steps
        assert executor._resumable_steps["production_md"] == "run_md_resume.sh"

    def test_step_with_cpt_and_gro_not_recoverable(self, tmp_path):
        from runtime.executor import RuntimeExecutor

        ws = _make_workspace(tmp_path)
        step_dir = ws / "steps" / "01_production_md"
        (step_dir / "md.cpt").touch()
        (step_dir / "md.gro").touch()   # gro present → simulation completed

        executor = RuntimeExecutor(ws, dry_run=False)
        executor.state = executor._initialize_state()
        plan = executor._plan_resume()

        assert "production_md" not in plan.recoverable

    def test_recoverable_step_gets_recoverable_status(self, tmp_path):
        from runtime.executor import RuntimeExecutor

        ws = _make_workspace(tmp_path)
        step_dir = ws / "steps" / "01_production_md"
        (step_dir / "md.cpt").touch()
        (step_dir / "run_md_resume.sh").touch()

        executor = RuntimeExecutor(ws, dry_run=False)
        executor.state = executor._initialize_state()
        executor._plan_resume()

        record = executor.state.get_step("production_md")
        assert record is not None
        assert record.status == StepStatus.RECOVERABLE

    def test_find_script_prefers_resume_for_recoverable(self, tmp_path):
        from runtime.executor import RuntimeExecutor

        ws = _make_workspace(tmp_path)
        step_dir = ws / "steps" / "01_production_md"
        (step_dir / "md.cpt").touch()
        (step_dir / "run_md.sh").touch()
        (step_dir / "run_md_resume.sh").touch()

        executor = RuntimeExecutor(ws, dry_run=False)
        executor._resumable_steps["production_md"] = "run_md_resume.sh"

        script = executor._find_script(step_dir, "production_md")
        assert script is not None
        assert script.name == "run_md_resume.sh"

    def test_find_script_falls_back_when_resume_missing(self, tmp_path):
        from runtime.executor import RuntimeExecutor

        ws = _make_workspace(tmp_path)
        step_dir = ws / "steps" / "01_production_md"
        (step_dir / "run_md.sh").touch()
        # run_md_resume.sh NOT created

        executor = RuntimeExecutor(ws, dry_run=False)
        executor._resumable_steps["production_md"] = "run_md_resume.sh"

        script = executor._find_script(step_dir, "production_md")
        assert script is not None
        assert script.name == "run_md.sh"

    def test_step_with_no_checkpoint_not_recoverable(self, tmp_path):
        from runtime.executor import RuntimeExecutor

        ws = _make_workspace(tmp_path)
        # No cpt file at all

        executor = RuntimeExecutor(ws, dry_run=False)
        executor.state = executor._initialize_state()
        plan = executor._plan_resume()

        assert "production_md" not in plan.recoverable
