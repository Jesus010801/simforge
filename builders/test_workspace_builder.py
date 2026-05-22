# builders/test_workspace_builder.py
from __future__ import annotations

from pathlib import Path
import json


# ── Workspace structure ───────────────────────────────────────────────────────

def test_workspace_path_exists(workspace):
    assert Path(workspace).exists()


def test_workspace_has_steps_directory(workspace):
    assert (Path(workspace) / "steps").is_dir()


def test_workspace_has_execution_manifest(workspace):
    assert (Path(workspace) / "metadata" / "execution_manifest.json").exists()


def test_manifest_has_correct_step_count(workspace, compilation_result):
    manifest = json.loads(
        (Path(workspace) / "metadata" / "execution_manifest.json").read_text()
    )
    assert manifest["n_steps"] == len(compilation_result.execution_order)


def test_manifest_steps_have_depends_on(workspace):
    manifest = json.loads(
        (Path(workspace) / "metadata" / "execution_manifest.json").read_text()
    )
    for entry in manifest["steps"]:
        assert "depends_on" in entry


def test_each_step_has_a_directory(workspace, compilation_result):
    steps_dir = Path(workspace) / "steps"
    dir_names = {d.name for d in steps_dir.iterdir() if d.is_dir()}
    for step in compilation_result.execution_order:
        assert any(step.step_id in name for name in dir_names)


def test_each_step_directory_has_metadata(workspace):
    steps_dir = Path(workspace) / "steps"
    for step_dir in sorted(steps_dir.iterdir()):
        if step_dir.is_dir():
            assert (step_dir / "metadata.json").exists(), \
                f"Missing metadata.json in {step_dir.name}"


# ── MDP files ─────────────────────────────────────────────────────────────────

def test_minimization_has_em_mdp(workspace):
    em_dir = next(
        (Path(workspace) / "steps" / d
         for d in _step_dirs(workspace)
         if "energy_minimization" in d),
        None,
    )
    assert em_dir is not None
    assert (em_dir / "em.mdp").exists()


def test_equilibration_has_nvt_and_npt_mdp(workspace):
    eq_dir = _find_step_dir(workspace, "equilibration")
    assert (eq_dir / "nvt.mdp").exists()
    assert (eq_dir / "npt.mdp").exists()


def test_production_has_md_mdp(workspace):
    prod_dir = _find_step_dir(workspace, "production_md")
    assert (prod_dir / "md.mdp").exists()


def test_mdp_temperature_matches_policy(workspace, plan):
    prod_dir = _find_step_dir(workspace, "production_md")
    mdp = (prod_dir / "md.mdp").read_text()
    expected_temp = str(plan.workflow_policy.temperature_K)
    assert expected_temp in mdp


# ── Run scripts — inter-step paths ────────────────────────────────────────────

def test_equilibration_nvt_script_has_em_dir_var(workspace):
    eq_dir = _find_step_dir(workspace, "equilibration")
    script = (eq_dir / "run_nvt.sh").read_text()
    assert "EM_DIR=" in script
    assert "TOPOL_DIR=" in script


def test_production_script_has_eq_dir_var(workspace):
    prod_dir = _find_step_dir(workspace, "production_md")
    script = (prod_dir / "run_md.sh").read_text()
    assert "EQ_DIR=" in script
    assert "TOPOL_DIR=" in script


def test_solvate_script_has_assemble_dir_var(workspace):
    sol_dir = _find_step_dir(workspace, "solvate_system")
    script = (sol_dir / "run.sh").read_text()
    assert "ASSEMBLE_DIR=" in script


# ── Analysis scripts ──────────────────────────────────────────────────────────

def test_rmsd_script_uses_gmx_rms(workspace):
    script = _read_analysis_script(workspace, "analysis_rmsd")
    assert "gmx rms" in script


def test_rmsd_script_uses_backbone_group(workspace):
    script = _read_analysis_script(workspace, "analysis_rmsd")
    assert '"4 4"' in script or "4 4" in script


def test_hbond_script_uses_gmx_hbond(workspace):
    script = _read_analysis_script(workspace, "analysis_hydrogen_bonds")
    assert "gmx hbond" in script


def test_distance_script_uses_gmx_distance(workspace):
    script = _read_analysis_script(workspace, "analysis_distance_analysis")
    assert "gmx distance" in script


def test_distance_step_has_make_ndx_script(workspace):
    dist_dir = _find_step_dir(workspace, "analysis_distance_analysis")
    assert (dist_dir / "make_ndx.sh").exists()


def test_analysis_scripts_reference_prod_dir(workspace):
    for name in ("analysis_rmsd", "analysis_hydrogen_bonds", "analysis_distance_analysis"):
        script = _read_analysis_script(workspace, name)
        assert "PROD_DIR=" in script


def test_analysis_metadata_has_gromacs_groups_field(workspace):
    for name in ("analysis_rmsd", "analysis_hydrogen_bonds"):
        d = _find_step_dir(workspace, name)
        meta = json.loads((d / "metadata.json").read_text())
        assert "gromacs_groups" in meta


# ── Helpers ───────────────────────────────────────────────────────────────────

def _step_dirs(workspace) -> list[str]:
    return [d.name for d in (Path(workspace) / "steps").iterdir() if d.is_dir()]


def _find_step_dir(workspace, step_id: str) -> Path:
    steps_dir = Path(workspace) / "steps"
    match = next(
        (d for d in steps_dir.iterdir() if d.is_dir() and step_id in d.name),
        None,
    )
    assert match is not None, f"Step directory for '{step_id}' not found"
    return match


def _read_analysis_script(workspace, step_id: str) -> str:
    d = _find_step_dir(workspace, step_id)
    return (d / "run_analysis.sh").read_text()
