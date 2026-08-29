"""Regression tests for MD-workflow trajectory-stage classification.

The motivating real bug: a directory containing ``nvt.xtc`` + ``npt.xtc`` +
``md.xtc`` was grouped as three production segments of one system.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from analysis.campaign.discovery import discover_study
from analysis.campaign.manifest import assigned_paths, build_manifest, load_manifest, write_manifest
from analysis.campaign.models import TrajectoryStage, ValidationState
from analysis.campaign.orchestration.study_analyzer import run_inspect
from analysis.campaign.trajectory.stage_classifier import classify_trajectory_stage
from analysis.campaign.validator import validate_study
from tests.analysis.campaign.conftest import write_pdb, write_top

_NVT_MDP = "integrator = md\ndt = 0.002\nnsteps = 25000\nnstxtcout = 10\ngen_vel = yes\n"
_NPT_MDP = ("integrator = md\ndt = 0.002\nnsteps = 50000\nnstxtcout = 100\n"
            "gen_vel = no\nPcoupl = Parrinello-Rahman\n")
_MD_MDP = ("integrator = md\ndt = 0.001\nnsteps = 25000000\n"
           "nstxout-compressed = 10000\ngen_vel = no\n"
           "tcoupl = V-rescale\npcoupl = Parrinello-Rahman\n")
_EM_MDP = "integrator = steep\nnsteps = 50000\nemtol = 1000\n"


def _standard_workflow_dir(tmp_path: Path) -> Path:
    d = tmp_path / "AA-A1"
    d.mkdir(parents=True)
    for stem, mdp in (("nvt", _NVT_MDP), ("npt", _NPT_MDP), ("md", _MD_MDP), ("em", _EM_MDP)):
        (d / f"{stem}.xtc").write_bytes(b"FAKE")
        (d / f"{stem}.tpr").write_bytes(b"FAKE")
        (d / f"{stem}.gro").write_bytes(b"FAKE")
        (d / f"{stem}.mdp").write_text(mdp)
    write_top(d / "topol.top")
    write_pdb(d / "md.gro" if False else d / "conf.pdb", {"A": 120})
    return d


# ── A. standard nvt/npt/md workflow directory ───────────────────────────────

def test_standard_workflow_one_production_two_equilibration(tmp_path):
    root = _standard_workflow_dir(tmp_path)
    res = discover_study(root)
    assert len(res.candidates) == 1
    cand = res.candidates[0]

    prod = [Path(p).name for p in cand.production_trajectory_paths]
    assert prod == ["md.xtc"], f"expected only md.xtc as production, got {prod}"

    stages = {Path(a.path).name: a.stage for a in cand.trajectory_artifacts}
    assert stages["md.xtc"] == TrajectoryStage.PRODUCTION
    assert stages["nvt.xtc"] == TrajectoryStage.EQUILIBRATION_NVT
    assert stages["npt.xtc"] == TrajectoryStage.EQUILIBRATION_NPT
    assert stages["em.xtc"] == TrajectoryStage.MINIMIZATION

    wf = {Path(p).name for p in cand.workflow_trajectory_paths}
    assert wf == {"nvt.xtc", "npt.xtc", "em.xtc"}

    assert not any("trajectory segments in one directory" in e
                   for e in cand.discovery_evidence)


def test_standard_workflow_stage_matched_topology(tmp_path):
    root = _standard_workflow_dir(tmp_path)
    cand = discover_study(root).candidates[0]
    assert Path(cand.topology_path).name == "md.tpr"
    assert not any(w.code.startswith("ambiguous") for w in cand.warnings)
    # each workflow artifact keeps its own stage-matched tpr for provenance
    by_name = {Path(a.path).name: a for a in cand.trajectory_artifacts}
    assert Path(by_name["npt.xtc"].topology_path).name == "npt.tpr"


def test_standard_workflow_no_inconsistent_timestep_warning(tmp_path):
    root = _standard_workflow_dir(tmp_path)
    m = build_manifest(root, inspect_trajectories=False, gmx="/nonexistent")
    codes = {w.code for s in m.systems for w in s.warnings}
    assert "inconsistent_timestep" not in codes
    assert "possibly_not_production" not in codes


# ── B. arbitrary filename, mdp decides ──────────────────────────────────────

def test_mdp_overrides_misleading_filename(tmp_path):
    d = tmp_path / "s"
    d.mkdir()
    # a file called "md.xtc" whose mdp is clearly NVT equilibration
    (d / "md.xtc").write_bytes(b"F")
    (d / "md.mdp").write_text(_NVT_MDP)
    r = classify_trajectory_stage(d / "md.xtc")
    assert r.stage == TrajectoryStage.EQUILIBRATION_NVT
    assert any("gen_vel" in e for e in r.evidence)


def test_short_production_not_misclassified_as_equilibration(tmp_path):
    """A short (<1 ns) production run named md.xtc, alone in its dir, with
    gen_vel=no — must stay production, not be downgraded to NPT equilibration."""
    d = tmp_path / "s"
    d.mkdir()
    (d / "md.xtc").write_bytes(b"F")
    (d / "md.mdp").write_text(
        "integrator = md\ndt = 0.002\nnsteps = 250000\nnstxout-compressed = 5000\n"
        "gen_vel = no\npcoupl = Parrinello-Rahman\n")
    r = classify_trajectory_stage(d / "md.xtc")
    assert r.stage == TrajectoryStage.PRODUCTION


def test_generic_name_no_mdp_is_production_by_elimination(tmp_path):
    root = tmp_path / "study" / "GLP1" / "TIP3P" / "rep1"
    root.mkdir(parents=True)
    (root / "final.xtc").write_bytes(b"F")
    write_top(root / "topol.top")
    write_pdb(root / "conf.pdb", {"A": 50})
    cand = discover_study(tmp_path / "study").candidates[0]
    assert [Path(p).name for p in cand.production_trajectory_paths] == ["final.xtc"]
    # single trajectory, no equilibration siblings -> benign INFO, not a blocking warning
    unverified = [w for w in cand.warnings if w.code == "production_stage_unverified"]
    assert unverified and unverified[0].severity == "info"


# ── C. true production continuation segments ────────────────────────────────

def test_true_production_parts_one_collection(tmp_path):
    d = tmp_path / "s"
    d.mkdir()
    for i in (1, 2, 3):
        (d / f"md.part{i:04d}.xtc").write_bytes(b"F")
    (d / "md.mdp").write_text(_MD_MDP)
    (d / "md.tpr").write_bytes(b"F")
    write_pdb(d / "conf.pdb", {"A": 30})
    cand = discover_study(tmp_path).candidates[0]
    names = [Path(p).name for p in cand.production_trajectory_paths]
    assert names == ["md.part0001.xtc", "md.part0002.xtc", "md.part0003.xtc"]
    assert any("continuation segments" in e for e in cand.discovery_evidence)


# ── D/E. multiple production trajectories WITHOUT part markers ──────────────

def test_multiple_unmarked_production_trajectories_not_merged(tmp_path):
    d = tmp_path / "s"
    d.mkdir()
    (d / "prod_run1.xtc").write_bytes(b"F")
    (d / "prod_run2.xtc").write_bytes(b"F")
    write_pdb(d / "conf.pdb", {"A": 30})
    cand = discover_study(tmp_path).candidates[0]
    assert any(w.code == "ambiguous_production_trajectories"
               and w.severity == "review_required" for w in cand.warnings)


# ── F. derived trajectory beside production ─────────────────────────────────

def test_derived_trajectory_never_production_segment(tmp_path):
    d = tmp_path / "s"
    d.mkdir()
    (d / "md.xtc").write_bytes(b"F")
    (d / "md.mdp").write_text(_MD_MDP)
    (d / "md_nojump.xtc").write_bytes(b"F")
    (d / "md_center.xtc").write_bytes(b"F")
    (d / "md_fit.xtc").write_bytes(b"F")
    write_pdb(d / "conf.pdb", {"A": 30})
    cand = discover_study(tmp_path).candidates[0]
    assert [Path(p).name for p in cand.production_trajectory_paths] == ["md.xtc"]
    derived = {Path(p).name for p in cand.derived_trajectory_paths}
    assert derived == {"md_nojump.xtc", "md_center.xtc", "md_fit.xtc"}
    stages = {a.stage for a in cand.trajectory_artifacts
              if Path(a.path).name != "md.xtc"}
    assert stages == {TrajectoryStage.DERIVED}


# ── G. warning propagation ─────────────────────────────────────────────────

def test_warnings_propagate_to_validation_report(tmp_path):
    d = tmp_path / "s"
    d.mkdir()
    (d / "prod_a.xtc").write_bytes(b"F")
    (d / "prod_b.xtc").write_bytes(b"F")
    write_pdb(d / "conf.pdb", {"A": 30})
    m = build_manifest(tmp_path, inspect_trajectories=False, gmx="/nonexistent")
    rep = validate_study(m, [])
    assert rep.study_state in (ValidationState.WARNING, ValidationState.INVALID)
    assert rep.warnings, "system warnings must appear in validation_report"
    assert any(w.code == "ambiguous_production_trajectories" for w in rep.warnings)


def test_equilibration_only_dir_is_invalid(tmp_path):
    d = tmp_path / "s"
    d.mkdir()
    (d / "nvt.xtc").write_bytes(b"F")
    (d / "nvt.mdp").write_text(_NVT_MDP)
    (d / "npt.xtc").write_bytes(b"F")
    (d / "npt.mdp").write_text(_NPT_MDP)
    write_pdb(d / "conf.pdb", {"A": 30})
    m = build_manifest(tmp_path, inspect_trajectories=False, gmx="/nonexistent")
    assert m.systems[0].production_trajectory_paths == []
    rep = validate_study(m, ["rmsd-receptor"])
    assert rep.system_states[m.systems[0].system_id] == ValidationState.INVALID
    assert any(w.code == "no_production_trajectory" for w in m.systems[0].warnings)


# ── H. assigned / unassigned disjoint invariant ───────────────────────────

def test_assigned_and_unassigned_are_disjoint(tmp_path):
    root = _standard_workflow_dir(tmp_path)
    m = build_manifest(root, inspect_trajectories=False, gmx="/nonexistent")
    claimed = assigned_paths(m)
    assert not (set(m.unassigned_files) & claimed), (
        f"overlap: {set(m.unassigned_files) & claimed}"
    )
    # the selected production topology/structure must NOT be unassigned
    s = m.systems[0]
    assert s.topology_path not in m.unassigned_files
    assert s.structure_path not in m.unassigned_files
    for p in s.production_trajectory_paths:
        assert p not in m.unassigned_files


# ── I. full run_inspect orchestration integration ─────────────────────────

def test_run_inspect_orchestration_classifies_stages(tmp_path):
    root = _standard_workflow_dir(tmp_path)
    out = tmp_path / "out"
    result = run_inspect(root, output_dir=out, inspect_trajectories=False, gmx="/nonexistent")
    sysrec = result.manifest.systems[0]
    assert [Path(p).name for p in sysrec.production_trajectory_paths] == ["md.xtc"]
    assert {Path(p).name for p in sysrec.workflow_trajectory_paths} == {"nvt.xtc", "npt.xtc", "em.xtc"}

    manifest_yaml = json.loads((out / "study_manifest.json").read_text())
    art = {a["path"].split("/")[-1]: a["stage"]
           for a in manifest_yaml["systems"][0]["trajectory_artifacts"]}
    assert art["md.xtc"] == "production"
    assert art["nvt.xtc"] == "equilibration_nvt"
    assert art["npt.xtc"] == "equilibration_npt"

    # reload round-trips the new structure
    reloaded = load_manifest(out / "study_manifest.yaml")
    assert [Path(p).name for p in reloaded.systems[0].production_trajectory_paths] == ["md.xtc"]
    assert len(reloaded.systems[0].trajectory_artifacts) == 4


# ── J. CLI integration ────────────────────────────────────────────────────

def test_cli_inspect_reports_production_only(tmp_path):
    import subprocess
    import sys
    root = _standard_workflow_dir(tmp_path)
    out = tmp_path / "cliout"
    proc = subprocess.run(
        [sys.executable, "-m", "cli", "study", "inspect", str(root),
         "--output", str(out), "--no-trajectory-inspection", "--json"],
        capture_output=True, text=True, cwd=Path(__file__).resolve().parents[3],
    )
    assert proc.returncode == 0, proc.stderr
    data = json.loads(proc.stdout)
    sysrec = data["manifest_stats"]
    assert sysrec["n_systems"] == 1
    payload = json.loads((out / "study_manifest.json").read_text())
    s = payload["systems"][0]
    assert [p.split("/")[-1] for p in s["production_trajectory_paths"]] == ["md.xtc"]
