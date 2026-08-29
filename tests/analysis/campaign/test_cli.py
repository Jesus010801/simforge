"""`simforge study inspect|analyze|analyses` CLI surface (subprocess-driven).

Click 8.2 in this repo broke ``CliRunner(mix_stderr=...)`` so every case shells
out to ``python -m cli`` from the repo root, exactly like the analysis.md tests.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from tests.analysis.campaign.conftest import (
    REAL_GRO, REAL_TPR, REAL_XTC, REPO_ROOT, file_stat, requires_gmx,
    requires_real_traj, write_pdb, write_top,
)


def _cli(*args: str, timeout: int = 180) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "cli", *args],
        cwd=REPO_ROOT, capture_output=True, text=True, timeout=timeout,
    )


def _fake_study(tmp_path: Path) -> Path:
    root = tmp_path / "study"
    d = root / "GLP1" / "TIP3P" / "rep1"
    d.mkdir(parents=True)
    (d / "final.xtc").write_bytes(b"FAKEXTCDATA")
    write_top(d / "topol.top")
    write_pdb(d / "conf.pdb", {"A": 120})
    return root


# ── analyses ───────────────────────────────────────────────────────────────

def test_analyses_lists_four_ids():
    r = _cli("study", "analyses")
    assert r.returncode == 0
    for aid in ("rmsd-receptor", "rmsd-complex", "rmsd-peptide-intrinsic",
                "rmsd-peptide-receptor-frame"):
        assert aid in r.stdout


def test_analyses_json():
    r = _cli("study", "analyses", "--json")
    assert r.returncode == 0
    data = json.loads(r.stdout)
    assert isinstance(data, list)
    assert {e["id"] for e in data} == {
        "rmsd-receptor", "rmsd-complex", "rmsd-peptide-intrinsic",
        "rmsd-peptide-receptor-frame",
    }


# ── analyze: nothing selected ──────────────────────────────────────────────

def test_analyze_no_selection_exits_2_and_writes_nothing(tmp_path):
    study = _fake_study(tmp_path)
    r = _cli("study", "analyze", str(study), "--non-interactive")
    assert r.returncode == 2
    assert "No analysis selected" in (r.stdout + r.stderr)
    assert not (study / "simforge_analysis" / "systems").exists()


def test_analyze_bogus_id_exits_1_names_valid_ids(tmp_path):
    study = _fake_study(tmp_path)
    r = _cli("study", "analyze", str(study), "--analysis", "bogus", "--non-interactive")
    assert r.returncode == 1
    out = r.stdout + r.stderr
    assert "rmsd-receptor" in out


# ── analyze: dry-run ───────────────────────────────────────────────────────

@requires_gmx
@requires_real_traj
def test_analyze_dry_run_plans_without_executing(tmp_path):
    root = tmp_path / "study"
    d = root / "GLP1R" / "rep1"
    d.mkdir(parents=True)
    (d / "trajectory_x.xtc").symlink_to(REAL_XTC)
    (d / "topol_x.tpr").symlink_to(REAL_TPR)
    (d / "ref_x.gro").symlink_to(REAL_GRO)

    r = _cli("study", "analyze", str(root), "--analysis", "rmsd-receptor",
             "--non-interactive", "--dry-run", "--json")
    assert r.returncode == 0, r.stderr
    data = json.loads(r.stdout)
    planned = [res for res in data["results"] if res["status"] == "planned"]
    assert planned
    assert planned[0]["data_summary"].get("planned_command")
    # no measurement output written
    assert list(root.rglob("rmsd.xvg")) == []


# ── inspect ────────────────────────────────────────────────────────────────

def test_inspect_writes_manifest_and_report(tmp_path):
    study = _fake_study(tmp_path)
    r = _cli("study", "inspect", str(study), "--no-trajectory-inspection")
    assert r.returncode == 0, r.stderr
    out = study / "simforge_analysis"
    for name in ("study_manifest.yaml", "study_manifest.json",
                 "validation_report.json", "validation_report.txt"):
        assert (out / name).is_file(), name
    json.loads((out / "study_manifest.json").read_text())


def test_inspect_json_has_manifest_stats(tmp_path):
    study = _fake_study(tmp_path)
    r = _cli("study", "inspect", str(study), "--no-trajectory-inspection", "--json")
    assert r.returncode == 0
    data = json.loads(r.stdout)
    assert data["manifest_stats"]["n_systems"] == 1


# ── legacy XVG mode still works ────────────────────────────────────────────

def test_legacy_xvg_mode_not_dispatched_to_campaign(tmp_path):
    empty = tmp_path / "xvgdir"
    empty.mkdir()
    r = _cli("study", str(empty))
    combined = r.stdout + r.stderr
    # legacy analyzer output, NOT the campaign "SimForge Study" manifest header
    assert "XVG" in combined or "No study structure" in combined
    assert "study_manifest.yaml" not in combined


# ── originals are never modified ───────────────────────────────────────────

def test_source_files_untouched_after_analyze(tmp_path):
    study = _fake_study(tmp_path)
    sources = list(study.rglob("*.xtc")) + list(study.rglob("*.top")) + list(study.rglob("*.pdb"))
    before = {p: file_stat(p) for p in sources}

    _cli("study", "analyze", str(study), "--analysis", "rmsd-receptor",
         "--non-interactive")

    after = {p: file_stat(p) for p in sources}
    assert before == after
