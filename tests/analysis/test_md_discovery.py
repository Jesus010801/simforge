"""Tests for analysis.md — MD run discovery, observable planning, and CLI routing.

All tests use synthetic temporary directories; no GROMACS or network access.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from analysis.md.discovery import discover_md_run
from analysis.md.models import DiscoveredMDRun
from analysis.md.observables import plan_observables

# Repo root — used to set cwd for subprocess calls so cli.py is importable
_REPO = Path(__file__).parent.parent.parent


# ── Helpers ───────────────────────────────────────────────────────────────────

def _fake_file(path: Path, content: str = "x") -> Path:
    """Write a tiny file so stat().st_size > 0."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    return path


def _complete_run(tmp_path: Path) -> Path:
    """Create a fake but structurally complete GROMACS run directory."""
    for name in ("md.tpr", "md.xtc", "md.gro", "md.edr", "topol.top", "index.ndx", "md.log"):
        _fake_file(tmp_path / name)
    return tmp_path


def _cli(*args: str) -> subprocess.CompletedProcess:
    """Run the simforge CLI via `python -m cli` from the repo root."""
    return subprocess.run(
        [sys.executable, "-m", "cli", *args],
        capture_output=True,
        text=True,
        cwd=str(_REPO),
    )


# ── Test 1: complete fake run detected correctly ───────────────────────────────

def test_discover_complete_run(tmp_path):
    _complete_run(tmp_path)
    run = discover_md_run(tmp_path)

    assert run.tpr        is not None and run.tpr.name        == "md.tpr"
    assert run.trajectory is not None and run.trajectory.name == "md.xtc"
    assert run.structure  is not None and run.structure.name  == "md.gro"
    assert run.edr        is not None and run.edr.name        == "md.edr"
    assert run.topology   is not None and run.topology.name   == "topol.top"
    assert run.index      is not None and run.index.name      == "index.ndx"
    assert run.log        is not None and run.log.name        == "md.log"
    assert len(run.all_files) == 7

    warn_codes = {w.code for w in run.warnings}
    assert "no_tpr"        not in warn_codes
    assert "no_trajectory" not in warn_codes
    assert "no_edr"        not in warn_codes


# ── Test 2: missing .tpr generates warning ────────────────────────────────────

def test_discover_warns_no_tpr(tmp_path):
    for name in ("md.xtc", "md.gro", "md.edr"):
        _fake_file(tmp_path / name)
    run = discover_md_run(tmp_path)

    assert run.tpr is None
    assert any(w.code == "no_tpr" for w in run.warnings)


# ── Test 3: .xtc preferred over .trr ─────────────────────────────────────────

def test_discover_prefers_xtc_over_trr(tmp_path):
    _fake_file(tmp_path / "md.xtc", "xtc_data")
    _fake_file(tmp_path / "md.trr", "trr_data")
    _fake_file(tmp_path / "md.tpr")
    run = discover_md_run(tmp_path)

    assert run.trajectory is not None
    assert run.trajectory.suffix == ".xtc"


# ── Test 4: multiple trajectories recorded ────────────────────────────────────

def test_discover_multiple_trajectories(tmp_path):
    _fake_file(tmp_path / "md.xtc",  "a")
    _fake_file(tmp_path / "md2.xtc", "b")
    _fake_file(tmp_path / "md.trr",  "c")
    _fake_file(tmp_path / "md.tpr")
    run = discover_md_run(tmp_path)

    assert len(run.trajectory_candidates) == 3
    assert any(w.code == "multiple_trajectories" for w in run.warnings)
    assert run.trajectory is not None
    assert run.trajectory.suffix == ".xtc"


# ── Test 5: RMSD/RMSF/Rg/SASA planned when tpr+traj present ─────────────────

def test_plan_trajectory_observables_when_complete(tmp_path):
    _complete_run(tmp_path)
    run  = discover_md_run(tmp_path)
    plan = plan_observables(run)

    status_by_name = {o.name: o.status for o in plan.observables}
    for name in ("rmsd", "rmsf", "radius_of_gyration", "sasa"):
        assert status_by_name[name] == "planned", f"{name} should be planned"


# ── Test 6: energy QC unavailable when .edr missing ──────────────────────────

def test_plan_energy_unavailable_without_edr(tmp_path):
    _fake_file(tmp_path / "md.tpr")
    _fake_file(tmp_path / "md.xtc")
    run  = discover_md_run(tmp_path)
    plan = plan_observables(run)

    status_by_name = {o.name: o.status for o in plan.observables}
    assert status_by_name["energy_temperature_pressure_qc"] == "unavailable"


# ── Test 7: CLI routing — `analyze md` creates all expected output files ──────

def test_cli_analyze_md_creates_output_files(tmp_path):
    run_dir = tmp_path / "run"
    out_dir = tmp_path / "analysis"
    _complete_run(run_dir)

    result = _cli("analyze", "md", str(run_dir), "--out", str(out_dir), "--dry-run")
    assert result.returncode == 0, f"stderr: {result.stderr}\nstdout: {result.stdout}"

    assert (out_dir / "metadata" / "run_discovery.json").exists()
    assert (out_dir / "metadata" / "analysis_plan.json").exists()
    assert (out_dir / "metadata" / "provenance.json").exists()
    assert (out_dir / "report.md").exists()


# ── Test 8: CLI routing — `analyze md --help` shows --out and --dry-run ───────

def test_cli_analyze_md_help(tmp_path):
    result = _cli("analyze", "md", "--help")
    assert result.returncode == 0, result.stderr
    assert "--out" in result.stdout or "--out" in result.stderr
    assert "--dry-run" in result.stdout or "--dry-run" in result.stderr


# ── Test 9: Legacy `analyze <path>` is still routed correctly ─────────────────

def test_cli_analyze_legacy_routing(tmp_path):
    """Verify `simforge analyze <path>` still invokes the XVG quality command, not the MD sub-command."""
    # The legacy command exits with an analysis error (no XVG files), not a routing error.
    # A routing error would give "No such option: --out" or "No such command: analyze".
    result = _cli("analyze", str(tmp_path))
    # It should NOT error with option/command routing issues
    assert "No such option: --out" not in result.stderr
    assert "No such command" not in result.stderr
    # May exit 1 (no XVG data) but not crash with usage error from wrong routing
    assert result.returncode in (0, 1)


# ── Test 10: CLI dry-run does not invoke GROMACS ──────────────────────────────

def test_cli_analyze_md_no_gromacs(tmp_path):
    """Verify no subprocess call to gmx occurs during a dry-run."""
    run_dir = tmp_path / "run"
    out_dir = tmp_path / "analysis"
    _complete_run(run_dir)

    result = _cli("analyze", "md", str(run_dir), "--out", str(out_dir), "--dry-run")
    assert result.returncode == 0, result.stderr
    combined = result.stdout + result.stderr
    assert "No such option" not in combined
    assert "No such command" not in combined


# ── Test 11: JSON outputs valid and contain warnings ──────────────────────────

def test_cli_analyze_md_json_valid_with_warnings(tmp_path):
    run_dir = tmp_path / "run"
    out_dir = tmp_path / "analysis"
    # Deliberately omit .tpr and .edr to force warnings
    _fake_file(run_dir / "md.xtc")
    _fake_file(run_dir / "md.gro")

    result = _cli("analyze", "md", str(run_dir), "--out", str(out_dir), "--dry-run")
    assert result.returncode == 0, f"stderr: {result.stderr}"

    discovery = json.loads((out_dir / "metadata" / "run_discovery.json").read_text())
    plan_data = json.loads((out_dir / "metadata" / "analysis_plan.json").read_text())
    prov_data = json.loads((out_dir / "metadata" / "provenance.json").read_text())

    assert "run_dir"     in discovery
    assert "warnings"    in discovery
    assert "observables" in plan_data
    assert "timestamp"   in prov_data
    assert prov_data["dry_run"] is True

    warn_codes = {w["code"] for w in discovery["warnings"]}
    assert "no_tpr" in warn_codes
    assert "no_edr" in warn_codes

    obs_status = {o["name"]: o["status"] for o in plan_data["observables"]}
    assert obs_status["rmsd"] == "unavailable"
    assert obs_status["energy_temperature_pressure_qc"] == "unavailable"


# ── Test 12: explicit --trajectory/--topology override selects user files ──────

def test_cli_analyze_md_explicit_overrides_top_level(tmp_path):
    """User-specified --trajectory and --topology must win over auto-discovered files."""
    run_dir = tmp_path / "run"
    out_dir = tmp_path / "analysis"

    # Populate: auto-discovery would pick R-1/md.xtc and R-1/md.tpr
    _fake_file(run_dir / "R-1" / "md.xtc")
    _fake_file(run_dir / "R-1" / "md.tpr")
    # Top-level files the user wants to use explicitly
    _fake_file(run_dir / "md_centered.xtc")
    _fake_file(run_dir / "md.tpr")
    _fake_file(run_dir / "md.edr")
    _fake_file(run_dir / "md.gro")

    result = _cli(
        "analyze", "md", str(run_dir),
        "--trajectory", "md_centered.xtc",
        "--topology", "md.tpr",
        "--out", str(out_dir),
        "--dry-run",
    )
    assert result.returncode == 0, f"stderr: {result.stderr}\nstdout: {result.stdout}"

    report = (out_dir / "report.md").read_text()
    plan_data = json.loads((out_dir / "metadata" / "analysis_plan.json").read_text())

    # Report shows user-specified files
    assert "md_centered.xtc" in report
    assert "md.tpr" in report
    assert "user-specified" in report

    # Auto-selected trajectory R-1/md.xtc must NOT appear as the selected trajectory
    assert "R-1/md.xtc" not in report.split("Primary trajectory")[1].split("\n")[0]

    # Proposed GROMACS commands must use the user-specified paths
    cmds = [o["gmx_command"] for o in plan_data["observables"] if o.get("gmx_command")]
    assert any("-f md_centered.xtc" in cmd for cmd in cmds), f"Commands: {cmds}"
    assert any("-s md.tpr" in cmd for cmd in cmds), f"Commands: {cmds}"
    # No command should use just the basename of the auto-discovered trajectory
    assert not any("-f md.xtc" in cmd for cmd in cmds), f"Commands leaked auto-discovered path: {cmds}"


# ── Test 13: --trajectory/--topology with subdirectory paths preserve them ─────

def test_cli_analyze_md_explicit_overrides_subdir(tmp_path):
    """When user specifies R-1/md.xtc explicitly, commands must keep R-1/md.xtc."""
    run_dir = tmp_path / "run"
    out_dir = tmp_path / "analysis"

    _fake_file(run_dir / "R-1" / "md.xtc")
    _fake_file(run_dir / "R-1" / "md.tpr")

    result = _cli(
        "analyze", "md", str(run_dir),
        "--trajectory", "R-1/md.xtc",
        "--topology", "R-1/md.tpr",
        "--out", str(out_dir),
        "--dry-run",
    )
    assert result.returncode == 0, f"stderr: {result.stderr}\nstdout: {result.stdout}"

    plan_data = json.loads((out_dir / "metadata" / "analysis_plan.json").read_text())
    cmds = [o["gmx_command"] for o in plan_data["observables"] if o.get("gmx_command")]

    # Paths must not be collapsed to bare basenames
    assert any("-f R-1/md.xtc" in cmd for cmd in cmds), f"Commands: {cmds}"
    assert any("-s R-1/md.tpr" in cmd for cmd in cmds), f"Commands: {cmds}"
    assert not any("-f md.xtc" in cmd for cmd in cmds), f"Path collapsed: {cmds}"
    assert not any("-s md.tpr" in cmd for cmd in cmds), f"Path collapsed: {cmds}"


# ── Test 14: multiple_trajectories warning suppressed when user specifies one ──

def test_discover_multiple_trajectories_user_override_changes_warning(tmp_path):
    """When user overrides trajectory, the warning should be info-level and mention the user file."""
    from analysis.md.discovery import discover_md_run

    _fake_file(tmp_path / "md.xtc",  "a")
    _fake_file(tmp_path / "md2.xtc", "b")
    _fake_file(tmp_path / "md.tpr")

    run = discover_md_run(tmp_path, overrides={"trajectory": "md2.xtc"})

    multi_warns = [w for w in run.warnings if w.code == "multiple_trajectories"]
    assert len(multi_warns) == 1
    assert multi_warns[0].severity == "info", "Warning should be downgraded to info"
    assert "md2.xtc" in multi_warns[0].message
    assert "user-specified" in multi_warns[0].message.lower()

    # The trajectory must be md2.xtc, not the auto-selected md.xtc
    assert run.trajectory is not None
    assert run.trajectory.name == "md2.xtc"
