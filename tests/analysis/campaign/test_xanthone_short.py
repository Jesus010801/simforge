"""Controlled execution of the xanthone_short thesis core.

Every GROMACS invocation is faked (``run_gmx`` monkeypatched): these tests
exercise eligibility gating, command construction, the missing-only state
machine, provenance, output-path safety and source integrity — never real MD.
The scenario matrix mirrors spec §8.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from analysis.campaign import xanthone_short as xs
from analysis.campaign.gmx import GmxResult
from analysis.campaign.manifest import build_manifest
from analysis.campaign.models import TrajectoryInspection
from tests.analysis.campaign.test_legacy import snapshot, system

WINDOW_END = 25000.0


# ── fixtures ─────────────────────────────────────────────────────────────────

def _ds(tmp_path):
    """Dataset root, kept separate from the (sibling) managed output dir."""
    root = tmp_path / "ds"
    root.mkdir(exist_ok=True)
    return root


def _out(tmp_path):
    return tmp_path / "out"


def _record(tmp_path, name="AA-A6", *, ligand_id=13, site_id=21,
            mdp_duration=25000, measured_end=25000.0, dt_ps=10.0,
            inspected=True):
    """A discovered SystemRecord for one synthetic A6 system with fake timing."""
    root = _ds(tmp_path)
    d, _ = system(root, name, ligand_id=ligand_id, site_id=site_id,
                  duration=mdp_duration)
    (d / "md.tpr").write_bytes(b"fake tpr header")
    rec = _rebuild(root, name, measured_end=measured_end, dt_ps=dt_ps,
                   inspected=inspected)
    return d, rec


def _rebuild(root, name, *, measured_end=25000.0, dt_ps=10.0, inspected=True):
    """Re-discover *name* under *root* (after an index edit) with fake timing."""
    manifest = build_manifest(root, inspect_trajectories=False)
    rec = next(s for s in manifest.systems
               if (s.legacy_study or {}).get("system") == name)
    rec.trajectory_inspection = TrajectoryInspection(
        trajectory_paths=list(rec.production_trajectory_paths),
        topology_path=rec.topology_path,
        inspected=inspected,
        end_time_ps=measured_end if inspected else None,
        dt_ps=dt_ps if inspected else None,
    )
    return rec


def _xvg(n_points=2501, *, title="RMSD", yunit="RMSD (nm)", start=0.0,
         end=WINDOW_END):
    step = (end - start) / (n_points - 1)
    head = (f'@    title "{title}"\n@    xaxis  label "Time (ps)"\n'
            f'@    yaxis  label "{yunit}"\n@TYPE xy\n')
    rows = "".join(f"{start + i * step:.3f}  {0.1 + 0.001 * (i % 7):.5f}\n"
                   for i in range(n_points))
    return head + rows


class FakeGmx:
    """Records calls and writes a valid XVG to every output flag it sees."""

    def __init__(self, *, fail=False, n_points=2501, end=WINDOW_END):
        self.calls: list[dict] = []
        self.fail = fail
        self.n_points = n_points
        self.end = end

    def __call__(self, args, *, stdin=None, gmx="gmx", cwd=None, timeout=0,
                 extra_env=None):
        argv = [gmx, *args]
        self.calls.append({"argv": argv, "stdin": stdin})
        if args[0] == "check":                      # trajectory probe
            return GmxResult(argv=argv, returncode=0, stdout="", stderr="")
        if args[0] == "dump":                        # TPR readability probe
            return GmxResult(argv=argv, returncode=0,
                             stdout="inputrec:\n   integrator = md\n", stderr="")
        if not self.fail:
            flags = {"-o", "-od", "-on", "-oall", "-ox"}
            it = iter(args)
            for a in it:
                if a in flags:
                    out = Path(next(it))
                    out.parent.mkdir(parents=True, exist_ok=True)
                    if a == "-on":
                        title, yunit = "Number of Contacts < 0.6 nm", "Number"
                    elif a in ("-od", "-oall", "-ox"):
                        title, yunit = "Distance", "Distance (nm)"
                    else:
                        title, yunit = "RMSD", "RMSD (nm)"
                    out.write_text(_xvg(self.n_points, title=title, yunit=yunit,
                                        end=self.end))
        rc = 1 if self.fail else 0
        return GmxResult(argv=argv, returncode=rc, stdout="",
                         stderr="fatal error" if self.fail else "",
                         gmx_version="2025.2")


@pytest.fixture
def fake_gmx(monkeypatch):
    fg = FakeGmx()
    monkeypatch.setattr(xs, "run_gmx", fg)
    monkeypatch.setattr(xs, "gmx_available", lambda _b="gmx": True)
    monkeypatch.setattr(xs, "gmx_version", lambda _b="gmx": "2025.2")
    return fg


def _run(rec, out_root, **kw):
    return xs.run_system(rec, out_root=Path(out_root), gmx="gmx",
                         force=kw.get("force", False),
                         dry_run=kw.get("dry_run", False))


def _by_analysis(run):
    return {a.analysis: a for a in run.analyses}


# ── command construction ─────────────────────────────────────────────────────

def test_window_and_cutoff_are_explicit(tmp_path):
    d, rec = _record(tmp_path)
    legacy = rec.legacy_study
    out = _out(tmp_path)
    for task in xs.TASKS:
        argv, _ = xs.build_command(task, tpr="md.tpr", traj="md.xtc",
                                   ndx=legacy["index"], legacy=legacy,
                                   out_dir=out, gmx="gmx")
        assert argv[argv.index("-b") + 1] == "0"
        assert argv[argv.index("-e") + 1] == "25000"
    mindist = next(t for t in xs.TASKS if t.tool == "mindist")
    argv, _ = xs.build_command(mindist, tpr="md.tpr", traj="md.xtc",
                               ndx=legacy["index"], legacy=legacy, out_dir=out)
    assert argv[argv.index("-d") + 1] == "0.6"


def test_groups_selected_by_name_not_id(tmp_path):
    # LIG at a non-standard position; command must still reference it by name.
    d, rec = _record(tmp_path, ligand_id=17)
    legacy = rec.legacy_study
    assert legacy["selections"]["ligand"]["group_id"] == 17
    lig_rmsd = next(t for t in xs.TASKS if t.key == "ligand_rmsd")
    argv, stdin = xs.build_command(lig_rmsd, tpr="s.tpr", traj="t.xtc",
                                   ndx=legacy["index"], legacy=legacy,
                                   out_dir=_out(tmp_path))
    assert stdin == "Protein\nLIG\n"
    assert "17" not in stdin
    cat = next(t for t in xs.TASKS if t.tool == "distance")
    argv, _ = xs.build_command(cat, tpr="s.tpr", traj="t.xtc",
                               ndx=legacy["index"], legacy=legacy,
                               out_dir=_out(tmp_path))
    sel = argv[argv.index("-select") + 1]
    assert sel == 'com of group "LIG" plus com of group "Catalytic_AA"'


def test_build_command_refuses_unresolved_group(tmp_path):
    """A None group id must never be turned into a GROMACS selection."""
    d, rec = _record(tmp_path)
    legacy = dict(rec.legacy_study)
    legacy["selections"] = dict(legacy["selections"])
    legacy["selections"]["catalytic_site"] = {"value": None, "group_id": None,
                                              "confidence": 0}
    cat = next(t for t in xs.TASKS if t.tool == "distance")
    with pytest.raises(ValueError, match="unresolved semantic group"):
        xs.build_command(cat, tpr="s.tpr", traj="t.xtc", ndx=legacy["index"],
                         legacy=legacy, out_dir=_out(tmp_path))


def test_equilibration_tpr_is_not_used_as_reference(tmp_path, fake_gmx):
    """Only md.tpr is a valid reference; em.tpr/npt.tpr must never substitute."""
    d, rec = _record(tmp_path)
    (d / "md.tpr").unlink()
    (d / "em.tpr").write_bytes(b"equilibration tpr")
    (d / "npt.tpr").write_bytes(b"npt tpr")
    rec = _rebuild(_ds(tmp_path), "AA-A6")
    run = _run(rec, _out(tmp_path))
    assert run.tpr is None
    assert {a.status for a in run.analyses} == {xs.REVIEW_REQUIRED}
    assert all("md.tpr not found" in a.message for a in run.analyses)
    # nothing was executed with the wrong reference
    assert not any(c["argv"][1] in ("rms", "mindist", "distance")
                   for c in fake_gmx.calls)


def test_catalytic_uses_com_distance_semantics(tmp_path):
    d, rec = _record(tmp_path)
    cat = next(t for t in xs.TASKS if t.tool == "distance")
    argv, stdin = xs.build_command(cat, tpr="s.tpr", traj="t.xtc",
                                   ndx=rec.legacy_study["index"],
                                   legacy=rec.legacy_study,
                                   out_dir=_out(tmp_path))
    assert "-oall" in argv and stdin == ""
    assert argv[1] == "distance"


# ── scenario matrix (spec §8) ────────────────────────────────────────────────

def test_trajectory_longer_than_window_executes_over_0_25(tmp_path, fake_gmx):
    # 200 ns production, folder claims 25 ns → override, window still 0-25 ns.
    d, rec = _record(tmp_path, "HMG-R-25ns-A6", ligand_id=13, site_id=21,
                     mdp_duration=200000, measured_end=200000.0)
    run = _run(rec, _out(tmp_path))
    assert run.eligibility["eligible"] is True
    assert run.eligibility["override_applied"] is True
    assert {a.status for a in run.analyses} == {xs.EXECUTED}
    for c in fake_gmx.calls:
        if c["argv"][1] in ("rms", "mindist", "distance"):
            assert c["argv"][c["argv"].index("-e") + 1] == "25000"


def test_exact_25ns_trajectory(tmp_path, fake_gmx):
    d, rec = _record(tmp_path, measured_end=25000.0, mdp_duration=25000)
    run = _run(rec, _out(tmp_path))
    assert {a.status for a in run.analyses} == {xs.EXECUTED}


def test_trajectory_too_short_is_blocked(tmp_path, fake_gmx):
    d, rec = _record(tmp_path, measured_end=12000.0, dt_ps=10.0)
    run = _run(rec, _out(tmp_path))
    assert {a.status for a in run.analyses} == {xs.REVIEW_REQUIRED}
    assert any("25" in a.message for a in run.analyses)


def test_missing_active_site_group_blocks_only_mindist(tmp_path, fake_gmx):
    d, rec = _record(tmp_path)
    idx = d / "index.ndx"
    idx.write_text(idx.read_text().replace("ActiveSite_AA", "SomethingElse"))
    rec = _rebuild(_ds(tmp_path), "AA-A6")
    run = _run(rec, _out(tmp_path))
    st = _by_analysis(run)
    assert st["active_site_mindist"].status == xs.REVIEW_REQUIRED
    assert st["active_site_contacts"].status == xs.REVIEW_REQUIRED
    # protein / ligand RMSD do not depend on the active-site group
    assert st["protein_rmsd"].status == xs.EXECUTED
    assert st["ligand_rmsd"].status == xs.EXECUTED


def test_missing_catalytic_group_blocks_only_catalytic(tmp_path, fake_gmx):
    d, rec = _record(tmp_path)
    idx = d / "index.ndx"
    idx.write_text(idx.read_text().replace("Catalytic_AA", "Nope"))
    rec = _rebuild(_ds(tmp_path), "AA-A6")
    run = _run(rec, _out(tmp_path))
    st = _by_analysis(run)
    assert st["catalytic_com_distance"].status == xs.REVIEW_REQUIRED
    assert st["protein_rmsd"].status == xs.EXECUTED
    assert st["active_site_mindist"].status == xs.EXECUTED


def test_lig_at_nonstandard_group_id_still_runs_and_is_recorded(tmp_path, fake_gmx):
    d, rec = _record(tmp_path, ligand_id=19)
    run = _run(rec, _out(tmp_path))
    st = _by_analysis(run)
    assert st["ligand_rmsd"].status == xs.EXECUTED
    prov = json.loads(Path(st["ligand_rmsd"].provenance_path).read_text())
    assert prov["numeric_group_ids"]["LIG"] == 19
    assert prov["semantic_groups"]["ligand"] == "LIG"


def test_existing_output_is_skipped(tmp_path, fake_gmx):
    d, rec = _record(tmp_path)
    out_root = _out(tmp_path)
    first = _run(rec, out_root)
    assert {a.status for a in first.analyses} == {xs.EXECUTED}
    # rerun: everything now present in the managed output tree
    second = _run(rec, out_root)
    assert {a.status for a in second.analyses} == {xs.SKIP_EXISTING}


def test_force_recomputes_existing(tmp_path, fake_gmx):
    d, rec = _record(tmp_path)
    out_root = _out(tmp_path)
    _run(rec, out_root)
    forced = xs.run_system(rec, out_root=out_root, gmx="gmx", force=True,
                           dry_run=False)
    assert {a.status for a in forced.analyses} == {xs.EXECUTED}


def test_existing_result_in_source_dir_is_skipped(tmp_path, fake_gmx):
    d, rec = _record(tmp_path)
    (d / "rmsd_protein.xvg").write_text(_xvg(title="RMSD"))
    run = _run(rec, _out(tmp_path))
    st = _by_analysis(run)
    assert st["protein_rmsd"].status == xs.SKIP_EXISTING
    assert st["ligand_rmsd"].status == xs.EXECUTED


def test_failed_gromacs_command_is_failed(tmp_path, monkeypatch):
    fg = FakeGmx(fail=True)
    monkeypatch.setattr(xs, "run_gmx", fg)
    monkeypatch.setattr(xs, "gmx_available", lambda _b="gmx": True)
    monkeypatch.setattr(xs, "gmx_version", lambda _b="gmx": "2025.2")
    d, rec = _record(tmp_path)
    run = _run(rec, _out(tmp_path))
    assert {a.status for a in run.analyses} == {xs.FAILED}
    # a failed run still writes provenance recording the non-zero return code
    prov = json.loads(Path(run.analyses[0].provenance_path).read_text())
    assert prov["gmx_returncode"] == 1


def test_output_path_collision_with_source_is_refused(tmp_path):
    d, rec = _record(tmp_path)
    with pytest.raises(ValueError, match="unsafe output path"):
        xs.assert_safe_output(tmp_path / "sub", tmp_path, [d])
    with pytest.raises(ValueError, match="unsafe output path"):
        xs.assert_safe_output(d, tmp_path, [d])


def test_rerun_determinism(tmp_path, fake_gmx):
    d, rec = _record(tmp_path)
    out_root = _out(tmp_path)
    run1 = _run(rec, out_root)
    argv1 = [a.argv for a in run1.analyses]
    # wipe managed outputs, run again
    import shutil
    shutil.rmtree(out_root)
    fg2 = FakeGmx()
    fake_gmx.calls.clear()
    run2 = _run(rec, out_root)
    argv2 = [a.argv for a in run2.analyses]
    assert argv1 == argv2


def test_source_files_unchanged(tmp_path, fake_gmx):
    d, rec = _record(tmp_path)
    before = snapshot(_ds(tmp_path))
    _run(rec, _out(tmp_path))
    assert snapshot(_ds(tmp_path)) == before


# ── provenance completeness (spec §6) ────────────────────────────────────────

def test_provenance_has_every_required_field(tmp_path, fake_gmx):
    d, rec = _record(tmp_path)
    run = _run(rec, _out(tmp_path))
    prov = json.loads(Path(_by_analysis(run)["protein_rmsd"].provenance_path).read_text())
    for key in ("source_trajectory", "source_tpr", "index_file",
                "semantic_groups", "numeric_group_ids", "analysis_window_ps",
                "command_argv", "gromacs_version", "timestamp", "outputs",
                "protocol_version", "profile", "source_fingerprints"):
        assert key in prov, key
    assert prov["analysis_window_ps"] == [0.0, 25000.0]
    assert prov["outputs"][0]["sha256"]
    assert prov["profile"] == "xanthone_short"


# ── post-run XVG validation (spec §10) ───────────────────────────────────────

def test_xvg_post_validation_flags_range_and_points(tmp_path, fake_gmx):
    d, rec = _record(tmp_path)
    run = _run(rec, _out(tmp_path))
    chk = _by_analysis(run)["protein_rmsd"].xvg_check
    assert chk["readable"] and chk["all_finite"]
    assert chk["n_points"] == 2501
    assert chk["time_start_ps"] == 0.0
    assert abs(chk["time_end_ps"] - 25000.0) < 1e-6
    assert chk["time_range_ok"] is True
    assert chk["points_match_10ps_sampling"] is True


def test_xvg_post_validation_notices_wrong_point_count(tmp_path, monkeypatch):
    fg = FakeGmx(n_points=1001, end=WINDOW_END)
    monkeypatch.setattr(xs, "run_gmx", fg)
    monkeypatch.setattr(xs, "gmx_available", lambda _b="gmx": True)
    monkeypatch.setattr(xs, "gmx_version", lambda _b="gmx": "2025.2")
    d, rec = _record(tmp_path)
    run = _run(rec, _out(tmp_path))
    chk = _by_analysis(run)["protein_rmsd"].xvg_check
    assert chk["points_match_10ps_sampling"] is False


# ── full study driver + integrity ───────────────────────────────────────────

def test_run_study_writes_reports_and_integrity(tmp_path, fake_gmx, monkeypatch):
    root = tmp_path / "dataset"
    root.mkdir()
    system(root, "AA-A6")
    (root / "AA-A6" / "md.tpr").write_bytes(b"tpr")

    # deterministic measured timing without touching real gmx check
    real_build = xs.build_manifest

    def patched(study_root, **kw):
        m = real_build(study_root, inspect_trajectories=False,
                       gmx=kw.get("gmx", "gmx"))
        for s in m.systems:
            s.trajectory_inspection = TrajectoryInspection(
                trajectory_paths=list(s.production_trajectory_paths),
                inspected=True, end_time_ps=25000.0, dt_ps=10.0)
        return m

    monkeypatch.setattr(xs, "build_manifest", patched)

    out = tmp_path / "results"
    report = xs.run_study(root, output_dir=out, gmx="gmx")
    assert (out / "run_report.json").is_file()
    assert (out / "run_report.md").is_file()
    assert (out / "source_integrity.json").is_file()
    assert report.integrity["verdict"]["ok"] is True
    assert report.counts()["EXECUTED"] == 5


def test_run_study_refuses_output_inside_dataset(tmp_path, monkeypatch):
    root = tmp_path / "dataset"
    root.mkdir()
    system(root, "AA-A6")
    with pytest.raises((ValueError,)):
        xs.run_study(root, output_dir=root / "results", gmx="gmx")
