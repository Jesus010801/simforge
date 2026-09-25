"""Phase 4 — read-only trajectory diagnostics."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from analysis.campaign.diagnostics import (
    ComponentMotionSeries, DiagnosticContext, DiagnosticParams, all_detectors,
    detector_identity, groups_from_semantic_index, run_diagnostics, write_report,
)
from analysis.campaign.diagnostics import registry as dreg
from analysis.campaign.gmx import run_gmx
from analysis.campaign.models import (
    CampaignWarning, DetectorRunStatus, DiagnosticReport, FrameTimeIndex, Severity,
    ValidationState,
)
from analysis.campaign.trajectory.time_index import validate_timeline
from tests.analysis.campaign.conftest import (
    REAL_TPR, REAL_XTC, requires_gmx, requires_real_traj,
)

BOX10 = [10.0, 10.0, 10.0, 0, 0, 0, 0, 0, 0]


def _index(times, boxes=None, steps=None, truncated=False) -> FrameTimeIndex:
    idx = FrameTimeIndex(trajectory_path="x.xtc", fingerprint=None, backend="fake",
                         backend_version="1", times_ps=[float(t) for t in times],
                         steps=steps if steps is not None else list(range(len(times))),
                         boxes=boxes if boxes is not None else [BOX10] * len(times))
    idx.truncated = truncated
    idx.truncation_evidence = "WARNING: Incomplete frame" if truncated else ""
    idx.backend_reported_frames = len(times)
    return validate_timeline(idx)


def _ctx(tmp_path, ti, motion=None, **kw) -> DiagnosticContext:
    traj = tmp_path / "md.xtc"
    if not traj.exists():
        traj.write_bytes(b"synthetic")
    return DiagnosticContext(trajectory_path=str(traj), time_index=ti,
                             motion_series=motion or {}, **kw)


def _series(role, xs, ys=None, zs=None) -> ComponentMotionSeries:
    n = len(xs)
    pos = np.column_stack([xs, ys if ys is not None else [5.0] * n,
                           zs if zs is not None else [5.0] * n]).astype(float)
    return ComponentMotionSeries(role=role, group=role.capitalize(), convention="center_of_geometry",
                                 positions=pos, group_evidence={"name": role, "atoms_sha256": role})


def _codes(report) -> dict[str, list]:
    out: dict[str, list] = {}
    for d in report.diagnostics:
        out.setdefault(d.code, []).append(d)
    return out


# ═══════════════════════════════════════════════════════════════════════════════
# Registry
# ═══════════════════════════════════════════════════════════════════════════════

def test_registry_levels_and_deferred():
    specs = {s.id: s for s in all_detectors()}
    assert {"timeline_integrity", "duplicate_timestamps", "non_monotonic_time", "time_gaps",
            "box_behaviour", "component_periodic_jump", "partner_receptor_separation"} <= set(specs)
    # Phase 13.5: the molecule-level detector is implemented (no longer reserved)
    assert "molecule_split" not in specs
    mol = specs["molecule_periodic_image_change"]
    assert not mol.deferred and mol.level == 3
    for mid in ("membrane_split_z", "leaflet_discontinuity", "membrane_normal_rotation"):
        assert specs[mid].deferred and specs[mid].level == 4
    assert [s.level for s in all_detectors()] == sorted(s.level for s in all_detectors())


def test_new_detector_needs_no_runner_change(tmp_path, monkeypatch):
    monkeypatch.setattr(dreg, "_DETECTORS", dict(dreg._DETECTORS))
    from analysis.campaign.models import Diagnostic

    @dreg.detector("test_custom", "1", 4, "custom")
    def _custom(ctx):
        return [Diagnostic(code="custom_hit", severity=Severity.INFO, scope="x", message="m",
                           measured={"value": 1})]

    rep = run_diagnostics(_ctx(tmp_path, _index([0, 20, 40])))
    assert "custom_hit" in _codes(rep)
    d = _codes(rep)["custom_hit"][0]
    assert d.detector == {"id": "test_custom", "version": "1"}


# ═══════════════════════════════════════════════════════════════════════════════
# Timeline
# ═══════════════════════════════════════════════════════════════════════════════

def test_regular_timeline_no_diagnostics(tmp_path):
    rep = run_diagnostics(_ctx(tmp_path, _index([0, 20, 40, 60])))
    assert rep.diagnostics == []
    assert rep.n_frames == 4
    ran = {r.detector_id: r for r in rep.detectors}
    assert ran["duplicate_timestamps"].status == DetectorRunStatus.RAN
    assert ran["duplicate_timestamps"].sampling == {
        "strategy": "metadata_only", "frames_examined": 4, "total_frames": 4}


def test_duplicate_timestamps_same_and_distinct_step(tmp_path):
    rep = run_diagnostics(_ctx(tmp_path, _index([0, 20, 40, 40, 60], steps=[0, 10, 20, 20, 30])))
    (d,) = _codes(rep)["duplicate_timestamps"]
    assert d.measured["frames"] == [2, 3] and d.measured["time_ps"] == 40.0
    assert d.measured["md_steps"] == [20, 20]
    assert d.interpretation == "duplicated_state_same_md_step" and d.severity == Severity.WARN

    rep = run_diagnostics(_ctx(tmp_path, _index([0, 20, 40, 40, 60], steps=[0, 10, 20, 21, 30])))
    (d,) = _codes(rep)["duplicate_timestamps"]
    assert d.interpretation == "time_precision_collapse_or_retimed"
    assert d.severity == Severity.REVIEW


def test_non_monotonic_transitions(tmp_path):
    rep = run_diagnostics(_ctx(tmp_path, _index([0, 20, 40, 10, 30])))
    (d,) = _codes(rep)["non_monotonic_time"]
    assert d.measured["n_transitions"] == 1
    assert d.measured["transitions"][0]["frames"] == [2, 3]
    assert d.measured["transitions"][0]["times_ps"] == [40.0, 10.0]
    assert d.severity == Severity.REVIEW and "segment_selection" in d.candidate_remediations


def test_large_gap_reported_not_invalid(tmp_path):
    rep = run_diagnostics(_ctx(tmp_path, _index([0, 20, 40, 60, 200, 220, 240])))
    (g,) = _codes(rep)["time_gap"]
    assert g.measured["typical_dt_ps"] == 20.0
    assert g.measured["gaps"][0] == {"frames": [3, 4], "times_ps": [60.0, 200.0],
                                     "gap_ps": 140.0, "ratio_to_typical": 7.0}
    assert g.severity == Severity.WARN
    assert _codes(rep)["irregular_sampling"][0].severity == Severity.INFO


def test_moderate_irregularity_is_not_a_gap(tmp_path):
    rep = run_diagnostics(_ctx(tmp_path, _index([0, 20, 40, 60, 100, 120])))
    assert "time_gap" not in _codes(rep)                     # 2× typical < 5× threshold
    assert "irregular_sampling" in _codes(rep)


def test_truncated_and_empty(tmp_path):
    rep = run_diagnostics(_ctx(tmp_path, _index([0, 20, 40], truncated=True)))
    (t,) = _codes(rep)["truncated_trajectory"]
    assert t.measured["n_complete_frames"] == 3 and t.severity == Severity.WARN
    empty = _index([])
    rep = run_diagnostics(_ctx(tmp_path, empty))
    (e,) = _codes(rep)["empty_trajectory"]
    assert e.severity == Severity.ERROR


# ═══════════════════════════════════════════════════════════════════════════════
# Box
# ═══════════════════════════════════════════════════════════════════════════════

def _box(a, b, c):
    return [a, b, c, 0, 0, 0, 0, 0, 0]


def test_small_box_fluctuations_not_flagged(tmp_path):
    rng = np.random.default_rng(0)
    boxes = [_box(*(10 + rng.normal(0, 0.01, 3))) for _ in range(50)]
    rep = run_diagnostics(_ctx(tmp_path, _index(range(0, 1000, 20), boxes=boxes)))
    assert "box_discontinuity" not in _codes(rep)


def test_abrupt_volume_jump_flagged(tmp_path):
    rng = np.random.default_rng(1)
    boxes = [_box(*(9.33 + rng.normal(0, 0.005, 3))) for _ in range(20)]
    boxes += [_box(*(9.9 + rng.normal(0, 0.005, 3))) for _ in range(20)]
    rep = run_diagnostics(_ctx(tmp_path, _index(range(0, 800, 20), boxes=boxes)))
    (d,) = _codes(rep)["box_discontinuity"]
    e = d.measured["events"][0]
    assert e["frames"] == [19, 20] and e["times_ps"] == [380.0, 400.0]
    assert e["relative_volume_change"] == pytest.approx(0.195, abs=0.01)
    assert "volume" in e["triggered_by"]
    assert d.severity == Severity.REVIEW and d.interpretation == "abrupt_box_change"


def test_steady_semi_isotropic_deformation_is_drift_not_discontinuity(tmp_path):
    # c grows ~1.5 %/frame, a,b shrink — like the bundled membrane run
    n = 30
    boxes = [_box(12.5 * (0.994 ** i), 12.5 * (0.994 ** i), 20.0 * (1.015 ** i)) for i in range(n)]
    rep = run_diagnostics(_ctx(tmp_path, _index(range(0, 20 * n, 20), boxes=boxes)))
    assert "box_discontinuity" not in _codes(rep)
    (d,) = _codes(rep)["box_shape_drift"]
    assert d.severity == Severity.INFO and d.measured["relative_edge_change"][2] > 0.4


def test_missing_box(tmp_path):
    rep = run_diagnostics(_ctx(tmp_path, _index([0, 20], boxes=[[0.0] * 9, [0.0] * 9])))
    assert _codes(rep)["box_missing"][0].severity == Severity.ERROR


# ═══════════════════════════════════════════════════════════════════════════════
# Component motion
# ═══════════════════════════════════════════════════════════════════════════════

def test_periodic_wrap_evidence(tmp_path):
    ctx = _ctx(tmp_path, _index([0, 20, 40, 60]),
               motion={"ligand": _series("ligand", [1.0, 0.2, 9.9, 9.8])},
               groups={"ligand": "Ligand"})
    rep = run_diagnostics(ctx)
    (d,) = _codes(rep)["component_periodic_wrap"]
    e = d.measured["events"][0]
    assert e["frames"] == [1, 2] and e["raw_magnitude_nm"] == pytest.approx(9.7)
    assert e["min_image_magnitude_nm"] == pytest.approx(0.3)
    assert e["box_translation"] == [1, 0, 0] and e["box_edges_nm"] == [10.0, 10.0, 10.0]
    assert d.interpretation == "likely_periodic_wrap" and d.confidence == "strong"
    assert d.scope == "component:ligand" and d.severity == Severity.WARN
    assert "nojump" in d.candidate_remediations
    assert d.sampling["strategy"] == "all_frames" and d.sampling["frames_examined"] == 4
    assert "component_large_displacement" not in _codes(rep)


def test_genuine_large_motion_not_labelled_wrap(tmp_path):
    ctx = _ctx(tmp_path, _index([0, 20, 40]),
               motion={"peptide": _series("peptide", [1.0, 1.1, 4.6])},
               groups={"peptide": "Peptide"})
    rep = run_diagnostics(ctx)
    assert "component_periodic_wrap" not in _codes(rep)
    (d,) = _codes(rep)["component_large_displacement"]
    assert d.interpretation == "large_motion_uncertain" and d.confidence == "weak"
    assert d.measured["events"][0]["min_image_magnitude_nm"] == pytest.approx(3.5)


def test_system_translation(tmp_path):
    ctx = _ctx(tmp_path, _index([0, 20, 40]),
               motion={"system": _series("system", [5.0, 9.8, 0.1])},
               groups={"system": "System"})
    rep = run_diagnostics(ctx)
    # 5.0→9.8 is 4.8 nm (0.48 box) with min-image 4.8: uncertain; 9.8→0.1 is a wrap
    (w,) = _codes(rep)["system_periodic_translation"]
    assert w.scope == "system" and w.candidate_remediations == ["center", "nojump"]


# ═══════════════════════════════════════════════════════════════════════════════
# Partner ↔ receptor
# ═══════════════════════════════════════════════════════════════════════════════

def _pair_ctx(tmp_path, lig_x):
    n = len(lig_x)
    return _ctx(tmp_path, _index(range(0, 20 * n, 20)),
                motion={"receptor": _series("receptor", [5.0] * n),
                        "ligand": _series("ligand", lig_x)},
                groups={"receptor": "Receptor", "ligand": "Ligand"})


def test_ligand_periodic_separation(tmp_path):
    rep = run_diagnostics(_pair_ctx(tmp_path, [5.8, 5.8, 14.1, 14.1]))
    (d,) = _codes(rep)["ligand_receptor_periodic_separation"]
    e = d.measured["events"][0]
    assert e["frames"] == [1, 2]
    assert e["raw_distance_nm"] == [0.8, 9.1]
    assert e["min_image_distance_nm"] == [0.8, 0.9]
    assert d.confidence == "strong" and d.severity == Severity.WARN
    assert d.candidate_remediations == ["complex_reconstruction", "minimum_image_analysis"]
    (disc,) = _codes(rep)["ligand_receptor_image_discrepancy"]
    assert disc.measured["frame_ranges"] == [[2, 3]]
    assert "ligand_receptor_separation_change" not in _codes(rep)


def test_ligand_genuine_separation_not_periodic(tmp_path):
    rep = run_diagnostics(_pair_ctx(tmp_path, [5.8, 5.8, 8.8, 8.9]))
    assert "ligand_receptor_periodic_separation" not in _codes(rep)
    (d,) = _codes(rep)["ligand_receptor_separation_change"]
    e = d.measured["events"][0]
    assert e["raw_distance_nm"] == [0.8, 3.8] and e["min_image_distance_nm"] == [0.8, 3.8]
    assert d.interpretation == "possible_physical_separation" and d.severity == Severity.INFO
    assert d.candidate_remediations == []


def test_partner_detector_not_applicable_without_partner(tmp_path):
    ctx = _ctx(tmp_path, _index([0, 20]), motion={"receptor": _series("receptor", [5, 5])},
               groups={"receptor": "Receptor"})
    rep = run_diagnostics(ctx)
    run = next(r for r in rep.detectors if r.detector_id == "partner_receptor_separation")
    assert run.status == DetectorRunStatus.NOT_APPLICABLE and "ligand" in run.reason


# ═══════════════════════════════════════════════════════════════════════════════
# Provenance / identity
# ═══════════════════════════════════════════════════════════════════════════════

_NDX = "[ System ]\n1 2 3 4 5 6\n[ Receptor ]\n1 2 3\n[ Ligand ]\n4\n[ Unused ]\n6\n"


def _id_ctx(tmp_path, ndx_text):
    ndx = tmp_path / "semantic.ndx"
    ndx.write_text(ndx_text)
    ctx = _ctx(tmp_path, _index([0, 20]), index_path=str(ndx),
               groups={"system": "System", "receptor": "Receptor", "ligand": "Ligand"})
    return ctx


def test_detector_identity_tracks_used_atoms_only(tmp_path):
    spec = dreg.get("partner_receptor_separation")
    base = detector_identity(spec, _id_ctx(tmp_path, _NDX))
    assert detector_identity(spec, _id_ctx(tmp_path, _NDX)) == base
    unused = detector_identity(spec, _id_ctx(tmp_path, _NDX.replace("[ Unused ]\n6", "[ Unused ]\n5 6")))
    assert unused == base
    lig = detector_identity(spec, _id_ctx(tmp_path, _NDX.replace("[ Ligand ]\n4", "[ Ligand ]\n4 5")))
    assert lig != base
    ctx = _id_ctx(tmp_path, _NDX)
    ctx.params = DiagnosticParams(jump_box_fraction=0.4)
    assert detector_identity(spec, ctx) != base


def test_motion_probe_key_tracks_group_atoms(tmp_path):
    k1 = _id_ctx(tmp_path, _NDX)._motion_key("ligand", "center_of_geometry")
    k2 = _id_ctx(tmp_path, _NDX.replace("[ Unused ]\n6", "[ Unused ]\n5 6"))._motion_key(
        "ligand", "center_of_geometry")
    k3 = _id_ctx(tmp_path, _NDX.replace("[ Ligand ]\n4", "[ Ligand ]\n4 5"))._motion_key(
        "ligand", "center_of_geometry")
    assert k1 == k2 != k3


def test_report_roundtrip_and_write(tmp_path):
    rep = run_diagnostics(_pair_ctx(tmp_path, [5.8, 5.8, 14.1, 14.1]))
    path = write_report(rep, tmp_path / "out")
    assert path.name == "trajectory_diagnostics.json"
    d = json.loads(path.read_text())
    assert d["schema_version"] == "simforge/trajectory-diagnostics/v1"
    back = DiagnosticReport.from_dict(d)
    assert back.to_dict() == rep.to_dict()
    diag = d["diagnostics"][0]
    assert diag["provenance"]["detector_identity"] and "params" in diag["provenance"]
    assert d["summary"]["by_severity"]["warn"] >= 1


def test_groups_from_semantic_index(tmp_path):
    ndx = tmp_path / "i.ndx"
    ndx.write_text("[ System ]\n1 2\n[ Receptor ]\n1\n[ Membrane ]\n2\n[ Other ]\n2\n")
    assert groups_from_semantic_index(str(ndx)) == {
        "system": "System", "receptor": "Receptor", "membrane": "Membrane"}


# ═══════════════════════════════════════════════════════════════════════════════
# Real GROMACS probes — read-only
# ═══════════════════════════════════════════════════════════════════════════════

def _sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _synthetic_wrap_traj(tmp_path: Path) -> tuple[Path, Path, Path]:
    """5 atoms: receptor = 1..3 near x≈9.5, ligand = 4 crossing x=10 → 0, one water."""
    lig_x = [9.80, 9.90, 0.20, 0.25, 0.30]
    lines = []
    for i, lx in enumerate(lig_x):
        lines += [f"wrap t= {20.0 * i:.5f} step= {i * 10}", "5"]
        for k, (rn, an, x) in enumerate([("REC", "CA", 9.4), ("REC", "CB", 9.5),
                                         ("REC", "CC", 9.6), ("LIG", "C1", lx),
                                         ("SOL", "OW", 3.0)]):
            lines.append(f"{1 if rn == 'REC' else (2 if rn == 'LIG' else 3):>5}{rn:<5}{an:>5}"
                         f"{k + 1:>5}{x:8.3f}{5.0:8.3f}{5.0:8.3f}")
        lines.append("  10.00000  10.00000  10.00000")
    src = tmp_path / "src"
    src.mkdir()
    gro_in = src / "frames.gro"
    gro_in.write_text("\n".join(lines) + "\n")
    xtc = src / "md.xtc"
    res = run_gmx(["trjconv", "-f", str(gro_in), "-s", str(gro_in), "-o", str(xtc)], stdin="0\n")
    assert res.ok, res.stderr[-400:]
    ref = src / "ref.gro"
    ref.write_text("\n".join(lines[:8]) + "\n")
    ndx = src / "semantic.ndx"
    ndx.write_text("[ Receptor ]\n1 2 3\n[ Ligand ]\n4\n[ System ]\n1 2 3 4 5\n")
    return xtc, ref, ndx


@requires_gmx
def test_gromacs_probe_read_only(tmp_path, monkeypatch):
    xtc, ref, ndx = _synthetic_wrap_traj(tmp_path)
    before = {p: _sha(p) for p in (xtc, ref, ndx)}
    listing = sorted(p.name for p in xtc.parent.iterdir())

    import analysis.campaign.diagnostics.context as dctx
    import analysis.campaign.trajectory.time_index as tix
    import analysis.campaign.trajectory.preprocessor as pre
    calls = []

    def spy(fn):
        def wrapped(args, **kw):
            calls.append(list(args))
            return fn(args, **kw)
        return wrapped
    monkeypatch.setattr(dctx, "run_gmx", spy(dctx.run_gmx))
    monkeypatch.setattr(tix, "run_gmx", spy(tix.run_gmx))

    def _no_views(**kw):
        raise AssertionError("diagnostics must not build views")
    monkeypatch.setattr(pre, "build_view", _no_views)

    ctx = DiagnosticContext(trajectory_path=str(xtc), structure_path=str(ref),
                            index_path=str(ndx),
                            groups=groups_from_semantic_index(str(ndx)),
                            cache_dir=str(tmp_path / "diag_cache"))
    rep = run_diagnostics(ctx)
    codes = _codes(rep)

    assert rep.n_frames == 5
    wrap = codes["component_periodic_wrap"][0]
    assert wrap.scope == "component:ligand"
    assert wrap.provenance["centre_convention"] == "center_of_geometry"   # no .tpr → COG, said so
    sep = codes["ligand_receptor_periodic_separation"][0]
    assert sep.measured["events"][0]["frames"] == [1, 2]
    assert "ligand_receptor_image_discrepancy" in codes

    # read-only: sources untouched, nothing written beside them, no corrective commands
    assert {p: _sha(p) for p in (xtc, ref, ndx)} == before
    assert sorted(p.name for p in xtc.parent.iterdir()) == listing
    for a in calls:
        assert a[0] in ("trjconv", "trajectory")
        assert not {"-pbc", "-center", "-fit", "-ur", "-cat"} & set(a)
        if a[0] == "trjconv":                                 # Phase 1 time-index probe only
            assert a[a.index("-o") + 1].endswith("frames.gro")
    n_calls = len(calls)

    # probe cache: a second run re-uses evidence (no new gmx trajectory pass)
    ctx2 = DiagnosticContext(trajectory_path=str(xtc), structure_path=str(ref),
                             index_path=str(ndx), groups=groups_from_semantic_index(str(ndx)),
                             cache_dir=str(tmp_path / "diag_cache"))
    rep2 = run_diagnostics(ctx2)
    assert len(calls) == n_calls
    assert [d.code for d in rep2.diagnostics] == [d.code for d in rep.diagnostics]


@requires_gmx
@requires_real_traj
def test_real_trajectory_diagnostics(tmp_path):
    from analysis.campaign.structure.component_detector import detect_components
    from analysis.campaign.structure.index_builder import build_semantic_index
    from tests.analysis.campaign.conftest import REAL_GRO

    before = _sha(REAL_XTC)
    det = detect_components(structure_path=REAL_GRO)
    idx = build_semantic_index(structure_path=REAL_GRO, components=det.components,
                               out_ndx=tmp_path / "semantic.ndx")
    groups = groups_from_semantic_index(idx.path, det.components)
    assert {"system", "receptor", "membrane"} <= set(groups)
    ctx = DiagnosticContext(trajectory_path=str(REAL_XTC), structure_path=str(REAL_TPR),
                            index_path=idx.path, groups=groups,
                            cache_dir=str(tmp_path / "cache"))
    rep = run_diagnostics(ctx)
    assert rep.n_frames == 51
    status = {r.detector_id: r.status for r in rep.detectors}
    assert status["component_periodic_jump"] == DetectorRunStatus.RAN
    assert status["partner_receptor_separation"] == DetectorRunStatus.NOT_APPLICABLE
    assert status["membrane_split_z"] == DetectorRunStatus.DEFERRED
    assert status["molecule_periodic_image_change"] == DetectorRunStatus.NOT_APPLICABLE
    motion = [p for p in rep.probes if p["probe"] == "component_motion"]
    assert motion and all(p["convention"] == "center_of_mass" for p in motion)
    assert {d.code for d in rep.diagnostics} <= {"box_shape_drift"}   # observed result
    assert _sha(REAL_XTC) == before
