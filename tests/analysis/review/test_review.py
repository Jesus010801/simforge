"""Phase 9 — ReviewDataset / review-session preparation (headless)."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from analysis.campaign.models import (
    AnalysisResult, AnalysisStatus, AnnotationRecord, Axis, AxisKind, ClassificationState,
    FrameAlignment, MolecularComponent, ResultArray, SourceFingerprint, StorageRef,
    SystemRecord,
)
from analysis.review import (
    Availability, CapabilityState, DisplayRequest, ObservableEntry, ObservableRequest,
    RequestError, ReviewDataset, ReviewRequest, ReviewResult, load_request_file,
    open_review_session, parse_show, prepare_review,
)
from analysis.review.prepare import _classify, session_identity_evidence
from analysis.review.store import put_file, verify
from analysis.review.sync import ReviewTimeline, SyncMapping, sync_metadata
from tests.analysis.campaign.conftest import (
    REAL_GRO, REAL_TPR, REAL_XTC, REPO_ROOT, requires_gmx, requires_real_traj,
)

# ═══════════════════════════════════════════════════════════════════════════════
# Requests
# ═══════════════════════════════════════════════════════════════════════════════


def test_observable_instances_parse_canonically():
    reqs = parse_show(["rmsd-receptor,rg(selection=component:receptor)",
                       "com-distance(selection_b=annotation:site, selection_a=component:ligand)",
                       "rg(selection=component:receptor)"])
    assert [r.instance_id for r in reqs] == [
        "rmsd-receptor", "rg(selection=component:receptor)",
        "com-distance(selection_a=component:ligand,selection_b=annotation:site)"]
    assert reqs[2].params == {"selection_a": "component:ligand", "selection_b": "annotation:site"}
    # same observable twice with different parameters = two instances
    two = parse_show(["rg(selection=component:receptor),rg(selection=annotation:sub)"])
    assert len(two) == 2 and two[0].observable == two[1].observable == "rg"
    for bad in ("rg(selection)", "rg(selection=a", "rg)", "rg(a=1,a=2)"):
        with pytest.raises(RequestError):
            parse_show([bad])


def test_request_file(tmp_path):
    f = tmp_path / "req.yaml"
    f.write_text("observables:\n  - observable: rg\n    parameters: {selection: 'annotation:x'}\n"
                 "  - rmsd-receptor\n")
    assert [r.instance_id for r in load_request_file(f)] == [
        "rg(selection=annotation:x)", "rmsd-receptor"]


def test_display_request_is_typed_and_intent_is_scoped():
    assert DisplayRequest.parse(None).requirements().view_kind() == "raw"
    assert DisplayRequest.parse("whole").requirements().view_kind() == "whole"
    c = DisplayRequest.parse("center:Receptor", intent=True)
    assert c.requirements().centering_target == "Receptor"
    intent = c.policy_intent()
    assert intent.source == "user_flag" and set(intent.operations) == {"make_whole", "center"}
    assert DisplayRequest.parse("fit:Receptor").policy_intent().source == "auto"
    for bad in ("center", "raw:X", "wobble"):
        with pytest.raises(RequestError):
            DisplayRequest.parse(bad)


# ═══════════════════════════════════════════════════════════════════════════════
# Session identity
# ═══════════════════════════════════════════════════════════════════════════════

def _rec(**over) -> SystemRecord:
    fp = lambda d: SourceFingerprint("/x", 1, 1, "strong", d)  # noqa: E731
    rec = SystemRecord("s", "c", "r", source_fingerprints={
        "trajectory": fp("traj-1"), "topology": fp("tpr-1")})
    rec.components = [MolecularComponent("receptor", "receptor",
                                         classification_state=ClassificationState.RESOLVED,
                                         chain_ids=["A"], atom_count=100)]
    rec.annotations = [AnnotationRecord("site", "residue_set", "binding_site", "user_yaml",
                                        "active", identity="id-1", atoms_sha256="atoms-1",
                                        description="old words")]
    for k, v in over.items():
        setattr(rec, k, v)
    return rec


def _req(*show, display="raw", intent=False):
    return ReviewRequest(parse_show(list(show) or ["rg(selection=annotation:site)"]),
                         DisplayRequest.parse(display, intent=intent))


def _sid(rec, req):
    from analysis.campaign.results import definition_token
    return definition_token(session_identity_evidence(rec, req))


def test_session_identity_scientific_inputs_only():
    base = _sid(_rec(), _req())
    assert base == _sid(_rec(), _req())                              # deterministic
    rec = _rec()
    rec.annotations[0].description = "new words"                    # not scientific
    rec.trajectory_paths.append("/elsewhere/md.xtc")                 # a path, not content
    assert _sid(rec, _req()) == base
    rec = _rec()
    rec.annotations[0].atoms_sha256, rec.annotations[0].identity = "atoms-2", "id-2"
    assert _sid(rec, _req()) != base                                 # annotation atoms
    assert _sid(_rec(), _req("rg(selection=component:receptor)")) != base      # observable
    assert _sid(_rec(), _req(display="center:Receptor")) != \
        _sid(_rec(), _req(display="center:Receptor", intent=True))  # display intent
    rec = _rec()
    rec.source_fingerprints["trajectory"] = SourceFingerprint("/x", 1, 1, "strong", "traj-2")
    assert _sid(rec, _req()) != base                                 # source content
    rec = _rec()
    rec.components[0].atom_count = 101
    assert _sid(rec, _req()) != base                                 # component definition


# ═══════════════════════════════════════════════════════════════════════════════
# Sync metadata (generic, by axes — never by observable name)
# ═══════════════════════════════════════════════════════════════════════════════

TIMES = [0.0, 20.0, 40.0, 60.0, 80.0]


def _xvg(path, rows):
    path.write_text("".join(" ".join(f"{v:.3f}" for v in r) + "\n" for r in rows))
    return path


def _time_array(tmp_path, name, times, axes_extra=(), quantity="q"):
    f = _xvg(tmp_path / f"{name}.xvg", [(t, 1.0) for t in times])
    axes = [Axis("time", AxisKind.TIME, unit="ps", values_ref=StorageRef(str(f), "xvg", column=0),
                 alignment=FrameAlignment("one_row_per_frame", trajectory_digest="traj"))]
    return ResultArray(name, quantity, "nm", axes=axes + list(axes_extra),
                       storage=StorageRef(str(f), "xvg", column=1))


def test_sync_metadata_by_axis_signature(tmp_path):
    rmsd = sync_metadata(_time_array(tmp_path, "rmsd", TIMES), TIMES, "traj")
    assert rmsd["syncable"] and rmsd["mapping"] == SyncMapping.ROW_IS_FRAME
    assert rmsd["mapping_policy"]["interpolation"] == "none"
    sub = sync_metadata(_time_array(tmp_path, "rg", TIMES[::2]), TIMES, "traj")
    assert sub["syncable"] and sub["mapping"] == SyncMapping.EXACT_TIMESTAMP   # frame != row
    rmsf = ResultArray("rmsf", "rmsf", "nm", axes=[Axis("residue", AxisKind.CATEGORY)])
    assert not sync_metadata(rmsf, TIMES, "traj")["syncable"]
    pore = _time_array(tmp_path, "pore", TIMES, axes_extra=[Axis("z", AxisKind.COORDINATE, "nm")])
    p = sync_metadata(pore, TIMES, "traj")
    assert p["syncable"] and p["time_axis_position"] == 0
    contact = ResultArray("map", "contact", None, axes=[Axis("residue_i", AxisKind.CATEGORY),
                                                        Axis("residue_j", AxisKind.CATEGORY)])
    c = sync_metadata(contact, TIMES, "traj")
    assert not c["syncable"] and "no time axis" in c["reason"]
    off = sync_metadata(_time_array(tmp_path, "off", [0.0, 21.0]), TIMES, "traj")
    assert not off["syncable"] and "match no frame" in off["reason"]
    declared = ResultArray("d", "q", axes=[Axis("time", AxisKind.TIME, unit="ps",
                                                alignment=FrameAlignment("one_row_per_frame",
                                                                         trajectory_digest="traj"))])
    assert sync_metadata(declared, TIMES, "traj")["syncable"]
    assert not sync_metadata(declared, TIMES, "other")["syncable"]


def test_review_timeline_mapping_with_duplicates(tmp_path):
    times = [0.0, 10.0, 10.0, 20.0]
    tl = ReviewTimeline(times)
    assert tl.display_frame_to_time(2) == 10.0
    look = tl.time_to_display_frame(10.0)
    assert look.frame == 1 and look.duplicate_frames == [1, 2]             # exposed, not hidden
    assert tl.time_to_display_frame(15.0).frame is None                     # exact by default
    cur = tl.time_to_display_frame(15.0, cursor=True)
    assert cur.frame == 1 and cur.tie                                       # earlier on tie
    assert tl.time_to_display_frame(99.0, cursor=True).frame is None        # no extrapolation
    arr = _time_array(tmp_path, "dup", times)
    res = ReviewResult(arr, sync=sync_metadata(arr, times, "traj"))
    assert res.sync["mapping"] == SyncMapping.ROW_IS_FRAME
    assert tl.sample_to_frames(res, 2) == [2]                               # row i is frame i
    sub = _time_array(tmp_path, "sub", [0.0, 20.0])
    r2 = ReviewResult(sub, sync=sync_metadata(sub, times, "traj"))
    assert tl.sample_to_frames(r2, 1) == [3]


# ═══════════════════════════════════════════════════════════════════════════════
# States and dataset round trip
# ═══════════════════════════════════════════════════════════════════════════════

def test_capability_states_reuse_existing_vocabulary():
    r = lambda st, msg="", **kw: AnalysisResult("x", "s", st, message=msg, **kw)  # noqa: E731
    assert _classify(r("success"), {}) == ("available", "computed")
    assert _classify(r("cached"), {}) == ("available", "cached")
    assert _classify(r("cached", reused_from={"user_selected": True}), {}) == \
        ("available", "explicit_import")
    assert _classify(r("skipped"), {})[0] == "not_applicable"
    none = {"reasons": ["selection (component:ligand): component 'ligand': 0 components of this type"]}
    assert _classify(r("review_required"), none)[0] == "not_applicable"       # nothing to measure
    assert _classify(r("review_required", "component 'ligand' is ambiguous, not RESOLVED"),
                     {"reasons": ["x"]})[0] == "ambiguous"
    unresolved = {"reasons": ["annotation 'site' is unresolved"]}
    assert _classify(r("review_required"), unresolved)[0] == "review_required"  # applicable, blocked
    assert _classify(r("failed"), {})[0] == "failed"
    assert _classify(r("error"), {})[0] == "unsupported"


def test_dataset_round_trip(tmp_path):
    arr = _time_array(tmp_path, "rg", TIMES)
    arr.view_ref = "raw-view"
    entries = [
        ObservableEntry("rg(selection=component:receptor)", "rg", {"selection": "component:receptor"},
                        availability=Availability.COMPUTED,
                        results=[ReviewResult(arr, sync=sync_metadata(arr, TIMES, "traj"))],
                        view_ref="raw-view"),
        ObservableEntry("rmsd-receptor", "rmsd-receptor", availability=Availability.CACHED,
                        reuse={"decision": "reused", "reused_from": {"provenance": "/p"}},
                        view_ref="whole-view"),
        ObservableEntry("rg(selection=annotation:tm_6)", "rg", {"selection": "annotation:tm_6"},
                        state=CapabilityState.REVIEW_REQUIRED, reason="tm_6 unresolved",
                        annotations_used=["tm_6"], blocking_annotations=["tm_6"]),
        ObservableEntry("rg(selection=annotation:x)", "rg", availability=Availability.EXPLICIT_IMPORT,
                        externally_supplied=True, independently_verified=False),
    ]
    ds = ReviewDataset(
        session_id="abc", identity_evidence={"schema": "x"},
        system={"system_id": "s", "membrane_present": True},
        display_view={"role": "display", "view_ref": "raw-view", "kind": "raw",
                      "purpose": "display", "status": "available", "decisions": []},
        analysis_views=[{"role": "analysis", "view_ref": "raw-view", "kind": "raw"},
                        {"role": "analysis", "view_ref": "whole-view", "kind": "whole",
                         "decisions": [{"operation": "make_whole", "classification": "valid"}]}],
        timeline={"n_frames": 5, "times_ref": {"path": "t.npy", "sha256": "0"}},
        annotations=[{"annotation_id": "tm_6", "state": "unresolved",
                      "blocks": ["rg(selection=annotation:tm_6)"]}],
        diagnostics={"execution": "auto", "status": "complete", "report_identity": "r"},
        observables=entries)
    back = ReviewDataset.from_dict(json.loads(json.dumps(ds.to_dict())))
    assert back.to_dict() == ds.to_dict()
    assert back.display_view["view_ref"] == "raw-view" and len(back.analysis_views) == 2
    assert [o.instance_id for o in back.blocked()] == ["rg(selection=annotation:tm_6)"]
    imp = back.observable("rg(selection=annotation:x)")
    assert imp.externally_supplied and not imp.independently_verified
    assert back.observables[0].results[0].sync["syncable"]
    s = back.summary()["observables"]
    assert (s["computed"], s["cached"], s["explicit_import"], s["blocked"]) == (1, 1, 1, 1)


def test_store_is_write_once_and_content_addressed(tmp_path):
    f = tmp_path / "r.xvg"
    f.write_text("0 1\n")
    p1, h1 = put_file(tmp_path / "store", f)
    f.write_text("0 2\n")                           # later "recomputation" of the source
    p2, h2 = put_file(tmp_path / "store", f)
    assert p1 != p2 and p1.read_text() == "0 1\n" and verify(p1, h1) and verify(p2, h2)


# ═══════════════════════════════════════════════════════════════════════════════
# Real GROMACS: bundled GLP-1R membrane system (build-spec annotations)
# ═══════════════════════════════════════════════════════════════════════════════

SHOW = ["rg(selection=component:receptor)", "rg(selection=annotation:tm_1)",
        "rg(selection=annotation:tm_6)", "rmsd-receptor"]


def _run_dir(root: Path) -> Path:
    run = root / "run"
    (run / "metadata").mkdir(parents=True)
    spec = json.loads((REAL_GRO.parents[2] / "metadata" / "run_info.json").read_text())["yaml_source"]
    (run / "metadata" / "run_info.json").write_text(json.dumps({"yaml_source": spec}))
    step = run / "steps" / "11_production_md"
    step.mkdir(parents=True)
    for f in (REAL_XTC, REAL_TPR, REAL_GRO):
        (step / f.name).symlink_to(f)
    return run


def _prepare(run, show=SHOW, **kw):
    display = kw.pop("display", "raw")
    intent = kw.pop("intent", False)
    return prepare_review(run, ReviewRequest(parse_show(show), DisplayRequest.parse(
        display, intent=intent)), **kw)


def _guard(monkeypatch):
    import analysis.campaign.observables.generic as gen
    import analysis.campaign.observables.rmsd as rmsd

    def boom(*a, **k):
        raise AssertionError("no observable GROMACS execution expected")
    monkeypatch.setattr(gen, "run_gmx", boom)
    monkeypatch.setattr(rmsd, "run_gmx", boom)


@pytest.fixture(scope="module")
def glp(tmp_path_factory):
    if not (REAL_XTC.is_file() and __import__("shutil").which("gmx")):
        pytest.skip("needs gmx and the bundled trajectory")
    from analysis.campaign.orchestration.study_analyzer import clear_view_cache
    clear_view_cache()
    run = _run_dir(tmp_path_factory.mktemp("glp"))
    return run, _prepare(run)


@requires_gmx
@requires_real_traj
def test_membrane_session_prepared(glp):
    run, prep = glp
    ds = prep.dataset
    assert prep.dataset_path.is_file() and ds.status == "prepared"
    by = {o.instance_id: o for o in ds.observables}
    for iid in ("rg(selection=component:receptor)", "rg(selection=annotation:tm_1)", "rmsd-receptor"):
        o = by[iid]
        assert o.state == "available" and o.availability == "computed", (iid, o.reason)
        (r,) = o.results
        assert r.sync["syncable"] and r.sync["mapping"] == "row_i_is_frame_i"
        assert r.array.storage.path.startswith("../../store/")            # frozen, relative
    blocked = by["rg(selection=annotation:tm_6)"]
    assert blocked.state == "review_required" and "tm_6" in blocked.reason
    assert blocked.blocking_annotations == ["tm_6"] and not blocked.results
    ann = {a["annotation_id"]: a for a in ds.annotations}
    assert ann["tm_6"]["state"] == "unresolved"
    assert ann["tm_6"]["blocks"] == ["rg(selection=annotation:tm_6)"]
    assert ann["tm_1"]["state"] == "active" and ann["tm_1"]["used_by"]
    # membrane-safe display: the source coordinates, planned for purpose DISPLAY
    dv = ds.display_view
    assert (dv["kind"], dv["purpose"], dv["status"], dv["decisions"]) == ("raw", "display",
                                                                        "available", [])
    # analysis views are separate records: RMSD's whole view never became the display
    kinds = {v["view_ref"]: v for v in ds.analysis_views}
    rmsd_view = kinds[by["rmsd-receptor"].view_ref]
    assert rmsd_view["kind"] == "whole" and rmsd_view["view_ref"] != dv["view_ref"]
    assert rmsd_view["purpose"] == "intramolecular_shape"
    tl = ds.timeline
    assert (tl["n_frames"], tl["start_time_ps"], tl["end_time_ps"]) == (51, 0.0, 1000.0)
    assert tl["display"]["relation"] == "identical_to_source"
    assert ds.diagnostics["execution"] == "auto" and "ran_found_nothing" in ds.diagnostics["detectors"]
    _, val = open_review_session(prep.dataset_path)
    assert val.valid, val.problems
    assert not list((prep.session_dir.parent).glob(".building-*"))


@requires_gmx
@requires_real_traj
def test_identical_request_reuses_session_and_results(glp, monkeypatch):
    run, prep = glp
    _guard(monkeypatch)
    again = _prepare(run)
    assert again.reused_session and again.dataset.session_id == prep.dataset.session_id
    forced = _prepare(run, force=True)                   # re-prepare: Phase 8 reuse, no gmx
    assert not forced.reused_session and forced.dataset.session_id == prep.dataset.session_id
    for o in forced.dataset.available():
        assert o.availability == "cached" and o.reuse["decision"] == "reused"
    # reopening needs no discovery / analysis at all
    ds, val = open_review_session(forced.dataset_path)
    assert val.valid and ds.session_id == prep.dataset.session_id


@requires_gmx
@requires_real_traj
def test_frozen_results_survive_recomputation_and_corruption_is_detected(glp):
    run, prep = glp
    _prepare(run, force=True)
    o = prep.dataset.observable("rg(selection=component:receptor)")
    origin = Path(o.results[0].origin["storage"]["path"])
    origin.write_text(origin.read_text() + "\n")          # campaign output rewritten later
    ds, val = open_review_session(prep.dataset_path)
    assert val.valid                                      # the session references frozen copies
    frozen = prep.session_dir / ds.observable(o.instance_id).results[0].array.storage.path
    frozen.chmod(0o644)
    frozen.write_text(frozen.read_text().replace("0.000", "9.000", 1))
    ds, val = open_review_session(prep.dataset_path)
    assert not val.valid and o.instance_id in val.unavailable
    assert "changed after the session was prepared" in val.unavailable[o.instance_id]
    again = _prepare(run)                                # stale session is never trusted
    assert not again.reused_session and again.replaced
    assert open_review_session(again.dataset_path)[1].valid
    # the campaign file was modified too: the result was recomputed, not reused
    rec = again.dataset.observable(o.instance_id)
    assert rec.availability == "computed" and rec.reuse["decision"] == "recomputed"


@requires_gmx
@requires_real_traj
def test_session_isolation_and_unrelated_changes(tmp_path):
    run = _run_dir(tmp_path)
    a = _prepare(run, show=["rg(selection=component:receptor)"])
    b = _prepare(run, show=["rg(selection=annotation:tm_1)"])
    c = _prepare(run, show=["rg(selection=component:receptor)"], display="center:Receptor",
                 intent=True)
    ids = {a.dataset.session_id, b.dataset.session_id, c.dataset.session_id}
    assert len(ids) == 3 and len({p.session_dir for p in (a, b, c)}) == 3
    assert open_review_session(a.dataset_path)[1].valid          # not overwritten by b / c
    # explicit display intent: a centred view, distinct from the raw analysis view
    dv = c.dataset.display_view
    assert dv["status"] == "available" and dv["kind"] == "centered"
    assert {d["intent_source"] for d in dv["decisions"]} == {"user_flag"}
    assert dv["view_ref"] not in {v["view_ref"] for v in c.dataset.analysis_views}
    assert c.dataset.timeline["display"]["relation"] == "identical_to_source"
    # without intent the same request is REVIEW_REQUIRED, never silently chosen
    d = _prepare(run, show=["rg(selection=component:receptor)"], display="center:Receptor")
    assert d.dataset.display_view["status"] == "review_required"
    assert "--display-intent" in d.dataset.display_view["reason"]
    # a different output directory does not change the scientific identity
    e = _prepare(run, show=["rg(selection=component:receptor)"], output=tmp_path / "elsewhere")
    assert e.dataset.session_id == a.dataset.session_id and e.session_dir != a.session_dir


@requires_gmx
@requires_real_traj
def test_annotation_change_invalidates_only_dependent_results(tmp_path, monkeypatch):
    run = _run_dir(tmp_path)
    y = run / "simforge_annotations.yaml"
    decl = ("schema_version: simforge/annotations/v1\nstructural_annotation:\n  residue_sets:\n"
            "  - id: site\n    kind: binding_site\n    description: {desc}\n"
            "    selection: {{numbering: topology, residues: \"{res}\", polymer_only: true}}\n")
    y.write_text(decl.format(desc="first", res="177-190"))
    show = ["rg(selection=component:receptor)", "rg(selection=annotation:site)"]
    first = _prepare(run, show=show)
    y.write_text(decl.format(desc="reworded", res="177-190"))          # description only
    assert _prepare(run, show=show).reused_session
    y.write_text(decl.format(desc="reworded", res="177-191"))          # atom set changed
    second = _prepare(run, show=show)
    assert second.dataset.session_id != first.dataset.session_id
    by = {o.instance_id: o for o in second.dataset.observables}
    assert by["rg(selection=component:receptor)"].availability == "cached"
    site = by["rg(selection=annotation:site)"]
    assert site.availability == "computed" and site.reuse["decision"] == "recomputed"
    assert "atoms_sha256" in site.reuse["evaluated"][0]["reason"]
    assert open_review_session(first.dataset_path)[1].valid           # older session intact


@requires_gmx
@requires_real_traj
def test_failed_preparation_leaves_no_session(tmp_path, monkeypatch):
    import analysis.review.validate as v
    from analysis.review.prepare import ReviewError
    run = _run_dir(tmp_path)
    monkeypatch.setattr(v, "validate_review_dataset",
                        lambda ds, d: v.ReviewValidation(False, ["injected"]))
    with pytest.raises(ReviewError, match="injected"):
        _prepare(run, show=["rmsd-receptor"])
    sessions = run / "simforge_analysis" / "review" / "sessions"
    assert not list(sessions.iterdir())                               # nothing half-written


# ═══════════════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════════════

def _cli(*args, timeout=600):
    return subprocess.run([sys.executable, "-m", "cli", *args], cwd=REPO_ROOT,
                          capture_output=True, text=True, timeout=timeout)


@requires_gmx
@requires_real_traj
def test_cli_dry_run_strict_json_and_trajectory_mode(tmp_path):
    run = _run_dir(tmp_path)
    show = ["--show", "rg(selection=component:receptor),rg(selection=annotation:tm_6)"]
    dry = _cli("trajectory", "review", str(run), *show, "--dry-run")
    assert dry.returncode == 0, dry.stderr
    assert "dry run" in dry.stdout and not (run / "simforge_analysis").exists()
    strict = _cli("trajectory", "review", str(run), *show, "--strict")
    assert strict.returncode == 1 and "review_required" in strict.stdout
    js = _cli("trajectory", "review", str(run), *show, "--json")
    assert js.returncode == 0
    out = json.loads(js.stdout)
    ds = ReviewDataset.from_dict(out["review_dataset"])              # same schema, no variant
    assert out["preparation"]["reused_session"] and ds.schema_version.startswith("simforge/review")
    # TRAJECTORY --topology mode: explicit files, own workspace, same machinery
    traj = _cli("trajectory", "review", str(REAL_XTC), "--topology", str(REAL_TPR),
                "--structure", str(REAL_GRO), "--show", "rg(selection=component:receptor)",
                "--output", str(tmp_path / "t"), "--json")
    assert traj.returncode == 0, traj.stderr
    tds = json.loads(traj.stdout)["review_dataset"]
    assert tds["sources"]["inputs"]["mode"] == "trajectory"
    assert tds["observables"][0]["state"] == "available"
    lst = _cli("trajectory", "review-observables", str(run), "--json")
    assert lst.returncode == 0
    listing = json.loads(lst.stdout)
    refs = {s["ref"]: s for s in listing["selections"]}
    assert refs["annotation:tm_1"]["usable"] and not refs["annotation:tm_6"]["usable"]
    assert {o["id"] for o in listing["observables"]} >= {"rg", "rmsd-receptor"}
