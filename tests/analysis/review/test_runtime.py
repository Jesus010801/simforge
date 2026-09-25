"""Phase 10 — sync hub, viewer boundary, local review server (presentation only)."""
from __future__ import annotations

import json
import os
import signal
import socket
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.request
from dataclasses import replace
from pathlib import Path

import pytest

from analysis.campaign.fingerprint import fingerprint_file
from analysis.campaign.models import Axis, AxisKind, FrameAlignment, ResultArray, StorageRef
from analysis.review.dataset import ObservableEntry, ReviewDataset, ReviewResult
from analysis.review.runtime import (
    Coalescer, FakeViewerAdapter, NullViewerAdapter, SyncError, SyncHub,
)
from analysis.review.runtime.results import (
    DIFFERENT_COORDINATES_NOTE, ResultCatalog, sample_for_frame,
)
from analysis.review.runtime.server import ReviewServeError, ReviewServer
from analysis.review.store import put_array, put_file
from analysis.review.sync import ReviewTimeline, sync_metadata
from tests.analysis.campaign.conftest import REPO_ROOT, requires_gmx, requires_real_traj

TIMES = [0.0, 10.0, 10.0, 20.0, 30.0]            # a duplicate timestamp, kept


# ═══════════════════════════════════════════════════════════════════════════════
# Synthetic prepared session (every renderer case + a display/analysis mismatch)
# ═══════════════════════════════════════════════════════════════════════════════

def _stored(store: Path, session: Path, text: str, name: str) -> tuple[str, str]:
    tmp = store.parent / name
    tmp.write_text(text)
    dest, sha = put_file(store, tmp)
    tmp.unlink()
    return os.path.relpath(dest, session), sha


def _ref(path, sha, column, fmt="xvg"):
    return StorageRef(path, fmt, column=column, attrs={"sha256": sha})


def make_session(root: Path) -> Path:
    review = root / "review"
    store, session = review / "store", review / "sessions" / "s1"
    session.mkdir(parents=True)
    traj = root / "md.xtc"
    traj.write_bytes(b"not really an xtc")
    tdest, tsha = put_array(store, TIMES, "float64")
    rows = "".join(f"{t:.3f} {1 + i * 0.1:.4f} {2 + i:.4f}\n" for i, t in enumerate(TIMES))
    ts_path, ts_sha = _stored(store, session, rows, "ts.xvg")
    sub_path, sub_sha = _stored(store, session, "0.000 5.0\n20.000 7.0\n", "sub.xvg")
    rmsf_path, rmsf_sha = _stored(store, session, "1 0.11\n2 0.22\n3 0.33\n", "rmsf.xvg")

    def tseries(name, quantity, view, path, sha, extra_axes=(), column=1):
        tax = Axis("time", AxisKind.TIME, "ps", values_ref=_ref(path, sha, 0),
                   alignment=FrameAlignment("one_row_per_frame", trajectory_digest="d", view_ref=view))
        return ResultArray(name, quantity, "nm", axes=[tax, *extra_axes],
                           storage=_ref(path, sha, column), view_ref=view)

    def with_sync(arr):
        absolute = replace(arr, axes=[replace(a, values_ref=replace(
            a.values_ref, path=str(session / a.values_ref.path))) if a.values_ref else a
            for a in arr.axes])
        return ReviewResult(arr, sync=sync_metadata(absolute, TIMES, "d"))

    rg = tseries("rg", "radius_of_gyration", "raw", ts_path, ts_sha)
    rmsd = tseries("rmsd", "rmsd", "fit", ts_path, ts_sha)
    multi = tseries("apl", "area_per_lipid", "raw", ts_path, ts_sha,
                    [Axis("leaflet", AxisKind.CATEGORY, values=["upper", "lower"])])
    sub = tseries("sub", "com_distance", "raw", sub_path, sub_sha)
    rmsf = ResultArray("rmsf", "rmsf", "nm",
                       axes=[Axis("residue", AxisKind.CATEGORY, values=["A:1", "A:2", "A:3"])],
                       storage=_ref(rmsf_path, rmsf_sha, 1), view_ref="fit")
    pore = tseries("pore", "pore_radius", "raw", ts_path, ts_sha,
                   [Axis("z", AxisKind.COORDINATE, "nm")])
    obs = [
        ObservableEntry("rg", "rg", display_name="Radius of gyration", availability="computed",
                        results=[with_sync(rg)], view_ref="raw"),
        ObservableEntry("rmsd", "rmsd-receptor", display_name="RMSD", availability="cached",
                        results=[with_sync(rmsd)], view_ref="fit"),
        ObservableEntry("apl", "apl", availability="computed", results=[with_sync(multi)]),
        ObservableEntry("sub", "com-distance", availability="computed", results=[with_sync(sub)]),
        ObservableEntry("rmsf", "rmsf", availability="computed",
                        results=[ReviewResult(rmsf, sync=sync_metadata(rmsf, TIMES, "d"))]),
        ObservableEntry("pore", "pore", availability="computed", results=[with_sync(pore)]),
        ObservableEntry("rg(selection=annotation:tm_6)", "rg", state="review_required",
                        reason="annotation 'tm_6' is unresolved (residue(s) [401, 402, 403, 404] "
                               "not found)", annotations_used=["tm_6"],
                        blocking_annotations=["tm_6"]),
    ]
    ds = ReviewDataset(
        session_id="synthetic-session", identity_evidence={"synthetic": True},
        system={"system_id": "sys", "membrane_present": True, "components": []},
        sources={"trajectory": {"path": str(traj), "fingerprint": fingerprint_file(traj).to_dict()}},
        display_view={"role": "display", "view_ref": "raw", "kind": "raw", "purpose": "display",
                      "status": "available", "path": str(traj), "decisions": [],
                      "request": {"mode": "raw"}},
        analysis_views=[{"view_ref": "raw", "kind": "raw", "purpose": "intramolecular_shape",
                         "decisions": []},
                        {"view_ref": "fit", "kind": "fitted", "purpose": "intramolecular_shape",
                         "decisions": [{"operation": "make_whole", "applied": True},
                                       {"operation": "fit", "applied": True}]}],
        timeline={"n_frames": len(TIMES), "start_time_ps": 0.0, "end_time_ps": 30.0,
                  "state": "valid", "n_duplicate_groups": 1, "duplicate_groups": [[1, 2]],
                  "times_ref": {"path": os.path.relpath(tdest, session), "sha256": tsha},
                  "display": {"relation": "identical_to_source"}},
        annotations=[{"annotation_id": "tm_1", "kind": "transmembrane_segment", "state": "active",
                      "origin": "build_spec", "reasons": [], "used_by": [], "blocks": []},
                     {"annotation_id": "tm_6", "kind": "transmembrane_segment",
                      "state": "unresolved", "origin": "build_spec",
                      "reasons": ["residue(s) [401, 402, 403, 404] not found"],
                      "used_by": ["rg(selection=annotation:tm_6)"],
                      "blocks": ["rg(selection=annotation:tm_6)"]},
                     {"annotation_id": "pocket", "kind": "binding_site", "state": "proposed",
                      "origin": "derived", "reasons": [], "used_by": [], "blocks": []}],
        diagnostics={"execution": "auto", "status": "complete",
                     "counts": {"info": 1, "warnings": 0, "review_required": 0, "errors": 0},
                     "detectors": {"ran_found_nothing": [{"detector": "timeline_integrity"}]},
                     "findings": [{"code": "box_ok", "severity": "info", "message": "stored text"}]},
        observables=obs)
    (session / "review_dataset.json").write_text(json.dumps(ds.to_dict(), indent=2))
    return session


@pytest.fixture
def session(tmp_path):
    return make_session(tmp_path)


def _hub(viewer=None):
    return SyncHub(ReviewTimeline(TIMES), "synthetic-session", viewer=viewer)


# ═══════════════════════════════════════════════════════════════════════════════
# SyncHub
# ═══════════════════════════════════════════════════════════════════════════════

def test_frame_and_time_mapping():
    hub = _hub()
    ev = hub.set_frame(3)
    assert (ev.frame, ev.time_ps, ev.mapping_status) == (3, 20.0, "frame")
    ev = hub.set_time(30.0)
    assert (ev.frame, ev.mapping_status, ev.delta_ps) == (4, "exact", 0.0)
    ev = hub.set_time(24.0)
    assert (ev.frame, ev.time_ps, ev.mapping_status) == (3, 20.0, "nearest")   # never "exact"
    assert ev.requested_time_ps == 24.0 and ev.delta_ps == 4.0
    ev = hub.set_time(25.0)
    assert ev.frame == 3 and ev.tie                                            # tie → earlier
    ev = hub.set_time(10.0)
    assert ev.frame == 1 and ev.duplicate_frames == [1, 2]                     # kept visible
    ev = hub.set_frame(2)
    assert ev.time_ps == 10.0 and ev.duplicate_frames == [1, 2]
    for bad in (-1, 5, 1.5, "2", True):
        with pytest.raises(SyncError):
            hub.set_frame(bad)
    for bad in (-0.1, 30.5, "x"):
        with pytest.raises(SyncError):
            hub.set_time(bad)                                                  # no extrapolation
    assert hub.frame == 2                                                      # rejected ≠ applied


def test_loop_suppression_with_and_without_acknowledgement():
    for ack in (True, False):
        viewer = FakeViewerAdapter(echo=True, ack=ack)
        hub = _hub(viewer)
        seen = []
        hub.subscribe(seen.append)
        hub.set_frame(3)                              # browser → hub → viewer → echoes 3
        assert [e.frame for e in seen] == [3] and viewer.commands == [(3, 1)]
        assert hub.suppressed == 1
        viewer.user_moves_to(4)                       # independent viewer move
        assert [(e.frame, e.origin) for e in seen] == [(3, "browser"), (4, "viewer")]
        assert len(viewer.commands) == 1              # never echoed back to the viewer
        viewer.user_moves_to(4)                       # re-report of the same frame
        assert len(seen) == 2


def test_stale_acknowledgement_is_consumed():
    viewer = FakeViewerAdapter(echo=False)
    hub = _hub(viewer)
    seen = []
    hub.subscribe(seen.append)
    hub.set_frame(1)
    hub.set_frame(3)                                  # supersedes the first command
    assert hub.viewer_frame_changed(1, ack_sequence=1) is None
    assert hub.viewer_frame_changed(3, ack_sequence=2) is None
    assert hub.frame == 3 and [e.frame for e in seen] == [1, 3]


def test_coalescing_bounds_ui_rate_but_keeps_latest_and_explicit():
    hub = _hub()
    c = Coalescer(min_interval=0.05)
    delivered = []
    now = 0.0
    for i in range(200):                              # 200 viewer frames in 0.2 s
        now = i * 0.001
        ev = hub.viewer_frame_changed(i % 5 if i % 5 != hub.frame else (i + 1) % 5) \
            or hub.last_event
        delivered += c.offer(ev, now)
    delivered += c.due(now + 1.0)                     # the latest state is always delivered
    assert len(delivered) <= 6
    assert delivered[-1].sequence == hub.last_event.sequence
    c2 = Coalescer(min_interval=0.05)
    out = c2.offer(hub.set_frame(0), 0.0)             # delivered immediately
    click = hub.set_frame(2)
    out += c2.offer(click, 0.01)                      # explicit, within the interval
    out += c2.offer(hub.viewer_frame_changed(4), 0.02)
    out += c2.due(0.06)
    assert [e.frame for e in out] == [0, 2, 4]        # the click was not coalesced away


# ═══════════════════════════════════════════════════════════════════════════════
# Result catalog (data contracts per axis signature)
# ═══════════════════════════════════════════════════════════════════════════════

def test_renderer_contracts(session):
    from analysis.review.validate import open_review_session
    ds, val = open_review_session(session)
    assert val.valid, val.problems
    cat = ResultCatalog(ds, session, ReviewTimeline.from_dataset(ds, session))
    by = {cat.describe(r)["instance_id"]: cat.data(r) for r in cat.ids()}
    rg = by["rg"]
    assert rg["renderer"] == "timeseries" and rg["data"]["time_ps"] == TIMES
    assert rg["data"]["series"][0]["values"] == [1.0, 1.1, 1.2, 1.3, 1.4]      # exactly as stored
    assert rg["cursor"]["mode"] == "identity" and sample_for_frame(rg["cursor"], 2) == 2
    apl = by["apl"]
    assert apl["renderer"] == "multiline"
    assert [s["label"] for s in apl["data"]["series"]] == ["upper", "lower"]
    assert apl["data"]["series"][1]["values"] == [2.0, 3.0, 4.0, 5.0, 6.0]
    rmsf = by["rmsf"]
    assert rmsf["renderer"] == "profile" and not rmsf["syncable"]
    assert rmsf["data"]["x"] == ["A:1", "A:2", "A:3"] and rmsf["cursor"]["mode"] == "none"
    pore = by["pore"]
    assert pore["renderer"] == "unsupported" and pore["data"] is None
    assert "(time, coordinate)" in pore["message"] and pore["quantity"] == "pore_radius"
    sub = by["sub"]                                   # 2 samples on a 5-frame timeline
    assert sub["cursor"]["mode"] == "frames"
    assert [sample_for_frame(sub["cursor"], f) for f in range(5)] == [0, None, None, 1, None]
    # display vs analysis coordinates
    assert rg["view"]["same_coordinates"] and not rg["view"]["note"]
    rmsd = by["rmsd"]
    assert not rmsd["view"]["same_coordinates"] and rmsd["view"]["note"] == DIFFERENT_COORDINATES_NOTE
    assert rmsd["view"]["analysis_kind"] == "fitted" and rmsd["view"]["display_kind"] == "raw"
    assert rmsd["view"]["synchronization"] == "simulation time"


# ═══════════════════════════════════════════════════════════════════════════════
# HTTP / SSE
# ═══════════════════════════════════════════════════════════════════════════════

def _req(server, path, body=None, token=True, headers=None):
    h = dict(headers or {})
    if token:
        h["X-SimForge-Token"] = server.token
    data = None
    if body is not None:
        data = body if isinstance(body, bytes) else json.dumps(body).encode()
        h["Content-Type"] = "application/json"
    r = urllib.request.Request(server.base_url + path, data=data, headers=h)
    try:
        with urllib.request.urlopen(r, timeout=10) as resp:
            return resp.status, json.loads(resp.read() or b"null")
    except urllib.error.HTTPError as e:
        return e.code, json.loads(e.read() or b"null")


class _SSE:
    def __init__(self, server, token=None):
        self.resp = urllib.request.urlopen(
            f"{server.base_url}/api/sync/events?token={token or server.token}", timeout=10)

    def next(self, name=None, pred=lambda d: True, timeout=5.0):
        deadline = time.monotonic() + timeout
        event, data = None, None
        while time.monotonic() < deadline:
            try:
                line = self.resp.readline().decode().rstrip("\n")
            except (TimeoutError, socket.timeout):
                break                                 # idle stream (heartbeat is 10 s)
            if line.startswith("event:"):
                event = line.split(":", 1)[1].strip()
            elif line.startswith("data:"):
                data = json.loads(line.split(":", 1)[1])
            elif line == "" and data is not None:
                if (name is None or event == name) and pred(data):
                    return event, data
                event, data = None, None
        raise AssertionError("no matching SSE event")

    def close(self):
        self.resp.close()


def _raw_get(server, path):
    with socket.create_connection((server.host, server.port), timeout=5) as s:
        s.sendall(f"GET {path} HTTP/1.0\r\nHost: 127.0.0.1:{server.port}\r\n"
                  f"X-SimForge-Token: {server.token}\r\n\r\n".encode())
        return s.recv(65536).decode(errors="replace").split("\r\n", 1)[0]


def test_http_sse_roundtrip_and_viewer_path(session):
    viewer = FakeViewerAdapter()
    with ReviewServer(session, viewer=viewer, min_interval=0.0) as srv:
        assert srv.host == "127.0.0.1" and srv.port > 0 and srv.token in srv.url
        st, s = _req(srv, "/api/session")
        assert st == 200 and s["session_id"] == "synthetic-session"
        assert s["timeline"]["times_ps"] == TIMES and s["timeline"]["n_duplicate_groups"] == 1
        assert s["blocked"][0]["instance_id"] == "rg(selection=annotation:tm_6)"
        assert {a["state"] for a in s["annotations"]} == {"active", "unresolved", "proposed"}
        rid = next(r["id"] for r in s["results"] if r["instance_id"] == "rg")
        st, d = _req(srv, f"/api/results/{rid}")
        assert st == 200 and d["data"]["n_samples"] == 5
        sse = _SSE(srv)
        assert sse.next("state")[1]["frame"] == 0
        st, r = _req(srv, "/api/sync/frame", {"frame": 3, "client_event_id": "c1"})
        assert st == 200 and r["event"]["client_event_id"] == "c1"
        _, ev = sse.next("sync", lambda d: d["frame"] == 3)
        assert ev["origin"] == "browser" and ev["time_ps"] == 20.0
        st, r = _req(srv, "/api/sync/time", {"time_ps": 24.0})
        assert r["event"]["frame"] == 3 and r["event"]["mapping_status"] == "nearest"
        viewer.user_moves_to(4)                       # viewer → hub → SSE, no browser POST
        _, ev = sse.next("sync", lambda d: d["origin"] == "viewer")
        assert ev["frame"] == 4 and ev["time_ps"] == 30.0
        assert viewer.commands == [(3, 1), (3, 2)]    # browser moves reached the viewer
        st, err = _req(srv, "/api/sync/frame", {"frame": 99})
        assert st == 400 and "outside" in err["error"]
        st, _ = _req(srv, "/api/sync/frame", b"x" * 5000)
        assert st == 400
        sse.close()


def test_security_boundary(session):
    with ReviewServer(session) as srv:
        assert _req(srv, "/api/session", token=False)[0] == 401
        assert _req(srv, "/api/session", headers={"X-SimForge-Token": "wrong"}, token=False)[0] == 401
        assert _req(srv, "/api/sync/frame", {"frame": 1}, token=False)[0] == 401
        with pytest.raises(urllib.error.HTTPError) as e:
            _SSE(srv, token="wrong")
        assert e.value.code == 401
        assert _req(srv, "/health", token=False)[0] == 200             # no data in it
        # valid token, forbidden paths: only catalog ids, only authorised store files
        for path in ("/api/results/../../../../etc/passwd", "/api/results/%2e%2e%2fetc%2fpasswd",
                     "/api/results//etc/passwd", "/api/results/r99-0",
                     "/static/../review_dataset.json", "/etc/passwd",
                     "/api/results/r0-0?path=/etc/passwd&x=..", "/../outside-file"):
            status = _req(srv, path)[0] if not path.startswith("/api/results/r0-0") else None
            if status is not None:
                assert status in (404,), path
        assert "404" in _raw_get(srv, "/../../etc/passwd")
        assert "404" in _raw_get(srv, "/api/results/../../../etc/passwd")
        # r0-0 ignores any query: it is resolved from the dataset, never a path
        st, d = _req(srv, "/api/results/r0-0?path=/etc/passwd")
        assert st == 200 and d["instance_id"] == "rg"
        # DNS-rebinding style Host
        assert _req(srv, "/api/session", headers={"Host": "evil.example:80"})[0] == 403
        assert "default-src 'self'" in urllib.request.urlopen(srv.base_url + "/").headers[
            "Content-Security-Policy"]


def test_changed_store_file_is_refused(session):
    ds = json.loads((session / "review_dataset.json").read_text())
    p = session / ds["observables"][0]["results"][0]["array"]["storage"]["path"]
    p.chmod(0o644)
    p.write_text(p.read_text() + "999 9\n")
    with pytest.raises(ReviewServeError, match="stale or invalid"):
        ReviewServer(session)


def test_result_changed_while_serving_becomes_unavailable(session):
    with ReviewServer(session) as srv:
        rid = "r1-0"
        d = srv.catalog.describe(rid)
        ds = json.loads((session / "review_dataset.json").read_text())
        p = session / ds["observables"][1]["results"][0]["array"]["storage"]["path"]
        p.write_text("tampered\n")
        st, err = _req(srv, f"/api/results/{rid}")
        assert st == 404 and "missing or changed" in err["error"] and d["id"] == rid


def test_serving_never_recomputes(session, monkeypatch):
    import analysis.campaign.diagnostics as diag
    import analysis.campaign.gmx as gmx
    import analysis.campaign.orchestration.study_analyzer as sa
    import analysis.campaign.trajectory.preprocessor as pre
    import analysis.review.prepare as prep

    def boom(*a, **k):
        raise AssertionError("serving must not run scientific backend code")
    for mod, name in ((sa, "run_analyze"), (diag, "diagnose_system"), (pre, "build_view"),
                      (gmx, "run_gmx"), (sa, "prepare_system"), (prep, "prepare_review"),
                      (sa, "plan_view")):
        monkeypatch.setattr(mod, name, boom)
    before = {p: p.stat().st_mtime_ns for p in session.parent.parent.rglob("*") if p.is_file()}
    with ReviewServer(session, viewer=FakeViewerAdapter()) as srv:
        _, s = _req(srv, "/api/session")
        for r in s["results"]:
            _req(srv, f"/api/results/{r['id']}")
        _req(srv, "/api/sync/frame", {"frame": 2})
        _req(srv, "/api/sync/time", {"time_ps": 20.0})
    after = {p: p.stat().st_mtime_ns for p in session.parent.parent.rglob("*") if p.is_file()}
    assert before == after                                  # nothing written while serving


def test_lifecycle_leaves_no_threads(session):
    base = set(threading.enumerate())
    srv = ReviewServer(session, viewer=FakeViewerAdapter()).start()
    sse = _SSE(srv)
    sse.next("state")
    srv.stop()
    sse.close()
    deadline = time.monotonic() + 3
    while set(threading.enumerate()) - base and time.monotonic() < deadline:
        time.sleep(0.05)
    assert not (set(threading.enumerate()) - base)
    assert not srv.viewer.running
    with pytest.raises(OSError):
        socket.create_connection((srv.host, srv.port), timeout=1)


def test_null_viewer_and_dry_run_refusal(session, tmp_path):
    with ReviewServer(session) as srv:
        assert isinstance(srv.viewer, NullViewerAdapter)
        assert _req(srv, "/api/sync/frame", {"frame": 1})[1]["event"]["frame"] == 1
    ds = json.loads((session / "review_dataset.json").read_text())
    ds["status"] = "dry_run"
    other = tmp_path / "dry"
    other.mkdir()
    (other / "review_dataset.json").write_text(json.dumps(ds))
    with pytest.raises(ReviewServeError, match="only prepared"):
        ReviewServer(other)


def test_static_assets_are_offline():
    static = REPO_ROOT / "analysis" / "review" / "static"
    text = "".join(p.read_text() for p in static.iterdir() if p.is_file())
    for needle in ("http://", "https://", "//cdn", "googleapis", "<script>"):
        assert needle not in text.replace("http://127.0.0.1", ""), needle


# ═══════════════════════════════════════════════════════════════════════════════
# CLI + real membrane session
# ═══════════════════════════════════════════════════════════════════════════════

def _serve_cli(target, *extra):
    p = subprocess.Popen([sys.executable, "-m", "cli", "trajectory", "serve", str(target),
                          "--no-browser", *extra], cwd=REPO_ROOT, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, text=True)
    import queue
    lines: "queue.Queue[str]" = queue.Queue()
    threading.Thread(target=lambda: [lines.put(ln) for ln in p.stdout], daemon=True).start()
    url, deadline = None, time.monotonic() + 30
    while url is None and time.monotonic() < deadline:
        try:
            line = lines.get(timeout=max(0.0, deadline - time.monotonic()))
        except queue.Empty:
            break
        if "token=" in line:
            url = line.strip()
    return p, url


def test_cli_serve_and_clean_interrupt(session):
    p, url = _serve_cli(session)
    try:
        assert url and url.startswith("http://127.0.0.1:")
        base, token = url.split("/?token=")
        with urllib.request.urlopen(base + "/health", timeout=5) as r:
            assert json.loads(r.read())["ok"]
        r = urllib.request.Request(base + "/api/session", headers={"X-SimForge-Token": token})
        with urllib.request.urlopen(r, timeout=5) as resp:
            assert json.loads(resp.read())["session_id"] == "synthetic-session"
    finally:
        p.send_signal(signal.SIGINT)
        try:
            code = p.wait(timeout=15)
        except subprocess.TimeoutExpired:
            p.kill()
            raise
        assert code == 0


def test_cli_refuses_invalid_session(tmp_path):
    r = subprocess.run([sys.executable, "-m", "cli", "trajectory", "serve", str(tmp_path),
                        "--no-browser"], cwd=REPO_ROOT, capture_output=True, text=True, timeout=60)
    assert r.returncode == 1 and "Refusing to serve" in r.stdout


@requires_gmx
@requires_real_traj
def test_real_membrane_session_served(tmp_path):
    from analysis.review import DisplayRequest, ReviewRequest, parse_show, prepare_review
    from tests.analysis.review.test_review import _run_dir
    run = _run_dir(tmp_path)
    prep = prepare_review(run, ReviewRequest(parse_show([
        "rg(selection=component:receptor)", "rg(selection=annotation:tm_1)",
        "rg(selection=annotation:tm_6)", "rmsd-receptor"]), DisplayRequest.parse("raw")))
    with ReviewServer(prep.dataset_path, viewer=FakeViewerAdapter()) as srv:
        _, s = _req(srv, "/api/session")
        assert s["display_view"]["kind"] == "raw" and s["display_view"]["decisions"] == []
        assert s["timeline"]["n_frames"] == 51
        assert [b["instance_id"] for b in s["blocked"]] == ["rg(selection=annotation:tm_6)"]
        assert "unresolved" in s["blocked"][0]["reason"]
        states = {a["annotation_id"]: a["state"] for a in s["annotations"]}
        assert states["tm_6"] == "unresolved" and states["tm_1"] == "active"
        panels = {r["instance_id"]: r for r in s["results"]}
        assert len(panels) == 3 and all(r["syncable"] for r in panels.values())
        assert panels["rmsd-receptor"]["view"]["note"] == DIFFERENT_COORDINATES_NOTE   # whole ≠ raw
        assert panels["rg(selection=component:receptor)"]["view"]["same_coordinates"]
        for r in panels.values():
            st, d = _req(srv, f"/api/results/{r['id']}")
            assert st == 200 and d["data"]["n_samples"] == 51
        _req(srv, "/api/sync/frame", {"frame": 25})
        assert srv.viewer.frame == 25
    assert not list(run.rglob("centered.xtc")) and not list(run.rglob("fitted.xtc"))
