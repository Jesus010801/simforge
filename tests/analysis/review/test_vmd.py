"""Phase 11 — VMD viewer client: protocol, Tcl bridge, adapter, sync, real VMD."""
from __future__ import annotations

import json
import os
import shutil
import socket
import subprocess
import sys
import time
from pathlib import Path

import pytest

from analysis.review.runtime import SyncHub
from analysis.review.runtime import vmd_protocol as proto
from analysis.review.runtime.vmd import (
    BRIDGE_TCL, VMDError, VMDViewerAdapter, ViewerState, display_inputs, find_vmd,
)
from analysis.review.sync import ReviewTimeline
from tests.analysis.campaign.conftest import (
    REAL_GRO, REAL_XTC, REPO_ROOT, requires_real_traj,
)

HAVE_VMD = shutil.which("vmd") is not None
HAVE_TCLSH = shutil.which("tclsh") is not None
requires_vmd = pytest.mark.skipif(not HAVE_VMD, reason="VMD not installed (integration test)")
TOKEN = "0123456789abcdef0123456789abcdef"


def wait(pred, timeout=5.0):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if pred():
            return True
        time.sleep(0.005)
    return False


# ═══════════════════════════════════════════════════════════════════════════════
# Protocol (Python side)
# ═══════════════════════════════════════════════════════════════════════════════

def test_protocol_parses_the_allow_list_only():
    P = proto.parse_bridge_line
    assert P(f"HELLO simforge-vmd/1 {TOKEN}\n").fields == ("simforge-vmd/1", TOKEN)
    assert P(b"READY 0 283331 51 2.0.0a5\n").fields == (0, 283331, 51, "2.0.0a5")
    assert P("FRAME 1000 137").fields == (1000, 137)
    assert P("FRAME 12 NONE").fields == (12, None)
    assert P("LOADING").kind == "LOADING" and P("BYE").kind == "BYE"
    assert P("ERROR load_failed cannot open file").fields == ("load_failed", "cannot open file")
    bad = [
        f"HELLO simforge-vmd/2 {TOKEN}",           # wrong protocol version
        "HELLO simforge-vmd/1 not-a-token",         # malformed token
        "FRAME -1 NONE", "FRAME 1.5 NONE", "FRAME 3", "FRAME 3 ack", "FRAME 01 NONE",
        "READY 0 10 x 2.0", "READY 0 10 5 ver;sion",
        "EVAL puts hi", "GOTO 1 2",                  # unknown / wrong direction
        "", "FRAME 1  NONE", "BYE now", "ERROR BadCode x",
        "FRAME 1 NONE\x00", "frame 1 NONE",
        "FRAME " + "9" * 600,                        # oversized
    ]
    for line in bad:
        with pytest.raises(proto.ProtocolError):
            P(line)
    with pytest.raises(proto.ProtocolError):
        P(b"FRAME 1 N\xc3\x93NE")                   # non-ASCII
    assert proto.goto(137, 1000) == b"GOTO 137 1000\n"
    for bad_args in ((-1, 0), (1, -2), (True, 1), (1, 2.0)):
        with pytest.raises(proto.ProtocolError):
            proto.goto(*bad_args)


# ═══════════════════════════════════════════════════════════════════════════════
# Tcl bridge helpers under plain tclsh (no VMD)
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.mark.skipif(not HAVE_TCLSH, reason="tclsh not available")
def test_tcl_bridge_parsing_and_acknowledgement_logic(tmp_path):
    script = tmp_path / "t.tcl"
    script.write_text(f"""
set ::simforge_test_mode 1
source {{{BRIDGE_TCL}}}
foreach line {{
    "WELCOME simforge-vmd/1" "WELCOME simforge-vmd/2" "GOTO 137 1000" "GOTO 1 -2"
    "GOTO 1 2 3" "GOTO a 2" "BYE" "BYE now" "EVAL exit" "GOTO 1 \\[exit\\]"
    "puts hi" ""
}} {{ puts "[list $line] -> [::simforge::parse_command $line]" }}
puts "long -> [::simforge::parse_command [string repeat x 600]]"
set ::simforge::pending_seq 5
set ::simforge::pending_frame 10
puts "r1 [::simforge::frame_report 9]"
puts "r2 [::simforge::frame_report 10]"
puts "r3 [::simforge::frame_report 10]"
puts "clean [::simforge::clean "a\\nb\\x01c"]"
""")
    r = subprocess.run(["tclsh", str(script)], capture_output=True, text=True, timeout=30)
    assert r.returncode == 0, r.stderr
    out = dict(line.split(" -> ", 1) for line in r.stdout.splitlines() if " -> " in line)
    assert out['{WELCOME simforge-vmd/1}'] == "WELCOME"
    assert out['{WELCOME simforge-vmd/2}'] == "ERROR bad_welcome"
    assert out['{GOTO 137 1000}'] == "GOTO 137 1000"
    for k in ('{GOTO 1 -2}', '{GOTO 1 2 3}', '{GOTO a 2}', '{GOTO 1 [exit]}'):
        assert out[k] == "ERROR bad_goto", k
    assert out["BYE"] == "BYE" and out['{BYE now}'] == "ERROR bad_bye"
    assert out['{EVAL exit}'] == out['{puts hi}'] == out["{}"] == "ERROR unknown_command"
    assert out["long"] == "ERROR line_too_long"
    lines = r.stdout.splitlines()
    assert "r1 FRAME 9 NONE" in lines             # not the pending frame: independent
    assert "r2 FRAME 10 5" in lines               # acknowledges the pending command…
    assert "r3 FRAME 10 NONE" in lines            # …exactly once
    assert "clean a b c" in lines


# ═══════════════════════════════════════════════════════════════════════════════
# Adapter against a fake bridge peer (no VMD)
# ═══════════════════════════════════════════════════════════════════════════════

class Peer:
    """Plays the Tcl bridge's side of the protocol."""

    def __init__(self, adapter, token=None):
        self.s = socket.create_connection(("127.0.0.1", adapter.port), timeout=5)
        self.f = self.s.makefile("rwb")
        self.frame = None
        self.pending = None
        if token is not False:
            self.send(f"HELLO {proto.PROTOCOL} {token or adapter.token}")

    def send(self, line):
        self.f.write(line.encode() + b"\n")
        self.f.flush()

    def recv(self):
        return self.f.readline().decode().strip()

    def ready(self, frames=5, atoms=10):
        assert self.recv() == f"WELCOME {proto.PROTOCOL}"
        self.send("LOADING")
        self.send(f"READY 0 {atoms} {frames} 2.0.0a5")

    def serve_goto(self):
        """Receive one GOTO and behave like VMD (move, ack)."""
        _, seq, frame = self.recv().split()
        self.frame = int(frame)
        self.send(f"FRAME {frame} {seq}")
        return int(seq), int(frame)

    def close(self):
        self.f.close()                               # the makefile holds the socket open
        self.s.close()


def adapter(frames=5, atoms=10, **kw):
    a = VMDViewerAdapter(structure="s.gro", structure_type="gro", trajectory="t.xtc",
                         trajectory_type="xtc", n_frames=frames, n_atoms=atoms, launch=False,
                         startup_timeout=10, **kw)
    return a


def test_adapter_handshake_goto_ack_and_independent_move():
    a = adapter()
    seen = []
    a.on_frame_changed(lambda f, ack_sequence=None: seen.append((f, ack_sequence)))
    a.start({"current_frame": lambda: 2})
    try:
        p = Peer(a)
        p.ready()
        assert a.wait_ready(5) and a.status()["numframes"] == 5
        assert p.serve_goto() == (0, 2)            # hub frame is authoritative on connect
        assert wait(lambda: a.current_frame == 2) and seen == []   # own sync swallowed
        a.goto_frame(4, 7)
        assert p.serve_goto() == (7, 4)
        assert wait(lambda: seen == [(4, 7)])
        p.send("FRAME 1 NONE")                      # user moved VMD
        assert wait(lambda: seen[-1] == (1, None))
        st = a.status()
        assert st["state"] == "ready" and "token" not in json.dumps(st).lower()
        assert a.token not in json.dumps(st)
    finally:
        a.stop()
    assert a.state == ViewerState.STOPPED


@pytest.mark.parametrize("token", ["f" * 32, False])
def test_bad_handshake_is_rejected(token):
    a = adapter()
    a.start({})
    try:
        p = Peer(a, token=token)
        if token is False:
            p.send("READY 0 10 5 2.0")             # skipping HELLO
        assert a.wait_ready(5) is False
        assert a.state == "error" and "handshake" in a.reason
        assert p.recv() == "BYE"                    # no WELCOME was ever sent
    finally:
        a.stop()


@pytest.mark.parametrize("line,needle", [
    ("EVAL exec rm -rf /", "unknown message"),
    ("FRAME x NONE", "frame"),
    ("A" * 700, "oversized"),
])
def test_malformed_peer_input_stops_the_viewer(line, needle):
    a = adapter()
    a.start({})
    try:
        p = Peer(a)
        p.ready()
        assert a.wait_ready(5)
        p.serve_goto()
        p.send(line)
        assert wait(lambda: a.state == "error") and needle in a.reason
    finally:
        a.stop()


@pytest.mark.parametrize("frames,atoms,needle", [(4, 10, "4 frames"), (5, 11, "11 atoms")])
def test_frame_or_atom_count_mismatch_refuses_sync(frames, atoms, needle):
    a = adapter(frames=5, atoms=10)
    a.start({})
    try:
        p = Peer(a)
        p.ready(frames=frames, atoms=atoms)
        assert a.wait_ready(5) is False
        assert a.state == "error" and needle in a.reason
        assert p.recv() == "BYE"                    # never truncated / offset
        a.goto_frame(1, 3)                          # not applied, not claimed
        assert a.status()["last_command_applied"] is False
    finally:
        a.stop()


def test_disconnect_and_commands_before_ready():
    a = adapter()
    a.start({"current_frame": lambda: 0})
    try:
        p = Peer(a)
        a.goto_frame(3, 9)                          # VMD still loading: not applied
        assert a.status()["last_command_applied"] is False
        p.ready()
        assert a.wait_ready(5)
        assert p.serve_goto() == (9, 3)             # the pending command, with its sequence
        p.close()                                   # VMD went away
        assert wait(lambda: a.state == "disconnected")
        a.goto_frame(1, 10)
        assert a.status()["last_command_applied"] is False
    finally:
        a.stop()


def test_hub_integration_ack_supersession_and_in_flight_reports():
    a = adapter(frames=5, atoms=10)
    hub = SyncHub(ReviewTimeline([0.0, 10.0, 20.0, 30.0, 40.0]), "s", viewer=a)
    events = []
    hub.subscribe(events.append)
    a.start({"current_frame": lambda: hub.frame})
    try:
        p = Peer(a)
        p.ready()
        assert a.wait_ready(5)
        p.serve_goto()                              # initial sync (sequence 0)
        hub.set_frame(3)                            # browser-like
        p.serve_goto()
        assert wait(lambda: hub.suppressed == 1)
        assert [e.frame for e in events] == [3]     # one logical transition
        p.send("FRAME 4 NONE")                      # independent VMD move
        assert wait(lambda: hub.frame == 4 and events[-1].origin == "viewer")
        # burst 1, 2, 3, 0 — acks arrive late and in order
        for f in (1, 2, 3, 0):
            hub.set_frame(f)
        p.send("FRAME 4 NONE")                      # in flight before VMD saw the burst
        for _ in range(4):
            p.serve_goto()
        assert wait(lambda: a.last_ack == hub.last_event.sequence)
        assert hub.frame == 0 and p.frame == 0      # converged on the latest command
        assert sum(e.origin == "viewer" for e in events) == 1
    finally:
        a.stop()


def test_display_inputs_refuse_stale_or_unprepared(tmp_path):
    from tests.analysis.review.test_runtime import make_session
    from analysis.review.validate import load_review_dataset
    session = make_session(tmp_path)
    ds = load_review_dataset(session)
    with pytest.raises(VMDError, match="structure"):
        display_inputs(ds, session)                 # the synthetic session has no structure
    gro = tmp_path / "md.gro"
    gro.write_text("x\n")
    from analysis.campaign.fingerprint import fingerprint_file
    ds.sources["structure"] = {"path": str(gro), "fingerprint": fingerprint_file(gro).to_dict()}
    inp = display_inputs(ds, session)
    assert (inp["n_frames"], inp["structure_type"], inp["trajectory_type"]) == (5, "gro", "xtc")
    Path(ds.sources["trajectory"]["path"]).write_bytes(b"changed")
    with pytest.raises(VMDError, match="changed since"):
        display_inputs(ds, session)
    ds.display_view["status"] = "review_required"
    with pytest.raises(VMDError, match="not available"):
        display_inputs(ds, session)


def test_find_vmd_is_conservative(tmp_path, monkeypatch):
    with pytest.raises(VMDError, match="not an executable"):
        find_vmd(str(tmp_path / "nope"))
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(VMDError, match="not found on PATH"):
        find_vmd()


def test_cli_reports_missing_vmd(tmp_path):
    from tests.analysis.review.test_runtime import make_session
    session = make_session(tmp_path)
    r = subprocess.run([sys.executable, "-m", "cli", "trajectory", "serve", str(session),
                        "--no-browser", "--viewer", "vmd", "--vmd", str(tmp_path / "no-vmd")],
                       cwd=REPO_ROOT, capture_output=True, text=True, timeout=60)
    assert r.returncode == 1 and "not an executable" in r.stdout


# ═══════════════════════════════════════════════════════════════════════════════
# Real VMD (text mode) on the bundled 51-frame membrane trajectory
# ═══════════════════════════════════════════════════════════════════════════════

def _real_adapter(root: Path, frames=51, atoms=283331):
    d = root / "dir with spaces & 'quotes'"
    d.mkdir()
    (d / "ref struct.gro").symlink_to(REAL_GRO)
    (d / "traj $x [1].xtc").symlink_to(REAL_XTC)
    return VMDViewerAdapter(structure=str(d / "ref struct.gro"), structure_type="gro",
                            trajectory=str(d / "traj $x [1].xtc"), trajectory_type="xtc",
                            n_frames=frames, n_atoms=atoms, headless=True, startup_timeout=120)


def _console(a, cmd):
    a._proc.stdin.write((cmd + "\n").encode())       # test-only: act as the VMD user
    a._proc.stdin.flush()


@requires_vmd
@requires_real_traj
def test_real_vmd_sync(tmp_path, monkeypatch):
    import analysis.campaign.gmx as gmx
    monkeypatch.setattr(gmx, "run_gmx", lambda *a, **k: pytest.fail("GROMACS called"))
    before = (REAL_XTC.stat().st_mtime_ns, REAL_GRO.stat().st_mtime_ns)
    a = _real_adapter(tmp_path)
    hub = SyncHub(ReviewTimeline([20.0 * i for i in range(51)]), "s", viewer=a)
    events = []
    hub.subscribe(events.append)
    a.start({"current_frame": lambda: hub.frame})
    try:
        assert a.wait_ready(120), a.status()
        st = a.status()
        assert (st["numframes"], st["numatoms"]) == (51, 283331) and st["version"]
        assert wait(lambda: a.current_frame == 0)    # VMD starts on its last frame; synced to 0
        for f in (50, 0):                            # last and first frame, zero-based
            ev = hub.set_frame(f)
            assert wait(lambda: a.current_frame == f and a.last_ack == ev.sequence)
        _console(a, "animate goto 12")               # independent VMD move
        assert wait(lambda: hub.frame == 12) and events[-1].origin == "viewer"
        for f in (40, 41, 42, 43):                   # burst
            hub.set_frame(f)
        assert wait(lambda: a.current_frame == 43 and a.last_ack == hub.last_event.sequence)
        assert hub.frame == 43 and sum(e.origin == "viewer" for e in events) == 1
        n0 = len(events)
        _console(a, "animate goto 0; animate speed 1.0; animate style once; animate forward")
        assert wait(lambda: hub.frame == 50, 20)     # playback reaches the hub
        assert len(events) - n0 > 40
        _console(a, "animate goto 0; animate forward")
        time.sleep(0.02)
        ev = hub.set_frame(25)                       # browser during playback
        assert wait(lambda: a.last_ack == ev.sequence, 5)
        time.sleep(0.5)
        assert hub.frame == a.current_frame          # converged, whatever VMD does next
    finally:
        a.stop()
    assert a.state == "stopped" and a._proc.returncode is not None
    assert not a.runtime_dir.exists()
    assert (REAL_XTC.stat().st_mtime_ns, REAL_GRO.stat().st_mtime_ns) == before


@requires_vmd
@requires_real_traj
def test_real_vmd_frame_count_mismatch_fails(tmp_path):
    a = _real_adapter(tmp_path, frames=50)
    a.start({})
    try:
        assert a.wait_ready(120) is False
        assert a.state == "error" and "51 frames" in a.reason and "50" in a.reason
    finally:
        a.stop()
    assert a._proc.returncode is not None


@requires_vmd
@requires_real_traj
def test_real_vmd_unexpected_exit_is_reported(tmp_path):
    a = _real_adapter(tmp_path)
    a.start({})
    try:
        assert a.wait_ready(120)
        _console(a, "quit")
        assert wait(lambda: a.state == "disconnected", 20)
        a.goto_frame(3, 99)
        assert a.status()["last_command_applied"] is False
    finally:
        a.stop()


# ═══════════════════════════════════════════════════════════════════════════════
# Dashboard server + real VMD (HTTP/SSE), prepared GLP-1R session
# ═══════════════════════════════════════════════════════════════════════════════

@requires_vmd
@requires_real_traj
@pytest.mark.skipif(shutil.which("gmx") is None, reason="needs gmx to prepare the session")
def test_dashboard_and_vmd_end_to_end(tmp_path, monkeypatch):
    from analysis.review import DisplayRequest, ReviewRequest, parse_show, prepare_review
    from analysis.review.runtime.server import ReviewServer
    from tests.analysis.review.test_review import _run_dir
    from tests.analysis.review.test_runtime import _SSE, _req
    run = _run_dir(tmp_path)
    prep = prepare_review(run, ReviewRequest(parse_show([
        "rg(selection=component:receptor)", "rg(selection=annotation:tm_6)"]),
        DisplayRequest.parse("raw")))
    import analysis.campaign.gmx as gmx
    import analysis.campaign.orchestration.study_analyzer as sa
    import analysis.campaign.trajectory.policy as pol
    import analysis.campaign.trajectory.preprocessor as pre
    import analysis.review.prepare as rp

    def boom(*a, **k):
        raise AssertionError("scientific backend called while serving with VMD")
    for mod, name in ((gmx, "run_gmx"), (sa, "run_analyze"), (sa, "plan_view"),
                      (sa, "prepare_system"), (pre, "build_view"), (rp, "prepare_review"),
                      (pol, "plan_preprocessing")):
        monkeypatch.setattr(mod, name, boom)
    with ReviewServer(prep.dataset_path, viewer="vmd", viewer_options={"headless": True},
                      min_interval=0.0) as srv:
        assert srv.viewer.wait_ready(120)
        _, s = _req(srv, "/api/session")
        assert s["viewer"]["kind"] == "vmd" and s["viewer"]["state"] == "ready"
        assert s["viewer"]["numframes"] == 51 and s["display_view"]["kind"] == "raw"
        assert [b["instance_id"] for b in s["blocked"]] == ["rg(selection=annotation:tm_6)"]
        sse = _SSE(srv)
        sse.next("state")
        st, r = _req(srv, "/api/sync/frame", {"frame": 33})          # browser → VMD
        assert wait(lambda: srv.viewer.current_frame == 33)
        assert r["viewer"]["kind"] == "vmd"
        _, ev = sse.next("sync", lambda d: d["frame"] == 33)
        _console(srv.viewer, "animate goto 7")                        # VMD → browser
        _, ev = sse.next("sync", lambda d: d["origin"] == "viewer")
        assert ev["frame"] == 7 and ev["time_ps"] == 140.0
        with pytest.raises(AssertionError):
            sse.next("sync", timeout=0.5)                             # no echo / duplicate
        sse.close()


def _vmd_bridge_pids() -> set[str]:
    r = subprocess.run(["pgrep", "-f", "simforge-vmd-.*/bridge.tcl"], capture_output=True, text=True)
    return set(r.stdout.split())


@requires_vmd
@requires_real_traj
@pytest.mark.skipif(shutil.which("gmx") is None, reason="needs gmx to prepare the session")
def test_cli_serve_with_vmd_and_clean_interrupt(tmp_path):
    import signal
    from analysis.review import DisplayRequest, ReviewRequest, parse_show, prepare_review
    from tests.analysis.review.test_review import _run_dir
    from tests.analysis.review.test_runtime import _serve_cli
    run = _run_dir(tmp_path)
    prep = prepare_review(run, ReviewRequest(parse_show(["rg(selection=component:receptor)"]),
                                             DisplayRequest.parse("raw")))
    before = _vmd_bridge_pids()
    p, url = _serve_cli(prep.dataset_path, "--viewer", "vmd", "--vmd-headless")
    try:
        assert url
        assert wait(lambda: len(_vmd_bridge_pids() - before) == 1, 60)   # one owned VMD
        import urllib.request
        base, token = url.split("/?token=")

        def viewer_state():
            req = urllib.request.Request(base + "/api/session", headers={"X-SimForge-Token": token})
            with urllib.request.urlopen(req, timeout=5) as r:
                return json.loads(r.read())["viewer"]
        assert wait(lambda: viewer_state()["state"] == "ready", 120)
        v = viewer_state()
        assert v["numframes"] == 51 and "token" not in json.dumps(v)
    finally:
        p.send_signal(signal.SIGINT)
        assert p.wait(timeout=30) == 0
    assert wait(lambda: not (_vmd_bridge_pids() - before), 10)          # no orphan VMD
