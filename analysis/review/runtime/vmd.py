"""VMD as a viewer client of the SyncHub (Phase 11).

VMD only *displays* the ReviewDataset's already prepared display trajectory
and reports its current frame.  This adapter never chooses, builds or
transforms a view, never runs GROMACS and never computes anything: it
validates that the display inputs are the ones the session describes,
launches a SimForge-owned VMD process with a static Tcl bridge, and relays
frames over an authenticated loopback protocol (``vmd_protocol``).

Lifecycle states: starting → loading → ready, or error / disconnected /
stopped.  Startup runs on a background thread; ``wait_ready()`` blocks.
"""
from __future__ import annotations

import hmac
import os
import secrets
import shutil
import socket
import subprocess
import tempfile
import threading
import time
from pathlib import Path
from typing import Callable, Optional

from analysis.review.runtime import vmd_protocol as proto

BRIDGE_TCL = Path(__file__).resolve().parent / "vmd_bridge.tcl"
STRUCTURE_TYPES = {".gro": "gro", ".pdb": "pdb"}
TRAJECTORY_TYPES = {".xtc": "xtc", ".trr": "trr"}


class VMDError(RuntimeError):
    """VMD cannot be used for this session (absent, stale inputs, bad load)."""


class ViewerState:
    STARTING = "starting"
    LOADING = "loading"
    READY = "ready"
    DISCONNECTED = "disconnected"
    ERROR = "error"
    STOPPED = "stopped"


def find_vmd(explicit: Optional[str] = None) -> str:
    """An explicit path, else ``vmd`` on PATH — no other directories are searched."""
    if explicit:
        p = Path(explicit)
        if p.is_file() and os.access(p, os.X_OK):
            return str(p)
        raise VMDError(f"--vmd {explicit}: not an executable file")
    found = shutil.which("vmd")
    if not found:
        raise VMDError("VMD was not found on PATH; install VMD or pass --vmd /path/to/vmd "
                       "(the dashboard works without it: --viewer null)")
    return found


def display_inputs(dataset, session_dir) -> dict:
    """The exact files VMD must show, re-validated against the dataset.

    Refuses (never prepares) when the display view is not an available,
    already existing file whose fingerprint still matches.
    """
    from analysis.campaign.fingerprint import fingerprint_file
    dv = dataset.display_view
    if dv.get("status") != "available":
        raise VMDError(f"display view is {dv.get('status')!r}, not available — "
                       f"{dv.get('reason', '')}".rstrip(" —"))
    base = Path(session_dir)
    path = Path(dv.get("path") or "")
    traj = path if path.is_absolute() else (base / path)
    if not traj.is_file():
        raise VMDError(f"display trajectory {traj} does not exist")
    expected_fp = ((dataset.sources.get("trajectory") or {}).get("fingerprint")
                   if dv.get("kind") == "raw" else dv.get("output_fingerprint"))
    if not expected_fp:
        raise VMDError("the display trajectory has no recorded fingerprint")
    if fingerprint_file(traj, strong=expected_fp.get("mode") == "strong").digest != \
            expected_fp.get("digest"):
        raise VMDError(f"display trajectory {traj.name} changed since the session was prepared")
    st = (dataset.sources.get("structure") or {})
    structure = Path(st.get("path") or "")
    if not structure.is_file():
        raise VMDError("the session has no reference structure file for VMD")
    sfp = st.get("fingerprint") or {}
    if sfp and fingerprint_file(structure, strong=sfp.get("mode") == "strong").digest != \
            sfp.get("digest"):
        raise VMDError(f"structure {structure.name} changed since the session was prepared")
    stype = STRUCTURE_TYPES.get(structure.suffix.lower())
    ttype = TRAJECTORY_TYPES.get(traj.suffix.lower())
    if not stype or not ttype:
        raise VMDError(f"unsupported file types for VMD: {structure.suffix} / {traj.suffix}")
    n_frames = dataset.timeline.get("n_frames")
    if (dataset.timeline.get("display") or {}).get("relation") == "own_timeline":
        n_frames = dataset.timeline["display"].get("n_frames")
    if not isinstance(n_frames, int) or n_frames <= 0:
        raise VMDError("the display frame count is unknown")
    groups = ((dataset.system.get("semantic_index") or {}).get("groups") or [])
    n_atoms = next((g.get("n_atoms") for g in groups if g.get("name") == "System"), None)
    return {"structure": str(structure), "structure_type": stype, "trajectory": str(traj),
            "trajectory_type": ttype, "n_frames": n_frames, "n_atoms": n_atoms,
            "view_ref": dv.get("view_ref"), "view_kind": dv.get("kind")}


class VMDViewerAdapter:
    kind = "vmd"
    acknowledges = True             # the bridge acks every GOTO and reports in order

    def __init__(self, *, structure: str, structure_type: str, trajectory: str,
                 trajectory_type: str, n_frames: int, n_atoms: Optional[int] = None,
                 vmd: Optional[str] = None, headless: bool = False,
                 startup_timeout: float = 300.0, runtime_parent: Optional[str] = None,
                 launch: bool = True, view_kind: Optional[str] = None):
        self.structure, self.structure_type = structure, structure_type
        self.trajectory, self.trajectory_type = trajectory, trajectory_type
        self.n_frames, self.n_atoms, self.view_kind = n_frames, n_atoms, view_kind
        self.vmd, self.headless, self.launch = vmd, headless, launch
        self.startup_timeout, self.runtime_parent = startup_timeout, runtime_parent
        self.token = secrets.token_hex(16)             # runtime only; never shown / persisted
        self.state = ViewerState.STOPPED
        self.reason = ""
        self.info: dict = {}
        self.current_frame: Optional[int] = None
        self.last_command: Optional[tuple[int, int]] = None
        self.last_ack: Optional[int] = None
        self._callback: Optional[Callable] = None
        self._current_frame_provider: Callable[[], int] = lambda: 0
        self._listener: Optional[socket.socket] = None
        self._conn: Optional[socket.socket] = None
        self._send_lock = threading.Lock()
        self._ready = threading.Event()
        self._done = threading.Event()
        self._thread: Optional[threading.Thread] = None
        self._proc: Optional[subprocess.Popen] = None
        self.runtime_dir: Optional[Path] = None
        self.log_path: Optional[Path] = None
        self.port: Optional[int] = None

    @classmethod
    def for_dataset(cls, dataset, session_dir, **kw) -> "VMDViewerAdapter":
        return cls(**{k: v for k, v in display_inputs(dataset, session_dir).items()
                      if k != "view_ref"}, **kw)

    # ── ViewerAdapter interface ────────────────────────────────────────────
    def on_frame_changed(self, callback: Callable) -> None:
        self._callback = callback

    def start(self, session: dict) -> None:
        provider = session.get("current_frame")
        if callable(provider):
            self._current_frame_provider = provider
        self.state, self.reason = ViewerState.STARTING, ""
        self._listener = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self._listener.bind(("127.0.0.1", 0))
        self._listener.listen(1)
        self.port = self._listener.getsockname()[1]
        self.runtime_dir = Path(tempfile.mkdtemp(prefix="simforge-vmd-", dir=self.runtime_parent))
        self.log_path = self.runtime_dir / "vmd.log"
        if self.launch:
            try:
                self._launch()
            except (VMDError, OSError) as exc:
                self._fail(ViewerState.ERROR, str(exc))
                return
        self._thread = threading.Thread(target=self._run, name="simforge-vmd-bridge", daemon=True)
        self._thread.start()

    def goto_frame(self, frame: int, event_sequence: int) -> None:
        self.last_command = (frame, event_sequence)
        if self.state != ViewerState.READY:
            return                       # not applied; status() says so (applied=False)
        self._send(proto.goto(event_sequence, frame))

    def status(self) -> dict:
        applied = (self.state == ViewerState.READY and self.last_command is not None
                   and self.current_frame == self.last_command[0])
        return {"kind": self.kind, "state": self.state, "connected": self.state == ViewerState.READY,
                "reason": self.reason, "version": self.info.get("version"),
                "molid": self.info.get("molid"), "numatoms": self.info.get("numatoms"),
                "numframes": self.info.get("numframes"), "expected_frames": self.n_frames,
                "expected_atoms": self.n_atoms, "current_frame": self.current_frame,
                "last_command": ({"frame": self.last_command[0], "sequence": self.last_command[1]}
                                 if self.last_command else None),
                "last_ack": self.last_ack, "last_command_applied": applied,
                "headless": self.headless, "display_view_kind": self.view_kind,
                "pid": self._proc.pid if self._proc else None,
                "log": str(self.log_path) if self.log_path else None}

    def stop(self) -> None:
        was = self.state
        self._done.set()
        if self._conn is not None:
            try:
                self._send(proto.bye())
            except OSError:
                pass
        self._close_sockets()
        if self._proc is not None:                     # only the process we launched
            try:
                if self._proc.stdin:
                    self._proc.stdin.close()
                self._proc.wait(timeout=10)
            except subprocess.TimeoutExpired:
                self._proc.terminate()
                try:
                    self._proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    self._proc.kill()
                    self._proc.wait(timeout=5)
        if self._thread is not None:
            self._thread.join(timeout=5)
        if self.runtime_dir is not None:
            bridge = self.runtime_dir / "bridge.tcl"
            bridge.unlink(missing_ok=True)
            if was != ViewerState.ERROR:               # keep logs of a failed start
                shutil.rmtree(self.runtime_dir, ignore_errors=True)
        if self.state != ViewerState.ERROR:
            self.state = ViewerState.STOPPED

    def wait_ready(self, timeout: Optional[float] = None) -> bool:
        self._ready.wait(timeout if timeout is not None else self.startup_timeout)
        return self.state == ViewerState.READY

    # ── internals ──────────────────────────────────────────────────────────
    def _launch(self) -> None:
        exe = find_vmd(self.vmd)
        bridge = self.runtime_dir / "bridge.tcl"
        shutil.copyfile(BRIDGE_TCL, bridge)
        env = dict(os.environ,
                   SIMFORGE_VMD_PORT=str(self.port), SIMFORGE_VMD_TOKEN=self.token,
                   SIMFORGE_VMD_STRUCTURE=self.structure,
                   SIMFORGE_VMD_STRUCTURE_TYPE=self.structure_type,
                   SIMFORGE_VMD_TRAJECTORY=self.trajectory,
                   SIMFORGE_VMD_TRAJECTORY_TYPE=self.trajectory_type)
        argv = [exe] + (["-dispdev", "text"] if self.headless else []) + ["-e", str(bridge)]
        log = open(self.log_path, "wb")                  # file redirection: no pipe to fill
        try:
            # stdin stays open: text-mode VMD exits on console EOF; start_new_session
            # keeps a terminal Ctrl-C from reaching VMD before SimForge stops it
            self._proc = subprocess.Popen(argv, env=env, stdin=subprocess.PIPE, stdout=log,
                                          stderr=subprocess.STDOUT, cwd=self.runtime_dir,
                                          start_new_session=True)
        finally:
            log.close()

    def _log_tail(self, n: int = 12) -> str:
        try:
            lines = self.log_path.read_text(errors="replace").splitlines()
        except (OSError, AttributeError):
            return ""
        keep = [l for l in lines if not l.startswith("Info)")] or lines
        return "\n".join(keep[-n:])

    def _fail(self, state: str, reason: str) -> None:
        tail = self._log_tail() if state == ViewerState.ERROR else ""
        self.state = state
        self.reason = reason + (f" (log: {self.log_path})" if tail else "")
        self._ready.set()

    def _send(self, data: bytes) -> None:
        conn = self._conn
        if conn is None:
            return
        with self._send_lock:
            conn.sendall(data)

    def _close_sockets(self) -> None:
        for s in (self._conn, self._listener):
            if s is not None:
                try:
                    s.shutdown(socket.SHUT_RDWR)          # wakes a reader blocked in recv
                except OSError:
                    pass
                try:
                    s.close()
                except OSError:
                    pass
        self._conn = self._listener = None

    def _accept(self) -> Optional[socket.socket]:
        deadline = time.monotonic() + self.startup_timeout
        self._listener.settimeout(0.25)
        while not self._done.is_set() and time.monotonic() < deadline:
            if self._proc is not None and self._proc.poll() is not None:
                self._fail(ViewerState.ERROR, f"VMD exited during startup (code "
                                              f"{self._proc.returncode})")
                return None
            try:
                conn, addr = self._listener.accept()
            except socket.timeout:
                continue
            except OSError:
                return None
            if addr[0] != "127.0.0.1":
                conn.close()
                continue
            return conn
        if not self._done.is_set():
            self._fail(ViewerState.ERROR, "VMD bridge did not connect in time")
        return None

    def _run(self) -> None:
        conn = self._accept()
        if conn is None:
            return
        self._listener.close()                          # exactly one bridge client
        self._listener = None
        self._conn = conn
        conn.settimeout(None)
        reader = conn.makefile("rb")
        try:
            self._serve(reader)
        except (OSError, ValueError):
            pass
        finally:
            if self.state in (ViewerState.READY, ViewerState.LOADING, ViewerState.STARTING) \
                    and not self._done.is_set():
                self._fail(ViewerState.DISCONNECTED, "VMD bridge connection closed"
                           + (f"; VMD exit code {self._proc.poll()}" if self._proc and
                              self._proc.poll() is not None else ""))
            self._ready.set()

    def _serve(self, reader) -> None:
        authed = False
        while not self._done.is_set():
            raw = reader.readline(proto.MAX_LINE + 1)
            if not raw:
                return
            if len(raw) > proto.MAX_LINE or not raw.endswith(b"\n"):
                self._protocol_error("oversized or unterminated line")
                return
            try:
                msg = proto.parse_bridge_line(raw)
            except proto.ProtocolError as exc:
                self._protocol_error(str(exc))
                return
            if not authed:
                if msg.kind != "HELLO" or not hmac.compare_digest(msg.fields[1], self.token):
                    self._protocol_error("bridge handshake failed (bad or missing token)")
                    return
                authed = True
                self._send(proto.welcome())
                continue
            if msg.kind == "LOADING":
                self.state = ViewerState.LOADING
            elif msg.kind == "READY":
                self._on_ready(*msg.fields)
            elif msg.kind == "FRAME":
                self._on_frame(*msg.fields)
            elif msg.kind == "ERROR":
                code, text = msg.fields
                if self.state != ViewerState.READY:
                    self._fail(ViewerState.ERROR, f"VMD bridge error {code}: {text}")
                    self._send(proto.bye())
                    return
                self.reason = f"last bridge error {code}: {text}"
            elif msg.kind == "BYE":
                return
            else:
                self._protocol_error(f"unexpected {msg.kind} after handshake")
                return

    def _protocol_error(self, why: str) -> None:
        self._fail(ViewerState.ERROR, f"bridge protocol violation: {why}")
        try:
            self._send(proto.bye())
        except OSError:
            pass

    def _on_ready(self, molid, numatoms, numframes, version) -> None:
        self.info = {"molid": molid, "numatoms": numatoms, "numframes": numframes,
                     "version": version}
        if numframes != self.n_frames:
            self._fail(ViewerState.ERROR, f"VMD loaded {numframes} frames but the display "
                                          f"timeline has {self.n_frames}; not synchronizing")
            self._send(proto.bye())
            return
        if self.n_atoms is not None and numatoms != self.n_atoms:
            self._fail(ViewerState.ERROR, f"VMD loaded {numatoms} atoms but the session's "
                                          f"System group has {self.n_atoms}; not synchronizing")
            self._send(proto.bye())
            return
        self.state = ViewerState.READY
        # the hub/session frame is authoritative on connection; a command issued
        # while VMD was loading is sent with its own sequence so the hub sees its ack
        if self.last_command is not None:
            frame, seq = self.last_command
        else:
            frame, seq = self._current_frame_provider(), proto.INITIAL_SYNC_SEQUENCE
        self._send(proto.goto(seq, frame))
        self._ready.set()

    def _on_frame(self, frame: int, ack: Optional[int]) -> None:
        if frame >= self.n_frames:
            self.reason = f"VMD reported frame {frame} outside the display timeline"
            return
        self.current_frame = frame
        if ack is not None:
            self.last_ack = ack
        if ack == proto.INITIAL_SYNC_SEQUENCE:
            return                                      # our own startup sync
        if self._callback is not None:
            self._callback(frame, ack_sequence=ack)
