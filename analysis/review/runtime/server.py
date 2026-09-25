"""Local review server — stdlib HTTP, SSE down, POST up.

Serves exactly one already prepared and *validated* ReviewDataset.  It never
discovers, prepares, analyses or preprocesses anything; runtime state (the
current frame) lives only in memory and nothing is written to the session.

Security boundary (local scientific viewer):

* binds to 127.0.0.1 by default, OS-chosen port;
* an ephemeral random token (runtime metadata only — never part of the
  session identity) is required for every API and the event stream, via the
  ``X-SimForge-Token`` header or ``?token=``;
* the Host header must name the loopback address the server listens on
  (DNS-rebinding guard);
* results are addressed by catalog id only; files are read solely from the
  dataset's authorised, re-hashed review-store entries;
* static assets come from a fixed whitelist; a CSP forbids any external
  resource.
"""
from __future__ import annotations

import hmac
import json
import re
import secrets
import threading
import time
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Optional
from urllib.parse import parse_qs, urlsplit

from analysis.review.runtime.hub import Coalescer, Origin, SyncError, SyncHub
from analysis.review.runtime.results import ResultCatalog, ResultUnavailable, session_payload
from analysis.review.runtime.viewer import NullViewerAdapter

STATIC_DIR = Path(__file__).resolve().parent.parent / "static"
_STATIC = {"/": ("index.html", "text/html; charset=utf-8"),
           "/index.html": ("index.html", "text/html; charset=utf-8"),
           "/static/app.js": ("app.js", "text/javascript; charset=utf-8"),
           "/static/app.css": ("app.css", "text/css; charset=utf-8")}
_RESULT_RE = re.compile(r"^/api/results/(r\d+-\d+)$")
_CSP = ("default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' data:; "
        "connect-src 'self'; frame-ancestors 'none'; base-uri 'none'; form-action 'none'")
MAX_BODY = 4096
HEARTBEAT_S = 10.0


class ReviewServeError(RuntimeError):
    """The dataset cannot be served (missing, invalid or stale)."""


class _Client:
    """One SSE subscriber with its own delivery rate control."""

    def __init__(self, min_interval: float, clock):
        self.coalescer = Coalescer(min_interval)
        self.clock = clock
        self.cond = threading.Condition()
        self.ready: list = []

    def push(self, event) -> None:
        with self.cond:
            self.ready.extend(self.coalescer.offer(event, self.clock()))
            self.cond.notify()

    def take(self, stopping: threading.Event, max_wait: float = 0.25) -> list:
        with self.cond:
            if not self.ready and not stopping.is_set():
                deadline = self.coalescer.next_deadline()
                wait = max_wait if deadline is None else max(0.0, min(max_wait,
                                                                     deadline - self.clock()))
                self.cond.wait(wait)
            self.ready.extend(self.coalescer.due(self.clock()))
            out, self.ready = self.ready, []
            return out

    def wake(self) -> None:
        with self.cond:
            self.cond.notify_all()


class _HTTPServer(ThreadingHTTPServer):
    daemon_threads = False                     # server_close() joins handler threads

    def handle_error(self, request, client_address):
        import sys
        if isinstance(sys.exc_info()[1], (BrokenPipeError, ConnectionResetError,
                                          ConnectionAbortedError)):
            return                             # the client went away; nothing to report
        super().handle_error(request, client_address)


class ReviewServer:
    def __init__(self, dataset_path: str | Path, *, host: str = "127.0.0.1", port: int = 0,
                 viewer=None, viewer_options: Optional[dict] = None,
                 min_interval: float = 0.05, clock=time.monotonic):
        from analysis.review.sync import ReviewTimeline
        from analysis.review.validate import _dataset_file, open_review_session
        f = _dataset_file(dataset_path)
        if not f.is_file():
            raise ReviewServeError(f"no review dataset at {dataset_path}")
        try:
            ds, val = open_review_session(f)
        except ValueError as exc:
            raise ReviewServeError(str(exc)) from exc
        if ds.status != "prepared":
            raise ReviewServeError(f"dataset status is {ds.status!r}; only prepared sessions are served")
        if not val.valid:
            raise ReviewServeError("the review session is stale or invalid and will not be "
                                   "served (re-prepare it with `simforge trajectory review`): "
                                   + "; ".join(val.problems))
        self.dataset, self.session_dir = ds, f.parent
        self.timeline = ReviewTimeline.from_dataset(ds, self.session_dir)
        self.catalog = ResultCatalog(ds, self.session_dir, self.timeline)
        if isinstance(viewer, str):
            from analysis.review.runtime.viewer import make_viewer
            viewer = make_viewer(viewer, dataset=ds, session_dir=self.session_dir,
                                 **(viewer_options or {}))
        self.viewer = viewer or NullViewerAdapter()
        self.hub = SyncHub(self.timeline, ds.session_id, viewer=self.viewer)
        self.host, self.port = host, port
        self.token = secrets.token_urlsafe(24)
        self.min_interval, self.clock = min_interval, clock
        self._stopping = threading.Event()
        self._clients: set[_Client] = set()
        self._clients_lock = threading.Lock()
        self._httpd: Optional[ThreadingHTTPServer] = None
        self._thread: Optional[threading.Thread] = None
        self._payload: Optional[dict] = None

    # ── lifecycle ──────────────────────────────────────────────────────────
    def start(self) -> "ReviewServer":
        handler = type("ReviewHandler", (_Handler,), {"review": self})
        self._httpd = _HTTPServer((self.host, self.port), handler)
        self.port = self._httpd.server_address[1]
        self._thread = threading.Thread(target=self._httpd.serve_forever,
                                        kwargs={"poll_interval": 0.1},
                                        name=f"simforge-review-{self.port}")
        self._thread.start()
        self.viewer.start({"session_id": self.dataset.session_id,
                           "n_frames": self.timeline.n_frames,
                           "current_frame": lambda: self.hub.frame})
        return self

    def stop(self) -> None:
        self._stopping.set()
        with self._clients_lock:
            for c in list(self._clients):
                c.wake()
        if self._httpd is not None:
            self._httpd.shutdown()
            self._httpd.server_close()
        if self._thread is not None:
            self._thread.join(timeout=5)
        self.viewer.stop()

    def __enter__(self):
        return self.start()

    def __exit__(self, *exc):
        self.stop()

    @property
    def base_url(self) -> str:
        return f"http://{self.host}:{self.port}"

    @property
    def url(self) -> str:
        return f"{self.base_url}/?token={self.token}"

    # ── data ───────────────────────────────────────────────────────────────
    def session(self) -> dict:
        if self._payload is None:
            self._payload = session_payload(self.dataset, self.timeline, self.catalog,
                                            self.viewer.status())
        return {**self._payload, "viewer": self.viewer.status(), "sync": self.hub.state()}

    def add_client(self) -> _Client:
        c = _Client(self.min_interval, self.clock)
        with self._clients_lock:
            self._clients.add(c)
        return c

    def remove_client(self, c: _Client) -> None:
        with self._clients_lock:
            self._clients.discard(c)


class _Handler(BaseHTTPRequestHandler):
    review: ReviewServer
    server_version = "SimForgeReview/1"

    def log_message(self, fmt, *args):          # operational noise stays off stdout
        pass

    # ── helpers ────────────────────────────────────────────────────────────
    def _send(self, status: int, body: bytes, ctype: str, extra: Optional[dict] = None):
        self.send_response(status)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Content-Security-Policy", _CSP)
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        self.wfile.write(body)

    def _json(self, status: int, obj) -> None:
        self._send(status, json.dumps(obj, allow_nan=False).encode(), "application/json")

    def _error(self, status: int, message: str) -> None:
        self._json(status, {"error": message})

    def _host_ok(self) -> bool:
        host = (self.headers.get("Host") or "").strip().lower()
        port = self.review.port
        allowed = {f"{self.review.host}:{port}", f"localhost:{port}", f"127.0.0.1:{port}"}
        return host in allowed

    def _authorized(self, query: dict) -> bool:
        given = self.headers.get("X-SimForge-Token") or (query.get("token") or [""])[0]
        return bool(given) and hmac.compare_digest(given, self.review.token)

    def _body(self) -> Optional[dict]:
        try:
            n = int(self.headers.get("Content-Length") or 0)
        except ValueError:
            return None
        if n <= 0 or n > MAX_BODY:
            return None
        try:
            obj = json.loads(self.rfile.read(n))
        except ValueError:
            return None
        return obj if isinstance(obj, dict) else None

    # ── routing ────────────────────────────────────────────────────────────
    def do_GET(self):
        url = urlsplit(self.path)
        query = parse_qs(url.query)
        if not self._host_ok():
            return self._error(HTTPStatus.FORBIDDEN, "unexpected Host header")
        if url.path in _STATIC:
            name, ctype = _STATIC[url.path]
            return self._send(HTTPStatus.OK, (STATIC_DIR / name).read_bytes(), ctype)
        if url.path == "/health":
            return self._json(HTTPStatus.OK, {"ok": True})
        if not self._authorized(query):
            return self._error(HTTPStatus.UNAUTHORIZED, "missing or invalid session token")
        if url.path == "/api/session":
            return self._json(HTTPStatus.OK, self.review.session())
        if url.path == "/api/sync/state":
            return self._json(HTTPStatus.OK, self.review.hub.state())
        if url.path == "/api/sync/events":
            return self._events()
        m = _RESULT_RE.match(url.path)
        if m:
            try:
                payload = self.review.catalog.data(m.group(1))
            except (ResultUnavailable, ValueError, OSError) as exc:
                return self._error(HTTPStatus.NOT_FOUND, f"result unavailable: {exc}")
            return self._json(HTTPStatus.OK, payload)       # a client hang-up is not "unavailable"
        return self._error(HTTPStatus.NOT_FOUND, "not found")

    def do_POST(self):
        url = urlsplit(self.path)
        if not self._host_ok():
            return self._error(HTTPStatus.FORBIDDEN, "unexpected Host header")
        if not self._authorized(parse_qs(url.query)):
            return self._error(HTTPStatus.UNAUTHORIZED, "missing or invalid session token")
        if url.path not in ("/api/sync/frame", "/api/sync/time"):
            return self._error(HTTPStatus.NOT_FOUND, "not found")
        body = self._body()
        if body is None:
            return self._error(HTTPStatus.BAD_REQUEST, "expected a small JSON object")
        cid = body.get("client_event_id")
        cid = str(cid)[:64] if cid is not None else None
        hub = self.review.hub
        try:
            if url.path == "/api/sync/frame":
                ev = hub.set_frame(body.get("frame"), origin=Origin.BROWSER, client_event_id=cid)
            else:
                ev = hub.set_time(body.get("time_ps"), origin=Origin.BROWSER, client_event_id=cid)
        except SyncError as exc:
            return self._error(HTTPStatus.BAD_REQUEST, str(exc))
        # the hub state changed; whether a viewer applied it is reported separately
        return self._json(HTTPStatus.OK, {"event": ev.to_dict() if ev else None,
                                          "state": hub.state(),
                                          "viewer": self.review.viewer.status()})

    def do_PUT(self):
        self._error(HTTPStatus.METHOD_NOT_ALLOWED, "method not allowed")

    do_DELETE = do_PATCH = do_PUT

    # ── server-sent events ─────────────────────────────────────────────────
    def _events(self):
        review = self.review
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", "text/event-stream")
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.end_headers()
        client = review.add_client()
        unsubscribe = review.hub.subscribe(client.push)
        last_beat = time.monotonic()
        try:
            self._write_event("state", review.hub.state())
            while not review._stopping.is_set():
                for ev in client.take(review._stopping):
                    self._write_event("sync", ev.to_dict(), ev.sequence)
                if time.monotonic() - last_beat > HEARTBEAT_S:
                    self.wfile.write(b": keep-alive\n\n")
                    self.wfile.flush()
                    last_beat = time.monotonic()
        except (BrokenPipeError, ConnectionResetError, ConnectionAbortedError, OSError):
            pass
        finally:
            unsubscribe()
            review.remove_client(client)
            self.close_connection = True

    def _write_event(self, name: str, data: dict, eid: Optional[int] = None) -> None:
        msg = (f"id: {eid}\n" if eid is not None else "") + f"event: {name}\n" + \
              f"data: {json.dumps(data, allow_nan=False)}\n\n"
        self.wfile.write(msg.encode())
        self.wfile.flush()


def serve_forever(server: ReviewServer, *, open_browser: bool = True, echo=print) -> None:
    """Blocking helper for the CLI: start, print the URL, optionally open a
    browser (convenience only — failure is reported, never fatal), wait for
    Ctrl-C, stop cleanly."""
    server.start()
    echo(f"Review session {server.dataset.session_id} served at:\n  {server.url}")
    echo("Local only (127.0.0.1); the token in the URL is required. Ctrl-C to stop.")
    if open_browser:
        try:
            import webbrowser
            if not webbrowser.open(server.url):
                echo("(no browser could be opened; open the URL manually)")
        except Exception as exc:  # noqa: BLE001 — convenience only
            echo(f"(browser launch failed: {exc}; open the URL manually)")
    last = None
    try:
        while True:
            st = server.viewer.status()
            key = (st.get("state"), st.get("reason"))
            if st.get("state") and key != last:          # viewer lifecycle is reported, the
                last = key                                # dashboard keeps running regardless
                detail = (f" ({st.get('numframes')} frames, {st.get('numatoms')} atoms, "
                          f"VMD {st.get('version')})" if st.get("state") == "ready" else "")
                echo(f"viewer {st.get('kind')}: {st.get('state')}{detail}"
                     + (f" — {st['reason']}" if st.get("reason") else ""))
            time.sleep(0.5)
    except KeyboardInterrupt:
        echo("stopping…")
    finally:
        server.stop()
