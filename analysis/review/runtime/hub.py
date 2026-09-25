"""SyncHub — the in-memory synchronization core of a review session.

Simulation time is the cross-view coordinate; the display frame is where a
viewer is.  The hub knows only the session's :class:`ReviewTimeline` (Phase 9
frame/time semantics) — never GROMACS, never the scientific backend.

Loop suppression (deterministic, no timing heuristics)
------------------------------------------------------
Every accepted state change gets a monotonic hub ``sequence``.  When a change
did not come from the viewer, the hub sends ``goto_frame(frame, sequence)``
and remembers it as the *pending viewer command*.  A later viewer report is:

* an **acknowledgement** — it carries that sequence, or (for viewers that
  cannot echo sequences) reports exactly the pending frame — → consumed, no
  new event;
* a **no-op** — it reports the frame the hub already holds → no new event;
* an **independent move** otherwise → a new ``viewer`` event, which is
  never sent back to the viewer.

Event delivery is rate-limited per subscriber by :class:`Coalescer`; the hub's
own state is never downsampled.
"""
from __future__ import annotations

import threading
from dataclasses import dataclass, field
from typing import Callable, Optional


class Origin:
    BROWSER = "browser"
    VIEWER = "viewer"
    SERVER = "server"


class MappingStatus:
    FRAME = "frame"           # a frame was requested; its time is the recorded one
    EXACT = "exact"           # a time was requested and a frame has exactly that time
    NEAREST = "nearest"       # a time was requested; the nearest frame was chosen


class SyncError(ValueError):
    """A frame / time request the session timeline cannot honour."""


@dataclass
class SyncEvent:
    sequence: int
    origin: str
    frame: int
    time_ps: float
    mapping_status: str
    session_id: str
    requested_frame: Optional[int] = None
    requested_time_ps: Optional[float] = None
    delta_ps: Optional[float] = None
    tie: bool = False
    duplicate_frames: list[int] = field(default_factory=list)
    client_event_id: Optional[str] = None
    #: UI-delivery hint: explicit user actions are never coalesced away
    explicit: bool = False

    def to_dict(self) -> dict:
        return {"sequence": self.sequence, "origin": self.origin, "frame": self.frame,
                "time_ps": self.time_ps, "mapping_status": self.mapping_status,
                "session_id": self.session_id, "requested_frame": self.requested_frame,
                "requested_time_ps": self.requested_time_ps, "delta_ps": self.delta_ps,
                "tie": self.tie, "duplicate_frames": self.duplicate_frames,
                "client_event_id": self.client_event_id, "explicit": self.explicit}


class SyncHub:
    def __init__(self, timeline, session_id: str, *, viewer=None, initial_frame: int = 0):
        self.timeline = timeline
        self.session_id = session_id
        self._lock = threading.RLock()
        self._subscribers: list[Callable[[SyncEvent], None]] = []
        self._sequence = 0
        self._pending_viewer: Optional[tuple[int, int]] = None     # (sequence, frame)
        self.frame = initial_frame
        self.time_ps = timeline.display_frame_to_time(initial_frame)
        self.last_event: Optional[SyncEvent] = None
        self.suppressed = 0                                         # echoes / no-ops consumed
        self.viewer = None
        if viewer is not None:
            self.attach_viewer(viewer)

    # ── wiring ─────────────────────────────────────────────────────────────
    def attach_viewer(self, viewer) -> None:
        self.viewer = viewer
        viewer.on_frame_changed(self.viewer_frame_changed)

    def subscribe(self, callback: Callable[[SyncEvent], None]) -> Callable[[], None]:
        with self._lock:
            self._subscribers.append(callback)

        def unsubscribe():
            with self._lock:
                if callback in self._subscribers:
                    self._subscribers.remove(callback)
        return unsubscribe

    def state(self) -> dict:
        with self._lock:
            return {"session_id": self.session_id, "sequence": self._sequence,
                    "frame": self.frame, "time_ps": self.time_ps,
                    "last_event": self.last_event.to_dict() if self.last_event else None}

    # ── resolution (pure) ──────────────────────────────────────────────────
    def resolve_frame(self, frame) -> tuple[float, list[int]]:
        if isinstance(frame, bool) or not isinstance(frame, int):
            raise SyncError(f"frame must be an integer, got {frame!r}")
        try:
            t = self.timeline.display_frame_to_time(frame)
        except IndexError as exc:
            raise SyncError(f"frame {frame} outside 0..{self.timeline.n_frames - 1}") from exc
        dups = self.timeline.display_frames_at(t)
        return t, dups if len(dups) > 1 else []

    def resolve_time(self, time_ps) -> tuple[int, float, str, float, bool, list[int]]:
        from analysis.campaign.trajectory.time_map import LookupStatus
        try:
            q = float(time_ps)
        except (TypeError, ValueError) as exc:
            raise SyncError(f"time must be a number, got {time_ps!r}") from exc
        look = self.timeline.time_to_display_frame(q, cursor=True)
        if look.status in (LookupStatus.BEFORE_START, LookupStatus.AFTER_END) or look.frame is None:
            raise SyncError(f"time {q} ps is outside the trajectory "
                            f"({self.timeline.display_frame_to_time(0)}..."
                            f"{self.timeline.display_frame_to_time(self.timeline.n_frames - 1)} ps);"
                            f" no extrapolation")
        status = MappingStatus.EXACT if look.status == LookupStatus.EXACT else MappingStatus.NEAREST
        return look.frame, look.time_ps, status, look.delta_ps, look.tie, look.duplicate_frames

    # ── requests ───────────────────────────────────────────────────────────
    def set_frame(self, frame, *, origin: str = Origin.BROWSER,
                  client_event_id: Optional[str] = None) -> Optional[SyncEvent]:
        t, dups = self.resolve_frame(frame)
        return self._apply(origin, frame, t, MappingStatus.FRAME, requested_frame=frame,
                           duplicate_frames=dups, client_event_id=client_event_id)

    def set_time(self, time_ps, *, origin: str = Origin.BROWSER,
                 client_event_id: Optional[str] = None) -> Optional[SyncEvent]:
        frame, t, status, delta, tie, dups = self.resolve_time(time_ps)
        return self._apply(origin, frame, t, status, requested_time_ps=float(time_ps),
                           delta_ps=delta, tie=tie, duplicate_frames=dups,
                           client_event_id=client_event_id)

    def viewer_frame_changed(self, frame, ack_sequence: Optional[int] = None) -> Optional[SyncEvent]:
        """A viewer reports its current frame (``ack_sequence``: the hub command
        it is executing, when the viewer can say so)."""
        with self._lock:
            pending = self._pending_viewer
            if pending is not None and (ack_sequence == pending[0] or
                                        (ack_sequence is None and frame == pending[1])):
                self._pending_viewer = None
                self.suppressed += 1
                return None                                  # our own command, echoed back
            if ack_sequence is not None and ack_sequence < self._sequence:
                self.suppressed += 1
                return None                                  # acknowledges a superseded command
            t, dups = self.resolve_frame(frame)
            return self._apply(Origin.VIEWER, frame, t, MappingStatus.FRAME,
                               requested_frame=frame, duplicate_frames=dups)

    # ── core ───────────────────────────────────────────────────────────────
    def _apply(self, origin, frame, t, status, **kw) -> Optional[SyncEvent]:
        with self._lock:
            if frame == self.frame and origin == Origin.VIEWER:
                self.suppressed += 1
                return None                                  # viewer re-reports the current frame
            self._sequence += 1
            ev = SyncEvent(sequence=self._sequence, origin=origin, frame=frame, time_ps=t,
                           mapping_status=status, session_id=self.session_id,
                           explicit=origin == Origin.BROWSER, **kw)
            self.frame, self.time_ps, self.last_event = frame, t, ev
            if origin != Origin.VIEWER and self.viewer is not None:
                self._pending_viewer = (ev.sequence, frame)
            subscribers = list(self._subscribers)
            viewer = self.viewer if origin != Origin.VIEWER else None
        for cb in subscribers:
            cb(ev)
        if viewer is not None:
            viewer.goto_frame(frame, ev.sequence)
        return ev


# ═══════════════════════════════════════════════════════════════════════════════
# UI delivery rate control (not scientific decimation)
# ═══════════════════════════════════════════════════════════════════════════════

class Coalescer:
    """At most one delivery per ``min_interval`` seconds per subscriber.

    Between deliveries only the *latest* implicit event is kept (intermediate
    playback frames are skipped for the UI); explicit events (user actions)
    are always kept, in order.  The latest state is always delivered when the
    interval elapses.  The clock is injected so this is testable without sleeps.
    """

    def __init__(self, min_interval: float = 0.05):
        self.min_interval = min_interval
        self._last_sent: Optional[float] = None
        self._pending: list[SyncEvent] = []

    def offer(self, event: SyncEvent, now: float) -> list[SyncEvent]:
        self._pending = [e for e in self._pending if e.explicit] + [event]
        return self.due(now)

    def due(self, now: float) -> list[SyncEvent]:
        if not self._pending:
            return []
        if self._last_sent is not None and now - self._last_sent < self.min_interval:
            return []
        out, self._pending = self._pending, []
        self._last_sent = now
        return out

    def next_deadline(self) -> Optional[float]:
        if not self._pending:
            return None
        return (self._last_sent or 0.0) + self.min_interval

    @property
    def pending(self) -> int:
        return len(self._pending)
