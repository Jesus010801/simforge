"""Viewer boundary.  Phase 10 ships only the null and fake adapters; a real
molecular viewer (VMD, Phase 11) implements the same interface.

Contract for implementations:

* ``goto_frame(frame, event_sequence)`` moves the viewer; when the viewer then
  reports the frame it must pass ``ack_sequence=event_sequence`` if it can —
  this is what lets the hub suppress the echo deterministically;
* user-driven moves are reported with ``ack_sequence=None``;
* ``acknowledges = True`` declares that the viewer reports in order and acks
  every command it receives (the hub then discards un-acked reports that
  were in flight before a command);
* ``stop()`` releases every thread / process the adapter owns;
* ``start(session)`` may receive ``session["current_frame"]`` (a callable):
  the hub's frame, authoritative when a viewer connects.
"""
from __future__ import annotations

import threading
from typing import Callable, Optional, Protocol, runtime_checkable

FrameCallback = Callable[..., object]          # (frame, ack_sequence=None)


@runtime_checkable
class ViewerAdapter(Protocol):
    def start(self, session: dict) -> None: ...
    def stop(self) -> None: ...
    def goto_frame(self, frame: int, event_sequence: int) -> None: ...
    def on_frame_changed(self, callback: FrameCallback) -> None: ...
    def status(self) -> dict: ...


class NullViewerAdapter:
    """No viewer: the dashboard alone drives the session."""

    kind = "null"
    acknowledges = False

    def __init__(self):
        self._callback: Optional[FrameCallback] = None
        self.commands: list[tuple[int, int]] = []

    def start(self, session: dict) -> None:
        pass

    def stop(self) -> None:
        pass

    def goto_frame(self, frame: int, event_sequence: int) -> None:
        self.commands.append((frame, event_sequence))   # nothing to move

    def on_frame_changed(self, callback: FrameCallback) -> None:
        self._callback = callback

    def status(self) -> dict:
        return {"kind": self.kind, "connected": False,
                "note": "no molecular viewer attached (Phase 10)"}


class FakeViewerAdapter:
    """Test / development stand-in for a molecular viewer.

    ``echo=True``: every ``goto_frame`` is reported back like a real viewer
    would (with the acknowledged sequence unless ``ack=False``, which imitates
    a viewer that cannot echo sequences).  ``user_moves_to(frame)`` imitates
    the user moving the viewer.
    """

    kind = "fake"

    def __init__(self, *, echo: bool = True, ack: bool = True):
        self.echo, self.ack = echo, ack
        self.acknowledges = echo and ack
        self._callback: Optional[FrameCallback] = None
        self._lock = threading.Lock()
        self.frame: Optional[int] = None
        self.commands: list[tuple[int, int]] = []
        self.running = False

    def start(self, session: dict) -> None:
        self.running = True

    def stop(self) -> None:
        self.running = False

    def goto_frame(self, frame: int, event_sequence: int) -> None:
        with self._lock:
            self.commands.append((frame, event_sequence))
            self.frame = frame
        if self.echo and self._callback is not None:
            self._callback(frame, ack_sequence=event_sequence if self.ack else None)

    def on_frame_changed(self, callback: FrameCallback) -> None:
        self._callback = callback

    def user_moves_to(self, frame: int):
        with self._lock:
            self.frame = frame
        return self._callback(frame, ack_sequence=None) if self._callback else None

    def status(self) -> dict:
        return {"kind": self.kind, "connected": self.running, "frame": self.frame,
                "commands": len(self.commands), "note": "fake viewer (tests / development)"}


VIEWERS = ("null", "fake", "vmd")


def make_viewer(kind: str, *, dataset=None, session_dir=None, **options):
    """``null`` / ``fake`` need nothing; ``vmd`` is built from the validated
    dataset's display view (``options``: vmd, headless, startup_timeout)."""
    if kind == "null":
        return NullViewerAdapter()
    if kind == "fake":
        return FakeViewerAdapter()
    if kind == "vmd":
        from analysis.review.runtime.vmd import VMDViewerAdapter
        if dataset is None:
            raise ValueError("the vmd viewer needs the served dataset")
        return VMDViewerAdapter.for_dataset(dataset, session_dir, **options)
    raise ValueError(f"unknown viewer {kind!r} (available: {', '.join(VIEWERS)})")
