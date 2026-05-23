from __future__ import annotations

import json
import threading
from pathlib import Path
from typing import TYPE_CHECKING

from runtime.events import EventBus, EventType, ExecutionEvent

if TYPE_CHECKING:
    pass


class JournalWriter:
    """
    Append-only JSONL journal — one ExecutionEvent per line.

    Writes to <workspace_dir>/metadata/execution_journal.jsonl.
    Thread-safe: a single lock guards all writes.
    """

    # Events too noisy for persistent storage (still dispatched on the bus)
    _SKIP_TYPES: frozenset[EventType] = frozenset({
        EventType.HEARTBEAT,
        EventType.STDOUT,
    })

    def __init__(self, workspace_dir: Path) -> None:
        journal_dir = workspace_dir / "metadata"
        journal_dir.mkdir(parents=True, exist_ok=True)
        self._path = journal_dir / "execution_journal.jsonl"
        self._lock = threading.Lock()

    # ── EventBus integration ──────────────────────────────────────────────────

    def register(self, bus: EventBus) -> None:
        """Subscribe to all events on the given bus."""
        bus.subscribe(self._handle)

    def _handle(self, event: ExecutionEvent) -> None:
        if event.event_type in self._SKIP_TYPES:
            return
        line = json.dumps(event.to_dict(), ensure_ascii=False) + "\n"
        with self._lock:
            with self._path.open("a", encoding="utf-8") as fh:
                fh.write(line)

    # ── Reading ───────────────────────────────────────────────────────────────

    def read_all(self) -> list[ExecutionEvent]:
        if not self._path.exists():
            return []
        events: list[ExecutionEvent] = []
        with self._lock:
            with self._path.open("r", encoding="utf-8") as fh:
                for raw in fh:
                    raw = raw.strip()
                    if raw:
                        try:
                            events.append(ExecutionEvent.from_dict(json.loads(raw)))
                        except Exception:
                            pass
        return events

    def path(self) -> Path:
        return self._path
