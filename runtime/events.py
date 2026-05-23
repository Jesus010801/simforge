from __future__ import annotations

import threading
import time
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Callable


class EventType(str, Enum):
    STEP_STARTED          = "STEP_STARTED"
    STEP_COMPLETED        = "STEP_COMPLETED"
    STEP_FAILED           = "STEP_FAILED"
    STEP_SKIPPED          = "STEP_SKIPPED"          # cache hit
    STDOUT                = "STDOUT"
    STDERR                = "STDERR"
    WARNING               = "WARNING"
    HEARTBEAT             = "HEARTBEAT"
    HANG_DETECTED         = "HANG_DETECTED"
    INTERACTIVE_PROMPT    = "INTERACTIVE_PROMPT"    # GROMACS "Select a group:"
    REMEDIATION_APPLIED   = "REMEDIATION_APPLIED"
    RETRY_STARTED         = "RETRY_STARTED"
    GPU_DETECTED          = "GPU_DETECTED"
    SCIENTIFIC_WARNING    = "SCIENTIFIC_WARNING"
    ARTIFACT_CREATED      = "ARTIFACT_CREATED"
    ARTIFACT_MODIFIED     = "ARTIFACT_MODIFIED"
    CACHE_HIT             = "CACHE_HIT"
    CACHE_MISS            = "CACHE_MISS"
    METRICS_SNAPSHOT      = "METRICS_SNAPSHOT"
    PERFORMANCE           = "PERFORMANCE"           # ns/day parsed from mdrun


class EventSeverity(str, Enum):
    DEBUG   = "DEBUG"
    INFO    = "INFO"
    WARNING = "WARNING"
    ERROR   = "ERROR"


@dataclass
class ExecutionEvent:
    event_type:   EventType
    workspace_id: str
    step_id:      str
    timestamp:    float         = field(default_factory=time.time)
    severity:     EventSeverity = EventSeverity.INFO
    message:      str           = ""
    data:         dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "event_type":   self.event_type.value,
            "workspace_id": self.workspace_id,
            "step_id":      self.step_id,
            "timestamp":    self.timestamp,
            "severity":     self.severity.value,
            "message":      self.message,
            "data":         self.data,
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> ExecutionEvent:
        return cls(
            event_type   = EventType(d["event_type"]),
            workspace_id = d["workspace_id"],
            step_id      = d["step_id"],
            timestamp    = d["timestamp"],
            severity     = EventSeverity(d["severity"]),
            message      = d.get("message", ""),
            data         = d.get("data", {}),
        )


Handler = Callable[[ExecutionEvent], None]


class EventBus:
    """Thread-safe pub/sub bus. Handlers run synchronously in the publisher's thread."""

    def __init__(self) -> None:
        self._lock:    threading.Lock              = threading.Lock()
        self._global:  list[Handler]               = []
        self._typed:   dict[EventType, list[Handler]] = {}

    def subscribe(self, handler: Handler, event_type: EventType | None = None) -> None:
        with self._lock:
            if event_type is None:
                self._global.append(handler)
            else:
                self._typed.setdefault(event_type, []).append(handler)

    def unsubscribe(self, handler: Handler, event_type: EventType | None = None) -> None:
        with self._lock:
            if event_type is None:
                self._global = [h for h in self._global if h is not handler]
            else:
                bucket = self._typed.get(event_type, [])
                self._typed[event_type] = [h for h in bucket if h is not handler]

    def publish(self, event: ExecutionEvent) -> None:
        with self._lock:
            handlers = list(self._global) + list(self._typed.get(event.event_type, []))
        for h in handlers:
            try:
                h(event)
            except Exception:
                pass  # never let a handler crash the runtime

    def make_publisher(self, workspace_id: str, step_id: str) -> BoundPublisher:
        return BoundPublisher(bus=self, workspace_id=workspace_id, step_id=step_id)


class BoundPublisher:
    """Pre-bound workspace_id/step_id so callers don't repeat them."""

    def __init__(self, bus: EventBus, workspace_id: str, step_id: str) -> None:
        self._bus          = bus
        self.workspace_id  = workspace_id
        self.step_id       = step_id

    def emit(
        self,
        event_type: EventType,
        message:    str = "",
        severity:   EventSeverity = EventSeverity.INFO,
        **data: Any,
    ) -> None:
        self._bus.publish(ExecutionEvent(
            event_type   = event_type,
            workspace_id = self.workspace_id,
            step_id      = self.step_id,
            message      = message,
            severity     = severity,
            data         = data,
        ))

    def rebind(self, step_id: str) -> BoundPublisher:
        return BoundPublisher(self._bus, self.workspace_id, step_id)
