"""Pure frame ↔ simulation-time mapping over a validated per-frame timeline.

Semantics (all times in ps):

* ``frame_to_time(i)`` — the recorded timestamp of frame ``i``; no negative
  indexing, out-of-range raises ``IndexError``.
* ``nearest_frame(t)`` — the frame whose timestamp is closest to ``t``:

  - exact match → ``EXACT``; if several frames share that timestamp the
    earliest frame is selected and all of them are listed in
    ``duplicate_frames`` (the ambiguity is exposed, not hidden);
  - between two frames → ``NEAREST``; an exact tie selects the **earlier**
    frame and sets ``tie=True``;
  - before the first / after the last timestamp → ``BEFORE_START`` /
    ``AFTER_END`` with ``frame=None`` (no extrapolation).  With
    ``out_of_range="clamp"`` the boundary frame is returned and ``clamped``
    is set — never silently;
  - ``max_delta_ps`` rejects matches farther than the tolerance
    (``BEYOND_TOLERANCE``, ``frame=None``, nearest candidate still reported).

A timeline whose times decrease anywhere is refused
(:class:`TimelineAmbiguityError`): time → frame is not a function there.
Duplicate timestamps are allowed but always reported per lookup.
"""
from __future__ import annotations

import bisect
import math
from dataclasses import dataclass, field
from typing import Optional, Sequence

from analysis.campaign.models import FrameTimeIndex, ValidationState


class TimelineAmbiguityError(ValueError):
    """The timeline cannot support an unambiguous time → frame mapping."""


class LookupStatus:
    EXACT = "exact"
    NEAREST = "nearest"
    BEFORE_START = "before_start"
    AFTER_END = "after_end"
    BEYOND_TOLERANCE = "beyond_tolerance"


@dataclass
class FrameLookup:
    query_ps: float
    status: str                               # LookupStatus.*
    frame: Optional[int] = None               # selected frame; None when no frame is valid
    time_ps: Optional[float] = None           # actual timestamp of the selected / nearest frame
    delta_ps: Optional[float] = None          # |time_ps - query_ps|
    tie: bool = False                         # equidistant between two timestamps
    clamped: bool = False                     # out-of-range query clamped on request
    duplicate_frames: list[int] = field(default_factory=list)  # all frames sharing time_ps (if >1)

    @property
    def ambiguous(self) -> bool:
        return self.tie or len(self.duplicate_frames) > 1

    def to_dict(self) -> dict:
        return {
            "query_ps": self.query_ps, "status": self.status, "frame": self.frame,
            "time_ps": self.time_ps, "delta_ps": self.delta_ps, "tie": self.tie,
            "clamped": self.clamped, "duplicate_frames": self.duplicate_frames,
        }


class FrameTimeMap:
    """Frame ↔ time lookups over a non-decreasing timeline (file order)."""

    def __init__(self, times_ps: Sequence[float]):
        times = [float(t) for t in times_ps]
        if not times:
            raise ValueError("empty timeline")
        bad = [i for i, t in enumerate(times) if not math.isfinite(t)]
        if bad:
            raise ValueError(f"non-finite timestamp at frame(s) {bad[:5]}")
        dec = [i for i in range(len(times) - 1) if times[i + 1] < times[i]]
        if dec:
            i = dec[0]
            raise TimelineAmbiguityError(
                f"timeline is non-monotonic (frame {i} {times[i]} ps -> frame {i + 1} "
                f"{times[i + 1]} ps); time→frame mapping is ambiguous")
        self._t = times

    @classmethod
    def from_index(cls, idx: FrameTimeIndex) -> "FrameTimeMap":
        if idx.state == ValidationState.INVALID:
            codes = ", ".join(w.code for w in idx.warnings) or "invalid"
            raise TimelineAmbiguityError(f"time index is invalid ({codes})")
        return cls(idx.times_ps)

    @property
    def n_frames(self) -> int:
        return len(self._t)

    @property
    def start_ps(self) -> float:
        return self._t[0]

    @property
    def end_ps(self) -> float:
        return self._t[-1]

    def frame_to_time(self, frame: int) -> float:
        if isinstance(frame, bool) or not isinstance(frame, int):
            raise TypeError(f"frame must be an int, got {type(frame).__name__}")
        if frame < 0 or frame >= len(self._t):
            raise IndexError(f"frame {frame} outside 0..{len(self._t) - 1}")
        return self._t[frame]

    def frames_at(self, time_ps: float) -> list[int]:
        """All frames whose timestamp equals ``time_ps`` exactly."""
        lo = bisect.bisect_left(self._t, time_ps)
        hi = bisect.bisect_right(self._t, time_ps)
        return list(range(lo, hi))

    def nearest_frame(
        self,
        time_ps: float,
        *,
        out_of_range: str = "none",
        max_delta_ps: Optional[float] = None,
    ) -> FrameLookup:
        if out_of_range not in ("none", "clamp"):
            raise ValueError("out_of_range must be 'none' or 'clamp'")
        q = float(time_ps)
        if not math.isfinite(q):
            raise ValueError(f"non-finite query time {time_ps!r}")
        t = self._t

        if q < t[0] or q > t[-1]:
            before = q < t[0]
            status = LookupStatus.BEFORE_START if before else LookupStatus.AFTER_END
            edge_time = t[0] if before else t[-1]
            look = FrameLookup(query_ps=q, status=status, time_ps=edge_time,
                               delta_ps=abs(q - edge_time))
            if out_of_range == "clamp":
                frames = self.frames_at(edge_time)
                look.frame = frames[0]
                look.clamped = True
                look.duplicate_frames = frames if len(frames) > 1 else []
            return look

        i = bisect.bisect_left(t, q)
        if t[i] == q:
            chosen_time, status, tie = q, LookupStatus.EXACT, False
        else:
            left, right = t[i - 1], t[i]          # t[i-1] < q < t[i]; i >= 1 here
            dl, dr = q - left, right - q
            tie = dl == dr
            chosen_time = left if dl <= dr else right   # tie -> earlier frame
            status = LookupStatus.NEAREST
        frames = self.frames_at(chosen_time)
        look = FrameLookup(query_ps=q, status=status, frame=frames[0], time_ps=chosen_time,
                           delta_ps=abs(chosen_time - q), tie=tie,
                           duplicate_frames=frames if len(frames) > 1 else [])
        if max_delta_ps is not None and look.delta_ps > max_delta_ps:
            look.status = LookupStatus.BEYOND_TOLERANCE
            look.frame = None
        return look
