"""Time-synchronisation metadata and frame ↔ time ↔ sample mapping.

Nothing is keyed by observable name: a result is *time-syncable* when it has
exactly one ``time`` axis whose samples can be related to the review
timeline.  Frame == row is never assumed — it is verified sample by sample
against the review trajectory's recorded frame times (Phase 1), or accepted
only from a ``one_row_per_frame`` alignment whose analysed-file digest *is*
the review trajectory.  Values are never interpolated or extrapolated.
"""
from __future__ import annotations

import bisect
from pathlib import Path
from typing import Optional, Sequence

from analysis.campaign.compatibility import same_time
from analysis.campaign.models import AxisKind, FrameAlignmentMode, ResultArray

MAPPING_POLICY = {
    "sample_to_time": "recorded time value of the sample (no interpolation)",
    "time_to_frame": "exact_timestamp",
    "cursor_to_frame": "nearest_frame (earlier frame on an exact tie; duplicates reported)",
    "interpolation": "none",
    "extrapolation": "none",
}


class SyncMapping:
    ROW_IS_FRAME = "row_i_is_frame_i"             # verified: sample i <-> frame i
    EXACT_TIMESTAMP = "exact_timestamp"           # every sample time is a frame time


def _time_values(axis, base: Optional[Path] = None) -> Optional[list[float]]:
    """Recorded time values; a relative ``values_ref`` resolves against ``base``
    (the session directory of a loaded dataset)."""
    from dataclasses import replace
    from analysis.campaign.results import read_column
    if axis.values is not None:
        return [float(v) for v in axis.values]
    if axis.values_ref is not None:
        ref = axis.values_ref
        if base is not None and not Path(ref.path).is_absolute():
            ref = replace(ref, path=str(Path(base) / ref.path))
        return read_column(ref)
    return None


def sync_metadata(array: ResultArray, timeline_times: Optional[Sequence[float]],
                  timeline_digest: Optional[str]) -> dict:
    """Whether (and how) ``array`` can follow the review trajectory cursor."""
    time_axes = [(i, a) for i, a in enumerate(array.axes) if a.kind == AxisKind.TIME]
    out = {"syncable": False, "time_axis": None, "time_axis_position": None,
           "alignment": None, "mapping": None, "mapping_policy": None,
           "n_samples": None, "reason": ""}
    if not time_axes:
        out["reason"] = f"no time axis (axes: {list(array.axis_signature()) or 'none'})"
        return out
    if len(time_axes) > 1:
        out["reason"] = "more than one time axis"
        return out
    pos, ax = time_axes[0]
    out.update(time_axis=ax.name, time_axis_position=pos,
               alignment=ax.alignment.to_dict() if ax.alignment else None)
    if ax.unit not in (None, "ps"):
        out["reason"] = f"time axis unit {ax.unit!r} is not ps"
        return out
    if timeline_times is None:
        out["reason"] = "review timeline unavailable"
        return out
    try:
        values = _time_values(ax)
    except (OSError, ValueError) as exc:
        out["reason"] = f"time axis unreadable: {exc}"
        return out
    frames = list(timeline_times)
    if values is None:
        al = ax.alignment
        if (al and al.mode == FrameAlignmentMode.ONE_ROW_PER_FRAME and timeline_digest
                and al.trajectory_digest == timeline_digest):
            out.update(syncable=True, mapping=SyncMapping.ROW_IS_FRAME,
                       mapping_policy=MAPPING_POLICY,
                       reason="one row per frame of the review trajectory itself (alignment digest)")
        else:
            out["reason"] = "time values not stored and the alignment does not name the review trajectory"
        return out
    out["n_samples"] = len(values)
    if len(values) == len(frames) and all(same_time(a, b) for a, b in zip(values, frames)):
        out.update(syncable=True, mapping=SyncMapping.ROW_IS_FRAME, mapping_policy=MAPPING_POLICY,
                   reason=f"all {len(values)} sample times equal the frame times in order")
        return out
    unmatched = [i for i, t in enumerate(values) if _frames_at(frames, t) == []]
    if unmatched:
        out["reason"] = (f"{len(unmatched)} of {len(values)} sample times match no frame "
                         f"(first: sample {unmatched[0]} at {values[unmatched[0]]} ps)")
        return out
    out.update(syncable=True, mapping=SyncMapping.EXACT_TIMESTAMP, mapping_policy=MAPPING_POLICY,
               reason="every sample time equals a recorded frame time")
    return out


def _frames_at(times: Sequence[float], t: float) -> list[int]:
    """Frames whose recorded time equals ``t`` (printed-precision tolerance).
    ``times`` must be non-decreasing; otherwise a linear scan is used."""
    lo = bisect.bisect_left(times, t - max(1e-3, 6e-6 * abs(t)))
    out = []
    for i in range(lo, len(times)):
        if times[i] > t + max(1e-3, 6e-6 * abs(t)):
            break
        if same_time(times[i], t):
            out.append(i)
    return out


# ═══════════════════════════════════════════════════════════════════════════════
# Mapping over a prepared session (used by later sync layers)
# ═══════════════════════════════════════════════════════════════════════════════

class ReviewTimeline:
    """source frame ↔ display frame ↔ simulation time ↔ result sample."""

    def __init__(self, source_times: Sequence[float], *,
                 display_times: Optional[Sequence[float]] = None,
                 base: Optional[Path] = None):
        from analysis.campaign.trajectory.time_map import FrameTimeMap
        self.base = Path(base) if base is not None else None
        self.source_times = [float(t) for t in source_times]
        self.display_times = ([float(t) for t in display_times]
                              if display_times is not None else None)
        self._source_map = FrameTimeMap(self.source_times)
        self._display_map = (FrameTimeMap(self.display_times)
                             if self.display_times is not None else self._source_map)

    @classmethod
    def from_dataset(cls, dataset, session_dir: str | Path) -> "ReviewTimeline":
        import numpy as np
        base = Path(session_dir)
        tl = dataset.timeline
        if not tl.get("times_ref"):
            raise ValueError("the session has no review timeline")
        src = np.load(base / tl["times_ref"]["path"]).tolist()
        disp = tl.get("display") or {}
        own = disp.get("times_ref")
        return cls(src, display_times=np.load(base / own["path"]).tolist() if own else None,
                   base=base)

    @property
    def n_frames(self) -> int:
        return len(self.display_times if self.display_times is not None else self.source_times)

    def display_frame_to_time(self, frame: int) -> float:
        return self._display_map.frame_to_time(frame)

    def display_frame_to_source_frames(self, frame: int) -> list[int]:
        """Identity when the display timeline is the source's; else by exact time."""
        if self.display_times is None:
            self._source_map.frame_to_time(frame)          # range check
            return [frame]
        return self._source_map.frames_at(self.display_frame_to_time(frame))

    def time_to_display_frame(self, time_ps: float, *, cursor: bool = False):
        """Exact lookup by default; ``cursor=True`` → nearest frame (Phase 1 ties)."""
        from analysis.campaign.trajectory.time_map import LookupStatus
        look = self._display_map.nearest_frame(time_ps)
        if not cursor and look.status != LookupStatus.EXACT:
            look.frame = None
        return look

    def sample_to_frames(self, result, sample: int) -> list[int]:
        """Display frames of one sample of a syncable ``ReviewResult``."""
        if not result.sync.get("syncable"):
            raise ValueError(f"result {result.array.name!r} is not time-syncable: "
                             f"{result.sync.get('reason')}")
        if result.sync.get("mapping") == SyncMapping.ROW_IS_FRAME and self.display_times is None:
            self._display_map.frame_to_time(sample)
            return [sample]
        return self._display_map.frames_at(self.sample_time(result, sample))

    def sample_time(self, result, sample: int) -> float:
        ax = result.array.axes[result.sync["time_axis_position"]]
        values = _time_values(ax, self.base)
        if values is None:
            raise ValueError("time values are not stored for this result")
        return values[sample]
