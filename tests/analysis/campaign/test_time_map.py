"""Pure frame ↔ time mapping (analysis/campaign/trajectory/time_map.py)."""
from __future__ import annotations

import math

import pytest

from analysis.campaign.models import CampaignWarning, FrameTimeIndex, Severity, ValidationState
from analysis.campaign.trajectory.time_map import (
    FrameTimeMap, LookupStatus, TimelineAmbiguityError,
)


def test_frame_to_time_uses_recorded_times_not_dt():
    m = FrameTimeMap([50000.0, 50020.0, 50040.0, 50100.0])
    assert m.frame_to_time(0) == 50000.0
    assert m.frame_to_time(3) == 50100.0          # irregular: not 50060
    assert (m.start_ps, m.end_ps, m.n_frames) == (50000.0, 50100.0, 4)


@pytest.mark.parametrize("bad", [-1, 4, 100])
def test_frame_to_time_out_of_range(bad):
    with pytest.raises(IndexError):
        FrameTimeMap([0.0, 20.0, 40.0, 60.0]).frame_to_time(bad)


def test_frame_to_time_rejects_non_int():
    m = FrameTimeMap([0.0, 20.0])
    with pytest.raises(TypeError):
        m.frame_to_time(1.0)
    with pytest.raises(TypeError):
        m.frame_to_time(True)


def test_exact_lookup():
    m = FrameTimeMap([0.0, 20.0, 40.0, 60.0])
    look = m.nearest_frame(40.0)
    assert look.status == LookupStatus.EXACT
    assert (look.frame, look.time_ps, look.delta_ps) == (2, 40.0, 0.0)
    assert not look.ambiguous


def test_between_frames_nearest():
    m = FrameTimeMap([0.0, 20.0, 40.0, 60.0])
    a = m.nearest_frame(27.0)
    assert (a.status, a.frame, a.time_ps, a.delta_ps, a.tie) == (LookupStatus.NEAREST, 1, 20.0, 7.0, False)
    b = m.nearest_frame(33.0)
    assert (b.frame, b.time_ps) == (2, 40.0)


def test_tie_prefers_earlier_frame():
    look = FrameTimeMap([0.0, 20.0, 40.0, 60.0]).nearest_frame(30.0)
    assert look.frame == 1 and look.time_ps == 20.0
    assert look.tie is True and look.ambiguous


def test_irregular_timeline_lookup():
    m = FrameTimeMap([0.0, 20.0, 40.0, 100.0, 120.0])
    assert m.nearest_frame(69.0).frame == 2
    assert m.nearest_frame(71.0).frame == 3
    tie = m.nearest_frame(70.0)
    assert tie.frame == 2 and tie.tie


def test_out_of_range_is_explicit():
    m = FrameTimeMap([50000.0, 50020.0])
    before = m.nearest_frame(0.0)
    assert before.status == LookupStatus.BEFORE_START and before.frame is None
    assert before.time_ps == 50000.0 and before.delta_ps == 50000.0
    after = m.nearest_frame(60000.0)
    assert after.status == LookupStatus.AFTER_END and after.frame is None


def test_out_of_range_clamp_is_flagged():
    m = FrameTimeMap([10.0, 20.0])
    look = m.nearest_frame(5.0, out_of_range="clamp")
    assert look.frame == 0 and look.clamped and look.status == LookupStatus.BEFORE_START
    look = m.nearest_frame(25.0, out_of_range="clamp")
    assert look.frame == 1 and look.clamped and look.status == LookupStatus.AFTER_END


def test_duplicate_timestamps_exposed():
    m = FrameTimeMap([0.0, 20.0, 40.0, 40.0, 60.0])
    exact = m.nearest_frame(40.0)
    assert exact.frame == 2 and exact.duplicate_frames == [2, 3] and exact.ambiguous
    near = m.nearest_frame(45.0)
    assert near.frame == 2 and near.duplicate_frames == [2, 3]
    # the frame *after* the duplicates maps normally
    assert m.nearest_frame(59.0).frame == 4
    assert m.frames_at(40.0) == [2, 3]


def test_non_monotonic_timeline_refused():
    with pytest.raises(TimelineAmbiguityError):
        FrameTimeMap([0.0, 20.0, 40.0, 10.0, 30.0])


def test_invalid_timelines():
    with pytest.raises(ValueError):
        FrameTimeMap([])
    with pytest.raises(ValueError):
        FrameTimeMap([0.0, math.nan])
    with pytest.raises(ValueError):
        FrameTimeMap([0.0, 1.0]).nearest_frame(math.inf)
    with pytest.raises(ValueError):
        FrameTimeMap([0.0, 1.0]).nearest_frame(0.5, out_of_range="extrapolate")


def test_tolerance():
    m = FrameTimeMap([0.0, 100.0])
    look = m.nearest_frame(30.0, max_delta_ps=10.0)
    assert look.status == LookupStatus.BEYOND_TOLERANCE and look.frame is None
    assert look.time_ps == 0.0 and look.delta_ps == 30.0
    assert m.nearest_frame(5.0, max_delta_ps=10.0).frame == 0


def test_from_index_refuses_invalid_index():
    idx = FrameTimeIndex(trajectory_path="x.xtc", fingerprint=None, backend="b",
                         backend_version="1", times_ps=[0.0, 1.0],
                         state=ValidationState.INVALID,
                         warnings=[CampaignWarning("frame_count_mismatch", "m", Severity.ERROR)])
    with pytest.raises(TimelineAmbiguityError):
        FrameTimeMap.from_index(idx)
    idx.state = ValidationState.WARNING
    assert FrameTimeMap.from_index(idx).n_frames == 2
