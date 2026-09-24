"""Per-frame trajectory time index (analysis/campaign/trajectory/time_index.py).

Two layers:

* backend-neutral validation / cache / segment logic, driven by a fake backend
  (no GROMACS needed);
* gmx-guarded tests that build tiny XTC/TRR files from multi-frame GRO input
  with chosen timestamps and index them with the real GROMACS backend.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from analysis.campaign.gmx import run_gmx
from analysis.campaign.models import (
    FrameTimeIndex, SegmentRelationKind, TrajectoryRequirements, TrajectoryView,
    ValidationState, ViewKind,
)
from analysis.campaign.trajectory import time_index as ti
from analysis.campaign.trajectory.inspector import _parse_gmx_check
from analysis.campaign.trajectory.time_map import FrameTimeMap
from tests.analysis.campaign.conftest import REAL_XTC, requires_gmx, requires_real_traj


# ═══════════════════════════════════════════════════════════════════════════════
# Fake backend
# ═══════════════════════════════════════════════════════════════════════════════

class FakeBackend:
    def __init__(self, times, *, id="fake/v1", version="1.0", truncated=False,
                 reported=None, ok=True):
        self.id = id
        self._version = version
        self.times = list(times)
        self.truncated = truncated
        self.reported = reported
        self.ok = ok
        self.scans = 0

    def available(self):
        return True

    def version(self):
        return self._version

    def scan(self, path):
        self.scans += 1
        idx = FrameTimeIndex(trajectory_path=str(path), fingerprint=None, backend=self.id,
                             backend_version=self._version)
        if not self.ok:
            idx.state = ValidationState.INVALID
            return idx, False
        idx.times_ps = list(self.times)
        idx.steps = list(range(len(self.times)))
        idx.boxes = [[3.0, 3.0, 3.0, 0, 0, 0, 0, 0, 0] for _ in self.times]
        idx.truncated = self.truncated
        idx.truncation_evidence = "WARNING: Incomplete frame: nr 9 time 9" if self.truncated else ""
        idx.backend_reported_frames = len(self.times) if self.reported is None else self.reported
        return idx, True


def _traj(tmp_path, name="md.xtc", payload=b"\x00" * 64) -> Path:
    p = tmp_path / name
    p.write_bytes(payload)
    return p


def _codes(idx):
    return {w.code for w in idx.warnings}


# ═══════════════════════════════════════════════════════════════════════════════
# Validation (fake backend)
# ═══════════════════════════════════════════════════════════════════════════════

def test_regular_timeline_valid(tmp_path):
    idx, cached = ti.index_trajectory(_traj(tmp_path), backend=FakeBackend([0, 20, 40, 60]))
    assert not cached
    assert idx.state == ValidationState.VALID
    assert idx.times_ps == [0, 20, 40, 60] and idx.n_frames == 4
    assert idx.strictly_increasing and idx.non_decreasing
    assert idx.duplicate_groups == [] and idx.non_monotonic_steps == []
    assert idx.fingerprint is not None and idx.backend == "fake/v1"


def test_non_zero_start_preserved(tmp_path):
    idx, _ = ti.index_trajectory(_traj(tmp_path), backend=FakeBackend([50000, 50020, 50040]))
    assert idx.start_time_ps == 50000 and idx.end_time_ps == 50040
    assert idx.state == ValidationState.VALID


def test_irregular_timeline_is_valid_not_resampled(tmp_path):
    idx, _ = ti.index_trajectory(_traj(tmp_path), backend=FakeBackend([0, 20, 40, 100, 120]))
    assert idx.times_ps == [0, 20, 40, 100, 120]
    assert idx.state == ValidationState.VALID


def test_duplicates_detected_not_removed(tmp_path):
    idx, _ = ti.index_trajectory(_traj(tmp_path), backend=FakeBackend([0, 20, 40, 40, 60]))
    assert idx.n_frames == 5 and idx.times_ps == [0, 20, 40, 40, 60]
    assert idx.duplicate_groups == [[2, 3]]
    assert idx.non_decreasing and not idx.strictly_increasing
    assert idx.state == ValidationState.WARNING
    assert "duplicate_timestamps" in _codes(idx)
    # distinct steps on identical times are flagged (precision or re-timing)
    assert "duplicate_time_distinct_step" in _codes(idx)


def test_non_monotonic_detected_not_reordered(tmp_path):
    idx, _ = ti.index_trajectory(_traj(tmp_path), backend=FakeBackend([0, 20, 40, 10, 30]))
    assert idx.times_ps == [0, 20, 40, 10, 30]
    assert idx.non_monotonic_steps == [[2, 3]]
    assert idx.non_decreasing is False and idx.strictly_increasing is False
    assert idx.state == ValidationState.WARNING
    assert "non_monotonic_timestamps" in _codes(idx)


def test_empty_trajectory_invalid(tmp_path):
    idx, _ = ti.index_trajectory(_traj(tmp_path), backend=FakeBackend([]))
    assert idx.state == ValidationState.INVALID and "empty_trajectory" in _codes(idx)


def test_non_finite_invalid(tmp_path):
    idx, _ = ti.index_trajectory(_traj(tmp_path), backend=FakeBackend([0.0, float("nan"), 40.0]))
    assert idx.state == ValidationState.INVALID
    assert idx.non_finite_frames == [1] and "non_finite_timestamps" in _codes(idx)


def test_backend_frame_count_mismatch_invalid(tmp_path):
    idx, _ = ti.index_trajectory(_traj(tmp_path), backend=FakeBackend([0, 20, 40], reported=4))
    assert idx.state == ValidationState.INVALID and "frame_count_mismatch" in _codes(idx)


def test_expected_frame_count_mismatch_reported(tmp_path):
    idx, _ = ti.index_trajectory(_traj(tmp_path), backend=FakeBackend([0, 20, 40]),
                                 expected_n_frames=5)
    assert idx.state == ValidationState.WARNING and "unexpected_frame_count" in _codes(idx)
    assert idx.expected_n_frames == 5


def test_truncation_reported(tmp_path):
    idx, _ = ti.index_trajectory(_traj(tmp_path), backend=FakeBackend([0, 20], truncated=True))
    assert idx.truncated and idx.state == ValidationState.WARNING
    assert "truncated_trajectory" in _codes(idx)


def test_backend_failure_not_cached(tmp_path):
    be = FakeBackend([0, 20], ok=False)
    cache = tmp_path / "cache"
    idx, _ = ti.index_trajectory(_traj(tmp_path), backend=be, cache_dir=cache)
    assert idx.state == ValidationState.INVALID
    assert not cache.exists() or not any(cache.rglob("index.json"))


def test_unsupported_format_and_missing_file(tmp_path):
    with pytest.raises(ValueError):
        ti.index_trajectory(_traj(tmp_path, "md.gro"), backend=FakeBackend([0]))
    with pytest.raises(FileNotFoundError):
        ti.index_trajectory(tmp_path / "missing.xtc", backend=FakeBackend([0]))


def test_unavailable_backend(tmp_path):
    idx, cached = ti.index_trajectory(_traj(tmp_path),
                                      backend=ti.GromacsTimeIndexBackend("/nonexistent/gmx"))
    assert idx.state == ValidationState.INVALID and not cached
    assert "time_index_backend_unavailable" in _codes(idx)


# ═══════════════════════════════════════════════════════════════════════════════
# Cache
# ═══════════════════════════════════════════════════════════════════════════════

def test_cache_hit_reuses_index(tmp_path):
    traj, cache = _traj(tmp_path), tmp_path / "cache"
    be = FakeBackend([50000, 50020, 50040, 50040])
    first, c1 = ti.index_trajectory(traj, backend=be, cache_dir=cache)
    second, c2 = ti.index_trajectory(traj, backend=be, cache_dir=cache)
    assert (c1, c2) == (False, True) and be.scans == 1
    assert second.times_ps == first.times_ps and second.steps == first.steps
    assert second.boxes == first.boxes
    assert second.duplicate_groups == [[2, 3]] and second.state == ValidationState.WARNING
    assert second.fingerprint.digest == first.fingerprint.digest


def test_cache_persists_backend_metadata(tmp_path):
    traj, cache = _traj(tmp_path), tmp_path / "cache"
    ti.index_trajectory(traj, backend=FakeBackend([0, 20], version="9.9"), cache_dir=cache)
    (meta_file,) = list(cache.rglob("index.json"))
    meta = json.loads(meta_file.read_text())
    assert meta["backend"] == "fake/v1" and meta["backend_version"] == "9.9"
    assert meta["schema_version"] == ti.TIME_INDEX_SCHEMA_VERSION
    assert meta["fingerprint"]["digest"]
    assert "times_ps" not in meta                       # arrays live in frames.npz
    assert (meta_file.parent / "frames.npz").is_file()


def test_cache_invalidated_by_content_change(tmp_path):
    traj, cache = _traj(tmp_path), tmp_path / "cache"
    be = FakeBackend([0, 20])
    ti.index_trajectory(traj, backend=be, cache_dir=cache)
    traj.write_bytes(b"\x01" * 64)                       # same size, new content
    os.utime(traj, ns=(1, 1))
    _, cached = ti.index_trajectory(traj, backend=be, cache_dir=cache)
    assert not cached and be.scans == 2


def test_cache_invalidated_by_backend_version_and_id(tmp_path):
    traj, cache = _traj(tmp_path), tmp_path / "cache"
    ti.index_trajectory(traj, backend=FakeBackend([0, 20], version="1"), cache_dir=cache)
    be2 = FakeBackend([0, 20], version="2")
    _, cached = ti.index_trajectory(traj, backend=be2, cache_dir=cache)
    assert not cached and be2.scans == 1
    be3 = FakeBackend([0, 20], id="other/v1", version="2")
    _, cached = ti.index_trajectory(traj, backend=be3, cache_dir=cache)
    assert not cached and be3.scans == 1


def test_cache_rejects_stale_schema(tmp_path):
    traj, cache = _traj(tmp_path), tmp_path / "cache"
    be = FakeBackend([0, 20])
    ti.index_trajectory(traj, backend=be, cache_dir=cache)
    (meta_file,) = list(cache.rglob("index.json"))
    meta = json.loads(meta_file.read_text())
    meta["schema_version"] = "simforge/trajectory-time-index/v0"
    meta_file.write_text(json.dumps(meta))
    _, cached = ti.index_trajectory(traj, backend=be, cache_dir=cache)
    assert not cached and be.scans == 2


def test_force_rescans(tmp_path):
    traj, cache = _traj(tmp_path), tmp_path / "cache"
    be = FakeBackend([0, 20])
    ti.index_trajectory(traj, backend=be, cache_dir=cache)
    _, cached = ti.index_trajectory(traj, backend=be, cache_dir=cache, force=True)
    assert not cached and be.scans == 2


def test_cached_expected_frames_hint_does_not_poison_cache(tmp_path):
    traj, cache = _traj(tmp_path), tmp_path / "cache"
    be = FakeBackend([0, 20])
    idx, _ = ti.index_trajectory(traj, backend=be, cache_dir=cache, expected_n_frames=3)
    assert "unexpected_frame_count" in _codes(idx)
    again, cached = ti.index_trajectory(traj, backend=be, cache_dir=cache)
    assert cached and again.state == ValidationState.VALID
    assert "unexpected_frame_count" not in _codes(again)


def test_no_cache_dir_writes_nothing(tmp_path):
    traj = _traj(tmp_path)
    before = set(tmp_path.rglob("*"))
    ti.index_trajectory(traj, backend=FakeBackend([0, 20]))
    assert set(tmp_path.rglob("*")) == before


def test_serialization_roundtrip():
    idx = FrameTimeIndex(trajectory_path="/x/md.xtc", fingerprint=None, backend="b",
                         backend_version="1", times_ps=[1.0, 2.0], steps=[10, 20],
                         boxes=[[1.0] * 9, [2.0] * 9])
    d = idx.to_dict(include_frames=True)
    back = FrameTimeIndex.from_dict(json.loads(json.dumps(d)))
    assert back.same_timeline(idx) and back.boxes == idx.boxes
    assert "times_ps" not in idx.to_dict()


# ═══════════════════════════════════════════════════════════════════════════════
# Segments and views
# ═══════════════════════════════════════════════════════════════════════════════

def _idx(path, times):
    return FrameTimeIndex(trajectory_path=path, fingerprint=None, backend="b",
                          backend_version="1", times_ps=list(times))


def test_segment_relations_descriptive_only():
    a = _idx("a.xtc", [0, 20, 40])
    b = _idx("b.xtc", [40, 60, 80])       # gmx continuation sharing a boundary frame
    c = _idx("c.xtc", [100, 120])
    d = _idx("d.xtc", [30, 50])           # overlaps a and b
    rels = {(r.first_path, r.second_path): r for r in ti.relate_segments([a, b, c, d])}
    assert rels[("a.xtc", "b.xtc")].relation == SegmentRelationKind.BOUNDARY_SHARED
    assert rels[("a.xtc", "b.xtc")].shared_timestamps == 1
    assert rels[("b.xtc", "c.xtc")].relation == SegmentRelationKind.SEQUENTIAL
    assert rels[("b.xtc", "c.xtc")].gap_ps == 20
    ad = rels[("a.xtc", "d.xtc")]
    assert ad.relation == SegmentRelationKind.OVERLAP and ad.overlap_ps == 10
    assert rels[("c.xtc", "d.xtc")].relation == SegmentRelationKind.REVERSED
    assert len(rels) == 6


def test_segment_relation_unknown_for_empty():
    (r,) = ti.relate_segments([_idx("a.xtc", []), _idx("b.xtc", [1.0])])
    assert r.relation == SegmentRelationKind.UNKNOWN


def test_index_segments_independent(tmp_path):
    p1, p2 = _traj(tmp_path, "part1.xtc"), _traj(tmp_path, "part2.xtc", b"\x02" * 64)
    out = ti.index_segments([p1, p2], backend=FakeBackend([0, 20]))
    assert [Path(i.trajectory_path).name for i in out] == ["part1.xtc", "part2.xtc"]


def test_same_timeline():
    assert _idx("a", [0, 20]).same_timeline(_idx("b", [0, 20]))
    assert not _idx("a", [0, 20]).same_timeline(_idx("b", [0, 40]))


def test_index_view_without_file_returns_none():
    view = TrajectoryView(kind=ViewKind.RAW, requirements=TrajectoryRequirements())
    assert ti.index_view(view) is None


# ═══════════════════════════════════════════════════════════════════════════════
# GRO probe parser
# ═══════════════════════════════════════════════════════════════════════════════

def test_parse_probe_gro(tmp_path):
    g = tmp_path / "f.gro"
    g.write_text(
        "probe t= 100000016.00000 step= 70\n1\n    1PRB     PX    1   0.1   0.2   0.3\n"
        "   3.00000   3.00000   3.00000\n"
        "probe t=   0.00100 step= 80\n1\n    1PRB     PX    1   0.1   0.2   0.3\n"
        "   3.0 3.0 3.0 0.0 0.0 1.0 0.0 1.0 1.0\n")
    times, steps, boxes = ti.parse_probe_gro(g)
    assert times == [100000016.0, 0.001] and steps == [70, 80]
    assert boxes[0] == [3.0, 3.0, 3.0, 0, 0, 0, 0, 0, 0] and boxes[1][5] == 1.0


def test_parse_probe_gro_rejects_malformed(tmp_path):
    g = tmp_path / "f.gro"
    g.write_text("no time here\n1\n    1PRB     PX    1   0.1   0.2   0.3\n   3 3 3\n")
    with pytest.raises(ValueError):
        ti.parse_probe_gro(g)
    g.write_text("p t= 1.0 step= 1\n2\n    1PRB     PX    1   0.1   0.2   0.3\n")
    with pytest.raises(ValueError):
        ti.parse_probe_gro(g)


# ═══════════════════════════════════════════════════════════════════════════════
# Real GROMACS backend on synthetic XTC / TRR
# ═══════════════════════════════════════════════════════════════════════════════

def _write_traj(tmp_path: Path, times, suffix: str, name: str = "traj", natoms: int = 3) -> Path:
    gro = tmp_path / f"{name}_in.gro"
    lines = []
    for i, t in enumerate(times):
        lines.append(f"synthetic t= {t:.5f} step= {i * 10}")
        lines.append(str(natoms))
        for a in range(natoms):
            lines.append(f"{1:>5}{'SOL':<5}{'OW':>5}{a + 1:>5}"
                         f"{0.1 * (a % 20):8.3f}{0.2 + 0.001 * i:8.3f}{0.3:8.3f}")
        lines.append(f"   3.00000   3.00000   {3.0 + 0.01 * i:.5f}")
    gro.write_text("\n".join(lines) + "\n")
    out = tmp_path / f"{name}{suffix}"
    res = run_gmx(["trjconv", "-f", str(gro), "-s", str(gro), "-o", str(out)], stdin="0\n")
    assert res.ok and out.is_file(), res.stderr[-500:]
    return out


@requires_gmx
@pytest.mark.parametrize("suffix", [".xtc", ".trr"])
@pytest.mark.parametrize("times,state", [
    ([0, 20, 40, 60], ValidationState.VALID),
    ([50000, 50020, 50040, 50060], ValidationState.VALID),
    ([0, 20, 40, 100, 120], ValidationState.VALID),
    ([0, 20, 40, 40, 60], ValidationState.WARNING),
    ([0, 20, 40, 10, 30], ValidationState.WARNING),
])
def test_gromacs_backend_times(tmp_path, suffix, times, state):
    traj = _write_traj(tmp_path, times, suffix)
    idx, _ = ti.index_trajectory(traj)
    assert idx.times_ps == [float(t) for t in times]      # file order, nothing removed
    assert idx.steps == [i * 10 for i in range(len(times))]
    assert idx.state == state
    assert idx.backend == ti.GromacsTimeIndexBackend.id and idx.backend_version
    assert idx.backend_reported_frames == len(times)
    assert idx.boxes[1][2] == pytest.approx(3.01)
    assert idx.command[1] == "trjconv" and idx.stdin == "0\n"


@requires_gmx
@pytest.mark.parametrize("suffix", [".xtc", ".trr"])
def test_gromacs_backend_large_time_precision(tmp_path, suffix):
    # values gmx traj -ob (%g) and gmx dump (%.7e) cannot resolve
    times = [1234567.5, 1234568.0, 20000020.0, 20000040.0, 100000000.0, 100000016.0]
    idx, _ = ti.index_trajectory(_write_traj(tmp_path, times, suffix))
    assert idx.times_ps == times
    assert idx.strictly_increasing


@requires_gmx
def test_gromacs_backend_float32_precision_warning(tmp_path):
    # XTC time is float32: 1e8 + 20 is stored as 100000016 and 1e8 + 4 collapses to 1e8
    idx, _ = ti.index_trajectory(_write_traj(tmp_path, [100000000.0, 100000004.0], ".xtc"))
    assert idx.duplicate_groups == [[0, 1]]
    assert {"duplicate_timestamps", "duplicate_time_distinct_step"} <= _codes(idx)


@requires_gmx
def test_gromacs_backend_strided_derived_trajectory(tmp_path):
    src = _write_traj(tmp_path, [50000 + 10 * i for i in range(11)], ".xtc", "src")
    strided = tmp_path / "strided.xtc"
    res = run_gmx(["trjconv", "-f", str(src), "-o", str(strided), "-dt", "40"], stdin="0\n")
    assert res.ok, res.stderr[-400:]
    a, _ = ti.index_trajectory(src)
    b, _ = ti.index_trajectory(strided)
    assert b.times_ps == [50000.0, 50040.0, 50080.0]
    assert not a.same_timeline(b)                    # derived view scanned, not inherited
    view = TrajectoryView(kind=ViewKind.WHOLE, requirements=TrajectoryRequirements(),
                          path=str(strided))
    assert ti.index_view(view).times_ps == b.times_ps


@requires_gmx
def test_gromacs_backend_truncated_xtc(tmp_path):
    traj = _write_traj(tmp_path, [0, 20, 40, 60, 80], ".xtc", natoms=50)
    data = traj.read_bytes()
    cut = tmp_path / "cut.xtc"
    cut.write_bytes(data[: len(data) - 40])
    idx, _ = ti.index_trajectory(cut)
    assert idx.truncated and "Incomplete" in idx.truncation_evidence
    assert idx.times_ps == [0.0, 20.0, 40.0, 60.0]
    assert idx.state == ValidationState.WARNING
    assert "truncated_trajectory" in _codes(idx)


@requires_gmx
def test_gromacs_backend_segments(tmp_path):
    a = _write_traj(tmp_path, [0, 20, 40], ".xtc", "part1")
    b = _write_traj(tmp_path, [40, 60, 80], ".xtc", "part2")
    ia, ib = ti.index_segments([a, b])
    (rel,) = ti.relate_segments([ia, ib])
    assert rel.relation == SegmentRelationKind.BOUNDARY_SHARED and rel.shared_timestamps == 1


@requires_gmx
def test_gromacs_backend_cache_roundtrip(tmp_path):
    traj = _write_traj(tmp_path, [0, 20, 40], ".trr")
    cache = tmp_path / "cache"
    first, c1 = ti.index_trajectory(traj, cache_dir=cache)
    second, c2 = ti.index_trajectory(traj, cache_dir=cache)
    assert (c1, c2) == (False, True)
    assert second.times_ps == first.times_ps and second.backend_version == first.backend_version
    assert FrameTimeMap.from_index(second).nearest_frame(30.0).frame == 1


# ═══════════════════════════════════════════════════════════════════════════════
# Real trajectory — compared with independent `gmx check` output
# ═══════════════════════════════════════════════════════════════════════════════

@requires_gmx
@requires_real_traj
def test_real_trajectory_matches_gmx_check():
    idx, _ = ti.index_trajectory(REAL_XTC)            # read-only; no cache written
    res = run_gmx(["check", "-f", str(REAL_XTC)])
    chk = _parse_gmx_check(res.stderr + "\n" + res.stdout)
    assert idx.state == ValidationState.VALID
    assert idx.n_frames == chk["n_frames"]
    assert idx.start_time_ps == chk["start_time_ps"]
    assert idx.end_time_ps == chk["end_time_ps"]
    assert idx.strictly_increasing
