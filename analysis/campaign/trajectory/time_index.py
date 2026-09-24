"""Per-frame trajectory time index — the actual timestamp of every frame.

Nothing here assumes ``time = start + frame * dt`` or a zero start.  Times are
*read* from the trajectory by a backend, validated, and cached by the
trajectory's content fingerprint.

Backend abstraction
-------------------
A :class:`TimeIndexBackend` extracts raw per-frame data; the backend-neutral
:func:`validate_timeline` then evaluates it.  Only the authoritative GROMACS
backend exists today; MDAnalysis / native-XDR backends can be added later by
implementing the same protocol (a native parser must first be validated
frame-by-frame against this backend).

GROMACS strategy (verified on GROMACS 2025.2)
---------------------------------------------
``gmx trjconv -f TRAJ -s probe.gro -n probe.ndx -o frames.gro`` with a
self-generated one-atom reference and index writes a one-atom GRO frame per
trajectory frame whose title ends in ``t= %.5f step= N``.  This route was
chosen over the alternatives because:

* ``gmx traj -ob`` prints times with ``%g`` (6 significant digits): at
  >= 1e6 ps (1 µs) distinct frames become indistinguishable
  (1234567.5 and 1234568.0 both print as ``1.23457e+06``);
* ``gmx dump -f`` prints ``%.7e`` (8 significant digits, not enough for every
  float32 value, e.g. 100000016 -> ``1.0000002e+08``) and dumps every
  coordinate;
* ``gmx check`` only reports a subset of frames.

The GRO title keeps 1e-5 ps resolution at any magnitude, carries the MD step,
preserves file order (duplicates / non-monotonic times are kept), needs no
user topology, and truncation is reported on stderr
(``WARNING: Incomplete frame`` / ``Incomplete header``).

Nothing here transforms or writes next to the source trajectory; the probe
files live in a temporary directory.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import re
import tempfile
from pathlib import Path
from typing import Optional, Protocol, Sequence, runtime_checkable

from analysis.campaign.fingerprint import fingerprint_file
from analysis.campaign.gmx import gmx_available, gmx_version, run_gmx
from analysis.campaign.models import (
    TIME_INDEX_SCHEMA_VERSION, CampaignWarning, FrameTimeIndex, SegmentRelation,
    SegmentRelationKind, Severity, TrajectoryView, ValidationState,
)

SUPPORTED_SUFFIXES = (".xtc", ".trr")

_GRO_TIME_RE = re.compile(r"\bt=\s*(\S+)\s+step=\s*(-?\d+)\s*$")
_LAST_WRITTEN_RE = re.compile(r"Last written:\s*frame\s+(\d+)\s+time\s+(\S+)")
_INCOMPLETE_RE = re.compile(r"WARNING:\s*Incomplete (frame|header)[^\n]*")

_INDEX_JSON = "index.json"
_FRAMES_NPZ = "frames.npz"


# ═══════════════════════════════════════════════════════════════════════════════
# Backend protocol
# ═══════════════════════════════════════════════════════════════════════════════

@runtime_checkable
class TimeIndexBackend(Protocol):
    id: str                               # stable id; part of the cache key

    def available(self) -> bool: ...

    def version(self) -> Optional[str]: ...

    def scan(self, path: Path) -> tuple[FrameTimeIndex, bool]:
        """Extract per-frame data (unvalidated).

        Returns ``(index, ran_ok)``; ``ran_ok`` is False when the tool itself
        failed, in which case the result must not be cached.
        """
        ...


class GromacsTimeIndexBackend:
    """Authoritative backend: one-atom GRO dump via ``gmx trjconv``."""

    id = "gromacs-trjconv-gro/v1"
    time_resolution_ps = 1e-5            # GRO title prints t with %.5f (ps)

    def __init__(self, gmx: str = "gmx", timeout: float = 4 * 3600.0):
        self.gmx = gmx
        self.timeout = timeout

    def available(self) -> bool:
        return gmx_available(self.gmx)

    def version(self) -> Optional[str]:
        return gmx_version(self.gmx)

    def scan(self, path: Path) -> tuple[FrameTimeIndex, bool]:
        path = Path(path)
        idx = FrameTimeIndex(
            trajectory_path=str(path.resolve()), fingerprint=None,
            backend=self.id, backend_version=self.version(),
            time_unit_reported="ps", time_resolution_ps=self.time_resolution_ps,
        )
        with tempfile.TemporaryDirectory(prefix="simforge_tindex_") as tmp:
            tmp_p = Path(tmp)
            ref = tmp_p / "probe.gro"
            ndx = tmp_p / "probe.ndx"
            out = tmp_p / "frames.gro"
            ref.write_text("simforge time-index probe\n1\n"
                           "    1PRB     PX    1   0.000   0.000   0.000\n"
                           "   1.00000   1.00000   1.00000\n")
            ndx.write_text("[ probe ]\n1\n")
            args = ["trjconv", "-f", str(path.resolve()), "-s", str(ref), "-n", str(ndx),
                    "-o", str(out)]
            stdin = "0\n"
            res = run_gmx(args, stdin=stdin, gmx=self.gmx, timeout=self.timeout)
            idx.command = list(res.argv)
            idx.stdin = stdin
            stderr = res.stderr or ""

            m = _INCOMPLETE_RE.search(stderr)
            if m:
                idx.truncated = True
                idx.truncation_evidence = m.group(0).strip()
            m = _LAST_WRITTEN_RE.search(stderr)
            if m:
                idx.backend_reported_frames = int(m.group(1)) + 1

            if not res.ok:
                idx.state = ValidationState.INVALID
                idx.warnings.append(CampaignWarning(
                    "time_index_backend_failed",
                    f"`gmx trjconv` failed (rc={res.returncode}): "
                    f"{stderr.strip()[-400:]}",
                    Severity.ERROR,
                ))
                return idx, False
            if not out.is_file():
                # rc 0 but nothing written: an empty / unreadable trajectory
                idx.backend_reported_frames = idx.backend_reported_frames or 0
                return idx, True

            try:
                times, steps, boxes = parse_probe_gro(out)
            except ValueError as exc:
                idx.state = ValidationState.INVALID
                idx.warnings.append(CampaignWarning(
                    "time_index_parse_failed", str(exc), Severity.ERROR))
                return idx, False
        idx.times_ps, idx.steps, idx.boxes = times, steps, boxes
        return idx, True


def parse_probe_gro(path: Path) -> tuple[list[float], list[int], list[list[float]]]:
    """Parse a multi-frame GRO written by trjconv: title (with ``t= step=``),
    atom count, atom lines, box line.  Raises ``ValueError`` on any
    structural inconsistency — a partially understood file is never used."""
    times: list[float] = []
    steps: list[int] = []
    boxes: list[list[float]] = []
    with Path(path).open(errors="replace") as fh:
        lineno = 0
        while True:
            title = fh.readline()
            if not title:
                break
            lineno += 1
            if not title.strip():
                continue
            m = _GRO_TIME_RE.search(title.rstrip("\n"))
            if not m:
                raise ValueError(f"frame title without 't= ... step= ...' at line {lineno}: {title!r}")
            times.append(float(m.group(1)))       # float('nan') / 'inf' preserved for validation
            steps.append(int(m.group(2)))
            count_line = fh.readline()
            lineno += 1
            try:
                natoms = int(count_line.strip())
            except ValueError:
                raise ValueError(f"bad atom-count line at line {lineno}: {count_line!r}") from None
            for _ in range(natoms):
                if not fh.readline():
                    raise ValueError(f"GRO frame {len(times) - 1} ends inside its atom block")
                lineno += 1
            box_line = fh.readline()
            lineno += 1
            try:
                vals = [float(v) for v in box_line.split()]
            except ValueError:
                raise ValueError(f"bad box line at line {lineno}: {box_line!r}") from None
            if len(vals) == 3:
                vals = vals + [0.0] * 6
            if len(vals) != 9:
                raise ValueError(f"box line with {len(vals)} values at line {lineno}")
            boxes.append(vals)
    return times, steps, boxes


# ═══════════════════════════════════════════════════════════════════════════════
# Timeline validation (backend-neutral)
# ═══════════════════════════════════════════════════════════════════════════════

def validate_timeline(idx: FrameTimeIndex, *, expected_n_frames: Optional[int] = None) -> FrameTimeIndex:
    """Evaluate the timeline in place and return it.

    Reports — never repairs — empty trajectories, non-finite times, duplicate
    and non-monotonic timestamps, truncation, and frame-count mismatches.
    """
    if idx.state == ValidationState.INVALID and not idx.times_ps:
        return idx                       # backend already failed; nothing to evaluate
    t = idx.times_ps
    n = len(t)
    invalid = False
    warn = False
    idx.expected_n_frames = expected_n_frames

    if n == 0:
        idx.warnings.append(CampaignWarning(
            "empty_trajectory", "no frames could be read from the trajectory", Severity.ERROR))
        invalid = True

    idx.non_finite_frames = [i for i, v in enumerate(t) if not math.isfinite(v)]
    if idx.non_finite_frames:
        idx.warnings.append(CampaignWarning(
            "non_finite_timestamps",
            f"{len(idx.non_finite_frames)} frame(s) have non-finite timestamps "
            f"(first: frame {idx.non_finite_frames[0]})", Severity.ERROR))
        invalid = True

    finite = not idx.non_finite_frames
    if n and finite:
        groups: dict[float, list[int]] = {}
        for i, v in enumerate(t):
            groups.setdefault(v, []).append(i)
        idx.duplicate_groups = sorted((g for g in groups.values() if len(g) > 1), key=lambda g: g[0])
        idx.non_monotonic_steps = [[i, i + 1] for i in range(n - 1) if t[i + 1] < t[i]]
        idx.non_decreasing = not idx.non_monotonic_steps
        idx.strictly_increasing = idx.non_decreasing and not idx.duplicate_groups

        if idx.duplicate_groups:
            n_dup = sum(len(g) for g in idx.duplicate_groups)
            first = idx.duplicate_groups[0]
            msg = (f"{len(idx.duplicate_groups)} timestamp(s) shared by {n_dup} frames "
                   f"(first: {t[first[0]]} ps at frames {first}); frames kept as-is")
            idx.warnings.append(CampaignWarning("duplicate_timestamps", msg, Severity.REVIEW))
            warn = True
            if idx.steps and any(len({idx.steps[i] for i in g}) > 1 for g in idx.duplicate_groups):
                idx.warnings.append(CampaignWarning(
                    "duplicate_time_distinct_step",
                    "frames with identical reported times carry different MD steps — "
                    "either the timeline is genuinely duplicated (e.g. re-timed "
                    "concatenation) or the stored time precision cannot resolve them; "
                    "not guessed", Severity.REVIEW))
        if idx.non_monotonic_steps:
            a, b = idx.non_monotonic_steps[0]
            idx.warnings.append(CampaignWarning(
                "non_monotonic_timestamps",
                f"time decreases at {len(idx.non_monotonic_steps)} position(s) "
                f"(first: frame {a} {t[a]} ps -> frame {b} {t[b]} ps); frames not reordered",
                Severity.REVIEW))
            warn = True
        _check_float32_precision(idx)

    if idx.truncated:
        idx.warnings.append(CampaignWarning(
            "truncated_trajectory",
            f"trajectory ends with an incomplete frame ({idx.truncation_evidence}); "
            f"only the {n} complete frame(s) are indexed", Severity.WARN))
        warn = True

    if idx.backend_reported_frames is not None and idx.backend_reported_frames != n:
        idx.warnings.append(CampaignWarning(
            "frame_count_mismatch",
            f"backend reported {idx.backend_reported_frames} frame(s) but "
            f"{n} were parsed", Severity.ERROR))
        invalid = True

    if expected_n_frames is not None and expected_n_frames != n:
        idx.warnings.append(CampaignWarning(
            "unexpected_frame_count",
            f"expected {expected_n_frames} frame(s), indexed {n}", Severity.REVIEW))
        warn = True

    if invalid:
        idx.state = ValidationState.INVALID
    elif warn:
        idx.state = ValidationState.WARNING
    else:
        idx.state = ValidationState.VALID
    return idx


def _check_float32_precision(idx: FrameTimeIndex) -> None:
    """XTC stores time as float32.  Warn when its spacing at the largest time
    approaches the frame spacing — timestamps may then no longer resolve frames."""
    if not idx.trajectory_path.lower().endswith(".xtc") or len(idx.times_ps) < 2:
        return
    import numpy as np
    t = sorted(set(idx.times_ps))
    if len(t) < 2:
        return
    min_dt = min(b - a for a, b in zip(t, t[1:]))
    ulp = float(np.spacing(np.float32(max(abs(t[0]), abs(t[-1])))))
    if ulp * 2 >= min_dt:
        idx.warnings.append(CampaignWarning(
            "time_precision_limit",
            f"XTC float32 time spacing at {t[-1]} ps is {ulp} ps, comparable to the "
            f"smallest frame interval {min_dt} ps — timestamps may not resolve frames",
            Severity.REVIEW))


# ═══════════════════════════════════════════════════════════════════════════════
# Cache
# ═══════════════════════════════════════════════════════════════════════════════

def default_cache_dir() -> Path:
    """``$SIMFORGE_CACHE_DIR/time_index`` or ``$XDG_CACHE_HOME/simforge/time_index``."""
    base = os.environ.get("SIMFORGE_CACHE_DIR")
    if base:
        return Path(base) / "time_index"
    xdg = os.environ.get("XDG_CACHE_HOME") or str(Path.home() / ".cache")
    return Path(xdg) / "simforge" / "time_index"


def cache_key(fingerprint_digest: str, backend_id: str, backend_version: Optional[str]) -> str:
    parts = [TIME_INDEX_SCHEMA_VERSION, backend_id, backend_version or "unknown",
             fingerprint_digest]
    return hashlib.sha256("\x00".join(parts).encode()).hexdigest()[:24]


def _save(idx: FrameTimeIndex, entry: Path) -> None:
    import numpy as np
    entry.mkdir(parents=True, exist_ok=True)
    npz_tmp = entry / (_FRAMES_NPZ + ".tmp.npz")
    arrays = {
        "times_ps": np.asarray(idx.times_ps, dtype=np.float64),
        "steps": np.asarray(idx.steps, dtype=np.int64),
    }
    if idx.boxes is not None:
        arrays["boxes"] = np.asarray(idx.boxes, dtype=np.float64).reshape(-1, 9)
    np.savez(npz_tmp, **arrays)
    os.replace(npz_tmp, entry / _FRAMES_NPZ)
    # index.json is written last: its presence marks a complete entry
    json_tmp = entry / (_INDEX_JSON + ".tmp")
    json_tmp.write_text(json.dumps(idx.to_dict(), indent=2) + "\n")
    os.replace(json_tmp, entry / _INDEX_JSON)


def _load(entry: Path, *, digest: str, backend_id: str,
          backend_version: Optional[str]) -> Optional[FrameTimeIndex]:
    jf, nf = entry / _INDEX_JSON, entry / _FRAMES_NPZ
    if not (jf.is_file() and nf.is_file()):
        return None
    try:
        import numpy as np
        meta = json.loads(jf.read_text())
        if (meta.get("schema_version") != TIME_INDEX_SCHEMA_VERSION
                or meta.get("backend") != backend_id
                or meta.get("backend_version") != backend_version
                or (meta.get("fingerprint") or {}).get("digest") != digest):
            return None
        with np.load(nf) as z:
            times = [float(v) for v in z["times_ps"]]
            steps = [int(v) for v in z["steps"]]
            boxes = z["boxes"].tolist() if "boxes" in z.files else None
        if len(times) != meta.get("n_frames"):
            return None
        return FrameTimeIndex.from_dict(meta, times_ps=times, steps=steps, boxes=boxes)
    except (OSError, ValueError, KeyError, json.JSONDecodeError):
        return None


# ═══════════════════════════════════════════════════════════════════════════════
# Public API
# ═══════════════════════════════════════════════════════════════════════════════

def index_trajectory(
    path: str | Path,
    *,
    backend: Optional[TimeIndexBackend] = None,
    cache_dir: Optional[str | Path] = None,
    force: bool = False,
    strong_fingerprint: bool = False,
    expected_n_frames: Optional[int] = None,
) -> tuple[FrameTimeIndex, bool]:
    """Index one trajectory file.  Returns ``(index, from_cache)``.

    ``cache_dir=None`` disables persistent caching (nothing is written).
    The cache entry is keyed by the trajectory content fingerprint (existing
    ``fingerprint_file`` policy), backend id, backend version and index schema.
    ``expected_n_frames`` is a validation hint, not part of the cache identity.
    """
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"trajectory not found: {p}")
    if p.suffix.lower() not in SUPPORTED_SUFFIXES:
        raise ValueError(f"unsupported trajectory format {p.suffix!r}; "
                         f"supported: {', '.join(SUPPORTED_SUFFIXES)}")
    backend = backend or GromacsTimeIndexBackend()
    fp = fingerprint_file(p, strong=strong_fingerprint)

    if not backend.available():
        idx = FrameTimeIndex(trajectory_path=str(p.resolve()), fingerprint=fp,
                             backend=backend.id, backend_version=None,
                             state=ValidationState.INVALID)
        idx.warnings.append(CampaignWarning(
            "time_index_backend_unavailable",
            f"time-index backend '{backend.id}' is not available", Severity.ERROR))
        return idx, False

    version = backend.version()
    entry = None
    if cache_dir is not None:
        entry = Path(cache_dir) / cache_key(fp.digest, backend.id, version)
        if not force:
            cached = _load(entry, digest=fp.digest, backend_id=backend.id,
                           backend_version=version)
            if cached is not None:
                cached.trajectory_path = str(p.resolve())
                if expected_n_frames is not None:
                    # re-evaluate with the caller's hint; the cached entry is untouched
                    validate_timeline(_reset_validation(cached), expected_n_frames=expected_n_frames)
                return cached, True

    idx, ran_ok = backend.scan(p)
    idx.fingerprint = fp
    idx.backend_version = version
    if ran_ok:
        validate_timeline(idx)             # cache the file's intrinsic validation
        if entry is not None:
            _save(idx, entry)
        if expected_n_frames is not None:
            validate_timeline(_reset_validation(idx), expected_n_frames=expected_n_frames)
    return idx, False


def _reset_validation(idx: FrameTimeIndex) -> FrameTimeIndex:
    """Drop derived validation so :func:`validate_timeline` can be rerun."""
    keep = {"time_index_backend_failed", "time_index_parse_failed"}
    idx.warnings = [w for w in idx.warnings if w.code in keep]
    idx.state = ValidationState.VALID
    idx.duplicate_groups, idx.non_monotonic_steps, idx.non_finite_frames = [], [], []
    return idx


def load_cached_index(
    path: str | Path,
    *,
    cache_dir: str | Path,
    backend_id: str,
    backend_version: Optional[str],
    strong_fingerprint: bool = False,
) -> Optional[FrameTimeIndex]:
    """Load a previously computed index by its recorded identity (current file
    content digest + backend id/version) without running any backend.
    Returns None when no matching, valid cache entry exists."""
    p = Path(path)
    if not p.is_file():
        return None
    fp = fingerprint_file(p, strong=strong_fingerprint)
    entry = Path(cache_dir) / cache_key(fp.digest, backend_id, backend_version)
    idx = _load(entry, digest=fp.digest, backend_id=backend_id, backend_version=backend_version)
    if idx is not None:
        idx.trajectory_path = str(p.resolve())
    return idx


def index_segments(paths: Sequence[str | Path], **kwargs) -> list[FrameTimeIndex]:
    """Index each segment independently, in the given order.  Never concatenates."""
    return [index_trajectory(p, **kwargs)[0] for p in paths]


def index_view(view: TrajectoryView, **kwargs) -> Optional[FrameTimeIndex]:
    """Index the actual file of a :class:`TrajectoryView` (None if it has no file).

    Derived views are scanned, never assumed to keep their source timeline;
    use :meth:`FrameTimeIndex.same_timeline` to demonstrate equality.
    """
    if not view.path or not Path(view.path).is_file():
        return None
    return index_trajectory(view.path, **kwargs)[0]


def relate_segments(indices: Sequence[FrameTimeIndex]) -> list[SegmentRelation]:
    """Describe every ordered pair (i < j, in the given order) in time.

    Purely descriptive: no precedence is chosen and nothing is merged.
    """
    rels: list[SegmentRelation] = []
    for i in range(len(indices)):
        for j in range(i + 1, len(indices)):
            rels.append(_relate(indices[i], indices[j]))
    return rels


def _relate(a: FrameTimeIndex, b: FrameTimeIndex) -> SegmentRelation:
    rel = SegmentRelation(first_path=a.trajectory_path, second_path=b.trajectory_path,
                          relation=SegmentRelationKind.UNKNOWN)
    fa = [v for v in a.times_ps if math.isfinite(v)]
    fb = [v for v in b.times_ps if math.isfinite(v)]
    if not fa or not fb:
        return rel
    a0, a1, b0, b1 = min(fa), max(fa), min(fb), max(fb)
    rel.first_start_ps, rel.first_end_ps = a0, a1
    rel.second_start_ps, rel.second_end_ps = b0, b1
    rel.gap_ps = b0 - a1
    rel.shared_timestamps = len(set(fa) & set(fb))
    overlap = min(a1, b1) - max(a0, b0)
    rel.overlap_ps = overlap if overlap > 0 else 0.0
    if b0 > a1:
        rel.relation = SegmentRelationKind.SEQUENTIAL
    elif b0 == a1 and b1 >= a1:
        rel.relation = SegmentRelationKind.BOUNDARY_SHARED
    elif b1 < a0:
        rel.relation = SegmentRelationKind.REVERSED
    else:
        rel.relation = SegmentRelationKind.OVERLAP
    return rel


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Developer entry point: ``python -m analysis.campaign.trajectory.time_index TRAJ``."""
    import argparse
    ap = argparse.ArgumentParser(description="Print the per-frame time index of a trajectory (JSON).")
    ap.add_argument("trajectory")
    ap.add_argument("--gmx", default="gmx")
    ap.add_argument("--cache-dir", default=None)
    ap.add_argument("--frames", action="store_true", help="include per-frame times/steps")
    ns = ap.parse_args(argv)
    idx, cached = index_trajectory(ns.trajectory, backend=GromacsTimeIndexBackend(ns.gmx),
                                   cache_dir=ns.cache_dir)
    out = idx.to_dict(include_frames=ns.frames)
    out["from_cache"] = cached
    print(json.dumps(out, indent=2))
    return 0 if idx.state != ValidationState.INVALID else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
