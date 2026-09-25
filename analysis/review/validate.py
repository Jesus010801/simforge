"""Reopen and validate a prepared ``ReviewDataset`` without rediscovering or
reanalysing the project.

Validation never trusts stored metadata: sources are re-fingerprinted (in the
mode they were recorded with), every store file is re-hashed, the timeline
is reloaded and checked, and every reference (views, annotations) must
resolve.  A result whose frozen file no longer matches is reported as
unavailable and makes the dataset invalid — stale science is never shown.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from analysis.review.dataset import (
    DATASET_FILENAME, SUPPORTED_SCHEMAS, CapabilityState, ReviewDataset,
)
from analysis.review.store import verify


class ReviewDatasetError(ValueError):
    """The file is not a readable review dataset."""


@dataclass
class ReviewValidation:
    valid: bool
    problems: list[str] = field(default_factory=list)
    unavailable: dict[str, str] = field(default_factory=dict)   # instance_id -> reason
    checked: dict[str, int] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {"valid": self.valid, "problems": self.problems,
                "unavailable": self.unavailable, "checked": self.checked}


def _dataset_file(path: str | Path) -> Path:
    p = Path(path)
    return p / DATASET_FILENAME if p.is_dir() else p


def load_review_dataset(path: str | Path) -> ReviewDataset:
    """Read ``review_dataset.json`` (or a session directory)."""
    f = _dataset_file(path)
    try:
        d = json.loads(f.read_text())
    except (OSError, ValueError) as exc:
        raise ReviewDatasetError(f"cannot read review dataset {f}: {exc}") from exc
    if not isinstance(d, dict) or "session_id" not in d:
        raise ReviewDatasetError(f"{f} is not a review dataset")
    return ReviewDataset.from_dict(d)


def _resolve(base: Path, path: str) -> Path:
    p = Path(path)
    return p if p.is_absolute() else base / p


def _store_ok(base: Path, ref: Optional[dict]) -> Optional[str]:
    """None if the frozen file is intact, else the reason."""
    if not ref:
        return None
    p = _resolve(base, ref["path"])
    if not p.is_file():
        return f"{ref['path']} is missing"
    if not verify(p, ref["sha256"]):
        return f"{ref['path']} changed after the session was prepared"
    return None


def _fingerprint_ok(entry: Optional[dict]) -> Optional[str]:
    from analysis.campaign.fingerprint import fingerprint_file
    if not entry or not entry.get("path"):
        return None
    p = Path(entry["path"])
    fp = entry.get("fingerprint") or {}
    if not p.is_file():
        return f"{p} is missing"
    if fp and fingerprint_file(p, strong=fp.get("mode") == "strong").digest != fp.get("digest"):
        return f"{p.name} changed since the session was prepared"
    return None


def validate_review_dataset(ds: ReviewDataset, session_dir: str | Path) -> ReviewValidation:
    base = Path(session_dir)
    problems: list[str] = []
    unavailable: dict[str, str] = {}
    checked = {"sources": 0, "store_files": 0, "results": 0, "views": 0, "annotations": 0}

    if ds.schema_version not in SUPPORTED_SCHEMAS:
        problems.append(f"unsupported schema {ds.schema_version!r}")
        return ReviewValidation(False, problems, unavailable, checked)
    if ds.status != "prepared":
        problems.append(f"dataset status is {ds.status!r}, not 'prepared'")

    # sources
    for role in ("trajectory", "topology", "structure", "index"):
        entry = ds.sources.get(role)
        if entry:
            checked["sources"] += 1
            why = _fingerprint_ok(entry)
            if why:
                problems.append(f"source {role}: {why}")

    # display view
    dv = ds.display_view
    if dv.get("status") == CapabilityState.AVAILABLE:
        checked["views"] += 1
        path = dv.get("path")
        if not path or not _resolve(base, path).is_file():
            problems.append("display view: trajectory file missing")
        elif dv.get("kind") != "raw":
            why = _fingerprint_ok({"path": str(_resolve(base, path)),
                                   "fingerprint": dv.get("output_fingerprint")})
            if why:
                problems.append(f"display view: {why}")
    elif not dv.get("status"):
        problems.append("display view: not described")

    # timeline
    tl = ds.timeline
    times = None
    for key in ("times_ref", "steps_ref", "index_ref"):
        if tl.get(key):
            checked["store_files"] += 1
            why = _store_ok(base, tl[key])
            if why:
                problems.append(f"timeline: {why}")
    if tl.get("times_ref") and not any(p.startswith("timeline") for p in problems):
        import numpy as np
        times = np.load(_resolve(base, tl["times_ref"]["path"])).tolist()
        if len(times) != tl.get("n_frames"):
            problems.append(f"timeline: {len(times)} stored frame times vs {tl.get('n_frames')} recorded")
        try:
            from analysis.review.sync import ReviewTimeline
            ReviewTimeline.from_dataset(ds, base)
        except ValueError as exc:
            if tl.get("state") != "invalid":
                problems.append(f"timeline: frame/time mapping unavailable ({exc})")
    elif ds.available():
        problems.append("timeline: no review timeline for a session with results")

    # diagnostics report snapshot
    if ds.diagnostics.get("report_ref"):
        checked["store_files"] += 1
        why = _store_ok(base, ds.diagnostics["report_ref"])
        if why:
            problems.append(f"diagnostics: {why}")

    # observables
    view_refs = {v.get("view_ref") for v in ds.analysis_views}
    ann_ids = {a["annotation_id"] for a in ds.annotations}
    for o in ds.observables:
        for a in o.annotations_used + o.blocking_annotations:
            checked["annotations"] += 1
            if a not in ann_ids and o.available:
                problems.append(f"{o.instance_id}: annotation '{a}' not described in the session")
        if not o.available:
            continue
        for ref in (o.provenance_ref, (o.reuse or {}).get("report_ref")):
            if ref:
                checked["store_files"] += 1
                why = _store_ok(base, ref)
                if why:
                    unavailable[o.instance_id] = why
        for r in o.results:
            checked["results"] += 1
            st = r.array.storage
            ref = {"path": st.path, "sha256": (st.attrs or {}).get("sha256")} if st else None
            why = ("no storage" if not ref or not ref["sha256"] else _store_ok(base, ref))
            if why:
                unavailable[o.instance_id] = f"result {r.array.name}: {why}"
            if r.array.view_ref not in view_refs:
                problems.append(f"{o.instance_id}: view_ref {r.array.view_ref} not described")
            if (r.sync.get("syncable") and r.sync.get("mapping") == "row_i_is_frame_i"
                    and times is not None and r.sync.get("n_samples") not in (None, len(times))):
                problems.append(f"{o.instance_id}: {r.sync.get('n_samples')} samples vs "
                                f"{len(times)} frames for a row-per-frame mapping")
    for iid, why in unavailable.items():
        problems.append(f"{iid}: {why}")
    return ReviewValidation(not problems, problems, unavailable, checked)


def open_review_session(path: str | Path) -> tuple[ReviewDataset, ReviewValidation]:
    """Load + validate a prepared session (what a viewer layer starts from)."""
    f = _dataset_file(path)
    ds = load_review_dataset(f)
    return ds, validate_review_dataset(ds, f.parent)
