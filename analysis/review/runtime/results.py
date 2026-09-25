"""Read-only presentation data for a validated ReviewDataset.

* Results are addressed by catalog ids (``r<observable>-<array>``), never by
  path.  Only files the dataset references inside its review store are read,
  and each is re-hashed against the sha256 recorded at preparation.
* Values are served exactly as stored (non-finite values become ``null`` and
  are counted) — no interpolation, smoothing, averaging or downsampling.
* Renderer choice is by axis signature only; an unknown signature stays a
  valid scientific result, described but not drawn.
"""
from __future__ import annotations

import math
import threading
from dataclasses import replace
from pathlib import Path
from typing import Optional

from analysis.review.store import sha256_file

DIFFERENT_COORDINATES_NOTE = ("Synchronized by simulation time; coordinates used for analysis "
                              "differ from displayed coordinates.")


class Renderer:
    TIMESERIES = "timeseries"        # [time]
    MULTILINE = "multiline"          # [time, category]
    PROFILE = "profile"              # [index] | [category]
    UNSUPPORTED = "unsupported"


def renderer_for(signature) -> str:
    sig = tuple(signature)
    if sig == ("time",):
        return Renderer.TIMESERIES
    if sig == ("time", "category"):
        return Renderer.MULTILINE
    if sig in (("index",), ("category",)):
        return Renderer.PROFILE
    return Renderer.UNSUPPORTED


class ResultUnavailable(LookupError):
    """Unknown id, or the referenced frozen file is missing / changed / not authorised."""


def _clean(values) -> tuple[list, int]:
    out, bad = [], 0
    for v in values:
        if isinstance(v, float) and not math.isfinite(v):
            out.append(None)
            bad += 1
        else:
            out.append(v)
    return out, bad


class ResultCatalog:
    def __init__(self, dataset, session_dir: str | Path, timeline):
        self.ds = dataset
        self.base = Path(session_dir).resolve()
        self.store = (self.base.parent.parent / "store").resolve()
        self.timeline = timeline
        self._entries: dict[str, tuple] = {}
        self._cache: dict[str, dict] = {}
        self._lock = threading.Lock()
        for i, o in enumerate(dataset.observables):
            for j, r in enumerate(o.results):
                self._entries[f"r{i}-{j}"] = (o, r)
        self._authorized = {self._resolve(p) for p in self._referenced_paths()}

    # ── authorisation ──────────────────────────────────────────────────────
    def _resolve(self, rel: str) -> Path:
        return (self.base / rel).resolve()

    def _referenced_paths(self):
        for o, r in self._entries.values():
            if r.array.storage:
                yield r.array.storage.path
            for ax in r.array.axes:
                if ax.values_ref:
                    yield ax.values_ref.path

    def _authorized_file(self, ref) -> Path:
        if Path(ref.path).is_absolute():
            raise ResultUnavailable("result storage is not inside the review store")
        p = self._resolve(ref.path)
        if p not in self._authorized or self.store not in p.parents:
            raise ResultUnavailable("result storage is not an authorised review-store file")
        sha = (ref.attrs or {}).get("sha256")
        if not p.is_file() or not sha or sha256_file(p) != sha:
            raise ResultUnavailable(f"frozen result file {p.name} is missing or changed")
        return p

    def _read(self, ref, column: Optional[int] = None) -> list:
        from analysis.campaign.results import read_column
        p = self._authorized_file(ref)
        if ref.format in ("xvg", "csv"):
            return read_column(replace(ref, path=str(p),
                                       column=ref.column if column is None else column))
        if ref.format == "npy":
            import numpy as np
            arr = np.load(p, allow_pickle=False)
            return arr.tolist() if column is None else arr[:, column].tolist()
        raise ResultUnavailable(f"storage format {ref.format!r} is not rendered in this version")

    # ── description ────────────────────────────────────────────────────────
    def ids(self) -> list[str]:
        return list(self._entries)

    def view_relation(self, view_ref: Optional[str]) -> dict:
        dv = self.ds.display_view
        av = self.ds.view(view_ref) or {}
        same = bool(view_ref) and view_ref == dv.get("view_ref")
        return {"display_view_ref": dv.get("view_ref"), "display_kind": dv.get("kind"),
                "analysis_view_ref": view_ref, "analysis_kind": av.get("kind"),
                "analysis_purpose": av.get("purpose"),
                "analysis_operations": [d.get("operation") for d in av.get("decisions", [])
                                        if d.get("applied")],
                "same_coordinates": same,
                "synchronization": "simulation time",
                "note": "" if same else DIFFERENT_COORDINATES_NOTE}

    def describe(self, rid: str) -> dict:
        if rid not in self._entries:
            raise ResultUnavailable(f"no result {rid!r} in this session")
        o, r = self._entries[rid]
        a = r.array
        return {"id": rid, "instance_id": o.instance_id, "observable": o.observable,
                "display_name": o.display_name, "availability": o.availability,
                "externally_supplied": o.externally_supplied,
                "independently_verified": o.independently_verified,
                "name": a.name, "quantity": a.quantity, "unit": a.unit,
                "axis_signature": list(a.axis_signature()),
                "axes": [{"name": x.name, "kind": x.kind, "unit": x.unit,
                          "n_values": len(x.values) if x.values is not None else None}
                         for x in a.axes],
                "storage_format": a.storage.format if a.storage else None,
                "renderer": renderer_for(a.axis_signature()),
                "syncable": bool(r.sync.get("syncable")), "sync": r.sync,
                "view": self.view_relation(a.view_ref)}

    # ── data ───────────────────────────────────────────────────────────────
    def data(self, rid: str) -> dict:
        with self._lock:
            if rid in self._cache:
                return self._cache[rid]
        out = self.describe(rid)
        _, r = self._entries[rid]
        a = r.array
        kind = out["renderer"]
        if kind == Renderer.UNSUPPORTED:
            out["data"] = None
            out["message"] = ("Scientific result available. Renderer not implemented for axis "
                              f"signature: ({', '.join(a.axis_signature())})")
        elif a.storage is None:
            raise ResultUnavailable("result has no storage")
        elif kind == Renderer.TIMESERIES:
            t, y = self._axis_values(a.axes[0]), self._read(a.storage)
            if len(t) != len(y):
                raise ResultUnavailable(f"{len(t)} time values vs {len(y)} samples")
            y, bad = _clean(y)
            out["data"] = {"time_ps": t, "series": [{"label": a.name, "values": y}],
                           "n_samples": len(t), "non_finite": bad}
            out["cursor"] = self.cursor(r, t)
        elif kind == Renderer.MULTILINE:
            t, cat = self._axis_values(a.axes[0]), a.axes[1]
            if cat.values is None:
                raise ResultUnavailable("category labels must be stored inline for rendering")
            cols = (a.storage.attrs or {}).get("columns") or \
                [(a.storage.column or 0) + k for k in range(len(cat.values))]
            series, bad = [], 0
            for label, c in zip(cat.values, cols):
                vals, b = _clean(self._read(a.storage, column=c))
                bad += b
                series.append({"label": str(label), "values": vals})
            if any(len(s["values"]) != len(t) for s in series):
                raise ResultUnavailable("category traces do not match the time axis length")
            out["data"] = {"time_ps": t, "series": series, "n_samples": len(t), "non_finite": bad}
            out["cursor"] = self.cursor(r, t)
        else:                                               # profile
            x = self._axis_values(a.axes[0])
            y, bad = _clean(self._read(a.storage))
            if x is None:
                x = list(range(len(y)))
            out["data"] = {"x": x, "x_label": a.axes[0].name, "x_kind": a.axes[0].kind,
                           "values": y, "n_samples": len(y), "non_finite": bad}
            out["cursor"] = {"mode": "none", "reason": "no time axis — not time-syncable"}
        with self._lock:
            self._cache[rid] = out
        return out

    def _axis_values(self, axis):
        if axis.values is not None:
            return list(axis.values)
        if axis.values_ref is not None:
            return self._read(axis.values_ref)
        return None

    # ── cursor → sample (exact only; never interpolated) ───────────────────
    def cursor(self, review_result, times: Optional[list]) -> dict:
        from analysis.review.sync import SyncMapping, frames_matching
        sync = review_result.sync
        if not sync.get("syncable"):
            return {"mode": "none", "reason": sync.get("reason", "not time-syncable")}
        tl = self.timeline
        if sync.get("mapping") == SyncMapping.ROW_IS_FRAME and tl.display_times is None:
            return {"mode": "identity", "rule": "sample i = display frame i (verified)"}
        if times is None:
            return {"mode": "none", "reason": "time values not stored"}
        display = tl.display_times if tl.display_times is not None else tl.source_times
        frame_samples = [-1] * len(display)
        shared = 0
        for i, t in enumerate(times):
            for f in frames_matching(display, t):
                if frame_samples[f] == -1:
                    frame_samples[f] = i
                else:
                    shared += 1
        return {"mode": "frames", "rule": "exact recorded time only (no interpolation)",
                "frame_samples": frame_samples, "frames_with_several_samples": shared}


def sample_for_frame(cursor: dict, frame: int) -> Optional[int]:
    """The sample shown at a display frame — the same rule the dashboard applies."""
    if cursor.get("mode") == "identity":
        return frame
    if cursor.get("mode") == "frames":
        fs = cursor["frame_samples"]
        s = fs[frame] if 0 <= frame < len(fs) else -1
        return s if s >= 0 else None
    return None


def session_payload(ds, timeline, catalog: ResultCatalog, viewer_status: dict) -> dict:
    """Everything the dashboard shows, from the dataset only (no data arrays)."""
    tl = ds.timeline
    display_times = timeline.display_times if timeline.display_times is not None \
        else timeline.source_times
    return {
        "session_id": ds.session_id, "schema_version": ds.schema_version, "status": ds.status,
        "system": {k: ds.system.get(k) for k in ("system_id", "condition_id", "replicate_id",
                                                  "membrane_present", "components")},
        "sources": {k: (ds.sources.get(k) or {}).get("path")
                    for k in ("trajectory", "topology", "structure")},
        "display_view": {k: ds.display_view.get(k) for k in (
            "view_ref", "kind", "purpose", "status", "reason", "request", "decisions")},
        "timeline": {"n_frames": timeline.n_frames, "start_time_ps": tl.get("start_time_ps"),
                     "end_time_ps": tl.get("end_time_ps"), "state": tl.get("state"),
                     "n_duplicate_groups": tl.get("n_duplicate_groups", 0),
                     "duplicate_groups": tl.get("duplicate_groups", []),
                     "display_relation": (tl.get("display") or {}).get("relation"),
                     "times_ps": display_times},
        "diagnostics": {k: ds.diagnostics.get(k) for k in (
            "execution", "status", "reason", "counts", "worst_severity", "detectors", "findings")},
        "annotations": [{k: a.get(k) for k in ("annotation_id", "kind", "state", "origin",
                                               "reasons", "used_by", "blocks")}
                        for a in ds.annotations],
        "observables": [{"instance_id": o.instance_id, "observable": o.observable,
                         "display_name": o.display_name, "state": o.state,
                         "availability": o.availability, "reason": o.reason,
                         "blocking_annotations": o.blocking_annotations,
                         "results": [rid for rid, (oo, _) in catalog._entries.items() if oo is o]}
                        for o in ds.observables],
        "results": [catalog.describe(rid) for rid in catalog.ids()],
        "blocked": [{"instance_id": o.instance_id, "state": o.state, "reason": o.reason,
                     "blocking_annotations": o.blocking_annotations}
                    for o in ds.observables if not o.available],
        "summary": ds.summary(),
        "viewer": viewer_status,
        # provenance only: why the session was configured this way (never drives rendering)
        "profile": ({k: (ds.provenance.get("profile") or {}).get(k) for k in (
            "id", "version", "source", "status", "definition_identity", "requirements",
            "optional", "overrides")} if (ds.provenance or {}).get("profile") else None),
    }
