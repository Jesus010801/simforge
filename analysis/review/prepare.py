"""Prepare a reproducible review session (headless; no viewer, no server).

    inputs → campaign discovery → one system → identity → (existing valid
    session? reuse it) → semantic index + diagnostics → display view (policy,
    purpose DISPLAY) → each requested observable instance through
    ``run_analyze`` (Phase 8 reuse first) → freeze referenced files in the
    content store → ReviewDataset → validate → atomic publish.

Layout under ``<analysis root>/review/`` (analysis root = ``--output`` or the
campaign default ``<study>/simforge_analysis``)::

    sessions/<session key>/review_dataset.json   the session (written last)
    instances/<instance key>/                    campaign output per observable instance
    systems/<system id>/                         semantic index, diagnostics, display views
    store/<sha[:2]>/<sha256><suffix>             immutable frozen results / timelines
    inputs/<name>/                               symlink workspace (TRAJECTORY mode)
    cache/                                       diagnostics probes + time index

Nothing in ``store/`` or a published session is ever rewritten; a later
recomputation only changes ``instances/`` — the session keeps its frozen copy.
"""
from __future__ import annotations

import datetime as _dt
import json
import os
import shutil
import tempfile
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Optional

from analysis.campaign.models import (
    AnalysisStatus, AnnotationState, DecisionClass, ObservablePurpose, StudyManifest,
    SystemRecord,
)
from analysis.review.dataset import (
    DATASET_FILENAME, REVIEW_SCHEMA_VERSION, Availability, CapabilityState, ObservableEntry,
    ReviewDataset, ReviewResult,
)
from analysis.review.request import ObservableRequest, ReviewRequest, request_problems
from analysis.review.store import STORE_DIRNAME, put_array, put_file

REVIEW_DIRNAME = "review"
SESSIONS_DIRNAME = "sessions"
_BUILD_PREFIX = ".building-"
DEFAULT_ANALYSIS_DIRNAME = "simforge_analysis"


class ReviewError(ValueError):
    """The session cannot be prepared at all (inputs, system choice)."""


@dataclass
class Preparation:
    dataset: ReviewDataset
    session_dir: Optional[Path] = None
    dataset_path: Optional[Path] = None
    reused_session: bool = False
    replaced: list[str] = field(default_factory=list)   # why an older session was re-prepared
    dry_run: bool = False

    def strict_failures(self) -> list[str]:
        out = [f"{o.instance_id}: {o.state} — {o.reason}" for o in self.dataset.observables
               if o.state not in (CapabilityState.AVAILABLE, CapabilityState.PLANNED)]
        dv = self.dataset.display_view
        if dv.get("status") not in (CapabilityState.AVAILABLE, CapabilityState.PLANNED):
            out.append(f"display view: {dv.get('status')} — {dv.get('reason', '')}")
        return out


# ═══════════════════════════════════════════════════════════════════════════════
# Inputs
# ═══════════════════════════════════════════════════════════════════════════════

def _link(dest: Path, src: Path) -> None:
    src = src.resolve()
    if dest.is_symlink() and Path(os.readlink(dest)) == src:
        return
    if dest.is_symlink() or dest.exists():
        dest.unlink()
    dest.symlink_to(src)


def resolve_inputs(target: str | Path, *, topology=None, structure=None, index=None,
                   output=None) -> tuple[Path, Path, dict]:
    """``(study root, analysis root, inputs record)``.

    RUN_DIR mode uses campaign discovery on the directory as is.  TRAJECTORY
    mode builds a symlink workspace (``review/inputs/<name>``) holding exactly
    the given files — explicit inputs are production by construction — plus a
    ``simforge_annotations.yaml`` found next to the trajectory.
    """
    from analysis.campaign.results import definition_token
    target = Path(target)
    if target.is_dir():
        if topology or structure or index:
            raise ReviewError("--topology/--structure/--index are for TRAJECTORY mode; "
                              "RUN_DIR mode discovers them (or use --manifest)")
        root = target.resolve()
        analysis_root = Path(output).resolve() if output else root / DEFAULT_ANALYSIS_DIRNAME
        return root, analysis_root, {"mode": "run_dir", "run_dir": str(root)}
    if not target.is_file():
        raise ReviewError(f"{target} is neither a run directory nor a trajectory file")
    if not topology:
        raise ReviewError("TRAJECTORY mode needs --topology (a .tpr is required by most observables)")
    files = {"trajectory": target, "topology": Path(topology)}
    if structure:
        files["structure"] = Path(structure)
    if index:
        files["index"] = Path(index)
    for role, p in files.items():
        if not p.is_file():
            raise ReviewError(f"--{role}: {p} does not exist")
    analysis_root = Path(output).resolve() if output else Path.cwd() / DEFAULT_ANALYSIS_DIRNAME
    key = definition_token({k: str(v.resolve()) for k, v in files.items()})[:8]
    ws = analysis_root / REVIEW_DIRNAME / "inputs" / f"{target.stem}-{key}"
    ws.mkdir(parents=True, exist_ok=True)
    _link(ws / f"md{target.suffix}", target)
    _link(ws / f"md{files['topology'].suffix}", files["topology"])
    if "structure" in files:
        _link(ws / f"md{files['structure'].suffix}", files["structure"])
    if "index" in files:
        _link(ws / "index.ndx", files["index"])
    ann = target.resolve().parent / "simforge_annotations.yaml"
    if ann.is_file():
        _link(ws / "simforge_annotations.yaml", ann)
    return ws, analysis_root, {"mode": "trajectory", "workspace": str(ws),
                               **{k: str(v.resolve()) for k, v in files.items()},
                               "annotations_file": str(ann) if ann.is_file() else None}


def _select_system(manifest: StudyManifest, system_id: Optional[str]) -> SystemRecord:
    ids = [s.system_id for s in manifest.systems]
    if system_id:
        rec = manifest.system(system_id)
        if rec is None:
            raise ReviewError(f"system '{system_id}' not found; discovered: {ids}")
        return rec
    if len(ids) != 1:
        raise ReviewError(f"{len(ids)} systems discovered — choose one with --system: {ids}"
                          if ids else "no simulation system was discovered")
    return manifest.systems[0]


# ═══════════════════════════════════════════════════════════════════════════════
# Identity
# ═══════════════════════════════════════════════════════════════════════════════

def session_identity_evidence(rec: SystemRecord, request: ReviewRequest,
                              profile_identity: Optional[str] = None) -> dict:
    """Scientific inputs only — never paths, output location, time or process.

    With a profile, its definition identity (content, not path or wording) is
    added so the session records exactly which profile version produced it."""
    from analysis.campaign.results import atom_set_hash
    from analysis.campaign.trajectory.policy import RULE_VERSION
    comps = sorted(
        [c.component_type, c.classification_state, c.selection, sorted(c.chain_ids),
         sorted(c.resnames), c.residue_count, c.atom_count,
         atom_set_hash(c.atom_ids) if c.atom_ids else None]
        for c in rec.components)
    return {
        "schema": REVIEW_SCHEMA_VERSION,
        "sources": {k: fp.digest for k, fp in sorted(rec.source_fingerprints.items())},
        "components": comps,
        "annotations": sorted([a.annotation_id, a.state, a.identity, a.atoms_sha256]
                              for a in rec.annotations),
        "display": {**request.display.to_dict(),
                    "requirements": request.display.requirements().cache_token()},
        "observables": sorted([o.observable, sorted(o.params.items())]
                              for o in request.observables),
        "policy_rule_version": RULE_VERSION,
        "window": None, "stride": None,
    } | ({"profile": profile_identity} if profile_identity else {})


def resolved_request_identity(request: ReviewRequest) -> str:
    from analysis.campaign.results import definition_token
    return definition_token({"observables": sorted([o.observable, sorted(o.params.items())]
                                                   for o in request.observables),
                             "display": request.display.to_dict()})


def instance_key(req: ObservableRequest) -> str:
    from analysis.campaign.results import definition_token
    return f"{req.observable}-" + definition_token(
        {"observable": req.observable, "parameters": req.params})[:12]


# ═══════════════════════════════════════════════════════════════════════════════
# Freezing references
# ═══════════════════════════════════════════════════════════════════════════════

class _Refs:
    """Store files + path policy relative to the session directory."""

    def __init__(self, review_root: Path, session_dir: Path, *, freeze: bool = True):
        self.review_root = review_root.resolve()
        self.store = self.review_root / STORE_DIRNAME
        self.session_dir = session_dir
        self.freeze = freeze
        self.external: set[str] = set()

    def rel(self, path: str | Path, kind: str) -> str:
        p = Path(path).resolve() if Path(path).exists() else Path(path)
        try:
            p.relative_to(self.review_root)
            return os.path.relpath(p, self.session_dir)
        except ValueError:
            self.external.add(kind)
            return str(p)

    def file(self, src: str | Path) -> Optional[dict]:
        if not self.freeze or not src or not Path(src).is_file():
            return None
        dest, sha = put_file(self.store, src)
        return {"path": os.path.relpath(dest, self.session_dir), "sha256": sha,
                "size_bytes": dest.stat().st_size}

    def array(self, values, dtype) -> Optional[dict]:
        if not self.freeze:
            return None
        dest, sha = put_array(self.store, values, dtype)
        return {"path": os.path.relpath(dest, self.session_dir), "sha256": sha}


def _freeze_result_array(arr, refs: _Refs):
    """A copy of ``arr`` whose storage / values_ref point at frozen store files."""
    from analysis.campaign.fingerprint import fingerprint_file
    frozen: dict[str, dict] = {}

    def repoint(ref):
        if ref is None:
            return None
        if ref.path not in frozen:
            frozen[ref.path] = refs.file(ref.path)
        f = frozen[ref.path]
        if f is None:
            return ref
        fp = fingerprint_file(refs.session_dir / f["path"], strong=True) \
            if (refs.session_dir / f["path"]).is_file() else None
        return replace(ref, path=f["path"], fingerprint=fp,
                       attrs={**(ref.attrs or {}), "sha256": f["sha256"]})

    from analysis.campaign.models import ResultArray
    out = ResultArray.from_dict(arr.to_dict())
    out.storage = repoint(out.storage)
    for ax in out.axes:
        ax.values_ref = repoint(ax.values_ref)
    return out


# ═══════════════════════════════════════════════════════════════════════════════
# Records
# ═══════════════════════════════════════════════════════════════════════════════

def _annotation_refs(params: dict) -> list[str]:
    return sorted({str(v).split(":", 1)[1] for v in params.values()
                   if isinstance(v, str) and v.startswith("annotation:")})


def _classify(ares, applicability: dict) -> tuple[str, str]:
    """(CapabilityState, Availability) from the campaign result — no new states."""
    st = ares.status if ares is not None else None
    if st == AnalysisStatus.SUCCESS:
        return CapabilityState.AVAILABLE, Availability.COMPUTED
    if st == AnalysisStatus.CACHED:
        user = bool((ares.reused_from or {}).get("user_selected"))
        return CapabilityState.AVAILABLE, (Availability.EXPLICIT_IMPORT if user
                                           else Availability.CACHED)
    if st == AnalysisStatus.PLANNED:
        return CapabilityState.PLANNED, Availability.PLANNED
    if st == AnalysisStatus.FAILED:
        return CapabilityState.FAILED, Availability.NONE
    if st == AnalysisStatus.ERROR:
        return CapabilityState.UNSUPPORTED, Availability.NONE
    reasons = applicability.get("reasons") or []
    if st == AnalysisStatus.SKIPPED or (
            reasons and all(": 0 components of this type" in r for r in reasons)):
        return CapabilityState.NOT_APPLICABLE, Availability.NONE
    text = " ".join([ares.message if ares else "", *reasons]).lower()
    if "ambiguous" in text:
        return CapabilityState.AMBIGUOUS, Availability.NONE
    if "unsupported" in text:
        return CapabilityState.UNSUPPORTED, Availability.NONE
    return CapabilityState.REVIEW_REQUIRED, Availability.NONE


def _view_record(view: dict, *, role: str, refs: _Refs) -> dict:
    """A ``TrajectoryView.to_dict()`` reference (paths per the portability policy)."""
    v = dict(view)
    rec = {"role": role, "view_ref": v.get("cache_key"), "kind": v.get("kind"),
           "purpose": (v.get("policy_ref") or {}).get("purpose"),
           "path": refs.rel(v["path"], "views") if v.get("path") else None,
           "output_fingerprint": v.get("output_fingerprint"),
           "time_index_ref": v.get("time_index_ref"),
           "decisions": v.get("decisions", []), "policy_ref": v.get("policy_ref"),
           "build_status": v.get("build_status"), "safe": v.get("safe"),
           "requirements": v.get("requirements")}
    return rec


def _diagnostics_record(report, ref: dict, refs: _Refs) -> dict:
    out = {"execution": ref.get("execution"), "status": ref.get("status"),
           "reason": ref.get("reason", ""), "report_identity": ref.get("report_identity")}
    if report is None:
        return out
    from analysis.campaign.models import Severity
    summ = report.summary()
    tmp = Path(tempfile.mkdtemp(prefix="simforge-diag-"))
    try:
        f = tmp / "trajectory_diagnostics.json"
        f.write_text(json.dumps(report.to_dict(), indent=2) + "\n")
        out["report_ref"] = refs.file(f)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    by = {}
    for r in report.detectors:
        key = ("ran_with_findings" if r.status == "ran" and r.n_diagnostics else
               "ran_found_nothing" if r.status == "ran" else r.status)
        by.setdefault(key, []).append({"detector": r.detector_id, "version": r.version,
                                       "reason": r.reason, "n_diagnostics": r.n_diagnostics})
    counts = summ["by_severity"]
    out.update({
        "counts": {"info": counts.get(Severity.INFO, 0), "warnings": counts.get(Severity.WARN, 0),
                   "review_required": counts.get(Severity.REVIEW, 0),
                   "errors": counts.get(Severity.ERROR, 0)},
        "worst_severity": summ["worst_severity"],
        "detectors": {k: by.get(k, []) for k in ("ran_with_findings", "ran_found_nothing",
                                                 "not_applicable", "deferred", "failed")},
        "findings": [{"code": d.code, "severity": d.severity, "scope": d.scope,
                      "message": d.message, "frames": d.frames, "times_ps": d.times_ps,
                      "detector": d.detector.get("id")} for d in report.diagnostics],
        "n_frames": report.n_frames,
    })
    return out


def _annotations_record(rec: SystemRecord, entries: list[ObservableEntry]) -> list[dict]:
    out = []
    for a in sorted(rec.annotations, key=lambda x: x.annotation_id):
        out.append({"annotation_id": a.annotation_id, "category": a.category, "kind": a.kind,
                    "origin": a.origin, "state": a.state, "identity": a.identity,
                    "atoms_sha256": a.atoms_sha256, "n_atoms": a.n_atoms,
                    "group_name": a.group_name, "numbering": a.numbering,
                    "reasons": a.reasons, "description": a.description,
                    "used_by": [e.instance_id for e in entries
                                if a.annotation_id in e.annotations_used],
                    "blocks": [e.instance_id for e in entries
                               if a.annotation_id in e.blocking_annotations]})
    declared = {a.annotation_id for a in rec.annotations}
    for e in entries:                     # referenced but never declared
        for aid in e.blocking_annotations:
            if aid not in declared and not any(x["annotation_id"] == aid for x in out):
                out.append({"annotation_id": aid, "state": "undeclared", "identity": None,
                            "reasons": ["not declared for this system"], "used_by": [],
                            "blocks": [x.instance_id for x in entries
                                       if aid in x.blocking_annotations]})
    return out


# ═══════════════════════════════════════════════════════════════════════════════
# Timeline
# ═══════════════════════════════════════════════════════════════════════════════

def _source_time_index(traj: str, tic: Path, gmx: str, *, allow_scan: bool):
    from analysis.campaign.gmx import gmx_version
    from analysis.campaign.trajectory.time_index import (
        GromacsTimeIndexBackend, index_trajectory, load_cached_index,
    )
    backend = GromacsTimeIndexBackend(gmx)
    idx = load_cached_index(traj, cache_dir=tic, backend_id=backend.id,
                            backend_version=gmx_version(gmx))
    if idx is not None or not allow_scan:
        return idx, idx is not None
    idx, _ = index_trajectory(traj, backend=backend, cache_dir=tic)
    return idx, False


def _timeline_record(idx, from_cache: bool, display_view, refs: _Refs) -> tuple[dict, Optional[list]]:
    if idx is None:
        return {"state": "not_indexed",
                "reason": "no cached time index; a trajectory pass is needed (not run in a dry run)"}, None
    rec = {"source_digest": idx.fingerprint.digest if idx.fingerprint else None,
           "backend": idx.backend, "backend_version": idx.backend_version,
           "schema_version": idx.schema_version, "n_frames": idx.n_frames,
           "start_time_ps": idx.start_time_ps, "end_time_ps": idx.end_time_ps,
           "state": idx.state, "strictly_increasing": idx.strictly_increasing,
           "non_decreasing": idx.non_decreasing,
           "duplicate_groups": idx.duplicate_groups[:50],
           "n_duplicate_groups": len(idx.duplicate_groups),
           "non_monotonic_steps": idx.non_monotonic_steps[:50],
           "time_resolution_ps": idx.time_resolution_ps,
           "warning_codes": [w.code for w in idx.warnings], "from_cache": from_cache,
           "times_ref": refs.array(idx.times_ps, "float64"),
           "steps_ref": refs.array([s if s is not None else -1 for s in idx.steps], "int64")}
    # display frames ↔ source frames
    disp = {"relation": "unknown"}
    if display_view is not None and display_view.safe:
        if display_view.kind == "raw":
            disp = {"relation": "identical_to_source",
                    "evidence": "raw passthrough of the source trajectory"}
        else:
            from analysis.campaign.trajectory.preprocessor import load_view_time_index
            vidx = load_view_time_index(display_view)
            if vidx is None:
                disp = {"relation": "unknown", "evidence": "display view has no time index"}
            elif vidx.same_timeline(idx):
                disp = {"relation": "identical_to_source",
                        "evidence": "frame-by-frame equal time index (Phase 1)"}
            else:
                disp = {"relation": "own_timeline", "n_frames": vidx.n_frames,
                        "times_ref": refs.array(vidx.times_ps, "float64")}
    rec["display"] = disp
    return rec, list(idx.times_ps)


# ═══════════════════════════════════════════════════════════════════════════════
# Main entry point
# ═══════════════════════════════════════════════════════════════════════════════

def prepare_review(
    target: str | Path,
    request: ReviewRequest,
    *,
    topology=None, structure=None, index=None, manifest: Optional[str] = None,
    system: Optional[str] = None, output=None, reuse: bool = True,
    reuse_from: tuple = (), cache_dir=None, dry_run: bool = False, force: bool = False,
    gmx: str = "gmx", profile=None, hide: tuple = (),
) -> Preparation:
    """``profile`` (a composed :class:`ReviewProfile`) is resolved against the
    selected system into ordinary observable instances, merged with the
    explicit ``request`` (user display wins; ``hide`` removes profile
    requests), then prepared exactly like any other request."""
    from analysis.campaign.manifest import build_manifest, load_manifest
    from analysis.campaign.results import definition_token
    from analysis.review.validate import load_review_dataset, validate_review_dataset

    study_root, analysis_root, inputs = resolve_inputs(
        target, topology=topology, structure=structure, index=index, output=output)
    mf = (load_manifest(manifest) if manifest
          else build_manifest(study_root, inspect_trajectories=False, gmx=gmx))
    rec = _select_system(mf, system)
    request, profile_ctx = _apply_profile(profile, rec, request, hide)
    evidence = session_identity_evidence(
        rec, request, profile.definition_identity if profile is not None else None)
    session_id = definition_token(evidence)
    review_root = analysis_root / REVIEW_DIRNAME
    sessions = review_root / SESSIONS_DIRNAME
    session_dir = sessions / session_id[:16]
    replaced: list[str] = []

    # ── an identical, still-valid session is reused as is ─────────────────────
    if not force and (session_dir / DATASET_FILENAME).is_file():
        try:
            old = load_review_dataset(session_dir)
            check = validate_review_dataset(old, session_dir)
            if old.session_id == session_id and check.valid:
                return Preparation(old, session_dir, session_dir / DATASET_FILENAME,
                                   reused_session=True, dry_run=dry_run)
            replaced = check.problems or ["session identity mismatch"]
        except ValueError as exc:
            replaced = [str(exc)]

    cache = Path(cache_dir).resolve() if cache_dir else review_root / "cache"
    scratch = Path(tempfile.mkdtemp(prefix="simforge-review-dry-")) if dry_run else None
    work = scratch if dry_run else review_root
    if not dry_run:
        sessions.mkdir(parents=True, exist_ok=True)
    build = (Path(tempfile.mkdtemp(prefix=f"{_BUILD_PREFIX}{session_id[:16]}-{os.getpid()}-",
                                   dir=sessions)) if not dry_run else scratch / "session")
    build.mkdir(parents=True, exist_ok=True)
    try:
        ds = _assemble(rec, mf, request, evidence, session_id, study_root, analysis_root,
                       inputs, work, review_root, build, cache, reuse, list(reuse_from),
                       dry_run, force, gmx, profile_ctx)
        if dry_run:
            return Preparation(ds, None, None, replaced=replaced, dry_run=True)
        (build / DATASET_FILENAME).write_text(json.dumps(ds.to_dict(), indent=2) + "\n")
        check = validate_review_dataset(ds, build)
        if not check.valid:
            # never publish a dataset pointing at inconsistent data
            (build / DATASET_FILENAME).unlink()
            raise ReviewError("prepared session failed validation: " + "; ".join(check.problems))
        _publish(build, session_dir)
        return Preparation(ds, session_dir, session_dir / DATASET_FILENAME, replaced=replaced)
    finally:
        if scratch is not None:
            shutil.rmtree(scratch, ignore_errors=True)
        elif build.exists():
            shutil.rmtree(build, ignore_errors=True)


def _publish(build: Path, final: Path) -> None:
    """Atomic: the dataset only ever appears complete at its final path."""
    stale = None
    if final.exists():
        stale = final.with_name(f".stale-{final.name}-{os.getpid()}")
        os.rename(final, stale)
    os.rename(build, final)
    if stale is not None:
        shutil.rmtree(stale, ignore_errors=True)


def _apply_profile(profile, rec, request: ReviewRequest, hide) -> tuple:
    """Profile data → the same ReviewRequest type the generic path prepares."""
    from dataclasses import replace as _replace
    from analysis.review.request import DisplayRequest
    if profile is None:
        if hide:
            raise ReviewError("--hide removes profile requests; it needs --profile")
        if request.display is None:
            request = _replace(request, display=DisplayRequest())
        return request, None
    from analysis.review.profiles import ProfileError, ProfileStatus, merge_request, resolve_profile
    resolution = resolve_profile(profile, rec)
    if resolution.status == ProfileStatus.NOT_APPLICABLE:
        raise ReviewError(f"profile {profile.id!r} is not applicable to {rec.system_id}: "
                          + "; ".join(resolution.reasons))
    try:
        merged, overrides = merge_request(resolution, request, hide)
    except ProfileError as exc:
        raise ReviewError(str(exc)) from exc
    return merged, {"resolution": resolution, "overrides": overrides}


def _profile_record(ctx, request: ReviewRequest, entries) -> Optional[dict]:
    """Why this session was configured this way (provenance; not the definition
    of any observable)."""
    if ctx is None:
        return None
    from analysis.review.profiles import PROFILE_SCHEMA, final_status
    res, p = ctx["resolution"], ctx["resolution"].profile
    return {"id": p.id, "version": p.version, "schema": PROFILE_SCHEMA, "source": p.source,
            "description": p.description, "definition_identity": p.definition_identity,
            "composition": p.composition, "duplicates_dropped": p.duplicates_dropped,
            "status": final_status(res, entries), "requirements": res.requirements,
            "requests": [r.instance_id for r in res.requests], "optional": res.optional,
            "display": p.display, "overrides": ctx["overrides"],
            "resolved_request_identity": resolved_request_identity(request),
            "policy_intent": "none — a profile never authorises transformations; its display "
                             "preference is judged like an unflagged request"}


def _assemble(rec, mf, request, evidence, session_id, study_root, analysis_root, inputs,
              work, review_root, build, cache, reuse, reuse_from, dry_run, force, gmx,
              profile_ctx=None) -> ReviewDataset:
    from analysis.campaign.compatibility import RESOLVER_VERSION
    from analysis.campaign.gmx import gmx_version
    from analysis.campaign.manifest import write_manifest
    from analysis.campaign.observables import registry
    from analysis.campaign.orchestration.study_analyzer import (
        plan_view, prepare_system, time_index_cache_dir,
    )
    from analysis.campaign.trajectory.policy import POLICY_ENGINE, RULE_VERSION
    registry.ensure_loaded()
    refs = _Refs(review_root, build, freeze=not dry_run)
    warnings: list[str] = []

    # semantic index + diagnostics once, shared by display and observables
    sys_dir = work / "systems" / rec.system_id
    report, diag_ref, w = prepare_system(rec, sys_dir, gmx=gmx, cache_dir=cache, dry_run=dry_run)
    warnings += [x.message for x in w]
    filtered = replace(mf, systems=[rec],
                       ambiguities=[a for a in mf.ambiguities if a.system_id == rec.system_id])
    manifest_path = write_manifest(filtered, build / "inputs")["yaml"]

    # display view: its own plan (purpose DISPLAY), never an analysis view
    dview, dw = plan_view(rec, request.display.requirements(), sys_dir,
                          purpose=ObservablePurpose.DISPLAY,
                          intent=request.display.policy_intent(), diagnostics=report,
                          diagnostics_ref=diag_ref, gmx=gmx, dry_run=dry_run)
    display = _display_record(dview, request, refs, dry_run)

    # observable instances
    entries: list[ObservableEntry] = []
    views: dict[str, dict] = {}
    for req in request.observables:
        entries.append(_run_instance(req, rec, study_root, analysis_root, review_root, work,
                                     manifest_path, cache, reuse, reuse_from, dry_run, force,
                                     gmx, refs, views))

    # timeline (after the observables: a fresh diagnostics pass may have indexed it)
    tic = time_index_cache_dir(sys_dir, diag_ref, cache)
    traj = rec.trajectory_paths[0] if len(rec.trajectory_paths) == 1 else None
    idx, from_cache = (_source_time_index(traj, tic, gmx, allow_scan=not dry_run)
                       if traj else (None, False))
    timeline, times = _timeline_record(idx, from_cache, dview if display["status"] ==
                                       CapabilityState.AVAILABLE else None, refs)
    if not traj:
        timeline = {"state": "unavailable",
                    "reason": f"{len(rec.trajectory_paths)} production segments; one is required"}
    digest = timeline.get("source_digest")
    for e in entries:
        for r in e.results:
            r.sync = _sync(r, times, digest)
        for r in e.results:
            r.array = _freeze_result_array(r.array, refs) if not dry_run else r.array

    fps = rec.source_fingerprints

    def src(role, path):
        fp = fps.get(role)
        return ({"path": str(Path(path).resolve()), "fingerprint": fp.to_dict() if fp else None}
                if path else None)

    ds = ReviewDataset(
        session_id=session_id, identity_evidence=evidence,
        status="dry_run" if dry_run else "prepared",
        system={"system_id": rec.system_id, "condition_id": rec.condition_id,
                "replicate_id": rec.replicate_id, "membrane_present": rec.membrane_present,
                "components": [{"type": c.component_type, "label": c.label,
                                "state": c.classification_state, "selection": c.selection,
                                "chain_ids": c.chain_ids, "resnames": c.resnames,
                                "residue_count": c.residue_count, "atom_count": c.atom_count}
                               for c in rec.components],
                "semantic_index": ({"groups": [g.to_dict() for g in rec.semantic_index.groups],
                                    "ndx_ref": refs.file(rec.semantic_index.path)}
                                   if rec.semantic_index else None),
                "manifest_ref": refs.file(manifest_path)},
        sources={"inputs": inputs,
                 "trajectory": src("trajectory", traj),
                 "topology": src("topology", rec.topology_path),
                 "structure": src("structure", rec.structure_path),
                 "index": src("index", rec.index_path)},
        display_view=display,
        analysis_views=sorted(views.values(), key=lambda v: v["view_ref"] or ""),
        timeline=timeline,
        annotations=_annotations_record(rec, entries),
        diagnostics=_diagnostics_record(report, diag_ref, refs),
        requested=[r.to_dict() for r in request.observables],
        observables=entries,
        provenance={"engine": "simforge/review-preparation/v1",
                    "simforge_version": mf.simforge_version, "gmx_version": gmx_version(gmx),
                    "compatibility_resolver": RESOLVER_VERSION, "policy_engine": POLICY_ENGINE,
                    "policy_rule_version": RULE_VERSION,
                    "profile": _profile_record(profile_ctx, request, entries),
                    "preparation": {"reuse": reuse, "reuse_from": [str(p) for p in reuse_from],
                                    "dry_run": dry_run, "force": force},
                    "prepared_at_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(
                        timespec="seconds"),
                    "note": "prepared_at_utc is informational; it is not part of the identity"},
        warnings=warnings + [f"display: {x.message}" for x in dw],
    )
    ds.portability = {
        "within_review_root": "results, provenance, timeline, diagnostics and manifest are "
                              "frozen in the store and referenced relative to the session",
        "external": sorted(refs.external | {"sources"}),
        "portable_without_sources": True,
        "note": "sources (and views outside the review root) are absolute paths + fingerprints; "
                "moving them invalidates the session until re-prepared",
    }
    return ds


def _display_record(view, request: ReviewRequest, refs: _Refs, dry_run: bool) -> dict:
    decisions = view.decisions
    if view.safe and view.build_status in ("raw", "built", "reused", "planned"):
        status = CapabilityState.PLANNED if view.build_status == "planned" else CapabilityState.AVAILABLE
        reason = ""
    else:
        blocking = [d for d in decisions if not d["applied"]]
        classes = {d["classification"] for d in blocking}
        status = (CapabilityState.REVIEW_REQUIRED if DecisionClass.VALID_WITH_INTENT in classes
                  else CapabilityState.UNSUPPORTED if DecisionClass.UNSUPPORTED in classes
                  else CapabilityState.REFUSED if DecisionClass.REFUSED in classes
                  else CapabilityState.FAILED)
        reason = "; ".join(f"{d['operation']} = {d['classification']} ({d['rule_id']}): "
                           f"{d['reason']}" for d in blocking) or \
            "; ".join(w.message for w in view.warnings)
        if status == CapabilityState.REVIEW_REQUIRED:
            reason += " — pass --display-intent to choose it explicitly"
    rec = _view_record(view.to_dict(), role="display", refs=refs)
    rec.update({"request": request.display.to_dict(), "status": status, "reason": reason,
                "purpose": ObservablePurpose.DISPLAY})
    return rec


def _sync(result: ReviewResult, times, digest) -> dict:
    from analysis.review.sync import sync_metadata
    return sync_metadata(result.array, times, digest)


def _run_instance(req, rec, study_root, analysis_root, review_root, work, manifest_path,
                  cache, reuse, reuse_from, dry_run, force, gmx, refs, views) -> ObservableEntry:
    from analysis.campaign.compatibility import definition_identity
    from analysis.campaign.observables import registry
    from analysis.campaign.orchestration.study_analyzer import run_analyze

    entry = ObservableEntry(instance_id=req.instance_id, observable=req.observable,
                            parameters=req.params,
                            annotations_used=_annotation_refs(req.params))
    problems = request_problems(req)
    if problems:
        entry.state, entry.reason = CapabilityState.UNSUPPORTED, "; ".join(problems)
        return entry
    spec = registry.get(req.observable)
    entry.display_name, entry.purpose = spec.display_name, spec.purpose
    params = {req.observable: req.params}
    appl = spec.applicability(rec, params, rec.semantic_index)
    entry.applicability = appl
    usable = {a.annotation_id for a in rec.annotations if a.state == AnnotationState.ACTIVE}
    entry.blocking_annotations = [a for a in entry.annotations_used if a not in usable]

    key = instance_key(req)
    real_inst = review_root / "instances" / key
    inst = (work / "instances" / key) if dry_run else real_inst
    roots = ([real_inst] if dry_run else []) + [analysis_root, *map(Path, reuse_from)]
    res = run_analyze(study_root, [req.observable], output_dir=inst, manifest_path=manifest_path,
                      # --force re-prepares the *session*; validated views and compatible
                      # results are still reused (--no-reuse is the way to recompute)
                      gmx=gmx, inspect_trajectories=False, dry_run=dry_run,
                      parameters=params, diagnostics="auto", diagnostics_cache_dir=cache,
                      reuse=reuse, reuse_roots=roots)
    ares = next((r for r in res.results if r.system_id == rec.system_id), None)
    if ares is None:
        entry.state = CapabilityState.UNSUPPORTED
        entry.reason = "; ".join(w.message for w in res.warnings) or "no result produced"
        return entry
    entry.state, entry.availability = _classify(ares, appl)
    entry.execution_status = ares.status
    entry.reason = ares.message if entry.state != CapabilityState.AVAILABLE else ""
    entry.definition_token = ares.definition_token
    entry.definition_identity = definition_identity(ares.definition_evidence)
    if entry.availability == Availability.EXPLICIT_IMPORT:
        entry.externally_supplied, entry.independently_verified = True, False
    comp = ares.compatibility or {}
    if comp:
        entry.reuse = {"decision": comp.get("decision"),
                       "reused_from": ares.reused_from or comp.get("reuse_from"),
                       "chosen": _compat_summary(comp.get("chosen")),
                       "evaluated": [_compat_summary(e) for e in comp.get("evaluated", [])]}
        if not dry_run:
            entry.reuse["report_ref"] = _json_ref(comp, refs)
    if entry.state != CapabilityState.AVAILABLE:
        return entry
    from analysis.campaign.compatibility import storage_integrity
    ok, why = storage_integrity(ares)
    if not ok:                                  # never freeze stale science
        entry.state, entry.availability = CapabilityState.FAILED, Availability.NONE
        entry.reason = f"result files not usable: {why}"
        return entry
    prov = {}
    if ares.provenance_path and Path(ares.provenance_path).is_file():
        prov = json.loads(Path(ares.provenance_path).read_text())
        entry.provenance_ref = refs.file(ares.provenance_path)
    tv = prov.get("trajectory_view")
    if tv and tv.get("cache_key"):
        views.setdefault(tv["cache_key"], _view_record(tv, role="analysis", refs=refs))
        views[tv["cache_key"]].setdefault("used_by", []).append(entry.instance_id)
    for arr in ares.arrays:
        entry.results.append(ReviewResult(
            array=arr, origin={"storage": arr.storage.to_dict() if arr.storage else None,
                               "provenance": ares.provenance_path}))
    entry.view_ref = next((a.view_ref for a in ares.arrays if a.view_ref), None)
    return entry


def _compat_summary(c: Optional[dict]) -> Optional[dict]:
    if not c:
        return None
    return {"status": c.get("status"), "reason": c.get("reason"),
            "candidate": (c.get("candidate") or {}).get("source"),
            "tier": (c.get("candidate") or {}).get("tier"),
            "checks": {x["dimension"]: x["outcome"] for x in c.get("checks", [])}}


def _json_ref(obj: dict, refs: _Refs) -> Optional[dict]:
    tmp = Path(tempfile.mkdtemp(prefix="simforge-review-"))
    try:
        f = tmp / "compatibility.json"
        f.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")
        return refs.file(f)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
