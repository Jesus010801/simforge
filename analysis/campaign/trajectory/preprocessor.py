"""Resolve :class:`TrajectoryRequirements` into a concrete :class:`TrajectoryView`.

One derived trajectory per *distinct* preprocessing definition.  This module
records **what exactly was done to the coordinates**; it does not decide
whether doing it was scientifically appropriate.

Preprocessing pipeline (only the steps a requirement asks for, in this order —
unchanged from the original implementation):

1. ``-pbc whole``   (make_whole)   — molecules made whole
2. ``-pbc nojump``  (nojump)       — unwrapped, diffusion-safe
3. ``-center``      (center)       — a semantic group centred in the box
4. ``-fit rot+trans`` (fit)        — least-squares fit to a semantic group

Each step is a separate ``gmx trjconv`` call fed the relevant *named* group
from the semantic index (never a numeric id).  If a required semantic group is
missing the view is marked ``safe=False`` and the caller must fail the analysis
rather than fall back to raw coordinates.

View identity (``TrajectoryView.cache_key``)
--------------------------------------------
A deterministic hash (``results.definition_token``) of *identity evidence*:

* content digests of the source trajectory segment(s) and of every file passed
  as ``-s`` (topology / reference coordinates);
* for every index group a step actually selects via ``-n``: name + atom-set
  hash (+ multiplicity hash when the group repeats atoms) — the ``.ndx`` file
  as a whole is *not* part of the identity, so editing an unused group keeps
  the view reusable while changing a used group's atoms does not;
* the ordered operations with their arguments, input/output paths replaced by
  role placeholders — so operation order and any timeline-changing argument
  (``-b -e -dt -t0 -settime -timestep``) are part of the identity;
* the requirement token and the GROMACS version.

Output / work directories are never part of the identity.

Cache entries
-------------
``<work_dir>/<cache_key>/`` holds exactly ``<kind>.xtc`` and
``view_manifest.json``.  Builds happen in ``<work_dir>/.building-<key>-<pid>-*``;
intermediates are deleted, the derived file is fingerprinted and time-indexed
(Phase 1), the manifest is written **last**, and the directory is renamed into
place atomically.  An entry is reused only if its manifest validates (schema,
identity, source digests, output digest); otherwise it is regenerated and the
rejection reason recorded on the view.  Entries without a manifest (pre-Phase-3
caches) are never trusted; legacy entries under old keys are left untouched.
"""
from __future__ import annotations

import datetime
import json
import os
import shutil
import tempfile
from pathlib import Path
from typing import Optional

from analysis.campaign.fingerprint import fingerprint_file
from analysis.campaign.gmx import gmx_available, gmx_version, run_gmx
from analysis.campaign.models import (
    CampaignWarning, PreprocessingOperation, SemanticIndex, Severity,
    SourceFingerprint, TrajectoryRequirements, TrajectoryView, ValidationState,
    ViewBuildStatus, ViewCacheStatus, ViewKind,
)

VIEW_MANIFEST_SCHEMA = "simforge/trajectory-view-manifest/v1"
VIEW_IDENTITY_SCHEMA = "simforge/trajectory-view-identity/v1"
MANIFEST_NAME = "view_manifest.json"
TIME_INDEX_DIRNAME = "time_index"
_BUILD_PREFIX = ".building-"

#: trjconv arguments that change the frame set or the time axis.
TIMELINE_ARGS = ("-b", "-e", "-dt", "-t0", "-settime", "-timestep", "-tu", "-skip", "-drop")
_PATH_ARGS = {"-f": "input", "-o": "output", "-s": "structure", "-n": "index"}


def _resolve_group(index: Optional[SemanticIndex], name: str) -> Optional[str]:
    if not index:
        return None
    return name if name in index.group_names() else None


# ═══════════════════════════════════════════════════════════════════════════════
# Planning (unchanged command construction)
# ═══════════════════════════════════════════════════════════════════════════════

def _plan_steps(requirements, src, topology_path, structure_path, semantic_index,
                out_dir: Path):
    """Return ``(steps, error_warning)``.  Commands are built exactly as the
    original implementation built them; ``out_dir`` only sets where outputs go."""
    ref_for_center = structure_path or topology_path
    current = src
    steps: list[dict] = []

    def add(op, args, stdin, reason, validation, uses_index):
        steps.append({"op": op, "args": args, "stdin": stdin, "reason": reason,
                      "validation": validation, "uses_index": uses_index})

    if requirements.requires_whole_molecules:
        out = out_dir / "whole.xtc"
        add("make_whole",
            ["trjconv", "-s", topology_path, "-f", current, "-o", str(out), "-pbc", "whole"],
            "System", "make molecules whole across periodic boundaries",
            "required before any fitting or geometric measurement", False)
        current = str(out)

    if requirements.requires_nojump:
        out = out_dir / "nojump.xtc"
        add("nojump",
            ["trjconv", "-s", topology_path, "-f", current, "-o", str(out), "-pbc", "nojump"],
            "System", "remove periodic jumps (unwrap) — preserves diffusion",
            "diffusion/MSD-safe unwrapping; no spatial fitting applied", False)
        current = str(out)

    if requirements.centering_target:
        grp = _resolve_group(semantic_index, requirements.centering_target)
        if grp is None:
            return None, CampaignWarning(
                "missing_semantic_group",
                f"centering target '{requirements.centering_target}' is not in the "
                f"semantic index; refusing to guess", Severity.ERROR)
        out = out_dir / "centered.xtc"
        add("center",
            ["trjconv", "-s", ref_for_center, "-f", current, "-o", str(out),
             "-pbc", "mol", "-center"],
            f"{grp}\nSystem", f"centre '{grp}' in the box",
            "centering only; no rotational fitting", True)
        current = str(out)

    if requirements.fit_selection:
        grp = _resolve_group(semantic_index, requirements.fit_selection)
        if grp is None:
            return None, CampaignWarning(
                "missing_semantic_group",
                f"fit selection '{requirements.fit_selection}' is not in the "
                f"semantic index; refusing to guess", Severity.ERROR)
        out = out_dir / "fitted.xtc"
        add("fit",
            ["trjconv", "-s", ref_for_center, "-f", current, "-o", str(out),
             "-fit", "rot+trans"],
            f"{grp}\nSystem", f"least-squares rot+trans fit to '{grp}'",
            f"reference frame = frame 0 of '{grp}'", True)
        current = str(out)

    index_args: list[str] = []
    if semantic_index and semantic_index.path:
        index_args = ["-n", semantic_index.path]
    for st in steps:
        st["full_args"] = list(st["args"]) + (index_args if st["uses_index"] else [])
    return steps, None


# ═══════════════════════════════════════════════════════════════════════════════
# Identity
# ═══════════════════════════════════════════════════════════════════════════════

def _group_evidence(index_path: Optional[str], name: str, uses_index: bool) -> dict:
    if not uses_index:
        return {"name": name, "source": "default groups of the -s file"}
    from analysis.campaign.results import index_group_evidence
    ev = index_group_evidence(index_path, name) if index_path else None
    return ev if ev is not None else {"name": name, "source": "not found in index"}


def _identity_evidence(requirements, kind, steps, src_fps, fp_cache, index_path,
                       gmx_ver) -> tuple[dict, list[dict]]:
    ops: list[dict] = []
    groups: list[dict] = []
    for st in steps:
        args = st["full_args"]
        templ: list[str] = []
        i = 0
        while i < len(args):
            a = args[i]
            templ.append(a)
            if a in _PATH_ARGS and i + 1 < len(args):
                role = _PATH_ARGS[a]
                val = args[i + 1]
                if role == "structure":
                    templ.append(f"<structure:{fp_cache[val].digest}>")
                elif role == "output":
                    templ.append(f"<output:{Path(val).name}>")
                else:
                    templ.append(f"<{role}>")
                i += 2
                continue
            i += 1
        sel = [g for g in (st["stdin"] or "").split("\n") if g]
        gev = [_group_evidence(index_path, g, st["uses_index"]) for g in sel]
        groups.extend(gev)
        ops.append({"operation": st["op"], "args": templ, "stdin": st["stdin"],
                    "groups": gev})
    evidence = {
        "schema": VIEW_IDENTITY_SCHEMA,
        "kind": kind,
        "requirements": requirements.cache_token(),
        "sources": [fp.digest for fp in src_fps],
        "operations": ops,
        "backend": {"tool": "gmx trjconv", "version": gmx_ver} if steps else None,
    }
    return evidence, groups


def view_identity(evidence: dict) -> str:
    from analysis.campaign.results import definition_token
    return definition_token(evidence)


def _timeline_args(steps) -> list[list[str]]:
    out: list[list[str]] = []
    for st in steps:
        a = st["full_args"]
        for i, tok in enumerate(a):
            if tok in TIMELINE_ARGS:
                out.append([st["op"], tok] + ([a[i + 1]] if i + 1 < len(a) else []))
    return out


# ═══════════════════════════════════════════════════════════════════════════════
# Manifest validation
# ═══════════════════════════════════════════════════════════════════════════════

def read_manifest(entry: Path) -> Optional[dict]:
    try:
        return json.loads((Path(entry) / MANIFEST_NAME).read_text())
    except (OSError, ValueError):
        return None


def validate_entry(entry: Path, *, identity: str, source_digests: list[str],
                   kind: str) -> tuple[bool, str]:
    """``(ok, reason)`` — whether a cache entry may be reused."""
    entry = Path(entry)
    if not entry.is_dir():
        return False, "no cache entry"
    mf = entry / MANIFEST_NAME
    if not mf.is_file():
        return False, "manifest missing (unverified / pre-manifest cache entry)"
    m = read_manifest(entry)
    if m is None:
        return False, "manifest unreadable or corrupt"
    if m.get("schema_version") != VIEW_MANIFEST_SCHEMA:
        return False, f"manifest schema {m.get('schema_version')!r} not supported"
    if m.get("status") != "complete":
        return False, f"manifest status {m.get('status')!r}"
    if m.get("view_identity") != identity:
        return False, "manifest identity does not match the requested view"
    recorded = [f.get("digest") for f in (m.get("sources") or {}).get("trajectories", [])]
    if recorded != list(source_digests):
        return False, "source trajectory fingerprint changed"
    out = (m.get("output") or {})
    rel = out.get("relative_path")
    if rel != f"{kind}.xtc":
        return False, "manifest output path inconsistent"
    out_path = entry / rel
    if not out_path.is_file():
        return False, "derived trajectory missing"
    fp = out.get("fingerprint") or {}
    strong = fp.get("mode") == "strong"
    if fingerprint_file(out_path, strong=strong).digest != fp.get("digest"):
        return False, "derived trajectory changed after the manifest was written"
    return True, "valid"


def _pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _sweep_stale_builds(work_dir: Path, key: str) -> None:
    """Remove abandoned build dirs *for this key only* whose process is gone."""
    for d in work_dir.glob(f"{_BUILD_PREFIX}{key}-*"):
        try:
            pid = int(d.name[len(_BUILD_PREFIX) + len(key) + 1:].split("-", 1)[0])
        except ValueError:
            continue
        if d.is_dir() and not d.is_symlink() and not _pid_alive(pid):
            shutil.rmtree(d, ignore_errors=True)


def _remove_entry(work_dir: Path, entry: Path) -> None:
    """Delete one invalid cache entry — only a real directory directly under work_dir."""
    if entry.parent.resolve() == work_dir.resolve() and entry.is_dir() and not entry.is_symlink():
        shutil.rmtree(entry)


# ═══════════════════════════════════════════════════════════════════════════════
# Time index of a view (Phase 1)
# ═══════════════════════════════════════════════════════════════════════════════

def _time_index_ref(idx, cache_dir: Path, same_as_source: Optional[bool]) -> dict:
    return {
        "backend": idx.backend,
        "backend_version": idx.backend_version,
        "schema_version": idx.schema_version,
        "trajectory_digest": idx.fingerprint.digest if idx.fingerprint else None,
        "cache_dir": str(cache_dir.resolve()),
        "n_frames": idx.n_frames,
        "start_time_ps": idx.start_time_ps,
        "end_time_ps": idx.end_time_ps,
        "state": idx.state,
        "warning_codes": [w.code for w in idx.warnings],
        "same_timeline_as_source": same_as_source,
    }


def load_view_time_index(view: TrajectoryView):
    """Return the view's authoritative ``FrameTimeIndex`` from the Phase 1 cache,
    identified by the digest/backend/version recorded on the view.  No backend
    is run; None if it is not retrievable (RAW views record none)."""
    from analysis.campaign.trajectory.time_index import load_cached_index
    ref = view.time_index_ref
    if not ref or not view.path:
        return None
    idx = load_cached_index(view.path, cache_dir=ref["cache_dir"],
                            backend_id=ref["backend"], backend_version=ref["backend_version"])
    if idx is None or idx.fingerprint is None or idx.fingerprint.digest != ref["trajectory_digest"]:
        return None
    return idx


# ═══════════════════════════════════════════════════════════════════════════════
# build_view
# ═══════════════════════════════════════════════════════════════════════════════

def build_view(
    *,
    requirements: TrajectoryRequirements,
    trajectory_paths: list[str],
    topology_path: str,
    structure_path: Optional[str],
    semantic_index: Optional[SemanticIndex],
    source_fingerprints: list[SourceFingerprint],
    topology_fingerprint: Optional[SourceFingerprint],
    work_dir: Path,
    gmx: str = "gmx",
    force: bool = False,
    dry_run: bool = False,
    time_index_cache_dir: Optional[Path] = None,
    decisions: Optional[list] = None,
    policy_ref: Optional[dict] = None,
) -> TrajectoryView:
    """Execute the requested operations.  Scientific admissibility is decided
    beforehand by ``trajectory.policy``; pass its decisions to mark the view as
    policy-planned.  Direct calls (``decisions=None``) stay legacy views."""
    kind = requirements.view_kind()
    view = TrajectoryView(
        kind=kind, requirements=requirements, topology_path=topology_path,
        source_fingerprints=source_fingerprints,
    )
    view.policy_planned = decisions is not None
    view.policy_ref = policy_ref
    view.decisions = [d.to_dict() if hasattr(d, "to_dict") else dict(d)
                      for d in (decisions or [])]
    if decisions is not None and any(not dd["applied"] for dd in view.decisions):
        view.safe = False
        view.build_status = ViewBuildStatus.FAILED
        view.warnings.append(CampaignWarning(
            "policy_blocked", "refusing to execute: not every requested operation was "
            "admitted by the preprocessing policy", Severity.ERROR))
        return view
    if requirements.fit_mode != "rot+trans" or requirements.reconstruction:
        view.safe = False
        view.build_status = ViewBuildStatus.FAILED
        view.warnings.append(CampaignWarning(
            "unsupported_operation",
            f"no executor for fit_mode={requirements.fit_mode!r} / "
            f"reconstruction={requirements.reconstruction!r}; nothing was built", Severity.ERROR))
        return view

    # Current content fingerprints of the sources (discovery-time ones may be stale).
    src_fps: list[SourceFingerprint] = []
    for i, p in enumerate(trajectory_paths):
        if Path(p).is_file():
            fp = fingerprint_file(Path(p))
            src_fps.append(fp)
            if i < len(source_fingerprints) and source_fingerprints[i].digest != fp.digest:
                view.warnings.append(CampaignWarning(
                    "source_changed_since_discovery",
                    f"{Path(p).name} changed since it was discovered; using its current "
                    f"content", Severity.WARN))
    if topology_fingerprint and topology_path and Path(topology_path).is_file():
        if fingerprint_file(Path(topology_path)).digest != topology_fingerprint.digest:
            view.warnings.append(CampaignWarning(
                "source_changed_since_discovery",
                f"{Path(topology_path).name} changed since it was discovered", Severity.WARN))

    if kind == ViewKind.RAW:
        view.build_status = ViewBuildStatus.RAW
        view.cache_status = ViewCacheStatus.BYPASSED
        evidence, _ = _identity_evidence(requirements, kind, [], src_fps, {}, None, None)
        view.cache_key = view_identity(evidence)
        if len(trajectory_paths) == 1:
            view.path = trajectory_paths[0]
        view.safe = len(trajectory_paths) == 1
        if not view.safe:
            view.warnings.append(CampaignWarning(
                "multi_segment_raw",
                "multiple trajectory segments and no preprocessing requested; "
                "an observable must handle segments explicitly",
                Severity.WARN,
            ))
        return view

    def fail(code, message, severity=Severity.ERROR):
        view.safe = False
        view.path = None
        view.build_status = ViewBuildStatus.FAILED
        view.warnings.append(CampaignWarning(code, message, severity))
        return view

    if not gmx_available(gmx):
        return fail("gmx_unavailable",
                    f"'{gmx}' not available; cannot build the '{kind}' trajectory view")

    if len(trajectory_paths) != 1:
        return fail("segments_not_concatenated",
                    "trajectory has multiple segments; automatic concatenation is not "
                    "performed. Provide a single trajectory or a validated concatenation.")

    src = trajectory_paths[0]
    work_dir = Path(work_dir)
    index_path = semantic_index.path if semantic_index and semantic_index.path else None

    # Plan against the entry location first (identity + dry run use it).
    steps, err = _plan_steps(requirements, src, topology_path, structure_path,
                             semantic_index, work_dir / "__entry__")
    if err is not None:
        view.safe = False
        view.build_status = ViewBuildStatus.FAILED
        view.warnings.append(err)
        return view

    fp_cache: dict[str, SourceFingerprint] = {}
    for st in steps:
        a = st["full_args"]
        s_val = a[a.index("-s") + 1]
        if s_val not in fp_cache:
            if not s_val or not Path(s_val).is_file():
                return fail("structure_missing", f"-s input {s_val!r} does not exist")
            fp_cache[s_val] = fingerprint_file(Path(s_val))
    gmx_ver = gmx_version(gmx)
    evidence, groups = _identity_evidence(requirements, kind, steps, src_fps, fp_cache,
                                          index_path, gmx_ver)
    key = view_identity(evidence)
    view.cache_key = key
    entry = work_dir / key
    final_name = f"{kind}.xtc"
    ti_cache = Path(time_index_cache_dir) if time_index_cache_dir else work_dir / TIME_INDEX_DIRNAME

    if dry_run:
        steps, _ = _plan_steps(requirements, src, topology_path, structure_path,
                               semantic_index, entry)
        for st in steps:
            args = st["full_args"]
            view.operations.append(_op_record(st, args, None, dry=True))
        view.path = str((entry / final_name).resolve())   # planned location; not created
        view.build_status = ViewBuildStatus.PLANNED
        view.cache_status = ViewCacheStatus.BYPASSED
        view.warnings.append(CampaignWarning(
            "dry_run", "trajectory preprocessing planned but not executed", Severity.INFO))
        return view

    # ── cache lookup ───────────────────────────────────────────────────────
    src_digests = [fp.digest for fp in src_fps]
    if entry.exists() and not force:
        ok, reason = validate_entry(entry, identity=key, source_digests=src_digests, kind=kind)
        if ok:
            return _reuse(view, entry, final_name, ti_cache, gmx)
        view.cache_status = ViewCacheStatus.INVALID
        view.cache_rejections.append(reason)
        view.warnings.append(CampaignWarning(
            "view_cache_rejected", f"cached view {key} not trusted: {reason}; regenerating",
            Severity.INFO))
    elif force:
        view.cache_status = ViewCacheStatus.BYPASSED
    else:
        view.cache_status = ViewCacheStatus.MISS

    # ── build in an isolated directory ─────────────────────────────────────
    work_dir.mkdir(parents=True, exist_ok=True)
    _sweep_stale_builds(work_dir, key)
    build = Path(tempfile.mkdtemp(prefix=f"{_BUILD_PREFIX}{key}-{os.getpid()}-", dir=work_dir))
    try:
        return _build(view, steps_plan=(requirements, src, topology_path, structure_path,
                                        semantic_index),
                      build=build, entry=entry, work_dir=work_dir, final_name=final_name,
                      gmx=gmx, gmx_ver=gmx_ver, evidence=evidence, groups=groups,
                      src_fps=src_fps, fp_cache=fp_cache, index_path=index_path,
                      ti_cache=ti_cache, src=src, fail=fail, force=force)
    finally:
        if build.exists():
            shutil.rmtree(build, ignore_errors=True)


def _op_record(st, args, res, *, dry=False) -> PreprocessingOperation:
    stdin = st["stdin"]
    return PreprocessingOperation(
        operation=st["op"], tool="gmx trjconv",
        command=(["gmx", *args] if dry else res.argv),
        stdin=stdin, semantic_selection=stdin.split("\n")[0] if stdin else None,
        input_trajectory=args[args.index("-f") + 1],
        input_topology=args[args.index("-s") + 1],
        output=args[args.index("-o") + 1],
        reason=st["reason"], validation_note=st["validation"],
        returncode=None if dry else res.returncode,
        stdout_tail="" if dry else res.tail("stdout", 800),
        stderr_tail="(dry run — not executed)" if dry else res.tail("stderr", 1500),
    )


def _reuse(view: TrajectoryView, entry: Path, final_name: str, ti_cache: Path,
           gmx: str) -> TrajectoryView:
    from analysis.campaign.models import SourceFingerprint as _SF
    m = read_manifest(entry)
    view.path = str((entry / final_name).resolve())
    view.reused = True
    view.build_status = ViewBuildStatus.REUSED
    view.cache_status = ViewCacheStatus.HIT
    view.manifest_path = str((entry / MANIFEST_NAME).resolve())
    view.output_fingerprint = _SF.from_dict(m["output"]["fingerprint"])
    view.operations = [PreprocessingOperation(**{
        k: o.get(k) for k in PreprocessingOperation.__dataclass_fields__ if k in o})
        for o in m.get("operations", [])]
    ref = dict(m.get("time_index") or {}) or None
    if ref:
        # re-establish the authoritative index under the current cache dir:
        # load it by its recorded identity, rescan only if it is gone
        from analysis.campaign.trajectory.time_index import (
            GromacsTimeIndexBackend, index_trajectory, load_cached_index,
        )
        ref["cache_dir"] = str(Path(ti_cache).resolve())
        idx = load_cached_index(view.path, cache_dir=ti_cache, backend_id=ref.get("backend"),
                                backend_version=ref.get("backend_version"))
        if idx is None:
            idx, _ = index_trajectory(view.path, backend=GromacsTimeIndexBackend(gmx),
                                      cache_dir=ti_cache)
            ref.update(_time_index_ref(idx, Path(ti_cache), ref.get("same_timeline_as_source")))
        if (idx.state == ValidationState.INVALID or idx.fingerprint is None
                or idx.fingerprint.digest != ref.get("trajectory_digest")):
            view.warnings.append(CampaignWarning(
                "time_index_unavailable",
                "cached view's time index could not be re-established", Severity.WARN))
    view.time_index_ref = ref
    return view


def _build(view, *, steps_plan, build, entry, work_dir, final_name, gmx, gmx_ver,
           evidence, groups, src_fps, fp_cache, index_path, ti_cache, src, fail, force):
    requirements, src_path, topology_path, structure_path, semantic_index = steps_plan
    steps, _ = _plan_steps(requirements, src_path, topology_path, structure_path,
                           semantic_index, build)

    protected = {str(Path(p).resolve()) for p in [src_path, topology_path, structure_path,
                                                   index_path] if p}
    for st in steps:
        out = str(Path(st["full_args"][st["full_args"].index("-o") + 1]).resolve())
        if out in protected or not out.startswith(str(build.resolve())):
            return fail("unsafe_output_path", f"refusing to write {out}")

    for st in steps:
        args = st["full_args"]
        res = run_gmx(args, stdin=(st["stdin"] + "\n") if st["stdin"] else None,
                      gmx=gmx, timeout=3600)
        view.operations.append(_op_record(st, args, res))
        if not res.ok:
            return fail("preprocessing_failed",
                        f"trajectory preprocessing step '{st['op']}' failed "
                        f"(rc={res.returncode}); analysis must not proceed on raw coordinates")

    last_out = Path(steps[-1]["full_args"][steps[-1]["full_args"].index("-o") + 1])
    if not last_out.is_file() or last_out.stat().st_size == 0:
        return fail("preprocessing_no_output", "expected derived trajectory was not produced")

    # keep only the final product
    removed: list[str] = []
    for st in steps[:-1]:
        inter = Path(st["full_args"][st["full_args"].index("-o") + 1])
        if inter.exists() and inter != last_out:
            inter.unlink()
            removed.append(inter.name)
    product = build / final_name
    if last_out != product:
        last_out.replace(product)
    for extra in build.iterdir():                 # e.g. GROMACS '#backup#' files
        if extra != product:
            if extra.is_dir():
                shutil.rmtree(extra, ignore_errors=True)
            else:
                extra.unlink(missing_ok=True)
                removed.append(extra.name)

    # authoritative per-view timeline (Phase 1) — scanned, never inherited
    from analysis.campaign.trajectory.time_index import GromacsTimeIndexBackend, index_trajectory
    backend = GromacsTimeIndexBackend(gmx)
    tidx, _ = index_trajectory(product, backend=backend, cache_dir=ti_cache)
    if tidx.state == ValidationState.INVALID:
        return fail("time_index_failed",
                    "derived trajectory could not be time-indexed: "
                    + "; ".join(w.message for w in tidx.warnings))
    if tidx.n_frames == 0:
        return fail("output_validation_failed", "derived trajectory contains no frames")
    try:
        sidx, _ = index_trajectory(src, backend=backend, cache_dir=ti_cache)
        same = tidx.same_timeline(sidx) if sidx.state != ValidationState.INVALID else None
    except (ValueError, OSError):
        same = None                               # source not indexable: equality unknown

    out_fp = fingerprint_file(product)
    final_path = (entry / final_name).resolve()
    out_fp.path = str(final_path)
    ti_ref = _time_index_ref(tidx, ti_cache, same)

    manifest = {
        "schema_version": VIEW_MANIFEST_SCHEMA,
        "status": "complete",
        "view_identity": view.cache_key,
        "identity_evidence": evidence,
        "kind": view.kind,
        "requirements": requirements.to_dict(),
        "sources": {
            "trajectories": [fp.to_dict() for fp in src_fps],
            "structures": {Path(k).name: v.to_dict() for k, v in fp_cache.items()},
            "index": fingerprint_file(Path(index_path)).to_dict() if index_path else None,
        },
        "groups": groups,
        "operations": [o.to_dict() for o in view.operations],
        "timeline_args": _timeline_args(steps),
        "gromacs_version": gmx_ver,
        "output": {"relative_path": final_name, "fingerprint": out_fp.to_dict()},
        "intermediates_removed": removed,
        "time_index": ti_ref,
        # scientific authorisation — not part of the identity
        "policy": {"planned": view.policy_planned, "decisions": view.decisions,
                   "context": view.policy_ref},
        # informational only — never part of the identity
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
    }
    tmp = build / (MANIFEST_NAME + ".tmp")
    tmp.write_text(json.dumps(manifest, indent=2) + "\n")
    os.replace(tmp, build / MANIFEST_NAME)        # manifest last

    # atomic promotion
    if entry.exists():
        ok, _ = validate_entry(entry, identity=view.cache_key,
                               source_digests=[fp.digest for fp in src_fps], kind=view.kind)
        if ok and not force:                      # a concurrent build finished first
            return _reuse(view, entry, final_name, ti_cache, gmx)
        _remove_entry(work_dir, entry)
    os.replace(build, entry)

    view.path = str(final_path)
    view.build_status = ViewBuildStatus.BUILT
    view.manifest_path = str((entry / MANIFEST_NAME).resolve())
    view.output_fingerprint = out_fp
    view.time_index_ref = ti_ref
    return view
