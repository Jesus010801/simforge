"""Top-level orchestration for the study-campaign layer.

``run_inspect``  – discovery -> manifest -> component inference -> trajectory
                   inspection -> validation.  Exports reports only when an
                   output directory is explicitly supplied.  Runs **no** scientific observables.

``run_analyze``  – reuses the exact same manifest/validation pipeline, then
                   executes **only** the explicitly requested analyses,
                   per system, with per-analysis output directories and full
                   provenance.

Nothing is analysed by default.  ``requested_analyses == []`` performs discovery
+ validation and returns (matching ``study inspect``).
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import yaml

from analysis.campaign.models import (
    AnalysisResult, AnalysisStatus, CampaignRunResult, CampaignWarning,
    ClassificationState, Severity, StudyManifest, SystemRecord, ValidationState,
)
from analysis.campaign.manifest import build_manifest, load_manifest, write_manifest
from analysis.campaign.observables import registry
from analysis.campaign.observables.base import AnalysisContext
from analysis.campaign.provenance import build_observable_provenance
from analysis.campaign.structure.index_builder import build_semantic_index
from analysis.campaign.trajectory.preprocessor import build_view
from analysis.campaign.validator import validate_study

DEFAULT_OUTPUT_DIRNAME = "simforge_analysis"


def _prepare_manifest(
    study_root: Path, *, manifest_path: Optional[str], gmx: str,
    strong_fingerprint: bool, inspect_trajectories: bool,
) -> StudyManifest:
    if manifest_path:
        return load_manifest(manifest_path)
    return build_manifest(
        study_root, inspect_trajectories=inspect_trajectories,
        strong_fingerprint=strong_fingerprint, gmx=gmx,
    )


def run_inspect(
    study_root: str | Path,
    *,
    output_dir: Optional[str | Path] = None,
    manifest_path: Optional[str] = None,
    gmx: str = "gmx",
    strong_fingerprint: bool = False,
    inspect_trajectories: bool = True,
) -> CampaignRunResult:
    study_root = Path(study_root).resolve()
    out_dir = Path(output_dir) if output_dir else study_root / DEFAULT_OUTPUT_DIRNAME
    manifest = _prepare_manifest(
        study_root, manifest_path=manifest_path, gmx=gmx,
        strong_fingerprint=strong_fingerprint, inspect_trajectories=inspect_trajectories,
    )
    validation = validate_study(manifest, [])
    output_files = []
    if output_dir is not None:
        destinations = [out_dir / n for n in (
            'study_manifest.yaml', 'study_manifest.json',
            'validation_report.json', 'validation_report.txt')]
        if any(p.exists() or p.is_symlink() for p in destinations):
            raise FileExistsError(f'inspection output already exists in {out_dir}')
        out_dir.mkdir(parents=True, exist_ok=True)
        paths = write_manifest(manifest, out_dir)
        destinations[2].write_text(json.dumps(validation.to_dict(), indent=2) + '\n')
        destinations[3].write_text('\n'.join(validation.lines) + '\n')
        output_files = [paths['yaml'], paths['json'], str(destinations[2]), str(destinations[3])]

    result = CampaignRunResult(
        study_root=str(study_root), output_dir=str(out_dir),
        manifest=manifest, validation=validation, requested_analyses=[],
        output_files=output_files,
        warnings=list(manifest.warnings),
    )
    return result


def run_analyze(
    study_root: str | Path,
    requested_analyses: list[str],
    *,
    output_dir: Optional[str | Path] = None,
    manifest_path: Optional[str] = None,
    gmx: str = "gmx",
    strong_fingerprint: bool = False,
    inspect_trajectories: bool = True,
    dry_run: bool = False,
    force: bool = False,
    parameters: Optional[dict] = None,
    policy_intent=None,
    diagnostics="auto",
    diagnostics_cache_dir: Optional[str | Path] = None,
    reuse: bool = True,
    reuse_roots: Optional[list] = None,
    use_imported: Optional[dict] = None,
) -> CampaignRunResult:
    """Discover/validate, then run the requested observables per system.

    ``diagnostics``: ``"auto"`` (default) runs the read-only Phase 4 detectors
    once per system and hands the report to the preprocessing policy;
    ``"off"`` supplies none (policy stays conservative); a ``DiagnosticReport``
    (or ``{system_id: report}``) is used after validating it against the
    system's trajectory and semantic groups.

    ``reuse`` (Phase 8): before computing, stored native results under this
    output root (and ``reuse_roots``) are checked for scientific compatibility;
    a fully compatible one is returned as CACHED without running GROMACS.
    ``use_imported``: ``{analysis_id: ResultArray}`` explicitly chosen
    external series — labelled externally supplied, never verified equivalent.
    """
    registry.ensure_loaded()
    study_root = Path(study_root).resolve()
    out_dir = Path(output_dir) if output_dir else study_root / DEFAULT_OUTPUT_DIRNAME
    out_dir.mkdir(parents=True, exist_ok=True)
    parameters = parameters or {}

    manifest = _prepare_manifest(
        study_root, manifest_path=manifest_path, gmx=gmx,
        strong_fingerprint=strong_fingerprint, inspect_trajectories=inspect_trajectories,
    )
    paths = write_manifest(manifest, out_dir)
    validation = validate_study(manifest, requested_analyses)
    (out_dir / "validation_report.json").write_text(
        json.dumps(validation.to_dict(), indent=2) + "\n")
    (out_dir / "validation_report.txt").write_text("\n".join(validation.lines) + "\n")

    result = CampaignRunResult(
        study_root=str(study_root), output_dir=str(out_dir),
        manifest=manifest, validation=validation,
        requested_analyses=list(requested_analyses),
        output_files=[paths["yaml"], paths["json"]],
        warnings=list(manifest.warnings),
    )

    unknown = [a for a in requested_analyses if not registry.is_registered(a)]
    for a in unknown:
        result.warnings.append(CampaignWarning(
            "unknown_analysis",
            f"'{a}' is not a registered analysis. Run `simforge study analyses`. "
            f"Valid: {registry.ids()}",
            Severity.ERROR))
    valid_analyses = [a for a in requested_analyses if registry.is_registered(a)]
    if not valid_analyses:
        return result

    status_lookup = {
        (o.system_id, o.analysis_id): o for o in validation.observable_statuses
    }

    systems_dir = out_dir / "systems"
    for rec in manifest.systems:
        sys_out = systems_dir / rec.system_id
        _ensure_semantic_index(rec, sys_out, gmx, result)
        diag_report, diag_ref = _system_diagnostics(
            rec, sys_out, gmx, result, diagnostics=diagnostics, dry_run=dry_run,
            cache_dir=diagnostics_cache_dir)

        for analysis_id in valid_analyses:
            spec = registry.get(analysis_id)
            obs_out = sys_out / "observables" / analysis_id
            st = status_lookup.get((rec.system_id, analysis_id))

            if st and st.state == ValidationState.INVALID:
                result.results.append(AnalysisResult(
                    analysis_id=analysis_id, system_id=rec.system_id,
                    status=AnalysisStatus.SKIPPED, message=st.reason))
                continue
            if st and st.state == ClassificationState.REVIEW_REQUIRED:
                result.results.append(AnalysisResult(
                    analysis_id=analysis_id, system_id=rec.system_id,
                    status=AnalysisStatus.REVIEW_REQUIRED,
                    message=st.reason + " | " + st.remediation))
                continue

            if use_imported and analysis_id in use_imported:
                result.results.append(_explicit_import(spec, rec, use_imported[analysis_id]))
                continue
            tic = _time_index_cache(sys_out, diag_ref, diagnostics_cache_dir)
            evaluated: list = []
            if reuse and not dry_run and not force:
                cached, evaluated = _try_reuse(
                    spec, rec, sys_out, out_dir, obs_out, gmx, result, parameters,
                    policy_intent, diag_report, diag_ref, tic, reuse_roots or [])
                if cached is not None:
                    result.results.append(cached)
                    result.output_files.extend(cached.output_files)
                    continue

            req = spec.trajectory_requirements(parameters)
            view = _resolve_view(rec, req, sys_out, gmx, force, result, dry_run=dry_run,
                                 purpose=spec.purpose, intent=policy_intent,
                                 diagnostics=diag_report, diagnostics_ref=diag_ref)
            ctx = AnalysisContext(
                system=rec,
                semantic_index=rec.semantic_index,
                trajectory_view=view,
                topology_path=rec.topology_path or (rec.structure_path or ""),
                structure_path=rec.structure_path,
                output_dir=obs_out,
                parameters=parameters,
                gmx=gmx,
                dry_run=dry_run,
                time_index_cache_dir=tic,
            )
            try:
                ares = spec.execute(ctx)
            except Exception as exc:  # noqa: BLE001 - observable must not crash the run
                ares = AnalysisResult(
                    analysis_id=analysis_id, system_id=rec.system_id,
                    status=AnalysisStatus.FAILED,
                    message=f"observable raised {type(exc).__name__}: {exc}")
                ares.warnings.append(CampaignWarning(
                    "observable_exception", ares.message, Severity.ERROR, scope=rec.system_id))

            ares.trajectory_view_kind = view.kind
            if evaluated:
                ares.compatibility = {"decision": "recomputed",
                                      "evaluated": [e.to_dict() for e in evaluated]}
            if ares.status in (AnalysisStatus.SUCCESS, AnalysisStatus.PLANNED,
                               AnalysisStatus.FAILED):
                prov = build_observable_provenance(
                    system=rec, result=ares,
                    semantic_index=rec.semantic_index, trajectory_view=view)
                obs_out.mkdir(parents=True, exist_ok=True)
                pj = obs_out / "provenance.json"
                pj.write_text(json.dumps(prov.to_dict(), indent=2) + "\n")
                ares.provenance_path = str(pj.resolve())
                ares.output_files.append(str(pj.resolve()))
            result.results.append(ares)
            result.output_files.extend(ares.output_files)

    # results summary file
    (out_dir / "analysis_results.json").write_text(
        json.dumps([r.to_dict() for r in result.results], indent=2) + "\n")
    result.output_files.append(str(out_dir / "analysis_results.json"))
    return result


def _ensure_semantic_index(
    rec: SystemRecord, sys_out: Path, gmx: str, result: CampaignRunResult,
) -> None:
    if rec.semantic_index is not None or not rec.structure_path:
        return
    ndx = sys_out / "indices" / "semantic_index.ndx"
    resolved = [c for c in rec.components
                if c.classification_state == ClassificationState.RESOLVED]
    if not resolved:
        return
    idx = build_semantic_index(
        structure_path=rec.structure_path, components=rec.components,
        out_ndx=ndx, gmx=gmx, existing_user_index=rec.index_path,
        annotations=rec.annotations,
    )
    rec.semantic_index = idx
    for w in idx.warnings:
        result.warnings.append(CampaignWarning(
            "semantic_index", w, Severity.INFO, scope=rec.system_id))


_VIEW_CACHE: dict[tuple, object] = {}


def _resolve_view(rec, req, sys_out: Path, gmx: str, force: bool,
                  result: CampaignRunResult, *, dry_run: bool = False,
                  purpose: Optional[str] = None, intent=None, diagnostics=None,
                  diagnostics_ref: Optional[dict] = None):
    # the semantic index content is part of the key: same group names with
    # different atoms must never share an in-process view
    idx = rec.semantic_index
    idx_digest = ""
    if idx is not None and idx.path and Path(idx.path).is_file():
        from analysis.campaign.fingerprint import fingerprint_file
        idx_digest = fingerprint_file(Path(idx.path)).digest
    from analysis.campaign.trajectory.policy import (
        PolicyContext, PolicyIntent, plan_preprocessing,
    )
    intent = intent or PolicyIntent()
    diag_key = (diagnostics_ref or {}).get("report_identity")
    key = (rec.system_id, req.cache_token(), dry_run, idx_digest, purpose,
           intent.source, intent.operations, diag_key)
    if key in _VIEW_CACHE and not force:
        return _VIEW_CACHE[key]
    # scientific admissibility first; build_view only executes admitted requests
    plan = plan_preprocessing(req, purpose=purpose,
                              context=PolicyContext.from_system(rec, diagnostics=diagnostics),
                              intent=intent)
    policy_ref = {"purpose": plan.purpose, "intent": intent.to_dict(),
                  "diagnostics": dict(diagnostics_ref or {"execution": "not_supplied",
                                                          "status": "unavailable"})}
    if not plan.executable:
        from analysis.campaign.models import TrajectoryView, ViewBuildStatus
        view = TrajectoryView(kind=req.view_kind(), requirements=req, safe=False,
                              build_status=ViewBuildStatus.FAILED, policy_planned=True,
                              decisions=[d.to_dict() for d in plan.decisions],
                              policy_ref=policy_ref)
        view.warnings.append(CampaignWarning(
            "policy_blocked",
            "preprocessing not admitted: " + "; ".join(
                f"{d.operation} = {d.classification} ({d.rule_id}): {d.reason}"
                for d in plan.decisions if not d.applied),
            Severity.REVIEW, scope=rec.system_id))
        _VIEW_CACHE[key] = view
        result.warnings.extend(view.warnings)
        return view
    topo_fp = rec.source_fingerprints.get("topology")
    view = build_view(
        requirements=req,
        trajectory_paths=rec.trajectory_paths,
        topology_path=rec.topology_path or (rec.structure_path or ""),
        structure_path=rec.structure_path,
        semantic_index=rec.semantic_index,
        source_fingerprints=[fp for k, fp in rec.source_fingerprints.items()
                             if k.startswith("trajectory")],
        topology_fingerprint=topo_fp,
        work_dir=sys_out / "trajectories",
        gmx=gmx, force=force, dry_run=dry_run,
        decisions=plan.decisions, policy_ref=policy_ref,
    )
    _VIEW_CACHE[key] = view
    for w in view.warnings:
        result.warnings.append(w)
    return view


def clear_view_cache() -> None:
    _VIEW_CACHE.clear()


# ═══════════════════════════════════════════════════════════════════════════════
# Diagnostics — once per system, handed to the (pure) policy
# ═══════════════════════════════════════════════════════════════════════════════

DIAGNOSTICS_DIRNAME = "diagnostics"


# ═══════════════════════════════════════════════════════════════════════════════
# Phase 8 — reuse of scientifically compatible results
# ═══════════════════════════════════════════════════════════════════════════════

def _requested_timeline(view, tic: Path, gmx: str) -> Optional[list[float]]:
    """The planned view's frame times from an existing Phase 1 cache (no scan)."""
    from analysis.campaign.gmx import gmx_version
    from analysis.campaign.trajectory.time_index import GromacsTimeIndexBackend, load_cached_index
    if not view.path or not Path(view.path).is_file():
        return None
    backend = GromacsTimeIndexBackend(gmx)
    for cache in (tic, Path(view.path).parent.parent / "time_index"):
        idx = load_cached_index(view.path, cache_dir=cache, backend_id=backend.id,
                                backend_version=gmx_version(gmx))
        if idx is not None:
            return list(idx.times_ps)
    return None


def _try_reuse(spec, rec, sys_out, out_dir, obs_out, gmx, result, parameters, intent,
               diag_report, diag_ref, tic, extra_roots):
    """(CACHED result | None, evaluations).  Order: current applicability ->
    candidates -> policy-planned view identity (no build) -> definition
    evidence -> per-dimension compatibility."""
    from analysis.campaign.compatibility import (
        definition_identity, discover_candidates, resolve_compatible_result,
    )
    if not spec.applicability(rec, parameters, rec.semantic_index)["applicable"]:
        return None, []                                  # current state decides first
    candidates = []
    for root in [out_dir, *extra_roots]:
        candidates += discover_candidates(root, rec.system_id, spec.id, current_dir=obs_out)
    if not candidates:
        return None, []
    planned = _resolve_view(rec, spec.trajectory_requirements(parameters), sys_out, gmx, False,
                            result, dry_run=True, purpose=spec.purpose, intent=intent,
                            diagnostics=diag_report, diagnostics_ref=diag_ref)
    if not planned.safe or not planned.cache_key:
        return None, []                                  # the current policy blocks it
    ctx = AnalysisContext(system=rec, semantic_index=rec.semantic_index, trajectory_view=planned,
                          topology_path=rec.topology_path or (rec.structure_path or ""),
                          structure_path=rec.structure_path, output_dir=obs_out,
                          parameters=parameters, gmx=gmx, dry_run=True, time_index_cache_dir=tic)
    evidence = spec.definition_evidence(ctx)
    if evidence is None:
        return None, []
    from analysis.campaign.observables.generic import _params
    chosen, reports = resolve_compatible_result(
        candidates, requested_evidence=evidence, requested_view_ref=planned.cache_key,
        requested_schema=spec.output_schema(_params(parameters, spec.id)),
        requested_timeline=_requested_timeline(planned, tic, gmx),
        view_invariant=spec.view_invariant)
    if chosen is None:
        return None, reports
    orig = chosen.candidate.result
    cached = AnalysisResult.from_dict(orig.to_dict())     # a new record; the stored one is untouched
    cached.status = AnalysisStatus.CACHED
    cached.cached = True
    cached.message = f"reused compatible result ({chosen.candidate.source})"
    cached.provenance_path = chosen.candidate.source
    cached.reused_from = {"provenance": chosen.candidate.source,
                          "definition_token": orig.definition_token,
                          "definition_identity": definition_identity(orig.definition_evidence),
                          "original_status": orig.status}
    cached.compatibility = {"decision": "reused", "chosen": chosen.to_dict(),
                            "evaluated": [r.to_dict() for r in reports]}
    obs_out.mkdir(parents=True, exist_ok=True)
    record = obs_out / "reuse.json"                        # never overwrites provenance.json
    record.write_text(json.dumps({"resolver": chosen.resolver, "reused": cached.reused_from,
                                  "compatibility": cached.compatibility,
                                  "requested_definition_identity": definition_identity(evidence)},
                                 indent=2) + "\n")
    cached.output_files = [a.storage.path for a in cached.arrays if a.storage] + [str(record.resolve())]
    return cached, reports


def _explicit_import(spec, rec, array) -> AnalysisResult:
    """A user-chosen external series: labelled, integrity-checked, never 'equivalent'."""
    from analysis.campaign.compatibility import Candidate, ProvenanceTier, evaluate, storage_integrity
    res = AnalysisResult(analysis_id=spec.id, system_id=rec.system_id,
                         status=AnalysisStatus.CACHED, cached=True, arrays=[array])
    ok, why = storage_integrity(res)
    if not ok:
        res.status = AnalysisStatus.REVIEW_REQUIRED
        res.arrays = []
        res.message = f"explicitly imported series not usable: {why}"
        return res
    report = evaluate(Candidate(source=array.storage.path, tier=ProvenanceTier.DECLARED, result=res),
                      requested_evidence=None, requested_view_ref=None,
                      requested_schema=spec.output_schema({}))
    res.message = ("EXPLICIT IMPORT: externally supplied series chosen by the user — "
                   "definition not independently verified")
    res.reused_from = {"external_file": array.storage.path,
                       "fingerprint": array.storage.fingerprint.digest, "user_selected": True}
    res.compatibility = {"decision": "explicit_import", "chosen": report.to_dict()}
    return res


def _time_index_cache(sys_out: Path, diag_ref: Optional[dict], cache_dir) -> Path:
    """The Phase 1 time-index cache the diagnostics pass used (so observables
    align against an already-indexed timeline instead of rescanning)."""
    probe = (diag_ref or {}).get("probe_cache")
    base = Path(probe) if probe else (Path(cache_dir) if cache_dir
                                      else sys_out / DIAGNOSTICS_DIRNAME / "cache")
    return base / "time_index"


def _system_diagnostics(rec, sys_out: Path, gmx: str, result: CampaignRunResult, *,
                        diagnostics, dry_run: bool, cache_dir) -> tuple:
    """Return ``(report | None, reference dict)``.  Read-only; never builds views.

    Precedence: a supplied report (validated) > an automatically generated one.
    """
    from analysis.campaign.diagnostics import (
        diagnose_system, diagnostics_status, report_identity, write_report,
    )
    from analysis.campaign.gmx import gmx_available

    def unavailable(execution, reason):
        return None, {"execution": execution, "status": "unavailable", "reason": reason}

    supplied = diagnostics
    if isinstance(diagnostics, dict) and not hasattr(diagnostics, "detectors"):
        supplied = diagnostics.get(rec.system_id)
        if supplied is None:
            return unavailable("not_supplied", "no report supplied for this system")
    if supplied is not None and not isinstance(supplied, str):
        problem = _validate_report(supplied, rec)
        if problem is None:
            return supplied, {"execution": "supplied", "status": diagnostics_status(supplied),
                              "report_identity": report_identity(supplied)}
        result.warnings.append(CampaignWarning(
            "diagnostics_report_rejected",
            f"supplied diagnostics report not used: {problem}; generating a fresh one",
            Severity.WARN, scope=rec.system_id))
        diagnostics = "auto"
    if diagnostics == "off":
        return unavailable("off", "diagnostics disabled by caller")
    if dry_run:
        return unavailable("skipped", "dry run: no probes executed")
    if not gmx_available(gmx):
        return unavailable("skipped", f"'{gmx}' not available")
    if len(rec.trajectory_paths) != 1:
        return unavailable("skipped", f"{len(rec.trajectory_paths)} production trajectories; "
                                      f"segments are not diagnosed as one timeline")
    cache = Path(cache_dir) if cache_dir else sys_out / DIAGNOSTICS_DIRNAME / "cache"
    try:
        report = diagnose_system(rec, gmx=gmx, cache_dir=cache)
    except Exception as exc:  # noqa: BLE001 — never silently drop: record it
        result.warnings.append(CampaignWarning(
            "diagnostics_failed", f"diagnostics could not run: {exc}", Severity.WARN,
            scope=rec.system_id))
        return None, {"execution": "auto", "status": "failed", "reason": str(exc)}
    path = write_report(report, sys_out / DIAGNOSTICS_DIRNAME)
    result.output_files.append(str(path))
    return report, {"execution": "auto", "status": diagnostics_status(report),
                    "report_identity": report_identity(report),
                    "report_path": str(path.resolve()),
                    "probe_cache": str(cache.resolve()),
                    "probes_from_cache": all(p.get("from_cache") for p in report.probes)
                    if report.probes else None}


def _validate_report(report, rec) -> Optional[str]:
    """A supplied report must describe this system's trajectory content and the
    same atoms for every group it used."""
    from analysis.campaign.fingerprint import fingerprint_file
    from analysis.campaign.results import index_group_evidence
    if len(rec.trajectory_paths) != 1:
        return "system has no single production trajectory"
    traj = Path(rec.trajectory_paths[0])
    if not traj.is_file():
        return "production trajectory missing"
    fp = report.trajectory_fingerprint
    if fp is None or fp.digest != fingerprint_file(traj).digest:
        return "trajectory fingerprint does not match"
    idx = rec.semantic_index
    for role, g in (report.groups or {}).items():
        name, ev = (g or {}).get("group"), (g or {}).get("evidence")
        if not name or ev is None:
            continue
        cur = index_group_evidence(idx.path, name) if idx is not None and idx.path else None
        if cur != ev:
            return f"group '{name}' ({role}) atoms differ from the current semantic index"
    return None
