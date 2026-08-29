"""Top-level orchestration for the study-campaign layer.

``run_inspect``  – discovery -> manifest -> component inference -> trajectory
                   inspection -> validation.  Writes the manifest + validation
                   report.  Runs **no** scientific observables.

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
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest = _prepare_manifest(
        study_root, manifest_path=manifest_path, gmx=gmx,
        strong_fingerprint=strong_fingerprint, inspect_trajectories=inspect_trajectories,
    )
    paths = write_manifest(manifest, out_dir)

    validation = validate_study(manifest, [])
    (out_dir / "validation_report.json").write_text(
        json.dumps(validation.to_dict(), indent=2) + "\n")
    (out_dir / "validation_report.txt").write_text("\n".join(validation.lines) + "\n")

    result = CampaignRunResult(
        study_root=str(study_root), output_dir=str(out_dir),
        manifest=manifest, validation=validation, requested_analyses=[],
        output_files=[paths["yaml"], paths["json"],
                      str(out_dir / "validation_report.json"),
                      str(out_dir / "validation_report.txt")],
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
) -> CampaignRunResult:
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

            req = spec.trajectory_requirements(parameters)
            view = _resolve_view(rec, req, sys_out, gmx, force, result, dry_run=dry_run)
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
    )
    rec.semantic_index = idx
    for w in idx.warnings:
        result.warnings.append(CampaignWarning(
            "semantic_index", w, Severity.INFO, scope=rec.system_id))


_VIEW_CACHE: dict[tuple[str, str], object] = {}


def _resolve_view(rec, req, sys_out: Path, gmx: str, force: bool,
                  result: CampaignRunResult, *, dry_run: bool = False):
    key = (rec.system_id, req.cache_token(), dry_run)
    if key in _VIEW_CACHE and not force:
        return _VIEW_CACHE[key]
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
    )
    _VIEW_CACHE[key] = view
    for w in view.warnings:
        result.warnings.append(w)
    return view


def clear_view_cache() -> None:
    _VIEW_CACHE.clear()
