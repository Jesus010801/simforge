"""Detector registry and the (read-only) diagnostic runner.

A detector is a plain function registered with :func:`detector`; adding one
never touches :func:`run_diagnostics`.  Levels:

* 0 — timeline / metadata          * 1 — box behaviour
* 2 — component centre motion      * 3 — sampled coordinate geometry
* 4 — domain-specific (membrane, channel, ...)

Reserved-but-unimplemented detectors are registered with ``deferred=...`` and
reported as such, never silently omitted.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

from analysis.campaign.diagnostics.context import DiagnosticContext, ProbeError
from analysis.campaign.models import (
    DetectorRun, DetectorRunStatus, Diagnostic, DiagnosticReport,
)

REPORT_FILENAME = "trajectory_diagnostics.json"


@dataclass(frozen=True)
class DetectorSpec:
    id: str
    version: str
    level: int
    description: str
    detect: Optional[Callable[[DiagnosticContext], list[Diagnostic]]] = None
    #: (ctx) -> (applicable, reason)
    applicable: Callable[[DiagnosticContext], tuple[bool, str]] = lambda ctx: (True, "")
    #: (ctx) -> roles whose group evidence enters this detector's identity
    roles: Callable[[DiagnosticContext], list[str]] = lambda ctx: []
    sampling: str = ""
    deferred: str = ""                       # non-empty => reserved, not run


_DETECTORS: dict[str, DetectorSpec] = {}


def detector(id: str, version: str, level: int, description: str, *,
             applicable=None, roles=None, sampling: str = ""):
    def wrap(fn):
        register(DetectorSpec(id=id, version=version, level=level, description=description,
                              detect=fn, applicable=applicable or (lambda ctx: (True, "")),
                              roles=roles or (lambda ctx: []), sampling=sampling))
        return fn
    return wrap


def register(spec: DetectorSpec) -> None:
    if spec.id in _DETECTORS:
        raise ValueError(f"detector already registered: {spec.id!r}")
    _DETECTORS[spec.id] = spec


def all_detectors() -> list[DetectorSpec]:
    _ensure_loaded()
    return sorted(_DETECTORS.values(), key=lambda s: (s.level, s.id))


def get(detector_id: str) -> DetectorSpec:
    _ensure_loaded()
    return _DETECTORS[detector_id]


def _ensure_loaded() -> None:
    # importing the modules registers the built-in detectors (idempotent)
    from analysis.campaign.diagnostics import box, deferred, motion, timeline  # noqa: F401


def detector_identity(spec: DetectorSpec, ctx: DiagnosticContext) -> str:
    """Deterministic identity of one detector run: detector id/version, params,
    trajectory content, view and the atom-set evidence of the groups it uses
    (unused groups never enter)."""
    from analysis.campaign.results import definition_token
    fp = ctx.trajectory_fingerprint()
    return definition_token({
        "detector": spec.id, "version": spec.version,
        "params": ctx.params.to_dict(),
        "trajectory": fp.digest if fp else None,
        "view_ref": ctx.view_ref,
        "groups": {r: ctx.group_evidence(r) for r in sorted(spec.roles(ctx))},
    })


def run_diagnostics(ctx: DiagnosticContext, *,
                    only: Optional[list[str]] = None) -> DiagnosticReport:
    """Run every applicable registered detector.  Read-only."""
    fp = ctx.trajectory_fingerprint()
    report = DiagnosticReport(trajectory_path=str(Path(ctx.trajectory_path).resolve()),
                              trajectory_fingerprint=fp, view_ref=ctx.view_ref,
                              params=ctx.params.to_dict())
    for spec in all_detectors():
        if only and spec.id not in only:
            continue
        run = DetectorRun(detector_id=spec.id, version=spec.version, level=spec.level,
                          status=DetectorRunStatus.RAN)
        if spec.deferred:
            run.status, run.reason = DetectorRunStatus.DEFERRED, spec.deferred
            report.detectors.append(run)
            continue
        try:
            ok, reason = spec.applicable(ctx)
        except ProbeError as exc:
            ok, reason = False, str(exc)
        if not ok:
            run.status, run.reason = DetectorRunStatus.NOT_APPLICABLE, reason
            report.detectors.append(run)
            continue
        run.identity = detector_identity(spec, ctx)
        try:
            found = spec.detect(ctx)
        except ProbeError as exc:
            run.status, run.reason = DetectorRunStatus.FAILED, f"probe failed: {exc}"
            report.detectors.append(run)
            continue
        for d in found:
            d.detector = {"id": spec.id, "version": spec.version}
            d.provenance.setdefault("detector_identity", run.identity)
            d.provenance.setdefault("trajectory_digest", fp.digest if fp else None)
            d.provenance.setdefault("view_ref", ctx.view_ref)
            d.provenance.setdefault("params", {k: v for k, v in ctx.params.to_dict().items()})
        run.n_diagnostics = len(found)
        run.sampling = found[0].sampling if found else _default_sampling(spec, ctx)
        report.detectors.append(run)
        report.diagnostics.extend(found)

    ti = ctx.time_index
    report.n_frames = ti.n_frames if ti is not None else None
    for role in sorted({r for s in all_detectors() if not s.deferred for r in _safe_roles(s, ctx)}):
        report.groups[role] = {"group": ctx.groups.get(role),
                               "evidence": ctx.group_evidence(role)}
    report.probes = list(ctx.probe_records)
    return report


def _safe_roles(spec: DetectorSpec, ctx: DiagnosticContext) -> list[str]:
    try:
        return spec.roles(ctx)
    except Exception:  # noqa: BLE001
        return []


def _default_sampling(spec: DetectorSpec, ctx: DiagnosticContext) -> dict:
    n = ctx.time_index.n_frames if ctx.time_index is not None else None
    return {"strategy": spec.sampling, "frames_examined": n, "total_frames": n}


def write_report(report: DiagnosticReport, out_dir: str | Path) -> Path:
    out = Path(out_dir) / REPORT_FILENAME
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(report.to_dict(), indent=2) + "\n")
    tmp.replace(out)
    return out


class DiagnosticsStatus:
    """Overall state of a diagnostics pass (DetectorRun statuses stay per detector)."""
    COMPLETE = "complete"        # no detector failed (not_applicable / deferred are normal)
    PARTIAL = "partial"          # some detectors ran, at least one failed
    FAILED = "failed"            # no detector ran
    UNAVAILABLE = "unavailable"  # diagnostics were not run at all


def diagnostics_status(report) -> str:
    if report is None:
        return DiagnosticsStatus.UNAVAILABLE
    active = [r for r in report.detectors if r.status != "deferred"]
    ran = [r for r in active if r.status == "ran"]
    failed = [r for r in active if r.status == "failed"]
    if not ran:
        return DiagnosticsStatus.FAILED
    return DiagnosticsStatus.PARTIAL if failed else DiagnosticsStatus.COMPLETE


def report_identity(report) -> str:
    """Deterministic identity of a report's evidence (not its file path or time)."""
    from analysis.campaign.results import definition_token
    fp = report.trajectory_fingerprint
    return definition_token({
        "schema": report.schema_version,
        "trajectory": fp.digest if fp else None,
        "view_ref": report.view_ref,
        "detectors": sorted((r.detector_id, r.version, r.status, r.identity or "")
                            for r in report.detectors),
        "diagnostics": sorted((d.code, d.scope, (d.provenance or {}).get("detector_identity") or "")
                              for d in report.diagnostics),
    })
