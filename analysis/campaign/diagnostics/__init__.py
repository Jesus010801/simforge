"""Read-only trajectory diagnostics (Trajectory Review Phase 4).

Observes and quantifies timeline, box and component-motion properties of a
trajectory; lists *candidate* remediations; never transforms coordinates,
never builds views, never runs corrective ``trjconv`` commands.
"""
from __future__ import annotations

from analysis.campaign.diagnostics.context import (
    ComponentMotionSeries, DiagnosticContext, ProbeError, groups_from_semantic_index,
)
from analysis.campaign.diagnostics.params import DiagnosticParams
from analysis.campaign.diagnostics.registry import (
    REPORT_FILENAME, DetectorSpec, DiagnosticsStatus, all_detectors, detector,
    detector_identity, diagnostics_status, report_identity, run_diagnostics, write_report,
)

__all__ = [
    "ComponentMotionSeries", "DiagnosticContext", "DiagnosticParams", "DetectorSpec",
    "DiagnosticsStatus", "diagnostics_status", "report_identity",
    "ProbeError", "REPORT_FILENAME", "all_detectors", "detector", "detector_identity",
    "diagnose_system", "groups_from_semantic_index", "run_diagnostics", "write_report",
]


def diagnose_system(rec, *, gmx: str = "gmx", cache_dir=None, params=None, view_ref=None):
    """Diagnose a campaign :class:`SystemRecord`'s single production trajectory.

    Uses the record's semantic index (resolved components only) and its .tpr
    (or structure) for probes.  Returns a :class:`DiagnosticReport`.
    """
    paths = rec.trajectory_paths
    if len(paths) != 1:
        raise ValueError("diagnose_system needs exactly one production trajectory "
                         f"(got {len(paths)}); segments are diagnosed individually")
    idx = rec.semantic_index
    groups = {}
    if idx is not None and idx.path:
        groups = groups_from_semantic_index(idx.path, rec.components)
    # -s for the read-only probes: a .tpr (real masses) when available, otherwise
    # the reference structure (COG) — a .top is not a GROMACS run-input/structure
    tpr = rec.topology_path if rec.topology_path and rec.topology_path.lower().endswith(".tpr") else None
    ctx = DiagnosticContext(
        trajectory_path=paths[0],
        structure_path=tpr or rec.structure_path,
        index_path=idx.path if idx is not None else None, groups=groups,
        view_ref=view_ref, membrane_present=rec.membrane_present,
        params=params or DiagnosticParams(), gmx=gmx,
        cache_dir=str(cache_dir) if cache_dir else None)
    return run_diagnostics(ctx)
