"""Provenance for one observable result.

Mirrors ``analysis/fel/provenance.py`` but also captures the GROMACS toolchain
(via ``core.provenance.collect_environment``) and the full semantic /
trajectory-view chain, so any result is reproducible from the manifest alone.
"""
from __future__ import annotations

import datetime
import sys
from dataclasses import dataclass, field
from typing import Any, Optional

from analysis.campaign.models import (
    AnalysisResult, SemanticIndex, SourceFingerprint, SystemRecord, TrajectoryView,
)


def _environment() -> dict:
    try:
        from core.provenance import collect_environment
        return collect_environment()
    except Exception:  # noqa: BLE001
        return {"note": "core.provenance.collect_environment unavailable"}


@dataclass
class ObservableProvenance:
    schema: str
    timestamp: str
    simforge_version: str
    argv: list[str]
    environment: dict
    system_id: str
    condition_id: str
    replicate_id: str
    analysis_id: str
    parameters: dict
    fit_selection: Optional[str]
    measure_selection: Optional[str]
    reference_frame: str
    source_files: dict
    source_fingerprints: dict
    semantic_index: Optional[dict]
    trajectory_view: Optional[dict]
    component_evidence: list[dict]
    trajectory_inspection: Optional[dict]
    result: dict
    warnings: list[dict] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "schema": self.schema,
            "timestamp": self.timestamp,
            "simforge_version": self.simforge_version,
            "argv": self.argv,
            "environment": self.environment,
            "system_id": self.system_id,
            "condition_id": self.condition_id,
            "replicate_id": self.replicate_id,
            "analysis_id": self.analysis_id,
            "parameters": self.parameters,
            "fit_selection": self.fit_selection,
            "measure_selection": self.measure_selection,
            "reference_frame": self.reference_frame,
            "source_files": self.source_files,
            "source_fingerprints": self.source_fingerprints,
            "semantic_index": self.semantic_index,
            "trajectory_view": self.trajectory_view,
            "component_evidence": self.component_evidence,
            "trajectory_inspection": self.trajectory_inspection,
            "result": self.result,
            "warnings": self.warnings,
        }


def build_observable_provenance(
    *,
    system: SystemRecord,
    result: AnalysisResult,
    semantic_index: Optional[SemanticIndex],
    trajectory_view: Optional[TrajectoryView],
) -> ObservableProvenance:
    try:
        from simforge import __version__ as version
    except Exception:  # noqa: BLE001
        version = "unknown"

    comp_ev: list[dict] = []
    for c in system.components:
        comp_ev.append({
            "component_type": c.component_type,
            "label": c.label,
            "classification_state": c.classification_state,
            "confidence": c.confidence,
            "chain_ids": c.chain_ids,
            "evidence": [e.to_dict() for e in c.evidence],
        })

    return ObservableProvenance(
        schema="simforge/study-campaign/observable-provenance/v1",
        timestamp=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        simforge_version=version,
        argv=sys.argv[:],
        environment=_environment(),
        system_id=system.system_id,
        condition_id=system.condition_id,
        replicate_id=system.replicate_id,
        analysis_id=result.analysis_id,
        parameters=result.parameters,
        fit_selection=result.fit_selection,
        measure_selection=result.measure_selection,
        reference_frame=result.reference_frame,
        source_files={
            "production_trajectory_paths": system.production_trajectory_paths,
            "topology_path": system.topology_path,
            "structure_path": system.structure_path,
            "index_path": system.index_path,
            "trajectory_artifacts": [
                {"path": a.path, "stage": a.stage, "stage_confidence": a.stage_confidence}
                for a in system.trajectory_artifacts
            ],
        },
        source_fingerprints={k: v.to_dict() for k, v in system.source_fingerprints.items()},
        semantic_index=semantic_index.to_dict() if semantic_index else None,
        trajectory_view=trajectory_view.to_dict() if trajectory_view else None,
        component_evidence=comp_ev,
        trajectory_inspection=(
            system.trajectory_inspection.to_dict() if system.trajectory_inspection else None
        ),
        result=result.to_dict(),
        warnings=[w.to_dict() for w in result.warnings],
    )
