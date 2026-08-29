"""Observable plugin interface for the study-campaign layer.

An observable receives a *validated* :class:`AnalysisContext` (system record,
semantic index, trajectory view already resolved) and produces an
:class:`AnalysisResult`.  It must NOT rediscover topology, receptor, peptide,
water model or trajectory — study semantics belong to the campaign layer
(architectural rule §65).

Trajectory preprocessing is declared, not performed, via
``trajectory_requirements()`` (architectural rule §66).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Protocol, runtime_checkable

from analysis.campaign.models import (
    AnalysisResult, SemanticIndex, SystemRecord, TrajectoryRequirements,
    TrajectoryView,
)


@dataclass
class AnalysisContext:
    system: SystemRecord
    semantic_index: Optional[SemanticIndex]
    trajectory_view: TrajectoryView
    topology_path: str
    structure_path: Optional[str]
    output_dir: Path                     # per (system, analysis) directory
    parameters: dict
    gmx: str = "gmx"
    dry_run: bool = False


@runtime_checkable
class Observable(Protocol):
    id: str
    display_name: str
    description: str
    category: str
    required_components: tuple[str, ...]        # ComponentType values that must be RESOLVED
    supports_replicate_aggregation: bool

    def parameters_schema(self) -> dict: ...

    def trajectory_requirements(self, ctx_params: dict) -> TrajectoryRequirements: ...

    def execute(self, ctx: AnalysisContext) -> AnalysisResult: ...


class ObservableSpec:
    """Concrete base class observables subclass (simpler than the Protocol).

    Subclasses set the metadata as plain class attributes.
    """
    id: str = ""
    display_name: str = ""
    description: str = ""
    category: str = "structural"
    required_components: tuple[str, ...] = ()
    supports_replicate_aggregation: bool = True

    def parameters_schema(self) -> dict:
        return {}

    def trajectory_requirements(self, ctx_params: dict) -> TrajectoryRequirements:  # pragma: no cover
        raise NotImplementedError

    def execute(self, ctx: AnalysisContext) -> AnalysisResult:  # pragma: no cover
        raise NotImplementedError

    def to_dict(self) -> dict:
        return {
            "id": self.id,
            "display_name": self.display_name,
            "description": self.description,
            "category": self.category,
            "required_components": list(self.required_components),
            "supports_replicate_aggregation": self.supports_replicate_aggregation,
            "parameters_schema": self.parameters_schema(),
        }
