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
    AnalysisResult, ResultArray, SemanticIndex, SystemRecord, TrajectoryRequirements,
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
    #: Phase 1 time-index cache shared with the diagnostics pass (timeline reuse)
    time_index_cache_dir: Optional[Path] = None


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
    #: Scientific coordinate purpose (``models.ObservablePurpose``; open string).
    #: Reserved for context-dependent preprocessing decisions — not acted on yet.
    purpose: Optional[str] = None
    #: ids of persistent annotations this observable is defined on (Phase 6);
    #: their identities belong in its definition evidence
    required_annotations: tuple[str, ...] = ()
    #: Phase 8: True only if results on *different* coordinate views are
    #: scientifically interchangeable for this observable.  Never assumed.
    view_invariant: bool = False

    def parameters_schema(self) -> dict:
        return {}

    def applicability(self, system: SystemRecord, params: Optional[dict] = None,
                      semantic_index: Optional[SemanticIndex] = None) -> dict:
        """Can this observable run on ``system``?  Never executes anything.

        ``{"applicable": bool, "reasons": [...], "missing_components": [...],
        "missing_annotations": [...]}`` — enough for a later capability listing.
        """
        from analysis.campaign.annotations import annotation_evidence
        from analysis.campaign.models import ClassificationState
        miss_c = [c for c in self.required_components
                  if not any(x.component_type == c and
                             x.classification_state == ClassificationState.RESOLVED
                             for x in system.components)]
        miss_a = [a for a in self.required_annotations
                  if annotation_evidence(system.annotations, a) is None]
        reasons = ([f"component '{c}' not RESOLVED" for c in miss_c]
                   + [f"annotation '{a}' not ACTIVE" for a in miss_a])
        return {"applicable": not reasons, "reasons": reasons,
                "missing_components": miss_c, "missing_annotations": miss_a}

    def annotation_evidence(self, ctx: "AnalysisContext") -> Optional[dict]:
        """``{id: evidence}`` for ``required_annotations`` — None if any is
        missing or not ACTIVE (an observable must not run on inactive intent).
        Include the result in :meth:`definition_evidence` so a changed
        annotation invalidates exactly the analyses that use it."""
        from analysis.campaign.annotations import annotation_evidence
        out: dict = {}
        for ann_id in self.required_annotations:
            ev = annotation_evidence(ctx.system.annotations, ann_id)
            if ev is None:
                return None
            out[ann_id] = ev
        return out

    def output_schema(self, params: Optional[dict] = None) -> list[ResultArray]:
        """The result arrays this observable will produce, declared before
        execution (same model as realised output, with ``storage=None``)."""
        return []

    def definition_evidence(self, ctx: "AnalysisContext") -> Optional[dict]:
        """Explicit evidence of what is being computed, for
        :meth:`definition_token`.  ``None`` = not enough evidence (no token).
        Must not contain output paths, timestamps or other run-local values."""
        return None

    def definition_token(self, ctx: "AnalysisContext") -> Optional[str]:
        evidence = self.definition_evidence(ctx)
        if evidence is None:
            return None
        from analysis.campaign.results import definition_token
        return definition_token(evidence)

    def dependencies(self, params: Optional[dict] = None) -> list[tuple[str, dict]]:
        """Reserved: ``(observable_id, parameters)`` this observable consumes.
        Not orchestrated yet; must stay empty until it is."""
        return []

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
            "purpose": self.purpose,
            "output_schema": [a.to_dict() for a in self.output_schema({})],
        }
