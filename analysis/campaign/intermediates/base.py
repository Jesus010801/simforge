"""Shared scientific intermediates — contracts (Phase 13).

An *intermediate* is a reproducible scientific artifact needed by one or more
consumers (observables or other intermediates) but not necessarily shown to
a user: a selection-centre series today; a membrane reference frame, a
leaflet assignment or a pore axis later.  It has an explicit identity and
provenance — never a scratch file.

Identity (all sha256 tokens of canonical JSON, like observable tokens):

* ``definition_identity`` — what is computed: id, version, parameters,
  selection / annotation / topology evidence, backend, *and the identities of
  the intermediates it depends on*; no paths, times or consumers;
* ``input_identity`` — which data it was computed on: the coordinate view
  identity (which encodes the source-trajectory fingerprints) for
  view/trajectory-scoped intermediates, the structure fingerprint for static
  ones;
* ``artifact_identity`` — both together; the key of the persistent store.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

from analysis.campaign.models import (
    ResultArray, SemanticIndex, StorageRef, SystemRecord, TrajectoryRequirements, TrajectoryView,
)

INTERMEDIATE_SCHEMA = "simforge/intermediate/v1"


class IntermediateScope:
    """What data an intermediate depends on (decides its input identity)."""
    STRUCTURE_STATIC = "structure_static"    # reference structure only
    VIEW = "view"                            # one coordinate view of the trajectory


class IntermediateStatus:
    COMPUTED = "computed"
    CACHED = "cached"
    PLANNED = "planned"                      # dry run: would be computed
    BLOCKED = "blocked"                      # applicability / policy / dependency blocked
    FAILED = "failed"

    USABLE = (COMPUTED, CACHED)


@dataclass(frozen=True)
class IntermediateRequest:
    """A parameterised dependency.  Its identity never mentions who asked:
    two consumers asking the same thing are one DAG node."""
    id: str
    parameters: tuple[tuple[str, Any], ...] = ()
    #: dependency contract: a consumer may pin the version it was written for
    version: Optional[int] = None

    @classmethod
    def make(cls, id: str, parameters: Optional[dict] = None,
             version: Optional[int] = None) -> "IntermediateRequest":
        return cls(id, tuple(sorted((parameters or {}).items())), version)

    @property
    def params(self) -> dict:
        return dict(self.parameters)

    @property
    def key(self) -> str:
        """Node key: id + parameters (a pinned version is a contract, checked
        by the planner, not a different request)."""
        if not self.parameters:
            return self.id
        return f"{self.id}(" + ",".join(f"{k}={v}" for k, v in self.parameters) + ")"

    def to_dict(self) -> dict:
        return {"id": self.id, "parameters": self.params, "version": self.version,
                "key": self.key}


@dataclass
class IntermediateContext:
    """What an intermediate is computed with; mirrors ``AnalysisContext``."""
    system: SystemRecord
    semantic_index: Optional[SemanticIndex]
    trajectory_view: Optional[TrajectoryView]
    topology_path: str
    structure_path: Optional[str]
    parameters: dict
    work_dir: Path                           # the build directory of this artifact
    gmx: str = "gmx"
    dry_run: bool = False
    time_index_cache_dir: Optional[Path] = None
    #: resolved dependencies of *this* intermediate, by request key
    dependencies: dict = field(default_factory=dict)


@dataclass
class IntermediateResult:
    intermediate_id: str
    version: int
    request: dict
    status: str
    reason: str = ""
    definition_identity: Optional[str] = None
    input_identity: Optional[str] = None
    artifact_identity: Optional[str] = None
    definition_evidence: dict = field(default_factory=dict)
    view_ref: Optional[str] = None
    arrays: list[ResultArray] = field(default_factory=list)
    #: non-array artifacts (structured mappings …): StorageRef with attrs["schema"]
    artifacts: list[StorageRef] = field(default_factory=list)
    dependencies: list[dict] = field(default_factory=list)
    provenance: dict = field(default_factory=dict)
    entry: Optional[str] = None              # published store entry (absolute path)

    @property
    def usable(self) -> bool:
        return self.status in IntermediateStatus.USABLE

    def reference(self) -> dict:
        """What a consumer records about this dependency (no payload)."""
        return {"id": self.intermediate_id, "version": self.version,
                "key": self.request.get("key"), "status": self.status,
                "definition_identity": self.definition_identity,
                "input_identity": self.input_identity,
                "artifact_identity": self.artifact_identity, "reason": self.reason}

    def array(self, name: str) -> Optional[ResultArray]:
        return next((a for a in self.arrays if a.name == name), None)

    def to_dict(self) -> dict:
        return {"intermediate_id": self.intermediate_id, "version": self.version,
                "request": self.request, "status": self.status, "reason": self.reason,
                "definition_identity": self.definition_identity,
                "input_identity": self.input_identity,
                "artifact_identity": self.artifact_identity,
                "definition_evidence": self.definition_evidence, "view_ref": self.view_ref,
                "arrays": [a.to_dict() for a in self.arrays],
                "artifacts": [a.to_dict() for a in self.artifacts],
                "dependencies": self.dependencies, "provenance": self.provenance,
                "entry": self.entry}


class IntermediateSpec:
    """Base contract.  Subclasses set the class attributes and implement
    ``definition_evidence`` / ``compute``; everything else has safe defaults."""
    id: str = ""
    version: int = 1
    purpose: str = ""                        # ObservablePurpose.* (policy purpose of its view)
    scope: str = IntermediateScope.VIEW
    description: str = ""

    def parameters_schema(self) -> dict:
        return {}

    def validate(self, params: dict) -> list[str]:
        """Parameter errors (cheap, no data access)."""
        unknown = sorted(set(params) - set(self.parameters_schema()))
        return [f"unknown parameter(s) {unknown}"] if unknown else []

    def dependencies(self, params: dict) -> list[IntermediateRequest]:
        return []

    def applicability(self, system: SystemRecord, params: dict,
                      semantic_index: Optional[SemanticIndex]) -> dict:
        return {"applicable": True, "reasons": []}

    def trajectory_requirements(self, params: dict) -> Optional[TrajectoryRequirements]:
        """The coordinate view it reads (judged by the Phase 5 policy with
        ``purpose``); None for structure-static intermediates."""
        return None

    def output_schema(self, params: dict) -> list[ResultArray]:
        return []

    def definition_evidence(self, ctx: IntermediateContext) -> Optional[dict]:
        raise NotImplementedError

    def compute(self, ctx: IntermediateContext) -> tuple[list[ResultArray], list[StorageRef], dict]:
        """Write artifacts under ``ctx.work_dir``; return (arrays, artifacts,
        provenance).  Must not publish anything itself."""
        raise NotImplementedError

    def to_dict(self) -> dict:
        return {"id": self.id, "version": self.version, "purpose": self.purpose,
                "scope": self.scope, "description": self.description,
                "parameters_schema": self.parameters_schema()}
