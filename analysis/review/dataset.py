"""``ReviewDataset`` — one reproducible, reopenable review-session description.

It *references* everything and embeds no trajectory or result arrays:

* sources (trajectory / topology / structure / index) by path + fingerprint;
* the display view and every analysis view (``TrajectoryView.to_dict()``);
* the review timeline (a Phase 1 ``FrameTimeIndex`` snapshot in the session
  store) and how the display frames relate to it;
* annotations, diagnostics and preprocessing decisions by identity;
* every *requested* observable instance — available ones with their
  ``ResultArray`` (storage re-pointed at an immutable store copy), blocked /
  not-applicable / unsupported / failed ones with their reason.

Scientific definitions stay in ``analysis.campaign``; this module only
sessionizes them.  Paths inside the review root are stored relative to the
session directory; external files keep absolute paths plus fingerprints.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

from analysis.campaign.models import (
    AnalysisStatus, ClassificationState, DecisionClass, ResultArray,
)

REVIEW_SCHEMA_VERSION = "simforge/review-dataset/v1"
SUPPORTED_SCHEMAS = (REVIEW_SCHEMA_VERSION,)
DATASET_FILENAME = "review_dataset.json"


class CapabilityState:
    """Why a requested capability is (not) available — existing vocabularies only.

    ``NOT_APPLICABLE`` (the system lacks what the observable is about) is kept
    apart from the blocked states (applicable in principle, but unresolved,
    ambiguous, needing review, or without an executor).
    """
    AVAILABLE = "available"
    REVIEW_REQUIRED = ClassificationState.REVIEW_REQUIRED
    AMBIGUOUS = ClassificationState.AMBIGUOUS
    UNSUPPORTED = DecisionClass.UNSUPPORTED
    REFUSED = DecisionClass.REFUSED                # policy refuses the transformation here
    NOT_APPLICABLE = "not_applicable"
    FAILED = AnalysisStatus.FAILED
    PLANNED = AnalysisStatus.PLANNED               # dry run only

    BLOCKED = (REVIEW_REQUIRED, AMBIGUOUS, UNSUPPORTED, REFUSED)


class Availability:
    """How an available result was obtained."""
    COMPUTED = "computed"                          # executed in this preparation
    CACHED = "cached"                              # Phase 8 compatible reuse
    EXPLICIT_IMPORT = "explicit_import"            # user-selected external series
    PLANNED = "planned"
    NONE = "none"


@dataclass
class ReviewResult:
    """One ``ResultArray`` of an observable, with storage in the session store."""
    array: ResultArray
    sync: dict[str, Any] = field(default_factory=dict)
    origin: dict[str, Any] = field(default_factory=dict)   # the campaign file it was frozen from

    def to_dict(self) -> dict:
        a = self.array
        return {"name": a.name, "quantity": a.quantity, "unit": a.unit,
                "axis_signature": list(a.axis_signature()), "view_ref": a.view_ref,
                "definition_token": a.definition_token, "sync": self.sync,
                "origin": self.origin, "array": a.to_dict()}

    @classmethod
    def from_dict(cls, d: dict) -> "ReviewResult":
        return cls(array=ResultArray.from_dict(d["array"]), sync=dict(d.get("sync", {})),
                   origin=dict(d.get("origin", {})))


@dataclass
class ObservableEntry:
    """One requested observable instance and what became of it."""
    instance_id: str
    observable: str
    parameters: dict[str, Any] = field(default_factory=dict)
    display_name: str = ""
    purpose: Optional[str] = None
    state: str = CapabilityState.AVAILABLE         # CapabilityState.*
    availability: str = Availability.NONE          # Availability.*
    execution_status: Optional[str] = None         # AnalysisStatus.* of the campaign result
    reason: str = ""
    applicability: dict[str, Any] = field(default_factory=dict)
    annotations_used: list[str] = field(default_factory=list)
    blocking_annotations: list[str] = field(default_factory=list)
    results: list[ReviewResult] = field(default_factory=list)
    definition_token: Optional[str] = None
    definition_identity: Optional[str] = None
    view_ref: Optional[str] = None
    reuse: Optional[dict[str, Any]] = None         # decision, reused_from, chosen summary, report ref
    provenance_ref: Optional[dict[str, Any]] = None
    externally_supplied: bool = False
    independently_verified: bool = True
    #: shared intermediates it was computed from (identities + status; no payload)
    dependencies: list[dict[str, Any]] = field(default_factory=list)

    @property
    def available(self) -> bool:
        return self.state == CapabilityState.AVAILABLE

    @property
    def blocked(self) -> bool:
        return self.state in CapabilityState.BLOCKED

    def to_dict(self) -> dict:
        return {
            "instance_id": self.instance_id, "observable": self.observable,
            "parameters": self.parameters, "display_name": self.display_name,
            "purpose": self.purpose, "state": self.state, "availability": self.availability,
            "execution_status": self.execution_status, "reason": self.reason,
            "applicability": self.applicability, "annotations_used": self.annotations_used,
            "blocking_annotations": self.blocking_annotations,
            "results": [r.to_dict() for r in self.results],
            "definition_token": self.definition_token,
            "definition_identity": self.definition_identity, "view_ref": self.view_ref,
            "reuse": self.reuse, "provenance_ref": self.provenance_ref,
            "externally_supplied": self.externally_supplied,
            "independently_verified": self.independently_verified,
            "dependencies": self.dependencies,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "ObservableEntry":
        return cls(
            instance_id=d["instance_id"], observable=d["observable"],
            parameters=dict(d.get("parameters", {})), display_name=d.get("display_name", ""),
            purpose=d.get("purpose"), state=d.get("state", CapabilityState.AVAILABLE),
            availability=d.get("availability", Availability.NONE),
            execution_status=d.get("execution_status"), reason=d.get("reason", ""),
            applicability=dict(d.get("applicability", {})),
            annotations_used=list(d.get("annotations_used", [])),
            blocking_annotations=list(d.get("blocking_annotations", [])),
            results=[ReviewResult.from_dict(r) for r in d.get("results", [])],
            definition_token=d.get("definition_token"),
            definition_identity=d.get("definition_identity"), view_ref=d.get("view_ref"),
            reuse=d.get("reuse"), provenance_ref=d.get("provenance_ref"),
            externally_supplied=d.get("externally_supplied", False),
            independently_verified=d.get("independently_verified", True),
            dependencies=list(d.get("dependencies", [])),
        )


@dataclass
class ReviewDataset:
    session_id: str
    identity_evidence: dict[str, Any]
    status: str = "prepared"                       # "prepared" | "dry_run"
    system: dict[str, Any] = field(default_factory=dict)
    sources: dict[str, Any] = field(default_factory=dict)
    display_view: dict[str, Any] = field(default_factory=dict)
    analysis_views: list[dict[str, Any]] = field(default_factory=list)
    timeline: dict[str, Any] = field(default_factory=dict)
    annotations: list[dict[str, Any]] = field(default_factory=list)
    diagnostics: dict[str, Any] = field(default_factory=dict)
    requested: list[dict[str, Any]] = field(default_factory=list)
    observables: list[ObservableEntry] = field(default_factory=list)
    provenance: dict[str, Any] = field(default_factory=dict)
    portability: dict[str, Any] = field(default_factory=dict)
    #: fields later layers will fill (viewer, window, stride) — never used here
    reserved: dict[str, Any] = field(default_factory=lambda: {
        "viewer": None, "window": None, "stride": None})
    warnings: list[str] = field(default_factory=list)
    schema_version: str = REVIEW_SCHEMA_VERSION

    # ── queries ────────────────────────────────────────────────────────────
    def observable(self, instance_id: str) -> Optional[ObservableEntry]:
        return next((o for o in self.observables if o.instance_id == instance_id), None)

    def available(self) -> list[ObservableEntry]:
        return [o for o in self.observables if o.available]

    def blocked(self) -> list[ObservableEntry]:
        return [o for o in self.observables if not o.available]

    def view(self, view_ref: Optional[str]) -> Optional[dict]:
        if view_ref and self.display_view.get("view_ref") == view_ref:
            return self.display_view
        return next((v for v in self.analysis_views if v.get("view_ref") == view_ref), None)

    def summary(self) -> dict:
        count = lambda pred: sum(1 for o in self.observables if pred(o))  # noqa: E731
        return {
            "session_id": self.session_id, "status": self.status,
            "system_id": self.system.get("system_id"),
            "frames": self.timeline.get("n_frames"),
            "time_range_ps": [self.timeline.get("start_time_ps"), self.timeline.get("end_time_ps")],
            "display": {"mode": (self.display_view.get("request") or {}).get("mode"),
                        "status": self.display_view.get("status")},
            "observables": {
                "requested": len(self.observables),
                "computed": count(lambda o: o.availability == Availability.COMPUTED),
                "cached": count(lambda o: o.availability == Availability.CACHED),
                "explicit_import": count(lambda o: o.availability == Availability.EXPLICIT_IMPORT),
                "planned": count(lambda o: o.state == CapabilityState.PLANNED),
                "blocked": count(lambda o: o.blocked),
                "not_applicable": count(lambda o: o.state == CapabilityState.NOT_APPLICABLE),
                "failed": count(lambda o: o.state == CapabilityState.FAILED),
            },
            "diagnostics": {k: self.diagnostics.get(k) for k in ("execution", "status", "counts")},
        }

    # ── serialisation ──────────────────────────────────────────────────────
    def to_dict(self) -> dict:
        return {
            "schema_version": self.schema_version, "session_id": self.session_id,
            "status": self.status, "identity_evidence": self.identity_evidence,
            "system": self.system, "sources": self.sources,
            "display_view": self.display_view, "analysis_views": self.analysis_views,
            "timeline": self.timeline, "annotations": self.annotations,
            "diagnostics": self.diagnostics, "requested": self.requested,
            "observables": [o.to_dict() for o in self.observables],
            "provenance": self.provenance, "portability": self.portability,
            "reserved": self.reserved, "warnings": self.warnings,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "ReviewDataset":
        return cls(
            session_id=d["session_id"], identity_evidence=dict(d.get("identity_evidence", {})),
            status=d.get("status", "prepared"), system=dict(d.get("system", {})),
            sources=dict(d.get("sources", {})), display_view=dict(d.get("display_view", {})),
            analysis_views=list(d.get("analysis_views", [])),
            timeline=dict(d.get("timeline", {})), annotations=list(d.get("annotations", [])),
            diagnostics=dict(d.get("diagnostics", {})), requested=list(d.get("requested", [])),
            observables=[ObservableEntry.from_dict(o) for o in d.get("observables", [])],
            provenance=dict(d.get("provenance", {})), portability=dict(d.get("portability", {})),
            reserved=dict(d.get("reserved", {"viewer": None, "window": None, "stride": None})),
            warnings=list(d.get("warnings", [])),
            schema_version=d.get("schema_version", ""),
        )
