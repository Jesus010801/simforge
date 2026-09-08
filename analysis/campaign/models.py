"""Data models for the comparative MD **study campaign** analysis layer.

Design rules (consistent with ``analysis/fel`` and ``analysis/md``):

* Plain ``@dataclass`` only — no third-party model framework, no import-time cost.
* Every model that is written to disk implements ``to_dict()`` returning
  JSON-serialisable primitives.  Models that are also *read back* implement a
  ``from_dict()`` classmethod.
* Enums are ``str``-valued so they serialise transparently and compare to plain
  strings in tests.
* Scientific ambiguity is represented explicitly (``ClassificationState``,
  ``Ambiguity``) — never silently resolved.

The :class:`StudyManifest` is the stable contract between discovery, structural
inference, trajectory validation, observable execution and (future) replicate
aggregation.  Treat changes to it as schema changes and bump
``MANIFEST_SCHEMA_VERSION``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

MANIFEST_SCHEMA_VERSION = "simforge/study-campaign/v1"


# ═══════════════════════════════════════════════════════════════════════════════
# Enumerations (str-valued)
# ═══════════════════════════════════════════════════════════════════════════════

class Severity:
    INFO = "info"                    # diagnostic only — does not affect usability
    WARN = "warn"                    # usable but noteworthy
    REVIEW = "review_required"       # scientific ambiguity blocks some/all analyses
    ERROR = "error"                  # cannot safely use the system

    ORDER = {INFO: 0, WARN: 1, REVIEW: 2, ERROR: 3}


class ClassificationState:
    """State of a semantic identity that SimForge tried to infer."""
    RESOLVED = "resolved"            # inferred with sufficient confidence
    AMBIGUOUS = "ambiguous"          # multiple interpretations are compatible
    UNKNOWN = "unknown"              # no evidence at all
    REVIEW_REQUIRED = "review_required"  # ambiguous AND blocks a requested analysis


class ValidationState:
    VALID = "valid"
    WARNING = "warning"
    INVALID = "invalid"


class AnalysisStatus:
    SUCCESS = "success"
    SKIPPED = "skipped"             # requirements not available for this system
    FAILED = "failed"              # execution error (tool failure etc.)
    CACHED = "cached"              # reused a prior result
    REVIEW_REQUIRED = "review_required"  # semantic identity unresolved
    PLANNED = "planned"           # dry-run: command constructed, not executed
    ERROR = "error"               # invalid request (unknown id, bad params)


class ComponentType:
    RECEPTOR = "receptor"
    PEPTIDE = "peptide"
    LIGAND = "ligand"
    PROTEIN_PARTNER = "protein_partner"
    NUCLEIC_ACID = "nucleic_acid"
    COMPLEX = "complex"
    MEMBRANE = "membrane"
    WATER = "water"
    IONS = "ions"
    COFACTOR = "cofactor"
    OTHER = "other"
    UNKNOWN = "unknown"


class TrajectoryStage:
    """MD-workflow stage a trajectory belongs to.

    Only ``PRODUCTION`` trajectories feed observables.  Equilibration /
    minimisation / preparation trajectories are recorded on the system for
    transparency but never analysed, and their frames/timing never contribute
    to system-level production metadata.
    """
    PRODUCTION = "production"
    EQUILIBRATION_NVT = "equilibration_nvt"
    EQUILIBRATION_NPT = "equilibration_npt"
    EQUILIBRATION_OTHER = "equilibration_other"
    MINIMIZATION = "minimization"
    PREPARATION = "preparation"
    DERIVED = "derived"              # a PBC/fit derivative of another trajectory
    UNKNOWN = "unknown"

    EQUILIBRATION = {EQUILIBRATION_NVT, EQUILIBRATION_NPT, EQUILIBRATION_OTHER}
    NON_PRODUCTION = EQUILIBRATION | {MINIMIZATION, PREPARATION, DERIVED}


class ViewKind:
    """First-class derived trajectory representations."""
    RAW = "raw"
    WHOLE = "whole"                 # molecules made whole across PBC
    NOJUMP = "nojump"              # unwrapped, no PBC jumps (diffusion-safe)
    CENTERED = "centered"          # a semantic target centred in the box
    FITTED = "fitted"             # rot+trans least-squares fit to a reference


# ═══════════════════════════════════════════════════════════════════════════════
# Primitive records
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class CampaignWarning:
    code: str
    message: str
    severity: str = Severity.WARN     # info | warn | error
    scope: Optional[str] = None       # e.g. system_id or "study"

    def to_dict(self) -> dict:
        return {
            "code": self.code,
            "message": self.message,
            "severity": self.severity,
            "scope": self.scope,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "CampaignWarning":
        return cls(
            code=d["code"], message=d["message"],
            severity=d.get("severity", Severity.WARN), scope=d.get("scope"),
        )


@dataclass
class SourceFingerprint:
    """Reproducibility fingerprint for one source file.

    ``mode == "fast"``  : ``size`` + ``mtime_ns`` + sha256 of the first & last
                          ``FAST_EDGE_BYTES`` — cheap, safe for multi-GB
                          trajectories, but *not* collision proof.
    ``mode == "strong"``: full-content sha256 (``digest``).
    """
    path: str
    size_bytes: int
    mtime_ns: int
    mode: str                        # "fast" | "strong"
    digest: str                      # hex; the primary identity token
    edge_bytes: Optional[int] = None
    note: str = ""

    def to_dict(self) -> dict:
        return {
            "path": self.path,
            "size_bytes": self.size_bytes,
            "mtime_ns": self.mtime_ns,
            "mode": self.mode,
            "digest": self.digest,
            "edge_bytes": self.edge_bytes,
            "note": self.note,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "SourceFingerprint":
        return cls(
            path=d["path"], size_bytes=d["size_bytes"], mtime_ns=d["mtime_ns"],
            mode=d["mode"], digest=d["digest"],
            edge_bytes=d.get("edge_bytes"), note=d.get("note", ""),
        )

    @property
    def short(self) -> str:
        return self.digest[:12] if self.digest else "0" * 12


@dataclass
class ComponentEvidence:
    """One weighted piece of evidence supporting a component classification."""
    kind: str                        # e.g. "topology_molecule_name", "independent_chain"
    detail: str
    weight: float                    # signed contribution, roughly 0..1
    value: Any = None

    def to_dict(self) -> dict:
        return {"kind": self.kind, "detail": self.detail,
                "weight": self.weight, "value": self.value}

    @classmethod
    def from_dict(cls, d: dict) -> "ComponentEvidence":
        return cls(kind=d["kind"], detail=d["detail"],
                   weight=d["weight"], value=d.get("value"))


@dataclass
class MolecularComponent:
    """A semantically meaningful part of the simulated system."""
    component_type: str              # ComponentType.*
    label: str                       # human label, e.g. "receptor", "peptide"
    selection: str = ""              # semantic selection string (see selections.py)
    chain_ids: list[str] = field(default_factory=list)
    resnames: list[str] = field(default_factory=list)
    residue_count: Optional[int] = None
    atom_count: Optional[int] = None
    molecule_names: list[str] = field(default_factory=list)
    confidence: float = 0.0
    classification_state: str = ClassificationState.UNKNOWN
    study_role: Optional[str] = None       # e.g. "receptor (assumed)"; distinct from biological identity
    evidence: list[ComponentEvidence] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "component_type": self.component_type,
            "label": self.label,
            "selection": self.selection,
            "chain_ids": self.chain_ids,
            "resnames": self.resnames,
            "residue_count": self.residue_count,
            "atom_count": self.atom_count,
            "molecule_names": self.molecule_names,
            "confidence": round(self.confidence, 4),
            "classification_state": self.classification_state,
            "study_role": self.study_role,
            "evidence": [e.to_dict() for e in self.evidence],
            "warnings": self.warnings,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "MolecularComponent":
        return cls(
            component_type=d["component_type"],
            label=d["label"],
            selection=d.get("selection", ""),
            chain_ids=list(d.get("chain_ids", [])),
            resnames=list(d.get("resnames", [])),
            residue_count=d.get("residue_count"),
            atom_count=d.get("atom_count"),
            molecule_names=list(d.get("molecule_names", [])),
            confidence=d.get("confidence", 0.0),
            classification_state=d.get("classification_state", ClassificationState.UNKNOWN),
            study_role=d.get("study_role"),
            evidence=[ComponentEvidence.from_dict(e) for e in d.get("evidence", [])],
            warnings=list(d.get("warnings", [])),
        )


# ═══════════════════════════════════════════════════════════════════════════════
# Semantic index
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class SemanticIndexGroup:
    name: str
    n_atoms: int
    n_residues: Optional[int] = None
    chain_ids: list[str] = field(default_factory=list)
    source: str = ""                 # how the selection was derived
    selection_logic: str = ""        # human-readable selection description

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "n_atoms": self.n_atoms,
            "n_residues": self.n_residues,
            "chain_ids": self.chain_ids,
            "source": self.source,
            "selection_logic": self.selection_logic,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "SemanticIndexGroup":
        return cls(
            name=d["name"], n_atoms=d["n_atoms"], n_residues=d.get("n_residues"),
            chain_ids=list(d.get("chain_ids", [])), source=d.get("source", ""),
            selection_logic=d.get("selection_logic", ""),
        )


@dataclass
class SemanticIndex:
    path: Optional[str] = None       # generated .ndx path
    source_structure: Optional[str] = None
    groups: list[SemanticIndexGroup] = field(default_factory=list)
    gmx_commands: list[list[str]] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def group_names(self) -> list[str]:
        return [g.name for g in self.groups]

    def to_dict(self) -> dict:
        return {
            "path": self.path,
            "source_structure": self.source_structure,
            "groups": [g.to_dict() for g in self.groups],
            "gmx_commands": self.gmx_commands,
            "warnings": self.warnings,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "SemanticIndex":
        return cls(
            path=d.get("path"), source_structure=d.get("source_structure"),
            groups=[SemanticIndexGroup.from_dict(g) for g in d.get("groups", [])],
            gmx_commands=[list(c) for c in d.get("gmx_commands", [])],
            warnings=list(d.get("warnings", [])),
        )


# ═══════════════════════════════════════════════════════════════════════════════
# Trajectory records
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class TrajectoryArtifact:
    """One trajectory file, classified by MD-workflow stage.

    Every trajectory found in a system's directory becomes one of these — the
    manifest represents the *whole* workflow, not just what will be analysed.
    """
    path: str
    stage: str = TrajectoryStage.UNKNOWN
    stage_confidence: float = 0.0
    stage_evidence: list[str] = field(default_factory=list)
    topology_path: Optional[str] = None
    structure_path: Optional[str] = None
    fingerprint: Optional[SourceFingerprint] = None
    n_frames: Optional[int] = None
    start_time_ps: Optional[float] = None
    end_time_ps: Optional[float] = None
    dt_ps: Optional[float] = None
    part_index: Optional[int] = None          # continuation ordering, when known
    derived_from: Optional[str] = None        # for DERIVED: the source trajectory

    def to_dict(self) -> dict:
        return {
            "path": self.path,
            "stage": self.stage,
            "stage_confidence": round(self.stage_confidence, 3),
            "stage_evidence": self.stage_evidence,
            "topology_path": self.topology_path,
            "structure_path": self.structure_path,
            "fingerprint": self.fingerprint.to_dict() if self.fingerprint else None,
            "n_frames": self.n_frames,
            "start_time_ps": self.start_time_ps,
            "end_time_ps": self.end_time_ps,
            "dt_ps": self.dt_ps,
            "part_index": self.part_index,
            "derived_from": self.derived_from,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "TrajectoryArtifact":
        fp = d.get("fingerprint")
        return cls(
            path=d["path"],
            stage=d.get("stage", TrajectoryStage.UNKNOWN),
            stage_confidence=d.get("stage_confidence", 0.0),
            stage_evidence=list(d.get("stage_evidence", [])),
            topology_path=d.get("topology_path"),
            structure_path=d.get("structure_path"),
            fingerprint=SourceFingerprint.from_dict(fp) if fp else None,
            n_frames=d.get("n_frames"),
            start_time_ps=d.get("start_time_ps"),
            end_time_ps=d.get("end_time_ps"),
            dt_ps=d.get("dt_ps"),
            part_index=d.get("part_index"),
            derived_from=d.get("derived_from"),
        )


@dataclass
class TrajectorySegment:
    path: str
    fingerprint: Optional[SourceFingerprint] = None
    n_frames: Optional[int] = None
    start_time_ps: Optional[float] = None
    end_time_ps: Optional[float] = None
    dt_ps: Optional[float] = None
    order_index: Optional[int] = None
    order_evidence: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "path": self.path,
            "fingerprint": self.fingerprint.to_dict() if self.fingerprint else None,
            "n_frames": self.n_frames,
            "start_time_ps": self.start_time_ps,
            "end_time_ps": self.end_time_ps,
            "dt_ps": self.dt_ps,
            "order_index": self.order_index,
            "order_evidence": self.order_evidence,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "TrajectorySegment":
        fp = d.get("fingerprint")
        return cls(
            path=d["path"],
            fingerprint=SourceFingerprint.from_dict(fp) if fp else None,
            n_frames=d.get("n_frames"),
            start_time_ps=d.get("start_time_ps"),
            end_time_ps=d.get("end_time_ps"),
            dt_ps=d.get("dt_ps"),
            order_index=d.get("order_index"),
            order_evidence=list(d.get("order_evidence", [])),
        )


@dataclass
class TrajectoryInspection:
    """Lightweight metadata gathered *before* any analysis runs.

    Populated by :mod:`analysis.campaign.trajectory.inspector` (``gmx check``).
    Fields left ``None`` mean "could not determine" — never guessed.
    """
    trajectory_paths: list[str] = field(default_factory=list)
    topology_path: Optional[str] = None
    structure_path: Optional[str] = None

    n_frames: Optional[int] = None
    start_time_ps: Optional[float] = None
    end_time_ps: Optional[float] = None
    dt_ps: Optional[float] = None
    total_duration_ps: Optional[float] = None
    tpr_duration_ps: Optional[float] = None

    box_present: Optional[bool] = None
    box_vectors: Optional[list[float]] = None
    precision: Optional[str] = None          # "single" | "double" | None

    atom_count_trajectory: Optional[int] = None
    atom_count_topology: Optional[int] = None
    atom_count_match: Optional[bool] = None

    membrane_present: Optional[bool] = None
    membrane_orientation_note: str = ""

    segments: list[TrajectorySegment] = field(default_factory=list)
    segment_overlap: Optional[bool] = None
    segment_ordering: str = "single"        # single | confident | uncertain

    coordinate_jump_warnings: list[str] = field(default_factory=list)
    tool: str = ""
    raw_tool_output: str = ""
    warnings: list[CampaignWarning] = field(default_factory=list)
    inspected: bool = False                 # False => tool unavailable / not run

    def to_dict(self) -> dict:
        return {
            "trajectory_paths": self.trajectory_paths,
            "topology_path": self.topology_path,
            "structure_path": self.structure_path,
            "n_frames": self.n_frames,
            "start_time_ps": self.start_time_ps,
            "end_time_ps": self.end_time_ps,
            "dt_ps": self.dt_ps,
            "total_duration_ps": self.total_duration_ps,
            "tpr_duration_ps": self.tpr_duration_ps,
            "box_present": self.box_present,
            "box_vectors": self.box_vectors,
            "precision": self.precision,
            "atom_count_trajectory": self.atom_count_trajectory,
            "atom_count_topology": self.atom_count_topology,
            "atom_count_match": self.atom_count_match,
            "membrane_present": self.membrane_present,
            "membrane_orientation_note": self.membrane_orientation_note,
            "segments": [s.to_dict() for s in self.segments],
            "segment_overlap": self.segment_overlap,
            "segment_ordering": self.segment_ordering,
            "coordinate_jump_warnings": self.coordinate_jump_warnings,
            "tool": self.tool,
            "inspected": self.inspected,
            "warnings": [w.to_dict() for w in self.warnings],
        }

    @classmethod
    def from_dict(cls, d: dict) -> "TrajectoryInspection":
        return cls(
            trajectory_paths=list(d.get("trajectory_paths", [])),
            topology_path=d.get("topology_path"),
            structure_path=d.get("structure_path"),
            n_frames=d.get("n_frames"),
            start_time_ps=d.get("start_time_ps"),
            end_time_ps=d.get("end_time_ps"),
            dt_ps=d.get("dt_ps"),
            total_duration_ps=d.get("total_duration_ps"),
            tpr_duration_ps=d.get("tpr_duration_ps"),
            box_present=d.get("box_present"),
            box_vectors=d.get("box_vectors"),
            precision=d.get("precision"),
            atom_count_trajectory=d.get("atom_count_trajectory"),
            atom_count_topology=d.get("atom_count_topology"),
            atom_count_match=d.get("atom_count_match"),
            membrane_present=d.get("membrane_present"),
            membrane_orientation_note=d.get("membrane_orientation_note", ""),
            segments=[TrajectorySegment.from_dict(s) for s in d.get("segments", [])],
            segment_overlap=d.get("segment_overlap"),
            segment_ordering=d.get("segment_ordering", "single"),
            coordinate_jump_warnings=list(d.get("coordinate_jump_warnings", [])),
            tool=d.get("tool", ""),
            inspected=d.get("inspected", False),
            warnings=[CampaignWarning.from_dict(w) for w in d.get("warnings", [])],
        )


@dataclass
class TrajectoryRequirements:
    """Typed, observable-declared preprocessing needs.

    The orchestration layer resolves these into a concrete
    :class:`TrajectoryView`.  Observables never embed shell recipes.
    """
    requires_whole_molecules: bool = False
    requires_nojump: bool = False
    centering_target: Optional[str] = None    # semantic group name, e.g. "Receptor"
    fit_selection: Optional[str] = None       # semantic group name for LSQ fit
    minimum_image_distances: bool = False     # observable does its own PBC minimum image
    rationale: str = ""

    def view_kind(self) -> str:
        if self.fit_selection:
            return ViewKind.FITTED
        if self.centering_target:
            return ViewKind.CENTERED
        if self.requires_nojump:
            return ViewKind.NOJUMP
        if self.requires_whole_molecules:
            return ViewKind.WHOLE
        return ViewKind.RAW

    def cache_token(self) -> str:
        return (
            f"whole={int(self.requires_whole_molecules)}"
            f":nojump={int(self.requires_nojump)}"
            f":center={self.centering_target or '-'}"
            f":fit={self.fit_selection or '-'}"
        )

    def to_dict(self) -> dict:
        return {
            "requires_whole_molecules": self.requires_whole_molecules,
            "requires_nojump": self.requires_nojump,
            "centering_target": self.centering_target,
            "fit_selection": self.fit_selection,
            "minimum_image_distances": self.minimum_image_distances,
            "view_kind": self.view_kind(),
            "rationale": self.rationale,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "TrajectoryRequirements":
        return cls(
            requires_whole_molecules=d.get("requires_whole_molecules", False),
            requires_nojump=d.get("requires_nojump", False),
            centering_target=d.get("centering_target"),
            fit_selection=d.get("fit_selection"),
            minimum_image_distances=d.get("minimum_image_distances", False),
            rationale=d.get("rationale", ""),
        )


@dataclass
class PreprocessingOperation:
    operation: str                   # "make_whole", "nojump", "center", "fit"
    tool: str
    command: list[str]
    stdin: Optional[str] = None
    semantic_selection: Optional[str] = None
    input_trajectory: Optional[str] = None
    input_topology: Optional[str] = None
    output: Optional[str] = None
    reason: str = ""
    validation_note: str = ""
    returncode: Optional[int] = None
    stdout_tail: str = ""
    stderr_tail: str = ""
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "operation": self.operation,
            "tool": self.tool,
            "command": self.command,
            "stdin": self.stdin,
            "semantic_selection": self.semantic_selection,
            "input_trajectory": self.input_trajectory,
            "input_topology": self.input_topology,
            "output": self.output,
            "reason": self.reason,
            "validation_note": self.validation_note,
            "returncode": self.returncode,
            "stdout_tail": self.stdout_tail,
            "stderr_tail": self.stderr_tail,
            "warnings": self.warnings,
        }


@dataclass
class TrajectoryView:
    kind: str                        # ViewKind.*
    requirements: TrajectoryRequirements
    path: Optional[str] = None       # derived trajectory; None for RAW passthrough
    topology_path: Optional[str] = None
    source_fingerprints: list[SourceFingerprint] = field(default_factory=list)
    operations: list[PreprocessingOperation] = field(default_factory=list)
    cache_key: str = ""
    reused: bool = False
    warnings: list[CampaignWarning] = field(default_factory=list)
    safe: bool = True                # False => could not build a scientifically safe path

    def to_dict(self) -> dict:
        return {
            "kind": self.kind,
            "requirements": self.requirements.to_dict(),
            "path": self.path,
            "topology_path": self.topology_path,
            "source_fingerprints": [f.to_dict() for f in self.source_fingerprints],
            "operations": [o.to_dict() for o in self.operations],
            "cache_key": self.cache_key,
            "reused": self.reused,
            "safe": self.safe,
            "warnings": [w.to_dict() for w in self.warnings],
        }


# ═══════════════════════════════════════════════════════════════════════════════
# System / condition / manifest
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class SystemRecord:
    """One independent simulated system (one replicate of one condition)."""
    system_id: str
    condition_id: str
    replicate_id: str

    # ── semantic identity (may be unknown) ──────────────────────────────────
    receptor: Optional[str] = None
    partner: Optional[str] = None
    partner_type: Optional[str] = None            # peptide | ligand | protein | ...
    water_model: Optional[str] = None
    forcefield: Optional[str] = None

    # ── arbitrary extra condition dimensions ────────────────────────────────
    condition_dimensions: dict[str, str] = field(default_factory=dict)

    # ── source files ───────────────────────────────────────────────────────
    # Only PRODUCTION trajectories are analysed.  ``trajectory_artifacts`` holds
    # every trajectory in the system directory with its stage classification.
    production_trajectory_paths: list[str] = field(default_factory=list)
    trajectory_artifacts: list[TrajectoryArtifact] = field(default_factory=list)
    topology_path: Optional[str] = None            # PRODUCTION topology
    structure_path: Optional[str] = None           # PRODUCTION reference structure
    index_path: Optional[str] = None
    source_fingerprints: dict[str, SourceFingerprint] = field(default_factory=dict)

    # ── derived structure / trajectory understanding ───────────────────────
    membrane_present: Optional[bool] = None
    components: list[MolecularComponent] = field(default_factory=list)
    semantic_index: Optional[SemanticIndex] = None
    trajectory_inspection: Optional[TrajectoryInspection] = None

    # ── provenance of the association itself ───────────────────────────────
    discovery_evidence: list[str] = field(default_factory=list)
    association_evidence: list[str] = field(default_factory=list)
    classification_state: str = ClassificationState.UNKNOWN
    warnings: list[CampaignWarning] = field(default_factory=list)
    user_overridden: bool = False
    legacy_study: dict = field(default_factory=dict)

    # -- helpers ----------------------------------------------------------------
    @property
    def trajectory_paths(self) -> list[str]:
        """The trajectories observables operate on — PRODUCTION only."""
        return self.production_trajectory_paths

    @property
    def workflow_trajectory_paths(self) -> list[str]:
        """Non-production trajectories (equilibration / minimisation / derived)."""
        return [a.path for a in self.trajectory_artifacts
                if a.stage in TrajectoryStage.NON_PRODUCTION]

    def artifact(self, path: str) -> Optional[TrajectoryArtifact]:
        for a in self.trajectory_artifacts:
            if a.path == path:
                return a
        return None

    def component(self, ctype: str) -> Optional[MolecularComponent]:
        for c in self.components:
            if c.component_type == ctype:
                return c
        return None

    def resolved_component(self, ctype: str) -> Optional[MolecularComponent]:
        c = self.component(ctype)
        if c and c.classification_state == ClassificationState.RESOLVED:
            return c
        return None

    def to_dict(self) -> dict:
        return {
            "system_id": self.system_id,
            "condition_id": self.condition_id,
            "replicate_id": self.replicate_id,
            "receptor": self.receptor,
            "partner": self.partner,
            "partner_type": self.partner_type,
            "water_model": self.water_model,
            "forcefield": self.forcefield,
            "condition_dimensions": self.condition_dimensions,
            "production_trajectory_paths": self.production_trajectory_paths,
            "workflow_trajectory_paths": self.workflow_trajectory_paths,
            "trajectory_artifacts": [a.to_dict() for a in self.trajectory_artifacts],
            "topology_path": self.topology_path,
            "structure_path": self.structure_path,
            "index_path": self.index_path,
            "source_fingerprints": {k: v.to_dict() for k, v in self.source_fingerprints.items()},
            "membrane_present": self.membrane_present,
            "components": [c.to_dict() for c in self.components],
            "semantic_index": self.semantic_index.to_dict() if self.semantic_index else None,
            "trajectory_inspection": (
                self.trajectory_inspection.to_dict() if self.trajectory_inspection else None
            ),
            "discovery_evidence": self.discovery_evidence,
            "association_evidence": self.association_evidence,
            "classification_state": self.classification_state,
            "user_overridden": self.user_overridden,
            "legacy_study": self.legacy_study,
            "warnings": [w.to_dict() for w in self.warnings],
        }

    @classmethod
    def from_dict(cls, d: dict) -> "SystemRecord":
        si = d.get("semantic_index")
        ti = d.get("trajectory_inspection")
        return cls(
            system_id=d["system_id"],
            condition_id=d["condition_id"],
            replicate_id=d["replicate_id"],
            receptor=d.get("receptor"),
            partner=d.get("partner"),
            partner_type=d.get("partner_type"),
            water_model=d.get("water_model"),
            forcefield=d.get("forcefield"),
            condition_dimensions=dict(d.get("condition_dimensions", {})),
            production_trajectory_paths=list(
                d.get("production_trajectory_paths", d.get("trajectory_paths", []))
            ),
            trajectory_artifacts=[
                TrajectoryArtifact.from_dict(a) for a in d.get("trajectory_artifacts", [])
            ],
            topology_path=d.get("topology_path"),
            structure_path=d.get("structure_path"),
            index_path=d.get("index_path"),
            source_fingerprints={
                k: SourceFingerprint.from_dict(v)
                for k, v in d.get("source_fingerprints", {}).items()
            },
            membrane_present=d.get("membrane_present"),
            components=[MolecularComponent.from_dict(c) for c in d.get("components", [])],
            semantic_index=SemanticIndex.from_dict(si) if si else None,
            trajectory_inspection=TrajectoryInspection.from_dict(ti) if ti else None,
            discovery_evidence=list(d.get("discovery_evidence", [])),
            association_evidence=list(d.get("association_evidence", [])),
            classification_state=d.get("classification_state", ClassificationState.UNKNOWN),
            warnings=[CampaignWarning.from_dict(w) for w in d.get("warnings", [])],
            user_overridden=d.get("user_overridden", False),
            legacy_study=dict(d.get("legacy_study", {})),
        )


@dataclass
class ConditionRecord:
    """A group of replicate systems that share every condition dimension."""
    condition_id: str
    dimensions: dict[str, str] = field(default_factory=dict)
    system_ids: list[str] = field(default_factory=list)
    replicate_ids: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "condition_id": self.condition_id,
            "dimensions": self.dimensions,
            "system_ids": self.system_ids,
            "replicate_ids": self.replicate_ids,
            "n_replicates": len(self.system_ids),
        }

    @classmethod
    def from_dict(cls, d: dict) -> "ConditionRecord":
        return cls(
            condition_id=d["condition_id"],
            dimensions=dict(d.get("dimensions", {})),
            system_ids=list(d.get("system_ids", [])),
            replicate_ids=list(d.get("replicate_ids", [])),
        )


@dataclass
class Ambiguity:
    system_id: str
    kind: str                        # e.g. "component_classification", "topology_pairing"
    message: str
    options: list[str] = field(default_factory=list)
    resolution: Optional[str] = None  # set when a manifest is manually corrected

    def to_dict(self) -> dict:
        return {
            "system_id": self.system_id,
            "kind": self.kind,
            "message": self.message,
            "options": self.options,
            "resolution": self.resolution,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "Ambiguity":
        return cls(
            system_id=d["system_id"], kind=d["kind"], message=d["message"],
            options=list(d.get("options", [])), resolution=d.get("resolution"),
        )


@dataclass
class StudyManifest:
    """The machine-readable record of every discovery decision for one study."""
    study_root: str
    generated_utc: str
    simforge_version: str
    schema_version: str = MANIFEST_SCHEMA_VERSION
    discovery_settings: dict[str, Any] = field(default_factory=dict)
    condition_dimensions: list[str] = field(default_factory=list)
    conditions: list[ConditionRecord] = field(default_factory=list)
    systems: list[SystemRecord] = field(default_factory=list)
    unassigned_files: list[str] = field(default_factory=list)
    ambiguities: list[Ambiguity] = field(default_factory=list)
    warnings: list[CampaignWarning] = field(default_factory=list)

    # -- helpers ----------------------------------------------------------------
    def system(self, system_id: str) -> Optional[SystemRecord]:
        for s in self.systems:
            if s.system_id == system_id:
                return s
        return None

    def stats(self) -> dict:
        partners = sorted({s.partner for s in self.systems if s.partner})
        waters = sorted({s.water_model for s in self.systems if s.water_model})
        replicate_counts = {c.condition_id: len(c.system_ids) for c in self.conditions}
        return {
            "n_systems": len(self.systems),
            "n_conditions": len(self.conditions),
            "n_partners": len(partners),
            "n_water_models": len(waters),
            "partners": partners,
            "water_models": waters,
            "replicate_counts": replicate_counts,
            "balanced": len(set(replicate_counts.values())) <= 1 if replicate_counts else True,
            "n_trajectory_segments": sum(len(s.trajectory_paths) for s in self.systems),
            "n_ambiguities": len(self.ambiguities),
        }

    def to_dict(self) -> dict:
        return {
            "schema_version": self.schema_version,
            "study_root": self.study_root,
            "generated_utc": self.generated_utc,
            "simforge_version": self.simforge_version,
            "discovery_settings": self.discovery_settings,
            "condition_dimensions": self.condition_dimensions,
            "stats": self.stats(),
            "conditions": [c.to_dict() for c in self.conditions],
            "systems": [s.to_dict() for s in self.systems],
            "unassigned_files": self.unassigned_files,
            "ambiguities": [a.to_dict() for a in self.ambiguities],
            "warnings": [w.to_dict() for w in self.warnings],
        }

    @classmethod
    def from_dict(cls, d: dict) -> "StudyManifest":
        return cls(
            study_root=d["study_root"],
            generated_utc=d.get("generated_utc", ""),
            simforge_version=d.get("simforge_version", "unknown"),
            schema_version=d.get("schema_version", MANIFEST_SCHEMA_VERSION),
            discovery_settings=dict(d.get("discovery_settings", {})),
            condition_dimensions=list(d.get("condition_dimensions", [])),
            conditions=[ConditionRecord.from_dict(c) for c in d.get("conditions", [])],
            systems=[SystemRecord.from_dict(s) for s in d.get("systems", [])],
            unassigned_files=list(d.get("unassigned_files", [])),
            ambiguities=[Ambiguity.from_dict(a) for a in d.get("ambiguities", [])],
            warnings=[CampaignWarning.from_dict(w) for w in d.get("warnings", [])],
        )


# ═══════════════════════════════════════════════════════════════════════════════
# Validation
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class ObservableSystemStatus:
    system_id: str
    analysis_id: str
    state: str                       # ValidationState.* or "review_required"
    reason: str = ""
    remediation: str = ""

    def to_dict(self) -> dict:
        return {
            "system_id": self.system_id,
            "analysis_id": self.analysis_id,
            "state": self.state,
            "reason": self.reason,
            "remediation": self.remediation,
        }


@dataclass
class ValidationReport:
    study_state: str = ValidationState.VALID
    system_states: dict[str, str] = field(default_factory=dict)
    observable_statuses: list[ObservableSystemStatus] = field(default_factory=list)
    lines: list[str] = field(default_factory=list)          # human-readable report
    warnings: list[CampaignWarning] = field(default_factory=list)
    counts: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "study_state": self.study_state,
            "system_states": self.system_states,
            "observable_statuses": [o.to_dict() for o in self.observable_statuses],
            "counts": self.counts,
            "lines": self.lines,
            "warnings": [w.to_dict() for w in self.warnings],
        }


# ═══════════════════════════════════════════════════════════════════════════════
# Analysis requests / results
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class AnalysisRequest:
    analysis_id: str
    parameters: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {"analysis_id": self.analysis_id, "parameters": self.parameters}


@dataclass
class AnalysisResult:
    analysis_id: str
    system_id: str
    status: str                      # AnalysisStatus.*
    message: str = ""
    fit_selection: Optional[str] = None
    measure_selection: Optional[str] = None
    reference_frame: str = ""
    trajectory_view_kind: str = ViewKind.RAW
    parameters: dict[str, Any] = field(default_factory=dict)
    output_files: list[str] = field(default_factory=list)
    data_summary: dict[str, Any] = field(default_factory=dict)
    cache_key: str = ""
    cached: bool = False
    provenance_path: Optional[str] = None
    warnings: list[CampaignWarning] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "analysis_id": self.analysis_id,
            "system_id": self.system_id,
            "status": self.status,
            "message": self.message,
            "fit_selection": self.fit_selection,
            "measure_selection": self.measure_selection,
            "reference_frame": self.reference_frame,
            "trajectory_view_kind": self.trajectory_view_kind,
            "parameters": self.parameters,
            "output_files": self.output_files,
            "data_summary": self.data_summary,
            "cache_key": self.cache_key,
            "cached": self.cached,
            "provenance_path": self.provenance_path,
            "warnings": [w.to_dict() for w in self.warnings],
        }


@dataclass
class CampaignRunResult:
    study_root: str
    output_dir: str
    manifest: Optional[StudyManifest] = None
    validation: Optional[ValidationReport] = None
    requested_analyses: list[str] = field(default_factory=list)
    results: list[AnalysisResult] = field(default_factory=list)
    output_files: list[str] = field(default_factory=list)
    warnings: list[CampaignWarning] = field(default_factory=list)

    def error_warnings(self) -> list[CampaignWarning]:
        return [w for w in self.warnings if w.severity == Severity.ERROR]

    def to_dict(self) -> dict:
        return {
            "study_root": self.study_root,
            "output_dir": self.output_dir,
            "manifest_stats": self.manifest.stats() if self.manifest else None,
            "legacy_systems": [s.legacy_study for s in self.manifest.systems] if self.manifest else [],
            "validation": self.validation.to_dict() if self.validation else None,
            "requested_analyses": self.requested_analyses,
            "results": [r.to_dict() for r in self.results],
            "output_files": self.output_files,
            "warnings": [w.to_dict() for w in self.warnings],
        }
