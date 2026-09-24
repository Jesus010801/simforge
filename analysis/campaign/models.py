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
    #: exact 1-based atom indices (reference-structure order) when membership is
    #: known atom-by-atom (resolved ligands, explicit resolutions); None otherwise
    atom_ids: Optional[list[int]] = None

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
            "atom_ids": self.atom_ids,
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
            atom_ids=list(d["atom_ids"]) if d.get("atom_ids") is not None else None,
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
    #: fit mode for ``fit_selection``; build_view executes only "rot+trans"
    fit_mode: str = "rot+trans"
    #: requested component reconstruction (e.g. "cluster"); no executor yet —
    #: representable so policy can evaluate (and refuse to fake) it
    reconstruction: Optional[str] = None

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
            + (f":fitmode={self.fit_mode}" if self.fit_mode != "rot+trans" else "")
            + (f":reconstruct={self.reconstruction}" if self.reconstruction else "")
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
            "fit_mode": self.fit_mode,
            "reconstruction": self.reconstruction,
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
            fit_mode=d.get("fit_mode", "rot+trans"),
            reconstruction=d.get("reconstruction"),
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
    cache_key: str = ""              # view identity (see trajectory/preprocessor.py)
    reused: bool = False
    warnings: list[CampaignWarning] = field(default_factory=list)
    safe: bool = True                # False => could not build a scientifically safe path
    # ── Phase 3: reproducible, validated derived views (additive) ──────────
    build_status: str = ""           # ViewBuildStatus.*
    cache_status: str = ""           # ViewCacheStatus.*
    cache_rejections: list[str] = field(default_factory=list)  # why a cached entry was not trusted
    manifest_path: Optional[str] = None
    output_fingerprint: Optional[SourceFingerprint] = None
    time_index_ref: Optional[dict] = None   # reference to the view's own FrameTimeIndex
    # ── Phase 5: scientific authorisation (additive) ───────────────────────
    #: True only when the requirements were evaluated by the policy engine;
    #: direct/legacy build_view calls stay False (never labelled as validated)
    policy_planned: bool = False
    decisions: list[dict] = field(default_factory=list)   # PreprocessingDecision.to_dict()
    #: purpose / intent / diagnostics report reference the plan was made with
    policy_ref: Optional[dict] = None

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
            "build_status": self.build_status,
            "cache_status": self.cache_status,
            "cache_rejections": self.cache_rejections,
            "manifest_path": self.manifest_path,
            "output_fingerprint": (self.output_fingerprint.to_dict()
                                   if self.output_fingerprint else None),
            "time_index_ref": self.time_index_ref,
            "policy_planned": self.policy_planned,
            "decisions": self.decisions,
            "policy_ref": self.policy_ref,
        }


class ViewBuildStatus:
    RAW = "raw"                      # passthrough of the source trajectory
    PLANNED = "planned"              # dry run: commands only
    BUILT = "built"                  # generated now, validated, manifest written
    REUSED = "reused"                # validated cache entry reused
    FAILED = "failed"                # no trustworthy view (see warnings)


class ViewCacheStatus:
    HIT = "hit"                      # entry validated against its manifest
    MISS = "miss"                    # no entry for this identity
    INVALID = "invalid"              # entry existed but failed validation → regenerated
    BYPASSED = "bypassed"            # force=True / dry run / raw


# ═══════════════════════════════════════════════════════════════════════════════
# Per-frame time index
# ═══════════════════════════════════════════════════════════════════════════════

TIME_INDEX_SCHEMA_VERSION = "simforge/trajectory-time-index/v1"


@dataclass
class FrameTimeIndex:
    """The actual simulation timestamp of every frame of *one* trajectory file.

    Frames are kept exactly in file order: duplicate and non-monotonic
    timestamps are *reported*, never removed or reordered.  Times are in ps.

    The per-frame arrays (``times_ps``, ``steps``, ``boxes``) can be large, so
    ``to_dict()`` omits them unless ``include_frames=True``; the cache stores
    them in a sidecar ``frames.npz`` (see ``trajectory/time_index.py``).
    """
    trajectory_path: str
    fingerprint: Optional[SourceFingerprint]
    backend: str                              # e.g. "gromacs-trjconv-gro"
    backend_version: Optional[str]            # e.g. GROMACS "2025.2"
    times_ps: list[float] = field(default_factory=list)
    steps: list[Optional[int]] = field(default_factory=list)
    boxes: Optional[list[list[float]]] = None  # per frame: 9 box elements (nm), GRO order
    schema_version: str = TIME_INDEX_SCHEMA_VERSION
    time_unit_reported: str = "ps"            # unit the backend reported the times in
    time_resolution_ps: Optional[float] = None  # textual resolution of the reported times

    # ── validation ─────────────────────────────────────────────────────────
    state: str = ValidationState.VALID        # ValidationState.*
    strictly_increasing: Optional[bool] = None
    non_decreasing: Optional[bool] = None
    duplicate_groups: list[list[int]] = field(default_factory=list)      # frames sharing one time
    non_monotonic_steps: list[list[int]] = field(default_factory=list)   # [i, i+1] where t[i+1] < t[i]
    non_finite_frames: list[int] = field(default_factory=list)
    truncated: bool = False
    truncation_evidence: str = ""
    backend_reported_frames: Optional[int] = None
    expected_n_frames: Optional[int] = None
    command: list[str] = field(default_factory=list)
    stdin: Optional[str] = None
    warnings: list[CampaignWarning] = field(default_factory=list)

    @property
    def n_frames(self) -> int:
        return len(self.times_ps)

    @property
    def start_time_ps(self) -> Optional[float]:
        return self.times_ps[0] if self.times_ps else None

    @property
    def end_time_ps(self) -> Optional[float]:
        return self.times_ps[-1] if self.times_ps else None

    def same_timeline(self, other: "FrameTimeIndex") -> bool:
        """Exact frame-by-frame equality of times (and steps when both have them).

        This is the only sanctioned way to claim that a derived trajectory
        keeps its source's timeline — timestamps are never inherited by assumption.
        """
        if self.times_ps != other.times_ps:
            return False
        if self.steps and other.steps and self.steps != other.steps:
            return False
        return True

    def to_dict(self, *, include_frames: bool = False) -> dict:
        d = {
            "schema_version": self.schema_version,
            "trajectory_path": self.trajectory_path,
            "fingerprint": self.fingerprint.to_dict() if self.fingerprint else None,
            "backend": self.backend,
            "backend_version": self.backend_version,
            "time_unit_reported": self.time_unit_reported,
            "time_resolution_ps": self.time_resolution_ps,
            "n_frames": self.n_frames,
            "start_time_ps": self.start_time_ps,
            "end_time_ps": self.end_time_ps,
            "has_boxes": self.boxes is not None,
            "state": self.state,
            "strictly_increasing": self.strictly_increasing,
            "non_decreasing": self.non_decreasing,
            "duplicate_groups": self.duplicate_groups,
            "non_monotonic_steps": self.non_monotonic_steps,
            "non_finite_frames": self.non_finite_frames,
            "truncated": self.truncated,
            "truncation_evidence": self.truncation_evidence,
            "backend_reported_frames": self.backend_reported_frames,
            "expected_n_frames": self.expected_n_frames,
            "command": self.command,
            "stdin": self.stdin,
            "warnings": [w.to_dict() for w in self.warnings],
        }
        if include_frames:
            d["times_ps"] = list(self.times_ps)
            d["steps"] = list(self.steps)
            d["boxes"] = [list(b) for b in self.boxes] if self.boxes is not None else None
        return d

    @classmethod
    def from_dict(cls, d: dict, *, times_ps: Optional[list[float]] = None,
                  steps: Optional[list[Optional[int]]] = None,
                  boxes: Optional[list[list[float]]] = None) -> "FrameTimeIndex":
        """Rebuild from ``to_dict()``; per-frame arrays come from ``d`` or the
        explicit arguments (sidecar storage)."""
        fp = d.get("fingerprint")
        return cls(
            trajectory_path=d["trajectory_path"],
            fingerprint=SourceFingerprint.from_dict(fp) if fp else None,
            backend=d["backend"],
            backend_version=d.get("backend_version"),
            times_ps=list(times_ps if times_ps is not None else d.get("times_ps", [])),
            steps=list(steps if steps is not None else d.get("steps", [])),
            boxes=boxes if boxes is not None else d.get("boxes"),
            schema_version=d.get("schema_version", TIME_INDEX_SCHEMA_VERSION),
            time_unit_reported=d.get("time_unit_reported", "ps"),
            time_resolution_ps=d.get("time_resolution_ps"),
            state=d.get("state", ValidationState.VALID),
            strictly_increasing=d.get("strictly_increasing"),
            non_decreasing=d.get("non_decreasing"),
            duplicate_groups=[list(g) for g in d.get("duplicate_groups", [])],
            non_monotonic_steps=[list(p) for p in d.get("non_monotonic_steps", [])],
            non_finite_frames=list(d.get("non_finite_frames", [])),
            truncated=d.get("truncated", False),
            truncation_evidence=d.get("truncation_evidence", ""),
            backend_reported_frames=d.get("backend_reported_frames"),
            expected_n_frames=d.get("expected_n_frames"),
            command=list(d.get("command", [])),
            stdin=d.get("stdin"),
            warnings=[CampaignWarning.from_dict(w) for w in d.get("warnings", [])],
        )


class SegmentRelationKind:
    """How two independently indexed segments relate in simulation time.

    Descriptive only: no precedence is implied and nothing is concatenated.
    """
    SEQUENTIAL = "sequential"            # second starts after first ends (gap_ps >= 0 apart)
    BOUNDARY_SHARED = "boundary_shared"  # second starts exactly at first's last time
    OVERLAP = "overlap"                  # time ranges overlap
    REVERSED = "reversed"                # second ends before first starts
    UNKNOWN = "unknown"                  # at least one segment has no frames


@dataclass
class SegmentRelation:
    first_path: str
    second_path: str
    relation: str                        # SegmentRelationKind.*
    first_start_ps: Optional[float] = None
    first_end_ps: Optional[float] = None
    second_start_ps: Optional[float] = None
    second_end_ps: Optional[float] = None
    gap_ps: Optional[float] = None       # second_start - first_end (negative => overlap)
    overlap_ps: Optional[float] = None
    shared_timestamps: int = 0           # exact timestamps present in both segments

    def to_dict(self) -> dict:
        return {
            "first_path": self.first_path,
            "second_path": self.second_path,
            "relation": self.relation,
            "first_start_ps": self.first_start_ps,
            "first_end_ps": self.first_end_ps,
            "second_start_ps": self.second_start_ps,
            "second_end_ps": self.second_end_ps,
            "gap_ps": self.gap_ps,
            "overlap_ps": self.overlap_ps,
            "shared_timestamps": self.shared_timestamps,
        }


# ═══════════════════════════════════════════════════════════════════════════════
# Scientific result contract — labelled N-dimensional arrays
# ═══════════════════════════════════════════════════════════════════════════════
#
# A result is described by its *axes*, never by a plot type: RMSD is
# ``[time]``, APL ``[time, leaflet]``, RMSF ``[residue]``, pore radius
# ``[time, z]``, a contact map ``[residue_i, residue_j]``.  How to render an
# axis signature is a presentation decision made elsewhere.
#
# Values stay in their existing files (XVG/CSV today, NPY/NPZ later); the
# models only *reference* them via :class:`StorageRef`.

RESULT_CONTRACT_VERSION = "simforge/result-array/v1"


class AxisKind:
    """What an axis coordinate *means* (open vocabulary; these are the core kinds)."""
    TIME = "time"                  # simulation time; values need NOT be unique
    COORDINATE = "coordinate"      # continuous spatial / physical coordinate (z, r, ...)
    INDEX = "index"                # ordinal position (frame, carbon number, bin, ...)
    CATEGORY = "category"          # discrete labels (leaflet, "A:45" residues, ...)

    ALL = (TIME, COORDINATE, INDEX, CATEGORY)


class StorageFormat:
    XVG = "xvg"
    CSV = "csv"
    NPY = "npy"
    NPZ = "npz"

    ALL = (XVG, CSV, NPY, NPZ)
    TABULAR = (XVG, CSV)           # addressed by column
    CONTAINER = (NPZ,)             # addressed by member


class MissingSemantics:
    """What a missing / NaN value means for this array (open vocabulary)."""
    UNSPECIFIED = "unspecified"            # not declared (legacy / unknown)
    NONE_EXPECTED = "none_expected"        # every sample is defined; NaN would be an error
    UNDEFINED_STATE = "undefined_state"    # NaN = quantity scientifically undefined there
                                           # (no pocket, outside the pore, ...)
    COMPUTATION_FAILED = "computation_failed"  # NaN = the calculation failed for that sample


class ObservablePurpose:
    """The scientific coordinate purpose an observable serves.

    Reserved for context-dependent preprocessing decisions (a later phase):
    a transformation is judged per (observable, system context, purpose), never
    as universally safe.  Open vocabulary — ``ObservableSpec.purpose`` is a
    plain string so new purposes need no schema change.
    """
    DISPLAY = "display"                                # a view for looking at, not measuring
    INTRAMOLECULAR_SHAPE = "intramolecular_shape"
    INTER_COMPONENT_GEOMETRY = "inter_component_geometry"
    MEMBRANE_FRAME_PROPERTY = "membrane_frame_property"
    DENSITY_PROFILE = "density_profile"
    DIFFUSION = "diffusion"
    SOLVENT_OCCUPANCY = "solvent_occupancy"

    ALL = (DISPLAY, INTRAMOLECULAR_SHAPE, INTER_COMPONENT_GEOMETRY, MEMBRANE_FRAME_PROPERTY,
           DENSITY_PROFILE, DIFFUSION, SOLVENT_OCCUPANCY)


# ═══════════════════════════════════════════════════════════════════════════════
# Preprocessing decisions (scientific admissibility, Phase 5)
# ═══════════════════════════════════════════════════════════════════════════════

class DecisionClass:
    """Admissibility of one operation for one (purpose, system context).

    VALID              acceptable here without further intent (not "universally safe")
    VALID_WITH_INTENT  meaningful, but interpretation-changing: applied only
                       with explicit (non-auto) intent
    REFUSED            conflicts with the purpose / system context
    UNSUPPORTED        meaningful in principle, cannot be done (missing
                       information, semantics or executor)
    """
    VALID = "valid"
    VALID_WITH_INTENT = "valid_with_intent"
    REFUSED = "refused"
    UNSUPPORTED = "unsupported"


class IntentSource:
    AUTO = "auto"                         # derived by SimForge — never explicit intent
    PROFILE = "profile"                   # a named preset the user chose
    USER_FLAG = "user_flag"
    MANIFEST_RESOLUTION = "manifest_resolution"
    EXPLICIT_API = "explicit_api"

    EXPLICIT = (PROFILE, USER_FLAG, MANIFEST_RESOLUTION, EXPLICIT_API)


@dataclass
class PreprocessingDecision:
    operation: str                        # make_whole | nojump | center | fit | minimum_image | reconstruction
    purpose: str
    classification: str                   # DecisionClass.*
    reason: str
    rule_id: str
    rule_version: str
    intent_source: str = IntentSource.AUTO
    parameters: dict[str, Any] = field(default_factory=dict)      # target / group / mode
    requested: bool = True
    applied: bool = False                 # admitted into the resolved requirements
    system_context: dict[str, Any] = field(default_factory=dict)
    diagnostic_evidence: list[dict] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    identity: str = ""
    #: which relevant detectors ran / found nothing / were not applicable / failed /
    #: were not supplied — "no finding" is evidence, "did not run" is not
    diagnostic_coverage: list[dict] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "operation": self.operation, "purpose": self.purpose,
            "classification": self.classification, "reason": self.reason,
            "rule_id": self.rule_id, "rule_version": self.rule_version,
            "intent_source": self.intent_source, "parameters": self.parameters,
            "requested": self.requested, "applied": self.applied,
            "system_context": self.system_context,
            "diagnostic_evidence": self.diagnostic_evidence,
            "notes": self.notes, "identity": self.identity,
            "diagnostic_coverage": self.diagnostic_coverage,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "PreprocessingDecision":
        return cls(
            operation=d["operation"], purpose=d["purpose"], classification=d["classification"],
            reason=d.get("reason", ""), rule_id=d["rule_id"], rule_version=d["rule_version"],
            intent_source=d.get("intent_source", IntentSource.AUTO),
            parameters=dict(d.get("parameters", {})), requested=d.get("requested", True),
            applied=d.get("applied", False), system_context=dict(d.get("system_context", {})),
            diagnostic_evidence=list(d.get("diagnostic_evidence", [])),
            notes=list(d.get("notes", [])), identity=d.get("identity", ""),
            diagnostic_coverage=list(d.get("diagnostic_coverage", [])),
        )


@dataclass
class StorageRef:
    """Where the numbers live and how to address them inside the file.

    ``column`` addresses a 0-based numeric column of a tabular file (XVG/CSV:
    for XVG, column 0 is the x column); ``member`` addresses an array inside a
    container (NPZ).  ``fingerprint`` identifies the file content when known.
    """
    path: str
    format: str                                   # StorageFormat.*
    column: Optional[int] = None
    member: Optional[str] = None
    fingerprint: Optional[SourceFingerprint] = None
    attrs: dict[str, Any] = field(default_factory=dict)

    def validate(self) -> list[str]:
        problems: list[str] = []
        if not self.path:
            problems.append("storage path is empty")
        if self.format not in StorageFormat.ALL:
            problems.append(f"unknown storage format {self.format!r}")
        if self.format in StorageFormat.TABULAR and self.column is None:
            problems.append(f"{self.format} storage needs a column")
        if self.column is not None and self.column < 0:
            problems.append("column must be >= 0")
        if self.format in StorageFormat.CONTAINER and not self.member:
            problems.append(f"{self.format} storage needs a member name")
        return problems

    def to_dict(self) -> dict:
        return {
            "path": self.path,
            "format": self.format,
            "column": self.column,
            "member": self.member,
            "fingerprint": self.fingerprint.to_dict() if self.fingerprint else None,
            "attrs": self.attrs,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "StorageRef":
        fp = d.get("fingerprint")
        return cls(
            path=d["path"], format=d["format"], column=d.get("column"),
            member=d.get("member"),
            fingerprint=SourceFingerprint.from_dict(fp) if fp else None,
            attrs=dict(d.get("attrs", {})),
        )


class FrameAlignmentMode:
    ONE_ROW_PER_FRAME = "one_row_per_frame"   # sample i <-> frame i of the analysed trajectory
    EXPLICIT = "explicit"                     # frame indices given by ``frame_index_ref``
    UNKNOWN = "unknown"


@dataclass
class FrameAlignment:
    """How the samples of a *time* axis correspond to trajectory frames.

    Time values are a physical coordinate, not an identity: they may repeat
    (float32 XTC time, concatenations).  Frame identity comes from the
    analysed trajectory itself — the per-frame ``FrameTimeIndex`` of the file
    whose content digest is ``trajectory_digest`` — so no timeline is copied
    into the result.  ``step_ref`` optionally points at MD-step evidence.
    """
    mode: str = FrameAlignmentMode.UNKNOWN
    trajectory_digest: Optional[str] = None       # SourceFingerprint.digest of the analysed file
    view_ref: Optional[str] = None                # TrajectoryView.cache_key
    frame_index_ref: Optional[StorageRef] = None  # EXPLICIT mode
    step_ref: Optional[StorageRef] = None

    def to_dict(self) -> dict:
        return {
            "mode": self.mode,
            "trajectory_digest": self.trajectory_digest,
            "view_ref": self.view_ref,
            "frame_index_ref": self.frame_index_ref.to_dict() if self.frame_index_ref else None,
            "step_ref": self.step_ref.to_dict() if self.step_ref else None,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "FrameAlignment":
        fr, sr = d.get("frame_index_ref"), d.get("step_ref")
        return cls(
            mode=d.get("mode", FrameAlignmentMode.UNKNOWN),
            trajectory_digest=d.get("trajectory_digest"),
            view_ref=d.get("view_ref"),
            frame_index_ref=StorageRef.from_dict(fr) if fr else None,
            step_ref=StorageRef.from_dict(sr) if sr else None,
        )


@dataclass
class Axis:
    """One labelled dimension of a result.

    Coordinates are either inline (``values`` — small, e.g. category labels)
    or referenced (``values_ref``), or absent when an axis is only declared.
    Values are kept exactly as given: never sorted, never de-duplicated.
    ``frame`` names the reference frame of a coordinate axis
    (e.g. ``"box"``, ``"membrane_normal"``, ``"pore_axis:<annotation-id>"``);
    it is an opaque string here.
    """
    name: str
    kind: str                                     # AxisKind.* (open vocabulary)
    unit: Optional[str] = None
    values: Optional[list[Any]] = None
    values_ref: Optional[StorageRef] = None
    frame: Optional[str] = None
    alignment: Optional[FrameAlignment] = None    # time axes: link to trajectory frames
    attrs: dict[str, Any] = field(default_factory=dict)

    def validate(self) -> list[str]:
        problems: list[str] = []
        if not self.name:
            problems.append("axis name is empty")
        if not self.kind:
            problems.append(f"axis {self.name!r} has no kind")
        if self.values is not None and self.values_ref is not None:
            problems.append(f"axis {self.name!r} has both inline values and values_ref")
        if self.values_ref is not None:
            problems += [f"axis {self.name!r}: {p}" for p in self.values_ref.validate()]
        if self.alignment is not None and self.kind != AxisKind.TIME:
            problems.append(f"axis {self.name!r}: frame alignment only applies to time axes")
        return problems

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "kind": self.kind,
            "unit": self.unit,
            "values": list(self.values) if self.values is not None else None,
            "values_ref": self.values_ref.to_dict() if self.values_ref else None,
            "frame": self.frame,
            "alignment": self.alignment.to_dict() if self.alignment else None,
            "attrs": self.attrs,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "Axis":
        vr, al = d.get("values_ref"), d.get("alignment")
        return cls(
            name=d["name"], kind=d["kind"], unit=d.get("unit"),
            values=list(d["values"]) if d.get("values") is not None else None,
            values_ref=StorageRef.from_dict(vr) if vr else None,
            frame=d.get("frame"),
            alignment=FrameAlignment.from_dict(al) if al else None,
            attrs=dict(d.get("attrs", {})),
        )


@dataclass
class ResultArray:
    """One scientific quantity over labelled axes.

    ``storage is None`` means *declared, not realised* — the same model is
    used by ``ObservableSpec.output_schema()`` before execution.
    """
    name: str                                     # unique within one AnalysisResult
    quantity: str                                 # semantic quantity id, e.g. "rmsd"
    unit: Optional[str] = None
    axes: list[Axis] = field(default_factory=list)
    storage: Optional[StorageRef] = None
    view_ref: Optional[str] = None                # TrajectoryView.cache_key of the analysed view
    definition_token: Optional[str] = None
    missing: str = MissingSemantics.UNSPECIFIED
    attrs: dict[str, Any] = field(default_factory=dict)
    contract_version: str = RESULT_CONTRACT_VERSION

    def axis_signature(self) -> tuple[str, ...]:
        """Axis kinds in order — e.g. ``("time",)``, ``("time", "category")``."""
        return tuple(a.kind for a in self.axes)

    def axis(self, name: str) -> Optional[Axis]:
        for a in self.axes:
            if a.name == name:
                return a
        return None

    @property
    def realised(self) -> bool:
        return self.storage is not None

    def validate(self) -> list[str]:
        problems: list[str] = []
        if not self.name:
            problems.append("result name is empty")
        if not self.quantity:
            problems.append(f"result {self.name!r} has no quantity")
        names = [a.name for a in self.axes]
        if len(set(names)) != len(names):
            problems.append(f"result {self.name!r} has duplicate axis names {names}")
        for a in self.axes:
            problems += a.validate()
        if self.storage is not None:
            problems += [f"result {self.name!r}: {p}" for p in self.storage.validate()]
        return problems

    def to_dict(self) -> dict:
        return {
            "contract_version": self.contract_version,
            "name": self.name,
            "quantity": self.quantity,
            "unit": self.unit,
            "axes": [a.to_dict() for a in self.axes],
            "storage": self.storage.to_dict() if self.storage else None,
            "view_ref": self.view_ref,
            "definition_token": self.definition_token,
            "missing": self.missing,
            "attrs": self.attrs,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "ResultArray":
        st = d.get("storage")
        return cls(
            name=d["name"], quantity=d["quantity"], unit=d.get("unit"),
            axes=[Axis.from_dict(a) for a in d.get("axes", [])],
            storage=StorageRef.from_dict(st) if st else None,
            view_ref=d.get("view_ref"),
            definition_token=d.get("definition_token"),
            missing=d.get("missing", MissingSemantics.UNSPECIFIED),
            attrs=dict(d.get("attrs", {})),
            contract_version=d.get("contract_version", RESULT_CONTRACT_VERSION),
        )


# ═══════════════════════════════════════════════════════════════════════════════
# Trajectory diagnostics (read-only observation; no remediation)
# ═══════════════════════════════════════════════════════════════════════════════

DIAGNOSTICS_SCHEMA_VERSION = "simforge/trajectory-diagnostics/v1"


class SamplingStrategy:
    METADATA_ONLY = "metadata_only"      # derived from the time index / box records only
    ALL_FRAMES = "all_frames"
    REGULAR_STRIDE = "regular_stride"
    ADAPTIVE = "adaptive"


class EvidenceConfidence:
    """Strength of the evidence behind a diagnostic's *interpretation*."""
    STRONG = "strong"
    MODERATE = "moderate"
    WEAK = "weak"


@dataclass
class Diagnostic:
    """One observed trajectory property or anomaly.

    Kept apart: ``measured`` (quantitative observation), ``interpretation``
    (inference), ``confidence`` / ``uncertainty`` (how sure), and
    ``candidate_remediations`` (possible future operations — descriptive ids,
    never ranked, never applied).  ``severity`` uses :class:`Severity`.
    """
    code: str
    severity: str
    scope: str                                   # "timeline", "box", "component:receptor", ...
    message: str
    frames: Optional[list[int]] = None           # [first, last] frame involved
    times_ps: Optional[list[float]] = None       # [first, last] time involved
    measured: dict[str, Any] = field(default_factory=dict)
    interpretation: str = ""
    confidence: str = ""                         # EvidenceConfidence.* ("" = not an inference)
    uncertainty: list[str] = field(default_factory=list)
    sampling: dict[str, Any] = field(default_factory=dict)   # strategy, frames_examined, total_frames
    detector: dict[str, str] = field(default_factory=dict)   # id, version
    candidate_remediations: list[str] = field(default_factory=list)
    provenance: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "code": self.code, "severity": self.severity, "scope": self.scope,
            "message": self.message, "frames": self.frames, "times_ps": self.times_ps,
            "measured": self.measured, "interpretation": self.interpretation,
            "confidence": self.confidence, "uncertainty": self.uncertainty,
            "sampling": self.sampling, "detector": self.detector,
            "candidate_remediations": self.candidate_remediations,
            "provenance": self.provenance,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "Diagnostic":
        return cls(
            code=d["code"], severity=d["severity"], scope=d.get("scope", ""),
            message=d.get("message", ""), frames=d.get("frames"), times_ps=d.get("times_ps"),
            measured=dict(d.get("measured", {})), interpretation=d.get("interpretation", ""),
            confidence=d.get("confidence", ""), uncertainty=list(d.get("uncertainty", [])),
            sampling=dict(d.get("sampling", {})), detector=dict(d.get("detector", {})),
            candidate_remediations=list(d.get("candidate_remediations", [])),
            provenance=dict(d.get("provenance", {})),
        )


class DetectorRunStatus:
    RAN = "ran"
    NOT_APPLICABLE = "not_applicable"
    DEFERRED = "deferred"                # reserved detector, not implemented yet
    FAILED = "failed"                    # probe/detector error (see reason)


@dataclass
class DetectorRun:
    detector_id: str
    version: str
    level: int
    status: str                          # DetectorRunStatus.*
    reason: str = ""
    identity: Optional[str] = None       # definition token of this detector run
    sampling: dict[str, Any] = field(default_factory=dict)
    n_diagnostics: int = 0

    def to_dict(self) -> dict:
        return {"detector_id": self.detector_id, "version": self.version, "level": self.level,
                "status": self.status, "reason": self.reason, "identity": self.identity,
                "sampling": self.sampling, "n_diagnostics": self.n_diagnostics}

    @classmethod
    def from_dict(cls, d: dict) -> "DetectorRun":
        return cls(detector_id=d["detector_id"], version=d["version"], level=d["level"],
                   status=d["status"], reason=d.get("reason", ""), identity=d.get("identity"),
                   sampling=dict(d.get("sampling", {})), n_diagnostics=d.get("n_diagnostics", 0))


@dataclass
class DiagnosticReport:
    trajectory_path: str
    trajectory_fingerprint: Optional[SourceFingerprint] = None
    view_ref: Optional[str] = None
    n_frames: Optional[int] = None
    params: dict[str, Any] = field(default_factory=dict)
    groups: dict[str, Any] = field(default_factory=dict)       # role -> group evidence
    detectors: list[DetectorRun] = field(default_factory=list)
    diagnostics: list[Diagnostic] = field(default_factory=list)
    probes: list[dict[str, Any]] = field(default_factory=list)
    schema_version: str = DIAGNOSTICS_SCHEMA_VERSION

    def summary(self) -> dict:
        counts = {s: 0 for s in (Severity.INFO, Severity.WARN, Severity.REVIEW, Severity.ERROR)}
        for d in self.diagnostics:
            counts[d.severity] = counts.get(d.severity, 0) + 1
        worst = max((d.severity for d in self.diagnostics),
                    key=lambda s: Severity.ORDER.get(s, 0), default=None)
        return {"n_diagnostics": len(self.diagnostics), "by_severity": counts,
                "worst_severity": worst,
                "detectors": {st: sum(1 for r in self.detectors if r.status == st)
                              for st in (DetectorRunStatus.RAN, DetectorRunStatus.NOT_APPLICABLE,
                                         DetectorRunStatus.DEFERRED, DetectorRunStatus.FAILED)}}

    def to_dict(self) -> dict:
        return {
            "schema_version": self.schema_version,
            "trajectory_path": self.trajectory_path,
            "trajectory_fingerprint": (self.trajectory_fingerprint.to_dict()
                                       if self.trajectory_fingerprint else None),
            "view_ref": self.view_ref, "n_frames": self.n_frames,
            "summary": self.summary(), "params": self.params, "groups": self.groups,
            "detectors": [r.to_dict() for r in self.detectors],
            "diagnostics": [d.to_dict() for d in self.diagnostics],
            "probes": self.probes,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "DiagnosticReport":
        fp = d.get("trajectory_fingerprint")
        return cls(
            trajectory_path=d["trajectory_path"],
            trajectory_fingerprint=SourceFingerprint.from_dict(fp) if fp else None,
            view_ref=d.get("view_ref"), n_frames=d.get("n_frames"),
            params=dict(d.get("params", {})), groups=dict(d.get("groups", {})),
            detectors=[DetectorRun.from_dict(r) for r in d.get("detectors", [])],
            diagnostics=[Diagnostic.from_dict(x) for x in d.get("diagnostics", [])],
            probes=list(d.get("probes", [])),
            schema_version=d.get("schema_version", DIAGNOSTICS_SCHEMA_VERSION),
        )


# ═══════════════════════════════════════════════════════════════════════════════
# Resolved scientific annotations (Phase 6)
# ═══════════════════════════════════════════════════════════════════════════════

class AnnotationState:
    """Resolution state of a persistent annotation.

    Only ACTIVE annotations are usable scientific selections (Ann_* groups,
    observable/policy evidence).  PROPOSED ones are visible but inert.
    """
    ACTIVE = "active"                     # explicit (or accepted) and resolved
    PROPOSED = "proposed"                 # derived, not accepted
    REVIEW_REQUIRED = "review_required"   # conflicting explicit declarations
    UNRESOLVED = "unresolved"             # selection matched nothing / missing referents
    UNSUPPORTED = "unsupported"           # numbering / reference cannot be honoured
    SUPERSEDED = "superseded"             # a higher-precedence declaration won
    REJECTED = "rejected"

    USABLE = (ACTIVE,)


class AnnotationCategory:
    RESIDUE_SET = "residue_set"
    AXIS = "axis"
    REFERENCE = "reference"


@dataclass
class AnnotationRecord:
    annotation_id: str
    category: str                         # AnnotationCategory.*
    kind: str
    origin: str
    state: str                            # AnnotationState.*
    definition: dict[str, Any] = field(default_factory=dict)   # the declared definition
    numbering: Optional[str] = None
    reference: dict[str, Any] = field(default_factory=dict)    # structure/file identity used
    residues: list[str] = field(default_factory=list)          # "chain:resname:resnum"
    atom_ids: list[int] = field(default_factory=list)
    n_atoms: int = 0
    atoms_sha256: Optional[str] = None
    identity: Optional[str] = None
    group_name: Optional[str] = None      # Ann_<id> when ACTIVE and atom-based
    depends_on: list[str] = field(default_factory=list)       # referenced annotation ids
    reasons: list[str] = field(default_factory=list)
    conflicts: list[dict] = field(default_factory=list)
    alternatives: list[dict] = field(default_factory=list)    # competing declarations
    provenance: dict[str, Any] = field(default_factory=dict)
    description: str = ""

    @property
    def usable(self) -> bool:
        return self.state in AnnotationState.USABLE

    def evidence(self) -> dict:
        """Compact, identity-bearing evidence for definitions / policy."""
        return {"id": self.annotation_id, "kind": self.kind, "category": self.category,
                "state": self.state, "identity": self.identity,
                "atoms_sha256": self.atoms_sha256, "n_atoms": self.n_atoms}

    def to_dict(self) -> dict:
        return {
            "annotation_id": self.annotation_id, "category": self.category, "kind": self.kind,
            "origin": self.origin, "state": self.state, "definition": self.definition,
            "numbering": self.numbering, "reference": self.reference,
            "residues": self.residues, "atom_ids": self.atom_ids, "n_atoms": self.n_atoms,
            "atoms_sha256": self.atoms_sha256, "identity": self.identity,
            "group_name": self.group_name, "depends_on": self.depends_on,
            "reasons": self.reasons, "conflicts": self.conflicts,
            "alternatives": self.alternatives, "provenance": self.provenance,
            "description": self.description,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "AnnotationRecord":
        return cls(
            annotation_id=d["annotation_id"], category=d["category"], kind=d["kind"],
            origin=d["origin"], state=d["state"], definition=dict(d.get("definition", {})),
            numbering=d.get("numbering"), reference=dict(d.get("reference", {})),
            residues=list(d.get("residues", [])), atom_ids=list(d.get("atom_ids", [])),
            n_atoms=d.get("n_atoms", 0), atoms_sha256=d.get("atoms_sha256"),
            identity=d.get("identity"), group_name=d.get("group_name"),
            depends_on=list(d.get("depends_on", [])), reasons=list(d.get("reasons", [])),
            conflicts=list(d.get("conflicts", [])), alternatives=list(d.get("alternatives", [])),
            provenance=dict(d.get("provenance", {})), description=d.get("description", ""),
        )


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
    #: resolved persistent scientific annotations (never molecular components)
    annotations: list["AnnotationRecord"] = field(default_factory=list)

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
            "annotations": [a.to_dict() for a in self.annotations],
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
            annotations=[AnnotationRecord.from_dict(a) for a in d.get("annotations", [])],
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
    # ── result contract (additive; empty for analyses not yet migrated) ────
    arrays: list[ResultArray] = field(default_factory=list)
    definition_token: Optional[str] = None
    definition_evidence: dict[str, Any] = field(default_factory=dict)

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
            "arrays": [a.to_dict() for a in self.arrays],
            "definition_token": self.definition_token,
            "definition_evidence": self.definition_evidence,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "AnalysisResult":
        """Tolerant of legacy dicts written before the result contract existed."""
        return cls(
            analysis_id=d["analysis_id"], system_id=d["system_id"], status=d["status"],
            message=d.get("message", ""),
            fit_selection=d.get("fit_selection"),
            measure_selection=d.get("measure_selection"),
            reference_frame=d.get("reference_frame", ""),
            trajectory_view_kind=d.get("trajectory_view_kind", ViewKind.RAW),
            parameters=dict(d.get("parameters", {})),
            output_files=list(d.get("output_files", [])),
            data_summary=dict(d.get("data_summary", {})),
            cache_key=d.get("cache_key", ""),
            cached=d.get("cached", False),
            provenance_path=d.get("provenance_path"),
            warnings=[CampaignWarning.from_dict(w) for w in d.get("warnings", [])],
            arrays=[ResultArray.from_dict(a) for a in d.get("arrays", [])],
            definition_token=d.get("definition_token"),
            definition_evidence=dict(d.get("definition_evidence", {})),
        )


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
