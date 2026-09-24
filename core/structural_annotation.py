# core/structural_annotation.py
"""
Structural Biology Annotation Layer.

Models structural knowledge about a molecular system — not execution parameters,
but what we know about the biology: topology, domains, orientation, and the
evidence behind those claims.

These models live at the top level of the YAML:

    structural_annotation:
      membrane_topology:
        extracellular_regions: ["1-50"]
        intracellular_regions: ["200-250"]
        transmembrane_segments: ["51-75", "90-112"]
      orientation:
        extracellular_side: "+z"
        source: "user_annotation"
        confidence: 0.9

Design decisions:
- Residue ranges support "1-50", "1-20,45-60", "5,10,15" formats
- TransmembraneSegment can be a plain string or a rich object in the YAML
- OrientationAnnotation is separated from MembraneTopologyAnnotation because
  orientation (which face is +Z) is a geometric claim, while topology
  (which residues are EC/IC/TM) is a biological claim
- StructuralAnnotation.is_complete_for_orient() decides AutomationLevel
"""

from __future__ import annotations

import re
from typing import Optional, Literal, Union
from pydantic import BaseModel, field_validator, model_validator


# ─────────────────────────────────────────────────────────────────────────────
# Residue range parsing
# ─────────────────────────────────────────────────────────────────────────────

_RANGE_PART = re.compile(r"^\s*(\d+)(?:\s*-\s*(\d+))?\s*$")


def parse_residue_range(s: str) -> list[tuple[int, int]]:
    """
    Parse a residue range string into (start, end) inclusive tuples.

    Supports:
        "1-50"        → [(1, 50)]
        "1-20,45-60"  → [(1, 20), (45, 60)]
        "5,10,15"     → [(5, 5), (10, 10), (15, 15)]
    """
    result: list[tuple[int, int]] = []
    for part in s.split(","):
        m = _RANGE_PART.match(part)
        if not m:
            raise ValueError(
                f"Formato de rango inválido: '{part.strip()}' "
                f"(esperados: '1-50', '1-20,45-60', o '5,10,15')"
            )
        start = int(m.group(1))
        end = int(m.group(2)) if m.group(2) else start
        if start < 1:
            raise ValueError(f"Número de residuo debe ser ≥ 1, se obtuvo {start}")
        if start > end:
            raise ValueError(f"Rango inválido: {start}-{end} (inicio > fin)")
        result.append((start, end))
    return result


def _validate_range_str(v: str) -> str:
    parse_residue_range(v)
    return v


def residues_in_range(range_str: str) -> set[int]:
    """Return all residue numbers covered by a range string."""
    result: set[int] = set()
    for start, end in parse_residue_range(range_str):
        result.update(range(start, end + 1))
    return result


# ─────────────────────────────────────────────────────────────────────────────
# Transmembrane segment
# ─────────────────────────────────────────────────────────────────────────────

class TransmembraneSegment(BaseModel):
    """A transmembrane helix or segment with optional metadata."""
    residues: str
    label: Optional[str] = None       # "TM1", "TM2", etc.
    helix_type: Optional[str] = None  # "alpha", "310", "pi"

    @field_validator("residues")
    @classmethod
    def validate_residues(cls, v: str) -> str:
        return _validate_range_str(v)

    def residue_set(self) -> set[int]:
        return residues_in_range(self.residues)


# ─────────────────────────────────────────────────────────────────────────────
# Biological domain
# ─────────────────────────────────────────────────────────────────────────────

DomainRole = Literal["extracellular", "intracellular", "transmembrane", "loop", "unknown"]


class BiologicalDomain(BaseModel):
    """A named structural region with a biological role."""
    residues: str
    role: DomainRole
    label: Optional[str] = None

    @field_validator("residues")
    @classmethod
    def validate_residues(cls, v: str) -> str:
        return _validate_range_str(v)

    def residue_set(self) -> set[int]:
        return residues_in_range(self.residues)


# ─────────────────────────────────────────────────────────────────────────────
# Annotation source and evidence
# ─────────────────────────────────────────────────────────────────────────────

AnnotationSource = Literal[
    "user_annotation",
    "opm_database",
    "pdbtm",
    "uniprot",
    "predicted",
    "inferred",
]


class OrientationEvidence(BaseModel):
    """Provenance record for an orientation or topology claim."""
    source: AnnotationSource = "user_annotation"
    confidence: float = 1.0
    reference: Optional[str] = None   # OPM accession, PubMed ID, URL, etc.
    notes: Optional[str] = None

    @field_validator("confidence")
    @classmethod
    def validate_confidence(cls, v: float) -> float:
        if not 0.0 <= v <= 1.0:
            raise ValueError(f"confidence debe estar en [0.0, 1.0], se obtuvo {v}")
        return v


# ─────────────────────────────────────────────────────────────────────────────
# Membrane topology
# ─────────────────────────────────────────────────────────────────────────────

class MembraneTopologyAnnotation(BaseModel):
    """
    Describes the disposition of protein regions relative to the membrane.

    Each entry in extracellular_regions, intracellular_regions, and
    transmembrane_segments can be:
        - a string: "1-50" or "1-20,45-60"
        - an object: {residues: "51-75", label: "TM1", helix_type: "alpha"}
    """
    extracellular_regions: list[str] = []
    intracellular_regions: list[str] = []
    transmembrane_segments: list[Union[str, TransmembraneSegment]] = []

    @field_validator("extracellular_regions", "intracellular_regions", mode="before")
    @classmethod
    def validate_region_list(cls, v: list) -> list[str]:
        validated = []
        for item in v:
            validated.append(_validate_range_str(str(item)))
        return validated

    @field_validator("transmembrane_segments", mode="before")
    @classmethod
    def validate_tm_list(cls, v: list) -> list:
        result = []
        for item in v:
            if isinstance(item, str):
                _validate_range_str(item)
                result.append(item)
            elif isinstance(item, dict):
                result.append(TransmembraneSegment(**item))
            elif isinstance(item, TransmembraneSegment):
                result.append(item)
            else:
                raise ValueError(f"Elemento TM inválido: {item!r}")
        return result

    def _ec_residues(self) -> set[int]:
        result: set[int] = set()
        for r in self.extracellular_regions:
            result |= residues_in_range(r)
        return result

    def _ic_residues(self) -> set[int]:
        result: set[int] = set()
        for r in self.intracellular_regions:
            result |= residues_in_range(r)
        return result

    def _tm_residues(self) -> set[int]:
        result: set[int] = set()
        for seg in self.transmembrane_segments:
            r = seg.residues if isinstance(seg, TransmembraneSegment) else seg
            result |= residues_in_range(r)
        return result

    def overlap_warnings(self) -> list[str]:
        """Return consistency warnings about overlapping regions."""
        warns: list[str] = []
        ec, ic, tm = self._ec_residues(), self._ic_residues(), self._tm_residues()

        def _fmt(s: set[int]) -> str:
            top = sorted(s)[:5]
            suffix = "..." if len(s) > 5 else ""
            return f"{top}{suffix}"

        if overlap := ec & ic:
            warns.append(f"Solapamiento EC↔IC en residuos: {_fmt(overlap)}")
        if overlap := ec & tm:
            warns.append(f"Solapamiento EC↔TM en residuos: {_fmt(overlap)}")
        if overlap := ic & tm:
            warns.append(f"Solapamiento IC↔TM en residuos: {_fmt(overlap)}")
        return warns

    def tm_segment_count(self) -> int:
        return len(self.transmembrane_segments)

    def is_sufficient_for_orient(self) -> bool:
        """True if EC and IC regions are defined (minimum for auto-orientation)."""
        return bool(self.extracellular_regions and self.intracellular_regions)


# ─────────────────────────────────────────────────────────────────────────────
# Orientation annotation
# ─────────────────────────────────────────────────────────────────────────────

ZAxis = Literal["+z", "-z"]


class OrientationAnnotation(BaseModel):
    """
    Declares which biological face maps to which geometric axis (+Z / -Z).

    This is a geometric claim, separate from the topological claim in
    MembraneTopologyAnnotation. The orient_protein step uses this to
    decide rotation angles.
    """
    extracellular_side: ZAxis = "+z"
    intracellular_side: ZAxis = "-z"
    source: AnnotationSource = "user_annotation"
    confidence: float = 0.9
    reference: Optional[str] = None

    @field_validator("confidence")
    @classmethod
    def validate_confidence(cls, v: float) -> float:
        if not 0.0 <= v <= 1.0:
            raise ValueError(f"confidence debe estar en [0.0, 1.0]")
        return v

    @model_validator(mode="after")
    def validate_opposite_sides(self) -> "OrientationAnnotation":
        if self.extracellular_side == self.intracellular_side:
            raise ValueError(
                f"extracellular_side e intracellular_side deben ser opuestos "
                f"(ambos son {self.extracellular_side!r})"
            )
        return self

    @property
    def is_standard(self) -> bool:
        """True when EC=+Z, the conventional GROMACS orientation."""
        return self.extracellular_side == "+z"


# ─────────────────────────────────────────────────────────────────────────────
# Persistent scientific annotations (analysis-side reusable intent)
# ─────────────────────────────────────────────────────────────────────────────
#
# An annotation declares scientific *intent* (a binding site, a gate, the
# subunit that matters, an axis, the fit reference).  It never proves the
# biology and it is never a molecular component.  Residue numbers are only
# meaningful together with their numbering scheme: "residue N" in one
# coordinate space is never assumed equal to "residue N" in another.

#: author_pdb — residue numbers/chain ids as deposited (PDB author numbering)
#: topology   — residue numbers of the prepared GROMACS system (.gro/.tpr), with
#:              chains = SimForge polymer segment ids for chain-less files
NumberingScheme = Literal["author_pdb", "topology"]

AnnotationOrigin = Literal[
    "user_cli", "explicit_api", "user_yaml", "study_manifest",
    "imported_ndx", "build_spec", "derived",
]

#: declared lifecycle state — derived annotations start as "proposed"
AnnotationDeclState = Literal["active", "proposed", "review_required", "rejected"]

ResidueSetKind = Literal[
    "binding_site", "catalytic_residues", "functional_residues", "gate",
    "pore_lining", "selected_subunit", "transmembrane_segment", "domain", "custom",
]

_ANNOTATION_ID = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")


def _validate_annotation_id(v: str) -> str:
    if not _ANNOTATION_ID.match(v or ""):
        raise ValueError(f"annotation id {v!r} must match [A-Za-z][A-Za-z0-9_]* "
                         f"(it becomes the index group name Ann_<id>)")
    return v


class ResidueSelection(BaseModel):
    """Residues (and optionally atoms) in an explicit numbering scheme."""
    numbering: NumberingScheme                 # required — never guessed
    chain: Optional[str] = None
    residues: Optional[str] = None             # "684-692", "10,12,15-20"
    resnames: list[str] = []
    atom_names: list[str] = []
    #: restrict matching to polymer (amino-acid / nucleotide) residues — residue
    #: numbering in prepared systems continues into lipids/solvent
    polymer_only: bool = False

    @field_validator("residues")
    @classmethod
    def validate_residues(cls, v: Optional[str]) -> Optional[str]:
        return _validate_range_str(v) if v is not None else v

    @model_validator(mode="after")
    def validate_not_empty(self) -> "ResidueSelection":
        if not (self.chain or self.residues or self.resnames):
            raise ValueError("a residue selection needs at least a chain, residues or resnames")
        return self

    def canonical(self) -> dict:
        """Order-independent definition used for identity."""
        return {"numbering": self.numbering, "chain": self.chain,
                "residues": sorted(residues_in_range(self.residues)) if self.residues else None,
                "resnames": sorted(set(self.resnames)), "atom_names": sorted(set(self.atom_names)),
                "polymer_only": self.polymer_only}


class IndexGroupRef(BaseModel):
    """An annotation defined by an existing .ndx group (identity = its atoms)."""
    path: str
    group: str


class AnnotationProvenance(BaseModel):
    origin: AnnotationOrigin
    source_id: Optional[str] = None            # e.g. a knowledge-bundle claim id (future)
    evidence_class: Optional[str] = None       # experimental | curated | homology_transferred | computational
    loaded_from: Optional[str] = None          # file the declaration came from
    notes: Optional[str] = None


class ResidueSetAnnotation(BaseModel):
    """A named, persistent residue/atom set with a scientific meaning."""
    id: str
    kind: ResidueSetKind
    selection: Optional[ResidueSelection] = None
    index_group: Optional[IndexGroupRef] = None
    label: Optional[str] = None
    description: str = ""                      # never part of the identity
    state: AnnotationDeclState = "active"
    provenance: AnnotationProvenance = AnnotationProvenance(origin="user_yaml")

    @field_validator("id")
    @classmethod
    def validate_id(cls, v: str) -> str:
        return _validate_annotation_id(v)

    @model_validator(mode="after")
    def validate_source(self) -> "ResidueSetAnnotation":
        if (self.selection is None) == (self.index_group is None):
            raise ValueError(f"annotation {self.id!r}: give exactly one of selection / index_group")
        return self


AxisType = Literal["explicit_vector", "com_to_com", "membrane_normal"]


class AxisDefinition(BaseModel):
    type: AxisType
    vector: Optional[tuple[float, float, float]] = None   # explicit_vector
    frame: Optional[str] = None                            # e.g. "box" for explicit vectors
    from_annotation: Optional[str] = None                  # com_to_com
    to_annotation: Optional[str] = None

    @model_validator(mode="after")
    def validate_fields(self) -> "AxisDefinition":
        if self.type == "explicit_vector":
            if self.vector is None or not any(self.vector):
                raise ValueError("explicit_vector needs a non-zero vector")
            if not self.frame:
                raise ValueError("explicit_vector needs its reference frame (e.g. 'box')")
        if self.type == "com_to_com":
            if not (self.from_annotation and self.to_annotation):
                raise ValueError("com_to_com needs from_annotation and to_annotation")
            if self.from_annotation == self.to_annotation:
                raise ValueError("com_to_com endpoints must differ")
        return self


class AxisAnnotation(BaseModel):
    """A persistent scientific axis *definition* (never computed here)."""
    id: str
    kind: str = "axis"                          # e.g. "pore_axis", "axis"
    definition: AxisDefinition
    description: str = ""
    state: AnnotationDeclState = "active"
    provenance: AnnotationProvenance = AnnotationProvenance(origin="user_yaml")

    @field_validator("id")
    @classmethod
    def validate_id(cls, v: str) -> str:
        return _validate_annotation_id(v)


class ReferenceIntent(BaseModel):
    """Which coordinates are the reference for a preprocessing purpose."""
    id: str
    purpose: Literal["fit", "center", "any"] = "fit"
    source: Literal["topology", "structure", "file"]
    file: Optional[str] = None
    description: str = ""
    state: AnnotationDeclState = "active"
    provenance: AnnotationProvenance = AnnotationProvenance(origin="user_yaml")

    @field_validator("id")
    @classmethod
    def validate_id(cls, v: str) -> str:
        return _validate_annotation_id(v)

    @model_validator(mode="after")
    def validate_file(self) -> "ReferenceIntent":
        if (self.source == "file") != bool(self.file):
            raise ValueError("source 'file' requires file (and only then)")
        return self


# ─────────────────────────────────────────────────────────────────────────────
# Top-level structural annotation
# ─────────────────────────────────────────────────────────────────────────────

class StructuralAnnotation(BaseModel):
    """
    First-class structural biology knowledge about a molecular system.

    This is not environment configuration. It represents what we know about
    the biological structure: membrane topology, domain assignments,
    geometric orientation, and the evidence behind those claims.

    Consumed by:
        - orient_protein builder  (decides AutomationLevel, computes rotation)
        - membrane dimension calculator (box sizing)
        - geometric validator (Phase 4)

    YAML location: top-level key `structural_annotation:`
    """
    membrane_topology: Optional[MembraneTopologyAnnotation] = None
    orientation: Optional[OrientationAnnotation] = None
    domains: list[BiologicalDomain] = []
    evidence: list[OrientationEvidence] = []
    #: persistent analysis-side annotations (Phase 6); ignored by build steps
    residue_sets: list[ResidueSetAnnotation] = []
    axes: list[AxisAnnotation] = []

    def is_complete_for_orient(self) -> bool:
        """
        True when enough data is present to fully automate orientation.
        Requires EC regions, IC regions, and an explicit orientation axis.
        """
        if self.membrane_topology is None:
            return False
        return (
            self.membrane_topology.is_sufficient_for_orient()
            and self.orientation is not None
        )

    def is_partial_for_orient(self) -> bool:
        """
        True when topology is present but orientation axis is not explicit.
        orient_protein can still run in GUIDED mode using topology alone.
        """
        if self.membrane_topology is None:
            return False
        return self.membrane_topology.is_sufficient_for_orient()

    def effective_confidence(self) -> float:
        """
        Effective confidence for orientation decisions.
        Uses orientation.confidence when set, otherwise infers from topology.
        """
        if self.orientation is not None:
            return self.orientation.confidence
        if self.is_partial_for_orient():
            return 0.6  # topology known, axis assumed standard (+Z)
        return 0.0

    def validation_warnings(self) -> list[str]:
        """Structural consistency warnings — not errors, but things to surface."""
        warns: list[str] = []
        if self.membrane_topology:
            warns.extend(self.membrane_topology.overlap_warnings())
            warns.extend(self._topology_order_warnings())
            mt = self.membrane_topology
            has_ec = bool(mt.extracellular_regions)
            has_ic = bool(mt.intracellular_regions)
            if has_ec and not has_ic:
                warns.append(
                    "orient_protein → GUIDED: hay regiones EC pero faltan intracellular_regions; "
                    "no se puede calcular el eje EC→IC. "
                    "Define IC regions para habilitar orientación automática."
                )
            elif has_ic and not has_ec:
                warns.append(
                    "orient_protein → GUIDED: hay regiones IC pero faltan extracellular_regions; "
                    "no se puede calcular el eje EC→IC. "
                    "Define EC regions para habilitar orientación automática."
                )
        if self.membrane_topology and not self.orientation:
            warns.append(
                "structural_annotation tiene topología pero orientation no está definida: "
                "se asumirá extracellular_side=+z (convención GROMACS estándar)"
            )
        return warns

    def _topology_order_warnings(self) -> list[str]:
        """
        Compile-time sanity check: for single-pass TM proteins, TM segment
        residue numbers should fall between EC and IC ranges (numerically).

        This is a heuristic — multi-pass TM proteins intentionally interleave
        EC/IC/TM regions, so the warning is suppressed when tm_segment_count > 1.
        """
        mt = self.membrane_topology
        if mt is None or mt.tm_segment_count() != 1:
            return []  # only check single-pass; multi-pass is always interleaved

        ec_nums = mt._ec_residues()
        ic_nums = mt._ic_residues()
        tm_nums = mt._tm_residues()

        if not (ec_nums and ic_nums and tm_nums):
            return []

        ec_mid = sum(ec_nums) / len(ec_nums)
        ic_mid = sum(ic_nums) / len(ic_nums)
        tm_mid = sum(tm_nums) / len(tm_nums)

        lo, hi = min(ec_mid, ic_mid), max(ec_mid, ic_mid)

        warns = []
        if not (lo <= tm_mid <= hi):
            warns.append(
                f"TM segment COM residue {tm_mid:.0f} no está entre el rango EC ({ec_mid:.0f}) "
                f"e IC ({ic_mid:.0f}). Para proteínas single-pass esto puede indicar residuos incorrectos."
            )
        return warns


# ─────────────────────────────────────────────────────────────────────────────
# Backward compat migration
# ─────────────────────────────────────────────────────────────────────────────

def migrate_from_legacy_orientation(
    extracellular_residues: Optional[str],
    intracellular_residues: Optional[str],
    tm_segments: Optional[str],
) -> StructuralAnnotation:
    """
    Build a StructuralAnnotation from the old environment.membrane.orientation fields.
    Called by the parser when structural_annotation is absent but the legacy
    fields are present, so existing YAML configs keep working.
    """
    ec = [extracellular_residues] if extracellular_residues else []
    ic = [intracellular_residues] if intracellular_residues else []
    tm = [tm_segments] if tm_segments else []

    topology = MembraneTopologyAnnotation(
        extracellular_regions=ec,
        intracellular_regions=ic,
        transmembrane_segments=tm,
    )
    evidence = [OrientationEvidence(
        source="user_annotation",
        confidence=0.9,
        notes="Migrado desde environment.membrane.orientation (formato legacy)",
    )]
    return StructuralAnnotation(membrane_topology=topology, evidence=evidence)


# ─────────────────────────────────────────────────────────────────────────────
# Study-level annotation file (simforge_annotations.yaml)
# ─────────────────────────────────────────────────────────────────────────────

ANNOTATION_SET_SCHEMA = "simforge/annotations/v1"
ANNOTATION_FILENAME = "simforge_annotations.yaml"


class AnnotationSet(BaseModel):
    """Versioned, human-editable persistence of annotations and intent.

    ``accept`` lists ids of derived/proposed annotations the user has accepted
    — the only way (besides a manifest ambiguity resolution) a proposal
    becomes active.
    """
    schema_version: str = ANNOTATION_SET_SCHEMA
    structural_annotation: StructuralAnnotation = StructuralAnnotation()
    reference_intents: list[ReferenceIntent] = []
    accept: list[str] = []
    #: numbering of a .pdb reference structure (.gro/.tpr are always "topology")
    reference_numbering: Optional[NumberingScheme] = None

    @field_validator("schema_version")
    @classmethod
    def validate_schema(cls, v: str) -> str:
        if v != ANNOTATION_SET_SCHEMA:
            raise ValueError(f"unsupported annotation schema {v!r} (expected {ANNOTATION_SET_SCHEMA})")
        return v
