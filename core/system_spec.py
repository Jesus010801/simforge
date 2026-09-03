"""
Unified declarative system specification for SimForge.

A ``SystemSpec`` describes *what* a molecular system is made of — a protein
(raw PDB or already-parameterized GRO + chain topologies) plus any number of
independently-parameterized non-protein components (ligands, cofactors,
substrates, inhibitors) — together with the force field and water model.

It is the high-level, file-based front end to the already-validated assembly
engine in ``ligand.integrate.assemble_system_multi``. This module ONLY parses,
resolves paths, validates, and (optionally) auto-discovers sibling topology
files. It never assembles anything itself — see ``core.system_build``.

Schema (YAML), either at the document root or nested under a ``system:`` key::

    system:
      name: glp1r_competitive          # optional label, used for output naming
      protein:
        structure: protein_only.gro    # .pdb/.ent  -> raw mode (runs pdb2gmx)
                                       # .gro        -> preparameterized mode
        topology:                      # preparameterized mode only
          - topol_Protein_chain_A.itp
          - topol_Protein_chain_B.itp
        restraints:                    # optional, preparameterized mode only
          - posre_Protein_chain_A.itp
        # topology: auto  /  restraints: auto   -> sibling-file auto-discovery
      components:
        - id: A6
          topology: A6.itp
          coordinates: A6.gro
          role: inhibitor
        - id: COA
          topology: COA_original.itp
          coordinates: COA.gro
          role: cofactor
      forcefield: oplsaa
      water_model: spce
      box:                             # optional, currently metadata only
        distance: 1.2
        type: triclinic

All relative paths are resolved relative to the directory containing the YAML
file (``load_system_spec``) — never the process working directory.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal, Optional, Union

import yaml
from pydantic import BaseModel, Field, model_validator


# Known values — used for *soft* validation (a warning-free happy path), not to
# reject inputs the downstream GROMACS tooling would actually accept.
KNOWN_FORCEFIELDS = {"oplsaa", "charmm36", "amber99sb", "amber99sb-ildn", "amber14sb", "gromos54a7"}
KNOWN_WATER_MODELS = {"spce", "spc", "tip3p", "tip4p", "tip4pew", "tip5p"}

_RAW_PROTEIN_SUFFIXES = {".pdb", ".ent"}
_PREPARAM_PROTEIN_SUFFIXES = {".gro"}

AUTO = "auto"


class SystemSpecError(ValueError):
    """Raised when a system specification is structurally or semantically invalid."""

    def __init__(self, messages: Union[str, list[str]]):
        if isinstance(messages, str):
            messages = [messages]
        self.messages = messages
        super().__init__("\n".join(f"  - {m}" for m in messages))


class BoxSpec(BaseModel):
    distance: float = 1.2
    type: str = "triclinic"


class ProteinSpec(BaseModel):
    structure: str
    topology: Union[list[str], Literal["auto"], None] = None
    restraints: Union[list[str], Literal["auto"], None] = None

    model_config = {"extra": "forbid"}

    @model_validator(mode="before")
    @classmethod
    def _coerce_scalars(cls, data: Any) -> Any:
        if not isinstance(data, dict):
            return data
        out = dict(data)
        for key in ("topology", "restraints"):
            val = out.get(key)
            if isinstance(val, str) and val != AUTO:
                out[key] = [val]
        return out


class ComponentSpec(BaseModel):
    id: str
    topology: str
    coordinates: str
    role: Optional[str] = None

    model_config = {"extra": "forbid"}


class SystemSpec(BaseModel):
    """A fully-parsed, path-resolved system specification.

    After ``load_system_spec`` all ``*_path`` attributes are absolute ``Path``
    objects; the original (as-written) strings are preserved on the nested
    ``ProteinSpec`` / ``ComponentSpec`` models for provenance.
    """

    name: Optional[str] = None
    protein: ProteinSpec
    components: list[ComponentSpec] = Field(default_factory=list)
    forcefield: str = "oplsaa"
    water_model: str = "spce"
    box: Optional[BoxSpec] = None

    model_config = {"extra": "forbid", "arbitrary_types_allowed": True}

    # ── Resolved state (populated by load_system_spec / resolve) ──────────────
    _base_dir: Optional[Path] = None
    _spec_path: Optional[Path] = None
    protein_mode: Literal["raw_pdb", "preparameterized", ""] = ""
    protein_structure_path: Optional[Path] = None
    protein_topology_paths: list[Path] = Field(default_factory=list)
    protein_restraint_paths: list[Path] = Field(default_factory=list)
    component_topology_paths: dict[str, Path] = Field(default_factory=dict)
    component_coordinate_paths: dict[str, Path] = Field(default_factory=dict)
    discovery_notes: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)

    # ── Construction helpers ─────────────────────────────────────────────────

    @staticmethod
    def _unwrap(doc: Any) -> dict:
        if not isinstance(doc, dict):
            raise SystemSpecError("System spec must be a YAML mapping.")
        if "system" in doc and isinstance(doc["system"], dict):
            return doc["system"]
        return doc

    def normalized_dict(self) -> dict:
        """A stable, JSON-serializable snapshot of the spec as declared."""
        return {
            "name": self.name,
            "protein": self.protein.model_dump(),
            "components": [c.model_dump() for c in self.components],
            "forcefield": self.forcefield,
            "water_model": self.water_model,
            "box": self.box.model_dump() if self.box else None,
        }

    # ── Path resolution + validation ─────────────────────────────────────────

    def resolve(self, base_dir: Path, spec_path: Optional[Path] = None) -> "SystemSpec":
        """Resolve every path relative to ``base_dir`` and validate the result.

        Raises ``SystemSpecError`` with an aggregated, actionable message list.
        """
        self._base_dir = Path(base_dir).resolve()
        self._spec_path = Path(spec_path).resolve() if spec_path else None
        errors: list[str] = []

        def _abs(p: str) -> Path:
            pp = Path(p)
            return pp if pp.is_absolute() else (self._base_dir / pp)

        # ── Protein input mode ───────────────────────────────────────────────
        struct = _abs(self.protein.structure)
        self.protein_structure_path = struct
        suffix = struct.suffix.lower()

        if suffix in _RAW_PROTEIN_SUFFIXES:
            self.protein_mode = "raw_pdb"
        elif suffix in _PREPARAM_PROTEIN_SUFFIXES:
            self.protein_mode = "preparameterized"
        else:
            errors.append(
                f"protein.structure '{self.protein.structure}' has unsupported "
                f"extension '{suffix or '(none)'}'. Use .pdb/.ent for a raw "
                f"protein (runs pdb2gmx) or .gro for an already-parameterized one."
            )
            self.protein_mode = ""

        if self.protein_mode == "raw_pdb":
            if self.protein.topology not in (None, []):
                errors.append(
                    "protein.topology must not be set for a raw .pdb protein — "
                    "pdb2gmx generates the topology. Give a .gro structure to "
                    "supply your own chain topologies."
                )
            if self.protein.restraints not in (None, []):
                errors.append(
                    "protein.restraints must not be set for a raw .pdb protein — "
                    "pdb2gmx generates position restraints."
                )
        elif self.protein_mode == "preparameterized":
            topo = self.protein.topology
            if topo is None or topo == []:
                errors.append(
                    "protein.topology is required for a .gro (preparameterized) "
                    "protein. List one topol_Protein_chain_*.itp per chain, or "
                    "use 'topology: auto' for sibling-file auto-discovery."
                )
            elif topo == AUTO:
                self._auto_discover_topology(errors)
            else:
                self.protein_topology_paths = [_abs(p) for p in topo]

            rest = self.protein.restraints
            if rest == AUTO:
                self._auto_discover_restraints()
            elif isinstance(rest, list):
                self.protein_restraint_paths = [_abs(p) for p in rest]

        # ── Existence checks: protein ────────────────────────────────────────
        if not struct.exists():
            errors.append(f"protein.structure not found: {struct}")
        for p in self.protein_topology_paths:
            if not p.exists():
                errors.append(f"protein topology not found: {p}")
        for p in self.protein_restraint_paths:
            if not p.exists():
                errors.append(f"protein restraint not found: {p}")

        # ── Components ───────────────────────────────────────────────────────
        seen_ids: set[str] = set()
        for comp in self.components:
            if comp.id in seen_ids:
                errors.append(f"duplicate component id '{comp.id}'.")
            seen_ids.add(comp.id)
            if not comp.topology:
                errors.append(f"component '{comp.id}': topology is required.")
            if not comp.coordinates:
                errors.append(
                    f"component '{comp.id}': coordinates are required "
                    f"(a parameterized component needs both an .itp and a .gro)."
                )
            t = _abs(comp.topology) if comp.topology else None
            c = _abs(comp.coordinates) if comp.coordinates else None
            if t is not None:
                self.component_topology_paths[comp.id] = t
                if not t.exists():
                    errors.append(f"component '{comp.id}': topology not found: {t}")
            if c is not None:
                self.component_coordinate_paths[comp.id] = c
                if not c.exists():
                    errors.append(f"component '{comp.id}': coordinates not found: {c}")

        if not self.components:
            errors.append(
                "at least one component is required "
                "(a protein-only system does not need `simforge build`)."
            )

        # ── Soft warnings ───────────────────────────────────────────────────
        if self.forcefield not in KNOWN_FORCEFIELDS:
            self.warnings.append(
                f"forcefield '{self.forcefield}' is not in the known set "
                f"{sorted(KNOWN_FORCEFIELDS)} — it will be passed through to GROMACS as-is."
            )
        if self.water_model not in KNOWN_WATER_MODELS:
            self.warnings.append(
                f"water_model '{self.water_model}' is not in the known set "
                f"{sorted(KNOWN_WATER_MODELS)} — it will be passed through to GROMACS as-is."
            )

        if errors:
            raise SystemSpecError(errors)
        return self

    # ── Auto-discovery (Phase 3) ────────────────────────────────────────────

    def _auto_discover_topology(self, errors: list[str]) -> None:
        d = self.protein_structure_path.parent
        chain = sorted(d.glob("topol_Protein_chain_*.itp"))
        chain = [p for p in chain if not p.name.startswith("posre")]
        single = sorted(p for p in d.glob("topol_Protein.itp"))
        generic = sorted(
            p for p in d.glob("topol*.itp")
            if p not in chain and p not in single and not p.name.startswith("posre")
        )

        if chain:
            if single or generic:
                errors.append(
                    "protein.topology: auto is ambiguous — found chain files "
                    f"{[p.name for p in chain]} alongside {[p.name for p in single + generic]}. "
                    "List the intended topology files explicitly."
                )
                return
            self.protein_topology_paths = chain
            self.discovery_notes.append(
                f"auto-discovered {len(chain)} protein chain topology file(s): "
                f"{[p.name for p in chain]}"
            )
        elif single:
            self.protein_topology_paths = single
            self.discovery_notes.append(
                f"auto-discovered protein topology: {single[0].name}"
            )
        elif generic:
            errors.append(
                "protein.topology: auto found no topol_Protein_chain_*.itp / "
                f"topol_Protein.itp, only {[p.name for p in generic]}. "
                "List the intended topology files explicitly."
            )
        else:
            errors.append(
                f"protein.topology: auto found no topol_Protein_chain_*.itp or "
                f"topol_Protein.itp next to {self.protein_structure_path.name} "
                f"(searched {d})."
            )

    def _auto_discover_restraints(self) -> None:
        d = self.protein_structure_path.parent
        posre = sorted(d.glob("posre_Protein_chain_*.itp"))
        if not posre:
            posre = sorted(d.glob("posre_Protein.itp"))
        if not posre:
            posre = sorted(d.glob("posre*.itp"))
        if posre:
            self.protein_restraint_paths = posre
            self.discovery_notes.append(
                f"auto-discovered {len(posre)} position-restraint file(s): "
                f"{[p.name for p in posre]}"
            )
        else:
            self.discovery_notes.append(
                "restraints: auto found no posre*.itp files — continuing without "
                "explicit restraint includes (pdb2gmx-embedded includes, if any, "
                "are still honored)."
            )


def load_system_spec(path: str | Path) -> SystemSpec:
    """Parse, path-resolve and validate a system spec YAML file.

    Relative paths inside the file are resolved relative to the file's own
    directory. Raises ``SystemSpecError`` (aggregated messages) on any problem.
    """
    p = Path(path).resolve()
    if not p.exists():
        raise SystemSpecError(f"system spec file not found: {p}")
    try:
        doc = yaml.safe_load(p.read_text())
    except yaml.YAMLError as exc:
        raise SystemSpecError(f"invalid YAML in {p}: {exc}")

    body = SystemSpec._unwrap(doc)
    try:
        spec = SystemSpec.model_validate(body)
    except Exception as exc:  # pydantic ValidationError -> flat message
        raise SystemSpecError(_flatten_pydantic_error(exc))
    return spec.resolve(base_dir=p.parent, spec_path=p)


def _flatten_pydantic_error(exc: Exception) -> list[str]:
    errs = getattr(exc, "errors", None)
    if not callable(errs):
        return [str(exc)]
    out: list[str] = []
    for e in exc.errors():
        loc = ".".join(str(x) for x in e.get("loc", ()))
        out.append(f"{loc}: {e.get('msg', 'invalid')}")
    return out or [str(exc)]
