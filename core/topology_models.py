# core/topology_models.py
"""
Data models for the Phase 11 protein-membrane topology workflow.

Architecture:
    generate_protein_topology  →  pdb2gmx on original protein PDB only
    generate_topology          →  assembles embed-time topol.top from components
    assemble_system_topology   →  final topol.top with post-refill molecule counts
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


class TopologyBuildError(RuntimeError):
    """Raised when pdb2gmx would be run on a mixed/embedded system file."""


@dataclass
class TopologyComponent:
    """One molecule type in the assembled system topology."""
    name: str
    kind: str             # "protein", "lipid", "ligand", "solvent", "ion"
    molecule_itp: str     # path to ITP defining this moleculetype
    count: int
    atomtypes_itp: Optional[str] = None  # separate atomtypes ITP (ligands only)


@dataclass
class TopologyAssemblyReport:
    passed: bool
    forcefield: str
    components: list[TopologyComponent] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    gro_source: str = ""
    molecule_counts: dict[str, int] = field(default_factory=dict)


# ── Residue name sets for mixed-system detection ──────────────────────────────

_LIPID_RESNAMES = frozenset({
    "DPP", "DPPC", "POPC", "POPE", "POPG", "POPS", "CHOL",
    "PALM", "OLEO", "PALC", "SM", "CER", "DLPC", "DMPC",
    "DOPC", "PLPC", "LYPC",
})

_SOLVENT_ION_RESNAMES = frozenset({
    "SOL", "HOH", "WAT", "TIP3", "TIP4", "TIP5",
    "NA", "CL", "SOD", "MG", "K", "CA", "ZN",
})


def is_mixed_system_gro(path: "str | Path") -> bool:
    """Return True if the GRO file contains both protein residues and lipid/solvent residues."""
    path = Path(path)
    if not path.exists():
        return False
    try:
        lines = path.read_text().splitlines()
        if len(lines) < 3:
            return False
        n_atoms = int(lines[1].strip())
        has_protein = False
        has_membrane = False
        for ln in lines[2: 2 + n_atoms]:
            if len(ln) < 10:
                continue
            resname = ln[5:10].strip()
            if resname in _LIPID_RESNAMES or resname in _SOLVENT_ION_RESNAMES:
                has_membrane = True
            elif resname:
                has_protein = True
            if has_protein and has_membrane:
                return True
    except Exception:
        pass
    return False
