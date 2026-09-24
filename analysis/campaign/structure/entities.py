"""Non-polymer molecular entities and explicit atom selections.

Structural evidence only (residue identity, chain, residue number, atom and
heavy-atom counts) — offline, no databases, no trajectory analysis, no
binding-site reasoning.

A *ligand candidate* is a residue whose name is not in any class SimForge
already understands (amino acid, nucleotide, water, ion, lipid, cofactor) and
which is not embedded inside a polymer sequence (a modified residue flanked by
polymer residues with consecutive numbering is part of that polymer).  Being a
candidate is not being a ligand: role resolution happens in the component
detector, conservatively.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from analysis.campaign.structure.residue_classes import classify_resname
from analysis.campaign.structure.spatial import _iter_atoms

#: A candidate needs at least this many heavy atoms to be a *plausible* ligand.
#: Smaller unknown molecules (glycerol 6, sulfate/phosphate 5, acetate or
#: ethylene glycol 4 …) are typical additives; they are reported, never resolved
#: automatically.  A conservative screening value, not a chemical definition.
LIGAND_MIN_HEAVY_ATOMS = 7

_POLYMER_CLASSES = ("amino_acid", "nucleotide")


@dataclass
class NonPolymerEntity:
    chain: str
    resseq: int
    resname: str
    atom_ids: list[int] = field(default_factory=list)   # 1-based, file order
    heavy_atoms: int = 0

    def to_dict(self) -> dict:
        return {"chain": self.chain, "resid": self.resseq, "resname": self.resname,
                "n_atoms": len(self.atom_ids), "heavy_atoms": self.heavy_atoms,
                "first_atom": self.atom_ids[0] if self.atom_ids else None}


@dataclass
class CandidateScan:
    candidates: list[NonPolymerEntity] = field(default_factory=list)
    excluded_classes: dict[str, int] = field(default_factory=dict)    # residues per known class
    embedded_nonstandard: list[dict] = field(default_factory=list)    # modified polymer residues


def _is_heavy(atom_name: str) -> bool:
    name = atom_name.lstrip("0123456789")
    return bool(name) and name[0].upper() != "H"


def _residues(structure_path: Path):
    """Residues in file order: (key, class, [atom ids], heavy count)."""
    order: list[tuple] = []
    data: dict[tuple, dict] = {}
    for a in _iter_atoms(Path(structure_path), first_model_only=True):
        key = (a.chain, a.resnum, a.resname)
        d = data.get(key)
        if d is None:
            d = data[key] = {"cls": classify_resname(a.resname), "atoms": [], "heavy": 0}
            order.append(key)
        d["atoms"].append(a.index)
        d["heavy"] += int(_is_heavy(a.atomname))
    return [(k, data[k]["cls"], data[k]["atoms"], data[k]["heavy"]) for k in order]


def scan_ligand_candidates(structure_path: str | Path) -> CandidateScan:
    res = _residues(Path(structure_path))
    scan = CandidateScan()
    excluded: Counter = Counter()
    for i, (key, cls, atoms, heavy) in enumerate(res):
        chain, resseq, resname = key
        if cls != "other":
            excluded[cls] += 1
            continue
        prev_ = res[i - 1] if i > 0 else None
        next_ = res[i + 1] if i + 1 < len(res) else None
        if (prev_ and next_ and prev_[1] in _POLYMER_CLASSES and next_[1] in _POLYMER_CLASSES
                and prev_[0][1] == resseq - 1 and next_[0][1] == resseq + 1):
            scan.embedded_nonstandard.append({"chain": chain, "resid": resseq, "resname": resname})
            continue
        scan.candidates.append(NonPolymerEntity(chain=chain, resseq=resseq, resname=resname,
                                                atom_ids=list(atoms), heavy_atoms=heavy))
    scan.excluded_classes = dict(sorted(excluded.items()))
    return scan


# ═══════════════════════════════════════════════════════════════════════════════
# Explicit selections (manifest resolution)
# ═══════════════════════════════════════════════════════════════════════════════

SELECTION_KEYS = ("chain", "resname", "resid", "group")


def parse_selection(text: str) -> dict:
    """``"chain:A,resname:LIG,resid:301"`` / ``"group:MyLig"`` / bare ``"A"`` (chain)."""
    text = text.strip()
    if ":" not in text:
        return {"chain": text}
    sel: dict = {}
    for part in text.split(","):
        k, _, v = part.partition(":")
        k, v = k.strip().lower(), v.strip()
        if k not in SELECTION_KEYS or not v:
            raise ValueError(f"unsupported selection term {part!r}; use {SELECTION_KEYS}")
        sel[k] = int(v) if k == "resid" else v
    if "group" in sel and len(sel) > 1:
        raise ValueError("'group' cannot be combined with other selection terms")
    return sel


def select_atoms(structure_path: str | Path, sel: dict,
                 index_path: Optional[str] = None) -> list[int]:
    """1-based atom ids (reference-structure order) matched by ``sel``."""
    if "group" in sel:
        from analysis.campaign.structure.index_groups import resolve_index_group
        if not index_path or not Path(index_path).is_file():
            raise ValueError(f"group {sel['group']!r} requested but no index file is available")
        grp = resolve_index_group(index_path, sel["group"])
        if grp is None:
            raise ValueError(f"index group {sel['group']!r} not found in {index_path}")
        return sorted(set(grp.atom_ids))
    out: list[int] = []
    for a in _iter_atoms(Path(structure_path), first_model_only=True):
        if "chain" in sel and a.chain != sel["chain"]:
            continue
        if "resname" in sel and a.resname != sel["resname"]:
            continue
        if "resid" in sel and a.resnum != sel["resid"]:
            continue
        out.append(a.index)
    return out
