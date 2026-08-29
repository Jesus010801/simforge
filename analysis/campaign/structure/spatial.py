"""Coordinate-level signals for component inference (membrane embedding).

Operates directly on the reference structure (``.gro`` or ``.pdb``).  Used as
*evidence*, never as a sole determinant.  All functions return ``None`` when the
signal cannot be computed rather than guessing.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from analysis.campaign.structure.residue_classes import (
    classify_resname, is_membrane_resname,
)

_LIPID_HEADGROUP_ATOMS = {"P", "P1", "P8", "PO4", "N", "N4", "NC3"}
_BACKBONE_ATOMS = {"CA", "C1'", "P"}   # CA for protein; C1'/P for nucleic


@dataclass
class _Atom:
    chain: str
    resname: str
    atomname: str
    z: float


def _iter_atoms(path: Path):
    p = Path(path)
    suffix = p.suffix.lower()
    text = p.read_text(errors="replace").splitlines()
    if suffix == ".pdb":
        for line in text:
            if line[:6].strip() not in ("ATOM", "HETATM"):
                continue
            try:
                chain = line[21:22].strip() or "_"
                resname = line[17:20].strip()
                atomname = line[12:16].strip()
                z = float(line[46:54])
            except (ValueError, IndexError):
                continue
            yield _Atom(chain, resname, atomname, z)
    elif suffix == ".gro":
        if len(text) < 3:
            return
        try:
            n = int(text[1].strip())
        except ValueError:
            return
        # Reconstruct polymer segment letters the SAME way parsers.parse_gro_structure
        # does (new segment on residue-number reset or polymer/non-polymer flip) so
        # chain ids line up between the structure model and this spatial signal.
        poly_letters = iter("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
        cur_chain = "_"
        prev_resnum: Optional[int] = None
        prev_is_poly: Optional[bool] = None
        for line in text[2:2 + n]:
            if len(line) < 44:
                continue
            try:
                resnum = int(line[0:5])
                resname = line[5:10].strip()
                atomname = line[10:15].strip()
                z = float(line[36:44])
            except (ValueError, IndexError):
                continue
            is_poly = classify_resname(resname) in ("amino_acid", "nucleotide")
            reset = prev_resnum is not None and resnum < prev_resnum - 1
            flip = prev_is_poly is not None and is_poly != prev_is_poly
            if prev_is_poly is None or reset or flip:
                cur_chain = next(poly_letters, "_") if is_poly else "_"
            prev_resnum, prev_is_poly = resnum, is_poly
            yield _Atom(cur_chain if is_poly else "_", resname, atomname, z)


@dataclass
class MembraneEmbedding:
    membrane_present: bool
    membrane_z_min: Optional[float] = None
    membrane_z_max: Optional[float] = None
    # chain_id -> fraction of the membrane slab spanned by that chain's backbone
    chain_span_fraction: dict[str, float] = None
    note: str = ""

    def transmembrane_chains(self, threshold: float = 0.6) -> list[str]:
        if not self.chain_span_fraction:
            return []
        return sorted(c for c, f in self.chain_span_fraction.items() if f >= threshold)


def analyse_membrane_embedding(structure_path: str | Path) -> MembraneEmbedding:
    """Estimate the membrane slab (from lipid headgroups) and how far each
    protein chain's backbone spans it.

    For ``.gro`` (no chain IDs) everything is chain ``"_"`` — the fraction still
    tells you whether *a* protein spans the bilayer.
    """
    lipid_head_z: list[float] = []
    chain_bb_z: dict[str, list[float]] = {}
    saw_lipid = False

    for atom in _iter_atoms(Path(structure_path)):
        cls = classify_resname(atom.resname)
        if cls == "lipid":
            saw_lipid = True
            if atom.atomname in _LIPID_HEADGROUP_ATOMS:
                lipid_head_z.append(atom.z)
        elif cls in ("amino_acid", "nucleotide"):
            if atom.atomname in _BACKBONE_ATOMS:
                chain_bb_z.setdefault(atom.chain, []).append(atom.z)

    if not saw_lipid:
        return MembraneEmbedding(membrane_present=False, chain_span_fraction={},
                                 note="no lipid residues in structure")
    if len(lipid_head_z) < 4:
        return MembraneEmbedding(
            membrane_present=True, chain_span_fraction={},
            note="lipids present but too few headgroup atoms to locate the slab",
        )

    lipid_head_z.sort()
    # robust slab bounds: 10th / 90th percentile of headgroup z
    lo = lipid_head_z[int(0.10 * (len(lipid_head_z) - 1))]
    hi = lipid_head_z[int(0.90 * (len(lipid_head_z) - 1))]
    thickness = hi - lo
    spans: dict[str, float] = {}
    if thickness > 1e-6:
        for chain, zs in chain_bb_z.items():
            if not zs:
                continue
            czlo, czhi = min(zs), max(zs)
            overlap = max(0.0, min(czhi, hi) - max(czlo, lo))
            spans[chain] = round(overlap / thickness, 3)

    return MembraneEmbedding(
        membrane_present=True, membrane_z_min=lo, membrane_z_max=hi,
        chain_span_fraction=spans,
        note=f"membrane slab z=[{lo:.2f},{hi:.2f}] nm from {len(lipid_head_z)} headgroup atoms",
    )
