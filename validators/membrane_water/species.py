"""
validators/membrane_water/species.py

Loads a .gro file (via utils.gro_parser) and splits its atoms into the
species groups the rest of the engine operates on: protein, lipid, water,
ion/other-solvent, and everything else (ligands etc.).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from utils.gro_parser import GroAtom, GroFile, parse_gro
from validators.membrane_water.constants import (
    DEFAULT_LIPID_RESNAMES,
    NON_PROTEIN_SOLVENT_LIKE,
    SOLVENT_RESNAMES,
)


@dataclass
class WaterMolecule:
    resid: int
    atom_indices: list[int]   # indices into SpeciesData.gro.atoms
    oxygen_index: int         # index into SpeciesData.gro.atoms


@dataclass
class SpeciesData:
    gro: GroFile
    lipid_resnames: frozenset[str]
    protein_idx: np.ndarray          # int array, indices of protein heavy atoms
    lipid_idx: np.ndarray            # int array, indices of all lipid atoms
    lipid_resid_to_idx: dict[int, list[int]]
    waters: list[WaterMolecule]
    other_idx: np.ndarray            # indices not protein/lipid/water (ligands, ions, ...)
    coords: np.ndarray               # (N, 3) float array, all atoms, nm

    def xyz(self, idx: "np.ndarray | list[int]") -> np.ndarray:
        return self.coords[np.asarray(idx, dtype=np.int64)]


def load_species(
    gro_path: "Path | str",
    lipid_resnames: "frozenset[str] | None" = None,
) -> SpeciesData:
    gro = parse_gro(gro_path)
    lip_set = frozenset(lipid_resnames) if lipid_resnames else DEFAULT_LIPID_RESNAMES

    n = len(gro.atoms)
    coords = np.empty((n, 3), dtype=np.float64)
    resnames = [""] * n
    atomnames = [""] * n
    resids = np.empty(n, dtype=np.int64)
    for i, atom in enumerate(gro.atoms):
        coords[i, 0] = atom.x
        coords[i, 1] = atom.y
        coords[i, 2] = atom.z
        resnames[i] = atom.residue_name
        atomnames[i] = atom.atom_name
        resids[i] = atom.residue_number

    is_lipid = np.array([r in lip_set for r in resnames], dtype=bool)
    is_solvent_like = np.array([r in NON_PROTEIN_SOLVENT_LIKE for r in resnames], dtype=bool)
    is_water = np.array([r in SOLVENT_RESNAMES for r in resnames], dtype=bool)
    is_hydrogen = np.array([a.strip().upper().lstrip("0123456789").startswith("H") for a in atomnames], dtype=bool)

    is_protein = (~is_lipid) & (~is_solvent_like) & (~is_hydrogen)
    protein_idx = np.nonzero(is_protein)[0]
    lipid_idx = np.nonzero(is_lipid)[0]
    other_idx = np.nonzero((~is_lipid) & (~is_water) & (~is_protein) & (~is_hydrogen))[0]

    lipid_resid_to_idx: dict[int, list[int]] = {}
    for i in lipid_idx:
        lipid_resid_to_idx.setdefault(int(resids[i]), []).append(int(i))

    waters = _group_waters(gro, is_water, resids, atomnames)

    return SpeciesData(
        gro=gro,
        lipid_resnames=lip_set,
        protein_idx=protein_idx,
        lipid_idx=lipid_idx,
        lipid_resid_to_idx=lipid_resid_to_idx,
        waters=waters,
        other_idx=other_idx,
        coords=coords,
    )


def _group_waters(
    gro: GroFile,
    is_water: np.ndarray,
    resids: np.ndarray,
    atomnames: list[str],
) -> list[WaterMolecule]:
    """Group water atoms into whole molecules using contiguous-run + resid,
    matching the wrap-safe convention already used by pore_hydration.py
    (resid alone is not reliable across the .gro 5-digit wraparound)."""
    waters: list[WaterMolecule] = []
    current: "WaterMolecule | None" = None
    last_index = -2
    for i in range(len(gro.atoms)):
        if not is_water[i]:
            current = None
            continue
        name = atomnames[i].strip().upper()
        is_oxygen = name.startswith("O")
        if current is None or int(resids[i]) != current.resid or is_oxygen or last_index != i - 1:
            current = WaterMolecule(resid=int(resids[i]), atom_indices=[], oxygen_index=-1)
            waters.append(current)
        current.atom_indices.append(i)
        if is_oxygen:
            current.oxygen_index = i
        last_index = i
    return [w for w in waters if w.oxygen_index >= 0]
