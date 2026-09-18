"""
validators/membrane_water/clash.py

Hard steric-clash test: water oxygen vs. nearest protein/lipid heavy atom,
using the real element-specific van der Waals radius sum (no solvent-probe
padding, no fixed cutoff) — independent of the void-graph voxel grid's
resolution, via the same CellList used elsewhere in this package.
"""
from __future__ import annotations

import numpy as np

from validators.membrane_water.constants import CLASH_VDW_SCALE, WATER_OXYGEN_VDW_RADIUS_NM, vdw_radius_nm
from validators.membrane_water.neighbors import CellList
from validators.membrane_water.species import SpeciesData


def build_obstacle_cell_list(species: SpeciesData) -> tuple["CellList | None", np.ndarray]:
    obstacle_idx = np.concatenate([species.protein_idx, species.lipid_idx])
    if len(obstacle_idx) == 0:
        return None, np.array([])
    coords = species.xyz(obstacle_idx)
    radii = np.array([
        vdw_radius_nm(species.gro.atoms[i].atom_name, species.gro.atoms[i].residue_name)
        for i in obstacle_idx.tolist()
    ])
    max_cell = max(0.4, float(radii.max()) + WATER_OXYGEN_VDW_RADIUS_NM)
    cell_list = CellList(coords, cell_size=max_cell)
    return cell_list, radii


def detect_clashing_waters(
    species: SpeciesData,
    oxygen_xyz: np.ndarray,
    cell_list: "CellList | None",
    radii: np.ndarray,
) -> np.ndarray:
    """Boolean array, one entry per row of oxygen_xyz: True if that water
    oxygen has a hard steric clash with any protein/lipid heavy atom, using
    CLASH_VDW_SCALE x the raw VDW radius sum (see constants.py)."""
    if cell_list is None or len(oxygen_xyz) == 0:
        return np.zeros(len(oxygen_xyz), dtype=bool)
    return cell_list.any_clash(
        oxygen_xyz, WATER_OXYGEN_VDW_RADIUS_NM * CLASH_VDW_SCALE, radii * CLASH_VDW_SCALE
    )
