"""
validators/membrane_water/excluded_volume.py

Vectorized excluded-volume voxel grid for the void graph. Unlike the legacy
`_build_occupancy_grid` (validators/membrane_protein_spatial_classifier.py),
this stamps BOTH protein and lipid atoms as obstacles (the root cause of the
lateral-defect misclassification bug — lipid was never blocking flood-fill
before), uses per-element van der Waals radii instead of one global padding
value, and adds a configurable solvent-probe radius on top.

Also emits a parallel "wall species" grid (0=void, 1=protein, 2=lipid) used
by regions.py to compute protein_wall_fraction / lipid_wall_fraction — the
signal that actually distinguishes a protein-lined pore from a lipid-lined
lateral defect.
"""
from __future__ import annotations

import math

import numpy as np

from validators.membrane_water.constants import vdw_radius_nm
from validators.membrane_water.species import SpeciesData

WALL_NONE = 0
WALL_PROTEIN = 1
WALL_LIPID = 2


def compute_grid_geometry(
    species: SpeciesData,
    grid_spacing_nm: float,
    buf_xy: float = 1.0,
    buf_z: float = 0.5,
) -> tuple[tuple[float, float, float], tuple[int, int, int]]:
    """X/Y span the full periodic simulation box (origin (0, 0)) so the grid
    can wrap correctly (see voidgraph.pbc_shift) — a padded atom bounding
    box in X/Y would treat the box edges as hard walls, creating a false
    top-to-bottom bulk shortcut around what is physically one periodic
    image of an infinite membrane, not a finite island. `buf_xy` is
    accepted for backward compatibility and ignored. Z is not periodic here
    (it is what distinguishes "above" from "below" the membrane), so it
    keeps the padded-atom-bounding-box treatment.
    """
    idx = np.concatenate([species.protein_idx, species.lipid_idx])
    if len(idx) == 0:
        idx = np.arange(len(species.gro.atoms))
    pts = species.xyz(idx)
    box = species.gro.box
    box_x = float(box[0]) if len(box) > 0 and box[0] > 0 else float(pts[:, 0].max() - pts[:, 0].min() + 2 * buf_xy)
    box_y = float(box[1]) if len(box) > 1 and box[1] > 0 else float(pts[:, 1].max() - pts[:, 1].min() + 2 * buf_xy)
    z_min, z_max = float(pts[:, 2].min()), float(pts[:, 2].max())
    origin = (0.0, 0.0, z_min - buf_z)
    z_span = z_max + buf_z - origin[2]
    shape = (
        max(4, int(math.ceil(box_x / grid_spacing_nm))),
        max(4, int(math.ceil(box_y / grid_spacing_nm))),
        max(4, int(math.ceil(z_span / grid_spacing_nm)) + 1),
    )
    return origin, shape  # type: ignore[return-value]


def build_excluded_volume(
    species: SpeciesData,
    grid_spacing_nm: float,
    solvent_probe_radius_nm: float,
    origin: tuple[float, float, float],
    shape: tuple[int, int, int],
) -> tuple[np.ndarray, np.ndarray]:
    """Returns (blocked, wall_species) grids of `shape`.

    blocked: bool, True = solvent-excluded voxel (protein or lipid VDW +
             solvent-probe radius reaches this voxel center).
    wall_species: int8, WALL_PROTEIN/WALL_LIPID for the *closest* atom that
             claims a blocked voxel, WALL_NONE for void voxels.
    """
    ox, oy, oz = origin
    Nx, Ny, Nz = shape
    blocked = np.zeros(shape, dtype=bool)
    wall_species = np.zeros(shape, dtype=np.int8)
    closest_dist = np.full(shape, np.inf, dtype=np.float64)

    def stamp(atom_indices: np.ndarray, species_code: int) -> None:
        for i in atom_indices.tolist():
            atom = species.gro.atoms[i]
            r_eff = vdw_radius_nm(atom.atom_name, atom.residue_name) + solvent_probe_radius_nm
            ax, ay, az = atom.x, atom.y, atom.z
            r_cells = int(math.ceil(r_eff / grid_spacing_nm))
            ix_c = int(round((ax - ox) / grid_spacing_nm))
            iy_c = int(round((ay - oy) / grid_spacing_nm))
            iz_c = int(round((az - oz) / grid_spacing_nm))
            iz_lo, iz_hi = max(0, iz_c - r_cells), min(Nz, iz_c + r_cells + 1)
            if iz_lo >= iz_hi:
                continue

            # X/Y neighborhood is computed in continuous (unwrapped) space —
            # an atom near a box edge has real neighbor voxels on the OTHER
            # side of the periodic boundary — then wrapped (mod Nx/Ny) only
            # for the actual array index, matching pbc_shift's convention
            # (X/Y periodic, Z a hard boundary).
            raw_ix = np.arange(ix_c - r_cells, ix_c + r_cells + 1)
            raw_iy = np.arange(iy_c - r_cells, iy_c + r_cells + 1)
            iz_range = np.arange(iz_lo, iz_hi)

            xs = ox + raw_ix * grid_spacing_nm
            ys = oy + raw_iy * grid_spacing_nm
            zs = oz + iz_range * grid_spacing_nm
            dist2 = (
                (xs[:, None, None] - ax) ** 2
                + (ys[None, :, None] - ay) ** 2
                + (zs[None, None, :] - az) ** 2
            )
            mask = dist2 <= r_eff * r_eff
            if not mask.any():
                continue
            dist = np.sqrt(dist2)

            idx3 = np.ix_(raw_ix % Nx, raw_iy % Ny, iz_range)

            sub_blocked = blocked[idx3]
            blocked[idx3] = sub_blocked | mask

            sub_dist = closest_dist[idx3]
            better = mask & (dist < sub_dist)
            closest_dist[idx3] = np.where(better, dist, sub_dist)

            sub_species = wall_species[idx3]
            wall_species[idx3] = np.where(better, species_code, sub_species)

    stamp(species.protein_idx, WALL_PROTEIN)
    stamp(species.lipid_idx, WALL_LIPID)
    return blocked, wall_species


def voxel_index(
    x: float, y: float, z: float,
    origin: tuple[float, float, float],
    spacing: float,
    shape: tuple[int, int, int],
) -> tuple[int, int, int]:
    ox, oy, oz = origin
    Nx, Ny, Nz = shape
    # X/Y wrap (periodic box, atoms may carry unwrapped coordinates outside
    # [0, box)); Z is a hard boundary (not periodic — see pbc_shift).
    ix = int(round((x - ox) / spacing)) % Nx
    iy = int(round((y - oy) / spacing)) % Ny
    iz = max(0, min(int(round((z - oz) / spacing)), Nz - 1))
    return ix, iy, iz


def voxel_indices_batch(
    xyz: np.ndarray,
    origin: tuple[float, float, float],
    spacing: float,
    shape: tuple[int, int, int],
) -> np.ndarray:
    """Vectorized voxel_index for an (N, 3) array; returns (N, 3) int array."""
    ox, oy, oz = origin
    Nx, Ny, Nz = shape
    ijk = np.round((xyz - np.array([ox, oy, oz])) / spacing).astype(np.int64)
    ijk[:, 0] = np.mod(ijk[:, 0], Nx)
    ijk[:, 1] = np.mod(ijk[:, 1], Ny)
    ijk[:, 2] = np.clip(ijk[:, 2], 0, Nz - 1)
    return ijk
