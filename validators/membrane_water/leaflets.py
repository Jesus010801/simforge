"""
validators/membrane_water/leaflets.py

Local leaflet-surface membrane model. Unlike the legacy
`_compute_membrane_core_z` (a single global (z_bot, z_top) scalar pair for
the whole system), this bins lipid headgroup atoms into an XY grid and fits
a per-cell upper/lower leaflet z-surface, so curvature and local thickness
variation are represented. Any XY cell without enough headgroup atoms on a
given leaflet falls back to the global slab value for that cell — the model
is always fully defined, never NaN, and degrades gracefully to the old
global-slab behavior when headgroup detection is sparse or absent.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from validators.membrane_water.constants import HEADGROUP_ATOM_NAMES, TAIL_ATOM_NAMES
from validators.membrane_water.species import SpeciesData

try:
    from core.membrane_knowledge import lipid_atom_names
except Exception:  # pragma: no cover - core module always present in this repo
    lipid_atom_names = None  # type: ignore[assignment]

_MIN_POINTS_PER_CELL = 3
_CURVATURE_STD_THRESHOLD_NM = 0.15


@dataclass
class MembraneModel:
    cell_size_xy: float
    x_edges: np.ndarray
    y_edges: np.ndarray
    upper_z: np.ndarray            # (nx, ny) per-cell upper leaflet surface
    lower_z: np.ndarray            # (nx, ny) per-cell lower leaflet surface
    is_fallback_cell: np.ndarray   # (nx, ny) bool, True where global slab was used
    global_upper_z: float
    global_lower_z: float
    is_curved: bool
    n_headgroup_atoms: int
    warnings: list[str]

    @property
    def global_core_thickness_nm(self) -> float:
        return max(self.global_upper_z - self.global_lower_z, 0.0)

    def _cell_indices(self, xy: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        ix = np.clip(np.searchsorted(self.x_edges, xy[:, 0], side="right") - 1, 0, self.upper_z.shape[0] - 1)
        iy = np.clip(np.searchsorted(self.y_edges, xy[:, 1], side="right") - 1, 0, self.upper_z.shape[1] - 1)
        return ix, iy

    def core_bounds_batch(self, xy: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Vectorized per-point (lower_z, upper_z) using the local leaflet
        surface (falling back to the global slab per-cell as built)."""
        if len(xy) == 0:
            return np.array([]), np.array([])
        ix, iy = self._cell_indices(xy)
        return self.lower_z[ix, iy], self.upper_z[ix, iy]

    def summary(self) -> dict:
        valid_upper = self.upper_z[~self.is_fallback_cell]
        valid_lower = self.lower_z[~self.is_fallback_cell]
        return {
            "type": "curved_local_leaflet" if self.is_curved else "flat_slab",
            "upper_leaflet_surface_summary": {
                "mean_z_nm": round(float(self.upper_z.mean()), 4),
                "std_z_nm": round(float(valid_upper.std()) if len(valid_upper) else 0.0, 4),
                "fallback_cell_fraction": round(float(self.is_fallback_cell.mean()), 4),
            },
            "lower_leaflet_surface_summary": {
                "mean_z_nm": round(float(self.lower_z.mean()), 4),
                "std_z_nm": round(float(valid_lower.std()) if len(valid_lower) else 0.0, 4),
            },
            "core_thickness_summary": {
                "global_thickness_nm": round(self.global_core_thickness_nm, 4),
                "local_mean_thickness_nm": round(float((self.upper_z - self.lower_z).mean()), 4),
            },
        }


def _headgroup_tail_names(lipid: "str | None", forcefield: "str | None") -> tuple[frozenset[str], frozenset[str]]:
    """Union the canonical registry names (when the lipid+forcefield pair is
    registered in core.membrane_knowledge) with the broader heuristic set,
    rather than replacing it — the registry adds precision for known lipids
    without narrowing detection for files using a different but still-common
    naming variant (e.g. a generic "P" instead of OPLS-AA DPPC's "O33")."""
    hg_names, tail_names = set(HEADGROUP_ATOM_NAMES), set(TAIL_ATOM_NAMES)
    if lipid and forcefield and lipid_atom_names is not None:
        try:
            names = lipid_atom_names(lipid, forcefield)
            hg_names.add(names.headgroup_ref)
            tail_names.add(names.tail_middle)
        except KeyError:
            pass
    return frozenset(hg_names), frozenset(tail_names)


def build_membrane_model(
    species: SpeciesData,
    lipid: "str | None" = None,
    forcefield: "str | None" = None,
    cell_size_xy: float = 1.5,
) -> MembraneModel:
    warnings: list[str] = []
    hg_names, tail_names = _headgroup_tail_names(lipid, forcefield)

    lipid_atoms = [species.gro.atoms[i] for i in species.lipid_idx.tolist()]
    hg = [a for a in lipid_atoms if a.atom_name.strip() in hg_names]
    tail = [a for a in lipid_atoms if a.atom_name.strip() in tail_names]

    if not hg:
        all_lip_z = [a.z for a in lipid_atoms]
        if not all_lip_z:
            warnings.append("no lipid atoms found; falling back to a fixed default membrane slab (4.0-6.0 nm)")
            global_lower, global_upper = 4.0, 6.0
        else:
            mid = sum(all_lip_z) / len(all_lip_z)
            global_lower, global_upper = mid - 1.0, mid + 1.0
            warnings.append("no recognized headgroup atoms found; using lipid-atom mean +/-1.0 nm as membrane slab")
        n_hg = 0
        hg_xyz = np.empty((0, 3))
    else:
        hg_xyz = np.array([[a.x, a.y, a.z] for a in hg])
        tail_xyz = np.array([[a.x, a.y, a.z] for a in tail]) if tail else np.empty((0, 3))
        z_mid = float(tail_xyz[:, 2].mean()) if len(tail_xyz) else float(hg_xyz[:, 2].mean())
        top_z = hg_xyz[hg_xyz[:, 2] > z_mid, 2]
        bot_z = hg_xyz[hg_xyz[:, 2] <= z_mid, 2]
        global_upper = float(top_z.mean()) if len(top_z) else z_mid + 1.0
        global_lower = float(bot_z.mean()) if len(bot_z) else z_mid - 1.0
        n_hg = len(hg)

    xy_source = species.xyz(np.concatenate([species.protein_idx, species.lipid_idx])) \
        if (len(species.protein_idx) or len(species.lipid_idx)) else species.coords
    x_min, y_min = xy_source[:, 0].min(), xy_source[:, 1].min()
    x_max, y_max = xy_source[:, 0].max(), xy_source[:, 1].max()
    nx = max(1, int(np.ceil((x_max - x_min) / cell_size_xy)))
    ny = max(1, int(np.ceil((y_max - y_min) / cell_size_xy)))
    x_edges = x_min + np.arange(nx + 1) * cell_size_xy
    y_edges = y_min + np.arange(ny + 1) * cell_size_xy

    upper_z = np.full((nx, ny), global_upper, dtype=np.float64)
    lower_z = np.full((nx, ny), global_lower, dtype=np.float64)
    is_fallback = np.ones((nx, ny), dtype=bool)

    if n_hg:
        z_mid = float(tail_xyz[:, 2].mean()) if len(tail_xyz) else float(hg_xyz[:, 2].mean())
        ix_all = np.clip(np.searchsorted(x_edges, hg_xyz[:, 0], side="right") - 1, 0, nx - 1)
        iy_all = np.clip(np.searchsorted(y_edges, hg_xyz[:, 1], side="right") - 1, 0, ny - 1)
        is_top = hg_xyz[:, 2] > z_mid
        for cell_ix in range(nx):
            for cell_iy in range(ny):
                sel_cell = (ix_all == cell_ix) & (iy_all == cell_iy)
                top_pts = hg_xyz[sel_cell & is_top, 2]
                bot_pts = hg_xyz[sel_cell & ~is_top, 2]
                if len(top_pts) >= _MIN_POINTS_PER_CELL:
                    upper_z[cell_ix, cell_iy] = float(top_pts.mean())
                    is_fallback[cell_ix, cell_iy] = False
                if len(bot_pts) >= _MIN_POINTS_PER_CELL:
                    lower_z[cell_ix, cell_iy] = float(bot_pts.mean())
                    is_fallback[cell_ix, cell_iy] = False

    fallback_fraction = float(is_fallback.mean())
    if fallback_fraction > 0.5:
        warnings.append(
            f"{fallback_fraction:.0%} of membrane-model XY cells lacked enough headgroup "
            "atoms for a local leaflet surface; using the global slab there"
        )

    valid_upper = upper_z[~is_fallback]
    is_curved = bool(len(valid_upper) >= 4 and float(valid_upper.std()) > _CURVATURE_STD_THRESHOLD_NM)

    return MembraneModel(
        cell_size_xy=cell_size_xy,
        x_edges=x_edges,
        y_edges=y_edges,
        upper_z=upper_z,
        lower_z=lower_z,
        is_fallback_cell=is_fallback,
        global_upper_z=global_upper,
        global_lower_z=global_lower,
        is_curved=is_curved,
        n_headgroup_atoms=n_hg,
        warnings=warnings,
    )
