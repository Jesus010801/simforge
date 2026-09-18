"""
validators/membrane_water/neighbors.py

Uniform-grid cell list for nearest-neighbor / min-distance queries over a
numpy coordinate array. Pure numpy (no scipy dependency) but vectorized by
grouping query points into their own occupied grid cells first, so the
Python-level loop runs over the number of *occupied cells* rather than over
every atom pair or every individual query point — this is what keeps clash
checks and excluded-volume construction tractable for ~250k-atom systems
without O(N_water * N_protein_atoms) brute force.

If scipy becomes available in the environment in the future, this class's
public interface (min_distance / query_radius_counts) is the seam where a
scipy.spatial.cKDTree-backed implementation could be swapped in without
touching any caller — deliberately not implemented here since it cannot be
exercised/tested in this environment (scipy is not an installed dependency).
"""
from __future__ import annotations

import numpy as np


class CellList:
    """Bucket-grid spatial index over an (N, 3) coordinate array."""

    def __init__(self, coords: np.ndarray, cell_size: float):
        self.coords = np.asarray(coords, dtype=np.float64)
        self.cell_size = float(cell_size)
        n = len(self.coords)
        if n == 0:
            self._keys = np.zeros((0, 3), dtype=np.int64)
            self._order = np.array([], dtype=np.int64)
            self._bucket_start: dict[tuple[int, int, int], int] = {}
            self._bucket_end: dict[tuple[int, int, int], int] = {}
            return

        self._keys = np.floor(self.coords / self.cell_size).astype(np.int64)
        order = np.lexsort((self._keys[:, 2], self._keys[:, 1], self._keys[:, 0]))
        self._order = order
        sorted_keys = self._keys[order]
        boundaries = np.nonzero(np.any(sorted_keys[1:] != sorted_keys[:-1], axis=1))[0] + 1
        starts = np.concatenate(([0], boundaries))
        ends = np.concatenate((boundaries, [n]))
        self._bucket_start = {}
        self._bucket_end = {}
        for s, e in zip(starts.tolist(), ends.tolist()):
            key = tuple(sorted_keys[s].tolist())
            self._bucket_start[key] = s
            self._bucket_end[key] = e

    def _candidates_near_cell(self, cell_key: tuple[int, int, int], radius_cells: int) -> np.ndarray:
        cx, cy, cz = cell_key
        pieces = []
        for dx in range(-radius_cells, radius_cells + 1):
            for dy in range(-radius_cells, radius_cells + 1):
                for dz in range(-radius_cells, radius_cells + 1):
                    key = (cx + dx, cy + dy, cz + dz)
                    s = self._bucket_start.get(key)
                    if s is None:
                        continue
                    e = self._bucket_end[key]
                    pieces.append(self._order[s:e])
        if not pieces:
            return np.array([], dtype=np.int64)
        return np.concatenate(pieces)

    def min_distance(self, query_coords: np.ndarray, cutoff: float) -> np.ndarray:
        """Return, for each query point, the distance to the nearest indexed
        coordinate that is within `cutoff` of the search cells scanned
        (values >= cutoff, including np.inf, mean "no neighbor within
        cutoff"). Queries are grouped by cell so the Python loop is bounded
        by the number of occupied query cells, not len(query_coords)."""
        query_coords = np.asarray(query_coords, dtype=np.float64)
        k = len(query_coords)
        result = np.full(k, np.inf, dtype=np.float64)
        if k == 0 or len(self.coords) == 0:
            return result

        radius_cells = max(1, int(np.ceil(cutoff / self.cell_size)))
        query_keys = np.floor(query_coords / self.cell_size).astype(np.int64)

        order = np.lexsort((query_keys[:, 2], query_keys[:, 1], query_keys[:, 0]))
        sorted_keys = query_keys[order]
        boundaries = np.nonzero(np.any(sorted_keys[1:] != sorted_keys[:-1], axis=1))[0] + 1
        starts = np.concatenate(([0], boundaries))
        ends = np.concatenate((boundaries, [k]))

        for s, e in zip(starts.tolist(), ends.tolist()):
            group_idx = order[s:e]
            key = tuple(sorted_keys[s].tolist())
            candidates = self._candidates_near_cell(key, radius_cells)
            if len(candidates) == 0:
                continue
            pts = query_coords[group_idx]                # (g, 3)
            cand = self.coords[candidates]                # (c, 3)
            diff = pts[:, None, :] - cand[None, :, :]      # (g, c, 3)
            dist = np.sqrt(np.einsum("gcd,gcd->gc", diff, diff))
            result[group_idx] = dist.min(axis=1)
        return result

    def any_within(self, query_coords: np.ndarray, cutoff: float) -> np.ndarray:
        """Vectorized boolean: is there an indexed point within `cutoff` of
        each query point."""
        return self.min_distance(query_coords, cutoff) < cutoff

    def nearest_index(self, query_coords: np.ndarray, cutoff: float) -> np.ndarray:
        """Index (into self.coords) of the nearest indexed point within
        `cutoff` for each query, or -1 if none found."""
        query_coords = np.asarray(query_coords, dtype=np.float64)
        k = len(query_coords)
        result = np.full(k, -1, dtype=np.int64)
        if k == 0 or len(self.coords) == 0:
            return result

        radius_cells = max(1, int(np.ceil(cutoff / self.cell_size)))
        query_keys = np.floor(query_coords / self.cell_size).astype(np.int64)
        order = np.lexsort((query_keys[:, 2], query_keys[:, 1], query_keys[:, 0]))
        sorted_keys = query_keys[order]
        boundaries = np.nonzero(np.any(sorted_keys[1:] != sorted_keys[:-1], axis=1))[0] + 1
        starts = np.concatenate(([0], boundaries))
        ends = np.concatenate((boundaries, [k]))

        for s, e in zip(starts.tolist(), ends.tolist()):
            group_idx = order[s:e]
            key = tuple(sorted_keys[s].tolist())
            candidates = self._candidates_near_cell(key, radius_cells)
            if len(candidates) == 0:
                continue
            pts = query_coords[group_idx]
            cand = self.coords[candidates]
            diff = pts[:, None, :] - cand[None, :, :]
            dist = np.sqrt(np.einsum("gcd,gcd->gc", diff, diff))
            best = dist.argmin(axis=1)
            best_dist = dist[np.arange(len(group_idx)), best]
            ok = best_dist < cutoff
            result[group_idx[ok]] = candidates[best[ok]]
        return result

    def any_clash(
        self,
        query_coords: np.ndarray,
        query_radius: float,
        atom_radii: np.ndarray,
    ) -> np.ndarray:
        """Boolean per query point: True if any indexed atom i satisfies
        dist(query, atom_i) < atom_radii[i] + query_radius — i.e. a real
        variable-radius steric clash, not a single fixed cutoff."""
        query_coords = np.asarray(query_coords, dtype=np.float64)
        atom_radii = np.asarray(atom_radii, dtype=np.float64)
        k = len(query_coords)
        result = np.zeros(k, dtype=bool)
        if k == 0 or len(self.coords) == 0:
            return result

        max_r = float(atom_radii.max())
        cutoff = max_r + float(query_radius)
        radius_cells = max(1, int(np.ceil(cutoff / self.cell_size)))
        query_keys = np.floor(query_coords / self.cell_size).astype(np.int64)
        order = np.lexsort((query_keys[:, 2], query_keys[:, 1], query_keys[:, 0]))
        sorted_keys = query_keys[order]
        boundaries = np.nonzero(np.any(sorted_keys[1:] != sorted_keys[:-1], axis=1))[0] + 1
        starts = np.concatenate(([0], boundaries))
        ends = np.concatenate((boundaries, [k]))

        for s, e in zip(starts.tolist(), ends.tolist()):
            group_idx = order[s:e]
            key = tuple(sorted_keys[s].tolist())
            candidates = self._candidates_near_cell(key, radius_cells)
            if len(candidates) == 0:
                continue
            pts = query_coords[group_idx]
            cand = self.coords[candidates]
            diff = pts[:, None, :] - cand[None, :, :]
            dist = np.sqrt(np.einsum("gcd,gcd->gc", diff, diff))
            thresh = atom_radii[candidates][None, :] + query_radius
            clash = (dist < thresh).any(axis=1)
            result[group_idx] = clash
        return result
