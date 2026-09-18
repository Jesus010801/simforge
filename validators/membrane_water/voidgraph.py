"""
validators/membrane_water/voidgraph.py

Solvent-accessible void voxelization and connectivity analysis.

Two things happen here, both generalizing the legacy `_flood_fill`/`_dilate_step`
primitives (validators/membrane_protein_spatial_classifier.py) rather than
replacing them with something unrelated:

1. Height classification: instead of one global (z_bot, z_top) slab, each
   voxel is classified TOP_BULK / BOT_BULK / CORE_BAND using the *local*
   leaflet surfaces from leaflets.MembraneModel (evaluated per XY cell, then
   broadcast along Z), so curvature and local thickness are respected.
2. Bulk reachability (top_reached / bottom_reached) is the same 6-connected
   dilation flood-fill as before, but now against a `blocked` grid that
   includes lipid atoms (see excluded_volume.py) — the direct fix for water
   reaching "pore" classification through a lateral lipid packing defect.
3. Connected-component labeling of the CORE_BAND interior (the only region
   that matters for pore/vestibule/cavity/defect classification — bulk
   volume is excluded from labeling entirely, keeping this cheap regardless
   of overall box/solvent size) uses integer-min label propagation: the same
   6-connected neighbor-shift step as `_dilate_step`, generalized from
   boolean OR to integer minimum. This converges to one label per connected
   component without scipy.ndimage.label.
"""
from __future__ import annotations

import numpy as np

from validators.membrane_water.leaflets import MembraneModel

HEIGHT_TOP_BULK = 0
HEIGHT_BOT_BULK = 1
HEIGHT_CORE_BAND = 2


def compute_height_classes(
    membrane_model: MembraneModel,
    origin: tuple[float, float, float],
    spacing: float,
    shape: tuple[int, int, int],
    margin_nm: float = 0.3,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (top_bulk_mask, bot_bulk_mask, core_band_mask), each shape."""
    ox, oy, oz = origin
    Nx, Ny, Nz = shape
    xs = ox + (np.arange(Nx) + 0.5) * spacing
    ys = oy + (np.arange(Ny) + 0.5) * spacing
    zs = oz + (np.arange(Nz) + 0.5) * spacing

    xy_grid = np.stack(np.meshgrid(xs, ys, indexing="ij"), axis=-1).reshape(-1, 2)
    lower_2d, upper_2d = membrane_model.core_bounds_batch(xy_grid)
    lower_2d = lower_2d.reshape(Nx, Ny)
    upper_2d = upper_2d.reshape(Nx, Ny)

    z3 = zs[None, None, :]
    upper3 = upper_2d[:, :, None]
    lower3 = lower_2d[:, :, None]

    top_bulk = z3 > (upper3 + margin_nm)
    bot_bulk = z3 < (lower3 - margin_nm)
    core_band = ~top_bulk & ~bot_bulk
    return top_bulk, bot_bulk, core_band


def pbc_shift(arr: np.ndarray, axis: int, step: int, fill=0):
    """Shift `arr` by one voxel along `axis`. X and Y (axes 0, 1) wrap
    around (np.roll) — the simulation box is periodic in X/Y, and a grid
    that treats its X/Y edges as hard walls creates a false "escape route"
    where bulk water above and below the membrane can connect directly
    around the edge of a lipid patch that (like nearly all real periodic
    membrane systems) is only ever built as one periodic image, not a
    finite island. Z (axis 2) never wraps — that axis is what actually
    distinguishes "above" from "below" the membrane, so it stays a hard
    boundary filled with `fill`."""
    if axis in (0, 1):
        return np.roll(arr, step, axis=axis)
    out = np.full_like(arr, fill)
    src = [slice(None)] * 3
    dst = [slice(None)] * 3
    if step == 1:
        src[axis] = slice(0, -1)
        dst[axis] = slice(1, None)
    else:
        src[axis] = slice(1, None)
        dst[axis] = slice(0, -1)
    out[tuple(dst)] = arr[tuple(src)]
    return out


def _dilate_step_bool(r: np.ndarray) -> np.ndarray:
    e = r.copy()
    for axis in range(3):
        e |= pbc_shift(r, axis, 1, fill=False)
        e |= pbc_shift(r, axis, -1, fill=False)
    return e


def flood_fill(blocked: np.ndarray, seeds: np.ndarray) -> np.ndarray:
    """Iterative 6-connected dilation flood-fill from a seed mask, blocked by
    excluded volume. Identical algorithm to the legacy `_flood_fill`."""
    valid = ~blocked
    reached = seeds & valid
    while True:
        expanded = _dilate_step_bool(reached) & valid
        if int(expanded.sum()) == int(reached.sum()):
            break
        reached = expanded
    return reached


def label_components(domain_mask: np.ndarray) -> np.ndarray:
    """Connected-component labels (6-connectivity) restricted to
    `domain_mask`. Returns an int64 array of the same shape: -1 outside the
    domain, otherwise a label shared by every voxel in the same component
    (the label value is the minimum flat index in that component — its
    numeric value carries no meaning beyond identity)."""
    if not domain_mask.any():
        return np.full(domain_mask.shape, -1, dtype=np.int64)

    flat_idx = np.arange(domain_mask.size, dtype=np.int64).reshape(domain_mask.shape)
    big = np.int64(domain_mask.size)
    labels = np.where(domain_mask, flat_idx, big)

    while True:
        candidate = labels
        for axis in range(3):
            candidate = np.minimum(candidate, pbc_shift(labels, axis, 1, fill=big))
            candidate = np.minimum(candidate, pbc_shift(labels, axis, -1, fill=big))
        candidate = np.where(domain_mask, candidate, big)
        if np.array_equal(candidate, labels):
            break
        labels = candidate

    return np.where(domain_mask, labels, -1)


def distance_to_blocked(blocked: np.ndarray, spacing: float) -> np.ndarray:
    """Discrete BFS distance (in nm) from every void voxel to the nearest
    blocked voxel, 6-connectivity. Not a true Euclidean distance transform,
    but a monotonic proxy for local solvent-accessible radius, used only as
    a secondary signal (bottleneck plausibility) alongside wall-composition
    connectivity — never the primary pore/defect discriminator."""
    dist_steps = np.where(blocked, 0, -1)
    frontier = blocked.copy()
    step = 0
    undetermined = dist_steps < 0
    while undetermined.any():
        step += 1
        frontier = _dilate_step_bool(frontier)
        newly_reached = frontier & undetermined
        if not newly_reached.any():
            break
        dist_steps[newly_reached] = step
        undetermined = dist_steps < 0
    dist_steps[dist_steps < 0] = step + 1
    return dist_steps.astype(np.float64) * spacing
