"""
validators/membrane_water/backbone.py

Backbone/branch decomposition for two-sided (membrane-crossing) void
components.

Why this exists: a single connected void component can contain BOTH the
true transmembrane lumen and a lateral, lipid-facing branch/crevice that
happens to link up with it underground (e.g. a gap between two subunits of
an oligomeric channel). regions.py's whole-component wall-fraction vote
would classify such a merged component as "pore" in aggregate if the lumen
portion dominates the vote count, which lets the lateral branch inherit the
lumen's classification and preserve its water. This module never assumes
one connected component == one biological pore: it finds the principal
upper-to-lower weighted-shortest-path corridor through the component
(favoring protein-lined, wide, solvent-accessible void; penalizing
lipid-facing space), grows a "lumen" region around that corridor only where
the surroundings remain protein-lined, and reclassifies whatever remains in
the component as one or more independent branches — each judged on its own
wall composition and connectivity, never on the classification of the
component it happened to be attached to.

No cylinder, no fixed radius, no fixed XY center, no residue-specific
hard-coding: only local wall composition, void topology and a
distance-to-wall field, all derived solely from protein/lipid atom
positions.
"""
from __future__ import annotations

import heapq
import math
from dataclasses import dataclass

import numpy as np

from validators.membrane_water.constants import (
    DEFAULT_PROTEIN_WALL_THRESHOLD,
    DEFAULT_VESTIBULE_WALL_THRESHOLD,
    REGION_ISOLATED_CAVITY,
    REGION_MEMBRANE_DEFECT,
    REGION_PORE,
    REGION_VESTIBULE,
)
from validators.membrane_water.excluded_volume import WALL_LIPID, WALL_PROTEIN
from validators.membrane_water.voidgraph import label_components, pbc_shift

_NEIGHBOR_OFFSETS = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))

# Path-cost weighting for the backbone search. These are generic, system-
# independent knobs (documented, defaulted) — not a fixed radius/axis for
# any particular channel.
LIPID_PATH_PENALTY = 4.0        # cost multiplier applied for fully lipid-facing voxels
RADIUS_BONUS_CAP_NM = 0.5       # distance-to-wall beyond this stops giving extra path discount
RADIUS_BONUS_MAX_FRACTION = 0.5  # widest voxels get up to this fractional cost discount
CORRIDOR_SLACK_FACTOR = 1.35    # backbone corridor = voxels within this factor of the shortest path cost

# Physical reach (nm) of the smoothed wall-adjacency field, converted to a
# hop count at call time from the actual grid spacing (a fixed hop count
# would give an inconsistent physical reach across different grid
# resolutions). Must comfortably cover a real pore's own radius — a pore
# whose center never sees ANY wall influence in this field would wrongly
# fail the "has_wall_nearby" enclosure requirement below — while staying
# well short of the many-nm scale of genuinely open bulk space.
LOCAL_CHARACTER_SMOOTHING_REACH_NM = 1.2
MIN_COMPONENT_VOXELS_FOR_DECOMPOSITION = 6

# Angular-coverage enclosure test (see _angular_coverage docstring): being
# "near protein" is not the same as being "enclosed by protein" — a void
# voxel skimming the outside of one solid object is near protein on one
# side only, with open space on the other. A real pore's interior is
# surrounded by wall material from most lateral directions, however many
# separate structural elements form it.
ANGULAR_COVERAGE_SEARCH_RADIUS_NM = 1.2
ANGULAR_COVERAGE_N_DIRECTIONS = 16
ANGULAR_COVERAGE_THRESHOLD = 0.65

# Medial-axis-style bound on lumen growth: a void voxel may join the lumen
# from an already-accepted neighbor only if its own distance-to-wall (local
# pore radius) does not exceed the neighbor's radius by more than this
# factor + absolute slack. This is what actually distinguishes "the channel
# gradually widens into a vestibule" (small, repeated radius increases,
# still admitted) from "the channel opens directly into wide lateral lipid
# space" (one large radius jump, rejected) — being merely closer to protein
# than lipid in a smoothed sense is NOT sufficient on its own, since an open
# lateral gap can still be nominally "protein-nearer" while being wide open.
LUMEN_RADIUS_GROWTH_FACTOR = 1.5
LUMEN_RADIUS_GROWTH_ABS_NM = 0.15


@dataclass
class BranchInfo:
    branch_id: int
    voxel_count: int
    connected_to_lumen: bool
    connected_to_upper_bulk: bool
    connected_to_lower_bulk: bool
    protein_wall_fraction: float
    lipid_wall_fraction: float
    wall_composition_defined: bool
    bottleneck_radius_nm: float
    mean_radius_nm: float
    max_radius_nm: float
    classification: str
    confidence: float


@dataclass
class DecompositionResult:
    attempted: bool
    split: bool                      # True if the component was actually split into lumen + branch(es)
    lumen_mask: "np.ndarray | None"  # boolean, same shape as the grid, subset of the component
    lumen_protein_wall_fraction: float
    lumen_classification: str        # REGION_PORE or REGION_MEMBRANE_DEFECT (no viable protein-lined path)
    lumen_confidence: float
    lumen_voxel_count: int
    lumen_bottleneck_radius_nm: float
    lumen_mean_radius_nm: float
    lumen_max_radius_nm: float
    backbone_voxel_count: int
    branch_labels: "np.ndarray | None"  # int, same shape as grid; -1 outside any branch
    branches: dict[int, BranchInfo]


def _shift(arr: np.ndarray, axis: int, step: int, fill=0):
    """X/Y-periodic, Z-bounded shift — see voidgraph.pbc_shift for why the
    periodic wrap matters (a non-periodic grid edge creates a false bulk
    top-to-bottom shortcut around a lipid patch that is, physically, one
    periodic image of an infinite membrane, not a finite island)."""
    return pbc_shift(arr, axis, step, fill=fill)


def _smoothed_wall_votes(blocked: np.ndarray, wall_species: np.ndarray, spacing: float) -> tuple[np.ndarray, np.ndarray]:
    """Spread immediate wall-adjacency votes outward by 6-connected steps (a
    discrete diffusion) reaching LOCAL_CHARACTER_SMOOTHING_REACH_NM,
    giving a soft 'local wall character' field usable at interior void
    voxels several voxels from any wall — not just the ones directly
    touching one. Nearby walls contribute more than distant ones; this is a
    smoothing kernel, not a hard classification."""
    hops = max(1, int(math.ceil(LOCAL_CHARACTER_SMOOTHING_REACH_NM / spacing)))
    protein_votes = (blocked & (wall_species == WALL_PROTEIN)).astype(np.float64)
    lipid_votes = (blocked & (wall_species == WALL_LIPID)).astype(np.float64)
    p, l = protein_votes.copy(), lipid_votes.copy()
    for _ in range(hops):
        p = p + sum(_shift(p, axis, step) for axis in range(3) for step in (1, -1))
        l = l + sum(_shift(l, axis, step) for axis in range(3) for step in (1, -1))
    return p, l


def _angular_coverage(
    blocked: np.ndarray,
    candidate_mask: np.ndarray,
    spacing: float,
    search_radius_nm: float,
    n_directions: int,
) -> np.ndarray:
    """For each True voxel in `candidate_mask`, the fraction of
    `n_directions` evenly-spaced lateral (XY, fixed Z) rays that hit a
    blocked (protein or lipid) voxel within `search_radius_nm`.

    This is the actual test for "enclosed by wall material" vs. "merely
    near one wall": a voxel skimming the outside surface of a single solid
    object has a wall in the direction facing that object, but the opposite
    direction escapes into open space — roughly half its rays miss. A
    voxel near the center of a real pore, surrounded by protein from most
    angles (however many separate helices form it), has most or all of its
    rays hit a wall within a short distance. Vectorized over all candidate
    voxels per (direction, step), not looped per voxel, so cost scales with
    n_directions * n_steps, not with the number of candidate voxels."""
    shape = blocked.shape
    Nx, Ny, Nz = shape
    idx = np.argwhere(candidate_mask)
    coverage_full = np.zeros(shape, dtype=np.float64)
    if len(idx) == 0:
        return coverage_full

    n_steps = max(1, int(round(search_radius_nm / spacing)))
    ix0, iy0, iz0 = idx[:, 0], idx[:, 1], idx[:, 2]
    hit_count = np.zeros(len(idx), dtype=np.int32)

    for d in range(n_directions):
        angle = 2 * math.pi * d / n_directions
        dx, dy = math.cos(angle), math.sin(angle)
        hit_this_dir = np.zeros(len(idx), dtype=bool)
        for step in range(1, n_steps + 1):
            pending = ~hit_this_dir
            if not pending.any():
                break
            ix = np.mod(np.round(ix0[pending] + dx * step).astype(np.int64), Nx)
            iy = np.mod(np.round(iy0[pending] + dy * step).astype(np.int64), Ny)
            iz = iz0[pending]
            vals = blocked[ix, iy, iz]
            pending_indices = np.nonzero(pending)[0]
            hit_this_dir[pending_indices[vals]] = True
        hit_count += hit_this_dir.astype(np.int32)

    coverage_full[ix0, iy0, iz0] = hit_count / n_directions
    return coverage_full


def _local_protein_fraction(local_protein: np.ndarray, local_lipid: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Returns (frac, has_wall_nearby). `frac` defaults to a neutral 0.5
    where no wall is within the smoothing radius at all — fine for the
    path-cost field (an open, wide-open stretch shouldn't be penalized just
    for being far from every wall). `has_wall_nearby` is False there, and
    MUST gate lumen-growth admission separately: "not lipid-dominated"
    (frac >= threshold) is not the same claim as "actually protein-lined" —
    a void voxel with no nearby wall at all is in open bulk-like space, not
    inside a bounded, protein-defined lumen, and must not be admitted just
    because the smoothed fraction defaults to a passing neutral value."""
    total = local_protein + local_lipid
    frac = np.full(total.shape, 0.5, dtype=np.float64)
    has_wall_nearby = total > 0
    frac[has_wall_nearby] = local_protein[has_wall_nearby] / total[has_wall_nearby]
    return frac, has_wall_nearby


def _voxel_cost(local_protein_frac: np.ndarray, dist_to_wall: np.ndarray, spacing: float) -> np.ndarray:
    lipid_penalty = 1.0 + LIPID_PATH_PENALTY * (1.0 - local_protein_frac)
    radius_bonus = np.clip(dist_to_wall / RADIUS_BONUS_CAP_NM, 0.0, 1.0) * RADIUS_BONUS_MAX_FRACTION
    return spacing * lipid_penalty * (1.0 - radius_bonus)


def _multi_source_dijkstra(
    sources: np.ndarray,
    component_set: set,
    cost_flat: np.ndarray,
    shape: tuple[int, int, int],
) -> dict:
    """Weighted shortest path (Dijkstra) from any of `sources` (flat voxel
    indices) through `component_set` (flat indices), edge cost = cost of
    entering the destination voxel. Returns {flat_idx: min_cost}."""
    Nx, Ny, Nz = shape
    dist: dict = {}
    heap = []
    for s in sources.tolist():
        if s in component_set and dist.get(s, math.inf) > 0:
            dist[s] = 0.0
            heapq.heappush(heap, (0.0, s))
    while heap:
        d, u = heapq.heappop(heap)
        if d > dist.get(u, math.inf):
            continue
        uz = u % Nz
        uy = (u // Nz) % Ny
        ux = u // (Ny * Nz)
        for dx, dy, dz in _NEIGHBOR_OFFSETS:
            # X/Y wrap (periodic simulation box); Z stays a hard boundary.
            vx, vy, vz = (ux + dx) % Nx, (uy + dy) % Ny, uz + dz
            if not (0 <= vz < Nz):
                continue
            v = (vx * Ny + vy) * Nz + vz
            if v not in component_set:
                continue
            nd = d + cost_flat[v]
            if nd < dist.get(v, math.inf):
                dist[v] = nd
                heapq.heappush(heap, (nd, v))
    return dist


def _wall_fraction_for_mask(mask: np.ndarray, protein_votes: np.ndarray, lipid_votes: np.ndarray) -> tuple[float, bool]:
    p = float(protein_votes[mask].sum())
    l = float(lipid_votes[mask].sum())
    total = p + l
    if total <= 0:
        return 0.5, False
    return p / total, True


def _no_split_result(whole_frac: float, attempted: bool) -> "DecompositionResult":
    return DecompositionResult(
        attempted=attempted, split=False, lumen_mask=None,
        lumen_protein_wall_fraction=whole_frac, lumen_classification=REGION_MEMBRANE_DEFECT,
        lumen_confidence=0.0, lumen_voxel_count=0,
        lumen_bottleneck_radius_nm=0.0, lumen_mean_radius_nm=0.0, lumen_max_radius_nm=0.0,
        backbone_voxel_count=0, branch_labels=None, branches={},
    )


def decompose_component(
    component_mask: np.ndarray,
    blocked: np.ndarray,
    wall_species: np.ndarray,
    top_bulk_void: np.ndarray,
    bot_bulk_void: np.ndarray,
    dist_to_wall: np.ndarray,
    spacing: float,
    protein_wall_threshold: float = DEFAULT_PROTEIN_WALL_THRESHOLD,
    vestibule_wall_threshold: float = DEFAULT_VESTIBULE_WALL_THRESHOLD,
) -> DecompositionResult:
    voxel_count = int(component_mask.sum())
    shape = component_mask.shape

    # Immediate wall-adjacency votes (same primitive as regions.py), used
    # for the actual pore/vestibule/defect wall-fraction decisions;
    # separately smoothed for the path-cost field only.
    is_protein_blocked = blocked & (wall_species == WALL_PROTEIN)
    is_lipid_blocked = blocked & (wall_species == WALL_LIPID)
    protein_votes = np.zeros(shape, dtype=np.int32)
    lipid_votes = np.zeros(shape, dtype=np.int32)
    for axis in range(3):
        for step in (1, -1):
            protein_votes += _shift(is_protein_blocked, axis, step).astype(np.int32)
            lipid_votes += _shift(is_lipid_blocked, axis, step).astype(np.int32)

    whole_frac, whole_defined = _wall_fraction_for_mask(component_mask, protein_votes, lipid_votes)

    if voxel_count < MIN_COMPONENT_VOXELS_FOR_DECOMPOSITION or not whole_defined or whole_frac <= 0.0:
        # Too small to meaningfully split, or zero protein contact anywhere
        # in the component — no protein-lined path can possibly exist, so
        # there is nothing to find (this is a cheap, principled short
        # circuit: not a radius/size cutoff on the *answer*, just skipping
        # a search that is mathematically guaranteed to fail).
        return _no_split_result(whole_frac, attempted=False)

    angular_coverage = _angular_coverage(
        blocked, component_mask, spacing,
        ANGULAR_COVERAGE_SEARCH_RADIUS_NM, ANGULAR_COVERAGE_N_DIRECTIONS,
    )
    local_p, local_l = _smoothed_wall_votes(blocked, wall_species, spacing)
    local_frac, local_has_wall = _local_protein_fraction(local_p, local_l)
    cost = _voxel_cost(local_frac, dist_to_wall, spacing)
    cost_flat = cost.reshape(-1)

    component_flat = np.flatnonzero(component_mask.reshape(-1))
    component_set = set(component_flat.tolist())

    portal_top = np.flatnonzero((component_mask & _touches(top_bulk_void)).reshape(-1))
    portal_bot = np.flatnonzero((component_mask & _touches(bot_bulk_void)).reshape(-1))

    if len(portal_top) == 0 or len(portal_bot) == 0:
        # Not actually reaching both sides at the voxel level (shouldn't
        # happen for a component already flagged traverses_membrane, but
        # fail safe rather than crash).
        return _no_split_result(whole_frac, attempted=False)

    dist_from_top = _multi_source_dijkstra(portal_top, component_set, cost_flat, shape)
    dist_from_bot = _multi_source_dijkstra(portal_bot, component_set, cost_flat, shape)

    common = set(dist_from_top) & set(dist_from_bot)
    if not common:
        return _no_split_result(whole_frac, attempted=True)

    totals = {v: dist_from_top[v] + dist_from_bot[v] for v in common}
    min_total = min(totals.values())
    threshold = min_total * CORRIDOR_SLACK_FACTOR if min_total > 0 else 0.0
    backbone_flat = np.array([v for v, t in totals.items() if t <= threshold + 1e-12], dtype=np.int64)
    backbone_mask = np.zeros(shape, dtype=bool)
    backbone_mask.reshape(-1)[backbone_flat] = True
    backbone_voxel_count = int(backbone_mask.sum())

    # Grow the lumen from the backbone corridor: a candidate voxel joins
    # only if it is protein-associated (local smoothed character) AND its
    # own local pore radius does not jump abruptly past an already-accepted
    # neighbor's radius (medial-axis growth-rate limit — see
    # LUMEN_RADIUS_GROWTH_FACTOR docstring above). The radius bound is what
    # actually separates the channel's own bore from open lateral space
    # that happens to be nominally protein-nearer in a smoothed sense.
    protein_associated = (
        local_has_wall
        & (local_frac >= vestibule_wall_threshold)
        & (angular_coverage >= ANGULAR_COVERAGE_THRESHOLD)
    )
    lumen_mask = _grow_radius_limited(
        backbone_mask, component_mask & protein_associated, dist_to_wall,
        LUMEN_RADIUS_GROWTH_FACTOR, LUMEN_RADIUS_GROWTH_ABS_NM,
    )

    lumen_frac, lumen_defined = _wall_fraction_for_mask(lumen_mask, protein_votes, lipid_votes)
    if lumen_defined and lumen_frac >= protein_wall_threshold:
        lumen_classification = REGION_PORE
        headroom = max(1e-9, 1.0 - protein_wall_threshold)
        lumen_confidence = float(min(1.0, 0.5 + 0.5 * min(1.0, (lumen_frac - protein_wall_threshold) / headroom)))
    else:
        lumen_classification = REGION_MEMBRANE_DEFECT
        lumen_confidence = 0.0
        lumen_mask = np.zeros(shape, dtype=bool)  # no viable protein-lined lumen found

    if lumen_mask.any():
        lumen_dist = dist_to_wall[lumen_mask]
        lumen_voxel_count = int(lumen_mask.sum())
        lumen_bottleneck = float(lumen_dist.min())
        lumen_mean = float(lumen_dist.mean())
        lumen_max = float(lumen_dist.max())
    else:
        lumen_voxel_count = 0
        lumen_bottleneck = lumen_mean = lumen_max = 0.0

    remainder = component_mask & ~lumen_mask
    if not remainder.any() or not lumen_mask.any():
        # Either the whole component IS the lumen (nothing to split off) or
        # no lumen was found at all (falls back to whole-component handling
        # by the caller) — either way, no split happened.
        return DecompositionResult(
            attempted=True, split=False, lumen_mask=lumen_mask if lumen_mask.any() else None,
            lumen_protein_wall_fraction=lumen_frac if lumen_defined else whole_frac,
            lumen_classification=lumen_classification, lumen_confidence=lumen_confidence,
            lumen_voxel_count=lumen_voxel_count,
            lumen_bottleneck_radius_nm=round(lumen_bottleneck, 4),
            lumen_mean_radius_nm=round(lumen_mean, 4),
            lumen_max_radius_nm=round(lumen_max, 4),
            backbone_voxel_count=backbone_voxel_count, branch_labels=None, branches={},
        )

    branch_labels = label_components(remainder)
    branch_ids = np.unique(branch_labels)
    branch_ids = branch_ids[branch_ids >= 0]

    dist_field = dist_to_wall
    branches: dict[int, BranchInfo] = {}
    lumen_dilated = _touches(lumen_mask)
    for bid in branch_ids.tolist():
        bmask = branch_labels == bid
        b_voxel_count = int(bmask.sum())
        b_frac, b_defined = _wall_fraction_for_mask(bmask, protein_votes, lipid_votes)
        connected_to_lumen = bool((bmask & lumen_dilated).any())
        connected_to_upper = bool((bmask & _touches(top_bulk_void)).any())
        connected_to_lower = bool((bmask & _touches(bot_bulk_void)).any())
        b_dist = dist_field[bmask]
        bottleneck = float(b_dist.min())
        mean_r = float(b_dist.mean())
        max_r = float(b_dist.max())

        if not b_defined:
            # No wall votes reachable at all (e.g. a tiny, fully-interior
            # grid-quantization sliver with no voxel touching excluded
            # volume) — fail safe like _classify_region does, not a
            # membrane_defect default. An undefined wall composition is not
            # evidence of a lipid defect; it is an absence of evidence.
            classification = "ambiguous"
            confidence = 0.3
        elif b_frac >= vestibule_wall_threshold and connected_to_lumen:
            classification = REGION_VESTIBULE if (connected_to_upper or connected_to_lower) else REGION_ISOLATED_CAVITY
            headroom = max(1e-9, 1.0 - vestibule_wall_threshold)
            confidence = float(min(1.0, 0.4 + 0.4 * min(1.0, (b_frac - vestibule_wall_threshold) / headroom)))
        else:
            classification = REGION_MEMBRANE_DEFECT
            confidence = float(min(1.0, 0.5 + 0.5 * min(1.0, (vestibule_wall_threshold - b_frac) / vestibule_wall_threshold)))

        branches[bid] = BranchInfo(
            branch_id=bid,
            voxel_count=b_voxel_count,
            connected_to_lumen=connected_to_lumen,
            connected_to_upper_bulk=connected_to_upper,
            connected_to_lower_bulk=connected_to_lower,
            protein_wall_fraction=round(b_frac, 4),
            lipid_wall_fraction=round(1.0 - b_frac, 4) if b_defined else 0.5,
            wall_composition_defined=b_defined,
            bottleneck_radius_nm=round(bottleneck, 4),
            mean_radius_nm=round(mean_r, 4),
            max_radius_nm=round(max_r, 4),
            classification=classification,
            confidence=round(confidence, 4),
        )

    return DecompositionResult(
        attempted=True,
        split=True,
        lumen_mask=lumen_mask,
        lumen_protein_wall_fraction=round(lumen_frac, 4) if lumen_defined else whole_frac,
        lumen_classification=lumen_classification,
        lumen_confidence=round(lumen_confidence, 4),
        lumen_voxel_count=lumen_voxel_count,
        lumen_bottleneck_radius_nm=round(lumen_bottleneck, 4),
        lumen_mean_radius_nm=round(lumen_mean, 4),
        lumen_max_radius_nm=round(lumen_max, 4),
        backbone_voxel_count=backbone_voxel_count,
        branch_labels=branch_labels,
        branches=branches,
    )


def _touches(mask: np.ndarray) -> np.ndarray:
    """Boolean array: True at voxels 6-adjacent to (or equal to) a True
    voxel of `mask`."""
    out = mask.copy()
    for axis in range(3):
        for step in (1, -1):
            out |= _shift(mask, axis, step, fill=False)
    return out


def _grow_within(seed: np.ndarray, allowed: np.ndarray) -> np.ndarray:
    """6-connected flood-fill from `seed`, restricted to `allowed` voxels."""
    reached = seed & allowed
    while True:
        expanded = _touches(reached) & allowed
        if int(expanded.sum()) == int(reached.sum()):
            break
        reached = expanded
    return reached


def _grow_radius_limited(
    seed: np.ndarray,
    allowed: np.ndarray,
    dist_to_wall: np.ndarray,
    growth_factor: float,
    growth_abs_nm: float,
) -> np.ndarray:
    """6-connected flood-fill from `seed`, restricted to `allowed` voxels,
    with an additional per-edge admission test: a not-yet-accepted voxel `v`
    may join from an already-accepted neighbor `u` only if
    dist_to_wall[v] <= dist_to_wall[u] * growth_factor + growth_abs_nm.
    Gradual widening (e.g. a vestibule opening up after a narrow
    constriction) is admitted one modest step at a time; a single sharp
    jump into much wider space (e.g. escaping into open lateral lipid-facing
    void) is not."""
    accepted = seed & allowed
    while True:
        newly = np.zeros_like(accepted)
        for axis in range(3):
            for step in (1, -1):
                neighbor_accepted = _shift(accepted, axis, step, fill=False)
                neighbor_radius = _shift(dist_to_wall, axis, step, fill=0.0)
                cap = neighbor_radius * growth_factor + growth_abs_nm
                candidate = allowed & ~accepted & neighbor_accepted & (dist_to_wall <= cap)
                newly |= candidate
        if not newly.any():
            break
        accepted = accepted | newly
    return accepted
