"""
validators/membrane_water/regions.py

Per-connected-component region analysis and the pore-vs-lateral-defect
discrimination that is the actual fix for the reported bug: a connected
void component that spans both bulk sides is only called a "pore" if its
membrane-core boundary is predominantly protein-walled. If it is
predominantly lipid-walled it is a "membrane_defect" — a lateral bulk-to-bulk
shortcut through a lipid packing gap — regardless of how large or
bulk-connected it is. This uses no cylinder, no fixed radius, no residue
ranges: only local wall composition (from the excluded-volume wall-species
grid), bulk connectivity (component-level, from the flood-fill masks), and a
secondary discrete-distance "bottleneck" plausibility check.
"""
from __future__ import annotations

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
from validators.membrane_water.voidgraph import distance_to_blocked, pbc_shift


@dataclass
class RegionInfo:
    region_id: int
    classification: str
    voxel_count: int
    connected_to_upper_bulk: bool
    connected_to_lower_bulk: bool
    traverses_membrane: bool
    protein_wall_fraction: float
    lipid_wall_fraction: float
    wall_composition_defined: bool
    bottleneck_radius_nm: float
    mean_radius_nm: float
    max_radius_nm: float
    confidence: float


def _neighbor_wall_votes(blocked: np.ndarray, wall_species: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """For every voxel, count how many of its 6 neighbors are blocked-protein
    vs. blocked-lipid voxels."""
    is_protein = blocked & (wall_species == WALL_PROTEIN)
    is_lipid = blocked & (wall_species == WALL_LIPID)
    protein_votes = np.zeros(blocked.shape, dtype=np.int32)
    lipid_votes = np.zeros(blocked.shape, dtype=np.int32)
    for axis in range(3):
        for step in (1, -1):
            protein_votes += pbc_shift(is_protein, axis, step, fill=False).astype(np.int32)
            lipid_votes += pbc_shift(is_lipid, axis, step, fill=False).astype(np.int32)
    return protein_votes, lipid_votes


def _classify_region(
    connected_upper: bool,
    connected_lower: bool,
    protein_frac: float,
    wall_defined: bool,
    bottleneck_nm: float,
    protein_wall_threshold: float,
    vestibule_wall_threshold: float,
) -> tuple[str, float]:
    if not wall_defined:
        return "ambiguous", 0.3

    if connected_upper and connected_lower:
        if protein_frac >= protein_wall_threshold:
            classification = REGION_PORE
            headroom = max(1e-9, 1.0 - protein_wall_threshold)
            confidence = 0.5 + 0.5 * min(1.0, (protein_frac - protein_wall_threshold) / headroom)
        else:
            classification = REGION_MEMBRANE_DEFECT
            confidence = 0.5 + 0.5 * min(1.0, (protein_wall_threshold - protein_frac) / protein_wall_threshold)
    elif connected_upper or connected_lower:
        if protein_frac >= vestibule_wall_threshold:
            classification = REGION_VESTIBULE
            headroom = max(1e-9, 1.0 - vestibule_wall_threshold)
            confidence = 0.4 + 0.4 * min(1.0, (protein_frac - vestibule_wall_threshold) / headroom)
        else:
            classification = REGION_MEMBRANE_DEFECT
            confidence = 0.4 + 0.4 * min(1.0, (vestibule_wall_threshold - protein_frac) / vestibule_wall_threshold)
    else:
        classification = REGION_ISOLATED_CAVITY
        confidence = 0.6

    if classification in (REGION_PORE, REGION_VESTIBULE) and bottleneck_nm < 0.10:
        confidence *= 0.7  # implausibly narrow for a water molecule; lower confidence, don't reclassify

    return classification, float(min(1.0, confidence))


def analyze_regions(
    labels: np.ndarray,
    blocked: np.ndarray,
    wall_species: np.ndarray,
    top_reached: np.ndarray,
    bottom_reached: np.ndarray,
    spacing: float,
    protein_wall_threshold: float = DEFAULT_PROTEIN_WALL_THRESHOLD,
    vestibule_wall_threshold: float = DEFAULT_VESTIBULE_WALL_THRESHOLD,
) -> dict[int, RegionInfo]:
    unique_labels = np.unique(labels)
    unique_labels = unique_labels[unique_labels >= 0]
    if len(unique_labels) == 0:
        return {}

    protein_votes, lipid_votes = _neighbor_wall_votes(blocked, wall_species)
    dist = distance_to_blocked(blocked, spacing)

    regions: dict[int, RegionInfo] = {}
    for label_id in unique_labels.tolist():
        mask = labels == label_id
        voxel_count = int(mask.sum())
        if voxel_count == 0:
            continue
        connected_upper = bool(top_reached[mask].any())
        connected_lower = bool(bottom_reached[mask].any())
        p_votes = int(protein_votes[mask].sum())
        l_votes = int(lipid_votes[mask].sum())
        total_votes = p_votes + l_votes
        wall_defined = total_votes > 0
        protein_frac = (p_votes / total_votes) if wall_defined else 0.5
        lipid_frac = 1.0 - protein_frac if wall_defined else 0.5
        region_dist = dist[mask]
        bottleneck_nm = float(region_dist.min())

        classification, confidence = _classify_region(
            connected_upper, connected_lower, protein_frac, wall_defined,
            bottleneck_nm, protein_wall_threshold, vestibule_wall_threshold,
        )

        regions[label_id] = RegionInfo(
            region_id=label_id,
            classification=classification,
            voxel_count=voxel_count,
            connected_to_upper_bulk=connected_upper,
            connected_to_lower_bulk=connected_lower,
            traverses_membrane=connected_upper and connected_lower,
            protein_wall_fraction=round(protein_frac, 4),
            lipid_wall_fraction=round(lipid_frac, 4),
            wall_composition_defined=wall_defined,
            bottleneck_radius_nm=round(bottleneck_nm, 4),
            # Lightweight medial-axis-adjacent stats (distance-to-nearest-wall
            # per voxel, aggregated over the region) — reusable later for a
            # pore-radius profile without needing to persist the full 3D
            # distance field. The full per-voxel field is available via
            # voidgraph.distance_to_blocked(blocked, spacing) for anyone who
            # needs it (e.g. a medial-axis/radius-profile tool).
            mean_radius_nm=round(float(region_dist.mean()), 4),
            max_radius_nm=round(float(region_dist.max()), 4),
            confidence=round(confidence, 4),
        )
    return regions
