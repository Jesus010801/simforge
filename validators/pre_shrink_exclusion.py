"""
validators/pre_shrink_exclusion.py

Cavity-aware pre-shrink lipid exclusion.

Removes complete lipid residues that are geometrically trapped inside the
protein TM cavity immediately after protein–bilayer merge (embed_in_bilayer)
and before topology generation.

Only removes lipids classified as geometrically surrounded by the protein TM
bundle — NOT surface-adjacent lipids that inflategro will resolve naturally.

No topology modification is performed here; the caller must ensure the
topology is regenerated from the filtered system.gro.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Optional

from validators.membrane_validators import detect_trapped_lipids


# Safety default: max lipids removable without an explicit safety-limit override.
# For a 512-lipid DPPC bilayer around a 7-TM GPCR, typically 20–40 cavity lipids
# are identified. 50 allows headroom while still blocking nonsensical deletions.
DEFAULT_SAFETY_LIMIT = 50
HARD_OVERLAP_NM = 0.2
_SOLVENT_NAMES = {"SOL", "HOH", "WAT", "NA", "CL", "K", "MG", "CA"}


def _parse_system(gro_path: Path, lipid_resname: str) -> tuple[list[dict], dict[int, list[dict]]]:
    protein_atoms: list[dict] = []
    lipid_atoms: dict[int, list[dict]] = {}
    for line in gro_path.read_text().splitlines()[2:-1]:
        if len(line) < 44:
            continue
        try:
            atom = {
                "resid": int(line[0:5]), "resname": line[5:10].strip(),
                "atomname": line[10:15].strip(), "x": float(line[20:28]),
                "y": float(line[28:36]), "z": float(line[36:44]),
            }
        except ValueError:
            continue
        if atom["resname"] == lipid_resname:
            lipid_atoms.setdefault(atom["resid"], []).append(atom)
        elif atom["resname"] not in _SOLVENT_NAMES:
            protein_atoms.append(atom)
    return protein_atoms, lipid_atoms


def _phosphorus_slab(lipid_atoms: dict[int, list[dict]]) -> tuple[float, float]:
    """Return mean lower/upper phosphorus planes, matching InflateGRO's slab."""
    phosphorus_z = [atom["z"] for atoms in lipid_atoms.values() for atom in atoms
                    if atom["atomname"].upper().startswith("P")]
    if len(phosphorus_z) < 2:
        raise ValueError("Cannot build membrane-exclusion mask: fewer than two lipid phosphorus atoms")
    middle = sum(phosphorus_z) / len(phosphorus_z)
    lower = [z for z in phosphorus_z if z < middle]
    upper = [z for z in phosphorus_z if z > middle]
    if not lower or not upper:
        raise ValueError("Cannot build membrane-exclusion mask: phosphorus leaflets could not be separated")
    return sum(lower) / len(lower), sum(upper) / len(upper)


def _spatial_index(atoms: list[dict], cell_size: float) -> dict[tuple[int, int, int], list[dict]]:
    index: dict[tuple[int, int, int], list[dict]] = {}
    for atom in atoms:
        cell = tuple(math.floor(atom[axis] / cell_size) for axis in ("x", "y", "z"))
        index.setdefault(cell, []).append(atom)
    return index


def _has_overlap(atoms: list[dict], index: dict, cutoff: float) -> bool:
    cutoff_sq = cutoff * cutoff
    for atom in atoms:
        cell = tuple(math.floor(atom[axis] / cutoff) for axis in ("x", "y", "z"))
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for other in index.get((cell[0] + dx, cell[1] + dy, cell[2] + dz), ()):
                        distance_sq = sum((atom[axis] - other[axis]) ** 2 for axis in ("x", "y", "z"))
                        if distance_sq < cutoff_sq:
                            return True
    return False


def _com_penetrates_mask(atoms: list[dict], index: dict, cutoff: float) -> bool:
    x, y, z = _lipid_com(atoms)
    probe = [{"x": x, "y": y, "z": z}]
    return _has_overlap(probe, index, cutoff)


def _convex_hull(points: list[tuple[float, float]]) -> list[tuple[float, float]]:
    points = sorted(set(points))
    if len(points) <= 2:
        return points
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    lower = []
    for point in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper = []
    for point in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return lower[:-1] + upper[:-1]


def _inside_hull(point: tuple[float, float], hull: list[tuple[float, float]]) -> bool:
    if len(hull) < 3:
        return False
    inside = False
    x, y = point
    for idx, current in enumerate(hull):
        previous = hull[idx - 1]
        if ((current[1] > y) != (previous[1] > y)) and (
            x < (previous[0] - current[0]) * (y - current[1])
            / (previous[1] - current[1]) + current[0]
        ):
            inside = not inside
    return inside


def _lipid_com(atoms: list[dict]) -> tuple[float, float, float]:
    count = len(atoms)
    return (sum(a["x"] for a in atoms) / count, sum(a["y"] for a in atoms) / count,
            sum(a["z"] for a in atoms) / count)


def _remove_residues_from_gro(
    gro_path: Path, resids: set[int], lipid_resname: str
) -> int:
    """
    Remove all atoms belonging to the given GRO residue numbers from gro_path.
    Rewrites the file in-place with the corrected atom-count header line.
    Returns the number of atoms removed.
    """
    lines = gro_path.read_text().splitlines()
    title   = lines[0]
    atom_lines = lines[2:-1]   # skip title, atom-count, box-vector
    box     = lines[-1]

    kept: list[str] = []
    n_removed = 0
    for line in atom_lines:
        if len(line) < 5:
            kept.append(line)
            continue
        try:
            resnum = int(line[0:5])
        except ValueError:
            kept.append(line)
            continue
        if resnum in resids and line[5:10].strip() == lipid_resname:
            n_removed += 1
        else:
            kept.append(line)

    new_text = f"{title}\n{len(kept)}\n" + "\n".join(kept) + f"\n{box}\n"
    gro_path.write_text(new_text)
    return n_removed


def apply_pre_shrink_exclusion(
    gro_path:         Path,
    lipid_resname:    str,
    tm_residues:      Optional[set[int]] = None,
    safety_limit:     int = DEFAULT_SAFETY_LIMIT,
    enabled:          bool = True,
    include_non_tm_slab_atoms: bool = False,
) -> dict:
    """
    Identify geometrically trapped lipid residues and remove them from gro_path.

    Criteria for removal (all required):
      1. Lipid COM lies between the mean lower/upper phosphorus planes.
      2. The lipid overlaps the TM-and-slab membrane mask.
      3. Its COM penetrates that mask within HARD_OVERLAP_NM.
      4. TM C-alpha atoms angularly surround the lipid COM (max gap < 185°).
    Full-protein contacts remain diagnostic-only, so soluble/out-of-slab atoms
    cannot authorize deletion of a normal annular neighbour.

    The GRO file is rewritten in-place with the corrected atom-count header.
    Residue renumbering (gmx editconf -resnr 1) is left to the caller.

    Args:
        gro_path:      Path to system.gro (modified in-place if enabled).
        lipid_resname: Residue name of lipid atoms (e.g. "DPP").
        tm_residues:   Optional set of protein residue numbers defining the TM
                       bundle. If None, all protein residues are used.
        safety_limit:  Maximum lipids removable in one pass. If the classifier
                       identifies more than this, raises ValueError.
        enabled:       If False, performs no action and returns a disabled report.
        include_non_tm_slab_atoms: Opt-in inclusion of non-TM protein atoms that
                       physically intersect the phosphorus-defined slab.

    Returns:
        dict with keys:
          enabled, n_lipids_checked, n_lipids_removed, removed_resids,
          removed_resnames, removed_reasons, n_hard_overlap_removed,
          n_low_overlap_geometrically_trapped_removed, safety_limit,
          safety_status

    Raises:
        ValueError if n_surrounded_candidates > safety_limit.
        FileNotFoundError if gro_path does not exist.
    """
    base_report = {
        "enabled":                               enabled,
        "n_lipids_checked":                      0,
        "n_lipids_removed":                      0,
        "removed_resids":                        [],
        "removed_resnames":                      [],
        "removed_reasons":                       [],
        "n_hard_overlap_removed":                0,
        "n_low_overlap_geometrically_trapped_removed": 0,
        "n_candidate_lipids":                    0,
        "n_removed_by_membrane_mask":            0,
        "n_rejected_out_of_slab":                0,
        "n_rejected_soluble_only_contact":       0,
        "membrane_mask_atom_count":              0,
        "membrane_mask_residue_count":           0,
        "slab_z_min":                            None,
        "slab_z_max":                            None,
        "dpp_before":                            0,
        "dpp_after":                             0,
        "annular_occupancy_before":              None,
        "annular_occupancy_after":               None,
        "safety_limit":                          safety_limit,
        "safety_status":                         "ok",
    }

    if not enabled:
        return base_report

    if not gro_path.exists():
        raise FileNotFoundError(f"GRO file not found: {gro_path}")

    protein_atoms, lipid_atoms = _parse_system(gro_path, lipid_resname)
    slab_z_min, slab_z_max = _phosphorus_slab(lipid_atoms)
    tm_set = set(tm_residues or ())
    # Z relevance is mandatory even for annotated TM residues: broad or
    # flanking TM annotations must not re-introduce soluble atoms into the
    # deletion mask. Optional non-TM atoms are admitted only where they truly
    # intersect the same physical membrane slab. The optional non-TM extension
    # is deliberately opt-in so soluble-domain contacts cannot delete annular
    # lipids by default.
    membrane_mask = [
        atom for atom in protein_atoms
        if slab_z_min <= atom["z"] <= slab_z_max
        and (atom["resid"] in tm_set or include_non_tm_slab_atoms or not tm_set)
    ]

    diag = detect_trapped_lipids(
        gro_path,
        lipid_resname=lipid_resname,
        tm_residues=tm_residues,
    )

    legacy_candidates = {c.resid: c for c in diag.lipid_residue_candidates}
    full_index = _spatial_index(protein_atoms, HARD_OVERLAP_NM)
    mask_index = _spatial_index(membrane_mask, HARD_OVERLAP_NM)
    to_remove: list[int] = []
    rejected_out_of_slab = 0
    rejected_soluble_only = 0
    candidate_lipids = 0
    for resid, atoms in lipid_atoms.items():
        full_overlap = _has_overlap(atoms, full_index, HARD_OVERLAP_NM)
        mask_overlap = _has_overlap(atoms, mask_index, HARD_OVERLAP_NM)
        com_mask_overlap = _com_penetrates_mask(atoms, mask_index, HARD_OVERLAP_NM)
        legacy_candidate = legacy_candidates.get(resid)
        is_surrounded = bool(legacy_candidate and legacy_candidate.is_surrounded_by_protein)
        if full_overlap or mask_overlap or legacy_candidate is not None:
            candidate_lipids += 1
        in_slab = slab_z_min <= _lipid_com(atoms)[2] <= slab_z_max
        if not in_slab:
            if full_overlap or mask_overlap or legacy_candidate is not None:
                rejected_out_of_slab += 1
            continue
        # Whole-lipid deletion requires molecular-centre penetration. A single
        # tail/headgroup atom contact is retained as a diagnostic candidate but
        # is not sufficient to delete an otherwise normal annular neighbour.
        if not mask_overlap or not com_mask_overlap or not is_surrounded:
            if full_overlap and not mask_overlap:
                rejected_soluble_only += 1
            continue
        to_remove.append(resid)

    n_to_remove = len(to_remove)
    mask_hull = _convex_hull([(a["x"], a["y"]) for a in membrane_mask])
    annular_before = sum(_inside_hull(_lipid_com(atoms)[:2], mask_hull)
                          for atoms in lipid_atoms.values())
    annular_removed = sum(_inside_hull(_lipid_com(lipid_atoms[r])[:2], mask_hull)
                           for r in to_remove)
    report_values = {
        "n_lipids_checked": len(lipid_atoms),
        "n_candidate_lipids": candidate_lipids,
        "n_removed_by_membrane_mask": n_to_remove,
        "n_rejected_out_of_slab": rejected_out_of_slab,
        "n_rejected_soluble_only_contact": rejected_soluble_only,
        "membrane_mask_atom_count": len(membrane_mask),
        "membrane_mask_residue_count": len({a["resid"] for a in membrane_mask}),
        "slab_z_min": round(slab_z_min, 4), "slab_z_max": round(slab_z_max, 4),
        "dpp_before": len(lipid_atoms), "dpp_after": len(lipid_atoms) - n_to_remove,
        "annular_occupancy_before": annular_before,
        "annular_occupancy_after": annular_before - annular_removed,
    }

    if n_to_remove > safety_limit:
        report = {
            **base_report,
            **report_values,
            "safety_status":    "exceeded",
        }
        raise ValueError(
            f"Pre-shrink exclusion would remove {n_to_remove} lipids "
            f"(exceeds safety_limit={safety_limit}). "
            "Too many lipids would be removed; check orientation, TM annotation, "
            "or bilayer selection."
        )

    if n_to_remove == 0:
        return {**base_report, **report_values}

    resids_to_remove = set(to_remove)
    removed_reasons = ["membrane_mask_com_penetration"] * n_to_remove
    n_hard = n_to_remove
    n_low_geo = 0

    _remove_residues_from_gro(gro_path, resids_to_remove, lipid_resname)

    return {
        **base_report,
        **report_values,
        "n_lipids_removed":                      n_to_remove,
        "removed_resids":                        sorted(resids_to_remove),
        "removed_resnames":                      [lipid_resname] * n_to_remove,
        "removed_reasons":                       removed_reasons,
        "n_hard_overlap_removed":                n_hard,
        "n_low_overlap_geometrically_trapped_removed": n_low_geo,
    }
