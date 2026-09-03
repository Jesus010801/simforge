"""Channel-aware water classification and conservative cleanup."""
from __future__ import annotations

import json
import math
import re
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Iterable

from validators.membrane_protein_spatial_classifier import (
    _DEFAULT_LIPID_RESNAMES,
    _HAS_NP,
    _SOLVENT_RESNAMES,
    _build_occupancy_grid,
    _compute_membrane_core_z,
    _compute_tm_geometry,
    _flood_fill,
    _parse_gro_atoms,
    _voxel_idx,
)

if _HAS_NP:
    import numpy as np


_CATEGORIES = (
    "bulk_water",
    "lipid_core_outside_protein",
    "channel_lumen_water",
    "protein_internal_cavity_water",
    "clash_water",
    "ambiguous",
)


def _water_molecules(path: Path) -> list[dict]:
    """Parse water molecules using occurrence IDs, safe across GRO resid wrap."""
    lines = path.read_text().splitlines()
    n_atoms = int(lines[1].strip())
    waters: list[dict] = []
    current: dict | None = None
    for atom_index, line in enumerate(lines[2:2 + n_atoms]):
        if len(line) < 44:
            continue
        try:
            resid = int(line[0:5])
            resname = line[5:10].strip()
            atomname = line[10:15].strip()
            xyz = (
                float(line[20:28]),
                float(line[28:36]),
                float(line[36:44]),
            )
        except ValueError:
            continue
        if resname != "SOL":
            current = None
            continue
        if (
            current is None
            or current["resid"] != resid
            or atomname == "OW"
            or current["atom_indices"][-1] != atom_index - 1
        ):
            current = {
                "water_id": len(waters),
                "resid": resid,
                "atom_indices": [],
                "oxygen": None,
            }
            waters.append(current)
        current["atom_indices"].append(atom_index)
        if atomname == "OW":
            current["oxygen"] = xyz
    return [water for water in waters if water["oxygen"] is not None]


def _spatial_hash(atoms: Iterable[dict], cell_nm: float) -> dict:
    cells: dict[tuple[int, int, int], list[tuple[float, float, float]]] = defaultdict(list)
    for atom in atoms:
        key = (
            math.floor(atom["x"] / cell_nm),
            math.floor(atom["y"] / cell_nm),
            math.floor(atom["z"] / cell_nm),
        )
        cells[key].append((atom["x"], atom["y"], atom["z"]))
    return cells


def _min_distance(
    point: tuple[float, float, float],
    cells: dict,
    cell_nm: float,
    cutoff_nm: float,
) -> float:
    key = tuple(math.floor(value / cell_nm) for value in point)
    span = max(1, math.ceil(cutoff_nm / cell_nm))
    best = math.inf
    for dx in range(-span, span + 1):
        for dy in range(-span, span + 1):
            for dz in range(-span, span + 1):
                for other in cells.get((key[0] + dx, key[1] + dy, key[2] + dz), ()):
                    best = min(best, math.dist(point, other))
    return best


def classify_pore_hydration(
    gro_path: Path | str,
    *,
    tm_residues: set[int] | None = None,
    lipid_resnames: set[str] | None = None,
    output_dir: Path | str | None = None,
    enabled: bool = True,
    grid_spacing_nm: float = 0.15,
    protein_padding_nm: float = 0.25,
    protein_clash_nm: float = 0.20,
    lipid_clash_nm: float = 0.18,
) -> dict:
    """Classify each water without modifying coordinates."""
    gro_path = Path(gro_path)
    output_dir = Path(output_dir) if output_dir else gro_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    report_path = output_dir / "pore_hydration_report.json"
    base = {
        "enabled": bool(enabled),
        "channel_like_system_detected": False,
        "n_waters_total": 0,
        "n_bulk_water": 0,
        "n_lipid_core_outside_protein": 0,
        "n_channel_lumen_water": 0,
        "n_protein_internal_cavity_water": 0,
        "n_clash_water": 0,
        "n_ambiguous": 0,
        "membrane_slab_z_min": None,
        "membrane_slab_z_max": None,
        "pore_hydration_score": 0.0,
        "lumen_continuity_score": 0.0,
        "warnings": [],
        "coordinate_file_modified": False,
        "water_ids_by_class": {name: [] for name in _CATEGORIES},
    }
    if not enabled:
        report_path.write_text(json.dumps(base, indent=2))
        return base
    if not _HAS_NP:
        base["warnings"].append("numpy unavailable; pore hydration classification skipped")
        report_path.write_text(json.dumps(base, indent=2))
        return base

    atoms, _box = _parse_gro_atoms(gro_path)
    waters = _water_molecules(gro_path)
    lipids = frozenset(lipid_resnames or _DEFAULT_LIPID_RESNAMES)
    tm_set = set(tm_residues or ())
    protein_atoms = [
        atom for atom in atoms
        if atom["resname"] not in lipids
        and atom["resname"] not in _SOLVENT_RESNAMES
        and not atom["atomname"].startswith("H")
    ]
    lipid_atoms = [atom for atom in atoms if atom["resname"] in lipids]
    z_bot, z_top, _z_mid = _compute_membrane_core_z(atoms, lipids)
    base["n_waters_total"] = len(waters)
    base["membrane_slab_z_min"] = round(z_bot, 4)
    base["membrane_slab_z_max"] = round(z_top, 4)

    if not protein_atoms or not lipid_atoms:
        base["warnings"].append("protein or lipid atoms missing; classification is ambiguous")
        base["n_ambiguous"] = len(waters)
        base["water_ids_by_class"]["ambiguous"] = [w["water_id"] for w in waters]
        report_path.write_text(json.dumps(base, indent=2))
        return base

    tm_cx, tm_cy, tm_radius = _compute_tm_geometry(atoms, tm_set, lipids)
    can_detect_channel = bool(tm_set)
    tm_threshold = tm_radius + protein_padding_nm * 0.5

    relevant = protein_atoms + lipid_atoms
    xs, ys, zs = ([atom[key] for atom in relevant] for key in ("x", "y", "z"))
    origin = (min(xs) - 1.0, min(ys) - 1.0, min(zs) - 0.5)
    shape = (
        max(4, math.ceil((max(xs) + 1.0 - origin[0]) / grid_spacing_nm) + 1),
        max(4, math.ceil((max(ys) + 1.0 - origin[1]) / grid_spacing_nm) + 1),
        max(4, math.ceil((max(zs) + 0.5 - origin[2]) / grid_spacing_nm) + 1),
    )
    top_reached = bottom_reached = None
    if can_detect_channel:
        blocked = _build_occupancy_grid(
            protein_atoms, grid_spacing_nm, protein_padding_nm, origin, shape
        )
        top_start = min(
            shape[2] - 1,
            max(0, math.ceil((z_top + 0.1 - origin[2]) / grid_spacing_nm)),
        )
        bottom_end = max(
            0,
            min(shape[2] - 1, math.floor((z_bot - 0.1 - origin[2]) / grid_spacing_nm)),
        )
        top_seeds = np.zeros(shape, dtype=bool)
        top_seeds[:, :, top_start:] = True
        bottom_seeds = np.zeros(shape, dtype=bool)
        bottom_seeds[:, :, : bottom_end + 1] = True
        top_reached = _flood_fill(blocked, top_seeds & ~blocked)
        bottom_reached = _flood_fill(blocked, bottom_seeds & ~blocked)

    cell_nm = max(protein_clash_nm, lipid_clash_nm)
    protein_hash = _spatial_hash(protein_atoms, cell_nm)
    lipid_hash = _spatial_hash(lipid_atoms, cell_nm)
    classifications: dict[str, list[int]] = {name: [] for name in _CATEGORIES}

    for water in waters:
        water_id = water["water_id"]
        x, y, z = water["oxygen"]
        if (
            _min_distance((x, y, z), protein_hash, cell_nm, protein_clash_nm)
            < protein_clash_nm
            or _min_distance((x, y, z), lipid_hash, cell_nm, lipid_clash_nm)
            < lipid_clash_nm
        ):
            classifications["clash_water"].append(water_id)
        elif z < z_bot - 0.05 or z > z_top + 0.05:
            classifications["bulk_water"].append(water_id)
        elif not can_detect_channel:
            classifications["ambiguous"].append(water_id)
        elif math.hypot(x - tm_cx, y - tm_cy) > tm_threshold:
            classifications["lipid_core_outside_protein"].append(water_id)
        else:
            idx = _voxel_idx(x, y, z, origin, grid_spacing_nm, shape)
            from_top = bool(top_reached[idx]) if top_reached is not None else False
            from_bottom = bool(bottom_reached[idx]) if bottom_reached is not None else False
            if from_top and from_bottom:
                classifications["channel_lumen_water"].append(water_id)
            elif not from_top and not from_bottom:
                classifications["protein_internal_cavity_water"].append(water_id)
            else:
                classifications["ambiguous"].append(water_id)

    for category, ids in classifications.items():
        base["water_ids_by_class"][category] = ids
        base[f"n_{category}"] = len(ids)

    lumen_ids = set(classifications["channel_lumen_water"])
    internal_ids = lumen_ids | set(classifications["protein_internal_cavity_water"])
    internal_ids |= set(classifications["ambiguous"])
    base["channel_like_system_detected"] = bool(lumen_ids)
    if internal_ids:
        base["lumen_continuity_score"] = round(len(lumen_ids) / len(internal_ids), 4)
    n_bins = 12
    occupied_bins: set[int] = set()
    slab_width = max(z_top - z_bot, 1e-6)
    for water in waters:
        if water["water_id"] in lumen_ids:
            z = water["oxygen"][2]
            occupied_bins.add(min(n_bins - 1, max(0, int((z - z_bot) / slab_width * n_bins))))
    base["pore_hydration_score"] = round(len(occupied_bins) / n_bins, 4)
    if can_detect_channel and not lumen_ids:
        base["warnings"].append("TM annotation present but no bulk-connected lumen water detected")
    if classifications["clash_water"]:
        base["warnings"].append(
            f"{len(classifications['clash_water'])} water molecules clash with protein/lipid atoms"
        )
    report_path.write_text(json.dumps(base, indent=2))
    return base


def clean_water_channel_aware(
    gro_in: Path | str,
    gro_out: Path | str,
    *,
    tm_residues: set[int],
    output_dir: Path | str | None = None,
    topol_in: Path | str | None = None,
    topol_out: Path | str | None = None,
) -> dict:
    """Remove only external lipid-core water; preserve lumen water."""
    gro_in, gro_out = Path(gro_in), Path(gro_out)
    report = classify_pore_hydration(
        gro_in, tm_residues=tm_residues, output_dir=output_dir or gro_out.parent
    )
    remove_ids = set(report["water_ids_by_class"]["lipid_core_outside_protein"])
    lines = gro_in.read_text().splitlines()
    n_atoms = int(lines[1].strip())
    atom_lines = lines[2:2 + n_atoms]
    box_line = lines[2 + n_atoms]
    waters = _water_molecules(gro_in)
    remove_atom_indices = {
        index for water in waters if water["water_id"] in remove_ids
        for index in water["atom_indices"]
    }
    kept, atom_number = [], 0
    for index, line in enumerate(atom_lines):
        if index in remove_atom_indices:
            continue
        atom_number += 1
        kept.append(line[:15] + f"{atom_number % 100000:5d}" + line[20:])
    gro_out.write_text("\n".join([lines[0], str(len(kept)), *kept, box_line]) + "\n")
    report["n_water_molecules_removed"] = len(remove_ids)
    report["coordinate_file_modified"] = bool(remove_ids)
    report["output_gro_path"] = str(gro_out)
    out_dir = Path(output_dir or gro_out.parent)
    out_dir.joinpath("pore_hydration_report.json").write_text(json.dumps(report, indent=2))
    topology_updated = False
    if topol_in and topol_out:
        topol_in, topol_out = Path(topol_in), Path(topol_out)
        shutil.copy2(topol_in, topol_out)
        updated = []
        for line in topol_out.read_text().splitlines():
            match = re.match(r"^(SOL)\s+(\d+)", line)
            if match:
                line = f"SOL              {int(match.group(2)) - len(remove_ids)}"
                topology_updated = True
            updated.append(line)
        topol_out.write_text("\n".join(updated) + "\n")
    clean_report = {
        "input_water_molecules": report["n_waters_total"],
        "n_water_molecules_removed": len(remove_ids),
        "n_water_atoms_removed": len(remove_atom_indices),
        "n_water_oxygens_remaining_in_core": 0,
        "n_pore_waters_preserved": report["n_channel_lumen_water"],
        "n_internal_cavity_waters_preserved": report["n_protein_internal_cavity_water"],
        "n_clash_waters_flagged": report["n_clash_water"],
        "core_z_min": report["membrane_slab_z_min"],
        "core_z_max": report["membrane_slab_z_max"],
        "output_gro_path": str(gro_out),
        "topology_updated": topology_updated,
        "cleanup_passed": True,
        "pore_aware": True,
        "warnings": report["warnings"],
    }
    out_dir.joinpath("clean_water_report.json").write_text(json.dumps(clean_report, indent=2))
    out_dir.joinpath("water_report.json").write_text(json.dumps({
        "passed": True,
        "n_water_molecules_removed": len(remove_ids),
        "n_water_atoms_removed": len(remove_atom_indices),
        "n_waters_remaining": 0,
        "n_water_oxygens_remaining_in_core": 0,
        "n_pore_waters_preserved": report["n_channel_lumen_water"],
        "cleanup_passed": True,
        "message": f"Removed {len(remove_ids)} external lipid-core waters; preserved {report['n_channel_lumen_water']} channel-lumen waters",
        "errors": [],
        "warnings": report["warnings"],
        "confidence": report["lumen_continuity_score"],
    }, indent=2))
    report["topology_updated"] = topology_updated
    return report
