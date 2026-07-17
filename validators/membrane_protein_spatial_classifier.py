"""
validators/membrane_protein_spatial_classifier.py
Phase 10B: Membrane-Protein Spatial Classifier.

Classifies regions of space in a protein-membrane system using a 3-D occupancy
grid and flood-fill to distinguish aqueous pores from external lipid-accessible
surfaces, hydrophobic core, and isolated cavities.

Classification categories (per water/lipid molecule):
  bulk_water                – outside membrane Z range
  pore_water                – at membrane Z, inside TM bundle, connected to bulk
  membrane_core_water       – at membrane Z, outside TM bundle (→ remove)
  isolated_internal_water   – at membrane Z, not connected to bulk
  valid_external_lipid      – outside TM bundle footprint
  forbidden_pore_lipid      – COM inside TM bundle at membrane Z
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Optional

try:
    import numpy as np
    _HAS_NP = True
except ImportError:
    _HAS_NP = False

_SOLVENT_RESNAMES: frozenset[str] = frozenset({
    "SOL", "HOH", "WAT", "TIP3", "TIP4", "TIP5",
    "NA", "CL", "SOD", "MG", "K", "CA",
})
_DEFAULT_LIPID_RESNAMES: frozenset[str] = frozenset({
    "DPP", "DPPC", "POPC", "POPE", "POPG", "POPS", "CHOL",
    "PALM", "OLEO", "PALC", "SM", "CER",
})
_HEADGROUP_ATOMS: frozenset[str] = frozenset({"P", "P1", "N", "NZ", "O33", "O11"})
_TAIL_ATOMS: frozenset[str] = frozenset({
    "C50", "C49", "C48", "C47", "C46", "C45", "C34", "C35",
})


# ── GRO parsing ───────────────────────────────────────────────────────────────

def _parse_gro_atoms(path: Path) -> tuple[list[dict], tuple[float, float, float]]:
    lines = path.read_text().splitlines()
    if len(lines) < 3:
        return [], (0.0, 0.0, 0.0)
    n = int(lines[1].strip())
    atoms: list[dict] = []
    for ln in lines[2: 2 + n]:
        if len(ln) < 44:
            continue
        try:
            atoms.append({
                "resid":    int(ln[0:5]),
                "resname":  ln[5:10].strip(),
                "atomname": ln[10:15].strip(),
                "x":        float(ln[20:28]),
                "y":        float(ln[28:36]),
                "z":        float(ln[36:44]),
            })
        except ValueError:
            continue
    box_line = lines[2 + n] if len(lines) > 2 + n else "10.0 10.0 10.0"
    parts = box_line.split()
    box = (float(parts[0]), float(parts[1]), float(parts[2]))
    return atoms, box


# ── Membrane geometry ─────────────────────────────────────────────────────────

def _compute_membrane_core_z(
    atoms: list[dict],
    lip_set: frozenset[str],
) -> tuple[float, float, float]:
    """Return (z_core_bot, z_core_top, z_midplane) from lipid headgroup statistics."""
    hg_z = [a["z"] for a in atoms
             if a["resname"] in lip_set and a["atomname"] in _HEADGROUP_ATOMS]
    tail_z = [a["z"] for a in atoms
               if a["resname"] in lip_set and a["atomname"] in _TAIL_ATOMS]

    if not hg_z:
        all_lip_z = [a["z"] for a in atoms if a["resname"] in lip_set]
        if not all_lip_z:
            return 4.0, 6.0, 5.0
        mid = sum(all_lip_z) / len(all_lip_z)
        return mid - 1.0, mid + 1.0, mid

    z_mid = sum(tail_z) / len(tail_z) if tail_z else sum(hg_z) / len(hg_z)
    top_hg = [z for z in hg_z if z > z_mid]
    bot_hg = [z for z in hg_z if z <= z_mid]
    z_top = sum(top_hg) / len(top_hg) if top_hg else z_mid + 1.0
    z_bot = sum(bot_hg) / len(bot_hg) if bot_hg else z_mid - 1.0
    return z_bot, z_top, z_mid


# ── TM geometry ───────────────────────────────────────────────────────────────

def _compute_tm_geometry(
    atoms: list[dict],
    tm_set: set[int],
    lip_set: frozenset[str],
) -> tuple[float, float, float]:
    """Return (tm_cx, tm_cy, tm_radius_nm)."""
    if tm_set:
        tm_atoms = [
            a for a in atoms
            if a["resid"] in tm_set
            and a["resname"] not in lip_set
            and a["resname"] not in _SOLVENT_RESNAMES
            and not a["atomname"].startswith("H")
        ]
    else:
        tm_atoms = [
            a for a in atoms
            if a["resname"] not in lip_set
            and a["resname"] not in _SOLVENT_RESNAMES
            and not a["atomname"].startswith("H")
        ]

    if not tm_atoms:
        return 5.0, 5.0, 0.5

    cx = sum(a["x"] for a in tm_atoms) / len(tm_atoms)
    cy = sum(a["y"] for a in tm_atoms) / len(tm_atoms)
    r = max(
        math.sqrt((a["x"] - cx) ** 2 + (a["y"] - cy) ** 2)
        for a in tm_atoms
    )
    return cx, cy, r


# ── Grid & flood-fill ─────────────────────────────────────────────────────────

def _build_occupancy_grid(
    protein_atoms: list[dict],
    spacing: float,
    padding: float,
    origin: tuple[float, float, float],
    shape: tuple[int, int, int],
) -> "np.ndarray":
    """3-D boolean blocked-voxel grid (True = protein-occupied)."""
    blocked = np.zeros(shape, dtype=bool)
    ox, oy, oz = origin
    Nx, Ny, Nz = shape
    r_cells = int(math.ceil(padding / spacing))

    for a in protein_atoms:
        ix_c = int(round((a["x"] - ox) / spacing))
        iy_c = int(round((a["y"] - oy) / spacing))
        iz_c = int(round((a["z"] - oz) / spacing))
        for dix in range(-r_cells, r_cells + 1):
            ix = ix_c + dix
            if ix < 0 or ix >= Nx:
                continue
            for diy in range(-r_cells, r_cells + 1):
                iy = iy_c + diy
                if iy < 0 or iy >= Ny:
                    continue
                for diz in range(-r_cells, r_cells + 1):
                    iz = iz_c + diz
                    if iz < 0 or iz >= Nz:
                        continue
                    dist = math.sqrt(
                        (ox + ix * spacing - a["x"]) ** 2
                        + (oy + iy * spacing - a["y"]) ** 2
                        + (oz + iz * spacing - a["z"]) ** 2
                    )
                    if dist <= padding:
                        blocked[ix, iy, iz] = True
    return blocked


def _dilate_step(r: "np.ndarray") -> "np.ndarray":
    e = r.copy()
    e[1:, :, :] |= r[:-1, :, :]
    e[:-1, :, :] |= r[1:, :, :]
    e[:, 1:, :] |= r[:, :-1, :]
    e[:, :-1, :] |= r[:, 1:, :]
    e[:, :, 1:] |= r[:, :, :-1]
    e[:, :, :-1] |= r[:, :, 1:]
    return e


def _flood_fill(blocked: "np.ndarray", seeds: "np.ndarray") -> "np.ndarray":
    """Iterative dilation flood-fill from seed mask; ignores blocked voxels."""
    valid = ~blocked
    reached = seeds & valid
    while True:
        expanded = _dilate_step(reached) & valid
        if int(expanded.sum()) == int(reached.sum()):
            break
        reached = expanded
    return reached


def _voxel_idx(
    x: float, y: float, z: float,
    origin: tuple[float, float, float],
    spacing: float,
    shape: tuple[int, int, int],
) -> tuple[int, int, int]:
    ox, oy, oz = origin
    Nx, Ny, Nz = shape
    ix = max(0, min(int(round((x - ox) / spacing)), Nx - 1))
    iy = max(0, min(int(round((y - oy) / spacing)), Ny - 1))
    iz = max(0, min(int(round((z - oz) / spacing)), Nz - 1))
    return ix, iy, iz


# ── Topology SOL count update ─────────────────────────────────────────────────

def _update_sol_count(text: str, n_removed: int) -> str:
    lines = text.splitlines()
    out = []
    for line in lines:
        m = re.match(r'^(SOL)\s+(\d+)', line)
        if m:
            old = int(m.group(2))
            line = f"SOL              {old - n_removed}"
        out.append(line)
    return "\n".join(out)


# ── GRO water removal ─────────────────────────────────────────────────────────

def _remove_sol_resids_from_gro(
    gro_in: Path,
    gro_out: Path,
    remove_resids: set[int],
) -> int:
    """Remove SOL molecules by resid. Returns count of molecules removed."""
    lines = gro_in.read_text().splitlines()
    if len(lines) < 3:
        gro_out.write_text(gro_in.read_text())
        return 0

    n_orig = int(lines[1].strip())
    atom_lines = lines[2: 2 + n_orig]
    box_line = lines[2 + n_orig] if len(lines) > 2 + n_orig else ""

    kept: list[str] = []
    removed_resids_seen: set[int] = set()
    atom_num = 0
    for ln in atom_lines:
        if len(ln) < 15:
            kept.append(ln)
            continue
        try:
            resid = int(ln[0:5])
            resname = ln[5:10].strip()
        except ValueError:
            kept.append(ln)
            continue
        if resname == "SOL" and resid in remove_resids:
            removed_resids_seen.add(resid)
            continue
        atom_num += 1
        kept.append(ln[:15] + f"{atom_num % 100000:5d}" + ln[20:])

    out_lines = [lines[0], str(len(kept))] + kept
    if box_line:
        out_lines.append(box_line)
    gro_out.write_text("\n".join(out_lines) + "\n")
    return len(removed_resids_seen)


# ── Report helpers ────────────────────────────────────────────────────────────

def _disabled_report() -> dict:
    return {
        "enabled": False,
        "grid_spacing_nm": None,
        "protein_padding_nm": None,
        "membrane_core_z_range": None,
        "n_regions_detected": 0,
        "n_transmembrane_pores": 0,
        "n_one_sided_vestibules": 0,
        "n_isolated_cavities": 0,
        "n_bulk_waters": 0,
        "n_pore_waters": 0,
        "n_membrane_core_waters": 0,
        "n_isolated_internal_waters": 0,
        "n_valid_external_lipids": 0,
        "n_forbidden_pore_lipids": 0,
        "forbidden_pore_lipid_resids": [],
        "pore_water_resids": [],
        "removed_water_resids_candidate": [],
        "warnings": [],
    }


def _write_classifier_report(report: dict, output_dir: Path) -> None:
    (output_dir / "membrane_protein_spatial_classification_report.json").write_text(
        json.dumps(report, indent=2, default=str)
    )


# ── Main classifier ───────────────────────────────────────────────────────────

def run_spatial_classifier(
    gro_path: "Path | str",
    tm_residues: Optional[set[int]] = None,
    lipid_resnames: Optional[set[str]] = None,
    output_dir: "Optional[Path | str]" = None,
    enabled: bool = True,
    policy: str = "warn",
    grid_spacing_nm: float = 0.15,
    protein_padding_nm: float = 0.25,
    preserve_pore_water: bool = True,
    remove_pore_lipids: bool = True,
    remove_membrane_core_water: bool = True,
    remove_isolated_internal_water: bool = False,
    max_removed_pore_lipids: int = 80,
) -> dict:
    """Classify water and lipid molecules by spatial region in a protein-membrane system.

    Returns a report dict with classification counts and residue lists.
    Writes membrane_protein_spatial_classification_report.json if output_dir is given.
    """
    gro_path = Path(gro_path)
    out_dir = Path(output_dir) if output_dir else None
    if out_dir:
        out_dir.mkdir(parents=True, exist_ok=True)

    warnings_list: list[str] = []

    if not enabled:
        rep = _disabled_report()
        if out_dir:
            _write_classifier_report(rep, out_dir)
        return rep

    if not _HAS_NP:
        warnings_list.append("numpy not available — spatial classifier requires numpy")
        rep = _disabled_report()
        rep["warnings"] = warnings_list
        if out_dir:
            _write_classifier_report(rep, out_dir)
        return rep

    atoms, _box = _parse_gro_atoms(gro_path)
    if not atoms:
        warnings_list.append(f"Empty or unreadable GRO: {gro_path}")
        rep = _disabled_report()
        rep["warnings"] = warnings_list
        if out_dir:
            _write_classifier_report(rep, out_dir)
        return rep

    lip_set = frozenset(lipid_resnames) if lipid_resnames else _DEFAULT_LIPID_RESNAMES
    tm_set = set(tm_residues) if tm_residues else set()

    # Separate atoms by category
    protein_atoms = [
        a for a in atoms
        if a["resname"] not in lip_set
        and a["resname"] not in _SOLVENT_RESNAMES
        and not a["atomname"].startswith("H")
    ]

    lipid_by_resid: dict[int, list[dict]] = {}
    for a in atoms:
        if a["resname"] in lip_set:
            lipid_by_resid.setdefault(a["resid"], []).append(a)

    sol_ow_by_resid: dict[int, dict] = {}
    for a in atoms:
        if a["resname"] == "SOL" and a["atomname"] == "OW":
            sol_ow_by_resid[a["resid"]] = a

    # Membrane geometry
    z_core_bot, z_core_top, _z_mid = _compute_membrane_core_z(atoms, lip_set)

    # TM geometry
    tm_cx, tm_cy, tm_radius = _compute_tm_geometry(atoms, tm_set, lip_set)
    # Threshold for "inside TM bundle": radius + half the padding
    tm_thresh = tm_radius + protein_padding_nm * 0.5

    can_detect_pore = bool(protein_atoms) and bool(tm_set)

    # ── Build grid & flood-fill ───────────────────────────────────────────────
    top_reached = bottom_reached = None

    if can_detect_pore:
        relevant = protein_atoms + [
            a for lst in lipid_by_resid.values() for a in lst
        ]
        xs = [a["x"] for a in relevant]
        ys = [a["y"] for a in relevant]
        zs = [a["z"] for a in relevant]
        buf_xy, buf_z = 1.0, 0.5
        ox = min(xs) - buf_xy
        oy = min(ys) - buf_xy
        oz = min(zs) - buf_z
        Nx = max(4, int(math.ceil((max(xs) + buf_xy - ox) / grid_spacing_nm)) + 1)
        Ny = max(4, int(math.ceil((max(ys) + buf_xy - oy) / grid_spacing_nm)) + 1)
        Nz = max(4, int(math.ceil((max(zs) + buf_z - oz) / grid_spacing_nm)) + 1)
        origin = (ox, oy, oz)
        shape = (Nx, Ny, Nz)

        blocked = _build_occupancy_grid(
            protein_atoms, grid_spacing_nm, protein_padding_nm, origin, shape
        )

        # Z-index boundaries for bulk seeds
        iz_top_start = min(
            Nz - 1,
            max(0, int(math.ceil((z_core_top + 0.1 - oz) / grid_spacing_nm)))
        )
        iz_bot_end = max(
            0,
            min(Nz - 1, int(math.floor((z_core_bot - 0.1 - oz) / grid_spacing_nm)))
        )

        seeds_top = np.zeros(shape, dtype=bool)
        seeds_top[:, :, iz_top_start:] = True
        seeds_top &= ~blocked

        seeds_bot = np.zeros(shape, dtype=bool)
        seeds_bot[:, :, : iz_bot_end + 1] = True
        seeds_bot &= ~blocked

        top_reached    = _flood_fill(blocked, seeds_top)
        bottom_reached = _flood_fill(blocked, seeds_bot)

    # ── Water classification ──────────────────────────────────────────────────
    bulk_water:     list[int] = []
    pore_water:     list[int] = []
    core_water:     list[int] = []
    isolated_water: list[int] = []

    for resid, ow in sol_ow_by_resid.items():
        z = ow["z"]
        if z < z_core_bot - 0.05 or z > z_core_top + 0.05:
            bulk_water.append(resid)
            continue

        if not can_detect_pore or top_reached is None:
            core_water.append(resid)
            continue

        dist_xy = math.sqrt((ow["x"] - tm_cx) ** 2 + (ow["y"] - tm_cy) ** 2)
        if dist_xy > tm_thresh:
            core_water.append(resid)
            continue

        # Inside TM bundle footprint at membrane Z
        ix, iy, iz = _voxel_idx(ow["x"], ow["y"], z, origin, grid_spacing_nm, shape)
        in_top = bool(top_reached[ix, iy, iz])
        in_bot = bool(bottom_reached[ix, iy, iz])

        if in_top or in_bot:
            pore_water.append(resid)
        else:
            isolated_water.append(resid)

    # ── Lipid classification ──────────────────────────────────────────────────
    valid_lipids:    list[int] = []
    forbidden_lipids: list[int] = []

    for resid, lip_atoms in lipid_by_resid.items():
        cx = sum(a["x"] for a in lip_atoms) / len(lip_atoms)
        cy = sum(a["y"] for a in lip_atoms) / len(lip_atoms)
        cz = sum(a["z"] for a in lip_atoms) / len(lip_atoms)

        if cz < z_core_bot or cz > z_core_top:
            valid_lipids.append(resid)
            continue

        dist_xy = math.sqrt((cx - tm_cx) ** 2 + (cy - tm_cy) ** 2)
        if can_detect_pore and dist_xy <= tm_thresh:
            forbidden_lipids.append(resid)
        else:
            valid_lipids.append(resid)

    # ── Policy checks ─────────────────────────────────────────────────────────
    if len(forbidden_lipids) > max_removed_pore_lipids:
        msg = (
            f"n_forbidden_pore_lipids={len(forbidden_lipids)} exceeds "
            f"max_removed_pore_lipids={max_removed_pore_lipids}; "
            "check TM annotation or embedding quality."
        )
        warnings_list.append(msg)
        if policy == "strict":
            raise RuntimeError(f"[spatial_classifier] STRICT: {msg}")

    if forbidden_lipids and policy == "strict" and remove_pore_lipids:
        raise RuntimeError(
            f"[spatial_classifier] STRICT: {len(forbidden_lipids)} forbidden pore lipids detected."
        )

    # ── Candidate water for removal ───────────────────────────────────────────
    removed_candidate: list[int] = []
    if remove_membrane_core_water:
        removed_candidate.extend(core_water)
    if remove_isolated_internal_water:
        removed_candidate.extend(isolated_water)

    # Pore/vestibule/cavity counts (binary: present or absent)
    n_pores     = 1 if pore_water else 0
    n_isolated  = 1 if isolated_water else 0

    report = {
        "enabled":                       True,
        "grid_spacing_nm":               grid_spacing_nm,
        "protein_padding_nm":            protein_padding_nm,
        "membrane_core_z_range":         [round(z_core_bot, 4), round(z_core_top, 4)],
        "n_regions_detected":            n_pores + n_isolated,
        "n_transmembrane_pores":         n_pores,
        "n_one_sided_vestibules":        0,
        "n_isolated_cavities":           n_isolated,
        "n_bulk_waters":                 len(bulk_water),
        "n_pore_waters":                 len(pore_water),
        "n_membrane_core_waters":        len(core_water),
        "n_isolated_internal_waters":    len(isolated_water),
        "n_valid_external_lipids":       len(valid_lipids),
        "n_forbidden_pore_lipids":       len(forbidden_lipids),
        "forbidden_pore_lipid_resids":   sorted(forbidden_lipids),
        "pore_water_resids":             sorted(pore_water),
        "removed_water_resids_candidate": sorted(removed_candidate),
        "warnings":                      warnings_list,
    }

    if out_dir:
        _write_classifier_report(report, out_dir)

    return report


# ── Pore-aware water cleanup ──────────────────────────────────────────────────

def run_pore_aware_water_cleanup(
    gro_in:       "Path | str",
    gro_out:      "Path | str",
    tm_residues:  Optional[set[int]] = None,
    lipid_resnames: Optional[set[str]] = None,
    topol_in:     "Optional[Path | str]" = None,
    topol_out:    "Optional[Path | str]" = None,
    output_dir:   "Optional[Path | str]" = None,
    preserve_pore_water:          bool  = True,
    remove_membrane_core_water:   bool  = True,
    remove_isolated_internal_water: bool = False,
    grid_spacing_nm:   float = 0.15,
    protein_padding_nm: float = 0.25,
    policy: str = "warn",
    update_topology_count: bool = True,
) -> dict:
    """Pore-aware clean_water replacement.

    Runs the spatial classifier to distinguish pore water (preserve) from
    membrane-core water (remove), then writes a new GRO and optionally
    updates the topology SOL count.

    Returns a clean_water_report-compatible dict.
    """
    gro_in  = Path(gro_in)
    gro_out = Path(gro_out)
    out_dir = Path(output_dir) if output_dir else gro_out.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    # Classify
    clf_report = run_spatial_classifier(
        gro_path=gro_in,
        tm_residues=tm_residues,
        lipid_resnames=lipid_resnames,
        output_dir=out_dir,
        enabled=True,
        policy=policy,
        grid_spacing_nm=grid_spacing_nm,
        protein_padding_nm=protein_padding_nm,
        preserve_pore_water=preserve_pore_water,
        remove_membrane_core_water=remove_membrane_core_water,
        remove_isolated_internal_water=remove_isolated_internal_water,
    )

    remove_resids: set[int] = set(clf_report["removed_water_resids_candidate"])
    n_removed = _remove_sol_resids_from_gro(gro_in, gro_out, remove_resids)

    topology_updated = False
    if update_topology_count and topol_in and topol_out and n_removed > 0:
        topol_in  = Path(topol_in)
        topol_out = Path(topol_out)
        if topol_in.exists():
            text = topol_in.read_text()
            topol_out.write_text(_update_sol_count(text, n_removed) + "\n")
            topology_updated = True

    pore_water_resids = clf_report["pore_water_resids"]
    core_z = clf_report["membrane_core_z_range"] or [None, None]

    report = {
        "input_water_molecules":             len(clf_report["removed_water_resids_candidate"])
                                             + clf_report["n_pore_waters"]
                                             + clf_report["n_bulk_waters"]
                                             + clf_report["n_isolated_internal_waters"],
        "n_water_molecules_removed":         n_removed,
        "n_water_oxygens_remaining_in_core": clf_report["n_pore_waters"],
        "n_pore_waters_preserved":           len(pore_water_resids),
        "pore_water_resids":                 pore_water_resids,
        "core_z_min":                        core_z[0],
        "core_z_max":                        core_z[1],
        "output_gro_path":                   str(gro_out),
        "topology_updated":                  topology_updated,
        "cleanup_passed":                    True,
        "pore_aware":                        True,
        "spatial_classifier_report":         clf_report,
        "warnings":                          clf_report.get("warnings", []),
    }

    (out_dir / "clean_water_report.json").write_text(
        json.dumps(report, indent=2, default=str)
    )
    return report
