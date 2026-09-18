"""
validators/membrane_protein_spatial_classifier.py
Phase 10B: Membrane-Protein Spatial Classifier.

`run_spatial_classifier` / `run_pore_aware_water_cleanup` (the public API of
this module) now delegate to validators.membrane_water — a real 3-D
solvent-topology engine (local leaflet surfaces, VDW excluded volume for
BOTH protein and lipid, connected-component region analysis, and
protein-vs-lipid wall-composition discrimination between a true
transmembrane pore and a lateral lipid-packing defect). See
validators/membrane_water/classify.py for the algorithm.

The private helpers below this docstring (_parse_gro_atoms,
_compute_membrane_core_z, _compute_tm_geometry, _build_occupancy_grid,
_dilate_step, _flood_fill, _voxel_idx, _update_sol_count,
_remove_sol_resids_from_gro) are the ORIGINAL (Phase 10B) global-Z-slab /
single-TM-circle implementation. They are kept, unmodified, only because
external code may still import them by name; run_spatial_classifier and
run_pore_aware_water_cleanup no longer call them. Do not extend this file's
own classification logic — extend validators/membrane_water instead.

Legacy classification categories, still reported for backward compatibility
(now computed for real by validators.membrane_water instead of the old
presence-flag / global-Z-slab approximations):
  bulk_water                – outside the local membrane core (leaflet model)
  pore_water                – protein-walled channel, connects both bulk sides
  membrane_core_water       – lipid-walled membrane-core water (→ remove)
  isolated_internal_water   – enclosed cavity, not connected to bulk
  valid_external_lipid      – lipid not intruding into a pore/vestibule
  forbidden_pore_lipid      – lipid COM inside a pore/vestibule region
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

from validators.membrane_water import classify_membrane_water as _classify_membrane_water
from validators.membrane_water import determine_removal_set as _determine_removal_set
from validators.membrane_water import run_cleanup as _run_membrane_water_cleanup
from validators.membrane_water.constants import CLEANUP_MODE_AGGRESSIVE, CLEANUP_MODE_CONSERVATIVE

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


# ── Main classifier (delegates to validators.membrane_water) ──────────────────

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
    lipid: str = "DPPC",
    forcefield: str = "opls-aa",
) -> dict:
    """Classify water and lipid molecules by spatial region in a
    protein-membrane system, using the validators.membrane_water engine
    (real leaflet-local membrane geometry, VDW excluded volume for both
    protein and lipid, connected-component pore/vestibule/cavity/defect
    classification — see validators/membrane_water/classify.py).

    Signature and legacy report fields are preserved for backward
    compatibility. `protein_padding_nm` is now interpreted as a
    solvent-probe-radius override (the primary excluded-volume definition is
    per-element van der Waals radii, not one global padding value) — pass
    the default (0.25) to get the engine's own default probe radius, or a
    different value to override it. `preserve_pore_water` and
    `remove_pore_lipids` are accepted for compatibility; pore/vestibule
    water and non-forbidden lipids are always preserved by the new engine
    regardless of these flags (there is no code path that removes them).
    Writes membrane_protein_spatial_classification_report.json if
    output_dir is given, extended with new fields (membrane_model,
    tm_annotation, regions, water_counts, cleanup_mode, confidence) on top
    of every legacy key.
    """
    gro_path = Path(gro_path)
    out_dir = Path(output_dir) if output_dir else None
    if out_dir:
        out_dir.mkdir(parents=True, exist_ok=True)

    if not enabled:
        rep = _disabled_report()
        if out_dir:
            _write_classifier_report(rep, out_dir)
        return rep

    if not _HAS_NP:
        rep = _disabled_report()
        rep["warnings"] = ["numpy not available — spatial classifier requires numpy"]
        if out_dir:
            _write_classifier_report(rep, out_dir)
        return rep

    cleanup_mode = CLEANUP_MODE_CONSERVATIVE if remove_membrane_core_water else CLEANUP_MODE_AGGRESSIVE
    # protein_padding_nm at its historical default (0.25) maps to the
    # engine's own default solvent-probe radius; any other value is treated
    # as an explicit override, preserving the old "bigger padding = bigger
    # excluded volume" knob without making it the primary VDW definition.
    probe_kwargs = {} if protein_padding_nm == 0.25 else {"solvent_probe_radius_nm": protein_padding_nm}

    try:
        report = _classify_membrane_water(
            gro_path,
            tm_residues=tm_residues,
            lipid_resnames=frozenset(lipid_resnames) if lipid_resnames else None,
            lipid=lipid,
            forcefield=forcefield,
            grid_spacing_nm=grid_spacing_nm,
            cleanup_mode=cleanup_mode,
            max_removed_pore_lipids=max_removed_pore_lipids,
            **probe_kwargs,
        )
    except Exception as exc:  # pragma: no cover - defensive fail-safe
        rep = _disabled_report()
        rep["warnings"] = [f"membrane_water engine failed: {exc}; falling back to disabled report"]
        if out_dir:
            _write_classifier_report(rep, out_dir)
        return rep

    removed_candidate = _determine_removal_set(
        report,
        cleanup_mode=report["cleanup_mode"],
        remove_membrane_core_water=remove_membrane_core_water,
        remove_isolated_internal_water=remove_isolated_internal_water,
    )

    if report["n_forbidden_pore_lipids"] > max_removed_pore_lipids and policy == "strict":
        raise RuntimeError(
            f"[spatial_classifier] STRICT: n_forbidden_pore_lipids="
            f"{report['n_forbidden_pore_lipids']} exceeds max_removed_pore_lipids="
            f"{max_removed_pore_lipids}; check TM annotation or embedding quality."
        )
    if report["forbidden_pore_lipid_resids"] and policy == "strict" and remove_pore_lipids:
        raise RuntimeError(
            f"[spatial_classifier] STRICT: {report['n_forbidden_pore_lipids']} forbidden "
            "pore lipids detected."
        )

    report["removed_water_resids_candidate"] = sorted(removed_candidate)

    if out_dir:
        _write_classifier_report(report, out_dir)
    return report


# ── Pore-aware water cleanup (delegates to validators.membrane_water) ─────────

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
    lipid: str = "DPPC",
    forcefield: str = "opls-aa",
) -> dict:
    """Pore-aware clean_water replacement, delegating to
    validators.membrane_water.run_cleanup() — see that module and
    validators/membrane_water/classify.py for the algorithm.

    Never deletes water solely because its Z coordinate falls inside the
    membrane core: only water confidently classified as a lipid/membrane
    core defect or a steric clash is removed by default (conservative mode);
    ambiguous water and any water inside a real pore/vestibule/isolated
    cavity is preserved.

    Returns a clean_water_report-compatible dict (every legacy key
    populated with real values, plus new fields under 'spatial_classifier_report').
    """
    gro_in  = Path(gro_in)
    gro_out = Path(gro_out)
    out_dir = Path(output_dir) if output_dir else gro_out.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    cleanup_mode = CLEANUP_MODE_CONSERVATIVE if remove_membrane_core_water else CLEANUP_MODE_AGGRESSIVE
    probe_kwargs = {} if protein_padding_nm == 0.25 else {"solvent_probe_radius_nm": protein_padding_nm}

    result = _run_membrane_water_cleanup(
        gro_in, gro_out,
        tm_residues=tm_residues,
        lipid_resnames=frozenset(lipid_resnames) if lipid_resnames else None,
        lipid=lipid,
        forcefield=forcefield,
        topol_in=topol_in if update_topology_count else None,
        topol_out=topol_out if update_topology_count else None,
        output_dir=out_dir,
        cleanup_mode=cleanup_mode,
        remove_membrane_core_water=remove_membrane_core_water,
        remove_isolated_internal_water=remove_isolated_internal_water,
        grid_spacing_nm=grid_spacing_nm,
        **probe_kwargs,
    )

    clf_report = result["classification"]
    if policy == "strict" and clf_report["forbidden_pore_lipid_resids"]:
        raise RuntimeError(
            f"[spatial_classifier] STRICT: {clf_report['n_forbidden_pore_lipids']} forbidden "
            "pore lipids detected."
        )

    core_z = clf_report["membrane_core_z_range"] or [None, None]
    report = {
        "input_water_molecules":             sum(len(v) for v in clf_report["water_resid_by_category"].values()),
        "n_water_molecules_removed":          result["n_water_molecules_removed"],
        "n_water_oxygens_remaining_in_core":  result["n_water_oxygens_remaining_in_core"],
        "n_pore_waters_preserved":            clf_report["n_pore_waters"],
        "pore_water_resids":                  clf_report["pore_water_resids"],
        "core_z_min":                         core_z[0],
        "core_z_max":                         core_z[1],
        "output_gro_path":                    result["output_gro_path"],
        "topology_updated":                   result["topology_updated"],
        "cleanup_passed":                     result["cleanup_passed"],
        "pore_aware":                         True,
        "spatial_classifier_report":          clf_report,
        "cleanup_mode":                       result["cleanup_mode"],
        "validation":                         result["validation"],
        "warnings":                           result["warnings"],
    }

    (out_dir / "clean_water_report.json").write_text(
        json.dumps(report, indent=2, default=str)
    )
    return report
