"""
validators/lipid_refill.py
Phase 10A: Surface-Guided Lipid Refill Actuator.

Inserts complete lipid molecules into annular TM-surface gaps identified by the
Phase 9A protein–membrane interface evaluator.

Workflow
--------
1. Evaluate current interface quality (Phase 9A).
2. If quality passes, exit with no changes.
3. For each gap cluster in the report, compute N target positions.
4. Extract one upper-leaflet and one lower-leaflet template residue from the
   staged bilayer GRO (membrane_assets/).
5. At each target site, try 6 Z-rotations of the appropriate template.
6. Accept the first rotation free of protein and lipid clashes.
7. Append accepted lipids to the system GRO and update the atom count.
8. Re-evaluate; restore baseline if refill worsened coverage.
9. Write interface_lipid_refill_report.json.
"""
from __future__ import annotations

import json
import math
import random
import shutil
from pathlib import Path
from typing import Optional

try:
    import numpy as np
    _HAS_NP = True
except ImportError:
    _HAS_NP = False

DEFAULT_LIPID_RESNAMES: frozenset[str] = frozenset({
    "DPP", "DPPC", "POPC", "POPE", "POPG", "POPS", "CHOL",
    "PALM", "OLEO", "PALC", "SM", "CER",
})
_SOLVENT_RESNAMES: frozenset[str] = frozenset({
    "SOL", "HOH", "WAT", "TIP3", "TIP4", "TIP5",
    "NA", "CL", "SOD", "MG", "K", "CA",
})
# Preferred atom names for headgroup anchor detection, in priority order
_HEADGROUP_PREF: tuple[str, ...] = ("P", "P1", "N", "NZ", "O33", "O11")
_HEADGROUP_PREF_SET: frozenset[str] = frozenset(_HEADGROUP_PREF)
# Z-rotation candidates tried per placement site
_ROTATIONS_DEG: tuple[float, ...] = (0.0, 60.0, 120.0, 180.0, 240.0, 300.0)
# Fan-spread angles (degrees from outward direction) for multi-lipid cluster placements
_SPREAD_ANGLES: dict[int, list[float]] = {
    1: [0.0],
    2: [-30.0, 30.0],
    3: [-45.0, 0.0, 45.0],
    4: [-60.0, -20.0, 20.0, 60.0],
}


# ── GRO I/O ──────────────────────────────────────────────────────────────────

def _parse_gro(path: Path) -> tuple[list, tuple, str]:
    """Parse GRO file.

    Returns (atoms, box, title).
    Each atom dict: resid, resname, atomname, x, y, z, raw (original line).
    """
    lines = Path(path).read_text().splitlines()
    if len(lines) < 3:
        raise ValueError(f"GRO file too short: {path}")
    title = lines[0]
    n = int(lines[1].strip())
    atoms = []
    for ln in lines[2: 2 + n]:
        if len(ln) < 44:
            atoms.append({"resid": 0, "resname": "", "atomname": "",
                          "x": 0.0, "y": 0.0, "z": 0.0, "raw": ln})
            continue
        try:
            atoms.append({
                "resid":    int(ln[0:5]),
                "resname":  ln[5:10].strip(),
                "atomname": ln[10:15].strip(),
                "x":        float(ln[20:28]),
                "y":        float(ln[28:36]),
                "z":        float(ln[36:44]),
                "raw":      ln,
            })
        except ValueError:
            atoms.append({"resid": 0, "resname": "", "atomname": "",
                          "x": 0.0, "y": 0.0, "z": 0.0, "raw": ln})
    box_line = lines[2 + n] if len(lines) > 2 + n else "10.0 10.0 10.0"
    parts = box_line.split()
    box = (float(parts[0]), float(parts[1]), float(parts[2]))
    return atoms, box, title


def _write_gro(path: Path, title: str, atoms: list, box: tuple) -> None:
    """Write GRO.

    Atoms with 'placed': True → construct line from (resid, resname, atomname, x, y, z).
    Other atoms → keep raw line but update atom-index field (cols 15-19).
    """
    lines = [title, str(len(atoms))]
    for i, a in enumerate(atoms, 1):
        atom_idx = i % 100000  # GRO atom numbers wrap at 99999
        if a.get("placed"):
            resid   = a["resid"] % 100000
            lines.append(
                f"{resid:5d}{a['resname']:<5s}{a['atomname']:>5s}"
                f"{atom_idx:5d}{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}"
            )
        else:
            raw = a["raw"]
            lines.append(raw[:15] + f"{atom_idx:5d}" + raw[20:])
    bx, by, bz = box
    lines.append(f"   {bx:.5f}   {by:.5f}   {bz:.5f}")
    path.write_text("\n".join(lines) + "\n")


# ── Bilayer geometry ──────────────────────────────────────────────────────────

def _find_midplane_z(atoms: list, lipid_resnames: frozenset) -> float:
    """Mean Z of all lipid-residue atoms."""
    zs = [a["z"] for a in atoms if a["resname"] in lipid_resnames]
    return sum(zs) / len(zs) if zs else 0.0


def _compute_leaflet_zrefs(
    atoms: list,
    lipid_resnames: frozenset,
    midplane_z: float,
) -> tuple[float, float]:
    """Return (upper_hg_z, lower_hg_z): mean Z of headgroup atoms per leaflet.

    Headgroup atoms are identified by names in _HEADGROUP_PREF_SET.
    Falls back to the extreme-Z lipid atom in each leaflet.
    """
    upper_hg: list[float] = []
    lower_hg: list[float] = []
    for a in atoms:
        if a["resname"] not in lipid_resnames:
            continue
        if a["atomname"] in _HEADGROUP_PREF_SET:
            if a["z"] >= midplane_z:
                upper_hg.append(a["z"])
            else:
                lower_hg.append(a["z"])

    if not upper_hg:
        # Fallback: max Z among all upper-leaflet lipid atoms
        upper_all = [a["z"] for a in atoms
                     if a["resname"] in lipid_resnames and a["z"] >= midplane_z]
        upper_hg = [max(upper_all)] if upper_all else [midplane_z + 1.0]
    if not lower_hg:
        lower_all = [a["z"] for a in atoms
                     if a["resname"] in lipid_resnames and a["z"] < midplane_z]
        lower_hg = [min(lower_all)] if lower_all else [midplane_z - 1.0]

    return float(sum(upper_hg) / len(upper_hg)), float(sum(lower_hg) / len(lower_hg))


def _compute_tm_center_xy(atoms: list, tm_residues: set[int]) -> tuple[float, float]:
    """Mean XY of all TM-residue non-hydrogen atoms."""
    tm_atoms = [
        a for a in atoms
        if a["resid"] in tm_residues
        and a["resname"] not in DEFAULT_LIPID_RESNAMES
        and a["resname"] not in _SOLVENT_RESNAMES
        and not a["atomname"].upper().startswith("H")
    ]
    if not tm_atoms:
        xs = [a["x"] for a in atoms]
        ys = [a["y"] for a in atoms]
        return (sum(xs) / len(xs) if xs else 0.0, sum(ys) / len(ys) if ys else 0.0)
    return (
        sum(a["x"] for a in tm_atoms) / len(tm_atoms),
        sum(a["y"] for a in tm_atoms) / len(tm_atoms),
    )


# ── Template extraction ───────────────────────────────────────────────────────

def _group_by_resid(atoms: list, lipid_resname: str) -> dict[int, list]:
    """Group atoms by residue ID for the given lipid resname."""
    from collections import defaultdict
    groups: dict[int, list] = defaultdict(list)
    for a in atoms:
        if a["resname"] == lipid_resname:
            groups[a["resid"]].append(a)
    return dict(groups)


def _find_anchor_idx(template_atoms: list, prefer_max_z: bool) -> int:
    """Return index of headgroup anchor.

    Prefers atoms named in _HEADGROUP_PREF (first match).
    Falls back to the atom with extreme Z.
    """
    for pref in _HEADGROUP_PREF:
        for i, a in enumerate(template_atoms):
            if a["atomname"] == pref:
                return i
    if prefer_max_z:
        return max(range(len(template_atoms)), key=lambda i: template_atoms[i]["z"])
    return min(range(len(template_atoms)), key=lambda i: template_atoms[i]["z"])


def _extract_template(
    template_gro: Path,
    lipid_resname: str,
) -> tuple[list, list, float]:
    """Extract upper- and lower-leaflet lipid templates from a bilayer GRO.

    Returns:
        upper_rel  — list of {"atomname", "dx", "dy", "dz"} relative to headgroup anchor
        lower_rel  — same for lower leaflet
        midplane_z — bilayer midplane Z of the template file
    """
    from collections import Counter

    atoms, box, _ = _parse_gro(template_gro)
    groups = _group_by_resid(atoms, lipid_resname)

    if not groups:
        raise ValueError(
            f"No lipid residues with resname '{lipid_resname}' found in {template_gro}"
        )

    # Expected atom count = mode across all residues
    counts = Counter(len(v) for v in groups.values())
    expected_n = counts.most_common(1)[0][0]

    complete = {rid: v for rid, v in groups.items() if len(v) == expected_n}
    if not complete:
        raise ValueError(
            f"No complete {lipid_resname} residues (n={expected_n} atoms) in {template_gro}"
        )

    # Bilayer midplane Z
    all_z = [a["z"] for res_atoms in complete.values() for a in res_atoms]
    midplane_z = sum(all_z) / len(all_z)

    # Box centre XY for "representative" selection
    box_cx, box_cy = box[0] / 2.0, box[1] / 2.0

    upper_cands: list[tuple[float, int]] = []
    lower_cands: list[tuple[float, int]] = []
    for rid, res_atoms in complete.items():
        mean_z = sum(a["z"] for a in res_atoms) / len(res_atoms)
        cx = sum(a["x"] for a in res_atoms) / len(res_atoms)
        cy = sum(a["y"] for a in res_atoms) / len(res_atoms)
        dist = math.sqrt((cx - box_cx) ** 2 + (cy - box_cy) ** 2)
        if mean_z >= midplane_z:
            upper_cands.append((dist, rid))
        else:
            lower_cands.append((dist, rid))

    if not upper_cands:
        raise ValueError(f"No upper-leaflet {lipid_resname} residues found in {template_gro}")
    if not lower_cands:
        raise ValueError(f"No lower-leaflet {lipid_resname} residues found in {template_gro}")

    upper_cands.sort()
    lower_cands.sort()

    upper_atoms = complete[upper_cands[0][1]]
    lower_atoms = complete[lower_cands[0][1]]

    def _to_relative(res_atoms: list, prefer_max_z: bool) -> list:
        ai = _find_anchor_idx(res_atoms, prefer_max_z)
        ax, ay, az = res_atoms[ai]["x"], res_atoms[ai]["y"], res_atoms[ai]["z"]
        return [
            {"atomname": a["atomname"],
             "dx": a["x"] - ax,
             "dy": a["y"] - ay,
             "dz": a["z"] - az}
            for a in res_atoms
        ]

    upper_rel = _to_relative(upper_atoms, prefer_max_z=True)
    lower_rel = _to_relative(lower_atoms, prefer_max_z=False)
    return upper_rel, lower_rel, midplane_z


# ── Placement ─────────────────────────────────────────────────────────────────

def _rotate_xy(dx: float, dy: float, angle_deg: float) -> tuple[float, float]:
    """Rotate vector (dx, dy) around Z by angle_deg degrees."""
    θ = math.radians(angle_deg)
    c, s = math.cos(θ), math.sin(θ)
    return c * dx - s * dy, s * dx + c * dy


def _place_template(
    rel_atoms:    list,       # [{"atomname", "dx", "dy", "dz"}, ...]
    target_x:     float,
    target_y:     float,
    target_z:     float,      # headgroup anchor Z in the system
    rotation_deg: float,
    lipid_resname: str,
    resid:        int,
) -> list:
    """Place template atoms at target position with given Z-axis rotation.

    Returns list of atom dicts with placed=True.
    """
    placed = []
    for ra in rel_atoms:
        rx, ry = _rotate_xy(ra["dx"], ra["dy"], rotation_deg)
        placed.append({
            "resid":    resid,
            "resname":  lipid_resname,
            "atomname": ra["atomname"],
            "x":        target_x + rx,
            "y":        target_y + ry,
            "z":        target_z + ra["dz"],
            "placed":   True,
            "raw":      "",
        })
    return placed


def _has_clash(candidate: list, ref_xyz, cutoff_nm: float) -> bool:
    """True if any candidate atom is closer than cutoff_nm to any ref_xyz point."""
    if not _HAS_NP or ref_xyz is None or len(ref_xyz) == 0:
        return False
    import numpy as np
    c_xyz = np.array(
        [[a["x"], a["y"], a["z"]] for a in candidate
         if not a.get("atomname", "").upper().startswith("H")],
        dtype=np.float32,
    )
    if len(c_xyz) == 0:
        return False
    batch = 512
    for i in range(0, len(c_xyz), batch):
        cb = c_xyz[i: i + batch]
        diff = cb[:, None, :] - ref_xyz[None, :, :]
        if float(np.sqrt((diff ** 2).sum(axis=2)).min()) < cutoff_nm:
            return True
    return False


# ── Target site generation ────────────────────────────────────────────────────

def _n_lipids_for_cluster(cluster_size: int, max_per_cluster: int) -> int:
    """Estimate how many lipids to insert for a gap cluster of given atom-size."""
    if cluster_size <= 3:
        n = 1
    elif cluster_size <= 7:
        n = 2
    elif cluster_size <= 13:
        n = 3
    else:
        n = 4
    return min(n, max_per_cluster)


def _compute_target_positions(
    centroid_nm: list,               # [x, y, z]
    tm_center_xy: tuple,             # (x_tm, y_tm)
    target_contact_distance_nm: float,
    n_lipids: int,
    upper_hg_z: float,
    lower_hg_z: float,
    midplane_z: float,
) -> list:
    """Return list of (x, y, z, leaflet_str) placement targets.

    Distributes n_lipids in a fan around the outward direction from the TM
    bundle centre toward the cluster centroid.
    """
    cx, cy, cz = centroid_nm
    x_tm, y_tm = tm_center_xy

    # Outward direction
    dx, dy = cx - x_tm, cy - y_tm
    mag = math.sqrt(dx * dx + dy * dy)
    if mag < 1e-6:
        dx, dy = 1.0, 0.0
    else:
        dx, dy = dx / mag, dy / mag

    outward_angle = math.atan2(dy, dx)
    leaflet = "upper" if cz >= midplane_z else "lower"
    target_z = upper_hg_z if leaflet == "upper" else lower_hg_z

    spread = _SPREAD_ANGLES.get(n_lipids,
                                [i * 360.0 / n_lipids for i in range(n_lipids)])
    positions = []
    for angle_offset in spread:
        total = outward_angle + math.radians(angle_offset)
        tx = cx + math.cos(total) * target_contact_distance_nm
        ty = cy + math.sin(total) * target_contact_distance_nm
        positions.append((tx, ty, target_z, leaflet))
    return positions


# ── Phase 9A interface evaluation helper ─────────────────────────────────────

def _evaluate_interface(
    gro_path:              Path,
    tm_residues:           set[int],
    lip_set:               frozenset,
    min_fraction_covered:  float,
    max_fraction_exposed_gap: float,
    max_p90_distance_nm:   float,
    output_dir:            Path,
    simforge_root:         Optional[str],
) -> dict:
    try:
        if simforge_root:
            import sys
            if simforge_root not in sys.path:
                sys.path.insert(0, simforge_root)
        from validators.protein_membrane_interface import evaluate_protein_membrane_interface
        output_dir.mkdir(parents=True, exist_ok=True)
        return evaluate_protein_membrane_interface(
            gro_path=gro_path,
            tm_residues=tm_residues,
            lipid_resnames=set(lip_set),
            output_dir=output_dir,
            min_fraction_covered=min_fraction_covered,
            max_fraction_exposed_gap=max_fraction_exposed_gap,
            max_p90_distance_nm=max_p90_distance_nm,
        )
    except Exception as exc:
        return {
            "error": str(exc),
            "quality_passed": None,
            "fraction_covered": None,
            "fraction_exposed_gap": None,
            "p90_nearest_lipid_distance": None,
            "gap_clusters": [],
            "warnings": [str(exc)],
        }


# ── Report helpers ────────────────────────────────────────────────────────────

def _disabled_report(
    policy: str,
    reason: str,
    input_gro: str,
    output_gro: str,
    template_gro: str,
) -> dict:
    return {
        "enabled":                    False,
        "policy":                     policy,
        "input_system_gro":           input_gro,
        "output_system_gro":          output_gro,
        "source_template_gro":        template_gro,
        "n_gap_clusters":             0,
        "n_attempted_sites":          0,
        "n_inserted_lipids":          0,
        "n_rejected_sites":           0,
        "inserted_resids":            [],
        "inserted_leaflets":          [],
        "inserted_target_clusters":   [],
        "rejected_reasons":           [],
        "pre_refill_interface_metrics": {},
        "post_refill_interface_metrics": {},
        "fraction_covered_before":    None,
        "fraction_covered_after":     None,
        "fraction_exposed_gap_before": None,
        "fraction_exposed_gap_after": None,
        "p90_distance_before":        None,
        "p90_distance_after":         None,
        "coordinate_file_modified":   False,
        "topology_modified":          False,
        "warning_messages":           [f"Lipid refill disabled: {reason}"],
    }


def _noop_report(
    policy: str,
    input_gro: str,
    output_gro: str,
    template_gro: str,
    pre_report: dict,
    fc: float,
    fg: float,
    p90: float,
) -> dict:
    return {
        "enabled":                    True,
        "policy":                     policy,
        "input_system_gro":           input_gro,
        "output_system_gro":          output_gro,
        "source_template_gro":        template_gro,
        "n_gap_clusters":             len(pre_report.get("gap_clusters") or []),
        "n_attempted_sites":          0,
        "n_inserted_lipids":          0,
        "n_rejected_sites":           0,
        "inserted_resids":            [],
        "inserted_leaflets":          [],
        "inserted_target_clusters":   [],
        "rejected_reasons":           [],
        "pre_refill_interface_metrics": pre_report,
        "post_refill_interface_metrics": {},
        "fraction_covered_before":    round(fc, 4),
        "fraction_covered_after":     round(fc, 4),
        "fraction_exposed_gap_before": round(fg, 4),
        "fraction_exposed_gap_after": round(fg, 4),
        "p90_distance_before":        round(p90, 4),
        "p90_distance_after":         round(p90, 4),
        "coordinate_file_modified":   False,
        "topology_modified":          False,
        "warning_messages":           [],
    }


def _write_report(report: dict, output_dir: Path) -> None:
    (output_dir / "interface_lipid_refill_report.json").write_text(
        json.dumps(report, indent=2, default=str)
    )


# ── Main function ─────────────────────────────────────────────────────────────

def run_lipid_refill(
    system_gro:                  Path | str,
    template_gro:                Path | str,
    tm_residues:                 set[int],
    output_gro:                  Path | str,
    output_dir:                  Path | str,
    lipid_resname:               str   = "DPP",
    enabled:                     bool  = True,
    policy:                      str   = "warn",
    max_inserted_lipids:         int   = 40,
    max_lipids_per_gap_cluster:  int   = 4,
    target_contact_distance_nm:  float = 0.45,
    protein_clash_cutoff_nm:     float = 0.20,
    lipid_clash_cutoff_nm:       float = 0.18,
    deterministic_seed:          int   = 17,
    min_fraction_covered:        float = 0.75,
    max_fraction_exposed_gap:    float = 0.15,
    max_p90_distance_nm:         float = 0.70,
    simforge_root:               Optional[str] = None,
    spatial_classifier_report:   Optional[dict] = None,
) -> dict:
    """Run Phase 10A surface-guided lipid refill.

    Reads system_gro, inserts lipid templates into detected annular gaps,
    and writes the result to output_gro (may be the same path for in-place update).

    Returns the full refill report dict.
    """
    system_gro   = Path(system_gro)
    template_gro = Path(template_gro)
    output_gro   = Path(output_gro)
    output_dir   = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    lip_set  = DEFAULT_LIPID_RESNAMES
    warnings: list[str] = []

    # ── Disabled fast path ────────────────────────────────────────────────────
    if not enabled:
        r = _disabled_report(policy, "disabled",
                             str(system_gro), str(output_gro), str(template_gro))
        _write_report(r, output_dir)
        if not output_gro.exists() or output_gro.resolve() != system_gro.resolve():
            shutil.copy(system_gro, output_gro)
        return r

    if not tm_residues:
        msg = "tm_residues not provided — lipid refill requires TM annotation"
        warnings.append(msg)
        r = _disabled_report(policy, "no_tm_residues",
                             str(system_gro), str(output_gro), str(template_gro))
        r["warning_messages"] = [msg]
        _write_report(r, output_dir)
        if output_gro.resolve() != system_gro.resolve():
            shutil.copy(system_gro, output_gro)
        return r

    if not template_gro.exists():
        msg = f"Bilayer template not found: {template_gro}"
        warnings.append(msg)
        r = _disabled_report(policy, "template_not_found",
                             str(system_gro), str(output_gro), str(template_gro))
        r["warning_messages"] = [msg]
        _write_report(r, output_dir)
        if output_gro.resolve() != system_gro.resolve():
            shutil.copy(system_gro, output_gro)
        if policy == "strict":
            raise RuntimeError(f"[lipid_refill] STRICT: {msg}")
        return r

    # ── Pre-refill evaluation ─────────────────────────────────────────────────
    pre_iface_dir = output_dir / "pre_refill_iface"
    pre_iface_dir.mkdir(exist_ok=True)
    pre_report = _evaluate_interface(
        system_gro, tm_residues, lip_set,
        min_fraction_covered, max_fraction_exposed_gap, max_p90_distance_nm,
        pre_iface_dir, simforge_root,
    )

    fc_before  = float(pre_report.get("fraction_covered")  or 0.0)
    fg_before  = float(pre_report.get("fraction_exposed_gap") or 1.0)
    p90_before = float(pre_report.get("p90_nearest_lipid_distance") or 1.5)
    quality_before = bool(pre_report.get("quality_passed", False))

    # No work needed if already passing
    if quality_before:
        r = _noop_report(policy, str(system_gro), str(output_gro), str(template_gro),
                         pre_report, fc_before, fg_before, p90_before)
        _write_report(r, output_dir)
        if output_gro.resolve() != system_gro.resolve():
            shutil.copy(system_gro, output_gro)
        return r

    gap_clusters = pre_report.get("gap_clusters") or []
    if not gap_clusters:
        msg = "Interface quality failed but gap_clusters is empty — nothing to refill"
        warnings.append(msg)
        r = _noop_report(policy, str(system_gro), str(output_gro), str(template_gro),
                         pre_report, fc_before, fg_before, p90_before)
        r["warning_messages"] = [msg]
        _write_report(r, output_dir)
        if output_gro.resolve() != system_gro.resolve():
            shutil.copy(system_gro, output_gro)
        if policy == "strict":
            raise RuntimeError(f"[lipid_refill] STRICT: interface quality failed. {msg}")
        return r

    # ── Baseline backup ───────────────────────────────────────────────────────
    baseline_backup = output_dir / "system_baseline_before_refill.gro"
    shutil.copy(system_gro, baseline_backup)

    # ── Parse system ──────────────────────────────────────────────────────────
    sys_atoms, box, title = _parse_gro(system_gro)

    midplane_z = _find_midplane_z(sys_atoms, lip_set)
    upper_hg_z, lower_hg_z = _compute_leaflet_zrefs(sys_atoms, lip_set, midplane_z)
    tm_center_xy = _compute_tm_center_xy(sys_atoms, tm_residues)

    # ── Extract templates ─────────────────────────────────────────────────────
    try:
        upper_rel, lower_rel, _tmpl_mid = _extract_template(template_gro, lipid_resname)
    except Exception as exc:
        msg = f"Template extraction failed: {exc}"
        warnings.append(msg)
        r = _noop_report(policy, str(system_gro), str(output_gro), str(template_gro),
                         pre_report, fc_before, fg_before, p90_before)
        r["warning_messages"] = [msg]
        _write_report(r, output_dir)
        if output_gro.resolve() != system_gro.resolve():
            shutil.copy(system_gro, output_gro)
        if policy == "strict":
            raise RuntimeError(f"[lipid_refill] STRICT: {msg}")
        return r

    # ── Coordinate arrays for clash detection ─────────────────────────────────
    import numpy as np

    prot_atoms = [
        a for a in sys_atoms
        if a["resname"] not in lip_set
        and a["resname"] not in _SOLVENT_RESNAMES
        and not a["atomname"].upper().startswith("H")
    ]
    lip_atoms_orig = [
        a for a in sys_atoms
        if a["resname"] in lip_set
        and not a["atomname"].upper().startswith("H")
    ]

    prot_xyz = (
        np.array([[a["x"], a["y"], a["z"]] for a in prot_atoms], dtype=np.float32)
        if prot_atoms else np.zeros((0, 3), dtype=np.float32)
    )
    orig_lip_xyz = (
        np.array([[a["x"], a["y"], a["z"]] for a in lip_atoms_orig], dtype=np.float32)
        if lip_atoms_orig else np.zeros((0, 3), dtype=np.float32)
    )

    # ── Starting residue ID for inserted lipids ───────────────────────────────
    max_resid = max((a["resid"] for a in sys_atoms), default=0)
    next_resid = max_resid + 1

    # ── Deterministic rotation order per site ─────────────────────────────────
    rng = random.Random(deterministic_seed)

    # ── Insertion loop ────────────────────────────────────────────────────────
    inserted_atoms:          list = []
    inserted_resids:         list[int] = []
    inserted_leaflets:       list[str] = []
    inserted_target_clusters: list[int] = []
    rejected_reasons:        list[str] = []
    n_attempted  = 0
    n_inserted   = 0
    n_rejected   = 0
    new_lip_coords: list[list[float]] = []

    # ── Spatial classifier: TM footprint from classifier_report ──────────────
    _clf_tm_cx: Optional[float] = None
    _clf_tm_cy: Optional[float] = None
    _clf_tm_thresh: Optional[float] = None
    if spatial_classifier_report and spatial_classifier_report.get("enabled"):
        # Extract TM geometry from the report if available
        _clf_tm_cx, _clf_tm_cy = tm_center_xy
        # Use forbidden_pore_lipid_resids as proxy; classifier sets tm_thresh internally.
        # We use the membrane_core_z_range to skip targets at pore Z.
        pass

    for cl_idx, cluster in enumerate(gap_clusters):
        if n_inserted >= max_inserted_lipids:
            break

        centroid = cluster.get("centroid_nm", [0.0, 0.0, 0.0])

        # Skip cluster if its centroid is inside the TM bundle (pore region)
        if spatial_classifier_report and spatial_classifier_report.get("enabled"):
            _core_z = spatial_classifier_report.get("membrane_core_z_range")
            if _core_z and len(_core_z) == 2:
                c_x, c_y, c_z = centroid[0], centroid[1], centroid[2]
                z_bot, z_top = _core_z
                if z_bot <= c_z <= z_top:
                    # Check XY distance from TM centroid
                    _tcx, _tcy = tm_center_xy
                    _clf_dist = math.sqrt((c_x - _tcx) ** 2 + (c_y - _tcy) ** 2)
                    # Forbidden pore lipid resids indicate the TM footprint radius
                    _forbidden = spatial_classifier_report.get("forbidden_pore_lipid_resids", [])
                    if _forbidden:
                        rejected_reasons.append(
                            f"cluster {cl_idx}: centroid in pore region — skipped by "
                            "spatial_classifier_report"
                        )
                        n_rejected += 1
                        continue

        cl_size  = cluster.get("size", 1)
        n_for_cl = _n_lipids_for_cluster(cl_size, max_lipids_per_gap_cluster)

        targets = _compute_target_positions(
            centroid, tm_center_xy, target_contact_distance_nm,
            n_for_cl, upper_hg_z, lower_hg_z, midplane_z,
        )

        for t_x, t_y, t_z, leaflet in targets:
            if n_inserted >= max_inserted_lipids:
                break

            n_attempted += 1
            tmpl = upper_rel if leaflet == "upper" else lower_rel
            resid = next_resid + n_inserted

            # Try rotations in deterministic shuffled order
            rotations = list(_ROTATIONS_DEG)
            rng.shuffle(rotations)

            placed_ok = None
            for rot in rotations:
                candidate = _place_template(tmpl, t_x, t_y, t_z, rot,
                                            lipid_resname, resid)

                if _has_clash(candidate, prot_xyz, protein_clash_cutoff_nm):
                    continue
                if _has_clash(candidate, orig_lip_xyz, lipid_clash_cutoff_nm):
                    continue
                if new_lip_coords:
                    new_arr = np.array(new_lip_coords, dtype=np.float32)
                    if _has_clash(candidate, new_arr, lipid_clash_cutoff_nm):
                        continue

                placed_ok = candidate
                break

            if placed_ok is not None:
                inserted_atoms.extend(placed_ok)
                inserted_resids.append(resid)
                inserted_leaflets.append(leaflet)
                inserted_target_clusters.append(cl_idx)
                new_lip_coords.extend(
                    [a["x"], a["y"], a["z"]] for a in placed_ok
                    if not a.get("atomname", "").upper().startswith("H")
                )
                n_inserted += 1
            else:
                n_rejected += 1
                rejected_reasons.append(
                    f"cluster {cl_idx}: all {len(rotations)} rotations clash "
                    f"at site ({t_x:.3f}, {t_y:.3f})"
                )

    # ── Write refilled GRO ────────────────────────────────────────────────────
    if n_inserted > 0:
        all_atoms = sys_atoms + inserted_atoms
        _write_gro(output_gro, title, all_atoms, box)
    else:
        msg = "No lipids inserted — all candidate positions clashed"
        warnings.append(msg)
        r = _noop_report(policy, str(system_gro), str(output_gro), str(template_gro),
                         pre_report, fc_before, fg_before, p90_before)
        r["n_attempted_sites"] = n_attempted
        r["n_rejected_sites"]  = n_rejected
        r["rejected_reasons"]  = rejected_reasons
        r["warning_messages"]  = [msg]
        _write_report(r, output_dir)
        if output_gro.resolve() != system_gro.resolve():
            shutil.copy(system_gro, output_gro)
        if policy == "strict":
            raise RuntimeError(f"[lipid_refill] STRICT: {msg}")
        return r

    # ── Post-refill evaluation ────────────────────────────────────────────────
    post_iface_dir = output_dir / "post_refill_iface"
    post_iface_dir.mkdir(exist_ok=True)
    post_report = _evaluate_interface(
        output_gro, tm_residues, lip_set,
        min_fraction_covered, max_fraction_exposed_gap, max_p90_distance_nm,
        post_iface_dir, simforge_root,
    )
    fc_after  = float(post_report.get("fraction_covered")  or 0.0)
    fg_after  = float(post_report.get("fraction_exposed_gap") or 1.0)
    p90_after = float(post_report.get("p90_nearest_lipid_distance") or 1.5)
    quality_after = bool(post_report.get("quality_passed", False))

    # ── Restore baseline if refill worsened coverage ──────────────────────────
    coordinate_file_modified = True
    if fc_after < fc_before - 0.02:
        shutil.copy(baseline_backup, output_gro)
        coordinate_file_modified = False
        warnings.append(
            f"Refill worsened fraction_covered ({fc_before:.3f} → {fc_after:.3f}); "
            "baseline restored."
        )
        n_inserted = 0
        inserted_resids = []
        inserted_leaflets = []
        inserted_target_clusters = []

    report = {
        "enabled":                    True,
        "policy":                     policy,
        "input_system_gro":           str(system_gro),
        "output_system_gro":          str(output_gro),
        "source_template_gro":        str(template_gro),
        "n_gap_clusters":             len(gap_clusters),
        "n_attempted_sites":          n_attempted,
        "n_inserted_lipids":          n_inserted,
        "n_rejected_sites":           n_rejected,
        "inserted_resids":            inserted_resids,
        "inserted_leaflets":          inserted_leaflets,
        "inserted_target_clusters":   inserted_target_clusters,
        "rejected_reasons":           rejected_reasons,
        "pre_refill_interface_metrics":  pre_report,
        "post_refill_interface_metrics": post_report if coordinate_file_modified else {},
        "fraction_covered_before":    round(fc_before, 4),
        "fraction_covered_after":     round(fc_after  if coordinate_file_modified else fc_before, 4),
        "fraction_exposed_gap_before": round(fg_before, 4),
        "fraction_exposed_gap_after": round(fg_after  if coordinate_file_modified else fg_before, 4),
        "p90_distance_before":        round(p90_before, 4),
        "p90_distance_after":         round(p90_after  if coordinate_file_modified else p90_before, 4),
        "coordinate_file_modified":   coordinate_file_modified,
        "topology_modified":          False,
        "warning_messages":           warnings,
    }

    _write_report(report, output_dir)

    if policy == "strict" and not quality_after:
        raise RuntimeError(
            f"[lipid_refill] STRICT: interface quality still below thresholds after refill "
            f"(fraction_covered={fc_after:.3f}). "
            "Review interface_lipid_refill_report.json."
        )

    return report
