"""
validators/tm_exclusion_mask.py

TM-aware exclusion mask prototype for protein-membrane setups.
Reconstructs a TM-only excluded volume using TM annotations and classifies
lipids without modifying coordinates.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Optional, Set, List, Dict, Any


def get_angular_gap(lx: float, ly: float, lz: float, tm_ca_coords: list[tuple[float, float, float]]) -> float:
    """
    Calculate the maximum angular gap in the XY plane of surrounding TM atoms
    relative to a lipid's center.
    """
    local_cas = [c for c in tm_ca_coords if abs(c[2] - lz) < 1.5]
    if len(local_cas) < 5:
        local_cas = tm_ca_coords

    angles = []
    for c in local_cas:
        dx = c[0] - lx
        dy = c[1] - ly
        angles.append(math.atan2(dy, dx))

    if not angles:
        return 360.0

    angles = sorted(angles)
    max_gap = 0.0
    for i in range(len(angles)):
        if i == len(angles) - 1:
            gap = (angles[0] + 2 * math.pi) - angles[i]
        else:
            gap = angles[i+1] - angles[i]
        if gap > max_gap:
            max_gap = gap

    return math.degrees(max_gap)


def classify_lipids_tm_aware(
    gro_path: Path | str,
    tm_residues: Optional[Set[int] | List[int]],
    lipid_resname: str,
    grid_spacing_nm: float = 0.1,
    mask_padding_nm: float = 0.4,
) -> dict[str, Any]:
    """
    Analyze lipid coordinates relative to the TM domain of the protein.
    Classifies each lipid residue and writes reports.
    """
    if not tm_residues:
        raise ValueError(
            "TM annotation is empty or not provided. TM-aware exclusion mask requires "
            "valid transmembrane residues. Fallback recommendation: use the legacy "
            "inflategro backend by setting membrane.embedding.backend to 'inflategro'."
        )

    tm_residues_set = set(tm_residues)
    path = Path(gro_path)
    if not path.exists():
        raise FileNotFoundError(f"GRO file not found: {path}")

    lines = path.read_text().splitlines()
    if len(lines) < 3:
        raise ValueError(f"Invalid GRO file format: {path}")

    atom_lines = lines[2:-1]

    tm_coords: list[tuple[float, float, float]] = []
    tm_ca_coords: list[tuple[float, float, float]] = []
    soluble_coords: list[tuple[float, float, float]] = []
    lipid_atoms: dict[int, list[tuple[float, float, float]]] = {}

    solvent_names = {"SOL", "HOH", "NA", "CL", "K", "MG", "CA"}

    # Parse GROMACS GRO file
    for line in atom_lines:
        if len(line) < 44:
            continue
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        try:
            resnum = int(line[0:5])
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
        except ValueError:
            continue

        if resname == lipid_resname:
            lipid_atoms.setdefault(resnum, []).append((x, y, z))
        elif resname not in solvent_names:
            # Protein atom
            if resnum in tm_residues_set:
                tm_coords.append((x, y, z))
                if atomname == "CA":
                    tm_ca_coords.append((x, y, z))
            else:
                soluble_coords.append((x, y, z))

    if not tm_coords:
        raise ValueError(
            "No protein atoms matching the TM residue annotation were found in the GRO file. "
            "TM-aware exclusion mask requires valid transmembrane residues present in the system. "
            "Fallback recommendation: use the legacy inflategro backend."
        )

    # Fallback to all TM atoms if CA atoms are not present
    if not tm_ca_coords:
        tm_ca_coords = tm_coords

    # TM geometries
    tm_z_min = min(c[2] for c in tm_coords)
    tm_z_max = max(c[2] for c in tm_coords)
    tm_cx = sum(c[0] for c in tm_coords) / len(tm_coords)
    tm_cy = sum(c[1] for c in tm_coords) / len(tm_coords)
    tm_cz = sum(c[2] for c in tm_coords) / len(tm_coords)
    tm_centroid = [round(tm_cx, 4), round(tm_cy, 4), round(tm_cz, 4)]

    classifications: dict[str, Any] = {}

    for resid, coords in lipid_atoms.items():
        lc_x = sum(c[0] for c in coords) / len(coords)
        lc_y = sum(c[1] for c in coords) / len(coords)
        lc_z = sum(c[2] for c in coords) / len(coords)
        lc = [round(lc_x, 4), round(lc_y, 4), round(lc_z, 4)]

        # 3D distance to closest TM atom
        d_min_tm = min(
            math.sqrt((a[0]-pa[0])**2 + (a[1]-pa[1])**2 + (a[2]-pa[2])**2)
            for a in coords for pa in tm_coords
        ) if tm_coords else 999.0

        # 3D distance to closest non-TM protein atom
        d_min_soluble = min(
            math.sqrt((a[0]-pa[0])**2 + (a[1]-pa[1])**2 + (a[2]-pa[2])**2)
            for a in coords for pa in soluble_coords
        ) if soluble_coords else 999.0

        max_gap = get_angular_gap(lc_x, lc_y, lc_z, tm_ca_coords)
        is_surrounded = max_gap < 185.0

        # Classification logic
        if tm_z_min <= lc_z <= tm_z_max:
            if is_surrounded:
                cls = "tm_cavity_trapped"
            elif d_min_tm < mask_padding_nm:
                cls = "tm_surface_overlap"
            else:
                cls = "bulk_lipid"
        else:
            if d_min_soluble < mask_padding_nm:
                cls = "soluble_domain_adjacent"
            else:
                cls = "bulk_lipid"

        classifications[str(resid)] = {
            "classification": cls,
            "center": lc,
            "min_dist_tm_nm": round(d_min_tm, 4),
            "min_dist_soluble_nm": round(d_min_soluble, 4),
            "max_gap_deg": round(max_gap, 2),
        }

    counts = {
        "bulk_lipid": sum(1 for c in classifications.values() if c["classification"] == "bulk_lipid"),
        "tm_surface_overlap": sum(1 for c in classifications.values() if c["classification"] == "tm_surface_overlap"),
        "tm_cavity_trapped": sum(1 for c in classifications.values() if c["classification"] == "tm_cavity_trapped"),
        "soluble_domain_adjacent": sum(1 for c in classifications.values() if c["classification"] == "soluble_domain_adjacent"),
    }

    mask_report = {
        "tm_z_min": round(tm_z_min, 4),
        "tm_z_max": round(tm_z_max, 4),
        "tm_centroid": tm_centroid,
        "n_tm_atoms": len(tm_coords),
        "grid_spacing_nm": grid_spacing_nm,
        "mask_padding_nm": mask_padding_nm,
    }

    classification_report = {
        "lipid_resname": lipid_resname,
        "total_lipids_checked": len(lipid_atoms),
        "counts": counts,
        "classifications": classifications,
    }

    # Write reports
    out_dir = path.parent
    (out_dir / "tm_mask_report.json").write_text(json.dumps(mask_report, indent=2))
    (out_dir / "tm_lipid_classification_report.json").write_text(json.dumps(classification_report, indent=2))

    return {
        "mask_report": mask_report,
        "classification_report": classification_report,
    }


def _remove_residues_from_gro(gro_path: Path, resids: set[int]) -> int:
    """
    Remove all atoms belonging to the given GRO residue numbers from gro_path.
    Rewrites the file in-place with the corrected atom-count header line.
    Returns the number of atoms removed.
    """
    lines = gro_path.read_text().splitlines()
    title = lines[0]
    atom_lines = lines[2:-1]   # skip title, atom-count, box-vector
    box = lines[-1]

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
        if resnum in resids:
            n_removed += 1
        else:
            kept.append(line)

    new_text = f"{title}\n{len(kept)}\n" + "\n".join(kept) + f"\n{box}\n"
    gro_path.write_text(new_text)
    return n_removed


def renumber_gro_residues_python(gro_path: Path) -> None:
    """Renumber residues sequentially starting from 1 in-place using Python."""
    lines = gro_path.read_text().splitlines()
    if len(lines) < 3:
        return
    title = lines[0]
    atom_count = lines[1]
    atom_lines = lines[2:-1]
    box = lines[-1]

    new_atom_lines = []
    last_resnum = -1
    current_resnum = 0

    for line in atom_lines:
        if len(line) < 5:
            new_atom_lines.append(line)
            continue
        try:
            resnum = int(line[0:5])
        except ValueError:
            new_atom_lines.append(line)
            continue

        if resnum != last_resnum:
            current_resnum += 1
            last_resnum = resnum

        new_resnum_str = f"{current_resnum % 100000:>5}"
        new_line = new_resnum_str + line[5:]
        new_atom_lines.append(new_line)

    new_text = f"{title}\n{atom_count}\n" + "\n".join(new_atom_lines) + f"\n{box}\n"
    gro_path.write_text(new_text)


def exclude_lipids_tm_aware(
    gro_path: Path | str,
    tm_residues: Optional[Set[int] | List[int] | set[int] | list[int]],
    lipid_resname: str,
    policy: str,
    safety_limit: int = 50,
    dry_run: bool = False,
) -> dict[str, Any]:
    """
    Controlled TM-aware lipid exclusion. Removes targeted TM lipids based on classification.
    """
    gro_path = Path(gro_path)
    if not gro_path.exists():
        raise FileNotFoundError(f"GRO file not found: {gro_path}")

    # Initialize default report
    report = {
        "enabled": policy != "off",
        "policy": policy,
        "input_gro": str(gro_path.name),
        "output_gro": str(gro_path.name),
        "n_total_lipids": 0,
        "n_bulk_lipids": 0,
        "n_soluble_domain_adjacent": 0,
        "n_tm_surface_overlap": 0,
        "n_tm_cavity_trapped": 0,
        # tm_cavity_trapped is a hard structural invalidity (lipid geometrically
        # surrounded by the protein) and is always removed, regardless of
        # safety_limit. tm_surface_overlap is a soft annular-density trimming
        # decision and remains capped by safety_limit. See
        # docs/audits/simforge_vs_protmemfiles_audit.md (Fix 2) — the original
        # protmemfiles inflategro-Jorge.pl overlap removal never capped removal
        # counts at all; safety_limit is a SimForge addition and must not block
        # the hard-invalidity class.
        "n_tm_cavity_trapped_removed": 0,
        "n_tm_surface_overlap_removed": 0,
        "n_lipids_removed": 0,
        "removed_resids": [],
        "removed_resnames": [],
        "removed_classes": [],
        "safety_limit": safety_limit,
        # surface_overlap_safety_status: whether safety_limit capped the
        # tm_surface_overlap subset. safety_status is a back-compat alias of
        # this field — cavity-trapped removal can no longer be blocked, so the
        # only thing safety_status can still report is the surface-overlap cap.
        "surface_overlap_safety_status": "ok",
        "safety_status": "ok",
        "post_exclusion_n_tm_surface_overlap": 0,
        "post_exclusion_n_tm_cavity_trapped": 0,
        "coordinate_file_modified": False,
        "topology_modified": False,
    }

    if policy == "off":
        report_file = gro_path.parent / "tm_aware_lipid_exclusion_report.json"
        report_file.write_text(json.dumps(report, indent=2))
        return report

    # 1. Run classifier on input gro
    res = classify_lipids_tm_aware(
        gro_path=gro_path,
        tm_residues=tm_residues,
        lipid_resname=lipid_resname,
    )
    class_report = res["classification_report"]
    counts = class_report["counts"]

    n_bulk = counts.get("bulk_lipid", 0)
    n_surface = counts.get("tm_surface_overlap", 0)
    n_cavity = counts.get("tm_cavity_trapped", 0)
    n_soluble = counts.get("soluble_domain_adjacent", 0)
    n_total = class_report["total_lipids_checked"]

    report.update({
        "n_total_lipids": n_total,
        "n_bulk_lipids": n_bulk,
        "n_soluble_domain_adjacent": n_soluble,
        "n_tm_surface_overlap": n_surface,
        "n_tm_cavity_trapped": n_cavity,
    })

    # Identify removal candidates, split by class.
    cavity_resids  = []
    surface_resids = []

    for resid_str, item in class_report["classifications"].items():
        cls = item["classification"]
        if cls == "tm_cavity_trapped":
            cavity_resids.append(int(resid_str))
        elif cls == "tm_surface_overlap":
            surface_resids.append(int(resid_str))

    # 2. Safety limit applies only to tm_surface_overlap. tm_cavity_trapped is
    # never capped — it is a hard structural invalidity (the lipid is
    # geometrically enclosed by the protein), and the original protmemfiles
    # methodology removed every such clash unconditionally.
    surface_exceeded = len(surface_resids) > safety_limit
    surface_status = "exceeded" if surface_exceeded else "ok"
    report.update({
        "surface_overlap_safety_status": surface_status,
        "safety_status": surface_status,   # back-compat alias
    })

    if dry_run:
        # Dry run never mutates the file and never enforces safety_limit —
        # it reports what a full (uncapped) apply would remove.
        apply_cavity_resids  = list(cavity_resids)
        apply_surface_resids = list(surface_resids)
    else:
        apply_cavity_resids  = list(cavity_resids)
        apply_surface_resids = [] if surface_exceeded else list(surface_resids)

    to_remove_resids   = apply_cavity_resids + apply_surface_resids
    to_remove_resnames = [lipid_resname] * len(to_remove_resids)
    to_remove_classes  = (
        ["tm_cavity_trapped"]  * len(apply_cavity_resids)
        + ["tm_surface_overlap"] * len(apply_surface_resids)
    )
    n_to_remove = len(to_remove_resids)

    # 3. Perform removal
    if policy in ("apply", "strict_apply") and n_to_remove > 0 and not dry_run:
        # Delete targeted residue IDs
        _remove_residues_from_gro(gro_path, set(to_remove_resids))

        # Renumber residues sequentially
        try:
            import subprocess
            ret = subprocess.run(
                ["gmx", "editconf", "-f", str(gro_path), "-o", str(gro_path), "-resnr", "1"],
                capture_output=True, text=True,
            )
            if ret.returncode != 0:
                renumber_gro_residues_python(gro_path)
        except Exception:
            renumber_gro_residues_python(gro_path)

        report.update({
            "n_lipids_removed": n_to_remove,
            "n_tm_cavity_trapped_removed": len(apply_cavity_resids),
            "n_tm_surface_overlap_removed": len(apply_surface_resids),
            "removed_resids": sorted(to_remove_resids),
            "removed_resnames": to_remove_resnames,
            "removed_classes": to_remove_classes,
            "coordinate_file_modified": True,
        })

        # Re-run classifier
        res_post = classify_lipids_tm_aware(
            gro_path=gro_path,
            tm_residues=tm_residues,
            lipid_resname=lipid_resname,
        )
        post_counts = res_post["classification_report"]["counts"]
        post_surface = post_counts.get("tm_surface_overlap", 0)
        post_cavity = post_counts.get("tm_cavity_trapped", 0)

        report.update({
            "post_exclusion_n_tm_surface_overlap": post_surface,
            "post_exclusion_n_tm_cavity_trapped": post_cavity,
        })

        # 4. Strict check — this is the only remaining path that can raise for
        # a surface_limit cap: if capped surface-overlap lipids (or any
        # cavity-trapped residual the classifier still finds) remain after
        # removal, strict_apply blocks. A plain "apply" policy does not raise;
        # it proceeds with whatever was safely removed, matching the original
        # methodology's unconditional-removal behavior for the cavity class.
        if policy == "strict_apply" and (post_surface > 0 or post_cavity > 0):
            report_file = gro_path.parent / "tm_aware_lipid_exclusion_report.json"
            report_file.write_text(json.dumps(report, indent=2))
            raise ValueError(
                f"TM-aware strict_apply validation failed: {post_surface} surface overlaps "
                f"and {post_cavity} cavity-trapped lipids remain after exclusion."
            )
    elif policy in ("apply", "strict_apply") and n_to_remove == 0 and not dry_run:
        # Nothing was safe to remove this pass (e.g. surface capped and no
        # cavity-trapped lipids found). Report the unchanged state.
        report.update({
            "post_exclusion_n_tm_surface_overlap": n_surface,
            "post_exclusion_n_tm_cavity_trapped": n_cavity,
        })
        if policy == "strict_apply" and (n_surface > 0 or n_cavity > 0):
            report_file = gro_path.parent / "tm_aware_lipid_exclusion_report.json"
            report_file.write_text(json.dumps(report, indent=2))
            raise ValueError(
                f"TM-aware strict_apply validation failed: {n_surface} surface overlaps "
                f"and {n_cavity} cavity-trapped lipids remain after exclusion."
            )
    else:
        # warn mode, off, or dry_run
        report.update({
            "n_lipids_removed": n_to_remove if dry_run else 0,
            "n_tm_cavity_trapped_removed": len(apply_cavity_resids) if dry_run else 0,
            "n_tm_surface_overlap_removed": len(apply_surface_resids) if dry_run else 0,
            "removed_resids": sorted(to_remove_resids) if dry_run else [],
            "removed_resnames": to_remove_resnames if dry_run else [],
            "removed_classes": to_remove_classes if dry_run else [],
            "post_exclusion_n_tm_surface_overlap": 0 if dry_run else n_surface,
            "post_exclusion_n_tm_cavity_trapped": 0 if dry_run else n_cavity,
            "coordinate_file_modified": False,
        })

    # Write final report
    report_file = gro_path.parent / "tm_aware_lipid_exclusion_report.json"
    report_file.write_text(json.dumps(report, indent=2))

    return report
