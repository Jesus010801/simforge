"""
validators/membrane_surface_registration.py
Phase 8B: Surface-interference adaptive membrane registration.

Extends Phase 8A by explicitly scoring how much EC/IC/loop protein regions
intersect the bilayer headgroup surfaces (the leaflet interface zones), not
only the hydrophobic core.  The objective is:

    maximise TM burial
    while minimising surface/interface interference by non-TM regions

Architecture
------------
- _parse_gro_atoms, _bilayer_midplane_from_atoms: same as membrane_registration.py
- score_surface_candidate():   per-δ scoring with bilayer-surface terms
- search_best_surface_registration():   exhaustive δ search
- surface_interference_adaptive_registration():  full pipeline → writes
  membrane_surface_registration_report.json
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

# ── GRO helpers ───────────────────────────────────────────────────────────────

def _parse_gro_atoms(gro_path: Path) -> tuple[list[tuple], tuple]:
    """Return (atom_list, box_xyz).  atom = (resid, resname, atomname, x, y, z)."""
    lines = gro_path.read_text().splitlines()
    n = int(lines[1].strip())
    atoms: list[tuple] = []
    for ln in lines[2 : 2 + n]:
        if len(ln) < 44:
            continue
        try:
            resid    = int(ln[0:5])
            resname  = ln[5:10].strip()
            atomname = ln[10:15].strip()
            x = float(ln[20:28])
            y = float(ln[28:36])
            z = float(ln[36:44])
            atoms.append((resid, resname, atomname, x, y, z))
        except ValueError:
            pass
    box_parts = lines[2 + n].split()
    box = (float(box_parts[0]), float(box_parts[1]), float(box_parts[2]))
    return atoms, box


def _bilayer_midplane_from_atoms(bil_atoms: list[tuple]) -> float:
    zs = [a[5] for a in bil_atoms]
    return (min(zs) + max(zs)) / 2.0


def _write_report(report: dict, output_dir: Path) -> None:
    (output_dir / "membrane_surface_registration_report.json").write_text(
        json.dumps(report, indent=2)
    )


# ── Candidate dataclass ───────────────────────────────────────────────────────

@dataclass
class SurfaceRegistrationCandidate:
    shift_z_nm:                          float
    actual_z_shift_nm:                   float
    new_bilayer_midplane_z:              float
    core_z_bot:                          float
    core_z_top:                          float
    upper_headgroup_surface_z:           float  # top of upper headgroup zone
    lower_headgroup_surface_z:           float  # bottom of lower headgroup zone
    tm_burial_score:                     float  # = tm_ca_burial_fraction
    tm_ca_burial_fraction:               float
    tm_burial_fraction:                  float  # all TM atoms
    non_tm_core_penetration:             int    # atom count
    non_tm_core_penetration_fraction:    float
    upper_surface_interference_score:    float  # EC atoms in upper HG zone / total EC
    lower_surface_interference_score:    float  # IC atoms in lower HG zone / total IC
    ec_core_atom_count:                  int
    ic_core_atom_count:                  int
    ec_interface_clashes:                int    # EC atoms in upper headgroup zone
    ic_interface_clashes:                int    # IC atoms in lower headgroup zone
    loop_interface_clashes:              int    # loop atoms in any headgroup zone
    registration_score:                  float
    warnings:                            list[str] = field(default_factory=list)


# ── Scoring ────────────────────────────────────────────────────────────────────

_MEMBRANE_RESNAMES = frozenset({
    'DPP', 'DPPC', 'POPC', 'POPE', 'CHOL', 'PALMITOYLOLEOYLPHOSPHATIDYLCHOLINE',
    'SOL', 'HOH', 'WAT', 'TIP', 'TIP3', 'TIP4',
    'NA', 'CL', 'SOD', 'MG', 'K', 'CA',
})


def score_surface_candidate(
    prot_atoms:               list[tuple],
    initial_bilayer_midplane: float,
    initial_tm_center:        float,
    shift_z_nm:               float,
    tm_residues:              set[int],
    ec_residues:              set[int],
    ic_residues:              set[int],
    hydrophobic_half_thickness: float,
    headgroup_thickness:        float,
    max_allowed_shift_nm:       float,
    w_burial:     float = 2.0,
    w_soluble:    float = 0.5,
    w_ec_core:    float = 0.4,
    w_ic_core:    float = 0.4,
    w_ec_surf:    float = 0.6,
    w_ic_surf:    float = 0.4,
    w_loop_surf:  float = 0.2,
    w_extreme:    float = 0.3,
) -> SurfaceRegistrationCandidate:
    """Score one candidate Z shift combining TM burial and surface interference."""
    new_midplane = initial_tm_center + shift_z_nm
    actual_shift = new_midplane - initial_bilayer_midplane

    core_bot   = new_midplane - hydrophobic_half_thickness
    core_top   = new_midplane + hydrophobic_half_thickness
    hg_low_bot = core_bot - headgroup_thickness   # bottom of lower headgroup zone
    hg_up_top  = core_top + headgroup_thickness   # top of upper headgroup zone

    tm_ca_in_core = tm_ca_total = 0
    tm_all_in_core = tm_all_total = 0
    non_tm_in_core = non_tm_total = 0
    ec_in_core = ec_total = ec_in_upper_hg = 0
    ic_in_core = ic_total = ic_in_lower_hg = 0
    loop_total = loop_in_any_hg = 0

    for (resid, resname, atomname, x, y, z) in prot_atoms:
        if resname in _MEMBRANE_RESNAMES:
            continue

        is_tm   = resid in tm_residues
        is_ec   = resid in ec_residues
        is_ic   = resid in ic_residues
        is_loop = not is_tm and not is_ec and not is_ic

        in_core   = core_bot   <= z <= core_top
        in_up_hg  = core_top   <  z <= hg_up_top
        in_low_hg = hg_low_bot <= z <  core_bot

        if is_tm:
            tm_all_total += 1
            if in_core:
                tm_all_in_core += 1
            if atomname == "CA":
                tm_ca_total += 1
                if in_core:
                    tm_ca_in_core += 1
        else:
            non_tm_total += 1
            if in_core:
                non_tm_in_core += 1

            if is_ec:
                ec_total += 1
                if in_core:
                    ec_in_core += 1
                if in_up_hg:
                    ec_in_upper_hg += 1
            elif is_ic:
                ic_total += 1
                if in_core:
                    ic_in_core += 1
                if in_low_hg:
                    ic_in_lower_hg += 1
            elif is_loop:
                loop_total += 1
                if in_up_hg or in_low_hg:
                    loop_in_any_hg += 1

    tm_ca_frac    = tm_ca_in_core  / max(1, tm_ca_total)
    tm_frac       = tm_all_in_core / max(1, tm_all_total)
    non_tm_frac   = non_tm_in_core / max(1, non_tm_total) if non_tm_total else 0.0
    ec_core_frac  = ec_in_core     / max(1, ec_total)     if ec_total  else 0.0
    ic_core_frac  = ic_in_core     / max(1, ic_total)     if ic_total  else 0.0
    ec_surf_frac  = ec_in_upper_hg / max(1, ec_total)     if ec_total  else 0.0
    ic_surf_frac  = ic_in_lower_hg / max(1, ic_total)     if ic_total  else 0.0
    loop_surf_frac = loop_in_any_hg / max(1, loop_total)  if loop_total else 0.0
    extreme_pen   = (shift_z_nm / max_allowed_shift_nm) ** 2 if max_allowed_shift_nm > 0 else 0.0

    score = (
          w_burial    * tm_ca_frac
        - w_soluble   * non_tm_frac
        - w_ec_core   * ec_core_frac
        - w_ic_core   * ic_core_frac
        - w_ec_surf   * ec_surf_frac
        - w_ic_surf   * ic_surf_frac
        - w_loop_surf * loop_surf_frac
        - w_extreme   * extreme_pen
    )

    return SurfaceRegistrationCandidate(
        shift_z_nm                       = round(shift_z_nm, 4),
        actual_z_shift_nm                = round(actual_shift, 4),
        new_bilayer_midplane_z           = round(new_midplane, 4),
        core_z_bot                       = round(core_bot, 4),
        core_z_top                       = round(core_top, 4),
        upper_headgroup_surface_z        = round(hg_up_top, 4),
        lower_headgroup_surface_z        = round(hg_low_bot, 4),
        tm_burial_score                  = round(tm_ca_frac, 4),
        tm_ca_burial_fraction            = round(tm_ca_frac, 4),
        tm_burial_fraction               = round(tm_frac, 4),
        non_tm_core_penetration          = non_tm_in_core,
        non_tm_core_penetration_fraction = round(non_tm_frac, 4),
        upper_surface_interference_score = round(ec_surf_frac, 4),
        lower_surface_interference_score = round(ic_surf_frac, 4),
        ec_core_atom_count               = ec_in_core,
        ic_core_atom_count               = ic_in_core,
        ec_interface_clashes             = ec_in_upper_hg,
        ic_interface_clashes             = ic_in_lower_hg,
        loop_interface_clashes           = loop_in_any_hg,
        registration_score               = round(score, 6),
        warnings                         = [],
    )


# ── Search ────────────────────────────────────────────────────────────────────

def search_best_surface_registration(
    protein_atoms:              list[tuple],
    initial_bilayer_midplane_z: float,
    initial_tm_center_z:        float,
    tm_residues:                set[int],
    ec_residues:                set[int],
    ic_residues:                set[int],
    search_min_nm:              float,
    search_max_nm:              float,
    search_step_nm:             float,
    hydrophobic_half_thickness: float,
    headgroup_thickness:        float,
    max_allowed_shift_nm:       float,
) -> tuple[SurfaceRegistrationCandidate, list[SurfaceRegistrationCandidate]]:
    """Exhaustive search over δ ∈ [search_min, search_max].

    δ=0 is always included even if it falls outside the search range.
    Returns (best_candidate, all_candidates).
    """
    import math

    n_steps  = round((search_max_nm - search_min_nm) / search_step_nm)
    shifts   = [round(search_min_nm + i * search_step_nm, 6) for i in range(n_steps + 1)]
    if not any(abs(s) < 1e-9 for s in shifts):
        shifts.append(0.0)

    kw = dict(
        prot_atoms               = protein_atoms,
        initial_bilayer_midplane = initial_bilayer_midplane_z,
        initial_tm_center        = initial_tm_center_z,
        tm_residues              = tm_residues,
        ec_residues              = ec_residues,
        ic_residues              = ic_residues,
        hydrophobic_half_thickness = hydrophobic_half_thickness,
        headgroup_thickness      = headgroup_thickness,
        max_allowed_shift_nm     = max_allowed_shift_nm,
    )

    candidates = [score_surface_candidate(shift_z_nm=δ, **kw) for δ in shifts]
    best = max(candidates, key=lambda c: c.registration_score)

    delta0 = next(c for c in candidates if abs(c.shift_z_nm) < 1e-9)
    if best.registration_score < delta0.registration_score:
        best = delta0

    return best, candidates


# ── Candidate → dict ──────────────────────────────────────────────────────────

def _candidate_to_dict(c: SurfaceRegistrationCandidate) -> dict:
    return {
        "shift_z_nm":                          c.shift_z_nm,
        "actual_z_shift_nm":                   c.actual_z_shift_nm,
        "new_bilayer_midplane_z":              c.new_bilayer_midplane_z,
        "core_z_bot":                          c.core_z_bot,
        "core_z_top":                          c.core_z_top,
        "upper_headgroup_surface_z":           c.upper_headgroup_surface_z,
        "lower_headgroup_surface_z":           c.lower_headgroup_surface_z,
        "tm_burial_score":                     c.tm_burial_score,
        "tm_ca_burial_fraction":               c.tm_ca_burial_fraction,
        "tm_burial_fraction":                  c.tm_burial_fraction,
        "non_tm_core_penetration":             c.non_tm_core_penetration,
        "non_tm_core_penetration_fraction":    c.non_tm_core_penetration_fraction,
        "upper_surface_interference_score":    c.upper_surface_interference_score,
        "lower_surface_interference_score":    c.lower_surface_interference_score,
        "ec_core_atom_count":                  c.ec_core_atom_count,
        "ic_core_atom_count":                  c.ic_core_atom_count,
        "ec_interface_clashes":                c.ec_interface_clashes,
        "ic_interface_clashes":                c.ic_interface_clashes,
        "loop_interface_clashes":              c.loop_interface_clashes,
        "registration_score":                  c.registration_score,
        "warnings":                            c.warnings,
    }


# ── Main pipeline ─────────────────────────────────────────────────────────────

def surface_interference_adaptive_registration(
    protein_gro:                  str | Path,
    bilayer_gro:                  str | Path,
    tm_residues:                  Optional[set[int]],
    ec_residues:                  Optional[set[int]] = None,
    ic_residues:                  Optional[set[int]] = None,
    search_min_nm:                float = -1.5,
    search_max_nm:                float = 1.5,
    search_step_nm:               float = 0.1,
    hydrophobic_half_thickness_nm: float = 1.25,
    headgroup_thickness_nm:        float = 0.5,
    policy:                        str   = "warn",
    max_allowed_shift_nm:          float = 1.5,
    min_tm_burial_fraction:        float = 0.30,
    output_dir:                    Optional[str | Path] = None,
) -> dict:
    """
    Full Phase 8B surface-interference registration pipeline.

    Reads protein_gro and bilayer_gro, searches Z shifts, selects the one
    that maximises TM burial while minimising EC/IC/loop surface interference,
    and writes membrane_surface_registration_report.json.

    Returns the report dict.
    """
    protein_gro = Path(protein_gro)
    bilayer_gro = Path(bilayer_gro)
    output_dir  = Path(output_dir) if output_dir else protein_gro.parent

    ec_set = ec_residues or set()
    ic_set = ic_residues or set()
    warnings_out: list[str] = []

    # ── No TM annotation → disabled fallback ──────────────────────────────────
    if not tm_residues:
        report = {
            "enabled":                         True,
            "protein_frame_used":              protein_gro.name,
            "phase":                           "8B",
            "selected_shift_z_nm":             0.0,
            "actual_z_shift_nm":               None,
            "initial_tm_center_z":             None,
            "tm_center_z_in_registration_frame": None,
            "bilayer_midplane_z_before":       None,
            "initial_bilayer_midplane_z":      None,
            "final_bilayer_midplane_z":        None,
            "search_range_nm":                 [search_min_nm, search_max_nm],
            "search_step_nm":                  search_step_nm,
            "candidates":                      [],
            "best_candidate":                  None,
            "candidate_scores":                [],
            "tm_burial_score":                 None,
            "non_tm_core_penetration":         None,
            "upper_surface_interference_score": None,
            "lower_surface_interference_score": None,
            "ec_interface_clashes":            None,
            "ic_interface_clashes":            None,
            "loop_interface_clashes":          None,
            "decision_reason":                 "No TM annotation — surface registration skipped",
            "warnings": ["No TM residues provided — Phase 8B requires transmembrane_segments annotation"],
        }
        _write_report(report, output_dir)
        return report

    # ── Parse GROs ────────────────────────────────────────────────────────────
    if not protein_gro.exists():
        raise FileNotFoundError(f"Protein GRO not found: {protein_gro}")
    if not bilayer_gro.exists():
        raise FileNotFoundError(f"Bilayer GRO not found: {bilayer_gro}")

    prot_atoms, _ = _parse_gro_atoms(protein_gro)
    bil_atoms,  _ = _parse_gro_atoms(bilayer_gro)
    initial_bil_midplane = _bilayer_midplane_from_atoms(bil_atoms)

    # ── TM CA centre ──────────────────────────────────────────────────────────
    tm_ca_zs = [a[5] for a in prot_atoms if a[0] in tm_residues and a[2] == "CA"]
    if not tm_ca_zs:
        warnings_out.append("No TM Cα atoms found — using full protein Z centre")
        all_zs = [a[5] for a in prot_atoms if a[1] not in _MEMBRANE_RESNAMES]
        initial_tm_center = (min(all_zs) + max(all_zs)) / 2.0
    else:
        initial_tm_center = sum(tm_ca_zs) / len(tm_ca_zs)

    # ── Search ────────────────────────────────────────────────────────────────
    best, candidates = search_best_surface_registration(
        protein_atoms              = prot_atoms,
        initial_bilayer_midplane_z = initial_bil_midplane,
        initial_tm_center_z        = initial_tm_center,
        tm_residues                = tm_residues,
        ec_residues                = ec_set,
        ic_residues                = ic_set,
        search_min_nm              = search_min_nm,
        search_max_nm              = search_max_nm,
        search_step_nm             = search_step_nm,
        hydrophobic_half_thickness = hydrophobic_half_thickness_nm,
        headgroup_thickness        = headgroup_thickness_nm,
        max_allowed_shift_nm       = max_allowed_shift_nm,
    )

    # ── TM burial safety check ─────────────────────────────────────────────────
    blocked = False
    block_reason: str | None = None

    if best.tm_ca_burial_fraction < min_tm_burial_fraction:
        msg = (
            f"Best candidate (δ={best.shift_z_nm:+.3f} nm) has "
            f"tm_ca_burial_fraction={best.tm_ca_burial_fraction:.3f} "
            f"below min_tm_burial_fraction={min_tm_burial_fraction}"
        )
        warnings_out.append(msg)
        if policy == "strict":
            blocked = True
            block_reason = msg

    if abs(best.shift_z_nm) > max_allowed_shift_nm:
        msg = (
            f"Best shift {best.shift_z_nm:+.3f} nm exceeds "
            f"max_allowed_shift_nm={max_allowed_shift_nm} nm"
        )
        warnings_out.append(msg)
        if policy == "strict":
            blocked = True
            block_reason = msg

    # ── Decision reason ────────────────────────────────────────────────────────
    if blocked:
        decision_reason = f"BLOCKED: {block_reason}"
    elif best.shift_z_nm == 0.0:
        decision_reason = (
            "Default TM-center alignment is optimal for surface interference scoring"
        )
    else:
        parts = []
        if best.ec_interface_clashes == 0 and best.ic_interface_clashes == 0:
            parts.append("zero EC/IC surface clashes")
        else:
            if best.ec_interface_clashes > 0:
                parts.append(f"ec_clashes={best.ec_interface_clashes}")
            if best.ic_interface_clashes > 0:
                parts.append(f"ic_clashes={best.ic_interface_clashes}")
        decision_reason = (
            f"Shift {best.shift_z_nm:+.3f} nm selected for surface interference "
            f"(score={best.registration_score:.4f}, "
            f"tm_ca_burial={best.tm_ca_burial_fraction:.3f}, "
            f"{', '.join(parts) if parts else 'minimised surface interference'})"
        )

    # ── Fallback to δ=0 when blocked ──────────────────────────────────────────
    if blocked:
        delta0 = next(c for c in candidates if abs(c.shift_z_nm) < 1e-9)
        effective = delta0
    else:
        effective = best

    # ── Build report ──────────────────────────────────────────────────────────
    best_dict = _candidate_to_dict(best)
    cand_list = [_candidate_to_dict(c) for c in candidates]

    report = {
        "enabled":                           True,
        "phase":                             "8B",
        "protein_frame_used":                protein_gro.name,
        "blocked":                           blocked,
        "selected_shift_z_nm":               effective.shift_z_nm,
        "actual_z_shift_nm":                 effective.actual_z_shift_nm,
        "initial_tm_center_z":               round(initial_tm_center, 4),
        "tm_center_z_in_registration_frame": round(initial_tm_center, 4),
        "bilayer_midplane_z_before":         round(initial_bil_midplane, 4),
        "initial_bilayer_midplane_z":        round(initial_bil_midplane, 4),
        "final_bilayer_midplane_z":          round(effective.new_bilayer_midplane_z, 4),
        "search_range_nm":                   [search_min_nm, search_max_nm],
        "search_step_nm":                    search_step_nm,
        "headgroup_thickness_nm":            headgroup_thickness_nm,
        "hydrophobic_half_thickness_nm":     hydrophobic_half_thickness_nm,
        "candidates":                        cand_list,
        "best_candidate":                    best_dict,
        "candidate_scores":                  [
            {"shift_z_nm": c["shift_z_nm"], "registration_score": c["registration_score"]}
            for c in cand_list
        ],
        "tm_burial_score":                   effective.tm_ca_burial_fraction,
        "tm_ca_burial_fraction":             effective.tm_ca_burial_fraction,
        "tm_burial_fraction":                effective.tm_burial_fraction,
        "non_tm_core_penetration":           effective.non_tm_core_penetration,
        "non_tm_core_penetration_fraction":  effective.non_tm_core_penetration_fraction,
        "upper_surface_interference_score":  effective.upper_surface_interference_score,
        "lower_surface_interference_score":  effective.lower_surface_interference_score,
        "ec_core_atom_count":                effective.ec_core_atom_count,
        "ic_core_atom_count":                effective.ic_core_atom_count,
        "ec_interface_clashes":              effective.ec_interface_clashes,
        "ic_interface_clashes":              effective.ic_interface_clashes,
        "loop_interface_clashes":            effective.loop_interface_clashes,
        "score_components": {
            "w_burial":    2.0,
            "w_soluble":   0.5,
            "w_ec_core":   0.4,
            "w_ic_core":   0.4,
            "w_ec_surf":   0.6,
            "w_ic_surf":   0.4,
            "w_loop_surf": 0.2,
            "w_extreme":   0.3,
        },
        "decision_reason": decision_reason,
        "warnings":        best.warnings + warnings_out,
    }

    _write_report(report, output_dir)
    return report
