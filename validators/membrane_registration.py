# validators/membrane_registration.py
"""
Adaptive Membrane Registration — Phase 8A.

Searches Z-shifts of the bilayer relative to the protein TM bundle to find
the placement that maximises TM burial while minimising soluble-domain
interference with the hydrophobic core.

Public API
----------
score_registration_candidate()    — score one shift value
search_best_registration()        — exhaustive search → best candidate
adaptive_membrane_registration()  — full pipeline: read GROs, search, write report
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


# ── GRO parsing ───────────────────────────────────────────────────────────────

def _parse_gro_atoms(gro_path: Path) -> tuple[list[tuple[int, str, str, float, float, float]], list[float]]:
    """
    Return (atom_tuples, box_xyz) from a GRO file.
    atom_tuples: list of (resnum, resname, atomname, x, y, z)
    """
    lines = gro_path.read_text().splitlines()
    n = int(lines[1].strip())
    atoms: list[tuple[int, str, str, float, float, float]] = []
    for ln in lines[2: 2 + n]:
        if len(ln) < 44:
            continue
        try:
            resnum   = int(ln[0:5].strip())
            resname  = ln[5:10].strip()
            atomname = ln[10:15].strip()
            x = float(ln[20:28])
            y = float(ln[28:36])
            z = float(ln[36:44])
            atoms.append((resnum, resname, atomname, x, y, z))
        except (ValueError, IndexError):
            continue
    box_parts = lines[2 + n].split()
    box = [float(v) for v in box_parts[:3]] if len(box_parts) >= 3 else [0.0, 0.0, 0.0]
    return atoms, box


def _bilayer_midplane_from_atoms(atoms: list[tuple]) -> float:
    """Estimate bilayer midplane Z as mean of all Z values (works for pure bilayer GRO)."""
    zs = [a[5] for a in atoms]
    return sum(zs) / len(zs) if zs else 0.0


# ── Scoring ───────────────────────────────────────────────────────────────────

@dataclass
class RegistrationCandidate:
    """Score record for one candidate Z shift."""
    shift_z_nm:                          float       # delta relative to TM-center alignment
    actual_z_shift_nm:                   float       # absolute shift to apply to bilayer atoms
    new_bilayer_midplane_z:              float       # bilayer midplane Z after shift
    core_z_bot:                          float       # hydrophobic core bottom boundary
    core_z_top:                          float       # hydrophobic core top boundary
    tm_burial_fraction:                  float       # fraction of all TM atoms in core
    tm_ca_burial_fraction:               float       # fraction of TM Cα atoms in core
    non_tm_core_atom_count:              int
    non_tm_core_residue_count:           int
    ec_core_atom_count:                  int
    ic_core_atom_count:                  int
    soluble_domain_core_penetration_fraction: float  # non_tm_core / total_non_tm
    estimated_lipid_interference_score:  float       # EC+IC atoms in core / total protein
    expected_void_risk_score:            float       # 1 - tm_ca_burial_fraction (proxy)
    registration_score:                  float
    warnings:                            list[str] = field(default_factory=list)


def score_registration_candidate(
    protein_atoms:              list[tuple[int, str, str, float, float, float]],
    shift_z_nm:                 float,
    initial_bilayer_midplane_z: float,
    initial_tm_center_z:        float,
    tm_residues:                set[int],
    ec_residues:                set[int],
    ic_residues:                set[int],
    hydrophobic_half_thickness: float,
    max_allowed_shift_nm:       float,
    w_burial:     float = 2.0,
    w_soluble:    float = 0.5,
    w_ec:         float = 0.3,
    w_ic:         float = 0.3,
    w_extreme:    float = 0.3,
) -> RegistrationCandidate:
    """
    Compute all metrics and registration score for a single candidate Z shift.

    shift_z_nm is the adjustment relative to the default TM-center alignment.
    actual_z_shift_nm = (initial_tm_center_z + shift_z_nm) - initial_bilayer_midplane_z.
    """
    actual_shift = (initial_tm_center_z + shift_z_nm) - initial_bilayer_midplane_z
    new_midplane = initial_bilayer_midplane_z + actual_shift
    core_bot     = new_midplane - hydrophobic_half_thickness
    core_top     = new_midplane + hydrophobic_half_thickness

    warnings: list[str] = []

    # — Classify protein atoms into TM / EC / IC / other ——————————————————————
    tm_atoms:      list[tuple] = []
    tm_ca_atoms:   list[tuple] = []
    ec_atoms:      list[tuple] = []
    ic_atoms:      list[tuple] = []
    soluble_atoms: list[tuple] = []     # non-TM protein atoms

    seen_tm_residues:  set[int] = set()
    seen_non_tm_res:   set[int] = set()

    for atom in protein_atoms:
        resnum, resname, atomname, x, y, z = atom
        in_tm = resnum in tm_residues
        in_ec = resnum in ec_residues
        in_ic = resnum in ic_residues
        if in_tm:
            tm_atoms.append(atom)
            seen_tm_residues.add(resnum)
            if atomname == "CA":
                tm_ca_atoms.append(atom)
        else:
            soluble_atoms.append(atom)
            seen_non_tm_res.add(resnum)
            if in_ec:
                ec_atoms.append(atom)
            elif in_ic:
                ic_atoms.append(atom)

    # — TM burial metrics —————————————————————————————————————————————————————
    n_tm_total      = len(tm_atoms)
    n_tm_ca_total   = len(tm_ca_atoms)
    n_tm_buried     = sum(1 for a in tm_atoms    if core_bot <= a[5] <= core_top)
    n_tm_ca_buried  = sum(1 for a in tm_ca_atoms if core_bot <= a[5] <= core_top)

    tm_burial_fraction    = n_tm_buried    / max(1, n_tm_total)
    tm_ca_burial_fraction = n_tm_ca_buried / max(1, n_tm_ca_total)

    if n_tm_total == 0:
        warnings.append("No TM atoms found for burial scoring — shift scoring unreliable")

    # — Soluble domain penetration ——————————————————————————————————————————————
    n_soluble_total    = len(soluble_atoms)
    n_non_tm_in_core   = sum(1 for a in soluble_atoms if core_bot <= a[5] <= core_top)
    non_tm_res_in_core = {a[0] for a in soluble_atoms if core_bot <= a[5] <= core_top}
    soluble_penetration_fraction = n_non_tm_in_core / max(1, n_soluble_total)

    # — EC/IC in core ——————————————————————————————————————————————————————————
    ec_in_core = sum(1 for a in ec_atoms if core_bot <= a[5] <= core_top)
    ic_in_core = sum(1 for a in ic_atoms if core_bot <= a[5] <= core_top)

    # — Lipid interference score ———————————————————————————————————————————————
    total_protein = max(1, len(protein_atoms))
    lipid_interference = (ec_in_core + ic_in_core) / total_protein

    # — Void risk proxy ————————————————————————————————————————————————————————
    void_risk = max(0.0, 1.0 - tm_ca_burial_fraction) ** 2

    # — Extreme shift penalty ———————————————————————————————————————————————————
    extreme_penalty = (shift_z_nm / max(max_allowed_shift_nm, 0.01)) ** 2

    # — EC/IC fractions for weighting ——————————————————————————————————————————
    ec_fraction = ec_in_core / max(1, len(ec_atoms))
    ic_fraction = ic_in_core / max(1, len(ic_atoms))

    # — Final score ————————————————————————————————————————————————————————————
    score = (
        w_burial  * tm_ca_burial_fraction
        - w_soluble * soluble_penetration_fraction
        - w_ec      * ec_fraction
        - w_ic      * ic_fraction
        - w_extreme * extreme_penalty
        - void_risk
    )

    return RegistrationCandidate(
        shift_z_nm                          = round(shift_z_nm, 4),
        actual_z_shift_nm                   = round(actual_shift, 4),
        new_bilayer_midplane_z              = round(new_midplane, 4),
        core_z_bot                          = round(core_bot, 4),
        core_z_top                          = round(core_top, 4),
        tm_burial_fraction                  = round(tm_burial_fraction, 4),
        tm_ca_burial_fraction               = round(tm_ca_burial_fraction, 4),
        non_tm_core_atom_count              = n_non_tm_in_core,
        non_tm_core_residue_count           = len(non_tm_res_in_core),
        ec_core_atom_count                  = ec_in_core,
        ic_core_atom_count                  = ic_in_core,
        soluble_domain_core_penetration_fraction = round(soluble_penetration_fraction, 4),
        estimated_lipid_interference_score  = round(lipid_interference, 4),
        expected_void_risk_score            = round(void_risk, 4),
        registration_score                  = round(score, 6),
        warnings                            = warnings,
    )


# ── Search ────────────────────────────────────────────────────────────────────

def search_best_registration(
    protein_atoms:              list[tuple[int, str, str, float, float, float]],
    initial_bilayer_midplane_z: float,
    initial_tm_center_z:        float,
    tm_residues:                set[int],
    ec_residues:                set[int],
    ic_residues:                set[int],
    search_min_nm:              float = -1.5,
    search_max_nm:              float = 1.5,
    search_step_nm:             float = 0.1,
    hydrophobic_half_thickness: float = 1.25,
    max_allowed_shift_nm:       float = 1.5,
) -> tuple[RegistrationCandidate, list[RegistrationCandidate]]:
    """
    Evaluate all candidate shifts and return (best_candidate, all_candidates).

    The search is over δ ∈ [search_min, search_max] where δ = 0 is the
    default TM-center alignment. The best candidate maximises registration_score.
    If no candidate improves over the default (δ=0), returns δ=0.
    """
    # Generate candidate shifts (rounded to avoid float noise)
    steps = round((search_max_nm - search_min_nm) / search_step_nm) + 1
    deltas = [
        round(search_min_nm + i * search_step_nm, 6)
        for i in range(steps)
    ]
    # Always include the default (δ=0) even if not in grid
    if 0.0 not in deltas:
        deltas.append(0.0)
    deltas.sort()

    candidates = [
        score_registration_candidate(
            protein_atoms=protein_atoms,
            shift_z_nm=delta,
            initial_bilayer_midplane_z=initial_bilayer_midplane_z,
            initial_tm_center_z=initial_tm_center_z,
            tm_residues=tm_residues,
            ec_residues=ec_residues,
            ic_residues=ic_residues,
            hydrophobic_half_thickness=hydrophobic_half_thickness,
            max_allowed_shift_nm=max_allowed_shift_nm,
        )
        for delta in deltas
    ]

    best = max(candidates, key=lambda c: c.registration_score)

    # If the best score is not better than δ=0, fall back to default alignment
    default_cands = [c for c in candidates if c.shift_z_nm == 0.0]
    if default_cands:
        default = default_cands[0]
        if best.registration_score <= default.registration_score:
            best = default  # no improvement — keep default

    return best, candidates


# ── Top-level entry point ─────────────────────────────────────────────────────

def adaptive_membrane_registration(
    protein_gro:             str | Path,
    bilayer_gro:             str | Path,
    tm_residues:             Optional[set[int]],
    ec_residues:             Optional[set[int]] = None,
    ic_residues:             Optional[set[int]] = None,
    search_min_nm:           float = -1.5,
    search_max_nm:           float = 1.5,
    search_step_nm:          float = 0.1,
    hydrophobic_half_thickness_nm: float = 1.25,
    policy:                  str = "warn",
    max_allowed_shift_nm:    float = 1.5,
    output_dir:              Optional[str | Path] = None,
) -> dict:
    """
    Full adaptive membrane registration pipeline.

    Reads protein_gro and bilayer_gro, searches Z shifts, selects the best,
    and writes membrane_registration_report.json to output_dir (or cwd).

    Returns the report dict.
    """
    protein_gro = Path(protein_gro)
    bilayer_gro = Path(bilayer_gro)
    output_dir  = Path(output_dir) if output_dir else protein_gro.parent

    warnings_out: list[str] = []

    # ── No TM annotation ──────────────────────────────────────────────────────
    if not tm_residues:
        report = {
            "enabled":                   True,
            "protein_frame_used":        protein_gro.name,
            "selected_shift_z_nm":       0.0,
            "actual_z_shift_nm":         None,
            "initial_tm_center_z":       None,
            "initial_bilayer_midplane_z": None,
            "final_bilayer_midplane_z":  None,
            "search_range_nm":           [search_min_nm, search_max_nm],
            "search_step_nm":            search_step_nm,
            "candidates":                [],
            "best_candidate":            None,
            "score_components":          {},
            "tm_burial_fraction":        None,
            "tm_ca_burial_fraction":     None,
            "soluble_domain_core_penetration_fraction": None,
            "non_tm_core_atom_count":    None,
            "ec_core_atom_count":        None,
            "ic_core_atom_count":        None,
            "expected_void_risk_score":  None,
            "decision_reason":           "No TM annotation — adaptive registration skipped, keeping default TM-center alignment",
            "warnings":                  ["No TM residues provided — adaptive registration requires transmembrane_segments annotation"],
        }
        _write_report(report, output_dir)
        return report

    # ── Parse GROs ────────────────────────────────────────────────────────────
    if not protein_gro.exists():
        raise FileNotFoundError(f"Protein GRO not found: {protein_gro}")
    if not bilayer_gro.exists():
        raise FileNotFoundError(f"Bilayer GRO not found: {bilayer_gro}")

    prot_atoms, _    = _parse_gro_atoms(protein_gro)
    bil_atoms,  _    = _parse_gro_atoms(bilayer_gro)

    initial_bil_midplane = _bilayer_midplane_from_atoms(bil_atoms)

    # ── TM CA center ──────────────────────────────────────────────────────────
    tm_ca_zs = [
        a[5] for a in prot_atoms
        if a[0] in tm_residues and a[2] == "CA"
    ]
    if not tm_ca_zs:
        warnings_out.append("No TM Cα atoms found in protein GRO — using full protein Z centre")
        all_zs = [a[5] for a in prot_atoms]
        initial_tm_center = (min(all_zs) + max(all_zs)) / 2.0
    else:
        initial_tm_center = sum(tm_ca_zs) / len(tm_ca_zs)

    ec_set = ec_residues or set()
    ic_set = ic_residues or set()

    # ── Search ────────────────────────────────────────────────────────────────
    best, candidates = search_best_registration(
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
        max_allowed_shift_nm       = max_allowed_shift_nm,
    )

    # ── Policy check ──────────────────────────────────────────────────────────
    blocked = False
    block_reason = None
    if abs(best.shift_z_nm) > max_allowed_shift_nm:
        msg = (
            f"Best registration shift {best.shift_z_nm:+.3f} nm exceeds "
            f"max_allowed_shift_nm ({max_allowed_shift_nm} nm)"
        )
        warnings_out.append(msg)
        if policy == "strict":
            blocked = True
            block_reason = msg

    # ── Decision reason ───────────────────────────────────────────────────────
    if best.shift_z_nm == 0.0:
        decision_reason = "Default TM-center alignment is optimal — no adjustment needed"
    else:
        decision_reason = (
            f"Shift {best.shift_z_nm:+.3f} nm from TM-center alignment maximises "
            f"tm_ca_burial_fraction={best.tm_ca_burial_fraction:.3f} "
            f"(score={best.registration_score:.4f})"
        )
    if blocked:
        decision_reason = f"BLOCKED: {block_reason}"

    # ── Build report ──────────────────────────────────────────────────────────
    best_dict = _candidate_to_dict(best)
    cand_list = [_candidate_to_dict(c) for c in candidates]

    report = {
        "enabled":                   True,
        "protein_frame_used":        protein_gro.name,
        "blocked":                   blocked,
        "selected_shift_z_nm":       best.shift_z_nm if not blocked else 0.0,
        "actual_z_shift_nm":         best.actual_z_shift_nm if not blocked else (initial_tm_center - initial_bil_midplane),
        "initial_tm_center_z":       round(initial_tm_center, 4),
        "tm_center_z_in_registration_frame": round(initial_tm_center, 4),
        "bilayer_midplane_z_before": round(initial_bil_midplane, 4),
        "initial_bilayer_midplane_z": round(initial_bil_midplane, 4),
        "final_bilayer_midplane_z":  round(best.new_bilayer_midplane_z if not blocked else initial_tm_center, 4),
        "search_range_nm":           [search_min_nm, search_max_nm],
        "search_step_nm":            search_step_nm,
        "candidates":                cand_list,
        "best_candidate":            best_dict,
        "score_components":          {
            "w_burial":   2.0,
            "w_soluble":  0.5,
            "w_ec":       0.3,
            "w_ic":       0.3,
            "w_extreme":  0.3,
        },
        "tm_burial_fraction":         best.tm_burial_fraction,
        "tm_ca_burial_fraction":      best.tm_ca_burial_fraction,
        "soluble_domain_core_penetration_fraction": best.soluble_domain_core_penetration_fraction,
        "non_tm_core_atom_count":     best.non_tm_core_atom_count,
        "ec_core_atom_count":         best.ec_core_atom_count,
        "ic_core_atom_count":         best.ic_core_atom_count,
        "expected_void_risk_score":   best.expected_void_risk_score,
        "decision_reason":            decision_reason,
        "warnings":                   best.warnings + warnings_out,
    }

    _write_report(report, output_dir)
    return report


def _candidate_to_dict(c: RegistrationCandidate) -> dict:
    return {
        "shift_z_nm":                          c.shift_z_nm,
        "actual_z_shift_nm":                   c.actual_z_shift_nm,
        "new_bilayer_midplane_z":              c.new_bilayer_midplane_z,
        "core_z_bot":                          c.core_z_bot,
        "core_z_top":                          c.core_z_top,
        "tm_burial_fraction":                  c.tm_burial_fraction,
        "tm_ca_burial_fraction":               c.tm_ca_burial_fraction,
        "non_tm_core_atom_count":              c.non_tm_core_atom_count,
        "non_tm_core_residue_count":           c.non_tm_core_residue_count,
        "ec_core_atom_count":                  c.ec_core_atom_count,
        "ic_core_atom_count":                  c.ic_core_atom_count,
        "soluble_domain_core_penetration_fraction": c.soluble_domain_core_penetration_fraction,
        "estimated_lipid_interference_score":  c.estimated_lipid_interference_score,
        "expected_void_risk_score":            c.expected_void_risk_score,
        "registration_score":                  c.registration_score,
        "warnings":                            c.warnings,
    }


def _write_report(report: dict, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "membrane_registration_report.json").write_text(
        json.dumps(report, indent=2)
    )
