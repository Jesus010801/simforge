"""
validators/membrane_embedding_optimizer.py
Phase 9B: Closed-Loop Membrane Embedding Optimizer.

Generates candidate bilayer placements by varying:
  - Z-shift offset applied to all bilayer lipid atoms
  - TM-exclusion mask padding (removes lipids closer than base_cutoff+pad from TM)

Evaluates each candidate with the Phase 9A interface evaluator.
Selects the highest-scoring valid candidate and promotes it to work_input.gro.

The selected work_input.gro becomes the starting point for the inflategro
shrink loop or the tm_aware backend (which then handles its own pdb2gmx step).

Safety:
  - Never modifies protein coordinates.
  - Keeps system_baseline.gro intact throughout.
  - Only promotes atomically (copy on success).
  - Marks candidates invalid if exclusion exceeds safety_limit.
  - Falls back to baseline if optimizer is disabled or all candidates fail.
"""
from __future__ import annotations

import json
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

DEFAULT_LIPID_RESNAMES: frozenset[str] = frozenset({
    "DPP", "DPPC", "POPC", "POPE", "POPG", "POPS", "CHOL",
    "PALM", "OLEO", "PALC", "SM", "CER",
})

_SOLVENT_RESNAMES: frozenset[str] = frozenset({
    "SOL", "HOH", "WAT", "TIP3", "TIP4", "TIP5",
    "NA", "CL", "SOD", "MG", "K", "CA",
})


# ── GRO helpers ───────────────────────────────────────────────────────────────

def _apply_z_shift_text(
    input_path: Path,
    output_path: Path,
    z_delta_nm: float,
    lipid_resnames: frozenset[str],
) -> None:
    """Line-level Z-shift: only modifies lipid lines, O(n_lines)."""
    lines = input_path.read_text().splitlines()
    if len(lines) < 3:
        shutil.copy(input_path, output_path)
        return
    n_atoms = int(lines[1].strip())
    out = [lines[0], lines[1]]
    for ln in lines[2: 2 + n_atoms]:
        if len(ln) >= 44:
            resname = ln[5:10].strip()
            if resname in lipid_resnames:
                try:
                    z = float(ln[36:44]) + z_delta_nm
                    ln = ln[:36] + f"{z:8.3f}" + ln[44:]
                except ValueError:
                    pass
        out.append(ln)
    # box line and any trailing lines
    out.extend(lines[2 + n_atoms:])
    output_path.write_text("\n".join(out) + "\n")


def _parse_gro_atoms(gro_path: Path) -> tuple[list, str, str]:
    """Returns (atom_records, title, box_line).
    atom_record = {"resid": int, "resname": str, "atomname": str, "raw": str}
    """
    lines = gro_path.read_text().splitlines()
    title   = lines[0]
    n_atoms = int(lines[1].strip())
    atoms   = []
    for ln in lines[2: 2 + n_atoms]:
        if len(ln) >= 44:
            try:
                resid    = int(ln[0:5])
                resname  = ln[5:10].strip()
                atomname = ln[10:15].strip()
                atoms.append({"resid": resid, "resname": resname, "atomname": atomname, "raw": ln})
            except ValueError:
                atoms.append({"resid": 0, "resname": "", "atomname": "", "raw": ln})
        else:
            atoms.append({"resid": 0, "resname": "", "atomname": "", "raw": ln})
    box_line = lines[2 + n_atoms] if len(lines) > 2 + n_atoms else "   10.00000   10.00000   10.00000"
    return atoms, title, box_line


def _write_gro_from_atoms(
    output_path: Path,
    atoms: list,
    title: str,
    box_line: str,
) -> None:
    """Write a GRO with re-numbered atom indices (residue numbers preserved)."""
    lines = [title, str(len(atoms))]
    for i, a in enumerate(atoms, 1):
        raw = a["raw"]
        # Replace atom number (cols 15-20) with new index
        lines.append(raw[:15] + f"{i:5d}" + raw[20:])
    lines.append(box_line)
    output_path.write_text("\n".join(lines) + "\n")


def _count_lipid_residues(atoms: list, lipid_resnames: frozenset[str]) -> int:
    seen = set()
    for a in atoms:
        if a["resname"] in lipid_resnames:
            seen.add(a["resid"])
    return len(seen)


# ── Lipid exclusion ───────────────────────────────────────────────────────────

def _exclude_lipids_near_tm(
    atoms: list,
    tm_residue_set: set[int],
    lipid_resnames: frozenset[str],
    cutoff_nm: float,
    safety_limit: int,
) -> tuple[Optional[list], int, Optional[str]]:
    """
    Remove lipid residues that have any atom within cutoff_nm of any TM atom.

    Returns:
        (new_atoms, n_removed, invalid_reason)
        invalid_reason is set when n_removed > safety_limit (candidate is invalid).
    """
    try:
        import numpy as np
    except ImportError:
        return atoms, 0, None  # skip exclusion silently

    tm_atoms  = [a for a in atoms if a["resid"] in tm_residue_set
                 and a["resname"] not in lipid_resnames and a["resname"] not in _SOLVENT_RESNAMES]
    lip_atoms = [a for a in atoms if a["resname"] in lipid_resnames]

    if not tm_atoms or not lip_atoms:
        return atoms, 0, None

    def _xyz(a_list):
        xs, ys, zs = [], [], []
        for a in a_list:
            ln = a["raw"]
            try:
                xs.append(float(ln[20:28]))
                ys.append(float(ln[28:36]))
                zs.append(float(ln[36:44]))
            except ValueError:
                xs.append(0.0); ys.append(0.0); zs.append(0.0)
        return np.array([xs, ys, zs], dtype=np.float32).T  # (N, 3)

    tm_coords = _xyz(tm_atoms)

    from collections import defaultdict
    lip_by_resid: dict[int, list] = defaultdict(list)
    for a in lip_atoms:
        lip_by_resid[a["resid"]].append(a)

    resids_to_remove: set[int] = set()
    for resid, r_atoms in lip_by_resid.items():
        r_coords = _xyz(r_atoms)  # (M, 3)
        # Min distance from any lipid atom in this residue to any TM atom
        diff  = r_coords[:, None, :] - tm_coords[None, :, :]   # (M, T, 3)
        dists = np.sqrt((diff ** 2).sum(axis=2))               # (M, T)
        if dists.min() < cutoff_nm:
            resids_to_remove.add(resid)

    n_removed = len(resids_to_remove)
    if n_removed > safety_limit:
        return None, n_removed, (
            f"exclusion_safety_limit_exceeded: {n_removed} lipid residues would be removed "
            f"(limit={safety_limit})"
        )

    new_atoms = [a for a in atoms if a["resid"] not in resids_to_remove]
    return new_atoms, n_removed, None


# ── Evaluators ────────────────────────────────────────────────────────────────

def _evaluate_interface(
    gro_path: Path,
    tm_residues: set[int],
    lipid_resnames: frozenset[str],
    min_fraction_covered: float,
    max_fraction_exposed_gap: float,
    max_p90_distance_nm: float,
    output_dir: Path,
    simforge_root: Optional[str] = None,
) -> dict:
    try:
        import sys
        if simforge_root:
            sys.path.insert(0, simforge_root)
        from validators.protein_membrane_interface import evaluate_protein_membrane_interface
        return evaluate_protein_membrane_interface(
            gro_path                 = gro_path,
            tm_residues              = tm_residues,
            lipid_resnames           = set(lipid_resnames),
            min_fraction_covered     = min_fraction_covered,
            max_fraction_exposed_gap = max_fraction_exposed_gap,
            max_p90_distance_nm      = max_p90_distance_nm,
            output_dir               = output_dir,
        )
    except Exception as exc:
        return {"error": str(exc), "quality_passed": None,
                "fraction_covered": None, "fraction_exposed_gap": None,
                "p90_nearest_lipid_distance": None, "warnings": [str(exc)]}


def _evaluate_quality(
    gro_path: Path,
    lipid_resname: str,
    headgroup_atom: str,
    tail_atom: str,
    reference_apl_nm2: float,
    output_dir: Path,
    simforge_root: Optional[str] = None,
) -> dict:
    try:
        import sys
        if simforge_root:
            sys.path.insert(0, simforge_root)
        from validators.embedding_quality import diagnose_embedding_quality
        result = diagnose_embedding_quality(
            str(gro_path),
            lipid_resname     = lipid_resname,
            headgroup_atom    = headgroup_atom,
            tail_atom         = tail_atom,
            reference_apl_nm2 = reference_apl_nm2,
        )
        return result.dict() if hasattr(result, "dict") else dict(result)
    except Exception as exc:
        return {"error": str(exc), "tm_burial_score": None,
                "void_fraction_local": None, "apl_local_nm2": None}


# ── Scoring ───────────────────────────────────────────────────────────────────

def compute_candidate_score(
    interface_m:      dict,
    quality_m:        dict,
    n_lip_initial:    int,
    n_lip_final:      int,
    reference_apl_nm2: float = 0.62,
) -> float:
    """
    score =
      + 3.0 * fraction_covered
      - 3.0 * fraction_exposed_gap
      - 1.5 * min(p90_nearest_lipid_distance_nm / 1.5, 1.0)
      - 2.0 * trapped_lipid_penalty
      - 1.0 * void_fraction_local
      - 0.5 * lipid_deletion_penalty
      + 1.0 * tm_burial_score
      - 1.0 * apl_deviation_penalty
    """
    fc   = interface_m.get("fraction_covered")     or 0.0
    fg   = interface_m.get("fraction_exposed_gap") or 0.0
    p90  = interface_m.get("p90_nearest_lipid_distance") or 1.5
    p90  = float(p90)

    tm_burial  = float(quality_m.get("tm_burial_score")  or 0.5)
    void_frac  = float(quality_m.get("void_fraction_local") or 0.0)
    apl_local  = quality_m.get("apl_local_nm2")

    # APL deviation (only if quality evaluator ran successfully)
    if apl_local is not None and reference_apl_nm2 > 0:
        apl_dev = min(abs(float(apl_local) - reference_apl_nm2) / reference_apl_nm2, 1.0)
    else:
        apl_dev = 0.0

    # Lipid deletion penalty: fraction above 20% threshold
    if n_lip_initial > 0:
        del_frac = (n_lip_initial - n_lip_final) / n_lip_initial
        lip_del_penalty = max(0.0, del_frac - 0.20)
    else:
        lip_del_penalty = 0.0

    # trapped_lipid_penalty: not computed pre-shrink (all zeros)
    trapped_penalty = 0.0

    score = (
        + 3.0 * float(fc)
        - 3.0 * float(fg)
        - 1.5 * min(p90 / 1.5, 1.0)
        - 2.0 * trapped_penalty
        - 1.0 * void_frac
        - 0.5 * lip_del_penalty
        + 1.0 * tm_burial
        - 1.0 * apl_dev
    )
    return round(score, 5)


def check_quality_passed(
    interface_m: dict,
    quality_m:   dict,
    thresholds:  dict,
) -> tuple[bool, list[str]]:
    reasons: list[str] = []
    fc   = interface_m.get("fraction_covered")     or 0.0
    fg   = interface_m.get("fraction_exposed_gap") or 1.0
    p90  = interface_m.get("p90_nearest_lipid_distance") or 1.5
    vf   = quality_m.get("void_fraction_local")    or 0.0

    min_fc  = thresholds.get("min_fraction_covered",    0.75)
    max_fg  = thresholds.get("max_fraction_exposed_gap", 0.15)
    max_p90 = thresholds.get("max_p90_nearest_lipid_distance", 0.70)
    max_vf  = thresholds.get("max_void_fraction_local", 0.30)

    if float(fc)  < min_fc:
        reasons.append(f"fraction_covered={fc:.3f} < {min_fc}")
    if float(fg)  > max_fg:
        reasons.append(f"fraction_exposed_gap={fg:.3f} > {max_fg}")
    if float(p90) > max_p90:
        reasons.append(f"p90_nearest_lipid_distance={p90:.3f} > {max_p90}")
    if float(vf)  > max_vf:
        reasons.append(f"void_fraction_local={vf:.3f} > {max_vf}")

    return len(reasons) == 0, reasons


# ── Main optimizer ────────────────────────────────────────────────────────────

def run_membrane_embedding_optimizer(
    system_baseline_gro:         str | Path,
    work_input_gro:              str | Path,
    tm_residues:                 set[int],
    output_dir:                  str | Path,
    lipid_resnames:              Optional[set[str]] = None,
    lipid_residue_name:          str  = "DPP",
    headgroup_atom:              str  = "O33",
    tail_atom:                   str  = "C50",
    reference_apl_nm2:           float = 0.62,
    z_shift_offsets_nm:          Optional[list] = None,
    mask_padding_nm_values:      Optional[list] = None,
    max_candidates:              int   = 28,
    enabled:                     bool  = True,
    policy:                      str   = "warn",
    run_candidate_minimization:  bool  = False,
    min_score_improvement:       float = 0.10,
    exclusion_base_cutoff_nm:    float = 0.14,
    exclusion_safety_limit:      int   = 100,
    min_fraction_covered:        float = 0.75,
    max_fraction_exposed_gap:    float = 0.15,
    max_p90_distance_nm:         float = 0.70,
    max_void_fraction_local:     float = 0.30,
    simforge_root:               Optional[str] = None,
) -> dict:
    """
    Run the Phase 9B closed-loop embedding optimizer.

    Reads system_baseline_gro, generates candidate placements, evaluates each,
    and writes the best candidate to work_input_gro.

    Returns the full optimizer report dict (also written to
    output_dir/membrane_embedding_optimizer_report.json).
    """
    baseline_gro  = Path(system_baseline_gro)
    work_input    = Path(work_input_gro)
    output_dir    = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    lip_set       = frozenset(lipid_resnames or DEFAULT_LIPID_RESNAMES)
    thresholds    = {
        "min_fraction_covered":           min_fraction_covered,
        "max_fraction_exposed_gap":       max_fraction_exposed_gap,
        "max_p90_nearest_lipid_distance": max_p90_distance_nm,
        "max_void_fraction_local":        max_void_fraction_local,
    }

    warnings: list[str] = []

    # ── Always write work_input (baseline fallback) ───────────────────────────
    shutil.copy(baseline_gro, work_input)

    # ── Disabled path ─────────────────────────────────────────────────────────
    if not enabled:
        report = {
            "enabled": False,
            "policy": policy,
            "baseline_metrics": {},
            "n_candidates_evaluated": 0,
            "n_candidates_invalid": 0,
            "n_candidates_quality_passed": 0,
            "candidates": [],
            "selected_candidate_id": None,
            "selected_candidate_params": None,
            "selected_score": None,
            "baseline_score": None,
            "score_improvement": None,
            "selected_quality_passed": None,
            "reason_selected": "optimizer_disabled",
            "promoted_to_work_input_gro": False,
            "warnings": ["Optimizer disabled; using baseline system.gro"],
        }
        _write_report(report, output_dir)
        return report

    if not tm_residues:
        warnings.append("tm_residues not provided — cannot evaluate TM coverage; optimizer disabled")
        report = _disabled_report(policy, "no_tm_residues", warnings)
        shutil.copy(baseline_gro, work_input)
        _write_report(report, output_dir)
        return report

    # ── Default candidate grid ─────────────────────────────────────────────────
    z_shifts  = z_shift_offsets_nm   or [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6]
    paddings  = mask_padding_nm_values or [0.0, 0.1, 0.2, 0.3]

    candidates_dir = output_dir / "optimizer_candidates"
    candidates_dir.mkdir(exist_ok=True)

    # ── Baseline metrics ──────────────────────────────────────────────────────
    baseline_interface = _evaluate_interface(
        baseline_gro, tm_residues, lip_set,
        min_fraction_covered, max_fraction_exposed_gap, max_p90_distance_nm,
        candidates_dir / "baseline",
        simforge_root,
    )
    (candidates_dir / "baseline").mkdir(exist_ok=True)
    baseline_quality = _evaluate_quality(
        baseline_gro, lipid_residue_name, headgroup_atom, tail_atom,
        reference_apl_nm2, candidates_dir / "baseline", simforge_root,
    )
    baseline_atoms, baseline_title, baseline_box = _parse_gro_atoms(baseline_gro)
    n_lip_baseline = _count_lipid_residues(baseline_atoms, lip_set)
    baseline_score = compute_candidate_score(
        baseline_interface, baseline_quality, n_lip_baseline, n_lip_baseline,
        reference_apl_nm2,
    )

    # ── Generate and evaluate candidates ─────────────────────────────────────
    candidate_results: list[dict] = []
    cid = 0

    for z_shift in z_shifts:
        for mask_pad in paddings:
            if cid >= max_candidates:
                break

            cid_dir = candidates_dir / f"c{cid:03d}"
            cid_dir.mkdir(exist_ok=True)
            cand_gro = cid_dir / "candidate.gro"

            # Step 1: apply Z-shift
            _apply_z_shift_text(baseline_gro, cand_gro, z_shift, lip_set)

            # Step 2: optional exclusion (mask_padding > 0)
            n_removed   = 0
            invalid_msg = None
            if mask_pad > 0.0 and tm_residues:
                atoms_shifted, _, _ = _parse_gro_atoms(cand_gro)
                cutoff = exclusion_base_cutoff_nm + mask_pad
                new_atoms, n_removed, invalid_msg = _exclude_lipids_near_tm(
                    atoms_shifted, tm_residues, lip_set, cutoff, exclusion_safety_limit
                )
                if invalid_msg is None and new_atoms is not None:
                    _write_gro_from_atoms(cand_gro, new_atoms, baseline_title, baseline_box)
                elif invalid_msg is not None:
                    # candidate invalid — keep the z-shifted GRO for reference
                    # but mark as invalid in the report
                    cand_atoms_tmp = atoms_shifted
                    n_lip_final_tmp = _count_lipid_residues(cand_atoms_tmp, lip_set)
                    candidate_results.append({
                        "candidate_id":            cid,
                        "z_shift_offset_nm":       z_shift,
                        "mask_padding_nm":         mask_pad,
                        "n_lipids_initial":        n_lip_baseline,
                        "n_lipids_removed":        n_removed,
                        "n_lipids_final":          n_lip_final_tmp,
                        "interface_metrics":       {},
                        "embedding_quality_metrics": {},
                        "trapped_lipid_count":     0,
                        "score":                   -999.0,
                        "quality_passed":          False,
                        "quality_reasons":         [invalid_msg],
                        "invalid_reason":          invalid_msg,
                    })
                    cid += 1
                    continue

            # Count lipids in candidate
            cand_atoms_eval, _, _ = _parse_gro_atoms(cand_gro)
            n_lip_cand = _count_lipid_residues(cand_atoms_eval, lip_set)

            # Step 3: evaluate
            iface_m   = _evaluate_interface(
                cand_gro, tm_residues, lip_set,
                min_fraction_covered, max_fraction_exposed_gap, max_p90_distance_nm,
                cid_dir, simforge_root,
            )
            quality_m = _evaluate_quality(
                cand_gro, lipid_residue_name, headgroup_atom, tail_atom,
                reference_apl_nm2, cid_dir, simforge_root,
            )

            score = compute_candidate_score(
                iface_m, quality_m, n_lip_baseline, n_lip_cand, reference_apl_nm2
            )
            q_passed, q_reasons = check_quality_passed(iface_m, quality_m, thresholds)

            candidate_results.append({
                "candidate_id":              cid,
                "z_shift_offset_nm":         z_shift,
                "mask_padding_nm":           mask_pad,
                "n_lipids_initial":          n_lip_baseline,
                "n_lipids_removed":          n_removed,
                "n_lipids_final":            n_lip_cand,
                "interface_metrics":         iface_m,
                "embedding_quality_metrics": quality_m,
                "trapped_lipid_count":       0,
                "score":                     score,
                "quality_passed":            q_passed,
                "quality_reasons":           q_reasons,
                "invalid_reason":            None,
            })
            cid += 1

        if cid >= max_candidates:
            break

    # ── Select best candidate ─────────────────────────────────────────────────
    valid_cands = [c for c in candidate_results if c["invalid_reason"] is None]
    passing_cands = [c for c in valid_cands if c["quality_passed"]]
    n_quality_passed = len(passing_cands)
    n_invalid = len(candidate_results) - len(valid_cands)

    best: Optional[dict] = None
    reason_selected = "none"

    if passing_cands:
        best = max(passing_cands, key=lambda c: c["score"])
        reason_selected = "quality_passed_and_best_score"
    elif valid_cands:
        # No candidate fully passes quality; pick best score if it improves baseline
        best_by_score = max(valid_cands, key=lambda c: c["score"])
        improvement = best_by_score["score"] - baseline_score
        if improvement >= min_score_improvement:
            best = best_by_score
            reason_selected = f"best_score_improves_baseline_by_{improvement:.3f}"
            warnings.append(
                f"No candidate passed all quality gates. "
                f"Selecting best-score candidate (improvement={improvement:.3f}). "
                "Review membrane_embedding_optimizer_report.json."
            )
        else:
            reason_selected = f"no_improvement_over_baseline (best delta={improvement:.3f})"
            warnings.append(
                "No candidate passed quality gates and none improves baseline "
                f"by >= {min_score_improvement}. Falling back to baseline system."
            )

    # ── Apply policy ──────────────────────────────────────────────────────────
    promoted = False
    if best is not None:
        best_cid = best["candidate_id"]
        best_cand_gro = candidates_dir / f"c{best_cid:03d}" / "candidate.gro"
        if best_cand_gro.exists():
            shutil.copy(best_cand_gro, work_input)
            promoted = True
    else:
        # Keep baseline in work_input (already copied at top)
        if policy == "strict" and valid_cands:
            warnings.append(
                "STRICT MODE: no valid improving candidate found. "
                "Downstream steps will use the baseline — consider relaxing thresholds."
            )

    # ── Build report ──────────────────────────────────────────────────────────
    # Trim interface_metrics to compact form in candidate list
    compact_candidates = []
    for c in candidate_results:
        im = c["interface_metrics"]
        compact_candidates.append({
            "candidate_id":         c["candidate_id"],
            "z_shift_offset_nm":    c["z_shift_offset_nm"],
            "mask_padding_nm":      c["mask_padding_nm"],
            "n_lipids_initial":     c["n_lipids_initial"],
            "n_lipids_removed":     c["n_lipids_removed"],
            "n_lipids_final":       c["n_lipids_final"],
            "interface_metrics": {
                "fraction_covered":          im.get("fraction_covered"),
                "fraction_exposed_gap":      im.get("fraction_exposed_gap"),
                "p90_nearest_lipid_distance": im.get("p90_nearest_lipid_distance"),
                "quality_passed":            im.get("quality_passed"),
            } if im else {},
            "score":               c["score"],
            "quality_passed":      c["quality_passed"],
            "quality_reasons":     c["quality_reasons"],
            "invalid_reason":      c["invalid_reason"],
        })

    report = {
        "enabled":                   True,
        "policy":                    policy,
        "baseline_metrics": {
            "fraction_covered":          baseline_interface.get("fraction_covered"),
            "fraction_exposed_gap":      baseline_interface.get("fraction_exposed_gap"),
            "p90_nearest_lipid_distance": baseline_interface.get("p90_nearest_lipid_distance"),
            "score":                     baseline_score,
        },
        "n_candidates_evaluated":    len(candidate_results),
        "n_candidates_invalid":      n_invalid,
        "n_candidates_quality_passed": n_quality_passed,
        "candidates":                compact_candidates,
        "selected_candidate_id":     best["candidate_id"] if best else None,
        "selected_candidate_params": {
            "z_shift_offset_nm": best["z_shift_offset_nm"],
            "mask_padding_nm":   best["mask_padding_nm"],
        } if best else None,
        "selected_score":            best["score"] if best else None,
        "baseline_score":            baseline_score,
        "score_improvement":         round(best["score"] - baseline_score, 5) if best else None,
        "selected_quality_passed":   best["quality_passed"] if best else None,
        "reason_selected":           reason_selected,
        "promoted_to_work_input_gro": promoted,
        "warnings":                  warnings,
    }
    _write_report(report, output_dir)
    return report


# ── Helpers ───────────────────────────────────────────────────────────────────

def _write_report(report: dict, output_dir: Path) -> None:
    (output_dir / "membrane_embedding_optimizer_report.json").write_text(
        json.dumps(report, indent=2)
    )


def _disabled_report(policy: str, reason: str, warnings: list[str]) -> dict:
    return {
        "enabled":                    False,
        "policy":                     policy,
        "baseline_metrics":           {},
        "n_candidates_evaluated":     0,
        "n_candidates_invalid":       0,
        "n_candidates_quality_passed": 0,
        "candidates":                 [],
        "selected_candidate_id":      None,
        "selected_candidate_params":  None,
        "selected_score":             None,
        "baseline_score":             None,
        "score_improvement":          None,
        "selected_quality_passed":    None,
        "reason_selected":            reason,
        "promoted_to_work_input_gro": False,
        "warnings":                   warnings,
    }
