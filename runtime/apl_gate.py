# runtime/apl_gate.py
"""
Gate: APL convergence check after membrane_embedding shrink loop.

Reads shrink_telemetry.json (written by the shrink loop bash script).
The shrink loop already exits 1 on non-convergence, so the gate is
informational — it emits a warning if APL is only barely within tolerance.
"""
from __future__ import annotations

import json
from pathlib import Path

from runtime.gate_runner import GateResult

_NEAR_TOLERANCE_FRACTION = 0.25   # warn if APL is within 25% of the tolerance band


def evaluate_apl_gate(step_dir: Path) -> GateResult | None:
    """
    Returns None when telemetry is absent (shrink loop did not run).
    Never blocks (script exits 1 on non-convergence, executor handles that).
    Warns when APL is very close to the upper tolerance limit.
    """
    telemetry_path = step_dir / "shrink_telemetry.json"
    if not telemetry_path.exists():
        return None
    try:
        data = json.loads(telemetry_path.read_text())
    except Exception:
        return None

    converged    = data.get("converged", False)
    final_apl    = data.get("final_apl_ang2")
    n_iterations = data.get("n_iterations", 0)

    # Read target + tolerance from metadata.json if available
    target_apl  = None
    tolerance   = None
    meta_path   = step_dir / "metadata.json"
    if meta_path.exists():
        try:
            meta = json.loads(meta_path.read_text())
            params = meta.get("params", {})
            target_apl = params.get("apl_target_ang2")
            tolerance  = params.get("apl_tolerance_ang2")
        except Exception:
            pass

    if not converged:
        # Script already failed, but if somehow we reach the gate, block.
        msg = f"Shrink loop did not converge (final APL={final_apl} Å², {n_iterations} iterations)"
        return GateResult(
            passed=False, blocked=True, confidence=1.0,
            errors=[msg],
        )

    warnings: list[str] = []
    if final_apl is not None and target_apl is not None and tolerance is not None:
        cutoff = target_apl + tolerance
        gap    = cutoff - final_apl
        if gap < tolerance * _NEAR_TOLERANCE_FRACTION:
            warnings.append(
                f"APL={final_apl:.1f} Å² is very close to tolerance limit "
                f"({cutoff:.1f} Å²) — consider additional deflation iterations"
            )

    # ── Trapped lipid diagnosis (Phase 4) ─────────────────────────────────────
    lipid_resname = "DPP"
    tm_residues = None
    meta_path = step_dir / "metadata.json"
    if meta_path.exists():
        try:
            meta = json.loads(meta_path.read_text())
            params = meta.get("params", {})
            lipid_resname = params.get("lipid_residue_name", "DPP")
        except Exception:
            pass

    embed_meta_path = step_dir.parent / "03_embed_in_bilayer" / "metadata.json"
    if embed_meta_path.exists():
        try:
            embed_meta = json.loads(embed_meta_path.read_text())
            tm_str = embed_meta.get("params", {}).get("tm_residues", "")
            if tm_str:
                tm_residues = set()
                for part in tm_str.split(','):
                    if '-' in part:
                        start, end = part.split('-')
                        tm_residues.update(range(int(start), int(end) + 1))
                    else:
                        tm_residues.add(int(part))
        except Exception:
            pass

    gro_path = step_dir / "converged.gro"
    if gro_path.exists():
        try:
            from validators.membrane_validators import detect_trapped_lipids
            diag = detect_trapped_lipids(gro_path, lipid_resname=lipid_resname, tm_residues=tm_residues)
            
            run_dir = step_dir.parent.parent
            gro_03 = run_dir / "steps/03_embed_in_bilayer/system.gro"
            gro_05 = step_dir / "converged.gro"
            gro_06 = run_dir / "steps/06_solvate_membrane/solvated.gro"
            gro_07 = run_dir / "steps/07_clean_water/system_clean.gro"
            
            files_inspected = [str(gro_03.resolve()) if gro_03.exists() else str(gro_03)]
            files_inspected.append(str(gro_05.resolve()))
            if gro_06.exists():
                files_inspected.append(str(gro_06.resolve()))
            if gro_07.exists():
                files_inspected.append(str(gro_07.resolve()))

            report_data = {
                "workspace_path": str(run_dir.resolve()),
                "files_inspected": files_inspected,
                "actual_coordinate_file_used_after_membrane_embedding": str(gro_05.resolve()),
                "actual_coordinate_file_passed_to_solvate_membrane": str(gro_05.resolve()),
                "bilayer_midplane_z": diag.bilayer_midplane_z,
                "protein_tm_center": diag.protein_tm_center,
                "protein_xy_footprint": diag.protein_xy_footprint,
                "lipid_residue_candidates": [
                    {
                        "resid": c.resid,
                        "resname": c.resname,
                        "stage_first_detected": c.stage_first_detected,
                        "present_since_embed": c.present_since_embed,
                        "introduced_by_shrink": c.introduced_by_shrink,
                        "inside_tm_footprint": c.inside_tm_footprint,
                        "inside_expanded_tm_footprint": c.inside_expanded_tm_footprint,
                        "low_overlap_but_geometrically_inside": c.low_overlap_but_geometrically_inside,
                        "min_distance_to_protein_nm": c.min_distance_to_protein_nm,
                        "atoms_within_0.25_nm": c.atoms_within_025,
                        "atoms_within_0.35_nm": c.atoms_within_035,
                        "atoms_within_0.45_nm": c.atoms_within_045,
                        "lipid_center": c.lipid_center,
                        "protein_tm_center": c.protein_tm_center,
                        "z_relative_to_midplane": c.z_relative_to_midplane,
                        "suspected_cause": c.suspected_cause,
                        "recommended_action": c.recommended_action
                    }
                    for c in diag.lipid_residue_candidates
                ],
                "summary": {
                    "n_lipids_checked": diag.n_lipids_checked,
                    "n_suspicious_lipids": diag.n_suspicious_lipids,
                    "n_hard_overlap_lipids": diag.n_hard_overlap_lipids,
                    "n_low_overlap_geometrically_trapped_lipids": diag.n_low_overlap_geometrically_trapped_lipids,
                    "likely_root_cause": diag.likely_root_cause,
                    "should_block_future_workflows": diag.should_block_future_workflows,
                    "remediation_not_applied": diag.remediation_not_applied
                }
            }

            report_path = step_dir / "trapped_lipid_diagnosis.json"
            report_path.write_text(json.dumps(report_data, indent=4))
            
            if diag.n_suspicious_lipids > 0:
                warnings.append(
                    f"[trapped_lipids] WARNING: {diag.n_suspicious_lipids} suspicious lipid(s) "
                    f"detected trapped in protein cavity or clashing after membrane embedding. "
                    f"See trapped_lipid_diagnosis.json for details."
                )
        except Exception as e:
            # Non-blocking: fail-safe
            warnings.append(f"[trapped_lipids] Failed to run diagnosis check: {e}")

    apl_str = f"{final_apl:.1f} Å²" if final_apl is not None else "unknown"
    msg = f"Shrink loop converged after {n_iterations} iterations (APL={apl_str})"

    return GateResult(
        passed=True, blocked=False, confidence=1.0,
        warnings=warnings,
    )

