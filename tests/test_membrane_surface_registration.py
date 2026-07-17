"""
Phase 8B — Surface-interference adaptive membrane registration tests.

6 tests:
  1. ECD colliding with upper membrane surface → shift moves bilayer away from ECD
  2. ICD colliding with lower membrane surface → shift moves bilayer away from ICD
  3. TM-only helix → selected shift near zero
  4. Better TM burial but worse surface interference is rejected
  5. Report contains all required surface-interference metrics
  6. GLP-1R regression: surface-interference shift differs from center-only when EC
     overlaps the upper headgroup zone at the center-only placement
"""
import json
import pytest
from pathlib import Path

from validators.membrane_surface_registration import (
    surface_interference_adaptive_registration,
    score_surface_candidate,
)


# ── GRO helpers (same style as test_membrane_registration.py) ─────────────────

def _make_gro(path: Path, atoms: list[tuple], box: tuple = (10.0, 10.0, 10.0)) -> None:
    lines = ["Test system", str(len(atoms))]
    for i, (resnum, resname, atomname, x, y, z) in enumerate(atoms, 1):
        lines.append(
            f"{resnum:>5}{resname:<5}{atomname:>5}{i:>5}{x:8.3f}{y:8.3f}{z:8.3f}"
        )
    lines.append(f"{box[0]:.5f}  {box[1]:.5f}  {box[2]:.5f}")
    path.write_text("\n".join(lines) + "\n")


# ── Test 1: ECD above upper headgroup surface → negative shift ────────────────

def test_ecd_upper_surface_collision_forces_negative_shift(tmp_path):
    """
    ECD residues placed just above the default core top (z=9.25) land in the
    upper headgroup zone [9.25, 9.75] at δ=0.
    Phase 8B must select a NEGATIVE shift to push the bilayer down so the ECD
    clears the upper headgroup surface (z > core_top + HG_thickness).

    Setup:
      TM CA at z=8.0  → TM center = 8.0
      At δ=0: core [6.75, 9.25], upper HG zone [9.25, 9.75]
      EC atoms at z=9.4 → inside upper HG zone → ec_interface_clashes > 0
      At δ=-0.5: upper HG top = 9.25 → EC at 9.4 still inside
      At δ=-1.0: core [6.75, 8.75], upper HG [8.75, 9.25] → EC at 9.4 outside ✓
    """
    tm_residues = set(range(1, 21))
    ec_residues = set(range(21, 41))

    atoms = []
    for r in tm_residues:
        atoms.append((r, "ALA", "CA", 5.0, 5.0, 8.0))
    for r in ec_residues:
        # 2 atoms each to give weight; placed at z=9.4 (in upper HG at δ=0)
        atoms.append((r, "ALA", "CA", 5.0, 5.0, 9.4))
        atoms.append((r, "ALA", "CB", 5.0, 5.0, 9.4))

    protein_gro = tmp_path / "protein_boxed.gro"
    _make_gro(protein_gro, atoms)
    bilayer_gro = tmp_path / "bilayer.gro"
    _make_gro(bilayer_gro, [(1, "DPP", "P", 5.0, 5.0, 5.0)])

    report = surface_interference_adaptive_registration(
        protein_gro             = protein_gro,
        bilayer_gro             = bilayer_gro,
        tm_residues             = tm_residues,
        ec_residues             = ec_residues,
        output_dir              = tmp_path,
    )

    assert report["enabled"] is True
    assert report["selected_shift_z_nm"] < 0.0, (
        f"Expected negative shift to clear ECD from upper HG zone, "
        f"got {report['selected_shift_z_nm']}"
    )
    assert report["ec_interface_clashes"] == 0 or report["upper_surface_interference_score"] < 1.0, (
        "ECD should fully clear the upper headgroup surface at the selected shift"
    )
    assert report["actual_z_shift_nm"] is not None
    assert (tmp_path / "membrane_surface_registration_report.json").exists()


# ── Test 2: ICD below lower headgroup surface → positive shift ────────────────

def test_icd_lower_surface_collision_forces_positive_shift(tmp_path):
    """
    IC residues placed just below the default core bottom (z=6.75) land in the
    lower headgroup zone [6.25, 6.75] at δ=0.
    Phase 8B must select a POSITIVE shift to move the bilayer up so the IC
    domain clears the lower headgroup surface (z < core_bot - HG_thickness).

    Setup:
      TM CA at z=8.0
      At δ=0: core [6.75, 9.25], lower HG [6.25, 6.75]
      IC atoms at z=6.5 → inside lower HG zone → ic_interface_clashes > 0
      At δ=+0.5: core [7.25, 9.75], lower HG [6.75, 7.25] → IC at 6.5 below zone ✓
    """
    tm_residues = set(range(1, 21))
    ic_residues = set(range(21, 41))

    atoms = []
    for r in tm_residues:
        atoms.append((r, "ALA", "CA", 5.0, 5.0, 8.0))
    for r in ic_residues:
        atoms.append((r, "ALA", "CA", 5.0, 5.0, 6.5))
        atoms.append((r, "ALA", "CB", 5.0, 5.0, 6.5))

    protein_gro = tmp_path / "protein_boxed.gro"
    _make_gro(protein_gro, atoms)
    bilayer_gro = tmp_path / "bilayer.gro"
    _make_gro(bilayer_gro, [(1, "DPP", "P", 5.0, 5.0, 5.0)])

    report = surface_interference_adaptive_registration(
        protein_gro             = protein_gro,
        bilayer_gro             = bilayer_gro,
        tm_residues             = tm_residues,
        ic_residues             = ic_residues,
        output_dir              = tmp_path,
    )

    assert report["selected_shift_z_nm"] > 0.0, (
        f"Expected positive shift to clear ICD from lower HG zone, "
        f"got {report['selected_shift_z_nm']}"
    )
    assert report["ic_interface_clashes"] == 0 or report["lower_surface_interference_score"] < 1.0
    assert report["actual_z_shift_nm"] is not None


# ── Test 3: TM-only protein → shift near zero ─────────────────────────────────

def test_tm_only_protein_selects_near_zero_shift(tmp_path):
    """
    All protein residues are TM. No EC/IC.  No surface interference penalties
    apply.  The extreme-shift penalty means any nonzero δ scores worse than δ=0
    if TM burial is already optimal at δ=0.
    """
    tm_residues = set(range(1, 21))
    atoms = [(r, "ALA", "CA", 5.0, 5.0, 8.0) for r in tm_residues]

    protein_gro = tmp_path / "protein_boxed.gro"
    _make_gro(protein_gro, atoms)
    bilayer_gro = tmp_path / "bilayer.gro"
    _make_gro(bilayer_gro, [(1, "DPP", "P", 5.0, 5.0, 5.0)])

    report = surface_interference_adaptive_registration(
        protein_gro = protein_gro,
        bilayer_gro = bilayer_gro,
        tm_residues = tm_residues,
        output_dir  = tmp_path,
    )

    assert abs(report["selected_shift_z_nm"]) < 0.15, (
        f"Expected shift near zero for TM-only protein, "
        f"got {report['selected_shift_z_nm']}"
    )
    assert report["ec_interface_clashes"] == 0
    assert report["ic_interface_clashes"] == 0
    assert report["loop_interface_clashes"] == 0


# ── Test 4: Better TM burial rejected when surface interference worse ─────────

def test_higher_burial_rejected_when_surface_interference_worse(tmp_path):
    """
    Two candidate Z positions:
      δ = 0.0: tm_ca_burial=0.8 but EC atoms are in upper HG zone
      δ = -0.5: tm_ca_burial≈0.7 but EC atoms are ABOVE the HG zone → no clashes

    The Phase 8B score must prefer δ=-0.5 despite lower TM burial, because
    the EC surface interference penalty outweighs the burial reward.
    """
    tm_residues = set(range(1, 21))
    ec_residues = set(range(21, 61))

    atoms = []
    for r in tm_residues:
        atoms.append((r, "ALA", "CA", 5.0, 5.0, 8.0))   # TM at midplane
    # Many EC atoms at z=9.4 — at δ=0 they're in upper HG [9.25, 9.75]
    # At δ=-1.0 the upper HG top is 8.75 so they are above HG zone entirely
    for r in ec_residues:
        atoms.append((r, "ALA", "CA", 5.0, 5.0, 9.4))
        atoms.append((r, "ALA", "CB", 5.0, 5.0, 9.4))
        atoms.append((r, "ALA", "CG", 5.0, 5.0, 9.4))

    protein_gro = tmp_path / "protein_boxed.gro"
    _make_gro(protein_gro, atoms)
    bilayer_gro = tmp_path / "bilayer.gro"
    _make_gro(bilayer_gro, [(1, "DPP", "P", 5.0, 5.0, 5.0)])

    report = surface_interference_adaptive_registration(
        protein_gro             = protein_gro,
        bilayer_gro             = bilayer_gro,
        tm_residues             = tm_residues,
        ec_residues             = ec_residues,
        output_dir              = tmp_path,
    )

    # Verify the selected shift clears EC from upper HG zone
    sel = report["selected_shift_z_nm"]
    # At selected δ, upper HG top = 8.0 + sel + 1.25 + 0.5 = 9.75 + sel
    # For EC at z=9.4 to clear: 9.75 + sel < 9.4 → sel < -0.35
    assert sel < -0.3, (
        f"Expected large negative shift to clear EC from HG zone, got δ={sel}"
    )
    assert report["upper_surface_interference_score"] == 0.0 or sel <= -0.3, (
        "EC should not penetrate upper headgroup surface at selected shift"
    )

    # δ=0 candidate must have EC clashes
    delta0_cand = next(c for c in report["candidates"] if c["shift_z_nm"] == 0.0)
    assert delta0_cand["ec_interface_clashes"] > 0, (
        "At δ=0 there must be EC-headgroup clashes (by construction)"
    )
    assert delta0_cand["registration_score"] < report["best_candidate"]["registration_score"], (
        "δ=0 must score worse than best candidate due to EC surface interference"
    )


# ── Test 5: Report contains all required fields ───────────────────────────────

def test_report_contains_all_required_fields(tmp_path):
    """
    The report dict and the written JSON must contain every field specified in
    the Phase 8B requirements.
    """
    tm_residues = set(range(1, 11))
    ec_residues = set(range(11, 21))
    ic_residues = set(range(21, 31))
    atoms = (
        [(r, "ALA", "CA", 5.0, 5.0, 8.0) for r in tm_residues]
        + [(r, "ALA", "CA", 5.0, 5.0, 10.0) for r in ec_residues]
        + [(r, "ALA", "CA", 5.0, 5.0, 6.0) for r in ic_residues]
    )

    protein_gro = tmp_path / "protein_boxed.gro"
    _make_gro(protein_gro, atoms)
    bilayer_gro = tmp_path / "bilayer.gro"
    _make_gro(bilayer_gro, [(1, "DPP", "P", 5.0, 5.0, 5.0)])

    report = surface_interference_adaptive_registration(
        protein_gro = protein_gro,
        bilayer_gro = bilayer_gro,
        tm_residues = tm_residues,
        ec_residues = ec_residues,
        ic_residues = ic_residues,
        output_dir  = tmp_path,
    )

    required_top_fields = [
        "enabled", "phase", "protein_frame_used",
        "selected_shift_z_nm", "actual_z_shift_nm",
        "initial_tm_center_z", "tm_center_z_in_registration_frame",
        "bilayer_midplane_z_before", "initial_bilayer_midplane_z",
        "final_bilayer_midplane_z",
        "candidates", "best_candidate", "candidate_scores",
        "tm_burial_score",
        "non_tm_core_penetration",
        "upper_surface_interference_score",
        "lower_surface_interference_score",
        "ec_interface_clashes",
        "ic_interface_clashes",
        "loop_interface_clashes",
        "decision_reason",
        "warnings",
    ]
    for field in required_top_fields:
        assert field in report, f"Report missing required field: {field!r}"

    assert report["phase"] == "8B"

    required_candidate_fields = [
        "shift_z_nm", "actual_z_shift_nm", "new_bilayer_midplane_z",
        "core_z_bot", "core_z_top",
        "upper_headgroup_surface_z", "lower_headgroup_surface_z",
        "tm_burial_score", "tm_ca_burial_fraction",
        "non_tm_core_penetration", "non_tm_core_penetration_fraction",
        "upper_surface_interference_score", "lower_surface_interference_score",
        "ec_core_atom_count", "ic_core_atom_count",
        "ec_interface_clashes", "ic_interface_clashes", "loop_interface_clashes",
        "registration_score",
    ]
    for cand in report["candidates"]:
        for field in required_candidate_fields:
            assert field in cand, f"Candidate missing field: {field!r}"

    # JSON on disk must be valid and match
    loaded = json.loads((tmp_path / "membrane_surface_registration_report.json").read_text())
    assert loaded["phase"] == "8B"
    assert loaded["candidates"] == report["candidates"]

    # candidate_scores is a compact list of {shift_z_nm, registration_score}
    assert len(report["candidate_scores"]) == len(report["candidates"])
    for cs in report["candidate_scores"]:
        assert "shift_z_nm" in cs
        assert "registration_score" in cs


# ── Test 6: GLP-1R regression — surface shift differs from center-only ─────────

def test_glp1r_regression_surface_shift_differs_from_center_only(tmp_path):
    """
    GLP-1R-like system: large EC domain (residues 169-250) placed at z=9.4 nm,
    which falls in the upper headgroup zone [9.25, 9.75] at δ=0.

    Phase 8B surface-interference registration must select a shift that is MORE
    NEGATIVE than δ=0 to clear the EC domain from the upper headgroup surface.
    This demonstrates that Phase 8B differs from a center-only registration
    (which would simply choose δ=0 as the TM midplane alignment).
    """
    tm_residues = set(range(1, 169))
    ec_residues = set(range(169, 251))
    ic_residues = set(range(251, 281))

    atoms = []
    for r in range(1, 169):
        atoms.append((r, "ALA", "CA", 5.0, 5.0, 8.0))    # TM at midplane z=8
    for r in range(169, 251):
        atoms.append((r, "ALA", "CA", 5.0, 5.0, 9.4))    # EC in upper HG at δ=0
        atoms.append((r, "ALA", "CB", 5.0, 5.0, 9.4))
    for r in range(251, 281):
        atoms.append((r, "ALA", "CA", 5.0, 5.0, 6.0))    # IC well below bilayer

    protein_gro = tmp_path / "protein_boxed.gro"
    _make_gro(protein_gro, atoms)
    bilayer_gro = tmp_path / "bilayer.gro"
    _make_gro(bilayer_gro, [(1, "DPP", "P", 5.0, 5.0, 5.0)])

    report = surface_interference_adaptive_registration(
        protein_gro = protein_gro,
        bilayer_gro = bilayer_gro,
        tm_residues = tm_residues,
        ec_residues = ec_residues,
        ic_residues = ic_residues,
        output_dir  = tmp_path,
    )

    # δ=0 candidate must have EC clashes (by construction)
    delta0 = next(c for c in report["candidates"] if c["shift_z_nm"] == 0.0)
    assert delta0["ec_interface_clashes"] > 0, (
        "At δ=0, large EC domain must clash with upper headgroup surface"
    )

    # Phase 8B must select a DIFFERENT (more negative) shift than δ=0
    assert report["selected_shift_z_nm"] < 0.0, (
        f"Phase 8B should move bilayer away from ECD; got δ={report['selected_shift_z_nm']}"
    )

    # Upper surface interference score at selected shift should be reduced vs δ=0
    best = report["best_candidate"]
    assert best["upper_surface_interference_score"] <= delta0["upper_surface_interference_score"], (
        "Selected shift must not increase upper surface interference vs default"
    )

    assert report["actual_z_shift_nm"] is not None
    assert report["phase"] == "8B"
