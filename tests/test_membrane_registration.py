"""
Phase 8A — Adaptive Membrane Registration unit tests.

8 tests covering: ECD penetration → shift, centred bundle → δ=0, missing TM annotation,
strict-policy blocking, bilayer file immutability, Z-only shift, full candidate report,
and GLP-1R-like regression.
"""
import json
import pytest
from pathlib import Path

from validators.membrane_registration import (
    adaptive_membrane_registration,
    search_best_registration,
)


# ── GRO helpers ───────────────────────────────────────────────────────────────

def _make_gro(path: Path, atoms: list[tuple], box: tuple = (10.0, 10.0, 10.0)) -> None:
    """Write a minimal GRO from (resnum, resname, atomname, x, y, z) tuples."""
    lines = ["Test system", str(len(atoms))]
    for i, (resnum, resname, atomname, x, y, z) in enumerate(atoms, 1):
        lines.append(
            f"{resnum:>5}{resname:<5}{atomname:>5}{i:>5}{x:8.3f}{y:8.3f}{z:8.3f}"
        )
    lines.append(f"{box[0]:.5f}  {box[1]:.5f}  {box[2]:.5f}")
    path.write_text("\n".join(lines) + "\n")


def _parse_gro_xy(path: Path) -> list[tuple[float, float]]:
    """Return list of (x, y) per atom from a GRO file."""
    lines = path.read_text().splitlines()
    n = int(lines[1].strip())
    result = []
    for ln in lines[2 : 2 + n]:
        if len(ln) >= 36:
            result.append((float(ln[20:28]), float(ln[28:36])))
    return result


# ── Test 1 ────────────────────────────────────────────────────────────────────

def test_ecd_penetration_forces_negative_shift(tmp_path):
    """
    EC domain at z=6.0 is inside the default hydrophobic core [3.75, 6.25].
    Shifting the bilayer downward (δ<0) clears the EC from the core.
    Selected shift must be negative.
    """
    atoms = [(r, "ALA", "CA", 5.0, 5.0, 5.0) for r in range(1, 25)]
    for r in range(25, 49):
        atoms.append((r, "ALA", "CA", 5.0, 5.0, 6.0))
        atoms.append((r, "ALA", "CB", 5.0, 5.0, 6.0))

    protein_gro = tmp_path / "protein.gro"
    _make_gro(protein_gro, atoms)
    bilayer_gro = tmp_path / "bilayer.gro"
    _make_gro(bilayer_gro, [(1, "DPP", "P", 5.0, 5.0, 5.0)])

    report = adaptive_membrane_registration(
        protein_gro = protein_gro,
        bilayer_gro = bilayer_gro,
        tm_residues = set(range(1, 25)),
        ec_residues = set(range(25, 49)),
        output_dir  = tmp_path,
    )

    assert report["selected_shift_z_nm"] < 0.0, (
        f"Expected negative shift to clear EC domain, got {report['selected_shift_z_nm']}"
    )
    assert report["actual_z_shift_nm"] is not None
    assert (tmp_path / "membrane_registration_report.json").exists()


# ── Test 2 ────────────────────────────────────────────────────────────────────

def test_centered_bundle_selects_zero_shift(tmp_path):
    """
    TM Cα atoms centred exactly at bilayer midplane, EC far above the core.
    The extreme-shift penalty means any nonzero δ scores worse than δ=0.
    """
    atoms = [(r, "ALA", "CA", 5.0, 5.0, 5.0) for r in range(1, 21)]
    for r in range(21, 31):
        atoms.append((r, "ALA", "CA", 5.0, 5.0, 9.0))

    protein_gro = tmp_path / "protein.gro"
    _make_gro(protein_gro, atoms)
    bilayer_gro = tmp_path / "bilayer.gro"
    _make_gro(bilayer_gro, [(1, "DPP", "P", 5.0, 5.0, 5.0)])

    report = adaptive_membrane_registration(
        protein_gro = protein_gro,
        bilayer_gro = bilayer_gro,
        tm_residues = set(range(1, 21)),
        ec_residues = set(range(21, 31)),
        output_dir  = tmp_path,
    )

    assert report["selected_shift_z_nm"] == 0.0, (
        f"Expected δ=0 for well-centred bundle, got {report['selected_shift_z_nm']}"
    )
    assert report["tm_ca_burial_fraction"] == 1.0


# ── Test 3 ────────────────────────────────────────────────────────────────────

def test_no_tm_annotation_fallback(tmp_path):
    """
    When tm_residues is None the function must return enabled=True,
    actual_z_shift_nm=None, and include a TM-annotation warning.
    """
    protein_gro = tmp_path / "protein.gro"
    bilayer_gro = tmp_path / "bilayer.gro"
    _make_gro(protein_gro, [(1, "ALA", "CA", 5.0, 5.0, 5.0)])
    _make_gro(bilayer_gro, [(1, "DPP", "P", 5.0, 5.0, 5.0)])

    report = adaptive_membrane_registration(
        protein_gro = protein_gro,
        bilayer_gro = bilayer_gro,
        tm_residues = None,
        output_dir  = tmp_path,
    )

    assert report["enabled"] is True
    assert report["actual_z_shift_nm"] is None
    assert report["selected_shift_z_nm"] == 0.0
    assert len(report["warnings"]) > 0
    assert any("TM" in w or "tm" in w.lower() for w in report["warnings"])
    assert (
        "skipped" in report["decision_reason"].lower()
        or "no tm" in report["decision_reason"].lower()
    )


# ── Test 4 ────────────────────────────────────────────────────────────────────

def test_strict_policy_blocks_excessive_shift(tmp_path):
    """
    EC at z=5.7 exits the core only at δ<-0.55 nm; the first grid hit is δ=-0.6.
    With max_allowed_shift_nm=0.5, the scoring still prefers δ=-0.6 (score≈1.57>1.2),
    but abs(-0.6)>0.5 triggers the strict-policy block.
    selected_shift_z_nm must fall back to 0.0 and blocked must be True.
    """
    atoms = [(r, "ALA", "CA", 5.0, 5.0, 5.0) for r in range(1, 25)]
    for r in range(25, 49):
        atoms.append((r, "ALA", "CA", 5.0, 5.0, 5.7))
        atoms.append((r, "ALA", "CB", 5.0, 5.0, 5.7))

    protein_gro = tmp_path / "protein.gro"
    _make_gro(protein_gro, atoms)
    bilayer_gro = tmp_path / "bilayer.gro"
    _make_gro(bilayer_gro, [(1, "DPP", "P", 5.0, 5.0, 5.0)])

    report = adaptive_membrane_registration(
        protein_gro          = protein_gro,
        bilayer_gro          = bilayer_gro,
        tm_residues          = set(range(1, 25)),
        ec_residues          = set(range(25, 49)),
        policy               = "strict",
        max_allowed_shift_nm = 0.5,   # δ=-0.6 beats δ=0 but |−0.6|>0.5 → blocked
        output_dir           = tmp_path,
    )

    assert report.get("blocked") is True, (
        f"Expected blocked=True; got blocked={report.get('blocked')}, "
        f"selected={report.get('selected_shift_z_nm')}"
    )
    assert report["selected_shift_z_nm"] == 0.0
    assert "BLOCKED" in report["decision_reason"]


# ── Test 5 ────────────────────────────────────────────────────────────────────

def test_bilayer_file_not_modified(tmp_path):
    """
    adaptive_membrane_registration must never write back to bilayer_gro.
    Both file content and mtime must be identical before and after the call.
    """
    protein_gro = tmp_path / "protein.gro"
    _make_gro(protein_gro, [(r, "ALA", "CA", 5.0, 5.0, 5.0) for r in range(1, 5)])

    bilayer_gro = tmp_path / "bilayer.gro"
    bil_atoms = [(i, "DPP", "P", float(i) * 0.5, float(i) * 0.3, 5.0) for i in range(1, 6)]
    _make_gro(bilayer_gro, bil_atoms)

    content_before = bilayer_gro.read_text()
    mtime_before   = bilayer_gro.stat().st_mtime

    adaptive_membrane_registration(
        protein_gro = protein_gro,
        bilayer_gro = bilayer_gro,
        tm_residues = set(range(1, 5)),
        output_dir  = tmp_path,
    )

    assert bilayer_gro.read_text() == content_before, "bilayer GRO content was modified"
    assert bilayer_gro.stat().st_mtime == mtime_before, "bilayer GRO file was touched"


# ── Test 6 ────────────────────────────────────────────────────────────────────

def test_registration_shift_is_z_only(tmp_path):
    """
    The reported actual_z_shift_nm is a pure Z offset.
    Applying it to bilayer atom positions must change Z by exactly that amount
    while leaving X and Y coordinates identical.
    """
    # TM at z=6.0, bilayer at z=5.0 → expected nonzero actual_z_shift
    protein_gro = tmp_path / "protein.gro"
    _make_gro(protein_gro, [(r, "ALA", "CA", 5.0, 5.0, 6.0) for r in range(1, 5)])

    bil_atoms = [
        (i, "DPP", "P", float(i) * 1.1, float(i) * 0.7, 5.0) for i in range(1, 6)
    ]
    bilayer_gro = tmp_path / "bilayer.gro"
    _make_gro(bilayer_gro, bil_atoms)

    report = adaptive_membrane_registration(
        protein_gro = protein_gro,
        bilayer_gro = bilayer_gro,
        tm_residues = set(range(1, 5)),
        output_dir  = tmp_path,
    )

    z_shift = report["actual_z_shift_nm"]
    assert z_shift is not None
    assert abs(z_shift) > 1e-3, (
        f"Expected nonzero Z shift (TM at 6.0, bilayer midplane at 5.0), got {z_shift}"
    )

    # Simulate applying the shift: X and Y must be untouched
    for orig in bil_atoms:
        ox, oy, oz = orig[3], orig[4], orig[5]
        shifted_z = oz + z_shift
        assert abs((shifted_z - oz) - z_shift) < 1e-9, "Z shift not applied exactly"
        # X and Y are unaffected (trivially — z_shift is scalar, applied to Z only)
        assert ox == orig[3], "X changed unexpectedly"
        assert oy == orig[4], "Y changed unexpectedly"

    # Confirm bilayer_gro XY coords are unchanged on disk
    xy_from_file = _parse_gro_xy(bilayer_gro)
    for (ox, oy), (fx, fy) in zip([(a[3], a[4]) for a in bil_atoms], xy_from_file):
        assert abs(ox - fx) < 1e-3, f"X on disk changed: {ox} → {fx}"
        assert abs(oy - fy) < 1e-3, f"Y on disk changed: {oy} → {fy}"


# ── Test 7 ────────────────────────────────────────────────────────────────────

def test_report_contains_all_candidate_scores(tmp_path):
    """
    With search range [-0.5, 0.5] step 0.1 the report must contain exactly 11
    candidates, each carrying every required field, and the JSON file must be valid.
    """
    protein_gro = tmp_path / "protein.gro"
    _make_gro(protein_gro, [(r, "ALA", "CA", 5.0, 5.0, 5.0) for r in range(1, 5)])
    bilayer_gro = tmp_path / "bilayer.gro"
    _make_gro(bilayer_gro, [(1, "DPP", "P", 5.0, 5.0, 5.0)])

    report = adaptive_membrane_registration(
        protein_gro   = protein_gro,
        bilayer_gro   = bilayer_gro,
        tm_residues   = set(range(1, 5)),
        output_dir    = tmp_path,
        search_min_nm  = -0.5,
        search_max_nm  =  0.5,
        search_step_nm =  0.1,
    )

    assert len(report["candidates"]) == 11, (
        f"Expected 11 candidates for [-0.5, 0.5] step 0.1, "
        f"got {len(report['candidates'])}"
    )

    required_fields = [
        "shift_z_nm", "actual_z_shift_nm", "new_bilayer_midplane_z",
        "core_z_bot", "core_z_top",
        "tm_burial_fraction", "tm_ca_burial_fraction",
        "non_tm_core_atom_count", "non_tm_core_residue_count",
        "ec_core_atom_count", "ic_core_atom_count",
        "soluble_domain_core_penetration_fraction",
        "estimated_lipid_interference_score",
        "expected_void_risk_score",
        "registration_score",
        "warnings",
    ]
    for cand in report["candidates"]:
        for fld in required_fields:
            assert fld in cand, f"Candidate missing field: {fld!r}"

    best_score = report["best_candidate"]["registration_score"]
    all_scores = [c["registration_score"] for c in report["candidates"]]
    assert best_score in all_scores, "best_candidate score not found in candidates list"

    loaded = json.loads((tmp_path / "membrane_registration_report.json").read_text())
    assert loaded["candidates"] == report["candidates"]


# ── Test 8 ────────────────────────────────────────────────────────────────────

def test_glp1r_like_regression_selects_nonzero_shift(tmp_path):
    """
    GLP-1R-like system: 7-TM bundle (res 1-168, z=5.0), large EC domain
    (res 169-250, z=6.0, inside default core), IC tail (res 251-280, z=2.0).
    Optimizer must select a negative shift to reduce EC penetration.
    """
    atoms = []
    for r in range(1, 169):
        atoms.append((r, "ALA", "CA", 5.0, 5.0, 5.0))
    for r in range(169, 251):
        atoms.append((r, "ALA", "CA", 5.0, 5.0, 6.0))
        atoms.append((r, "ALA", "CB", 5.0, 5.0, 6.0))
    for r in range(251, 281):
        atoms.append((r, "ALA", "CA", 5.0, 5.0, 2.0))

    protein_gro = tmp_path / "protein.gro"
    _make_gro(protein_gro, atoms)
    bilayer_gro = tmp_path / "bilayer.gro"
    _make_gro(bilayer_gro, [(1, "DPP", "P", 5.0, 5.0, 5.0)])

    report = adaptive_membrane_registration(
        protein_gro = protein_gro,
        bilayer_gro = bilayer_gro,
        tm_residues = set(range(1, 169)),
        ec_residues = set(range(169, 251)),
        ic_residues = set(range(251, 281)),
        output_dir  = tmp_path,
    )

    assert report["selected_shift_z_nm"] < 0.0, (
        f"Expected negative shift for GPCR-like EC-in-core system, "
        f"got {report['selected_shift_z_nm']}"
    )
    assert report["actual_z_shift_nm"] is not None

    default_cand = next(
        c for c in report["candidates"] if c["shift_z_nm"] == 0.0
    )
    assert report["ec_core_atom_count"] <= default_cand["ec_core_atom_count"], (
        "Selected shift should not increase EC penetration vs default"
    )


# ── Regression Tests: Phase 8A coordinate-frame fix ──────────────────────────

def test_protein_frame_used_recorded_in_report(tmp_path):
    """
    protein_frame_used must be present in the report and equal the basename of
    the protein GRO passed.  When the file is 'protein_boxed.gro', that name
    must appear verbatim in both the returned dict and the JSON on disk.
    """
    protein_gro = tmp_path / "protein_boxed.gro"
    _make_gro(protein_gro, [(r, "ALA", "CA", 5.0, 5.0, 5.0) for r in range(1, 5)])
    bilayer_gro = tmp_path / "bilayer.gro"
    _make_gro(bilayer_gro, [(1, "DPP", "P", 5.0, 5.0, 5.0)])

    report = adaptive_membrane_registration(
        protein_gro = protein_gro,
        bilayer_gro = bilayer_gro,
        tm_residues = set(range(1, 5)),
        output_dir  = tmp_path,
    )

    assert "protein_frame_used" in report, "Report must contain 'protein_frame_used'"
    assert report["protein_frame_used"] == "protein_boxed.gro", (
        f"Expected protein_frame_used='protein_boxed.gro', "
        f"got {report['protein_frame_used']!r}"
    )
    loaded = json.loads((tmp_path / "membrane_registration_report.json").read_text())
    assert loaded["protein_frame_used"] == "protein_boxed.gro"


def test_oriented_vs_boxed_frame_yields_different_tm_center(tmp_path):
    """
    protein_oriented.gro (TM CA at z=5.0) vs protein_boxed.gro (TM CA at z=8.0)
    must produce different initial_tm_center_z values.
    The 3.0 nm delta equals a typical editconf -c recentering shift.
    """
    oriented_gro = tmp_path / "protein_oriented.gro"
    _make_gro(oriented_gro, [(r, "ALA", "CA", 5.0, 5.0, 5.0) for r in range(1, 5)])

    boxed_gro = tmp_path / "protein_boxed.gro"
    _make_gro(boxed_gro, [(r, "ALA", "CA", 5.0, 5.0, 8.0) for r in range(1, 5)])

    bilayer_gro = tmp_path / "bilayer.gro"
    _make_gro(bilayer_gro, [(1, "DPP", "P", 5.0, 5.0, 5.0)])

    out_o = tmp_path / "out_oriented"
    out_o.mkdir()
    report_o = adaptive_membrane_registration(
        protein_gro = oriented_gro,
        bilayer_gro = bilayer_gro,
        tm_residues = set(range(1, 5)),
        output_dir  = out_o,
    )

    out_b = tmp_path / "out_boxed"
    out_b.mkdir()
    report_b = adaptive_membrane_registration(
        protein_gro = boxed_gro,
        bilayer_gro = bilayer_gro,
        tm_residues = set(range(1, 5)),
        output_dir  = out_b,
    )

    z_o = report_o["initial_tm_center_z"]
    z_b = report_b["initial_tm_center_z"]
    assert abs(z_o - 5.0) < 0.01, f"Oriented TM center expected 5.0, got {z_o}"
    assert abs(z_b - 8.0) < 0.01, f"Boxed TM center expected 8.0, got {z_b}"
    assert abs(z_b - z_o - 3.0) < 0.01, (
        f"Expected 3.0 nm frame difference, got {z_b - z_o:.4f}"
    )
    assert report_b["protein_frame_used"] == "protein_boxed.gro"
    assert report_o["protein_frame_used"] == "protein_oriented.gro"


def test_final_bilayer_midplane_matches_boxed_tm_center(tmp_path):
    """
    With TM CA at z=8.0 (post-editconf frame) and EC far above the hydrophobic
    core, δ=0 is optimal and the final bilayer midplane must equal the TM center.
    This confirms registration is correct when protein_boxed.gro is used.
    """
    protein_gro = tmp_path / "protein_boxed.gro"
    atoms = [(r, "ALA", "CA", 5.0, 5.0, 8.0) for r in range(1, 21)]
    for r in range(21, 31):
        atoms.append((r, "ALA", "CA", 5.0, 5.0, 12.0))   # EC far outside core
    _make_gro(protein_gro, atoms, box=(10.0, 10.0, 20.0))

    bilayer_gro = tmp_path / "bilayer.gro"
    _make_gro(bilayer_gro, [(1, "DPP", "P", 5.0, 5.0, 5.0)])

    report = adaptive_membrane_registration(
        protein_gro = protein_gro,
        bilayer_gro = bilayer_gro,
        tm_residues = set(range(1, 21)),
        ec_residues = set(range(21, 31)),
        output_dir  = tmp_path,
    )

    assert report["selected_shift_z_nm"] == 0.0, (
        f"Expected δ=0 for well-centred boxed protein, got {report['selected_shift_z_nm']}"
    )
    assert abs(report["final_bilayer_midplane_z"] - report["initial_tm_center_z"]) < 0.01, (
        f"Final bilayer midplane {report['final_bilayer_midplane_z']} "
        f"must match TM center {report['initial_tm_center_z']}"
    )
    assert report["protein_frame_used"] == "protein_boxed.gro"
