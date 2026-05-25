# core/test_orientation_validation.py
"""
Tests for OrientationValidationLayer.

Coverage:
    - EC/IC correctly separated → high confidence, passed
    - EC/IC inverted → error, low confidence
    - EC/IC small Z-separation → warning, penalty
    - EC atoms not found in GRO → error
    - IC atoms not found in GRO → error
    - EC→IC vector too short (degenerate) → warning, penalty
    - Excessive tilt → warning, penalty
    - TM COM outside bilayer zone → warning, penalty
    - TM COM inside bilayer zone → pass
    - No TM residues → TM check skipped
    - to_dict() shape
    - to_markdown() shape
    - extracellular_side="-z" convention
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from core.orientation_validation import (
    validate_orientation,
    OrientationValidationReport,
    CheckResult,
)


# ─────────────────────────────────────────────────────────────────────────────
# GRO builder helper (same format as test_membrane_geometry.py)
# ─────────────────────────────────────────────────────────────────────────────

def _make_gro(atoms: list[tuple[int, str, str, float, float, float]]) -> str:
    """Build a minimal GRO string. atoms = [(resnum, resname, atomname, x, y, z)]."""
    lines = ["Test protein", str(len(atoms))]
    for i, (resnum, resname, atomname, x, y, z) in enumerate(atoms, 1):
        lines.append(
            f"{resnum:5d}{resname:<5s}{atomname:>5s}{i:5d}"
            f"{x:8.3f}{y:8.3f}{z:8.3f}"
        )
    lines.append("10.0 10.0 10.0")
    return "\n".join(lines) + "\n"


def _write_gro(tmp_path: Path, atoms) -> Path:
    p = tmp_path / "protein_oriented.gro"
    p.write_text(_make_gro(atoms))
    return p


# ─────────────────────────────────────────────────────────────────────────────
# Happy-path: correct orientation
# ─────────────────────────────────────────────────────────────────────────────

class TestCorrectOrientation:
    def test_ec_plus_z_ic_minus_z(self, tmp_path):
        gro = _write_gro(tmp_path, [
            (1, "ALA", "CA", 0.0, 0.0,  3.0),   # EC residue, high Z
            (2, "ALA", "CA", 0.0, 0.0,  3.5),   # EC residue
            (5, "ALA", "CA", 0.0, 0.0, -2.5),   # IC residue, low Z
            (6, "ALA", "CA", 0.0, 0.0, -3.0),   # IC residue
        ])
        result = validate_orientation(
            gro,
            ec_residues={1, 2},
            ic_residues={5, 6},
            extracellular_side="+z",
        )
        assert result.passed
        assert result.confidence >= 0.9
        assert result.ec_com_z is not None and result.ec_com_z > 0
        assert result.ic_com_z is not None and result.ic_com_z < 0
        ec_ic_check = next(c for c in result.checks if c.name == "ec_ic_side")
        assert ec_ic_check.passed

    def test_confidence_one_when_all_pass(self, tmp_path):
        """Penalty=0 when all checks pass → confidence=1.0."""
        gro = _write_gro(tmp_path, [
            (1, "ALA", "CA", 0.0, 0.0,  3.0),
            (5, "ALA", "CA", 0.0, 0.0, -3.0),
        ])
        result = validate_orientation(
            gro, ec_residues={1}, ic_residues={5}
        )
        assert result.confidence == 1.0
        assert not result.warnings
        assert not result.errors

    def test_minus_z_convention(self, tmp_path):
        """extracellular_side='-z' → EC at -Z is correct."""
        gro = _write_gro(tmp_path, [
            (1, "ALA", "CA", 0.0, 0.0, -3.0),   # EC at -Z
            (5, "ALA", "CA", 0.0, 0.0,  3.0),   # IC at +Z
        ])
        result = validate_orientation(
            gro,
            ec_residues={1},
            ic_residues={5},
            extracellular_side="-z",
        )
        assert result.passed
        ec_ic_check = next(c for c in result.checks if c.name == "ec_ic_side")
        assert ec_ic_check.passed


# ─────────────────────────────────────────────────────────────────────────────
# Inverted orientation
# ─────────────────────────────────────────────────────────────────────────────

class TestInvertedOrientation:
    def test_ec_at_minus_z_when_plus_z_expected(self, tmp_path):
        gro = _write_gro(tmp_path, [
            (1, "ALA", "CA", 0.0, 0.0, -3.0),   # EC at -Z (wrong!)
            (5, "ALA", "CA", 0.0, 0.0,  3.0),   # IC at +Z (wrong!)
        ])
        result = validate_orientation(
            gro,
            ec_residues={1},
            ic_residues={5},
            extracellular_side="+z",
        )
        assert not result.passed
        assert any("inverted" in e.lower() or "wrong" in e.lower() for e in result.errors)
        assert result.confidence < 0.7

    def test_inversion_detected_via_check_result(self, tmp_path):
        gro = _write_gro(tmp_path, [
            (1, "ALA", "CA", 0.0, 0.0, -2.0),
            (5, "ALA", "CA", 0.0, 0.0,  2.0),
        ])
        result = validate_orientation(
            gro, ec_residues={1}, ic_residues={5}, extracellular_side="+z"
        )
        ec_ic_check = next(c for c in result.checks if c.name == "ec_ic_side")
        assert not ec_ic_check.passed


# ─────────────────────────────────────────────────────────────────────────────
# Small Z-separation
# ─────────────────────────────────────────────────────────────────────────────

class TestSmallSeparation:
    def test_separation_below_threshold(self, tmp_path):
        gro = _write_gro(tmp_path, [
            (1, "ALA", "CA", 0.0, 0.0, 0.2),    # EC barely above IC
            (5, "ALA", "CA", 0.0, 0.0, 0.0),
        ])
        result = validate_orientation(
            gro,
            ec_residues={1},
            ic_residues={5},
            min_ec_ic_separation_nm=0.5,
        )
        assert result.warnings
        assert any("separation" in w.lower() for w in result.warnings)
        assert result.confidence < 1.0

    def test_confidence_reduced_by_small_sep(self, tmp_path):
        gro = _write_gro(tmp_path, [
            (1, "ALA", "CA", 0.0, 0.0, 0.1),
            (5, "ALA", "CA", 0.0, 0.0, 0.0),
        ])
        result_tight = validate_orientation(
            gro, ec_residues={1}, ic_residues={5}
        )
        gro2 = _write_gro(tmp_path, [
            (1, "ALA", "CA", 0.0, 0.0, 3.0),
            (5, "ALA", "CA", 0.0, 0.0, -3.0),
        ])
        result_good = validate_orientation(
            gro2, ec_residues={1}, ic_residues={5}
        )
        assert result_tight.confidence < result_good.confidence


# ─────────────────────────────────────────────────────────────────────────────
# Missing atoms
# ─────────────────────────────────────────────────────────────────────────────

class TestMissingAtoms:
    def test_ec_atoms_not_found(self, tmp_path):
        gro = _write_gro(tmp_path, [
            (5, "ALA", "CA", 0.0, 0.0, -2.0),
        ])
        result = validate_orientation(
            gro, ec_residues={999}, ic_residues={5}
        )
        assert not result.passed
        assert any("EC" in e for e in result.errors)
        assert result.n_ec_atoms == 0

    def test_ic_atoms_not_found(self, tmp_path):
        gro = _write_gro(tmp_path, [
            (1, "ALA", "CA", 0.0, 0.0, 2.0),
        ])
        result = validate_orientation(
            gro, ec_residues={1}, ic_residues={999}
        )
        assert not result.passed
        assert any("IC" in e for e in result.errors)
        assert result.n_ic_atoms == 0

    def test_both_missing_large_penalty(self, tmp_path):
        gro = _write_gro(tmp_path, [
            (1, "ALA", "CB", 0.0, 0.0, 2.0),    # CB, not CA → skipped
        ])
        result = validate_orientation(
            gro, ec_residues={1}, ic_residues={2}
        )
        # Two missing-atom penalties (0.35 each) → confidence ≤ 0.35
        assert result.confidence <= 0.35
        assert len(result.errors) == 2


# ─────────────────────────────────────────────────────────────────────────────
# EC-IC distance (degenerate check)
# ─────────────────────────────────────────────────────────────────────────────

class TestEcIcDistance:
    def test_vector_too_short(self, tmp_path):
        gro = _write_gro(tmp_path, [
            (1, "ALA", "CA", 0.0, 0.0, 0.8),    # EC
            (5, "ALA", "CA", 0.0, 0.0, 0.0),    # IC — 0.8 nm apart (< 1.0)
        ])
        result = validate_orientation(
            gro,
            ec_residues={1},
            ic_residues={5},
            min_ec_ic_distance_nm=1.0,
        )
        assert result.warnings
        dist_check = next(c for c in result.checks if c.name == "ec_ic_distance")
        assert not dist_check.passed
        assert result.ec_ic_distance_nm is not None
        assert result.ec_ic_distance_nm < 1.0

    def test_vector_sufficient(self, tmp_path):
        gro = _write_gro(tmp_path, [
            (1, "ALA", "CA", 0.0, 0.0,  2.0),
            (5, "ALA", "CA", 0.0, 0.0, -2.0),
        ])
        result = validate_orientation(
            gro, ec_residues={1}, ic_residues={5}
        )
        dist_check = next(c for c in result.checks if c.name == "ec_ic_distance")
        assert dist_check.passed
        assert abs(result.ec_ic_distance_nm - 4.0) < 0.01


# ─────────────────────────────────────────────────────────────────────────────
# Tilt angle
# ─────────────────────────────────────────────────────────────────────────────

class TestTiltAngle:
    def test_zero_tilt_when_aligned_to_z(self, tmp_path):
        gro = _write_gro(tmp_path, [
            (1, "ALA", "CA", 0.0, 0.0,  3.0),   # EC directly above IC
            (5, "ALA", "CA", 0.0, 0.0, -3.0),
        ])
        result = validate_orientation(
            gro, ec_residues={1}, ic_residues={5}
        )
        assert result.tilt_angle_deg is not None
        assert result.tilt_angle_deg < 1.0   # ≈ 0°

    def test_large_tilt_warns(self, tmp_path):
        # EC and IC offset mostly in XY plane → high tilt
        gro = _write_gro(tmp_path, [
            (1, "ALA", "CA",  5.0, 0.0, 0.3),   # EC — large XY offset
            (5, "ALA", "CA", -5.0, 0.0, 0.1),   # IC
        ])
        result = validate_orientation(
            gro,
            ec_residues={1},
            ic_residues={5},
            max_tilt_angle_deg=45.0,
        )
        assert result.tilt_angle_deg is not None
        assert result.tilt_angle_deg > 80.0   # nearly horizontal
        tilt_check = next(c for c in result.checks if c.name == "tilt_angle")
        assert not tilt_check.passed
        assert result.warnings

    def test_tilt_exactly_at_threshold_passes(self, tmp_path):
        # Build atoms where tilt ≈ 30° (below 45° threshold)
        # EC at z=2, IC at z=-2, x offset = 2*tan(30°) ≈ 2.309
        import math as _m
        x_off = 2.0 * _m.tan(_m.radians(30))
        gro = _write_gro(tmp_path, [
            (1, "ALA", "CA",  x_off, 0.0,  2.0),
            (5, "ALA", "CA", 0.0,   0.0, -2.0),
        ])
        result = validate_orientation(
            gro, ec_residues={1}, ic_residues={5}, max_tilt_angle_deg=45.0
        )
        assert result.tilt_angle_deg is not None
        assert result.tilt_angle_deg < 35.0
        tilt_check = next(c for c in result.checks if c.name == "tilt_angle")
        assert tilt_check.passed


# ─────────────────────────────────────────────────────────────────────────────
# TM centering
# ─────────────────────────────────────────────────────────────────────────────

class TestTmCentering:
    def test_tm_inside_bilayer_zone(self, tmp_path):
        gro = _write_gro(tmp_path, [
            (1, "ALA", "CA", 0.0, 0.0,  3.0),   # EC
            (5, "ALA", "CA", 0.0, 0.0, -3.0),   # IC
            (3, "ALA", "CA", 0.0, 0.0,  0.2),   # TM — near Z=0
        ])
        result = validate_orientation(
            gro,
            ec_residues={1},
            ic_residues={5},
            tm_residues={3},
            bilayer_half_thickness_nm=2.0,
        )
        tm_check = next(c for c in result.checks if c.name == "tm_centered")
        assert tm_check.passed

    def test_tm_outside_bilayer_zone_warns(self, tmp_path):
        gro = _write_gro(tmp_path, [
            (1, "ALA", "CA", 0.0, 0.0,  4.0),   # EC
            (5, "ALA", "CA", 0.0, 0.0, -1.0),   # IC
            (3, "ALA", "CA", 0.0, 0.0,  3.5),   # TM — far from Z=0
        ])
        result = validate_orientation(
            gro,
            ec_residues={1},
            ic_residues={5},
            tm_residues={3},
            bilayer_half_thickness_nm=2.0,
        )
        tm_check = next(c for c in result.checks if c.name == "tm_centered")
        assert not tm_check.passed
        assert any("TM" in w for w in result.warnings)

    def test_no_tm_residues_skips_check(self, tmp_path):
        gro = _write_gro(tmp_path, [
            (1, "ALA", "CA", 0.0, 0.0,  3.0),
            (5, "ALA", "CA", 0.0, 0.0, -3.0),
        ])
        result = validate_orientation(
            gro, ec_residues={1}, ic_residues={5}, tm_residues=None
        )
        tm_checks = [c for c in result.checks if c.name == "tm_centered"]
        assert not tm_checks   # check not run


# ─────────────────────────────────────────────────────────────────────────────
# Geometry fields
# ─────────────────────────────────────────────────────────────────────────────

class TestGeometryFields:
    def test_n_atoms_counted(self, tmp_path):
        gro = _write_gro(tmp_path, [
            (1, "ALA", "CA", 0.0, 0.0,  3.0),
            (2, "ALA", "CA", 0.0, 0.0,  2.5),   # second EC atom
            (5, "ALA", "CA", 0.0, 0.0, -2.5),
            (3, "ALA", "CA", 0.0, 0.0,  0.1),   # TM
        ])
        result = validate_orientation(
            gro,
            ec_residues={1, 2},
            ic_residues={5},
            tm_residues={3},
        )
        assert result.n_ec_atoms == 2
        assert result.n_ic_atoms == 1
        assert result.n_tm_atoms == 1

    def test_com_values_correct(self, tmp_path):
        gro = _write_gro(tmp_path, [
            (1, "ALA", "CA", 0.0, 0.0, 2.0),
            (2, "ALA", "CA", 0.0, 0.0, 4.0),   # EC COM = 3.0
            (5, "ALA", "CA", 0.0, 0.0, -3.0),  # IC COM = -3.0
        ])
        result = validate_orientation(
            gro, ec_residues={1, 2}, ic_residues={5}
        )
        assert result.ec_com_z is not None
        assert abs(result.ec_com_z - 3.0) < 0.01
        assert result.ic_com_z is not None
        assert abs(result.ic_com_z - (-3.0)) < 0.01
        assert result.ec_ic_distance_nm is not None
        assert abs(result.ec_ic_distance_nm - 6.0) < 0.01


# ─────────────────────────────────────────────────────────────────────────────
# Report serialization
# ─────────────────────────────────────────────────────────────────────────────

class TestReportSerialization:
    def _make_result(self, tmp_path) -> OrientationValidationReport:
        gro = _write_gro(tmp_path, [
            (1, "ALA", "CA", 0.0, 0.0,  3.0),
            (5, "ALA", "CA", 0.0, 0.0, -3.0),
            (3, "ALA", "CA", 0.0, 0.0,  0.5),
        ])
        return validate_orientation(
            gro,
            ec_residues={1},
            ic_residues={5},
            tm_residues={3},
        )

    def test_to_dict_has_required_keys(self, tmp_path):
        d = self._make_result(tmp_path).to_dict()
        assert "confidence" in d
        assert "passed" in d
        assert "geometry" in d
        assert "checks" in d
        assert "warnings" in d
        assert "errors" in d

    def test_to_dict_geometry_fields(self, tmp_path):
        d = self._make_result(tmp_path).to_dict()
        g = d["geometry"]
        assert "ec_com_z_nm" in g
        assert "ic_com_z_nm" in g
        assert "tm_com_z_nm" in g
        assert "ec_ic_distance_nm" in g
        assert "tilt_angle_deg" in g
        assert "n_ec_atoms" in g

    def test_to_dict_confidence_rounded(self, tmp_path):
        d = self._make_result(tmp_path).to_dict()
        assert isinstance(d["confidence"], float)
        assert d["confidence"] <= 1.0

    def test_to_markdown_contains_headers(self, tmp_path):
        md = self._make_result(tmp_path).to_markdown()
        assert "# Orientation Validation Report" in md
        assert "## Geometry" in md
        assert "## Checks" in md
        assert "## Warnings" in md
        assert "## Errors" in md

    def test_to_markdown_contains_pass(self, tmp_path):
        md = self._make_result(tmp_path).to_markdown()
        assert "PASS" in md

    def test_to_markdown_contains_confidence(self, tmp_path):
        md = self._make_result(tmp_path).to_markdown()
        assert "Confidence" in md

    def test_to_markdown_fail_on_inverted(self, tmp_path):
        gro = _write_gro(tmp_path, [
            (1, "ALA", "CA", 0.0, 0.0, -3.0),  # EC at -Z (wrong for +z config)
            (5, "ALA", "CA", 0.0, 0.0,  3.0),
        ])
        result = validate_orientation(
            gro, ec_residues={1}, ic_residues={5}, extracellular_side="+z"
        )
        md = result.to_markdown()
        assert "FAIL" in md


# ─────────────────────────────────────────────────────────────────────────────
# Passed property
# ─────────────────────────────────────────────────────────────────────────────

class TestPassedProperty:
    def test_passed_when_no_errors(self, tmp_path):
        gro = _write_gro(tmp_path, [
            (1, "ALA", "CA", 0.0, 0.0,  3.0),
            (5, "ALA", "CA", 0.0, 0.0, -3.0),
        ])
        result = validate_orientation(
            gro, ec_residues={1}, ic_residues={5}
        )
        assert result.passed

    def test_not_passed_when_errors(self, tmp_path):
        gro = _write_gro(tmp_path, [
            (1, "ALA", "CA", 0.0, 0.0, -3.0),  # inverted
            (5, "ALA", "CA", 0.0, 0.0,  3.0),
        ])
        result = validate_orientation(
            gro, ec_residues={1}, ic_residues={5}, extracellular_side="+z"
        )
        assert not result.passed
