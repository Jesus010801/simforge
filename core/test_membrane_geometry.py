# core/test_membrane_geometry.py
"""Tests for membrane orientation geometry — no GROMACS dependency."""

from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from core.membrane_geometry import (
    parse_residue_range,
    compute_orient_rotation,
    validate_oriented_gro,
    GeometricValidationResult,
    OrientRotation,
)


# ── parse_residue_range ───────────────────────────────────────────────────────

def test_parse_simple_range():
    assert parse_residue_range("1-5") == [1, 2, 3, 4, 5]


def test_parse_single_residues():
    assert parse_residue_range("5,10,15") == [5, 10, 15]


def test_parse_mixed():
    result = parse_residue_range("1-3,7,10-12")
    assert result == [1, 2, 3, 7, 10, 11, 12]


def test_parse_dedup_and_sorted():
    result = parse_residue_range("5,5,3-5")
    assert result == [3, 4, 5]


def test_parse_empty():
    assert parse_residue_range("") == []
    assert parse_residue_range("  ") == []


def test_parse_single_value():
    assert parse_residue_range("42") == [42]


# ── Helper: build minimal GRO content ────────────────────────────────────────

def _make_gro(atoms: list[tuple[int, str, str, float, float, float]]) -> str:
    """
    Build a minimal GRO string.
    atoms: [(resnum, resname, atomname, x, y, z), ...]
    """
    lines = ["Test protein", str(len(atoms))]
    for i, (resnum, resname, atomname, x, y, z) in enumerate(atoms, 1):
        lines.append(
            f"{resnum:5d}{resname:<5s}{atomname:>5s}{i:5d}"
            f"{x:8.3f}{y:8.3f}{z:8.3f}"
        )
    lines.append("10.0 10.0 10.0")
    return "\n".join(lines) + "\n"


# ── compute_orient_rotation ───────────────────────────────────────────────────

def test_tm_along_z_ec_at_plus_z(tmp_path):
    """EC residues at +Z, IC at -Z → no rotation needed."""
    gro = tmp_path / "prot.gro"
    gro.write_text(_make_gro([
        (1, "ALA", "CA", 5.0, 5.0, 8.0),  # EC residue, high Z
        (2, "ALA", "CA", 5.0, 5.0, 2.0),  # IC residue, low Z
    ]))
    result = compute_orient_rotation(gro, ec_resids=[1], ic_resids=[2])
    assert result.rx == 0.0
    assert result.ry == 0.0
    assert result.rz == 0.0
    assert "no rotation" in result.description


def test_tm_along_z_ec_at_minus_z(tmp_path):
    """EC residues at -Z, IC at +Z → flip 180 around X."""
    gro = tmp_path / "prot.gro"
    gro.write_text(_make_gro([
        (1, "ALA", "CA", 5.0, 5.0, 2.0),  # EC at -Z (low Z)
        (2, "ALA", "CA", 5.0, 5.0, 8.0),  # IC at +Z (high Z)
    ]))
    result = compute_orient_rotation(gro, ec_resids=[1], ic_resids=[2])
    assert result.rx == 180.0
    assert result.ry == 0.0
    assert result.rz == 0.0
    assert "180" in result.description


def test_tm_along_x_ec_at_plus_x(tmp_path):
    """TM axis along X, EC at +X → rotate 0 270 0."""
    gro = tmp_path / "prot.gro"
    gro.write_text(_make_gro([
        (1, "ALA", "CA", 8.0, 5.0, 5.0),  # EC at +X
        (2, "ALA", "CA", 2.0, 5.0, 5.0),  # IC at -X
    ]))
    result = compute_orient_rotation(gro, ec_resids=[1], ic_resids=[2])
    assert result.rx == 0.0
    assert result.ry == 270.0
    assert result.rz == 0.0


def test_tm_along_x_ec_at_minus_x(tmp_path):
    """TM axis along X, EC at -X → rotate 0 90 0."""
    gro = tmp_path / "prot.gro"
    gro.write_text(_make_gro([
        (1, "ALA", "CA", 2.0, 5.0, 5.0),  # EC at -X
        (2, "ALA", "CA", 8.0, 5.0, 5.0),  # IC at +X
    ]))
    result = compute_orient_rotation(gro, ec_resids=[1], ic_resids=[2])
    assert result.rx == 0.0
    assert result.ry == 90.0
    assert result.rz == 0.0


def test_tm_along_y_ec_at_plus_y(tmp_path):
    """TM axis along Y, EC at +Y → rotate -90 0 0."""
    gro = tmp_path / "prot.gro"
    gro.write_text(_make_gro([
        (1, "ALA", "CA", 5.0, 8.0, 5.0),  # EC at +Y
        (2, "ALA", "CA", 5.0, 2.0, 5.0),  # IC at -Y
    ]))
    result = compute_orient_rotation(gro, ec_resids=[1], ic_resids=[2])
    assert result.rx == -90.0
    assert result.ry == 0.0
    assert result.rz == 0.0


def test_tm_along_y_ec_at_minus_y(tmp_path):
    """TM axis along Y, EC at -Y → rotate 90 0 0."""
    gro = tmp_path / "prot.gro"
    gro.write_text(_make_gro([
        (1, "ALA", "CA", 5.0, 2.0, 5.0),  # EC at -Y
        (2, "ALA", "CA", 5.0, 8.0, 5.0),  # IC at +Y
    ]))
    result = compute_orient_rotation(gro, ec_resids=[1], ic_resids=[2])
    assert result.rx == 90.0
    assert result.ry == 0.0
    assert result.rz == 0.0


def test_multi_residue_com(tmp_path):
    """COM uses mean of all Cα atoms, not just one."""
    gro = tmp_path / "prot.gro"
    gro.write_text(_make_gro([
        # EC group: two residues, mean Z = 7.5
        (1, "ALA", "CA", 5.0, 5.0, 7.0),
        (2, "ALA", "CA", 5.0, 5.0, 8.0),
        # IC group: two residues, mean Z = 2.5
        (3, "ALA", "CA", 5.0, 5.0, 2.0),
        (4, "ALA", "CA", 5.0, 5.0, 3.0),
    ]))
    result = compute_orient_rotation(gro, ec_resids=[1, 2], ic_resids=[3, 4])
    # EC mean Z > IC mean Z → TM along +Z → no rotation
    assert result.rx == 0.0
    assert result.ry == 0.0
    assert result.rz == 0.0


def test_non_ca_atoms_ignored(tmp_path):
    """Only CA atoms are used; CB, N, O are ignored."""
    gro = tmp_path / "prot.gro"
    gro.write_text(_make_gro([
        (1, "ALA", "CB", 5.0, 5.0, 9.0),  # CB at very high Z — should be ignored
        (1, "ALA", "CA", 5.0, 5.0, 7.0),  # CA at Z=7.0 — used
        (2, "ALA", "CB", 5.0, 5.0, 0.0),  # CB at very low Z — ignored
        (2, "ALA", "CA", 5.0, 5.0, 3.0),  # CA at Z=3.0 — used
    ]))
    result = compute_orient_rotation(gro, ec_resids=[1], ic_resids=[2])
    # EC CA at 7.0, IC CA at 3.0 → TM along +Z
    assert result.rx == 0.0
    assert result.ry == 0.0
    assert result.rz == 0.0


def test_missing_ec_residues_raises(tmp_path):
    """ValueError when EC residues have no Cα in GRO."""
    gro = tmp_path / "prot.gro"
    gro.write_text(_make_gro([
        (1, "ALA", "CA", 5.0, 5.0, 7.0),
    ]))
    with pytest.raises(ValueError, match="EC residues"):
        compute_orient_rotation(gro, ec_resids=[99], ic_resids=[1])


def test_missing_ic_residues_raises(tmp_path):
    """ValueError when IC residues have no Cα in GRO."""
    gro = tmp_path / "prot.gro"
    gro.write_text(_make_gro([
        (1, "ALA", "CA", 5.0, 5.0, 7.0),
    ]))
    with pytest.raises(ValueError, match="IC residues"):
        compute_orient_rotation(gro, ec_resids=[1], ic_resids=[99])


def test_orient_rotation_is_namedtuple():
    r = OrientRotation(0.0, 270.0, 0.0, "test")
    assert r.rx == 0.0
    assert r.ry == 270.0
    assert r.rz == 0.0
    assert r.description == "test"


# ═══════════════════════════════════════════════════════════════════════════════
# validate_oriented_gro — Phase 4 geometric validation
# ═══════════════════════════════════════════════════════════════════════════════

class TestValidateOrientedGro:
    """Geometric validation of an oriented GRO against structural annotation."""

    def _make_tm_gro(self, tmp_path, ec_z: float, ic_z: float, tm_z: float) -> Path:
        """GRO with EC at ec_z, IC at ic_z, TM at tm_z."""
        gro = tmp_path / "oriented.gro"
        gro.write_text(_make_gro([
            (1,  "GLY", "CA", 5.0, 5.0, ec_z),   # EC residue 1
            (2,  "GLY", "CA", 5.0, 5.0, ec_z),   # EC residue 2
            (50, "ALA", "CA", 5.0, 5.0, tm_z),   # TM residue 50
            (99, "LEU", "CA", 5.0, 5.0, ic_z),   # IC residue 99
        ]))
        return gro

    def test_correct_orientation_plus_z(self, tmp_path):
        """EC at +Z, IC at -Z, target +Z → orientation_correct."""
        gro = self._make_tm_gro(tmp_path, ec_z=7.0, ic_z=3.0, tm_z=5.0)
        results = validate_oriented_gro(
            gro,
            ec_resids=[1, 2],
            ic_resids=[99],
            tm_resids=[50],
            extracellular_side="+z",
        )
        codes = [r.code for r in results]
        assert "orientation_correct" in codes
        assert "orientation_inverted" not in codes
        assert "tm_in_membrane_zone" in codes

    def test_inverted_orientation_detected(self, tmp_path):
        """EC at -Z, IC at +Z, but target is +Z → orientation_inverted error."""
        gro = self._make_tm_gro(tmp_path, ec_z=3.0, ic_z=7.0, tm_z=5.0)
        results = validate_oriented_gro(
            gro,
            ec_resids=[1, 2],
            ic_resids=[99],
            extracellular_side="+z",
        )
        codes = [r.code for r in results]
        assert "orientation_inverted" in codes
        levels = {r.code: r.level for r in results}
        assert levels["orientation_inverted"] == "error"

    def test_minus_z_target_correct(self, tmp_path):
        """EC at -Z, IC at +Z, target -Z → orientation_correct."""
        gro = self._make_tm_gro(tmp_path, ec_z=3.0, ic_z=7.0, tm_z=5.0)
        results = validate_oriented_gro(
            gro,
            ec_resids=[1, 2],
            ic_resids=[99],
            extracellular_side="-z",
        )
        codes = [r.code for r in results]
        assert "orientation_correct" in codes
        assert "orientation_inverted" not in codes

    def test_tm_outside_zone_warns(self, tmp_path):
        """TM COM outside EC-IC span → tm_outside_membrane_zone warning."""
        gro = self._make_tm_gro(tmp_path, ec_z=8.0, ic_z=6.0, tm_z=2.0)
        results = validate_oriented_gro(
            gro,
            ec_resids=[1, 2],
            ic_resids=[99],
            tm_resids=[50],
            extracellular_side="+z",
        )
        codes = [r.code for r in results]
        assert "tm_outside_membrane_zone" in codes

    def test_tm_inside_zone_ok(self, tmp_path):
        """TM COM inside EC-IC span → tm_in_membrane_zone ok."""
        gro = self._make_tm_gro(tmp_path, ec_z=7.0, ic_z=3.0, tm_z=5.0)
        results = validate_oriented_gro(
            gro,
            ec_resids=[1, 2],
            ic_resids=[99],
            tm_resids=[50],
            extracellular_side="+z",
        )
        codes = [r.code for r in results]
        assert "tm_in_membrane_zone" in codes

    def test_missing_ec_atoms_warns(self, tmp_path):
        """No Cα for EC residues → warning, not crash."""
        gro = tmp_path / "prot.gro"
        gro.write_text(_make_gro([(99, "LEU", "CA", 5.0, 5.0, 3.0)]))
        results = validate_oriented_gro(
            gro,
            ec_resids=[1, 2],  # not in GRO
            ic_resids=[99],
            extracellular_side="+z",
        )
        codes = [r.code for r in results]
        assert "ec_no_ca" in codes

    def test_no_tm_residues_skips_tm_check(self, tmp_path):
        """No tm_resids → no tm check, still validates EC/IC orientation."""
        gro = self._make_tm_gro(tmp_path, ec_z=7.0, ic_z=3.0, tm_z=5.0)
        results = validate_oriented_gro(
            gro,
            ec_resids=[1, 2],
            ic_resids=[99],
            tm_resids=None,
            extracellular_side="+z",
        )
        codes = [r.code for r in results]
        assert "tm_in_membrane_zone" not in codes
        assert "tm_outside_membrane_zone" not in codes
        assert "orientation_correct" in codes


# ═══════════════════════════════════════════════════════════════════════════════
# Topology order warnings — compile-time heuristic
# ═══════════════════════════════════════════════════════════════════════════════

class TestTopologyOrderWarnings:
    def test_single_pass_correct_order(self):
        from core.structural_annotation import StructuralAnnotation, MembraneTopologyAnnotation
        ann = StructuralAnnotation(
            membrane_topology=MembraneTopologyAnnotation(
                extracellular_regions=["1-30"],
                intracellular_regions=["70-100"],
                transmembrane_segments=["31-69"],
            ),
        )
        warns = ann._topology_order_warnings()
        assert warns == []

    def test_single_pass_tm_out_of_order_warns(self):
        from core.structural_annotation import StructuralAnnotation, MembraneTopologyAnnotation
        ann = StructuralAnnotation(
            membrane_topology=MembraneTopologyAnnotation(
                extracellular_regions=["1-30"],
                intracellular_regions=["70-100"],
                transmembrane_segments=["200-230"],  # outside EC-IC span
            ),
        )
        warns = ann._topology_order_warnings()
        assert len(warns) == 1
        assert "single-pass" in warns[0].lower() or "tm segment" in warns[0].lower()

    def test_multi_pass_skips_order_check(self):
        from core.structural_annotation import StructuralAnnotation, MembraneTopologyAnnotation
        ann = StructuralAnnotation(
            membrane_topology=MembraneTopologyAnnotation(
                extracellular_regions=["1-30"],
                intracellular_regions=["70-100"],
                transmembrane_segments=["31-50", "55-69"],  # 2 TM → no check
            ),
        )
        warns = ann._topology_order_warnings()
        assert warns == []
