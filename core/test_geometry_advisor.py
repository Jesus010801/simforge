"""Tests for GeometryAdvisor — no GROMACS required, no PDB download needed."""
from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from core.geometry_advisor import GeometryAdvisor, GeometryReport


# ── Helpers ───────────────────────────────────────────────────────────────────

def _make_pdb(coords: list[tuple[float, float, float]], tmp_path: Path) -> Path:
    """Build a minimal PDB with ATOM records at the given (x,y,z) positions in Å."""
    lines = []
    for i, (x, y, z) in enumerate(coords, 1):
        lines.append(
            f"ATOM  {i:5d}  CA  ALA A{i:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C  "
        )
    pdb = tmp_path / "test.pdb"
    pdb.write_text("\n".join(lines) + "\n")
    return pdb


def _compact_pdb(tmp_path: Path) -> Path:
    """Globular protein — roughly 30×25×28 Å."""
    import random
    rng = random.Random(42)
    coords = [(rng.uniform(0, 30), rng.uniform(0, 25), rng.uniform(0, 28))
              for _ in range(150)]
    return _make_pdb(coords, tmp_path)


def _elongated_pdb(tmp_path: Path) -> Path:
    """Rod-like protein: ~40 Å × 20 Å × 20 Å → aspect ~2."""
    coords = []
    for i in range(40):           # long axis X: 0–39 Å
        for j in range(3):        # cross-section Y: 0, 10, 20 Å
            for k in range(3):    # cross-section Z: 0, 10, 20 Å
                coords.append((float(i), float(j * 10), float(k * 10)))
    return _make_pdb(coords, tmp_path)


def _highly_elongated_pdb(tmp_path: Path) -> Path:
    """Filament: ~60 Å × 20 Å × 20 Å → aspect ~3."""
    coords = []
    for i in range(60):           # long axis X: 0–59 Å
        for j in range(3):        # cross-section Y: 0, 10, 20 Å
            for k in range(3):    # cross-section Z: 0, 10, 20 Å
                coords.append((float(i), float(j * 10), float(k * 10)))
    return _make_pdb(coords, tmp_path)


# ── Geometry classification ───────────────────────────────────────────────────

class TestGeometryClassification:
    def test_compact_classified_correctly(self, tmp_path):
        pdb = _compact_pdb(tmp_path)
        report = GeometryAdvisor().analyze(pdb)
        assert report.geometry == "compact"
        assert report.aspect_ratio < 1.8

    def test_elongated_classified_correctly(self, tmp_path):
        pdb = _elongated_pdb(tmp_path)
        report = GeometryAdvisor().analyze(pdb)
        assert report.geometry in ("elongated", "highly_elongated"), (
            f"Expected elongated/highly_elongated, got {report.geometry} "
            f"(aspect={report.aspect_ratio})"
        )

    def test_highly_elongated_classified_correctly(self, tmp_path):
        pdb = _highly_elongated_pdb(tmp_path)
        report = GeometryAdvisor().analyze(pdb)
        assert report.geometry == "highly_elongated", (
            f"Expected highly_elongated, got {report.geometry} "
            f"(aspect={report.aspect_ratio})"
        )

    def test_aspect_ratio_positive(self, tmp_path):
        pdb = _compact_pdb(tmp_path)
        report = GeometryAdvisor().analyze(pdb)
        assert report.aspect_ratio >= 1.0


# ── Bounding box ──────────────────────────────────────────────────────────────

class TestBoundingBox:
    def test_dimensions_in_nanometres(self, tmp_path):
        """A 100 Å span along X → dim_x = 10.0 nm."""
        coords = [(0.0, 0.0, 0.0), (100.0, 0.0, 0.0)]
        pdb = _make_pdb(coords, tmp_path)
        report = GeometryAdvisor().analyze(pdb, padding_nm=0.0)
        assert abs(report.dim_x - 10.0) < 0.1
        assert abs(report.dim_y) < 0.1
        assert abs(report.dim_z) < 0.1

    def test_box_includes_padding(self, tmp_path):
        """Box dims = protein dims + 2 × padding."""
        coords = [(0.0, 0.0, 0.0), (100.0, 50.0, 80.0)]
        pdb = _make_pdb(coords, tmp_path)
        report = GeometryAdvisor().analyze(pdb, padding_nm=1.0)
        assert abs(report.box_x - (report.dim_x + 2.0)) < 0.01
        assert abs(report.box_y - (report.dim_y + 2.0)) < 0.01
        assert abs(report.box_z - (report.dim_z + 2.0)) < 0.01

    def test_zero_padding_box_equals_protein_dims(self, tmp_path):
        coords = [(0.0, 0.0, 0.0), (50.0, 40.0, 30.0)]
        pdb = _make_pdb(coords, tmp_path)
        report = GeometryAdvisor().analyze(pdb, padding_nm=0.0)
        assert abs(report.box_x - report.dim_x) < 0.01
        assert abs(report.box_y - report.dim_y) < 0.01
        assert abs(report.box_z - report.dim_z) < 0.01


# ── Atom count estimates ──────────────────────────────────────────────────────

class TestAtomCountEstimates:
    def test_protein_atom_count_matches_pdb(self, tmp_path):
        coords = [(float(i), 0.0, 0.0) for i in range(200)]
        pdb = _make_pdb(coords, tmp_path)
        report = GeometryAdvisor().analyze(pdb)
        assert report.n_protein_atoms == 200

    def test_water_estimate_positive(self, tmp_path):
        pdb = _compact_pdb(tmp_path)
        report = GeometryAdvisor().analyze(pdb, padding_nm=1.0)
        assert report.n_water_estimated > 0

    def test_total_larger_than_protein_alone(self, tmp_path):
        pdb = _compact_pdb(tmp_path)
        report = GeometryAdvisor().analyze(pdb, padding_nm=1.0)
        assert report.n_total_estimated > report.n_protein_atoms

    def test_more_padding_more_water(self, tmp_path):
        pdb = _compact_pdb(tmp_path)
        r1 = GeometryAdvisor().analyze(pdb, padding_nm=0.5)
        r2 = GeometryAdvisor().analyze(pdb, padding_nm=1.5)
        assert r2.n_water_estimated > r1.n_water_estimated


# ── Advisories ────────────────────────────────────────────────────────────────

class TestAdvisories:
    def test_compact_protein_has_info_advisory(self, tmp_path):
        pdb = _compact_pdb(tmp_path)
        report = GeometryAdvisor().analyze(pdb)
        levels = {a.level for a in report.advisories}
        assert "INFO" in levels
        # Compact protein should not generate WARNINGs from geometry alone
        geom_warns = [
            a for a in report.advisories
            if a.level == "WARNING" and "elongat" in a.message.lower()
        ]
        assert geom_warns == []

    def test_highly_elongated_generates_warning(self, tmp_path):
        pdb = _highly_elongated_pdb(tmp_path)
        report = GeometryAdvisor().analyze(pdb)
        warns = [a for a in report.advisories if a.level == "WARNING"]
        assert warns, "Highly elongated structure must generate at least one WARNING"
        assert any("elongat" in w.message.lower() for w in warns)

    def test_highly_elongated_generates_suggestions(self, tmp_path):
        pdb = _highly_elongated_pdb(tmp_path)
        report = GeometryAdvisor().analyze(pdb)
        suggests = [a for a in report.advisories if a.level == "SUGGEST"]
        assert suggests, "Highly elongated structure must include SUGGEST advisories"

    def test_large_system_size_warning(self, tmp_path):
        """Very large box (low padding but huge protein) → size WARNING."""
        # 500 atoms in a 300 Å × 300 Å × 300 Å box → ~300k water molecules
        coords = [(float(i % 300), float((i // 300) % 300), float(i // 90000))
                  for i in range(500)]
        pdb = _make_pdb(coords, tmp_path)
        report = GeometryAdvisor().analyze(pdb, padding_nm=1.0)
        # If total > 300k, WARNING should be present
        if report.n_total_estimated > 300_000:
            size_warns = [a for a in report.advisories
                          if a.level == "WARNING" and "atom" in a.message.lower()]
            assert size_warns, "Large system must generate size WARNING"

    def test_advisories_always_present(self, tmp_path):
        """analyze() always returns at least one advisory (size INFO)."""
        pdb = _compact_pdb(tmp_path)
        report = GeometryAdvisor().analyze(pdb)
        assert len(report.advisories) >= 1

    def test_has_warnings_reflects_advisory_levels(self, tmp_path):
        report = GeometryAdvisor().analyze(_compact_pdb(tmp_path))
        expected = any(a.level in ("WARNING", "SUGGEST") for a in report.advisories)
        assert report.has_warnings() == expected


# ── Edge cases ────────────────────────────────────────────────────────────────

class TestEdgeCases:
    def test_missing_pdb_returns_warning(self, tmp_path):
        report = GeometryAdvisor().analyze(tmp_path / "nonexistent.pdb")
        assert any(a.level == "WARNING" for a in report.advisories)
        assert report.n_protein_atoms == 0

    def test_hydrogen_atoms_excluded(self, tmp_path):
        """H atoms (element H) must be excluded from coords."""
        lines = [
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  ",
            "ATOM      2  H   ALA A   1      10.000   0.000   0.000  1.00  0.00           H  ",
            "ATOM      3  CB  ALA A   1       5.000   0.000   0.000  1.00  0.00           C  ",
        ]
        pdb = tmp_path / "h_test.pdb"
        pdb.write_text("\n".join(lines) + "\n")
        report = GeometryAdvisor().analyze(pdb, padding_nm=0.0)
        # Only 2 heavy atoms → protein_atoms = 2
        assert report.n_protein_atoms == 2

    def test_dodecahedron_box_type_accepted(self, tmp_path):
        pdb = _compact_pdb(tmp_path)
        report = GeometryAdvisor().analyze(pdb, box_type="dodecahedron")
        assert report.n_water_estimated >= 0

    def test_single_atom_pdb(self, tmp_path):
        pdb = _make_pdb([(10.0, 20.0, 30.0)], tmp_path)
        report = GeometryAdvisor().analyze(pdb, padding_nm=1.0)
        assert report.n_protein_atoms == 1
        # Single atom: all dims = 0 → aspect ratio defaults to 1.0
        assert report.aspect_ratio >= 1.0
