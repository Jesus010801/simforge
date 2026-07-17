"""
tests/test_protein_membrane_interface.py
Phase 9A: Protein–Membrane Interface Evaluator — unit tests.
"""
from __future__ import annotations

import json
import math
import tempfile
from pathlib import Path

import pytest

from validators.protein_membrane_interface import evaluate_protein_membrane_interface


# ── GRO builder ───────────────────────────────────────────────────────────────

def _write_gro(path, atoms, box=(10.0, 10.0, 10.0)):
    """atoms: list of (resid, resname, atomname, x, y, z)"""
    lines = ["test system", str(len(atoms))]
    for i, (resid, resname, atomname, x, y, z) in enumerate(atoms, 1):
        lines.append(f"{resid:5d}{resname:<5s}{atomname:>5s}{i:5d}{x:8.3f}{y:8.3f}{z:8.3f}")
    bx, by, bz = box
    lines.append(f"   {bx:.5f}   {by:.5f}   {bz:.5f}")
    Path(path).write_text("\n".join(lines) + "\n")


# ── Test 1: fully covered TM bundle → quality_passed True ─────────────────────

def test_covered_tm_bundle_passes_quality():
    """10-atom TM chain with lipids at 0.3 nm from each atom → quality_passed=True."""
    with tempfile.TemporaryDirectory() as tmpdir:
        gro = Path(tmpdir) / "system.gro"
        atoms = []
        tm_residues = set(range(1, 11))

        # 10 TM CA atoms in a vertical line
        for ri in range(10):
            z = 4.0 + ri * 0.2
            atoms.append((ri + 1, "ALA", "CA", 5.0, 5.0, z))

        # 4 lipid P atoms per TM atom, all at 0.3 nm → covered
        lipid_resid = 100
        for ri in range(10):
            z = 4.0 + ri * 0.2
            for dx, dy in [(0.3, 0.0), (-0.3, 0.0), (0.0, 0.3), (0.0, -0.3)]:
                atoms.append((lipid_resid, "DPPC", "P", 5.0 + dx, 5.0 + dy, z))
                lipid_resid += 1

        _write_gro(gro, atoms)
        report = evaluate_protein_membrane_interface(
            gro_path=gro,
            tm_residues=tm_residues,
            output_dir=Path(tmpdir),
        )

        assert report["quality_passed"] is True, report["warnings"]
        assert report["fraction_covered"] >= 0.75
        assert report["fraction_exposed_gap"] <= 0.15
        assert report["n_tm_surface_atoms"] == 10


# ── Test 2: TM bundle with ring-shaped lipid gap → quality_passed False ──────

def test_tm_bundle_with_lipid_gap_fails_quality():
    """Half the TM ring has no nearby lipids → fraction_exposed_gap > 0.15 → quality_passed=False."""
    with tempfile.TemporaryDirectory() as tmpdir:
        gro = Path(tmpdir) / "system.gro"
        atoms = []
        tm_residues = set(range(1, 21))

        # 20 TM CA atoms in a ring of radius 0.5 nm at z=5.0
        for ri in range(20):
            angle = ri * 2.0 * math.pi / 20.0
            x = 5.0 + 0.5 * math.cos(angle)
            y = 5.0 + 0.5 * math.sin(angle)
            atoms.append((ri + 1, "ALA", "CA", x, y, 5.0))

        # Lipids covering only the first half of the ring (ri 0-9)
        lipid_resid = 100
        for ri in range(10):
            angle = ri * 2.0 * math.pi / 20.0
            lx = 5.0 + 0.8 * math.cos(angle)
            ly = 5.0 + 0.8 * math.sin(angle)
            atoms.append((lipid_resid, "DPPC", "P", lx, ly, 5.0))
            lipid_resid += 1

        _write_gro(gro, atoms)
        report = evaluate_protein_membrane_interface(
            gro_path=gro,
            tm_residues=tm_residues,
            output_dir=Path(tmpdir),
        )

        assert report["quality_passed"] is False
        # Boundary lipids cover adjacent ring atoms, but the gap centre atoms remain exposed.
        # Ring radius 0.5 nm, lipid radius 0.8 nm → atoms ~90-120° from nearest lipid
        # are > 0.7 nm away.  Expect at least 3 genuinely exposed atoms.
        assert report["fraction_exposed_gap"] > 0.15
        assert report["n_exposed_gap"] >= 3
        # At least one gap cluster must be found
        assert len(report["gap_clusters"]) >= 1
        assert report["largest_contiguous_gap_cluster_size"] >= 1


# ── Test 3: lipids near ECD but not TM → no TM coverage credit ───────────────

def test_lipids_near_ecd_do_not_cover_tm():
    """Lipid atoms at z=8–9 (near ECD) give no coverage to TM atoms at z=4–5."""
    with tempfile.TemporaryDirectory() as tmpdir:
        gro = Path(tmpdir) / "system.gro"
        atoms = []
        tm_residues = set(range(1, 11))

        # 10 TM CA atoms at z=4.0–4.9
        for ri in range(10):
            atoms.append((ri + 1, "ALA", "CA", 5.0, 5.0, 4.0 + ri * 0.1))

        # 10 ECD atoms at z=8.0–8.9 (not in tm_residues)
        for ri in range(10):
            atoms.append((50 + ri, "ALA", "CA", 5.0, 5.0, 8.0 + ri * 0.1))

        # Lipid P atoms near ECD (z=8.0–8.9), far from TM (>3 nm)
        lipid_resid = 200
        for ri in range(10):
            atoms.append((lipid_resid, "DPPC", "P", 5.3, 5.0, 8.0 + ri * 0.1))
            lipid_resid += 1

        _write_gro(gro, atoms)
        report = evaluate_protein_membrane_interface(
            gro_path=gro,
            tm_residues=tm_residues,
            output_dir=Path(tmpdir),
        )

        # All TM atoms are exposed (nearest lipid > 3 nm away)
        assert report["quality_passed"] is False
        assert report["fraction_exposed_gap"] > 0.15
        assert report["mean_nearest_lipid_distance"] > 2.0


# ── Test 4: per-segment worst TM segment identified correctly ─────────────────

def test_worst_tm_segment_identified():
    """Two TM segments: seg1 covered, seg2 uncovered → worst_tm_segment = seg2."""
    with tempfile.TemporaryDirectory() as tmpdir:
        gro = Path(tmpdir) / "system.gro"
        atoms = []
        # Seg1: res 1-10 at z=4.0–4.9
        # Seg2: res 20-29 at z=7.0–7.9 (gap of 10 res → separate segment)
        tm_residues = set(range(1, 11)) | set(range(20, 30))

        for ri in range(10):
            atoms.append((ri + 1, "ALA", "CA", 5.0, 5.0, 4.0 + ri * 0.1))
        for ri in range(10):
            atoms.append((20 + ri, "ALA", "CA", 5.0, 5.0, 7.0 + ri * 0.1))

        # Lipids covering seg1 (z=4.0–4.9) at 0.3 nm, nothing near seg2
        lipid_resid = 100
        for ri in range(10):
            atoms.append((lipid_resid, "DPPC", "P", 5.3, 5.0, 4.0 + ri * 0.1))
            lipid_resid += 1

        _write_gro(gro, atoms)
        report = evaluate_protein_membrane_interface(
            gro_path=gro,
            tm_residues=tm_residues,
            output_dir=Path(tmpdir),
        )

        seg_metrics = report["coverage_fraction_by_tm_segment"]
        assert len(seg_metrics) == 2, f"Expected 2 TM segments, got {len(seg_metrics)}"

        worst = report["worst_tm_segment"]
        assert worst is not None
        # Worst segment must span seg2 residue range
        seg2_start, seg2_end = worst["residue_range"]
        assert seg2_start >= 20, f"Expected worst segment in seg2 (res 20-29), got {worst}"

        # Seg1 must have better coverage than seg2
        seg1_cov = next(
            s["coverage_fraction"] for s in seg_metrics if s["residue_range"][0] == 1
        )
        seg2_cov = next(
            s["coverage_fraction"] for s in seg_metrics if s["residue_range"][0] == 20
        )
        assert seg1_cov > seg2_cov, f"seg1={seg1_cov:.3f} should exceed seg2={seg2_cov:.3f}"


# ── Test 5: report file contains all required fields ─────────────────────────

def test_report_contains_all_required_fields():
    """All 23 expected fields must be present in the report and on disk."""
    required_fields = [
        "source_file",
        "tm_residue_ranges",
        "n_tm_surface_atoms",
        "fraction_covered",
        "fraction_weakly_covered",
        "fraction_exposed_gap",
        "mean_nearest_lipid_distance",
        "p90_nearest_lipid_distance",
        "max_nearest_lipid_distance",
        "mean_nearest_headgroup_distance",
        "mean_nearest_tail_distance",
        "n_covered",
        "n_weakly_covered",
        "n_exposed_gap",
        "largest_contiguous_gap_cluster_size",
        "gap_clusters",
        "coverage_fraction_by_residue",
        "mean_distance_by_residue",
        "coverage_fraction_by_tm_segment",
        "worst_residues",
        "worst_tm_segment",
        "worst_tm_segments",
        "quality_passed",
        "quality_thresholds",
        "warnings",
    ]
    with tempfile.TemporaryDirectory() as tmpdir:
        gro = Path(tmpdir) / "system.gro"
        atoms = []
        tm_residues = set(range(1, 6))

        for ri in range(5):
            atoms.append((ri + 1, "ALA", "CA", 5.0, 5.0, 4.0 + ri * 0.2))
        lipid_resid = 50
        for ri in range(5):
            atoms.append((lipid_resid + ri, "DPPC", "P", 5.3, 5.0, 4.0 + ri * 0.2))
            # tail C atoms
            atoms.append((lipid_resid + 20 + ri, "DPPC", "C16", 5.4, 5.0, 4.0 + ri * 0.2))

        _write_gro(gro, atoms)
        report = evaluate_protein_membrane_interface(
            gro_path=gro,
            tm_residues=tm_residues,
            output_dir=Path(tmpdir),
        )

        missing = [f for f in required_fields if f not in report]
        assert missing == [], f"Missing fields: {missing}"

        # JSON on disk
        report_path = Path(tmpdir) / "protein_membrane_interface_report.json"
        assert report_path.exists(), "JSON report not written to disk"
        on_disk = json.loads(report_path.read_text())
        for f in required_fields:
            assert f in on_disk, f"Field '{f}' missing from JSON on disk"

        # quality_thresholds sub-keys
        qt = report["quality_thresholds"]
        assert "min_fraction_covered" in qt
        assert "max_fraction_exposed_gap" in qt
        assert "max_p90_nearest_lipid_distance_nm" in qt


# ── Test 6: coordinate file not modified ─────────────────────────────────────

def test_evaluator_does_not_modify_coordinates():
    """evaluate_protein_membrane_interface must not alter the GRO file."""
    with tempfile.TemporaryDirectory() as tmpdir:
        gro = Path(tmpdir) / "system.gro"
        atoms = [
            (1, "ALA", "CA", 5.0, 5.0, 4.0),
            (2, "ALA", "CA", 5.0, 5.0, 4.2),
            (100, "DPPC", "P", 5.3, 5.0, 4.0),
            (101, "DPPC", "P", 5.3, 5.0, 4.2),
        ]
        _write_gro(gro, atoms)

        content_before = gro.read_bytes()

        evaluate_protein_membrane_interface(
            gro_path=gro,
            tm_residues={1, 2},
            output_dir=Path(tmpdir),
        )

        content_after = gro.read_bytes()
        assert content_before == content_after, "GRO file was modified by the evaluator"
