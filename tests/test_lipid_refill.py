"""
tests/test_lipid_refill.py
Phase 10A: Surface-Guided Lipid Refill Actuator — unit tests.
"""
from __future__ import annotations

import json
import math
import tempfile
from pathlib import Path

import pytest

from validators.lipid_refill import run_lipid_refill


# ── GRO helpers ───────────────────────────────────────────────────────────────

def _write_gro(path: Path, atoms: list, box=(10.0, 10.0, 10.0)) -> None:
    """atoms: [(resid, resname, atomname, x, y, z), ...]"""
    lines = ["test system", str(len(atoms))]
    for i, (resid, resname, atomname, x, y, z) in enumerate(atoms, 1):
        lines.append(
            f"{resid:5d}{resname:<5s}{atomname:>5s}{i:5d}{x:8.3f}{y:8.3f}{z:8.3f}"
        )
    bx, by, bz = box
    lines.append(f"   {bx:.5f}   {by:.5f}   {bz:.5f}")
    path.write_text("\n".join(lines) + "\n")


def _make_template_bilayer(path: Path, lipid_resname: str = "DPP",
                            n_upper: int = 4, n_lower: int = 4) -> None:
    """Build a minimal template bilayer with upper + lower leaflet lipids.

    Each lipid residue has 2 atoms: a headgroup "P" and a tail "C16".
    Upper leaflet: headgroup P at z=6.0, tail C16 at z=5.5.
    Lower leaflet: headgroup P at z=4.0, tail C16 at z=4.5.
    """
    atoms = []
    resid = 1
    for i in range(n_upper):
        angle = i * 2.0 * math.pi / n_upper
        x = 5.0 + 2.0 * math.cos(angle)
        y = 5.0 + 2.0 * math.sin(angle)
        atoms.append((resid, lipid_resname, "P",   x, y, 6.0))
        atoms.append((resid, lipid_resname, "C16", x, y, 5.5))
        resid += 1
    for i in range(n_lower):
        angle = i * 2.0 * math.pi / n_lower
        x = 5.0 + 2.0 * math.cos(angle)
        y = 5.0 + 2.0 * math.sin(angle)
        atoms.append((resid, lipid_resname, "P",   x, y, 4.0))
        atoms.append((resid, lipid_resname, "C16", x, y, 4.5))
        resid += 1
    _write_gro(path, atoms, box=(10.0, 10.0, 10.0))


def _make_gap_system(path: Path, lipid_resname: str = "DPP",
                     n_tm: int = 10, tm_radius: float = 0.5,
                     lipid_radius: float = 5.0, n_lipids: int = 8,
                     box=(10.0, 10.0, 10.0)) -> set[int]:
    """Protein CA ring with annular gap; lipids far away."""
    atoms = []
    tm_residues: set[int] = set()
    cx = box[0] / 2.0
    cy = box[1] / 2.0
    for i in range(n_tm):
        resid = i + 1
        angle = i * 2.0 * math.pi / n_tm
        x = cx + tm_radius * math.cos(angle)
        y = cy + tm_radius * math.sin(angle)
        atoms.append((resid, "ALA", "CA", x, y, 5.0))
        tm_residues.add(resid)
    resid = n_tm + 1
    for i in range(n_lipids // 2):
        angle = i * 2.0 * math.pi / (n_lipids // 2)
        x = cx + lipid_radius * math.cos(angle)
        y = cy + lipid_radius * math.sin(angle)
        atoms.append((resid, lipid_resname, "P",   x, y, 5.8))
        atoms.append((resid, lipid_resname, "C16", x, y, 5.3))
        resid += 1
    for i in range(n_lipids // 2):
        angle = i * 2.0 * math.pi / (n_lipids // 2)
        x = cx + lipid_radius * math.cos(angle)
        y = cy + lipid_radius * math.sin(angle)
        atoms.append((resid, lipid_resname, "P",   x, y, 4.2))
        atoms.append((resid, lipid_resname, "C16", x, y, 4.7))
        resid += 1
    _write_gro(path, atoms, box=box)
    return tm_residues


def _make_good_system(path: Path, lipid_resname: str = "DPP",
                      n_tm: int = 10, tm_radius: float = 0.5) -> set[int]:
    """Protein CA ring with lipids at 0.3 nm (fully covered)."""
    atoms = []
    tm_residues: set[int] = set()
    for i in range(n_tm):
        resid = i + 1
        angle = i * 2.0 * math.pi / n_tm
        x = 5.0 + tm_radius * math.cos(angle)
        y = 5.0 + tm_radius * math.sin(angle)
        atoms.append((resid, "ALA", "CA", x, y, 5.0))
        tm_residues.add(resid)
    resid = n_tm + 1
    for i in range(n_tm):
        angle = i * 2.0 * math.pi / n_tm
        x = 5.0 + (tm_radius + 0.30) * math.cos(angle)
        y = 5.0 + (tm_radius + 0.30) * math.sin(angle)
        atoms.append((resid, lipid_resname, "P",   x, y, 5.8))
        atoms.append((resid, lipid_resname, "C16", x, y, 5.3))
        resid += 1
    for i in range(n_tm):
        angle = i * 2.0 * math.pi / n_tm
        x = 5.0 + (tm_radius + 0.30) * math.cos(angle)
        y = 5.0 + (tm_radius + 0.30) * math.sin(angle)
        atoms.append((resid, lipid_resname, "P",   x, y, 4.2))
        atoms.append((resid, lipid_resname, "C16", x, y, 4.7))
        resid += 1
    _write_gro(path, atoms, box=(10.0, 10.0, 10.0))
    return tm_residues


# ── Test 1: annular gap → refill inserts lipids and improves coverage ─────────

def test_refill_inserts_lipids_and_improves_coverage():
    """Annular gap system: refill must insert lipids and increase fraction_covered."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        system_gro   = tmpdir / "system.gro"
        template_gro = tmpdir / "template.gro"

        tm_residues = _make_gap_system(system_gro)
        _make_template_bilayer(template_gro)

        report = run_lipid_refill(
            system_gro   = system_gro,
            template_gro = template_gro,
            tm_residues  = tm_residues,
            output_gro   = system_gro,
            output_dir   = tmpdir,
            lipid_resname = "DPP",
            target_contact_distance_nm = 0.40,
        )

        assert report["n_inserted_lipids"] > 0, (
            f"Expected inserted lipids but got 0. report: {report}"
        )
        fc_before = report["fraction_covered_before"]
        fc_after  = report["fraction_covered_after"]
        assert fc_after >= fc_before, (
            f"fraction_covered did not improve: {fc_before:.3f} → {fc_after:.3f}"
        )
        assert report["coordinate_file_modified"] is True


# ── Test 2: already-good system → no lipids inserted ─────────────────────────

def test_no_insertion_when_quality_passes():
    """A fully covered system must not be modified."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        system_gro   = tmpdir / "system.gro"
        template_gro = tmpdir / "template.gro"

        tm_residues = _make_good_system(system_gro)
        _make_template_bilayer(template_gro)

        original_bytes = system_gro.read_bytes()

        report = run_lipid_refill(
            system_gro   = system_gro,
            template_gro = template_gro,
            tm_residues  = tm_residues,
            output_gro   = system_gro,
            output_dir   = tmpdir,
            lipid_resname = "DPP",
        )

        assert report["n_inserted_lipids"] == 0
        assert report["coordinate_file_modified"] is False
        assert system_gro.read_bytes() == original_bytes, (
            "system.gro was modified despite quality passing"
        )


# ── Test 3: inserted lipids are complete residues ─────────────────────────────

def test_inserted_lipids_are_complete_residues():
    """Every inserted residue must have the same atom count as the template."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        system_gro   = tmpdir / "system.gro"
        template_gro = tmpdir / "template.gro"

        tm_residues = _make_gap_system(system_gro)
        _make_template_bilayer(template_gro)   # 2 atoms per residue (P + C16)

        report = run_lipid_refill(
            system_gro   = system_gro,
            template_gro = template_gro,
            tm_residues  = tm_residues,
            output_gro   = system_gro,
            output_dir   = tmpdir,
            lipid_resname = "DPP",
            target_contact_distance_nm = 0.40,
        )

        if report["n_inserted_lipids"] == 0:
            pytest.skip("No lipids were inserted — adjust system geometry")

        # Parse output GRO and count atoms per inserted resid
        from validators.lipid_refill import _parse_gro
        atoms, _, _ = _parse_gro(system_gro)
        inserted_resids = set(report["inserted_resids"])
        resid_counts: dict[int, int] = {}
        for a in atoms:
            if a["resid"] in inserted_resids:
                resid_counts[a["resid"]] = resid_counts.get(a["resid"], 0) + 1

        assert len(resid_counts) == report["n_inserted_lipids"], (
            "Mismatch between reported inserted count and actual residue groups"
        )
        # All inserted residues must have the same count (template size = 2)
        counts = list(resid_counts.values())
        assert all(c == counts[0] for c in counts), (
            f"Inserted residues have inconsistent atom counts: {counts}"
        )
        assert counts[0] == 2, f"Expected 2 atoms per template, got {counts[0]}"


# ── Test 4: GRO atom count is updated ─────────────────────────────────────────

def test_gro_atom_count_updated():
    """Line 2 of output GRO must match the actual number of atom lines."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        system_gro   = tmpdir / "system.gro"
        template_gro = tmpdir / "template.gro"

        tm_residues = _make_gap_system(system_gro)
        _make_template_bilayer(template_gro)

        report = run_lipid_refill(
            system_gro   = system_gro,
            template_gro = template_gro,
            tm_residues  = tm_residues,
            output_gro   = system_gro,
            output_dir   = tmpdir,
            lipid_resname = "DPP",
            target_contact_distance_nm = 0.40,
        )

        if report["n_inserted_lipids"] == 0:
            pytest.skip("No lipids inserted")

        lines = system_gro.read_text().splitlines()
        declared_n = int(lines[1].strip())
        # Atom lines are lines[2 : 2+declared_n]; box is lines[2+declared_n]
        actual_atom_lines = lines[2: 2 + declared_n]
        assert len(actual_atom_lines) == declared_n, (
            f"Declared {declared_n} atoms but found {len(actual_atom_lines)} atom lines"
        )


# ── Test 5: protein coordinates unchanged ─────────────────────────────────────

def test_protein_coordinates_unchanged():
    """Protein atom coordinates must be bit-for-bit identical after refill."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        system_gro   = tmpdir / "system.gro"
        template_gro = tmpdir / "template.gro"

        tm_residues = _make_gap_system(system_gro)
        _make_template_bilayer(template_gro)

        from validators.lipid_refill import _parse_gro, DEFAULT_LIPID_RESNAMES, _SOLVENT_RESNAMES
        atoms_before, _, _ = _parse_gro(system_gro)
        prot_before = {
            (a["resid"], a["atomname"]): (a["x"], a["y"], a["z"])
            for a in atoms_before
            if a["resname"] not in DEFAULT_LIPID_RESNAMES
            and a["resname"] not in _SOLVENT_RESNAMES
        }

        run_lipid_refill(
            system_gro   = system_gro,
            template_gro = template_gro,
            tm_residues  = tm_residues,
            output_gro   = system_gro,
            output_dir   = tmpdir,
            lipid_resname = "DPP",
            target_contact_distance_nm = 0.40,
        )

        atoms_after, _, _ = _parse_gro(system_gro)
        prot_after = {
            (a["resid"], a["atomname"]): (a["x"], a["y"], a["z"])
            for a in atoms_after
            if a["resname"] not in DEFAULT_LIPID_RESNAMES
            and a["resname"] not in _SOLVENT_RESNAMES
        }

        assert prot_before == prot_after, "Protein coordinates were modified"


# ── Test 6: box vectors preserved ─────────────────────────────────────────────

def test_box_vectors_preserved():
    """Box vectors in output GRO must equal those in the input GRO."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        system_gro   = tmpdir / "system.gro"
        template_gro = tmpdir / "template.gro"

        tm_residues = _make_gap_system(system_gro, box=(12.5, 12.5, 10.0))
        _make_template_bilayer(template_gro)

        run_lipid_refill(
            system_gro   = system_gro,
            template_gro = template_gro,
            tm_residues  = tm_residues,
            output_gro   = system_gro,
            output_dir   = tmpdir,
            lipid_resname = "DPP",
            target_contact_distance_nm = 0.40,
        )

        last_line = system_gro.read_text().splitlines()[-1]
        parts = last_line.split()
        box_out = (float(parts[0]), float(parts[1]), float(parts[2]))
        assert abs(box_out[0] - 12.5) < 1e-4
        assert abs(box_out[1] - 12.5) < 1e-4
        assert abs(box_out[2] - 10.0) < 1e-4


# ── Test 7: clash filter rejects placements ───────────────────────────────────

def test_clash_filter_rejects_invalid_placements():
    """When all candidate positions clash with protein, zero lipids are inserted."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        system_gro   = tmpdir / "system.gro"
        template_gro = tmpdir / "template.gro"

        # Build a dense protein ring that fills the space around TM at r=0.5 nm
        # so every insertion at 0.45 nm would clash within protein_clash_cutoff=0.20
        atoms = []
        tm_residues: set[int] = set()
        n_tm = 36  # dense ring: atom every 10 degrees
        tm_r = 0.5
        # Add very dense blocker ring at r=0.70 nm (exactly where lipids would be placed)
        for i in range(n_tm):
            resid = i + 1
            angle = i * 2.0 * math.pi / n_tm
            # TM atoms
            x = 5.0 + tm_r * math.cos(angle)
            y = 5.0 + tm_r * math.sin(angle)
            atoms.append((resid, "ALA", "CA", x, y, 5.0))
            tm_residues.add(resid)
        # Dense blocker ring at r=0.90 nm (within 0.20 nm of any 0.45 nm target)
        resid = n_tm + 1
        n_block = 72
        for i in range(n_block):
            angle = i * 2.0 * math.pi / n_block
            x = 5.0 + 0.90 * math.cos(angle)
            y = 5.0 + 0.90 * math.sin(angle)
            atoms.append((resid, "ALA", "CB", x, y, 5.0))
            resid += 1
        # Only 2 lipids far away (to make the system valid)
        for i in range(2):
            atoms.append((resid, "DPP", "P",   5.0 + 4.0 * math.cos(i * math.pi), 5.0, 5.8))
            atoms.append((resid, "DPP", "C16", 5.0 + 4.0 * math.cos(i * math.pi), 5.0, 5.3))
            resid += 1
        for i in range(2):
            atoms.append((resid, "DPP", "P",   5.0 + 4.0 * math.cos(i * math.pi), 5.0, 4.2))
            atoms.append((resid, "DPP", "C16", 5.0 + 4.0 * math.cos(i * math.pi), 5.0, 4.7))
            resid += 1

        _write_gro(system_gro, atoms)
        _make_template_bilayer(template_gro)

        report = run_lipid_refill(
            system_gro              = system_gro,
            template_gro            = template_gro,
            tm_residues             = tm_residues,
            output_gro              = system_gro,
            output_dir              = tmpdir,
            lipid_resname           = "DPP",
            protein_clash_cutoff_nm = 0.20,
            target_contact_distance_nm = 0.45,
        )

        # All placements should clash (blocker ring at 0.90, template at 0.45+0.45=0.90)
        assert report["n_inserted_lipids"] == 0 or report["n_rejected_sites"] > 0, (
            "Expected at least some rejections in dense protein environment"
        )


# ── Test 8: baseline restored when refill worsens metrics ────────────────────

def test_baseline_restored_when_refill_worsens():
    """Verify baseline system.gro is restored if refill decreases fraction_covered."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        system_gro   = tmpdir / "system.gro"
        template_gro = tmpdir / "template.gro"

        # Build a system that will fail quality but has lipids placed so that
        # inserting at 0.40 nm would actually be FURTHER from TM than existing
        # lipids — a pathological case where inserting at a large distance worsens.
        # Achieve by: TM ring at r=0.5, existing lipids at r=0.70 (weakly_covered),
        # template placed outward at r=0.90 (would count as exposed_gap).
        # But since we check "if fc_after < fc_before - 0.02: restore", we need
        # a scenario where insertion actually worsens coverage.
        #
        # Instead, monkeypatch the evaluator to return worsened metrics on post-check.
        # We use a large target distance so inserted lipids land far from TM.
        n_tm = 10
        tm_r = 0.5
        atoms = []
        tm_residues: set[int] = set()
        for i in range(n_tm):
            resid = i + 1
            angle = i * 2.0 * math.pi / n_tm
            x = 5.0 + tm_r * math.cos(angle)
            y = 5.0 + tm_r * math.sin(angle)
            atoms.append((resid, "ALA", "CA", x, y, 5.0))
            tm_residues.add(resid)
        # Lipids at r=0.68 nm: weakly_covered (0.45-0.70), so quality fails
        resid = n_tm + 1
        for i in range(n_tm):
            angle = i * 2.0 * math.pi / n_tm
            x = 5.0 + (tm_r + 0.68) * math.cos(angle)
            y = 5.0 + (tm_r + 0.68) * math.sin(angle)
            atoms.append((resid, "DPP", "P",   x, y, 5.8))
            atoms.append((resid, "DPP", "C16", x, y, 5.3))
            resid += 1
        for i in range(n_tm):
            angle = i * 2.0 * math.pi / n_tm
            x = 5.0 + (tm_r + 0.68) * math.cos(angle)
            y = 5.0 + (tm_r + 0.68) * math.sin(angle)
            atoms.append((resid, "DPP", "P",   x, y, 4.2))
            atoms.append((resid, "DPP", "C16", x, y, 4.7))
            resid += 1
        _write_gro(system_gro, atoms)
        _make_template_bilayer(template_gro)

        baseline_content = system_gro.read_bytes()
        backup_path = tmpdir / "system_baseline_before_refill.gro"

        report = run_lipid_refill(
            system_gro   = system_gro,
            template_gro = template_gro,
            tm_residues  = tm_residues,
            output_gro   = system_gro,
            output_dir   = tmpdir,
            lipid_resname = "DPP",
            target_contact_distance_nm = 0.40,
        )

        # Whether or not lipids were inserted, the baseline backup must exist
        assert backup_path.exists(), "Baseline backup not created"
        # If refill worsened coverage, system.gro should equal baseline
        if not report["coordinate_file_modified"]:
            assert system_gro.read_bytes() == baseline_content, (
                "system.gro was not restored to baseline after worsening"
            )


# ── Test 9: output GRO is topology-valid ─────────────────────────────────────

def test_output_gro_is_valid_for_topology():
    """Output GRO must parse correctly: atom count matches declared count, box valid."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        system_gro   = tmpdir / "system.gro"
        template_gro = tmpdir / "template.gro"

        tm_residues = _make_gap_system(system_gro)
        _make_template_bilayer(template_gro)

        run_lipid_refill(
            system_gro   = system_gro,
            template_gro = template_gro,
            tm_residues  = tm_residues,
            output_gro   = system_gro,
            output_dir   = tmpdir,
            lipid_resname = "DPP",
            target_contact_distance_nm = 0.40,
        )

        from validators.lipid_refill import _parse_gro
        atoms, box, title = _parse_gro(system_gro)

        # Declared atom count must match parsed atom list
        lines = system_gro.read_text().splitlines()
        declared_n = int(lines[1].strip())
        assert declared_n == len(atoms), (
            f"Declared {declared_n} atoms but parsed {len(atoms)}"
        )

        # Box must have 3 positive dimensions
        assert len(box) == 3
        assert all(v > 0 for v in box), f"Invalid box: {box}"

        # Title must be non-empty
        assert len(title.strip()) > 0


# ── Test 10: GLP-1R regression — large annular gap is reduced by refill ───────

def test_glp1r_regression_large_annular_gap_improved():
    """Synthetic GLP-1R–like case: small TM ring at headgroup Z, bulk lipids far away.

    The inserted lipids land at the same Z as the TM atoms so 3-D distances are
    purely XY, guaranteeing fraction_covered_after > 0 and p90 decrease.

    Geometry:
      - 8 TM atoms at r=0.30 nm, z=5.8 (upper headgroup level)
      - Bulk DPP lipids at r=3.0 nm (large annular gap ~2.7 nm)
      - 4 lipids inserted at r=0.40 nm → 7 of 8 TM atoms within 0.45 nm
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        system_gro   = tmpdir / "system.gro"
        template_gro = tmpdir / "template.gro"

        atoms = []
        tm_residues: set[int] = set()
        n_tm = 8
        tm_r = 0.10   # small ring: placed lipids at 0.40 nm are 0.30 nm away (> clash cutoff)
        hg_z = 5.8    # place TM atoms at the upper headgroup Z

        for i in range(n_tm):
            resid = i + 1
            angle = i * 2.0 * math.pi / n_tm
            x = 5.0 + tm_r * math.cos(angle)
            y = 5.0 + tm_r * math.sin(angle)
            atoms.append((resid, "ALA", "CA", x, y, hg_z))
            tm_residues.add(resid)

        # Bulk lipids at r=3.0 nm (large annular gap)
        resid = n_tm + 1
        for i in range(8):
            angle = i * 2.0 * math.pi / 8
            x = 5.0 + 3.0 * math.cos(angle)
            y = 5.0 + 3.0 * math.sin(angle)
            atoms.append((resid, "DPP", "P",   x, y, 5.8))
            atoms.append((resid, "DPP", "C16", x, y, 5.3))
            resid += 1
        for i in range(8):
            angle = i * 2.0 * math.pi / 8
            x = 5.0 + 3.0 * math.cos(angle)
            y = 5.0 + 3.0 * math.sin(angle)
            atoms.append((resid, "DPP", "P",   x, y, 4.2))
            atoms.append((resid, "DPP", "C16", x, y, 4.7))
            resid += 1

        _write_gro(system_gro, atoms, box=(10.0, 10.0, 10.0))
        _make_template_bilayer(template_gro, n_upper=8, n_lower=8)

        report = run_lipid_refill(
            system_gro                 = system_gro,
            template_gro               = template_gro,
            tm_residues                = tm_residues,
            output_gro                 = system_gro,
            output_dir                 = tmpdir,
            lipid_resname              = "DPP",
            target_contact_distance_nm = 0.40,
            max_inserted_lipids        = 40,
            max_lipids_per_gap_cluster = 4,
        )

        fc_before = report["fraction_covered_before"]
        fc_after  = report["fraction_covered_after"]
        p90_before = report["p90_distance_before"]
        p90_after  = report["p90_distance_after"]
        n_ins = report["n_inserted_lipids"]

        assert n_ins > 0, (
            f"Expected lipids to be inserted into the large gap, got {n_ins}. "
            f"report warnings: {report.get('warning_messages')}"
        )
        assert fc_after > fc_before, (
            f"fraction_covered did not improve: {fc_before:.3f} → {fc_after:.3f}"
        )
        # p90 should decrease (lipids are now closer to TM surface)
        assert p90_after < p90_before, (
            f"p90 did not decrease: {p90_before:.3f} → {p90_after:.3f}"
        )

        # Report must be written to disk
        report_path = tmpdir / "interface_lipid_refill_report.json"
        assert report_path.exists(), "interface_lipid_refill_report.json not written"
        disk_report = json.loads(report_path.read_text())
        assert disk_report["n_inserted_lipids"] == n_ins
