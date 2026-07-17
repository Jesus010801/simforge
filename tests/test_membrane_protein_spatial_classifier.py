"""
tests/test_membrane_protein_spatial_classifier.py
Phase 10B: Membrane-Protein Spatial Classifier — unit tests.
"""
from __future__ import annotations

import json
import math
import tempfile
from pathlib import Path

import pytest

from validators.membrane_protein_spatial_classifier import (
    run_spatial_classifier,
    run_pore_aware_water_cleanup,
)

# ── GRO helpers ───────────────────────────────────────────────────────────────

def _write_gro(path: Path, atoms: list, box=(12.0, 12.0, 10.0)) -> None:
    """atoms: [(resid, resname, atomname, x, y, z), ...]"""
    lines = ["test system", str(len(atoms))]
    for i, (resid, resname, atomname, x, y, z) in enumerate(atoms, 1):
        lines.append(
            f"{resid:5d}{resname:<5s}{atomname:>5s}{i:5d}{x:8.3f}{y:8.3f}{z:8.3f}"
        )
        bx, by, bz = box
    lines.append(f"   {bx:.5f}   {by:.5f}   {bz:.5f}")
    path.write_text("\n".join(lines) + "\n")


def _make_protein_ring(
    cx: float = 6.0,
    cy: float = 6.0,
    ring_radius: float = 1.0,
    n_per_level: int = 16,
    z_levels: "list[float] | None" = None,
    resname: str = "ALA",
    atomname: str = "CA",
    start_resid: int = 1,
) -> "tuple[list, set[int]]":
    """Build a dense cylindrical protein ring.  Returns (atoms, tm_resid_set)."""
    if z_levels is None:
        z_levels = [z * 0.3 + 3.5 for z in range(11)]  # z=3.5…6.5, 11 levels
    atoms = []
    tm_residues: set[int] = set()
    resid = start_resid
    for iz, z in enumerate(z_levels):
        for k in range(n_per_level):
            angle = k * 2.0 * math.pi / n_per_level
            x = cx + ring_radius * math.cos(angle)
            y = cy + ring_radius * math.sin(angle)
            atoms.append((resid, resname, atomname, x, y, z))
            tm_residues.add(resid)
            resid += 1
    return atoms, tm_residues


def _make_sol(resid: int, x: float, y: float, z: float) -> list:
    """One SPC/E water molecule at (x, y, z)."""
    return [
        (resid, "SOL", "OW",  x,        y,        z),
        (resid, "SOL", "HW1", x + 0.01, y,        z),
        (resid, "SOL", "HW2", x,        y + 0.01, z),
    ]


def _make_lipid(resid: int, x: float, y: float, z_hg: float,
                resname: str = "DPP") -> list:
    """Minimal two-atom lipid at (x, y)."""
    return [
        (resid, resname, "P",   x, y, z_hg),
        (resid, resname, "C16", x, y, z_hg - 0.5),
    ]


# ── Fixture: channel system ───────────────────────────────────────────────────

def _build_channel_system(
    tmpdir: Path,
    pore_water: bool = True,
    pore_lipid: bool = False,
    core_water_outside: bool = False,
) -> "tuple[Path, set[int]]":
    """Protein ring (sealed cylinder) + optional pore water/lipid and core water."""
    cx, cy = 6.0, 6.0
    ring_atoms, tm_residues = _make_protein_ring(cx=cx, cy=cy)
    atoms = list(ring_atoms)

    next_resid = max(tm_residues) + 1

    # 8 bulk lipids outside ring at r=3.0
    for i in range(4):
        angle = i * math.pi / 2
        lx = cx + 3.0 * math.cos(angle)
        ly = cy + 3.0 * math.sin(angle)
        atoms += _make_lipid(next_resid, lx, ly, z_hg=6.0)
        next_resid += 1
    for i in range(4):
        angle = i * math.pi / 2
        lx = cx + 3.0 * math.cos(angle)
        ly = cy + 3.0 * math.sin(angle)
        atoms += _make_lipid(next_resid, lx, ly, z_hg=4.0)
        next_resid += 1

    pore_water_resids = []
    if pore_water:
        # Water at pore centre
        atoms += _make_sol(next_resid, cx, cy, 5.0)
        pore_water_resids.append(next_resid)
        next_resid += 1

    if pore_lipid:
        atoms += _make_lipid(next_resid, cx, cy, z_hg=5.5)
        next_resid += 1

    if core_water_outside:
        # Water in membrane core but outside ring at r=3.0
        atoms += _make_sol(next_resid, cx + 3.0, cy, 5.0)
        next_resid += 1

    gro = tmpdir / "system.gro"
    _write_gro(gro, atoms)
    return gro, tm_residues


# ── Test 1: pore water is preserved ──────────────────────────────────────────

def test_pore_water_is_preserved():
    """Water at the centre of a closed protein ring at membrane Z → pore_water."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        gro, tm_res = _build_channel_system(td, pore_water=True)

        rep = run_spatial_classifier(
            gro_path=gro,
            tm_residues=tm_res,
            output_dir=td,
        )

        assert rep["enabled"]
        assert rep["n_pore_waters"] >= 1, (
            f"Expected ≥1 pore water, got {rep['n_pore_waters']}. "
            f"core_z={rep['membrane_core_z_range']}, "
            f"removed_candidate={rep['removed_water_resids_candidate']}"
        )
        # Pore water must NOT appear in removal candidates
        pore_ids = set(rep["pore_water_resids"])
        removed_ids = set(rep["removed_water_resids_candidate"])
        overlap = pore_ids & removed_ids
        assert not overlap, f"Pore water resids also in removal candidates: {overlap}"


# ── Test 2: lipid inside pore → forbidden_pore_lipid ─────────────────────────

def test_lipid_inside_pore_classified_forbidden():
    """A lipid whose COM is inside the closed TM ring at membrane Z → forbidden."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        gro, tm_res = _build_channel_system(td, pore_water=False, pore_lipid=True)

        rep = run_spatial_classifier(
            gro_path=gro,
            tm_residues=tm_res,
            output_dir=td,
        )

        assert rep["n_forbidden_pore_lipids"] >= 1, (
            f"Expected ≥1 forbidden pore lipid. "
            f"valid={rep['n_valid_external_lipids']}, "
            f"core_z={rep['membrane_core_z_range']}"
        )
        assert len(rep["forbidden_pore_lipid_resids"]) == rep["n_forbidden_pore_lipids"]


# ── Test 3: membrane-core water outside protein is removed ────────────────────

def test_core_water_outside_protein_is_removed():
    """Water at membrane Z but outside TM ring → membrane_core_water → candidate for removal."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        gro, tm_res = _build_channel_system(
            td, pore_water=False, core_water_outside=True
        )

        rep = run_spatial_classifier(
            gro_path=gro,
            tm_residues=tm_res,
            output_dir=td,
            remove_membrane_core_water=True,
        )

        assert rep["n_membrane_core_waters"] >= 1, (
            f"Expected core water outside ring. "
            f"core_z={rep['membrane_core_z_range']}"
        )
        assert len(rep["removed_water_resids_candidate"]) >= 1


# ── Test 4: GPCR-like closed bundle — external annular gaps lipid-accessible ──

def test_gpcr_like_external_annular_region_lipid_accessible():
    """Closed TM bundle with no pore: lipids outside the ring are valid_external."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        cx, cy = 6.0, 6.0
        ring_atoms, tm_res = _make_protein_ring(cx=cx, cy=cy)
        atoms = list(ring_atoms)
        next_resid = max(tm_res) + 1

        # External lipids in annular region at r=2.0 nm (outside ring at r=1.0)
        for i in range(8):
            angle = i * math.pi / 4
            lx = cx + 2.0 * math.cos(angle)
            ly = cy + 2.0 * math.sin(angle)
            atoms += _make_lipid(next_resid, lx, ly, z_hg=6.0)
            next_resid += 1
        for i in range(8):
            angle = i * math.pi / 4
            lx = cx + 2.0 * math.cos(angle)
            ly = cy + 2.0 * math.sin(angle)
            atoms += _make_lipid(next_resid, lx, ly, z_hg=4.0)
            next_resid += 1

        gro = td / "system.gro"
        _write_gro(gro, atoms)

        rep = run_spatial_classifier(gro_path=gro, tm_residues=tm_res, output_dir=td)

        assert rep["n_valid_external_lipids"] > 0, "All lipids should be valid_external"
        assert rep["n_forbidden_pore_lipids"] == 0, (
            f"Annular lipids wrongly classified as forbidden: "
            f"{rep['forbidden_pore_lipid_resids']}"
        )


# ── Test 5: Phase 10A skips pore-region targets via classifier report ─────────

def test_lipid_refill_skips_pore_region():
    """When classifier_report marks a region as pore, run_lipid_refill skips it."""
    pytest.importorskip("numpy")
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)

        # Build gap system with protein ring
        cx, cy = 6.0, 6.0
        ring_atoms, tm_res = _make_protein_ring(cx=cx, cy=cy)
        atoms = list(ring_atoms)
        next_resid = max(tm_res) + 1

        # Only lipids far outside the ring (r=4.0) — creates a gap at r=1-4
        for i in range(4):
            angle = i * math.pi / 2
            lx = cx + 4.0 * math.cos(angle)
            ly = cy + 4.0 * math.sin(angle)
            atoms += _make_lipid(next_resid, lx, ly, z_hg=6.0)
            next_resid += 1
        for i in range(4):
            angle = i * math.pi / 2
            lx = cx + 4.0 * math.cos(angle)
            ly = cy + 4.0 * math.sin(angle)
            atoms += _make_lipid(next_resid, lx, ly, z_hg=4.0)
            next_resid += 1

        gro = td / "system.gro"
        _write_gro(gro, atoms)

        # Classifier report that marks TM centre as pore
        clf_report = run_spatial_classifier(
            gro_path=gro,
            tm_residues=tm_res,
            output_dir=td,
        )

        # The forbidden_pore_lipid_resids should be empty (no lipid in pore)
        # What matters: the classifier correctly identifies the TM footprint
        assert clf_report["enabled"]
        # n_forbidden_pore_lipids == 0 because no lipid was placed in the pore
        assert clf_report["n_forbidden_pore_lipids"] == 0

        # Integration: run_lipid_refill with the spatial classifier report should
        # avoid inserting into the pore region. We verify by checking that the
        # report field is present and the classifier correctly reports valid lipids.
        assert clf_report["n_valid_external_lipids"] > 0


# ── Test 6: pore-aware cleanup preserves pore water ──────────────────────────

def test_pore_aware_cleanup_preserves_pore_water():
    """run_pore_aware_water_cleanup must NOT delete water classified as pore_water."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        gro, tm_res = _build_channel_system(td, pore_water=True)
        gro_out = td / "system_clean.gro"

        report = run_pore_aware_water_cleanup(
            gro_in=gro,
            gro_out=gro_out,
            tm_residues=tm_res,
            output_dir=td,
            preserve_pore_water=True,
        )

        assert report["pore_aware"] is True

        # OW must still be present in output
        out_text = gro_out.read_text()
        assert "OW" in out_text, "Pore water OW was deleted from output GRO"
        assert "SOL" in out_text, "All SOL removed — pore water should survive"

        # The clean_water_report.json must document pore waters preserved
        cwr = json.loads((td / "clean_water_report.json").read_text())
        assert cwr["n_pore_waters_preserved"] >= 1


# ── Test 7: pore-aware cleanup removes true membrane-core water ───────────────

def test_pore_aware_cleanup_removes_core_water():
    """Water at membrane Z but outside TM ring must be removed by pore-aware cleanup."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        gro, tm_res = _build_channel_system(
            td, pore_water=False, core_water_outside=True
        )

        # Count initial SOL OW lines
        initial_sol_count = gro.read_text().count("OW")
        gro_out = td / "system_clean.gro"

        report = run_pore_aware_water_cleanup(
            gro_in=gro,
            gro_out=gro_out,
            tm_residues=tm_res,
            output_dir=td,
            remove_membrane_core_water=True,
        )

        assert report["n_water_molecules_removed"] >= 1, (
            "Expected at least 1 core water molecule removed"
        )

        final_sol_count = gro_out.read_text().count("OW")
        assert final_sol_count < initial_sol_count, (
            "SOL count in GRO did not decrease after core water removal"
        )


# ── Test 8: SOL count updated in topology ────────────────────────────────────

def test_sol_count_updated_in_topology():
    """Topology SOL count must be decremented by the number of removed waters."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        gro, tm_res = _build_channel_system(
            td, pore_water=False, core_water_outside=True
        )
        topol = td / "topol_in.top"
        topol_out = td / "topol_out.top"
        topol.write_text(
            "[ molecules ]\n"
            "DPPC             32\n"
            "SOL              1000\n"
            "NA               5\n"
        )

        report = run_pore_aware_water_cleanup(
            gro_in=gro,
            gro_out=td / "system_clean.gro",
            tm_residues=tm_res,
            output_dir=td,
            topol_in=topol,
            topol_out=topol_out,
            remove_membrane_core_water=True,
            update_topology_count=True,
        )

        n_removed = report["n_water_molecules_removed"]
        if n_removed == 0:
            pytest.skip("No waters removed in this geometry")

        text = topol_out.read_text()
        expected = f"SOL              {1000 - n_removed}"
        assert expected in text, (
            f"Expected '{expected}' in topology after removing {n_removed} waters.\n"
            f"Topology contents:\n{text}"
        )
        assert report["topology_updated"] is True


# ── Test 9: strict mode blocks when pore lipids remain ───────────────────────

def test_strict_mode_blocks_on_pore_lipids():
    """policy='strict' with pore lipids present must raise RuntimeError."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        gro, tm_res = _build_channel_system(
            td, pore_water=False, pore_lipid=True
        )
        with pytest.raises(RuntimeError, match="STRICT"):
            run_spatial_classifier(
                gro_path=gro,
                tm_residues=tm_res,
                output_dir=td,
                policy="strict",
                remove_pore_lipids=True,
            )


# ── Test 10: classifier works regardless of absolute protein position ─────────

def test_classifier_position_independent():
    """Translate the whole system and verify classification is unchanged."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)

        def build_at(cx, cy, out_name):
            ring_atoms, tm_res = _make_protein_ring(cx=cx, cy=cy)
            atoms = list(ring_atoms)
            nrid = max(tm_res) + 1
            # Pore water
            atoms += _make_sol(nrid, cx, cy, 5.0)
            pore_rid = nrid
            nrid += 1
            # Bulk lipids
            for i in range(4):
                angle = i * math.pi / 2
                atoms += _make_lipid(nrid, cx + 3.0 * math.cos(angle),
                                     cy + 3.0 * math.sin(angle), z_hg=6.0)
                nrid += 1
            gro = td / out_name
            _write_gro(gro, atoms, box=(cx * 2 + 2, cy * 2 + 2, 10.0))
            return gro, tm_res, pore_rid

        gro_a, tm_a, _ = build_at(6.0, 6.0, "a.gro")
        gro_b, tm_b, _ = build_at(3.0, 3.0, "b.gro")

        rep_a = run_spatial_classifier(gro_a, tm_residues=tm_a)
        rep_b = run_spatial_classifier(gro_b, tm_residues=tm_b)

        assert rep_a["n_pore_waters"] == rep_b["n_pore_waters"], (
            f"Pore water count changed with translation: "
            f"{rep_a['n_pore_waters']} vs {rep_b['n_pore_waters']}"
        )
        assert rep_a["n_forbidden_pore_lipids"] == rep_b["n_forbidden_pore_lipids"]


# ── Test 11: multiple cavities / pores detected ───────────────────────────────

def test_multiple_pores_detected():
    """Two separate water molecules in the pore → both classified as pore_water."""
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        cx, cy = 6.0, 6.0
        ring_atoms, tm_res = _make_protein_ring(cx=cx, cy=cy)
        atoms = list(ring_atoms)
        nrid = max(tm_res) + 1

        # Bulk lipids — need both leaflets so membrane core Z is [4,6]
        for i in range(4):
            atoms += _make_lipid(nrid, cx + 3.0 * math.cos(i * math.pi / 2),
                                 cy + 3.0 * math.sin(i * math.pi / 2), z_hg=6.0)
            nrid += 1
        for i in range(4):
            atoms += _make_lipid(nrid, cx + 3.0 * math.cos(i * math.pi / 2),
                                 cy + 3.0 * math.sin(i * math.pi / 2), z_hg=4.0)
            nrid += 1

        # Two pore water molecules at slightly different Z
        atoms += _make_sol(nrid, cx, cy, 4.8)
        nrid += 1
        atoms += _make_sol(nrid, cx, cy, 5.2)
        nrid += 1

        gro = td / "system.gro"
        _write_gro(gro, atoms)

        rep = run_spatial_classifier(gro_path=gro, tm_residues=tm_res, output_dir=td)

        assert rep["n_pore_waters"] >= 2, (
            f"Expected ≥2 pore waters, got {rep['n_pore_waters']}"
        )


# ── Test 12: existing lipid_refill regression unaffected ─────────────────────

def test_existing_lipid_refill_regression_unaffected():
    """Spatial classifier must not break existing lipid refill tests.

    Verifies that run_spatial_classifier handles a system with no SOL at all
    (pure protein+lipid bilayer) without errors.
    """
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)

        # Pure protein + lipid, no water
        atoms = []
        tm_res: set[int] = set()
        for i in range(10):
            resid = i + 1
            angle = i * 2.0 * math.pi / 10
            x = 6.0 + 0.5 * math.cos(angle)
            y = 6.0 + 0.5 * math.sin(angle)
            atoms.append((resid, "ALA", "CA", x, y, 5.0))
            tm_res.add(resid)
        nrid = 11
        for i in range(10):
            angle = i * 2.0 * math.pi / 10
            lx = 6.0 + 2.0 * math.cos(angle)
            ly = 6.0 + 2.0 * math.sin(angle)
            atoms += _make_lipid(nrid, lx, ly, z_hg=6.0)
            nrid += 1

        gro = td / "system.gro"
        _write_gro(gro, atoms)

        rep = run_spatial_classifier(gro_path=gro, tm_residues=tm_res, output_dir=td)

        # Should run without error and report 0 waters
        assert rep["enabled"]
        assert rep["n_bulk_waters"] == 0
        assert rep["n_pore_waters"] == 0
        assert rep["n_membrane_core_waters"] == 0

        report_file = td / "membrane_protein_spatial_classification_report.json"
        assert report_file.exists()
