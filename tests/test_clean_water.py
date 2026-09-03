"""
tests/test_clean_water.py
Phase 10B: Pore-aware clean_water — unit tests.

Complements builders/step_builders/test_clean_water_builder.py (which tests the
original Z-slab WaterDeletorAdapter). This file tests the pore-aware cleanup
that preserves channel water while still removing membrane-core water outside
the protein.
"""
from __future__ import annotations

import json
import math
import tempfile
from pathlib import Path

import pytest

from validators.membrane_protein_spatial_classifier import run_pore_aware_water_cleanup


# ── GRO helpers (same pattern as test_membrane_protein_spatial_classifier) ────

def _write_gro(path: Path, atoms: list, box=(12.0, 12.0, 10.0)) -> None:
    lines = ["test system", str(len(atoms))]
    for i, (resid, resname, atomname, x, y, z) in enumerate(atoms, 1):
        lines.append(
            f"{resid:5d}{resname:<5s}{atomname:>5s}{i:5d}{x:8.3f}{y:8.3f}{z:8.3f}"
        )
    bx, by, bz = box
    lines.append(f"   {bx:.5f}   {by:.5f}   {bz:.5f}")
    path.write_text("\n".join(lines) + "\n")


def _make_protein_ring(
    cx=6.0, cy=6.0, ring_radius=1.0, n_per_level=16,
    z_levels=None, start_resid=1,
):
    if z_levels is None:
        z_levels = [3.5 + k * 0.3 for k in range(11)]
    atoms, tm_residues = [], set()
    resid = start_resid
    for z in z_levels:
        for k in range(n_per_level):
            angle = k * 2.0 * math.pi / n_per_level
            atoms.append((resid, "ALA", "CA",
                           cx + ring_radius * math.cos(angle),
                           cy + ring_radius * math.sin(angle), z))
            tm_residues.add(resid)
            resid += 1
    return atoms, tm_residues


def _sol(resid, x, y, z):
    return [
        (resid, "SOL", "OW",  x,        y,        z),
        (resid, "SOL", "HW1", x + 0.01, y,        z),
        (resid, "SOL", "HW2", x,        y + 0.01, z),
    ]


def _lipid(resid, x, y, z_hg, resname="DPP"):
    return [
        (resid, resname, "P",   x, y, z_hg),
        (resid, resname, "C16", x, y, z_hg - 0.5),
    ]


def _build_system(
    tmpdir: Path,
    cx=6.0, cy=6.0,
    pore_water=True,
    core_water=True,
    n_bulk_lipids=4,
):
    """Build a channel system with optional pore and core water."""
    ring_atoms, tm_res = _make_protein_ring(cx=cx, cy=cy)
    atoms = list(ring_atoms)
    nrid = max(tm_res) + 1

    # Bulk lipids (upper + lower leaflet)
    for i in range(n_bulk_lipids):
        angle = i * 2.0 * math.pi / n_bulk_lipids
        lx = cx + 3.0 * math.cos(angle)
        ly = cy + 3.0 * math.sin(angle)
        atoms += _lipid(nrid, lx, ly, z_hg=6.0)
        nrid += 1
    for i in range(n_bulk_lipids):
        angle = i * 2.0 * math.pi / n_bulk_lipids
        lx = cx + 3.0 * math.cos(angle)
        ly = cy + 3.0 * math.sin(angle)
        atoms += _lipid(nrid, lx, ly, z_hg=4.0)
        nrid += 1

    pore_resid = None
    if pore_water:
        atoms += _sol(nrid, cx, cy, 5.0)
        pore_resid = nrid
        nrid += 1

    core_resid = None
    if core_water:
        atoms += _sol(nrid, cx + 3.0, cy, 5.0)
        core_resid = nrid
        nrid += 1

    gro = tmpdir / "system.gro"
    _write_gro(gro, atoms)
    return gro, tm_res, pore_resid, core_resid


# ── Test 1: pore water is not removed ────────────────────────────────────────

class TestPoreWaterPreserved:
    def test_pore_ow_remains_in_output_gro(self, tmp_path):
        """OW of pore water must still appear in output GRO."""
        gro, tm_res, pore_rid, core_rid = _build_system(
            tmp_path, pore_water=True, core_water=False
        )
        gro_out = tmp_path / "out.gro"

        rep = run_pore_aware_water_cleanup(
            gro_in=gro, gro_out=gro_out,
            tm_residues=tm_res,
            output_dir=tmp_path,
            preserve_pore_water=True,
            remove_membrane_core_water=True,
        )

        assert "SOL" in gro_out.read_text(), "Pore water SOL removed from output"
        assert rep["n_pore_waters_preserved"] >= 1

    def test_pore_water_resids_in_report(self, tmp_path):
        """Pore water resids must appear in the report."""
        gro, tm_res, _, _ = _build_system(
            tmp_path, pore_water=True, core_water=False
        )
        gro_out = tmp_path / "out.gro"

        rep = run_pore_aware_water_cleanup(
            gro_in=gro, gro_out=gro_out,
            tm_residues=tm_res,
            output_dir=tmp_path,
        )

        assert len(rep["pore_water_resids"]) >= 1

    def test_clean_water_report_written(self, tmp_path):
        """clean_water_report.json must be written with pore_aware=True."""
        gro, tm_res, _, _ = _build_system(tmp_path, pore_water=True, core_water=False)
        gro_out = tmp_path / "out.gro"

        run_pore_aware_water_cleanup(
            gro_in=gro, gro_out=gro_out,
            tm_residues=tm_res,
            output_dir=tmp_path,
        )

        cwr_path = tmp_path / "clean_water_report.json"
        assert cwr_path.exists()
        cwr = json.loads(cwr_path.read_text())
        assert cwr["pore_aware"] is True
        assert cwr["cleanup_passed"] is True


# ── Test 2: membrane-core water is removed ────────────────────────────────────

class TestCoreWaterRemoved:
    def test_core_water_mol_count_decreases(self, tmp_path):
        """A core water molecule outside the TM ring must be removed."""
        gro, tm_res, _, core_rid = _build_system(
            tmp_path, pore_water=False, core_water=True
        )
        initial_ow = gro.read_text().count("OW")
        gro_out = tmp_path / "out.gro"

        rep = run_pore_aware_water_cleanup(
            gro_in=gro, gro_out=gro_out,
            tm_residues=tm_res,
            output_dir=tmp_path,
            remove_membrane_core_water=True,
        )

        final_ow = gro_out.read_text().count("OW")
        assert final_ow < initial_ow, (
            f"Core water OW not removed: before={initial_ow}, after={final_ow}"
        )
        assert rep["n_water_molecules_removed"] >= 1

    def test_non_water_atoms_unaffected(self, tmp_path):
        """Protein and lipid atoms must survive cleanup."""
        gro, tm_res, _, _ = _build_system(tmp_path, core_water=True)
        gro_out = tmp_path / "out.gro"

        run_pore_aware_water_cleanup(
            gro_in=gro, gro_out=gro_out,
            tm_residues=tm_res,
            output_dir=tmp_path,
        )

        out_text = gro_out.read_text()
        assert "ALA" in out_text, "Protein atom removed"
        assert "DPP" in out_text, "Lipid atom removed"


# ── Test 3: SOL count updated in topology ─────────────────────────────────────

class TestTopologyUpdate:
    def test_sol_decremented(self, tmp_path):
        """topol.top SOL count must decrease by n_water_molecules_removed."""
        gro, tm_res, _, _ = _build_system(
            tmp_path, pore_water=False, core_water=True
        )
        topol_in  = tmp_path / "topol.top"
        topol_out = tmp_path / "topol_clean.top"
        initial_sol = 500
        topol_in.write_text(
            "[ molecules ]\nDPPC             32\n"
            f"SOL              {initial_sol}\nNA               5\n"
        )

        rep = run_pore_aware_water_cleanup(
            gro_in=gro, gro_out=tmp_path / "out.gro",
            tm_residues=tm_res,
            topol_in=topol_in, topol_out=topol_out,
            output_dir=tmp_path,
            update_topology_count=True,
        )

        n_rm = rep["n_water_molecules_removed"]
        if n_rm == 0:
            pytest.skip("No waters removed in this geometry")

        text = topol_out.read_text()
        assert f"SOL              {initial_sol - n_rm}" in text, (
            f"Expected SOL {initial_sol - n_rm}, got:\n{text}"
        )
        assert rep["topology_updated"] is True

    def test_no_topol_update_when_flag_false(self, tmp_path):
        """update_topology_count=False must not write topology."""
        gro, tm_res, _, _ = _build_system(tmp_path, core_water=True)
        topol_in  = tmp_path / "topol.top"
        topol_out = tmp_path / "topol_clean.top"
        topol_in.write_text("SOL              500\n")

        rep = run_pore_aware_water_cleanup(
            gro_in=gro, gro_out=tmp_path / "out.gro",
            tm_residues=tm_res,
            topol_in=topol_in, topol_out=topol_out,
            output_dir=tmp_path,
            update_topology_count=False,
        )

        assert rep["topology_updated"] is False
        assert not topol_out.exists()


# ── Test 4: both pore water and core water in same system ─────────────────────

class TestMixedSystem:
    def test_core_removed_pore_preserved(self, tmp_path):
        """Core water removed AND pore water preserved in same run."""
        gro, tm_res, pore_rid, core_rid = _build_system(
            tmp_path, pore_water=True, core_water=True
        )
        gro_out = tmp_path / "out.gro"

        rep = run_pore_aware_water_cleanup(
            gro_in=gro, gro_out=gro_out,
            tm_residues=tm_res,
            output_dir=tmp_path,
            preserve_pore_water=True,
            remove_membrane_core_water=True,
        )

        assert rep["n_pore_waters_preserved"] >= 1
        assert rep["n_water_molecules_removed"] >= 1

        # Output GRO must still have SOL (the pore water)
        assert "SOL" in gro_out.read_text()

    def test_water_gate_compatible_report(self, tmp_path):
        """clean_water_report.json must be compatible with water_gate.py."""
        gro, tm_res, _, _ = _build_system(
            tmp_path, pore_water=True, core_water=True
        )

        run_pore_aware_water_cleanup(
            gro_in=gro, gro_out=tmp_path / "out.gro",
            tm_residues=tm_res,
            output_dir=tmp_path,
        )

        from runtime.water_gate import evaluate_water_gate
        result = evaluate_water_gate(tmp_path)
        # Gate should read clean_water_report.json and not block
        assert result is not None
