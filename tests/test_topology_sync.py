"""
tests/test_topology_sync.py

Regression tests for the membrane_embedding coordinate/topology sync fix.

Real-run evidence this fixes (see task description / docs/audits):
    work.gro:            30421 atoms, 477 DPP residues
    bootstrap_topol.top:  478 DPP  ->  implies 30471 atoms
    grompp fatal error: coordinate file has 30421 atoms, topology expects 30471

The mismatch is exactly 50 atoms = 1 DPP molecule (OPLS-AA DPPC has 50 atoms
per residue).
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from validators.topology_sync import sync_topology_lipid_count


# ── GRO/topology fixture helpers ───────────────────────────────────────────────

def _make_gro_line(resid: int, resname: str, atomname: str, atomnum: int) -> str:
    return f"{resid:>5}{resname:<5}{atomname:>5}{atomnum:>5}{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}"


def _write_synthetic_system(
    path: Path,
    n_protein_atoms: int = 6571,
    n_dpp_residues: int = 477,
    atoms_per_dpp: int = 50,
    lipid_resname: str = "DPP",
) -> int:
    """Write a synthetic work.gro with n_protein_atoms protein atoms and
    n_dpp_residues lipid residues (atoms_per_dpp atoms each). Returns total
    atom count written."""
    lines: list[str] = []
    atomnum = 0
    for _ in range(n_protein_atoms):
        atomnum += 1
        lines.append(_make_gro_line(1, "PRO", "CA", atomnum))
    for lip in range(n_dpp_residues):
        resid = 7000 + lip
        for a in range(atoms_per_dpp):
            atomnum += 1
            lines.append(_make_gro_line(resid, lipid_resname, f"C{a}", atomnum))

    out = ["synthetic sync-test system", str(len(lines))] + lines + ["  10.00000  10.00000  10.00000"]
    path.write_text("\n".join(out) + "\n")
    return len(lines)


def _write_topol(path: Path, protein_count: int = 1, dpp_count: int = 478,
                  lipid_resname: str = "DPP", extra_lines: list[str] | None = None) -> None:
    lines = [
        "; synthetic bootstrap topology",
        "[ system ]",
        "Protein-membrane system",
        "",
        "[ molecules ]",
        f"Protein_chain_A      {protein_count}",
        f"{lipid_resname:<20s} {dpp_count}",
    ]
    if extra_lines:
        lines.extend(extra_lines)
    path.write_text("\n".join(lines) + "\n")


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Core sync algorithm — matches the exact real-run evidence numbers
# ═══════════════════════════════════════════════════════════════════════════════

class TestSyncTopologyLipidCount:

    def test_sync_matches_real_run_evidence_numbers(self, tmp_path):
        """
        478 DPP in topology, 477 DPP in work.gro (as in the real failing run).
        Sync must update topology to 477 and the resulting expected atom
        count must match the coordinate file's actual atom count.
        """
        gro = tmp_path / "work.gro"
        n_atoms = _write_synthetic_system(gro, n_protein_atoms=6571, n_dpp_residues=477)
        assert n_atoms == 30421, "fixture sanity check: expected 30421 atoms"

        topol = tmp_path / "bootstrap_topol.top"
        _write_topol(topol, protein_count=1, dpp_count=478)

        report = sync_topology_lipid_count(gro, topol, "DPP")

        assert report["old_dpp_count"] == 478
        assert report["new_dpp_count"] == 477
        assert report["removed_dpp_count"] == 1
        assert report["coordinate_atoms"] == 30421
        assert report["old_expected_atoms"] == 30471
        assert report["new_expected_atoms"] == 30421
        assert report["new_expected_atoms"] == report["coordinate_atoms"], (
            "after sync, topology-implied atom count must match the coordinate file"
        )
        assert report["synchronized"] is True

        # Topology file itself was rewritten.
        updated_text = topol.read_text()
        assert "DPP             477" in updated_text or "DPP              477" in updated_text

    def test_protein_count_preserved(self, tmp_path):
        """Only the lipid line is touched — Protein count must survive untouched."""
        gro = tmp_path / "work.gro"
        _write_synthetic_system(gro, n_protein_atoms=6571, n_dpp_residues=477)
        topol = tmp_path / "bootstrap_topol.top"
        _write_topol(topol, protein_count=1, dpp_count=478)

        sync_topology_lipid_count(gro, topol, "DPP")

        updated_lines = topol.read_text().splitlines()
        protein_line = next(l for l in updated_lines if l.strip().startswith("Protein_chain_A"))
        assert protein_line.split()[1] == "1"

    def test_sol_and_ions_not_invented(self, tmp_path):
        """If SOL/ions are absent from the topology, sync must not add them."""
        gro = tmp_path / "work.gro"
        _write_synthetic_system(gro, n_protein_atoms=6571, n_dpp_residues=477)
        topol = tmp_path / "bootstrap_topol.top"
        _write_topol(topol, protein_count=1, dpp_count=478)  # no SOL/ion lines

        sync_topology_lipid_count(gro, topol, "DPP")

        updated_text = topol.read_text()
        assert "SOL" not in updated_text
        assert "NA" not in updated_text
        assert "CL" not in updated_text

    def test_no_change_needed_is_still_synchronized(self, tmp_path):
        """When topology already matches work.gro, no rewrite is needed but
        the report must still say synchronized=True."""
        gro = tmp_path / "work.gro"
        _write_synthetic_system(gro, n_protein_atoms=6571, n_dpp_residues=478)
        topol = tmp_path / "bootstrap_topol.top"
        original_text = None
        _write_topol(topol, protein_count=1, dpp_count=478)
        original_text = topol.read_text()

        report = sync_topology_lipid_count(gro, topol, "DPP")

        assert report["old_dpp_count"] == 478
        assert report["new_dpp_count"] == 478
        assert report["removed_dpp_count"] == 0
        assert report["synchronized"] is True
        assert topol.read_text() == original_text

    def test_lipid_entry_absent_is_not_invented(self, tmp_path):
        """If the lipid isn't in [ molecules ] at all, do not invent an entry."""
        gro = tmp_path / "work.gro"
        _write_synthetic_system(gro, n_protein_atoms=6571, n_dpp_residues=477)
        topol = tmp_path / "bootstrap_topol.top"
        topol.write_text(
            "[ molecules ]\n"
            "Protein_chain_A      1\n"
        )

        report = sync_topology_lipid_count(gro, topol, "DPP")

        assert report["old_dpp_count"] is None
        assert report["synchronized"] is False
        assert "DPP" not in topol.read_text()

    def test_moleculetype_section_not_confused_with_molecules_section(self, tmp_path):
        """
        A '[ moleculetype ]' header line like 'DPP   3' (nrexcl) inlined in the
        same topology file must not be mistaken for the [ molecules ] count.
        """
        gro = tmp_path / "work.gro"
        _write_synthetic_system(gro, n_protein_atoms=6571, n_dpp_residues=477)
        topol = tmp_path / "bootstrap_topol.top"
        topol.write_text(
            "[ moleculetype ]\n"
            "; name  nrexcl\n"
            "DPP   3\n"
            "\n"
            "[ molecules ]\n"
            "Protein_chain_A      1\n"
            "DPP                  478\n"
        )

        report = sync_topology_lipid_count(gro, topol, "DPP")

        assert report["old_dpp_count"] == 478, (
            "must read the [ molecules ] count (478), not the moleculetype nrexcl (3)"
        )
        assert report["new_dpp_count"] == 477

        updated_lines = topol.read_text().splitlines()
        # The moleculetype nrexcl line must be untouched.
        assert "DPP   3" in updated_lines

    def test_missing_coordinate_file_reports_unsynchronized(self, tmp_path):
        topol = tmp_path / "bootstrap_topol.top"
        _write_topol(topol, dpp_count=478)
        report = sync_topology_lipid_count(tmp_path / "missing.gro", topol, "DPP")
        assert report["synchronized"] is False
        assert report["warnings"]

    def test_missing_topology_file_reports_unsynchronized(self, tmp_path):
        gro = tmp_path / "work.gro"
        _write_synthetic_system(gro, n_dpp_residues=477)
        report = sync_topology_lipid_count(gro, tmp_path / "missing_topol.top", "DPP")
        assert report["synchronized"] is False
        assert report["warnings"]


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Report contract — required fields from the task spec
# ═══════════════════════════════════════════════════════════════════════════════

class TestSyncReportContract:

    REQUIRED_KEYS = {
        "coordinate_file", "topology_file", "old_dpp_count", "new_dpp_count",
        "removed_dpp_count", "old_expected_atoms", "new_expected_atoms",
        "coordinate_atoms", "synchronized",
    }

    def test_report_has_all_required_fields(self, tmp_path):
        gro = tmp_path / "work.gro"
        _write_synthetic_system(gro, n_dpp_residues=477)
        topol = tmp_path / "bootstrap_topol.top"
        _write_topol(topol, dpp_count=478)

        report = sync_topology_lipid_count(gro, topol, "DPP")
        assert self.REQUIRED_KEYS.issubset(set(report.keys())), (
            f"missing keys: {self.REQUIRED_KEYS - set(report.keys())}"
        )

    def test_report_written_to_disk(self, tmp_path):
        gro = tmp_path / "work.gro"
        _write_synthetic_system(gro, n_dpp_residues=477)
        topol = tmp_path / "bootstrap_topol.top"
        _write_topol(topol, dpp_count=478)

        sync_topology_lipid_count(gro, topol, "DPP")

        report_path = topol.parent / "topology_sync_report.json"
        assert report_path.exists()
        disk_report = json.loads(report_path.read_text())
        assert disk_report["old_dpp_count"] == 478
        assert disk_report["new_dpp_count"] == 477

    def test_report_path_override(self, tmp_path):
        gro = tmp_path / "work.gro"
        _write_synthetic_system(gro, n_dpp_residues=477)
        topol = tmp_path / "bootstrap_topol.top"
        _write_topol(topol, dpp_count=478)
        custom_report = tmp_path / "custom_sync_report.json"

        sync_topology_lipid_count(gro, topol, "DPP", report_path=custom_report)
        assert custom_report.exists()
