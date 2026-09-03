"""
Tests for ligand/campaign.py — multi-complex ligand campaign discovery,
pose extraction, identity/dedup, and manifest generation.

Hydrogenation is exercised through the real ``ligand.prepare`` /
``ligand.hydrogenation`` pipeline with RDKit/Open Babel availability mocked
(same convention as ligand/test_integrate.py), so these tests are fast and
deterministic regardless of what's installed in the environment.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from ligand.campaign import (
    CampaignDiscoveryError,
    CampaignManifest,
    ChemicalIdentity,
    ComplexEntry,
    ExtractedComplex,
    LigandSignature,
    compute_ligand_signature,
    discover_complexes,
    extract_complex_components,
    group_ligand_occurrences,
    prepare_campaign,
)

_RDKIT_AVAIL = "ligand.hydrogenation._rdkit_available"
_OBABEL_AVAIL = "ligand.hydrogenation._obabel_available"

# No hydrogenation backend available -> ManualRequiredBackend -> status stays
# "missing" (ligand fixtures below have no H atoms). Deterministic across envs.
_NO_HYDROGENATION_BACKEND = [
    patch(_RDKIT_AVAIL, return_value=False),
    patch(_OBABEL_AVAIL, return_value=False),
]


def _no_backend():
    from contextlib import ExitStack

    stack = ExitStack()
    for p in _NO_HYDROGENATION_BACKEND:
        stack.enter_context(p)
    return stack


# ═══════════════════════════════════════════════════════════════════════════
# PDB fixture builders
# ═══════════════════════════════════════════════════════════════════════════

def _pdb_line(
    record: str, serial: int, atom_name: str, resname: str, chain: str,
    resseq: int, x: float, y: float, z: float, element: str,
) -> str:
    return (
        f"{record:<6}{serial:>5} {atom_name:>4} {resname:>3} {chain:1}{resseq:>4}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {element:>2}"
    )


def _protein_lines(n_res: int = 2) -> list[str]:
    lines = []
    serial = 1
    for r in range(1, n_res + 1):
        for name, element in (("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")):
            lines.append(_pdb_line("ATOM", serial, name, "ALA", "A", r, r * 1.0, 0.0, 0.0, element))
            serial += 1
    return lines


def _hexagon_lines(
    resname: str, serial_start: int = 100, resseq: int = 900,
    center: tuple[float, float, float] = (10.0, 0.0, 0.0),
    elements: tuple[str, ...] | None = None,
) -> list[str]:
    """Six-membered ring HETATM block (heavy atoms only), no hydrogens."""
    import math

    cx, cy, cz = center
    n = 6
    elements = elements or tuple("C" for _ in range(n))
    lines = []
    for i in range(n):
        angle = 2 * math.pi * i / n
        x = cx + 1.2 * math.cos(angle)
        y = cy + 1.2 * math.sin(angle)
        z = cz
        lines.append(_pdb_line(
            "HETATM", serial_start + i, f"C{i+1}", resname, "A", resseq,
            x, y, z, elements[i],
        ))
    return lines


def _write_complex(path: Path, resname: str, center=(10.0, 0.0, 0.0), n_res=2, elements=None) -> Path:
    lines = _protein_lines(n_res) + ["TER"] + _hexagon_lines(resname, center=center, elements=elements) + ["END"]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


# ═══════════════════════════════════════════════════════════════════════════
# Phase 1 — discovery
# ═══════════════════════════════════════════════════════════════════════════

class TestDiscoverComplexes:
    def test_four_directories_four_pdbs(self, tmp_path):
        for sid in ("AA", "AG", "HMG-R", "LP"):
            _write_complex(tmp_path / sid / f"{sid}.pdb", "A6")

        entries = discover_complexes(tmp_path)
        assert [e.system_id for e in entries] == ["AA", "AG", "HMG-R", "LP"]

    def test_deterministic_order(self, tmp_path):
        for sid in ("LP", "AA", "HMG-R", "AG"):
            _write_complex(tmp_path / sid / f"{sid}.pdb", "A6")

        first = [e.system_id for e in discover_complexes(tmp_path)]
        second = [e.system_id for e in discover_complexes(tmp_path)]
        assert first == second == sorted(first)

    def test_ambiguous_directory_raises(self, tmp_path):
        _write_complex(tmp_path / "AA" / "one.pdb", "A6")
        _write_complex(tmp_path / "AA" / "two.pdb", "A6")

        with pytest.raises(CampaignDiscoveryError, match="Ambiguous"):
            discover_complexes(tmp_path)

    def test_ignores_simforge_output_dirs(self, tmp_path):
        entry = _write_complex(tmp_path / "AA" / "AA.pdb", "A6")
        # Simulate a prior run's output directory containing a stray PDB.
        (tmp_path / "AA" / "simforge").mkdir()
        (tmp_path / "AA" / "simforge" / "protein_only.pdb").write_text("END\n")
        (tmp_path / "simforge_campaign").mkdir()
        (tmp_path / "simforge_campaign" / "stray.pdb").write_text("END\n")

        entries = discover_complexes(tmp_path)
        assert [e.system_id for e in entries] == ["AA"]
        assert entries[0].complex_pdb == entry

    def test_missing_complexes_returns_empty(self, tmp_path):
        (tmp_path / "empty_dir").mkdir()
        assert discover_complexes(tmp_path) == []

    def test_missing_root_raises(self, tmp_path):
        with pytest.raises(CampaignDiscoveryError):
            discover_complexes(tmp_path / "does_not_exist")


# ═══════════════════════════════════════════════════════════════════════════
# Phase 2 — pose extraction
# ═══════════════════════════════════════════════════════════════════════════

class TestExtractComplexComponents:
    def test_extracts_protein_and_ligand(self, tmp_path):
        complex_pdb = _write_complex(tmp_path / "complex.pdb", "A6")
        result = extract_complex_components(complex_pdb, "A6", tmp_path / "out")

        assert result.success
        assert result.protein_pdb.exists()
        assert result.ligand_pose_pdb.exists()
        assert result.ligand_atom_count == 6
        assert result.protein_atom_count == 8  # 2 residues * 4 atoms

    def test_preserves_original_heavy_atom_coordinates(self, tmp_path):
        complex_pdb = _write_complex(tmp_path / "complex.pdb", "A6", center=(30.0, 5.0, 2.0))
        result = extract_complex_components(complex_pdb, "A6", tmp_path / "out")

        original_coords = []
        for line in complex_pdb.read_text().splitlines():
            if line.startswith("HETATM"):
                original_coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))

        extracted_coords = []
        for line in result.ligand_pose_pdb.read_text().splitlines():
            if line.startswith("HETATM"):
                extracted_coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))

        assert extracted_coords == original_coords

    def test_no_hydrogenation_performed(self, tmp_path):
        """extract_complex_components must not modify the ligand chemically."""
        complex_pdb = _write_complex(tmp_path / "complex.pdb", "A6")
        result = extract_complex_components(complex_pdb, "A6", tmp_path / "out")
        assert result.n_hydrogen_atoms == 0
        assert result.n_heavy_atoms == 6

    def test_missing_resname_fails(self, tmp_path):
        complex_pdb = _write_complex(tmp_path / "complex.pdb", "A6")
        result = extract_complex_components(complex_pdb, "ZZZ", tmp_path / "out")
        assert not result.success
        assert "ZZZ" in result.error

    def test_missing_complex_file_fails(self, tmp_path):
        result = extract_complex_components(tmp_path / "nonexistent.pdb", "A6", tmp_path / "out")
        assert not result.success
        assert "not found" in result.error


# ═══════════════════════════════════════════════════════════════════════════
# Phase 3 — ligand identity / signatures
# ═══════════════════════════════════════════════════════════════════════════

class TestLigandSignature:
    def test_same_signature_matches(self, tmp_path):
        c1 = _write_complex(tmp_path / "c1.pdb", "A6", center=(10, 0, 0))
        c2 = _write_complex(tmp_path / "c2.pdb", "A6", center=(30, 5, 2))
        e1 = extract_complex_components(c1, "A6", tmp_path / "o1")
        e2 = extract_complex_components(c2, "A6", tmp_path / "o2")

        sig1 = compute_ligand_signature(e1.ligand_pose_pdb, "A6")
        sig2 = compute_ligand_signature(e2.ligand_pose_pdb, "A6")
        assert sig1.matches(sig2)

    def test_different_heavy_atom_count_does_not_match(self, tmp_path):
        c1 = _write_complex(tmp_path / "c1.pdb", "A6")
        # Seven heavy atoms instead of six.
        lines = _protein_lines() + ["TER"] + _hexagon_lines("A6", center=(10, 0, 0)) + [
            _pdb_line("HETATM", 200, "C7", "A6", "A", 900, 12.0, 0.0, 0.0, "C"),
            "END",
        ]
        c2 = tmp_path / "c2.pdb"
        c2.write_text("\n".join(lines) + "\n")

        e1 = extract_complex_components(c1, "A6", tmp_path / "o1")
        e2 = extract_complex_components(c2, "A6", tmp_path / "o2")
        sig1 = compute_ligand_signature(e1.ligand_pose_pdb, "A6")
        sig2 = compute_ligand_signature(e2.ligand_pose_pdb, "A6")
        assert not sig1.matches(sig2)

    def test_different_element_composition_does_not_match(self, tmp_path):
        c1 = _write_complex(tmp_path / "c1.pdb", "A6")
        c2 = _write_complex(
            tmp_path / "c2.pdb", "A6",
            elements=("C", "C", "C", "C", "C", "N"),  # one C -> N
        )
        e1 = extract_complex_components(c1, "A6", tmp_path / "o1")
        e2 = extract_complex_components(c2, "A6", tmp_path / "o2")
        sig1 = compute_ligand_signature(e1.ligand_pose_pdb, "A6")
        sig2 = compute_ligand_signature(e2.ligand_pose_pdb, "A6")
        assert not sig1.matches(sig2)


class TestGroupLigandOccurrences:
    def test_same_signature_one_group(self):
        sig = LigandSignature("A6", 6, ("C1", "C2", "C3", "C4", "C5", "C6"),
                               ("C", "C", "C", "C", "C", "C"))
        occurrences = [("AA", sig), ("AG", sig), ("LP", sig)]
        groups, warnings = group_ligand_occurrences(occurrences, "A6")
        assert list(groups.keys()) == ["A6"]
        assert groups["A6"].system_ids == ["AA", "AG", "LP"]
        assert not warnings

    def test_divergent_signature_splits_group(self):
        sig_a = LigandSignature("A6", 6, ("C1",) * 1 + ("C2", "C3", "C4", "C5", "C6"),
                                 ("C", "C", "C", "C", "C", "C"))
        sig_b = LigandSignature("A6", 7, ("C1", "C2", "C3", "C4", "C5", "C6", "C7"),
                                 ("C", "C", "C", "C", "C", "C", "C"))
        occurrences = [("AA", sig_a), ("AG", sig_b)]
        groups, warnings = group_ligand_occurrences(occurrences, "A6")
        assert set(groups.keys()) == {"A6", "A6_2"}
        assert groups["A6"].system_ids == ["AA"]
        assert groups["A6_2"].system_ids == ["AG"]
        assert warnings

    def test_explicit_chemical_reference_groups_all_matching(self):
        sig = LigandSignature("A6", 6, (), ("C", "C", "C", "C", "C", "C"))
        chem = ChemicalIdentity(source_path=Path("A6.sdf"), heavy_atom_count=6, elements=("C",) * 6)
        occurrences = [("AA", sig), ("AG", sig)]
        groups, warnings = group_ligand_occurrences(occurrences, "A6", chemical_identity=chem)
        assert list(groups.keys()) == ["A6"]
        assert groups["A6"].identity_method == "chemical_reference"
        assert groups["A6"].identity_confidence == "high"
        assert not warnings

    def test_explicit_chemical_reference_excludes_inconsistent_occurrence(self):
        good = LigandSignature("A6", 6, (), ("C",) * 6)
        bad = LigandSignature("A6", 5, (), ("C",) * 5)
        chem = ChemicalIdentity(source_path=Path("A6.sdf"), heavy_atom_count=6, elements=("C",) * 6)
        occurrences = [("AA", good), ("AG", bad)]
        groups, warnings = group_ligand_occurrences(occurrences, "A6", chemical_identity=chem)
        assert groups["A6"].system_ids == ["AA"]
        assert any("AG" in w for w in warnings)


# ═══════════════════════════════════════════════════════════════════════════
# compute_chemical_identity — graceful RDKit-absent behavior
# ═══════════════════════════════════════════════════════════════════════════

def test_compute_chemical_identity_raises_importerror_without_rdkit():
    from ligand.campaign import compute_chemical_identity

    with patch("ligand.rdkit_reader.load_mol", side_effect=ImportError("RDKit is required")):
        with pytest.raises(ImportError):
            compute_chemical_identity(Path("A6.sdf"))


# ═══════════════════════════════════════════════════════════════════════════
# Phase 4-7 — prepare_campaign orchestration
# ═══════════════════════════════════════════════════════════════════════════

class TestPrepareCampaign:
    def _campaign_root(self, tmp_path, system_ids=("AA", "AG", "HMG-R", "LP")):
        for sid in system_ids:
            _write_complex(tmp_path / sid / f"{sid}.pdb", "A6")
        return tmp_path

    def test_extracts_all_proteins_and_poses(self, tmp_path):
        root = self._campaign_root(tmp_path)
        with _no_backend():
            result = prepare_campaign(root, "A6")

        assert result.success, result.error
        assert result.systems_discovered == 4
        for sid in ("AA", "AG", "HMG-R", "LP"):
            assert (root / sid / "simforge" / "protein_only.pdb").exists()
            assert (root / sid / "simforge" / "ligand_pose.pdb").exists()

    def test_only_one_hydrogenation_call_for_four_identical_ligands(self, tmp_path):
        root = self._campaign_root(tmp_path)
        with _no_backend(), patch(
            "ligand.campaign.extract_ligand_from_complex",
            side_effect=__import__("ligand.prepare", fromlist=["extract_ligand_from_complex"]).extract_ligand_from_complex,
        ) as mocked:
            result = prepare_campaign(root, "A6")

        assert result.success, result.error
        assert mocked.call_count == 1

    def test_only_one_ligpargen_input_written(self, tmp_path):
        root = self._campaign_root(tmp_path)
        with _no_backend():
            result = prepare_campaign(root, "A6")

        param_files = list((root / "simforge_campaign" / "ligands").glob("*/parameterization/ligand_for_ligpargen.pdb"))
        assert len(param_files) == 1

    def test_manifest_contains_all_systems(self, tmp_path):
        root = self._campaign_root(tmp_path)
        with _no_backend():
            result = prepare_campaign(root, "A6")

        assert len(result.manifest.systems) == 4
        assert {s.id for s in result.manifest.systems} == {"AA", "AG", "HMG-R", "LP"}
        assert result.unique_ligands == 1

    def test_manifest_uses_relative_paths(self, tmp_path):
        root = self._campaign_root(tmp_path)
        with _no_backend():
            result = prepare_campaign(root, "A6")

        for s in result.manifest.systems:
            assert not Path(s.protein).is_absolute()
            assert not Path(s.ligand_pose).is_absolute()
            assert not Path(s.complex_pdb).is_absolute()

    def test_missing_hydrogens_recorded(self, tmp_path):
        root = self._campaign_root(tmp_path, system_ids=("AA",))
        with _no_backend():
            result = prepare_campaign(root, "A6")

        lig = result.manifest.ligands["A6"]
        assert lig.hydrogenation_status == "missing"
        assert lig.hydrogenation_performed is False
        assert lig.parameterization_status == "chemical_validation_required"

    def test_manifest_json_round_trip(self, tmp_path):
        root = self._campaign_root(tmp_path, system_ids=("AA", "AG"))
        with _no_backend():
            result = prepare_campaign(root, "A6")

        loaded = CampaignManifest.load(result.manifest_path)
        assert loaded.campaign_id == "A6"
        assert len(loaded.systems) == 2

    def test_report_md_written(self, tmp_path):
        root = self._campaign_root(tmp_path, system_ids=("AA", "AG"))
        with _no_backend():
            result = prepare_campaign(root, "A6")

        assert result.report_path.exists()
        text = result.report_path.read_text()
        assert "A6" in text

    def test_resume_reuses_existing_manifest(self, tmp_path):
        root = self._campaign_root(tmp_path, system_ids=("AA", "AG"))
        with _no_backend():
            first = prepare_campaign(root, "A6")
            second = prepare_campaign(root, "A6")

        assert first.success and second.success
        assert second.resumed is True
        assert second.manifest.systems == first.manifest.systems

    def test_resume_detects_missing_files(self, tmp_path):
        root = self._campaign_root(tmp_path, system_ids=("AA", "AG"))
        with _no_backend():
            prepare_campaign(root, "A6")

        (root / "AA" / "simforge" / "protein_only.pdb").unlink()
        result = prepare_campaign(root, "A6")
        assert not result.success
        assert "missing" in result.error.lower()

    def test_no_complexes_discovered_fails_clearly(self, tmp_path):
        (tmp_path / "empty").mkdir()
        result = prepare_campaign(tmp_path, "A6")
        assert not result.success
        assert "No protein-ligand complex" in result.error

    def test_ambiguous_directory_propagates_as_failure(self, tmp_path):
        _write_complex(tmp_path / "AA" / "one.pdb", "A6")
        _write_complex(tmp_path / "AA" / "two.pdb", "A6")
        result = prepare_campaign(tmp_path, "A6")
        assert not result.success
        assert "Ambiguous" in result.error
