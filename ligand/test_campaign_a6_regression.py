"""
Campaign regression test using the real A6 ligand across four synthetic
receptor systems -- confirms the chemical-perception fix (element
normalization, geometry-based bond-order perception, reference-guided
mapping) composes correctly with the campaign architecture from the
previous phase, without requiring any change to that architecture beyond
consuming the new ``ligand.prepare.extract_ligand_from_complex`` fields
(already wired in ``ligand.campaign.prepare_campaign``).

Each of the four "receptors" gets the real A6 heavy-atom pose (reproducing
the reported docking defect: atom named "CL" with PDB element column "C"),
translated to a different location, so system-specific pose preservation is
verified against genuinely different coordinates per system.

RDKit >= 2022.09 (rdkit.Chem.rdDetermineBonds) is required; skipped
automatically otherwise -- see ligand/test_chemical_perception.py for the
same convention.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from ligand.campaign import prepare_campaign

FIXTURE_DIR = Path(__file__).parent.parent / "tests" / "fixtures" / "ligands" / "a6"
A6_REFERENCE = FIXTURE_DIR / "A6_reference_correct.pdb"
A6_HEAVY_ONLY_BUGGY = FIXTURE_DIR / "A6_docked_heavy_only_buggy_element.pdb"

_HAS_RDKIT_DETERMINE_BONDS = True
try:
    from rdkit.Chem import rdDetermineBonds  # noqa: F401
except ImportError:
    _HAS_RDKIT_DETERMINE_BONDS = False

requires_rdkit_determine_bonds = pytest.mark.skipif(
    not _HAS_RDKIT_DETERMINE_BONDS,
    reason="requires RDKit >= 2022.09 (rdkit.Chem.rdDetermineBonds)",
)

_SYSTEM_CENTERS = {
    "AA": (0.0, 0.0, 0.0),
    "AG": (30.0, 5.0, 2.0),
    "HMG-R": (-15.0, 20.0, -8.0),
    "LP": (10.0, -25.0, 12.0),
}


def _pdb_line(record, serial, atom_name, resname, chain, resseq, x, y, z, element):
    return (
        f"{record:<6}{serial:>5} {atom_name:>4} {resname:>3} {chain:1}{resseq:>4}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {element:>2}"
    )


def _protein_lines() -> list[str]:
    return [_pdb_line("ATOM", 1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N")]


def _a6_heavy_lines(center: tuple[float, float, float]) -> list[str]:
    """Real A6 heavy atoms (with the reported "CL"/element="C" defect),
    translated so each system has a distinct receptor-specific pose."""
    raw = [l for l in A6_HEAVY_ONLY_BUGGY.read_text().splitlines() if l.startswith("HETATM")]
    xs = [float(l[30:38]) for l in raw]
    ys = [float(l[38:46]) for l in raw]
    zs = [float(l[46:54]) for l in raw]
    cx, cy, cz = sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)
    dx, dy, dz = center[0] - cx, center[1] - cy, center[2] - cz

    out = []
    for l in raw:
        name = l[12:16].strip()
        element = l[76:78]
        x = float(l[30:38]) + dx
        y = float(l[38:46]) + dy
        z = float(l[46:54]) + dz
        out.append(_pdb_line("HETATM", int(l[6:11]), name, "A6", "A", 900, x, y, z, element.strip()))
    return out


def _write_complex(path: Path, center: tuple[float, float, float]) -> Path:
    lines = _protein_lines() + ["TER"] + _a6_heavy_lines(center) + ["END"]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


def _campaign_root(tmp_path) -> Path:
    for sid, center in _SYSTEM_CENTERS.items():
        _write_complex(tmp_path / sid / f"{sid}.pdb", center)
    return tmp_path


@requires_rdkit_determine_bonds
class TestA6CampaignRegression:
    def test_four_systems_deduplicate_to_one_ligand(self, tmp_path):
        root = _campaign_root(tmp_path)
        result = prepare_campaign(root, "A6")
        assert result.success, result.error
        assert result.systems_discovered == 4
        assert result.unique_ligands == 1
        assert set(result.manifest.ligands.keys()) == {"A6"}
        assert {s.id for s in result.manifest.systems} == set(_SYSTEM_CENTERS)

    def test_only_one_ligpargen_submission_written(self, tmp_path):
        root = _campaign_root(tmp_path)
        result = prepare_campaign(root, "A6")
        assert result.success, result.error
        param_files = list((root / "simforge_campaign" / "ligands").glob("*/parameterization/ligand_for_ligpargen*.pdb"))
        assert len(param_files) == 1

    def test_element_correction_recorded_once_not_per_system(self, tmp_path):
        """The Cl/C element defect is corrected during the single shared
        preparation pass, not once per receptor."""
        root = _campaign_root(tmp_path)
        result = prepare_campaign(root, "A6")
        assert result.success, result.error
        lig = result.manifest.ligands["A6"]
        assert len(lig.element_corrections) == 1
        assert lig.element_corrections[0]["atom_name"] == "CL"
        assert lig.element_corrections[0]["corrected_element"] == "Cl"

    def test_all_four_receptor_poses_preserved_unchanged(self, tmp_path):
        """Pose extraction must never be affected by chemical
        preparation/element correction -- every system's ligand_pose.pdb
        keeps its own original (uncorrected-element, receptor-specific)
        heavy-atom coordinates."""
        root = _campaign_root(tmp_path)
        result = prepare_campaign(root, "A6")
        assert result.success, result.error

        for sid, center in _SYSTEM_CENTERS.items():
            expected_lines = _a6_heavy_lines(center)
            expected_coords = [
                (round(float(l[30:38]), 3), round(float(l[38:46]), 3), round(float(l[46:54]), 3))
                for l in expected_lines
            ]
            pose_path = root / sid / "simforge" / "ligand_pose.pdb"
            actual_coords = [
                (round(float(l[30:38]), 3), round(float(l[38:46]), 3), round(float(l[46:54]), 3))
                for l in pose_path.read_text().splitlines()
                if l.startswith("HETATM")
            ]
            assert actual_coords == expected_coords, f"system {sid} pose was modified"

    def test_without_reference_heavy_only_pose_is_blocked(self, tmp_path):
        """No --ligand-reference and no explicit hydrogens in the docked
        pose: geometry alone cannot resolve bond orders (see
        ligand.chemical_perception). Must block, not guess."""
        root = _campaign_root(tmp_path)
        result = prepare_campaign(root, "A6")
        assert result.success, result.error
        lig = result.manifest.ligands["A6"]
        assert lig.parameterization_status == "chemical_validation_required"
        assert lig.readiness_block_reasons

    def test_with_reference_ligand_is_ready_with_correct_formula(self, tmp_path):
        """With --ligand-reference, the same heavy-atom-only, element-buggy
        campaign correctly reconstructs C18H15ClO3 and reports
        ready_for_ligpargen -- acceptance-criteria Option A applied to the
        full campaign workflow."""
        from collections import Counter

        root = _campaign_root(tmp_path)
        result = prepare_campaign(root, "A6", ligand_reference=A6_REFERENCE)
        assert result.success, result.error

        lig = result.manifest.ligands["A6"]
        assert lig.parameterization_status == "ready_for_ligpargen"
        assert lig.chemical_identity_method == "chemical_reference"
        assert not lig.readiness_block_reasons

        param_pdb = root / lig.ligand_for_ligpargen
        formula = Counter(
            l[76:78].strip().upper() for l in param_pdb.read_text().splitlines()
            if l.startswith(("ATOM", "HETATM"))
        )
        assert dict(formula) == {"C": 18, "O": 3, "CL": 1, "H": 15}

    def test_resname_and_chain_consistent_after_reference_guided_prep(self, tmp_path):
        root = _campaign_root(tmp_path)
        result = prepare_campaign(root, "A6", ligand_reference=A6_REFERENCE)
        assert result.success, result.error
        lig = result.manifest.ligands["A6"]
        param_pdb = root / lig.ligand_for_ligpargen
        lines = [l for l in param_pdb.read_text().splitlines() if l.startswith(("ATOM", "HETATM"))]
        assert {l[17:20].strip() for l in lines} == {"A6"}
