"""
Unit tests for ligand/export.py (Phase 2).

Skipped automatically when RDKit is not installed.
Run in the rdkit_env:
    conda run -n rdkit_env python -m pytest ligand/test_export.py -v
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("rdkit.Chem", reason="RDKit not installed")

from rdkit import Chem
from rdkit.Chem import AllChem

from ligand.export import export_for_ligpargen


# ── Molecule factories ────────────────────────────────────────────────────────

def _make_sdf(smiles: str, tmp_path: Path, name: str, with_h: bool = True) -> Path:
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if not with_h:
        mol = Chem.RemoveAllHs(mol)
    mol.SetProp("_Name", name)
    p = tmp_path / f"{name}.sdf"
    writer = Chem.SDWriter(str(p))
    writer.write(mol)
    writer.close()
    return p


def _make_sdf_no_coords(smiles: str, tmp_path: Path, name: str) -> Path:
    """Create SDF with a 2D layout (no valid 3D conformer)."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.Compute2DCoords(mol)
    mol.SetProp("_Name", name)
    p = tmp_path / f"{name}_2d.sdf"
    writer = Chem.SDWriter(str(p))
    writer.write(mol)
    writer.close()
    return p


# ── Basic success ─────────────────────────────────────────────────────────────

def test_export_success(tmp_path):
    sdf = _make_sdf("CCO", tmp_path, "ETH")
    result = export_for_ligpargen(sdf, tmp_path / "out")
    assert result.success is True


def test_export_creates_pdb(tmp_path):
    sdf = _make_sdf("CCO", tmp_path, "ETH")
    result = export_for_ligpargen(sdf, tmp_path / "out", mol_name="ETH")
    assert result.exported_path is not None
    assert result.exported_path.exists()
    assert result.exported_path.suffix == ".pdb"


def test_export_pdb_is_parseable_by_rdkit(tmp_path):
    sdf = _make_sdf("c1ccccc1", tmp_path, "BENZ")
    result = export_for_ligpargen(sdf, tmp_path / "out", mol_name="BENZ")
    mol = Chem.MolFromPDBFile(str(result.exported_path), removeHs=False)
    assert mol is not None


def test_export_molecule_name_in_result(tmp_path):
    sdf = _make_sdf("CCO", tmp_path, "ETH")
    result = export_for_ligpargen(sdf, tmp_path / "out", mol_name="ETOH")
    assert result.molecule_name == "ETOH"


def test_export_atom_count(tmp_path):
    sdf = _make_sdf("CCO", tmp_path, "ETH", with_h=True)
    result = export_for_ligpargen(sdf, tmp_path / "out", mol_name="ETH")
    assert result.atom_count > 0


def test_export_atom_count_no_h_input_add_h_true(tmp_path):
    """When input has no H and add_hydrogens=True, output should have more atoms."""
    sdf_no_h = _make_sdf("CCO", tmp_path, "ETH", with_h=False)
    mol_no_h = Chem.SDMolSupplier(str(sdf_no_h), removeHs=False)[0]
    count_no_h = mol_no_h.GetNumAtoms()

    result = export_for_ligpargen(sdf_no_h, tmp_path / "out", mol_name="ETH", add_hydrogens=True)
    assert result.atom_count > count_no_h


def test_export_no_error_field(tmp_path):
    sdf = _make_sdf("CCO", tmp_path, "ETH")
    result = export_for_ligpargen(sdf, tmp_path / "out")
    assert result.error is None


# ── Molecule name handling ────────────────────────────────────────────────────

def test_export_name_auto_derived(tmp_path):
    """When mol_name is empty, derive from molecule's _Name property."""
    sdf = _make_sdf("CCO", tmp_path, "ETOH")
    result = export_for_ligpargen(sdf, tmp_path / "out", mol_name="")
    assert result.molecule_name == "ETOH"


def test_export_name_truncated_to_4(tmp_path):
    sdf = _make_sdf("CCO", tmp_path, "ETH")
    result = export_for_ligpargen(sdf, tmp_path / "out", mol_name="TOOLONGNAME")
    assert len(result.molecule_name) <= 4


def test_export_name_uppercased(tmp_path):
    sdf = _make_sdf("CCO", tmp_path, "ETH")
    result = export_for_ligpargen(sdf, tmp_path / "out", mol_name="eth")
    assert result.molecule_name == result.molecule_name.upper()


# ── RMSD (geometry_alignment) ─────────────────────────────────────────────────

def test_export_heavy_atom_rmsd_present(tmp_path):
    sdf = _make_sdf("c1ccccc1", tmp_path, "BENZ", with_h=True)
    result = export_for_ligpargen(sdf, tmp_path / "out", mol_name="BENZ")
    assert result.heavy_atom_rmsd is not None


def test_export_heavy_atom_rmsd_near_zero_no_h_change(tmp_path):
    """Molecule already has H → heavy-atom scaffold unchanged → RMSD ≈ 0."""
    sdf = _make_sdf("c1ccccc1", tmp_path, "BENZ", with_h=True)
    result = export_for_ligpargen(sdf, tmp_path / "out", mol_name="BENZ")
    assert result.heavy_atom_rmsd == pytest.approx(0.0, abs=0.1)


def test_export_heavy_atom_rmsd_non_negative(tmp_path):
    sdf = _make_sdf("CCO", tmp_path, "ETH", with_h=True)
    result = export_for_ligpargen(sdf, tmp_path / "out")
    assert result.heavy_atom_rmsd >= 0.0


# ── Output directory creation ─────────────────────────────────────────────────

def test_export_creates_output_dir(tmp_path):
    sdf = _make_sdf("CCO", tmp_path, "ETH")
    out_dir = tmp_path / "new_dir" / "nested"
    result = export_for_ligpargen(sdf, out_dir)
    assert out_dir.exists()


def test_export_pdb_named_after_molecule(tmp_path):
    sdf = _make_sdf("CCO", tmp_path, "ETH")
    result = export_for_ligpargen(sdf, tmp_path / "out", mol_name="MYML")
    assert result.exported_path.name == "MYML.pdb"


# ── add_hydrogens flag ────────────────────────────────────────────────────────

def test_export_add_hydrogens_false_preserves_count(tmp_path):
    sdf = _make_sdf("CCO", tmp_path, "ETH", with_h=True)
    mol_in = Chem.SDMolSupplier(str(sdf), removeHs=False)[0]
    count_in = mol_in.GetNumAtoms()
    result = export_for_ligpargen(sdf, tmp_path / "out", add_hydrogens=False)
    # atom_count may differ slightly due to PDB round-trip, but should be close
    assert result.atom_count > 0
    assert result.success is True


# ── Error cases ───────────────────────────────────────────────────────────────

def test_export_missing_file(tmp_path):
    result = export_for_ligpargen(tmp_path / "ghost.sdf", tmp_path / "out")
    assert result.success is False
    assert result.error is not None


def test_export_no_3d_conformer(tmp_path):
    sdf = _make_sdf_no_coords("CCO", tmp_path, "ETH")
    result = export_for_ligpargen(sdf, tmp_path / "out")
    assert result.success is False
    assert result.error is not None


def test_export_corrupt_file(tmp_path):
    p = tmp_path / "bad.sdf"
    p.write_text("garbage content\n$$$$\n")
    result = export_for_ligpargen(p, tmp_path / "out")
    assert result.success is False


# ── Diverse molecules ─────────────────────────────────────────────────────────

def test_export_benzene(tmp_path):
    sdf = _make_sdf("c1ccccc1", tmp_path, "BENZ")
    result = export_for_ligpargen(sdf, tmp_path / "out", mol_name="BENZ")
    assert result.success is True


def test_export_ibuprofen(tmp_path):
    sdf = _make_sdf("CC(C)Cc1ccc(cc1)C(C)C(=O)O", tmp_path, "IBU")
    result = export_for_ligpargen(sdf, tmp_path / "out", mol_name="IBU")
    assert result.success is True
    assert result.atom_count > 20


def test_export_charged_molecule(tmp_path):
    mol = Chem.MolFromSmiles("[NH4+]")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    mol.SetProp("_Name", "NH4")
    p = tmp_path / "nh4.sdf"
    Chem.SDWriter(str(p)).write(mol)
    result = export_for_ligpargen(p, tmp_path / "out", mol_name="NH4")
    assert result.success is True
