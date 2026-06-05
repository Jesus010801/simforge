"""
Unit tests for ligand/preparation.py (Phase 2).

Skipped automatically when RDKit is not installed.
Run in the rdkit_env:
    conda run -n rdkit_env python -m pytest ligand/test_preparation.py -v
"""

from __future__ import annotations

from pathlib import Path

import pytest

# Skip this entire module if RDKit is not available
pytest.importorskip("rdkit.Chem", reason="RDKit not installed")

from rdkit import Chem
from rdkit.Chem import AllChem

from ligand.preparation import validate_ligand_for_parameterization


# ── Molecule factories ────────────────────────────────────────────────────────

def _sdf_with_h(smiles: str, tmp_path: Path, name: str = "MOL") -> Path:
    """Create an SDF file with 3D coords and explicit Hs for a SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    mol.SetProp("_Name", name)
    p = tmp_path / f"{name}.sdf"
    writer = Chem.SDWriter(str(p))
    writer.write(mol)
    writer.close()
    return p


def _sdf_no_h(smiles: str, tmp_path: Path, name: str = "MOL") -> Path:
    """Create an SDF with 3D coords but NO explicit Hs."""
    mol = Chem.MolFromSmiles(smiles)
    mol_h = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_h, AllChem.ETKDGv3())
    mol_no_h = Chem.RemoveAllHs(mol_h)  # heavy-atom 3D coords preserved
    mol_no_h.SetProp("_Name", name)
    p = tmp_path / f"{name}_noh.sdf"
    writer = Chem.SDWriter(str(p))
    writer.write(mol_no_h)
    writer.close()
    return p


def _sdf_2d(smiles: str, tmp_path: Path, name: str = "MOL") -> Path:
    """Create an SDF file with NO 3D conformer (only 2D layout)."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.Compute2DCoords(mol)
    mol.SetProp("_Name", name)
    p = tmp_path / f"{name}_2d.sdf"
    writer = Chem.SDWriter(str(p))
    writer.write(mol)
    writer.close()
    return p


# ── Happy path: valid molecule ────────────────────────────────────────────────

def test_valid_molecule_returns_valid(tmp_path):
    p = _sdf_with_h("c1ccccc1", tmp_path, "BENZ")
    result = validate_ligand_for_parameterization(p)
    assert result.valid is True


def test_valid_molecule_no_errors(tmp_path):
    p = _sdf_with_h("CCO", tmp_path, "ETH")
    result = validate_ligand_for_parameterization(p)
    assert result.errors == []


def test_valid_molecule_has_file_path(tmp_path):
    p = _sdf_with_h("CCO", tmp_path, "ETH")
    result = validate_ligand_for_parameterization(p)
    assert result.file_path == p


def test_valid_molecule_atom_count(tmp_path):
    # Benzene with H: 6C + 6H = 12 atoms
    p = _sdf_with_h("c1ccccc1", tmp_path, "BENZ")
    result = validate_ligand_for_parameterization(p)
    assert result.atom_count == 12


def test_valid_molecule_has_hydrogens(tmp_path):
    p = _sdf_with_h("CCO", tmp_path, "ETH")
    result = validate_ligand_for_parameterization(p)
    assert result.has_hydrogens is True


def test_molecule_name_derived(tmp_path):
    p = _sdf_with_h("CCO", tmp_path, "ETOH")
    result = validate_ligand_for_parameterization(p)
    assert result.molecule_name == "ETOH"


def test_molecule_name_truncated_to_4(tmp_path):
    p = _sdf_with_h("CCO", tmp_path, "LONGNAME")
    result = validate_ligand_for_parameterization(p)
    assert len(result.molecule_name) <= 4


# ── Warnings: no explicit H ───────────────────────────────────────────────────

def test_no_hydrogens_generates_warning(tmp_path):
    p = _sdf_no_h("c1ccccc1", tmp_path, "BENZ")
    result = validate_ligand_for_parameterization(p)
    assert result.has_hydrogens is False
    assert any("hydrogen" in w.lower() for w in result.warnings)


def test_no_hydrogens_still_valid(tmp_path):
    """Missing H is a warning, not a hard error."""
    p = _sdf_no_h("CCO", tmp_path, "ETH")
    result = validate_ligand_for_parameterization(p)
    assert result.valid is True


# ── Errors: no 3D conformer ───────────────────────────────────────────────────

def test_no_3d_conformer_invalid(tmp_path):
    p = _sdf_2d("CCO", tmp_path, "ETH")
    result = validate_ligand_for_parameterization(p)
    assert result.valid is False


def test_no_3d_conformer_has_error(tmp_path):
    p = _sdf_2d("c1ccccc1", tmp_path, "BENZ")
    result = validate_ligand_for_parameterization(p)
    assert any("3D" in e or "conformer" in e.lower() for e in result.errors)


# ── Errors: file not found ────────────────────────────────────────────────────

def test_missing_file_invalid(tmp_path):
    result = validate_ligand_for_parameterization(tmp_path / "ghost.sdf")
    assert result.valid is False


def test_missing_file_has_error(tmp_path):
    result = validate_ligand_for_parameterization(tmp_path / "ghost.sdf")
    assert len(result.errors) > 0


# ── Errors: corrupted file ────────────────────────────────────────────────────

def test_corrupt_sdf_invalid(tmp_path):
    p = tmp_path / "bad.sdf"
    p.write_text("this is not an sdf file\n$$$$\n")
    result = validate_ligand_for_parameterization(p)
    assert result.valid is False


# ── Warning: charged molecule ─────────────────────────────────────────────────

def test_charged_molecule_warning(tmp_path):
    # Ammonium cation NH4+
    mol = Chem.MolFromSmiles("[NH4+]")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    mol.SetProp("_Name", "NH4")
    p = tmp_path / "nh4.sdf"
    Chem.SDWriter(str(p)).write(mol)
    result = validate_ligand_for_parameterization(p)
    assert any("charge" in w.lower() or "+1" in w for w in result.warnings)


# ── Warning: unusual elements ─────────────────────────────────────────────────

def test_unusual_element_warning(tmp_path):
    # Trimethylsilanol: silicon
    mol = Chem.MolFromSmiles("C[Si](C)(C)O")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    mol.SetProp("_Name", "TMSOL")
    p = tmp_path / "si.sdf"
    Chem.SDWriter(str(p)).write(mol)
    result = validate_ligand_for_parameterization(p)
    assert any("Si" in w or "unusual" in w.lower() for w in result.warnings)


# ── Error: metal atom ─────────────────────────────────────────────────────────

def test_metal_atom_is_error(tmp_path):
    mol = Chem.MolFromSmiles("[Fe]")
    mol.SetProp("_Name", "FE")
    # No 3D needed to test the metal-blocking logic — but we need to give it a conformer
    # so the "no 3D" error doesn't shadow the metal error.
    mol3d = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol3d, AllChem.ETKDGv3())
    mol3d.SetProp("_Name", "FE")
    p = tmp_path / "fe.sdf"
    Chem.SDWriter(str(p)).write(mol3d)
    result = validate_ligand_for_parameterization(p)
    assert any("Fe" in e or "metal" in e.lower() for e in result.errors)
    assert result.valid is False


# ── Drug-like molecule end-to-end ─────────────────────────────────────────────

def test_ibuprofen_valid(tmp_path):
    # Ibuprofen SMILES
    p = _sdf_with_h("CC(C)Cc1ccc(cc1)C(C)C(=O)O", tmp_path, "IBU")
    result = validate_ligand_for_parameterization(p)
    assert result.valid is True
    assert result.atom_count > 0


def test_ibuprofen_molecule_name(tmp_path):
    p = _sdf_with_h("CC(C)Cc1ccc(cc1)C(C)C(=O)O", tmp_path, "IBU")
    result = validate_ligand_for_parameterization(p)
    assert result.molecule_name == "IBU"
