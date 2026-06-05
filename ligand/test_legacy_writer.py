"""
Regression tests for LigParGenLegacyWriter and export_for_ligpargen_legacy.

All tests are skipped automatically when RDKit is not installed.
Run in rdkit_env:
    conda run -n rdkit_env python -m pytest ligand/test_legacy_writer.py -v

Design notes
------------
The "accepted format" from the real LigParGen validation has these properties:
  - ATOM records (not HETATM)
  - Simple element-only atom names: " C  ", " O  ", " H  ", " N  "
  - No trailing digits in atom names
  - CONECT records for every bond (including H–heavy bonds)
  - Coordinates preserved to 3 decimal places (8.3f)
  - Heavy-atom RMSD ≈ 0 (scaffold unchanged)

The golden-format tests pin the exact byte layout for a hand-built molecule
with known coordinates; any change to the writer breaks them immediately.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("rdkit.Chem", reason="RDKit not installed")

from rdkit import Chem
from rdkit.Chem import AllChem

from ligand.legacy_writer import LigParGenLegacyWriter, _pdb_atom_name, write_legacy_pdb
from ligand.export import export_for_ligpargen_legacy


# ── Molecule helpers ──────────────────────────────────────────────────────────

def _mol_with_known_coords():
    """
    Return a methanol molecule (CH3OH, SMILES "CO") with manually placed coordinates.

    Atom order after AddHs: C(0) O(1) H(2) H(3) H(4) H(5)
    Coordinates are chosen so the expected PDB output can be written by hand.
    """
    from rdkit.Chem import RWMol, Conformer

    # AddHs returns a Mol; wrap in RWMol to call AddConformer
    mol_h = Chem.AddHs(Chem.MolFromSmiles("CO"))
    rw = RWMol(mol_h)
    n = rw.GetNumAtoms()
    conf = Conformer(n)
    positions = [
        (0.000,  0.000,  0.000),   # C
        (1.000,  0.000,  0.000),   # O
        (-0.500, 0.500,  0.500),   # H on C
        (-0.500, 0.500, -0.500),   # H on C
        (-0.500,-0.500,  0.000),   # H on C
        ( 1.500, 0.800,  0.000),   # H on O
    ]
    for i, (x, y, z) in enumerate(positions):
        conf.SetAtomPosition(i, (x, y, z))
    rw.AddConformer(conf, assignId=True)
    return rw.GetMol()


def _make_sdf(smiles: str, tmp_path: Path, name: str, with_h: bool = True) -> Path:
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if not with_h:
        mol = Chem.RemoveAllHs(mol)
    mol.SetProp("_Name", name)
    p = tmp_path / f"{name}.sdf"
    Chem.SDWriter(str(p)).write(mol)
    return p


# ══════════════════════════════════════════════════════════════════════════════
# _pdb_atom_name
# ══════════════════════════════════════════════════════════════════════════════

class TestPdbAtomName:
    def test_carbon(self):
        assert _pdb_atom_name("C") == " C  "

    def test_oxygen(self):
        assert _pdb_atom_name("O") == " O  "

    def test_hydrogen(self):
        assert _pdb_atom_name("H") == " H  "

    def test_nitrogen(self):
        assert _pdb_atom_name("N") == " N  "

    def test_sulfur(self):
        assert _pdb_atom_name("S") == " S  "

    def test_phosphorus(self):
        assert _pdb_atom_name("P") == " P  "

    def test_chlorine(self):
        assert _pdb_atom_name("Cl") == "CL  "

    def test_bromine(self):
        assert _pdb_atom_name("Br") == "BR  "

    def test_fluorine(self):
        assert _pdb_atom_name("F") == " F  "

    def test_iodine(self):
        assert _pdb_atom_name("I") == " I  "

    def test_always_4_chars(self):
        for element in ["C", "O", "H", "N", "S", "Cl", "Br", "F", "I", "P"]:
            assert len(_pdb_atom_name(element)) == 4, f"Failed for {element}"

    def test_uppercase(self):
        assert _pdb_atom_name("c") == " C  "
        assert _pdb_atom_name("cl") == "CL  "

    def test_no_digits_in_name(self):
        name = _pdb_atom_name("C")
        assert not any(ch.isdigit() for ch in name)


# ══════════════════════════════════════════════════════════════════════════════
# LigParGenLegacyWriter — format regression (golden tests)
# ══════════════════════════════════════════════════════════════════════════════

class TestLegacyWriterFormat:
    """
    Golden tests that pin the exact PDB output for a known molecule.
    These detect any accidental change to the column layout.
    """

    def test_atom_record_type(self, tmp_path):
        mol = _mol_with_known_coords()
        out = tmp_path / "mol.pdb"
        LigParGenLegacyWriter().write(mol, out, "LIG")
        lines = out.read_text().splitlines()
        atom_lines = [l for l in lines if l[:4] == "ATOM"]
        assert len(atom_lines) > 0
        for line in atom_lines:
            assert line.startswith("ATOM  "), f"Expected ATOM record, got: {line!r}"

    def test_no_hetatm_records(self, tmp_path):
        mol = _mol_with_known_coords()
        out = tmp_path / "mol.pdb"
        LigParGenLegacyWriter().write(mol, out, "LIG")
        lines = out.read_text().splitlines()
        hetatm = [l for l in lines if l.startswith("HETATM")]
        assert hetatm == [], f"Unexpected HETATM records: {hetatm}"

    def test_golden_first_atom_line(self, tmp_path):
        """Pin the exact byte layout of the first ATOM record (C at origin)."""
        mol = _mol_with_known_coords()
        out = tmp_path / "mol.pdb"
        LigParGenLegacyWriter().write(mol, out, "LIG")
        lines = out.read_text().splitlines()
        first_atom = next(l for l in lines if l.startswith("ATOM"))
        expected = (
            "ATOM      1  C   LIG     1       0.000   0.000   0.000  1.00  0.00"
        )
        assert first_atom == expected, (
            f"Golden mismatch:\n  got:      {first_atom!r}\n  expected: {expected!r}"
        )

    def test_golden_second_atom_line(self, tmp_path):
        """Pin the oxygen atom record."""
        mol = _mol_with_known_coords()
        out = tmp_path / "mol.pdb"
        LigParGenLegacyWriter().write(mol, out, "LIG")
        lines = out.read_text().splitlines()
        atom_lines = [l for l in lines if l.startswith("ATOM")]
        expected = (
            "ATOM      2  O   LIG     1       1.000   0.000   0.000  1.00  0.00"
        )
        assert atom_lines[1] == expected

    def test_golden_hydrogen_atom_line(self, tmp_path):
        """Pin a hydrogen atom record."""
        mol = _mol_with_known_coords()
        out = tmp_path / "mol.pdb"
        LigParGenLegacyWriter().write(mol, out, "LIG")
        lines = out.read_text().splitlines()
        atom_lines = [l for l in lines if l.startswith("ATOM")]
        # Atom 3 (index 2) is the first H, coords (-0.500, 0.500, 0.500)
        expected = (
            "ATOM      3  H   LIG     1      -0.500   0.500   0.500  1.00  0.00"
        )
        assert atom_lines[2] == expected

    def test_atom_name_col_positions(self, tmp_path):
        """Atom name occupies cols 12–15 (0-indexed); must be 4 chars."""
        mol = _mol_with_known_coords()
        out = tmp_path / "mol.pdb"
        LigParGenLegacyWriter().write(mol, out, "LIG")
        lines = out.read_text().splitlines()
        for line in [l for l in lines if l.startswith("ATOM")]:
            name_field = line[12:16]
            assert len(name_field) == 4
            assert not any(ch.isdigit() for ch in name_field), (
                f"Digit in atom name field: {name_field!r} in {line!r}"
            )

    def test_residue_name_col_positions(self, tmp_path):
        """Residue name occupies cols 17–19 (0-indexed)."""
        mol = _mol_with_known_coords()
        out = tmp_path / "mol.pdb"
        LigParGenLegacyWriter().write(mol, out, "MOL")
        lines = out.read_text().splitlines()
        for line in [l for l in lines if l.startswith("ATOM")]:
            resname_field = line[17:20]
            assert resname_field == "MOL", f"Residue name wrong: {resname_field!r}"

    def test_x_col_positions(self, tmp_path):
        """x coordinate occupies cols 30–37 (0-indexed)."""
        mol = _mol_with_known_coords()
        out = tmp_path / "mol.pdb"
        LigParGenLegacyWriter().write(mol, out, "LIG")
        lines = out.read_text().splitlines()
        first_atom = next(l for l in lines if l.startswith("ATOM"))
        x_field = first_atom[30:38]
        assert float(x_field.strip()) == pytest.approx(0.000, abs=1e-3)

    def test_line_length(self, tmp_path):
        """Each ATOM record must be exactly 66 chars (no element column appended)."""
        mol = _mol_with_known_coords()
        out = tmp_path / "mol.pdb"
        LigParGenLegacyWriter().write(mol, out, "LIG")
        for line in [l for l in out.read_text().splitlines() if l.startswith("ATOM")]:
            assert len(line) == 66, f"Line length {len(line)} != 66: {line!r}"

    def test_ends_with_end_record(self, tmp_path):
        mol = _mol_with_known_coords()
        out = tmp_path / "mol.pdb"
        LigParGenLegacyWriter().write(mol, out, "LIG")
        lines = out.read_text().strip().splitlines()
        assert lines[-1] == "END"


# ══════════════════════════════════════════════════════════════════════════════
# LigParGenLegacyWriter — CONECT records
# ══════════════════════════════════════════════════════════════════════════════

class TestConectRecords:
    def test_conect_records_present(self, tmp_path):
        mol = _mol_with_known_coords()
        out = tmp_path / "mol.pdb"
        LigParGenLegacyWriter().write(mol, out, "LIG")
        lines = out.read_text().splitlines()
        conect_lines = [l for l in lines if l.startswith("CONECT")]
        assert len(conect_lines) > 0

    def test_conect_covers_all_atoms_with_bonds(self, tmp_path):
        mol = _mol_with_known_coords()
        out = tmp_path / "mol.pdb"
        LigParGenLegacyWriter().write(mol, out, "LIG")
        lines = out.read_text().splitlines()
        conect_lines = [l for l in lines if l.startswith("CONECT")]
        # All 6 atoms have bonds → 6 CONECT lines (one per atom)
        serials_with_conect = {int(l[6:11]) for l in conect_lines}
        expected = set(range(1, mol.GetNumAtoms() + 1))
        assert serials_with_conect == expected

    def test_conect_format(self, tmp_path):
        """Each CONECT line: 'CONECT' (6) + serials (5 chars each)."""
        mol = _mol_with_known_coords()
        out = tmp_path / "mol.pdb"
        LigParGenLegacyWriter().write(mol, out, "LIG")
        for line in [l for l in out.read_text().splitlines() if l.startswith("CONECT")]:
            assert line.startswith("CONECT")
            remainder = line[6:]
            assert len(remainder) % 5 == 0, f"CONECT tail not multiple of 5: {line!r}"
            for i in range(0, len(remainder), 5):
                assert remainder[i:i+5].strip().isdigit(), (
                    f"Non-integer in CONECT: {remainder[i:i+5]!r}"
                )

    def test_conect_carbon_bonded_to_oxygen_and_hydrogens(self, tmp_path):
        """Carbon (serial 1) must list O (serial 2) and its 3 H atoms."""
        mol = _mol_with_known_coords()
        out = tmp_path / "mol.pdb"
        LigParGenLegacyWriter().write(mol, out, "LIG")
        lines = out.read_text().splitlines()
        carbon_conect = next(
            l for l in lines if l.startswith("CONECT") and int(l[6:11]) == 1
        )
        # Cols 6-10 are the source atom serial; bonded atoms start at col 11
        bonded = [int(carbon_conect[11 + 5*i : 11 + 5*(i+1)])
                  for i in range((len(carbon_conect) - 11) // 5)]
        # C bonded to O and 3 H atoms → 4 neighbours
        assert 2 in bonded   # O is atom 2
        assert len(bonded) == 4


# ══════════════════════════════════════════════════════════════════════════════
# LigParGenLegacyWriter — coordinate fidelity
# ══════════════════════════════════════════════════════════════════════════════

class TestCoordinateFidelity:
    def test_coordinates_preserved_to_3dp(self, tmp_path):
        mol = _mol_with_known_coords()
        out = tmp_path / "mol.pdb"
        LigParGenLegacyWriter().write(mol, out, "LIG")
        conf = mol.GetConformer()
        lines = out.read_text().splitlines()
        atom_lines = [l for l in lines if l.startswith("ATOM")]
        for i, line in enumerate(atom_lines):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            pos = conf.GetAtomPosition(i)
            assert x == pytest.approx(pos.x, abs=5e-4)
            assert y == pytest.approx(pos.y, abs=5e-4)
            assert z == pytest.approx(pos.z, abs=5e-4)

    def test_negative_coordinate_format(self, tmp_path):
        """Negative coordinates must fit in 8.3f without overflowing the field."""
        mol = _mol_with_known_coords()
        out = tmp_path / "mol.pdb"
        LigParGenLegacyWriter().write(mol, out, "LIG")
        lines = [l for l in out.read_text().splitlines() if l.startswith("ATOM")]
        # Atom 3 (index 2) has x=-0.500
        x_field = lines[2][30:38]
        assert float(x_field) == pytest.approx(-0.500, abs=5e-4)
        assert len(x_field) == 8

    def test_atom_count_matches_mol(self, tmp_path):
        mol = _mol_with_known_coords()
        out = tmp_path / "mol.pdb"
        LigParGenLegacyWriter().write(mol, out, "LIG")
        atom_lines = [l for l in out.read_text().splitlines() if l.startswith("ATOM")]
        assert len(atom_lines) == mol.GetNumAtoms()


# ══════════════════════════════════════════════════════════════════════════════
# LigParGenLegacyWriter — molecule name handling
# ══════════════════════════════════════════════════════════════════════════════

class TestMoleculeNameHandling:
    def test_resname_truncated_to_3(self, tmp_path):
        mol = _mol_with_known_coords()
        out = tmp_path / "mol.pdb"
        LigParGenLegacyWriter().write(mol, out, "LONGNAME")
        lines = [l for l in out.read_text().splitlines() if l.startswith("ATOM")]
        resname = lines[0][17:20]
        assert resname == "LON"

    def test_resname_uppercased(self, tmp_path):
        mol = _mol_with_known_coords()
        out = tmp_path / "mol.pdb"
        LigParGenLegacyWriter().write(mol, out, "eth")
        lines = [l for l in out.read_text().splitlines() if l.startswith("ATOM")]
        resname = lines[0][17:20]
        assert resname == "ETH"

    def test_resname_padded_to_3(self, tmp_path):
        """Short names (≤2 chars) are right-padded with spaces."""
        mol = _mol_with_known_coords()
        out = tmp_path / "mol.pdb"
        LigParGenLegacyWriter().write(mol, out, "C")
        lines = [l for l in out.read_text().splitlines() if l.startswith("ATOM")]
        resname = lines[0][17:20]
        assert len(resname) == 3

    def test_default_resname_is_LIG(self, tmp_path):
        mol = _mol_with_known_coords()
        out = tmp_path / "mol.pdb"
        LigParGenLegacyWriter().write(mol, out)  # no mol_name
        lines = [l for l in out.read_text().splitlines() if l.startswith("ATOM")]
        assert lines[0][17:20] == "LIG"


# ══════════════════════════════════════════════════════════════════════════════
# write_legacy_pdb convenience alias
# ══════════════════════════════════════════════════════════════════════════════

def test_write_legacy_pdb_alias(tmp_path):
    mol = _mol_with_known_coords()
    out = tmp_path / "alias.pdb"
    write_legacy_pdb(mol, out, "LIG")
    content = out.read_text()
    assert "ATOM" in content
    assert "CONECT" in content
    assert "END" in content


# ══════════════════════════════════════════════════════════════════════════════
# export_for_ligpargen_legacy — integration tests
# ══════════════════════════════════════════════════════════════════════════════

class TestExportLegacy:
    def test_success(self, tmp_path):
        sdf = _make_sdf("CCO", tmp_path, "ETH")
        result = export_for_ligpargen_legacy(sdf, tmp_path / "out", mol_name="ETH")
        assert result.success is True

    def test_filename_contains_legacy(self, tmp_path):
        sdf = _make_sdf("CCO", tmp_path, "ETH")
        result = export_for_ligpargen_legacy(sdf, tmp_path / "out", mol_name="ETH")
        assert "legacy" in result.exported_path.name.lower()

    def test_output_is_pdb(self, tmp_path):
        sdf = _make_sdf("c1ccccc1", tmp_path, "BENZ")
        result = export_for_ligpargen_legacy(sdf, tmp_path / "out", mol_name="BENZ")
        assert result.exported_path.suffix == ".pdb"

    def test_output_uses_atom_records(self, tmp_path):
        sdf = _make_sdf("CCO", tmp_path, "ETH")
        result = export_for_ligpargen_legacy(sdf, tmp_path / "out", mol_name="ETH")
        content = result.exported_path.read_text()
        assert "ATOM  " in content
        assert "HETATM" not in content

    def test_output_has_no_numbered_atom_names(self, tmp_path):
        """Atom name field (cols 12-15) must not contain digits."""
        sdf = _make_sdf("CCO", tmp_path, "ETH")
        result = export_for_ligpargen_legacy(sdf, tmp_path / "out", mol_name="ETH")
        for line in result.exported_path.read_text().splitlines():
            if line.startswith("ATOM"):
                name_field = line[12:16]
                assert not any(ch.isdigit() for ch in name_field), (
                    f"Digit in atom name: {name_field!r}"
                )

    def test_output_has_conect_records(self, tmp_path):
        sdf = _make_sdf("CCO", tmp_path, "ETH")
        result = export_for_ligpargen_legacy(sdf, tmp_path / "out", mol_name="ETH")
        content = result.exported_path.read_text()
        assert "CONECT" in content

    def test_heavy_atom_rmsd_near_zero(self, tmp_path):
        """Heavy-atom scaffold must be preserved: RMSD should be ~0."""
        sdf = _make_sdf("c1ccccc1", tmp_path, "BENZ")
        result = export_for_ligpargen_legacy(sdf, tmp_path / "out", mol_name="BENZ")
        assert result.heavy_atom_rmsd is not None
        assert result.heavy_atom_rmsd == pytest.approx(0.0, abs=0.1)

    def test_missing_file_error(self, tmp_path):
        result = export_for_ligpargen_legacy(tmp_path / "ghost.sdf", tmp_path / "out")
        assert result.success is False
        assert result.error is not None

    def test_no_3d_conformer_error(self, tmp_path):
        mol = Chem.MolFromSmiles("CCO")
        mol = Chem.AddHs(mol)
        AllChem.Compute2DCoords(mol)
        mol.SetProp("_Name", "ETH")
        sdf = tmp_path / "eth_2d.sdf"
        Chem.SDWriter(str(sdf)).write(mol)
        result = export_for_ligpargen_legacy(sdf, tmp_path / "out")
        assert result.success is False

    def test_atom_count_in_result(self, tmp_path):
        sdf = _make_sdf("CCO", tmp_path, "ETH")
        result = export_for_ligpargen_legacy(sdf, tmp_path / "out", mol_name="ETH")
        assert result.atom_count > 0

    def test_ibuprofen_accepted_format(self, tmp_path):
        """Drug-like molecule: verify the accepted format end-to-end."""
        sdf = _make_sdf("CC(C)Cc1ccc(cc1)C(C)C(=O)O", tmp_path, "IBU")
        result = export_for_ligpargen_legacy(sdf, tmp_path / "out", mol_name="IBU")
        assert result.success is True
        lines = result.exported_path.read_text().splitlines()
        atom_lines = [l for l in lines if l.startswith("ATOM")]
        conect_lines = [l for l in lines if l.startswith("CONECT")]
        # Ibuprofen has many atoms and bonds
        assert len(atom_lines) > 20
        assert len(conect_lines) > 20
        # All ATOM records use element-only names
        for line in atom_lines:
            name_field = line[12:16]
            assert not any(ch.isdigit() for ch in name_field)

    def test_creates_output_directory(self, tmp_path):
        sdf = _make_sdf("CCO", tmp_path, "ETH")
        out_dir = tmp_path / "new" / "nested"
        result = export_for_ligpargen_legacy(sdf, out_dir, mol_name="ETH")
        assert out_dir.exists()

    def test_rdkit_can_parse_output(self, tmp_path):
        """The output PDB must be parseable by RDKit (round-trip sanity)."""
        sdf = _make_sdf("CCO", tmp_path, "ETH")
        result = export_for_ligpargen_legacy(sdf, tmp_path / "out", mol_name="ETH")
        mol_back = Chem.MolFromPDBFile(str(result.exported_path), removeHs=False)
        assert mol_back is not None
