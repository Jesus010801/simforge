"""
Tests for ligand/ligpargen_import_validator.py.

Uses the real LigParGen fixture from tests/fixtures/ligpargen/a1/ for the
happy path; synthetic minimal files for error/edge-case paths.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from ligand.ligpargen_import_validator import LigParGenImportValidator
from ligand.normalization import LigandIdentity
from utils.gro_parser import GroAtom, GroFile, parse_gro, write_gro
from utils.itp_parser import parse_itp

FIXTURE_DIR = Path(__file__).parent.parent / "tests" / "fixtures" / "ligpargen" / "a1"
FIXTURE_GRO = FIXTURE_DIR / "A1_ligpargen.gro"
FIXTURE_ITP = FIXTURE_DIR / "A1_ligpargen.itp"

_FIXTURE_ATOM_COUNT = 35
_FIXTURE_SOURCE_MOLECULETYPE = "H"


# ── Shared helpers ─────────────────────────────────────────────────────────────

def _make_identity(internal_id: str = "L01") -> LigandIdentity:
    return LigandIdentity(
        component_id="a1",
        display_name="A1",
        source_filename="A1_ligpargen.itp",
        internal_id=internal_id,
        residue_name=internal_id,
        moleculetype=internal_id,
    )


def _minimal_gro(path: Path, n_atoms: int = 1, residue_name: str = "LIG") -> Path:
    """Write a syntactically valid GRO with n_atoms atoms."""
    atoms = [
        GroAtom(
            residue_number=1,
            residue_name=residue_name,
            atom_name=f"C{i + 1}",
            atom_number=i + 1,
            x=1.0, y=1.0, z=1.0,
        )
        for i in range(n_atoms)
    ]
    write_gro(GroFile(title="Minimal", atoms=atoms, box=[1.0, 1.0, 1.0]), path)
    return path


def _minimal_itp(
    path: Path,
    n_atoms: int = 1,
    residue_name: str = "LIG",
    moleculetype_name: str = "LIG",
    charge_per_atom: float = 0.0,
    include_moleculetype: bool = True,
    include_atoms: bool = True,
) -> Path:
    """Write a minimal but parseable ITP."""
    lines: list[str] = []
    if include_moleculetype:
        lines += [
            "[ moleculetype ]",
            f"; Name  nrexcl",
            f"{moleculetype_name}   3",
        ]
    if include_atoms:
        lines += [
            "[ atoms ]",
            ";   nr  type  resnr  residue  atom  cgnr  charge  mass",
        ]
        for i in range(n_atoms):
            lines.append(
                f"     {i + 1}   CT    1    {residue_name}     C{i + 1}    1"
                f"   {charge_per_atom:.4f}  12.0110"
            )
    path.write_text("\n".join(lines) + "\n")
    return path


# ── Reusable result fixture for the happy path ────────────────────────────────

@pytest.fixture(scope="module")
def happy_result(tmp_path_factory):
    work = tmp_path_factory.mktemp("validator_happy")
    v = LigParGenImportValidator()
    return v.validate(FIXTURE_GRO, FIXTURE_ITP, _make_identity(), work)


# ── Happy path ─────────────────────────────────────────────────────────────────

def test_valid_fixture_is_valid(happy_result):
    assert happy_result.valid is True


def test_valid_fixture_has_no_errors(happy_result):
    assert happy_result.errors == []


# ── File existence checks ──────────────────────────────────────────────────────

def test_missing_gro_fails(tmp_path):
    v = LigParGenImportValidator()
    r = v.validate(tmp_path / "ghost.gro", FIXTURE_ITP, _make_identity(), tmp_path / "out")
    assert r.valid is False
    assert any("GRO file not found" in e for e in r.errors)


def test_missing_itp_fails(tmp_path):
    v = LigParGenImportValidator()
    r = v.validate(FIXTURE_GRO, tmp_path / "ghost.itp", _make_identity(), tmp_path / "out")
    assert r.valid is False
    assert any("ITP file not found" in e for e in r.errors)


# ── Required ITP sections ─────────────────────────────────────────────────────

def test_missing_moleculetype_section_fails(tmp_path):
    gro = _minimal_gro(tmp_path / "a.gro")
    itp = _minimal_itp(tmp_path / "a.itp", include_moleculetype=False)
    v = LigParGenImportValidator()
    r = v.validate(gro, itp, _make_identity(), tmp_path / "out")
    assert r.valid is False
    assert any("moleculetype" in e for e in r.errors)


def test_missing_atoms_section_fails(tmp_path):
    gro = _minimal_gro(tmp_path / "a.gro")
    itp = _minimal_itp(tmp_path / "a.itp", include_atoms=False)
    v = LigParGenImportValidator()
    r = v.validate(gro, itp, _make_identity(), tmp_path / "out")
    assert r.valid is False
    assert any("atoms" in e for e in r.errors)


# ── Atom count consistency ────────────────────────────────────────────────────

def test_atom_count_mismatch_fails(tmp_path):
    # GRO has 1 atom, real ITP has 35 atoms
    gro = _minimal_gro(tmp_path / "short.gro", n_atoms=1)
    v = LigParGenImportValidator()
    r = v.validate(gro, FIXTURE_ITP, _make_identity(), tmp_path / "out")
    assert r.valid is False
    assert any("mismatch" in e.lower() for e in r.errors)


# ── Total charge ──────────────────────────────────────────────────────────────

def test_non_integer_charge_produces_warning(tmp_path):
    # charge_per_atom=0.3 → total = 0.3, not near-integer
    gro = _minimal_gro(tmp_path / "a.gro")
    itp = _minimal_itp(tmp_path / "a.itp", charge_per_atom=0.3)
    v = LigParGenImportValidator()
    r = v.validate(gro, itp, _make_identity(), tmp_path / "out")
    assert any("Non-integer" in w or "non-integer" in w.lower() for w in r.warnings)
    assert r.charge_integer is False


def test_integer_charge_produces_no_charge_warning(happy_result):
    assert not any("charge" in w.lower() for w in happy_result.warnings)
    assert happy_result.charge_integer is True


# ── Generic name detection ────────────────────────────────────────────────────

def test_generic_source_name_h_produces_warning(happy_result):
    # Fixture has moleculetype "H" which is in GENERIC_LIGPARGEN_NAMES
    assert any(_FIXTURE_SOURCE_MOLECULETYPE in w for w in happy_result.warnings)


def test_non_generic_source_name_no_generic_warning(tmp_path):
    gro = _minimal_gro(tmp_path / "a.gro", residue_name="MYLIG")
    itp = _minimal_itp(tmp_path / "a.itp", residue_name="MYLIG", moleculetype_name="MYLIG")
    v = LigParGenImportValidator()
    r = v.validate(gro, itp, _make_identity(), tmp_path / "out")
    assert not any("Generic" in w or "generic" in w for w in r.warnings)


# ── Normalized output files ───────────────────────────────────────────────────

def test_normalized_gro_exists(happy_result):
    assert happy_result.gro_path is not None
    assert happy_result.gro_path.exists()
    assert happy_result.gro_path.name == "L01.gro"


def test_normalized_itp_exists(happy_result):
    assert happy_result.itp_path is not None
    assert happy_result.itp_path.exists()
    assert happy_result.itp_path.name == "L01.itp"


# ── Normalized GRO / ITP consistency ─────────────────────────────────────────

def test_normalized_gro_residue_names_are_l01(happy_result):
    gro = parse_gro(happy_result.gro_path)
    assert all(a.residue_name == "L01" for a in gro.atoms)


def test_normalized_itp_moleculetype_is_l01(happy_result):
    itp = parse_itp(happy_result.itp_path)
    assert itp.moleculetype is not None
    assert itp.moleculetype.name == "L01"


def test_normalized_itp_atoms_residue_names_are_l01(happy_result):
    itp = parse_itp(happy_result.itp_path)
    assert all(a.residue_name == "L01" for a in itp.atoms)


# ── Result fields ─────────────────────────────────────────────────────────────

def test_result_atom_count_matches_fixture(happy_result):
    assert happy_result.atom_count == _FIXTURE_ATOM_COUNT


def test_result_molecule_name_is_source_moleculetype(happy_result):
    assert happy_result.molecule_name == _FIXTURE_SOURCE_MOLECULETYPE


def test_result_total_charge_is_populated(happy_result):
    # Total charge should be near 0.0 for this neutral molecule
    assert abs(happy_result.total_charge) < 0.01
