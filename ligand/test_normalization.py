"""
Tests for ligand/normalization.py — LigParGen identity normalization layer.

Uses real LigParGen output from tests/fixtures/ligpargen/a1/.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from ligand.normalization import (
    LigandIdentity,
    LigandNormalizationResult,
    normalize_ligpargen_outputs,
)
from utils.gro_parser import parse_gro
from utils.itp_parser import parse_itp

FIXTURE_DIR = Path(__file__).parent.parent / "tests" / "fixtures" / "ligpargen" / "a1"
FIXTURE_GRO = FIXTURE_DIR / "A1_ligpargen.gro"
FIXTURE_ITP = FIXTURE_DIR / "A1_ligpargen.itp"


@pytest.fixture(scope="module")
def identity_l01() -> LigandIdentity:
    return LigandIdentity(
        component_id="a1",
        display_name="A1",
        source_filename="A1_ligpargen.itp",
        internal_id="L01",
        residue_name="L01",
        moleculetype="L01",
    )


@pytest.fixture(scope="module")
def result(identity_l01, tmp_path_factory) -> LigandNormalizationResult:
    out = tmp_path_factory.mktemp("normalization_a1")
    return normalize_ligpargen_outputs(FIXTURE_GRO, FIXTURE_ITP, identity_l01, out)


# ── Moleculetype normalization ────────────────────────────────────────────────

def test_normalize_moleculetype_h_to_l01(result):
    itp = parse_itp(result.normalized_itp)
    assert itp.moleculetype is not None
    assert itp.moleculetype.name == "L01"


def test_normalize_atoms_residue_name_h_to_l01(result):
    itp = parse_itp(result.normalized_itp)
    assert all(a.residue_name == "L01" for a in itp.atoms)


def test_normalize_gro_residue_names_to_l01(result):
    gro = parse_gro(result.normalized_gro)
    assert all(a.residue_name == "L01" for a in gro.atoms)


# ── Atom count preservation ───────────────────────────────────────────────────

def test_atom_count_unchanged_gro(result):
    original = parse_gro(FIXTURE_GRO)
    normalized = parse_gro(result.normalized_gro)
    assert normalized.atom_count == original.atom_count


def test_atom_count_unchanged_itp(result):
    original = parse_itp(FIXTURE_ITP)
    normalized = parse_itp(result.normalized_itp)
    assert len(normalized.atoms) == len(original.atoms)


# ── Coordinate preservation ───────────────────────────────────────────────────

def test_coordinates_unchanged(result):
    original = parse_gro(FIXTURE_GRO)
    normalized = parse_gro(result.normalized_gro)
    for orig, norm in zip(original.atoms, normalized.atoms):
        assert norm.x == pytest.approx(orig.x, abs=1e-4)
        assert norm.y == pytest.approx(orig.y, abs=1e-4)
        assert norm.z == pytest.approx(orig.z, abs=1e-4)


# ── Atom name preservation ────────────────────────────────────────────────────

def test_atom_names_unchanged_gro(result):
    original = parse_gro(FIXTURE_GRO)
    normalized = parse_gro(result.normalized_gro)
    assert [a.atom_name for a in normalized.atoms] == [a.atom_name for a in original.atoms]


def test_atom_names_unchanged_itp(result):
    original = parse_itp(FIXTURE_ITP)
    normalized = parse_itp(result.normalized_itp)
    assert [a.atom_name for a in normalized.atoms] == [a.atom_name for a in original.atoms]


# ── Charge preservation ───────────────────────────────────────────────────────

def test_total_charge_unchanged(result):
    original = parse_itp(FIXTURE_ITP)
    normalized = parse_itp(result.normalized_itp)
    assert normalized.total_charge == pytest.approx(original.total_charge, abs=1e-4)


# ── Warning behavior ──────────────────────────────────────────────────────────

def test_generic_source_name_produces_warning(result):
    assert result.source_was_generic is True
    assert len(result.warnings) > 0
    assert any("H" in w for w in result.warnings)


def test_non_generic_source_normalized_and_records_original(tmp_path):
    """Non-generic moleculetype is still normalized but source name is preserved."""
    minimal_itp = tmp_path / "nongeneric.itp"
    minimal_itp.write_text(
        "[ moleculetype ]\n"
        "; Name  nrexcl\n"
        "MYLIG   3\n"
        "[ atoms ]\n"
        ";   nr  type  resnr  residue  atom  cgnr  charge  mass\n"
        "     1   CT    1    MYLIG     C1    1   0.0000  12.0110\n"
    )
    minimal_gro = tmp_path / "nongeneric.gro"
    # GRO fixed-width: resnum(5) resname(5) atomname(5) atomnum(5) x(8.3f) y(8.3f) z(8.3f)
    minimal_gro.write_text(
        "Minimal ligand\n"
        "   1\n"
        "    1MYLIG   C1    1   1.000   1.000   1.000\n"
        "   1.00000   1.00000   1.00000\n"
    )
    identity = LigandIdentity(
        component_id="a1",
        display_name="A1",
        source_filename="nongeneric.itp",
        internal_id="L01",
        residue_name="L01",
        moleculetype="L01",
    )
    r = normalize_ligpargen_outputs(minimal_gro, minimal_itp, identity, tmp_path / "out_ng")

    assert r.source_was_generic is False
    assert r.identity.source_moleculetype == "MYLIG"
    assert parse_itp(r.normalized_itp).moleculetype.name == "L01"


# ── Output files ──────────────────────────────────────────────────────────────

def test_output_files_exist(result):
    assert result.normalized_gro.exists()
    assert result.normalized_itp.exists()
    assert result.identity_json.exists()
    assert result.normalization_report.exists()


def test_output_filenames(result):
    assert result.normalized_gro.name == "L01.gro"
    assert result.normalized_itp.name == "L01.itp"


# ── identity.json ─────────────────────────────────────────────────────────────

def test_identity_json_has_source_and_normalized_names(result):
    data = json.loads(result.identity_json.read_text())
    assert data["source_moleculetype"] == "H"
    assert data["moleculetype"] == "L01"
    assert data["residue_name"] == "L01"
    assert data["internal_id"] == "L01"
