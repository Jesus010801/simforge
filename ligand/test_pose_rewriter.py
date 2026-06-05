"""
Tests for ligand/pose_rewriter.py — LigandPoseRewriter.

Uses real A1 LigParGen fixtures:
  - tests/fixtures/ligpargen/a1/A1_ligpargen_input.pdb  (pose PDB, Å)
  - tests/fixtures/ligpargen/a1/A1_ligpargen.gro        (raw LigParGen GRO)
  - tests/fixtures/ligpargen/a1/A1_ligpargen.itp        (raw LigParGen ITP)

Normalization to L01 is applied via normalize_ligpargen_outputs before the
pose rewriter runs, mirroring the real workflow order.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from core.ligand_workflow_models import LigandPoseRewriteResult
from ligand.normalization import LigandIdentity, normalize_ligpargen_outputs
from ligand.pose_rewriter import LigandPoseRewriter
from utils.gro_parser import parse_gro
from utils.itp_parser import parse_itp

FIXTURE_DIR = Path(__file__).parent.parent / "tests" / "fixtures" / "ligpargen" / "a1"
FIXTURE_GRO = FIXTURE_DIR / "A1_ligpargen.gro"
FIXTURE_ITP = FIXTURE_DIR / "A1_ligpargen.itp"
FIXTURE_PDB = FIXTURE_DIR / "A1_ligpargen_input.pdb"

_ATOM_COUNT = 35


# ── Module-scoped fixtures ────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def normalized(tmp_path_factory):
    """Normalize A1 GRO/ITP outputs to the L01 identity."""
    out = tmp_path_factory.mktemp("norm_a1_pose")
    identity = LigandIdentity(
        component_id="a1",
        display_name="A1",
        source_filename="A1_ligpargen.itp",
        internal_id="L01",
        residue_name="L01",
        moleculetype="L01",
    )
    return normalize_ligpargen_outputs(FIXTURE_GRO, FIXTURE_ITP, identity, out)


@pytest.fixture(scope="module")
def rewrite_result(normalized, tmp_path_factory) -> LigandPoseRewriteResult:
    """Run pose rewriter with A1 input PDB as pose source."""
    out = tmp_path_factory.mktemp("pose_rewrite_a1")
    return LigandPoseRewriter().rewrite(
        pose_pdb=FIXTURE_PDB,
        reference_gro=normalized.normalized_gro,
        reference_itp=normalized.normalized_itp,
        output_dir=out,
    )


# ── Basic success ─────────────────────────────────────────────────────────────

def test_rewrite_succeeds(rewrite_result):
    assert rewrite_result.success is True
    assert rewrite_result.error is None


def test_output_path_exists(rewrite_result):
    assert rewrite_result.output_path is not None
    assert rewrite_result.output_path.exists()


def test_output_filename_contains_residue_name(rewrite_result):
    assert rewrite_result.output_path.name == "L01_pose.gro"


# ── Atom count preservation ───────────────────────────────────────────────────

def test_preserves_atom_count(rewrite_result):
    gro = parse_gro(rewrite_result.output_path)
    assert gro.atom_count == _ATOM_COUNT


def test_atoms_written_field_matches_count(rewrite_result):
    assert rewrite_result.atoms_written == _ATOM_COUNT


# ── GRO atom order and names ──────────────────────────────────────────────────

def test_preserves_gro_atom_names(rewrite_result, normalized):
    ref = parse_gro(normalized.normalized_gro)
    out = parse_gro(rewrite_result.output_path)
    assert [a.atom_name for a in out.atoms] == [a.atom_name for a in ref.atoms]


def test_preserves_gro_atom_numbers(rewrite_result, normalized):
    ref = parse_gro(normalized.normalized_gro)
    out = parse_gro(rewrite_result.output_path)
    assert [a.atom_number for a in out.atoms] == [a.atom_number for a in ref.atoms]


def test_preserves_gro_residue_numbers(rewrite_result, normalized):
    ref = parse_gro(normalized.normalized_gro)
    out = parse_gro(rewrite_result.output_path)
    assert [a.residue_number for a in out.atoms] == [a.residue_number for a in ref.atoms]


# ── Residue name ──────────────────────────────────────────────────────────────

def test_all_atoms_have_l01_residue_name(rewrite_result):
    gro = parse_gro(rewrite_result.output_path)
    assert all(a.residue_name == "L01" for a in gro.atoms)


def test_result_ligand_residue_name_is_l01(rewrite_result):
    assert rewrite_result.ligand_residue_name == "L01"


# ── Coordinate transfer from PDB ──────────────────────────────────────────────

def test_coordinates_transferred_from_pdb(rewrite_result, normalized):
    """Output coordinates should differ from reference GRO (they come from PDB frame)."""
    ref = parse_gro(normalized.normalized_gro)
    out = parse_gro(rewrite_result.output_path)
    max_diff = max(
        abs(r.x - o.x) + abs(r.y - o.y) + abs(r.z - o.z)
        for r, o in zip(ref.atoms, out.atoms)
    )
    # PDB pose is in a different coordinate frame than the reference GRO
    assert max_diff > 0.001  # nm


# ── RMSD ─────────────────────────────────────────────────────────────────────

def test_rmsd_is_computed(rewrite_result):
    assert rewrite_result.rmsd_from_reference is not None


def test_heavy_atom_rmsd_near_zero(rewrite_result):
    """A1 input PDB and LigParGen GRO represent the same conformation.

    After element-based matching + Kabsch alignment the heavy-atom RMSD
    must be small (< 0.05 nm = 0.5 Å).
    """
    assert rewrite_result.rmsd_from_reference < 0.05


# ── Valid GRO output ──────────────────────────────────────────────────────────

def test_output_parseable_by_gro_parser(rewrite_result):
    gro = parse_gro(rewrite_result.output_path)
    assert gro.atom_count == _ATOM_COUNT


def test_output_box_matches_reference(rewrite_result, normalized):
    ref = parse_gro(normalized.normalized_gro)
    out = parse_gro(rewrite_result.output_path)
    assert out.box == pytest.approx(ref.box, abs=1e-6)


# ── Reference files not modified ──────────────────────────────────────────────

def test_normalized_gro_not_modified(normalized):
    gro = parse_gro(normalized.normalized_gro)
    assert all(a.residue_name == "L01" for a in gro.atoms)
    assert gro.atom_count == _ATOM_COUNT


def test_normalized_itp_not_modified(normalized):
    itp = parse_itp(normalized.normalized_itp)
    assert itp.moleculetype is not None
    assert itp.moleculetype.name == "L01"
    assert len(itp.atoms) == _ATOM_COUNT


def test_raw_fixture_gro_not_modified():
    """Raw fixture must keep its original LigParGen residue name."""
    gro = parse_gro(FIXTURE_GRO)
    assert all(a.residue_name == "H" for a in gro.atoms)


# ── Error handling ────────────────────────────────────────────────────────────

def test_missing_pose_pdb_returns_failure(normalized, tmp_path):
    result = LigandPoseRewriter().rewrite(
        pose_pdb=tmp_path / "nonexistent.pdb",
        reference_gro=normalized.normalized_gro,
        reference_itp=normalized.normalized_itp,
        output_dir=tmp_path,
    )
    assert result.success is False
    assert result.error is not None


def test_missing_reference_gro_returns_failure(normalized, tmp_path):
    result = LigandPoseRewriter().rewrite(
        pose_pdb=FIXTURE_PDB,
        reference_gro=tmp_path / "nonexistent.gro",
        reference_itp=normalized.normalized_itp,
        output_dir=tmp_path,
    )
    assert result.success is False
    assert result.error is not None


def test_missing_reference_itp_returns_failure(normalized, tmp_path):
    result = LigandPoseRewriter().rewrite(
        pose_pdb=FIXTURE_PDB,
        reference_gro=normalized.normalized_gro,
        reference_itp=tmp_path / "nonexistent.itp",
        output_dir=tmp_path,
    )
    assert result.success is False
    assert result.error is not None
