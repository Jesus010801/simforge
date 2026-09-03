"""
Tests for LigandPoseRewriter.rewrite_heavy_atom_transfer() — receptor-specific
pose reconstruction for the ligand campaign workflow (Phase 10-12).

Two fixture sources are used:
  - A synthetic 5-heavy-atom chain (C1-C2-C3-C4-C5, one H each) built in this
    file, used to test correctness of the per-atom local-frame hydrogen
    reconstruction under a genuinely non-rigid conformational change (the
    two chain "arms" are bent independently around the fixed C3 pivot —
    something a single whole-molecule Kabsch rotation cannot satisfy).
  - The real A1 LigParGen fixture (tests/fixtures/ligpargen/a1/), with its
    pose PDB's hydrogens stripped to simulate a receptor pose extracted
    straight from a docked complex (heavy atoms only, as described in the
    campaign workflow).
"""

from __future__ import annotations

from pathlib import Path
from textwrap import dedent

import numpy as np
import pytest

from ligand.normalization import LigandIdentity, normalize_ligpargen_outputs
from ligand.pose_rewriter import LigandPoseRewriter
from utils.gro_parser import parse_gro

FIXTURE_DIR = Path(__file__).parent.parent / "tests" / "fixtures" / "ligpargen" / "a1"
FIXTURE_GRO = FIXTURE_DIR / "A1_ligpargen.gro"
FIXTURE_ITP = FIXTURE_DIR / "A1_ligpargen.itp"
FIXTURE_PDB = FIXTURE_DIR / "A1_ligpargen_input.pdb"


# ═══════════════════════════════════════════════════════════════════════════
# Synthetic 5-carbon chain fixture
# ═══════════════════════════════════════════════════════════════════════════
#
# Heavy skeleton: C1-C2-C3-C4-C5 (straight-ish zigzag), one H per carbon.
# Chosen so each hydrogen's local reconstruction frame (center + 1/2-bond
# heavy neighbours) differs between the two ends:
#   H1 (on C1) frame -> {C1, C2, C3}
#   H5 (on C5) frame -> {C3, C4, C5}
# which lets a test apply two *different* rotations to each half of the
# chain (pivoting on C3) and verify each hydrogen follows its own local arm.

_HEAVY_REF = {
    "C1": np.array([-2.0, 0.0, 0.0]),
    "C2": np.array([-1.0, 0.3, 0.0]),
    "C3": np.array([0.0, 0.0, 0.0]),
    "C4": np.array([1.0, 0.3, 0.0]),
    "C5": np.array([2.0, 0.0, 0.0]),
}
_H_REF = {
    "H1": np.array([-2.5, 0.7, 0.3]),
    "H2": np.array([-1.0, 1.0, -0.5]),
    "H3": np.array([0.0, -0.8, 0.4]),
    "H4": np.array([1.0, 1.0, -0.5]),
    "H5": np.array([2.5, 0.7, 0.3]),
}
_ATOM_ORDER = ["C1", "C2", "C3", "C4", "C5", "H1", "H2", "H3", "H4", "H5"]
_ATOM_TYPE = {n: ("opls_c" if n.startswith("C") else "opls_h") for n in _ATOM_ORDER}
_MASS = {"C1": 12.011, "C2": 12.011, "C3": 12.011, "C4": 12.011, "C5": 12.011,
         "H1": 1.008, "H2": 1.008, "H3": 1.008, "H4": 1.008, "H5": 1.008}
_BONDS = [
    ("C1", "C2"), ("C2", "C3"), ("C3", "C4"), ("C4", "C5"),
    ("C1", "H1"), ("C2", "H2"), ("C3", "H3"), ("C4", "H4"), ("C5", "H5"),
]


def _rotz(theta: float) -> np.ndarray:
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _write_chain_itp(path: Path, bonds=None) -> None:
    idx = {name: i + 1 for i, name in enumerate(_ATOM_ORDER)}
    lines = [
        "[ atomtypes ]",
        "  opls_c  C  12.0110  0.000  A  3.50000E-01  2.76144E-01",
        "  opls_h  H   1.0080  0.000  A  2.50000E-01  1.25520E-01",
        "",
        "[ moleculetype ]",
        "; Name  nrexcl",
        "MOL     3",
        "",
        "[ atoms ]",
        ";  nr type resnr residue atom cgnr charge mass",
    ]
    for name in _ATOM_ORDER:
        lines.append(
            f"  {idx[name]:>2} {_ATOM_TYPE[name]:<7} 1  MOL {name:<4} 1  0.000  {_MASS[name]:.4f}"
        )
    lines.append("")
    lines.append("[ bonds ]")
    for a, b in (bonds if bonds is not None else _BONDS):
        lines.append(f"  {idx[a]:>2} {idx[b]:>2}  1  0.1400  300000.0")
    path.write_text("\n".join(lines) + "\n")


def _write_chain_gro(path: Path, coords: dict[str, np.ndarray]) -> None:
    lines = ["Synthetic C5 chain", f"{len(_ATOM_ORDER):5d}"]
    for i, name in enumerate(_ATOM_ORDER, start=1):
        x, y, z = coords[name]
        lines.append(f"{1:5d}{'MOL':<5s}{name:>5s}{i:5d}{x:8.3f}{y:8.3f}{z:8.3f}")
    lines.append(f"{5.00000:8.5f}{5.00000:8.5f}{5.00000:8.5f}")
    path.write_text("\n".join(lines) + "\n")


def _write_pose_pdb(path: Path, heavy_coords_nm: dict[str, np.ndarray], resname="MOL") -> None:
    """Heavy-atom-only PDB (Å), as produced by campaign pose extraction."""
    lines = []
    for i, name in enumerate(["C1", "C2", "C3", "C4", "C5"], start=1):
        x, y, z = heavy_coords_nm[name] * 10.0
        lines.append(
            f"HETATM{i:5d} {name:<4} {resname:<3} A 900    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
        )
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


@pytest.fixture()
def chain_reference(tmp_path) -> tuple[Path, Path]:
    itp = tmp_path / "chain.itp"
    gro = tmp_path / "chain.gro"
    all_coords = {**_HEAVY_REF, **_H_REF}
    _write_chain_itp(itp)
    _write_chain_gro(gro, all_coords)
    return gro, itp


# ═══════════════════════════════════════════════════════════════════════════
# Rigid case: identity transform (heavy atoms unchanged)
# ═══════════════════════════════════════════════════════════════════════════

def test_rigid_identity_pose_preserves_heavy_atoms(chain_reference, tmp_path):
    ref_gro, ref_itp = chain_reference
    pose_pdb = tmp_path / "pose.pdb"
    _write_pose_pdb(pose_pdb, _HEAVY_REF)

    result = LigandPoseRewriter().rewrite_heavy_atom_transfer(
        reference_gro=ref_gro, reference_itp=ref_itp,
        target_pose_pdb=pose_pdb, output_dir=tmp_path / "out",
    )

    assert result.success, result.error
    assert result.heavy_atom_mapping_method == "atom_order"
    assert result.hydrogens_reconstructed is True

    out = parse_gro(result.output_path)
    for i, name in enumerate(["C1", "C2", "C3", "C4", "C5"]):
        got = np.array([out.atoms[i].x, out.atoms[i].y, out.atoms[i].z])
        assert got == pytest.approx(_HEAVY_REF[name], abs=1e-3)


def test_atoms_written_matches_reference_count(chain_reference, tmp_path):
    ref_gro, ref_itp = chain_reference
    pose_pdb = tmp_path / "pose.pdb"
    _write_pose_pdb(pose_pdb, _HEAVY_REF)

    result = LigandPoseRewriter().rewrite_heavy_atom_transfer(
        reference_gro=ref_gro, reference_itp=ref_itp,
        target_pose_pdb=pose_pdb, output_dir=tmp_path / "out",
    )
    assert result.atoms_written == 10


def test_output_preserves_parameterized_atom_ordering(chain_reference, tmp_path):
    ref_gro, ref_itp = chain_reference
    pose_pdb = tmp_path / "pose.pdb"
    _write_pose_pdb(pose_pdb, _HEAVY_REF)

    result = LigandPoseRewriter().rewrite_heavy_atom_transfer(
        reference_gro=ref_gro, reference_itp=ref_itp,
        target_pose_pdb=pose_pdb, output_dir=tmp_path / "out",
    )
    out = parse_gro(result.output_path)
    ref = parse_gro(ref_gro)
    assert [a.atom_name for a in out.atoms] == [a.atom_name for a in ref.atoms]
    assert [a.atom_number for a in out.atoms] == [a.atom_number for a in ref.atoms]


# ═══════════════════════════════════════════════════════════════════════════
# Flexible case: two different local rotations on either side of a fixed
# pivot (C3) — a genuine conformational change no single global rotation
# can satisfy exactly for both arms simultaneously.
# ═══════════════════════════════════════════════════════════════════════════

def test_flexible_conformation_heavy_atoms_transferred_exactly(chain_reference, tmp_path):
    ref_gro, ref_itp = chain_reference
    pivot = _HEAVY_REF["C3"]
    R_left = _rotz(np.deg2rad(40))
    R_right = _rotz(np.deg2rad(-40))

    target_heavy = {
        "C1": pivot + R_left @ (_HEAVY_REF["C1"] - pivot),
        "C2": pivot + R_left @ (_HEAVY_REF["C2"] - pivot),
        "C3": pivot,
        "C4": pivot + R_right @ (_HEAVY_REF["C4"] - pivot),
        "C5": pivot + R_right @ (_HEAVY_REF["C5"] - pivot),
    }
    pose_pdb = tmp_path / "pose.pdb"
    _write_pose_pdb(pose_pdb, target_heavy)

    result = LigandPoseRewriter().rewrite_heavy_atom_transfer(
        reference_gro=ref_gro, reference_itp=ref_itp,
        target_pose_pdb=pose_pdb, output_dir=tmp_path / "out",
    )
    assert result.success, result.error

    out = parse_gro(result.output_path)
    name_to_atom = {a.atom_name: a for a in out.atoms}
    for name, expected in target_heavy.items():
        got = np.array([name_to_atom[name].x, name_to_atom[name].y, name_to_atom[name].z])
        assert got == pytest.approx(expected, abs=1e-3), f"{name} heavy atom was not transferred exactly"


def test_flexible_conformation_hydrogens_follow_correct_local_arm(chain_reference, tmp_path):
    """H1 (near C1) must follow the LEFT-arm rotation; H5 (near C5) the RIGHT-arm one.

    A single whole-molecule Kabsch rotation would place both hydrogens
    using the same (incorrect, averaged) rotation -- this test fails under
    that naive approach and only passes with per-atom local-frame repositioning.
    """
    ref_gro, ref_itp = chain_reference
    pivot = _HEAVY_REF["C3"]
    R_left = _rotz(np.deg2rad(40))
    R_right = _rotz(np.deg2rad(-40))

    target_heavy = {
        "C1": pivot + R_left @ (_HEAVY_REF["C1"] - pivot),
        "C2": pivot + R_left @ (_HEAVY_REF["C2"] - pivot),
        "C3": pivot,
        "C4": pivot + R_right @ (_HEAVY_REF["C4"] - pivot),
        "C5": pivot + R_right @ (_HEAVY_REF["C5"] - pivot),
    }
    pose_pdb = tmp_path / "pose.pdb"
    _write_pose_pdb(pose_pdb, target_heavy)

    result = LigandPoseRewriter().rewrite_heavy_atom_transfer(
        reference_gro=ref_gro, reference_itp=ref_itp,
        target_pose_pdb=pose_pdb, output_dir=tmp_path / "out",
    )
    assert result.success, result.error
    out = parse_gro(result.output_path)
    name_to_atom = {a.atom_name: a for a in out.atoms}

    # Expected H1: C1's local frame is {C1,C2,C3}, entirely within the LEFT
    # arm (which rotates rigidly by R_left about the pivot) -> exact.
    expected_h1 = target_heavy["C1"] + R_left @ (_H_REF["H1"] - _HEAVY_REF["C1"])
    got_h1 = np.array([name_to_atom["H1"].x, name_to_atom["H1"].y, name_to_atom["H1"].z])
    assert got_h1 == pytest.approx(expected_h1, abs=1e-3)

    # Expected H5: C5's local frame is {C3,C4,C5}, entirely within the RIGHT arm.
    expected_h5 = target_heavy["C5"] + R_right @ (_H_REF["H5"] - _HEAVY_REF["C5"])
    got_h5 = np.array([name_to_atom["H5"].x, name_to_atom["H5"].y, name_to_atom["H5"].z])
    assert got_h5 == pytest.approx(expected_h5, abs=1e-3)

    # And they must differ from each other's rotation applied to the wrong arm,
    # proving the two hydrogens were NOT repositioned via one shared rotation.
    wrong_h1 = target_heavy["C1"] + R_right @ (_H_REF["H1"] - _HEAVY_REF["C1"])
    assert got_h1 != pytest.approx(wrong_h1, abs=1e-3)


def test_heavy_atoms_never_touched_by_hydrogen_placement(chain_reference, tmp_path):
    """Re-run the flexible case and confirm heavy coordinates are bit-identical
    to the target pose regardless of hydrogen reconstruction — i.e. hydrogen
    placement never perturbs heavy atoms."""
    ref_gro, ref_itp = chain_reference
    pivot = _HEAVY_REF["C3"]
    R = _rotz(np.deg2rad(75))
    target_heavy = {name: pivot + R @ (coord - pivot) for name, coord in _HEAVY_REF.items()}

    pose_pdb = tmp_path / "pose.pdb"
    _write_pose_pdb(pose_pdb, target_heavy)

    result = LigandPoseRewriter().rewrite_heavy_atom_transfer(
        reference_gro=ref_gro, reference_itp=ref_itp,
        target_pose_pdb=pose_pdb, output_dir=tmp_path / "out",
    )
    assert result.success, result.error
    out = parse_gro(result.output_path)
    name_to_atom = {a.atom_name: a for a in out.atoms}
    for name, expected in target_heavy.items():
        got = np.array([name_to_atom[name].x, name_to_atom[name].y, name_to_atom[name].z])
        # .gro format stores 3 decimal places (nm) -- tolerance matches that precision.
        assert got == pytest.approx(expected, abs=1e-3)


# ═══════════════════════════════════════════════════════════════════════════
# Ambiguous / failure paths
# ═══════════════════════════════════════════════════════════════════════════

def test_heavy_atom_count_mismatch_fails_safely(chain_reference, tmp_path):
    ref_gro, ref_itp = chain_reference
    # Drop C5 from the pose -> only 4 heavy atoms instead of 5.
    partial = {k: v for k, v in _HEAVY_REF.items() if k != "C5"}
    pose_pdb = tmp_path / "pose.pdb"
    lines = []
    for i, name in enumerate(["C1", "C2", "C3", "C4"], start=1):
        x, y, z = partial[name] * 10.0
        lines.append(
            f"HETATM{i:5d} {name:<4} MOL A 900    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
        )
    lines.append("END")
    pose_pdb.write_text("\n".join(lines) + "\n")

    result = LigandPoseRewriter().rewrite_heavy_atom_transfer(
        reference_gro=ref_gro, reference_itp=ref_itp,
        target_pose_pdb=pose_pdb, output_dir=tmp_path / "out",
    )
    assert result.success is False
    assert "count mismatch" in result.error.lower()


def test_missing_bonds_section_fails_explicitly(tmp_path):
    itp = tmp_path / "chain.itp"
    gro = tmp_path / "chain.gro"
    _write_chain_itp(itp, bonds=[])  # no [ bonds ] entries written -> empty section is absent
    # Remove the (empty) [ bonds ] header entirely so parse_itp finds no bonds section.
    text = itp.read_text().split("[ bonds ]")[0]
    itp.write_text(text)
    _write_chain_gro(gro, {**_HEAVY_REF, **_H_REF})

    pose_pdb = tmp_path / "pose.pdb"
    _write_pose_pdb(pose_pdb, _HEAVY_REF)

    result = LigandPoseRewriter().rewrite_heavy_atom_transfer(
        reference_gro=gro, reference_itp=itp,
        target_pose_pdb=pose_pdb, output_dir=tmp_path / "out",
    )
    assert result.success is False
    assert "bonds" in result.error.lower()


def test_hydrogen_without_bonded_parent_fails_explicitly(tmp_path):
    itp = tmp_path / "chain.itp"
    gro = tmp_path / "chain.gro"
    # Omit H3's bond entirely -- H3 becomes unreconstructable.
    bonds_without_h3 = [b for b in _BONDS if b != ("C3", "H3")]
    _write_chain_itp(itp, bonds=bonds_without_h3)
    _write_chain_gro(gro, {**_HEAVY_REF, **_H_REF})

    pose_pdb = tmp_path / "pose.pdb"
    _write_pose_pdb(pose_pdb, _HEAVY_REF)

    result = LigandPoseRewriter().rewrite_heavy_atom_transfer(
        reference_gro=gro, reference_itp=itp,
        target_pose_pdb=pose_pdb, output_dir=tmp_path / "out",
    )
    assert result.success is False
    assert "H3" in result.error


def test_missing_target_pose_fails(chain_reference, tmp_path):
    ref_gro, ref_itp = chain_reference
    result = LigandPoseRewriter().rewrite_heavy_atom_transfer(
        reference_gro=ref_gro, reference_itp=ref_itp,
        target_pose_pdb=tmp_path / "nonexistent.pdb", output_dir=tmp_path / "out",
    )
    assert result.success is False
    assert result.error is not None


# ═══════════════════════════════════════════════════════════════════════════
# Real LigParGen fixture (A1) — heavy-atom-only pose, geometric fallback path
# ═══════════════════════════════════════════════════════════════════════════

@pytest.fixture(scope="module")
def a1_normalized(tmp_path_factory):
    out = tmp_path_factory.mktemp("norm_a1_heavy_transfer")
    identity = LigandIdentity(
        component_id="a1", display_name="A1", source_filename="A1_ligpargen.itp",
        internal_id="L01", residue_name="L01", moleculetype="L01",
    )
    return normalize_ligpargen_outputs(FIXTURE_GRO, FIXTURE_ITP, identity, out)


@pytest.fixture(scope="module")
def a1_heavy_only_pose(tmp_path_factory):
    lines = FIXTURE_PDB.read_text().splitlines()
    heavy = [
        l for l in lines
        if not (l[:6].strip() in ("ATOM", "HETATM") and len(l) >= 78 and l[76:78].strip().upper() == "H")
    ]
    out_dir = tmp_path_factory.mktemp("a1_heavy_pose")
    path = out_dir / "a1_heavy_only.pdb"
    path.write_text("\n".join(heavy) + "\n")
    return path


def test_real_fixture_heavy_only_pose_succeeds(a1_normalized, a1_heavy_only_pose, tmp_path):
    result = LigandPoseRewriter().rewrite_heavy_atom_transfer(
        reference_gro=a1_normalized.normalized_gro,
        reference_itp=a1_normalized.normalized_itp,
        target_pose_pdb=a1_heavy_only_pose,
        output_dir=tmp_path,
    )
    assert result.success, result.error
    assert result.atoms_written == 35
    assert result.hydrogens_reconstructed is True


def test_real_fixture_rmsd_small(a1_normalized, a1_heavy_only_pose, tmp_path):
    """Same conformation, hydrogens only stripped -> reconstruction should be
    near-exact (small heavy-atom RMSD after Kabsch alignment)."""
    result = LigandPoseRewriter().rewrite_heavy_atom_transfer(
        reference_gro=a1_normalized.normalized_gro,
        reference_itp=a1_normalized.normalized_itp,
        target_pose_pdb=a1_heavy_only_pose,
        output_dir=tmp_path,
    )
    assert result.success, result.error
    assert result.rmsd_from_reference is not None
    assert result.rmsd_from_reference < 0.05
