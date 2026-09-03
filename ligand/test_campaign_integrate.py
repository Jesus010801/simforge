"""
Tests for ligand/campaign_integrate.py — batch-integrate orchestration.

assemble_system() is patched (same convention as the rest of the ligand test
suite -- see ligand/test_integrate.py) so these tests run without GROMACS.
Pose reconstruction runs for real against the actual A1 LigParGen fixture
plus a heavy-atom-only pose derived from it, so the ligand-topology-reuse
and per-system-pose paths are exercised end to end up to assembly.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from ligand.campaign import prepare_campaign
from ligand.campaign_integrate import integrate_campaign
from ligand.integrate import AssemblyResult

FIXTURE_DIR = Path(__file__).parent.parent / "tests" / "fixtures" / "ligpargen" / "a1"
FIXTURE_GRO = FIXTURE_DIR / "A1_ligpargen.gro"
FIXTURE_ITP = FIXTURE_DIR / "A1_ligpargen.itp"
FIXTURE_PDB = FIXTURE_DIR / "A1_ligpargen_input.pdb"

_RDKIT_AVAIL = "ligand.hydrogenation._rdkit_available"
_OBABEL_AVAIL = "ligand.hydrogenation._obabel_available"


def _pdb_line(record, serial, atom_name, resname, chain, resseq, x, y, z, element):
    return (
        f"{record:<6}{serial:>5} {atom_name:>4} {resname:>3} {chain:1}{resseq:>4}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {element:>2}"
    )


def _protein_lines() -> list[str]:
    lines = []
    serial = 1
    for name, element in (("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")):
        lines.append(_pdb_line("ATOM", serial, name, "ALA", "A", 1, 1.0, 0.0, 0.0, element))
        serial += 1
    return lines


def _a1_heavy_hetatm_lines(center=(30.0, 0.0, 0.0), resname="A1") -> list[str]:
    """Heavy-atom HETATM lines for the real A1 ligand, translated to `center`."""
    heavy = []
    for line in FIXTURE_PDB.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        if len(line) >= 78 and line[76:78].strip().upper() == "H":
            continue
        heavy.append(line)

    # Recenter around the fixture's own centroid, then translate to `center`.
    xs = [float(l[30:38]) for l in heavy]
    ys = [float(l[38:46]) for l in heavy]
    zs = [float(l[46:54]) for l in heavy]
    cx, cy, cz = sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)
    dx, dy, dz = center[0] - cx, center[1] - cy, center[2] - cz

    out = []
    for i, line in enumerate(heavy, start=1):
        x = float(line[30:38]) + dx
        y = float(line[38:46]) + dy
        z = float(line[46:54]) + dz
        name = line[12:16].strip()
        element = line[76:78].strip() if len(line) >= 78 else (name[0] if name else "C")
        out.append(_pdb_line("HETATM", 100 + i, name, resname, "A", 900, x, y, z, element))
    return out


def _write_complex(path: Path, center) -> Path:
    lines = _protein_lines() + ["TER"] + _a1_heavy_hetatm_lines(center=center) + ["END"]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


def _prepared_campaign_root(tmp_path, system_centers: dict[str, tuple]) -> Path:
    """Run the real prepare_campaign() (hydrogenation mocked off) to get a
    valid manifest + directory layout, matching how batch-integrate expects
    its input."""
    for sid, center in system_centers.items():
        _write_complex(tmp_path / sid / f"{sid}.pdb", center)

    with patch(_RDKIT_AVAIL, return_value=False), patch(_OBABEL_AVAIL, return_value=False):
        result = prepare_campaign(tmp_path, "A1")
    assert result.success, result.error
    return tmp_path


def _place_fake_ligpargen_outputs(root: Path, ligand_id: str = "A1") -> None:
    ligpargen_dir = root / "simforge_campaign" / "ligands" / ligand_id / "ligpargen"
    ligpargen_dir.mkdir(parents=True, exist_ok=True)
    (ligpargen_dir / f"{ligand_id}.gro").write_bytes(FIXTURE_GRO.read_bytes())
    (ligpargen_dir / f"{ligand_id}.itp").write_bytes(FIXTURE_ITP.read_bytes())


def _fake_assembly_result(out_dir: Path) -> AssemblyResult:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    complex_gro = out_dir / "complex.gro"
    complex_gro.write_text("fake\n")
    return AssemblyResult(
        success=True,
        complex_gro=complex_gro,
        topol_top=out_dir / "topol.top",
        protein_atom_count=8,
        ligand_atom_count=35,
        total_atom_count=43,
    )


# ═══════════════════════════════════════════════════════════════════════════

class TestIntegrateCampaign:
    def test_no_manifest_fails_clearly(self, tmp_path):
        (tmp_path / "AA").mkdir()
        result = integrate_campaign(tmp_path)
        assert not result.success
        assert "batch-prepare" in result.error

    def test_single_itp_reused_across_systems(self, tmp_path):
        root = _prepared_campaign_root(tmp_path, {"AA": (30, 0, 0), "AG": (60, 5, 2), "LP": (10, -8, 3)})
        _place_fake_ligpargen_outputs(root)

        with patch("ligand.campaign_integrate.assemble_system") as mocked:
            mocked.side_effect = lambda **kw: _fake_assembly_result(kw["out_dir"])
            result = integrate_campaign(root)

        assert result.success, result.error
        assert mocked.call_count == 3
        # Every call must reuse the exact same normalized ligand .itp file.
        itp_paths = {call.kwargs["ligand_itp"] for call in mocked.call_args_list}
        assert len(itp_paths) == 1

    def test_receptor_specific_gro_generated_per_system(self, tmp_path):
        root = _prepared_campaign_root(tmp_path, {"AA": (30, 0, 0), "AG": (60, 5, 2)})
        _place_fake_ligpargen_outputs(root)

        with patch("ligand.campaign_integrate.assemble_system") as mocked:
            mocked.side_effect = lambda **kw: _fake_assembly_result(kw["out_dir"])
            result = integrate_campaign(root)

        assert result.success, result.error
        gro_paths = [call.kwargs["ligand_gro"] for call in mocked.call_args_list]
        assert len(set(gro_paths)) == 2  # distinct per-system pose files
        for p in gro_paths:
            assert Path(p).name == "A1_pose_parameterized.gro"
            assert Path(p).exists()

    def test_assemble_system_invoked_once_per_system(self, tmp_path):
        root = _prepared_campaign_root(
            tmp_path, {"AA": (30, 0, 0), "AG": (60, 5, 2), "HMG-R": (0, 20, 0), "LP": (10, -8, 3)}
        )
        _place_fake_ligpargen_outputs(root)

        with patch("ligand.campaign_integrate.assemble_system") as mocked:
            mocked.side_effect = lambda **kw: _fake_assembly_result(kw["out_dir"])
            result = integrate_campaign(root)

        assert mocked.call_count == 4
        assert result.systems_total == 4
        assert result.systems_succeeded == 4

    def test_partial_failure_reported_per_system(self, tmp_path):
        root = _prepared_campaign_root(tmp_path, {"AA": (30, 0, 0), "AG": (60, 5, 2)})
        _place_fake_ligpargen_outputs(root)

        def _side_effect(**kw):
            if "AG" in str(kw["out_dir"]):
                return AssemblyResult(success=False, error="grompp failed for AG")
            return _fake_assembly_result(kw["out_dir"])

        with patch("ligand.campaign_integrate.assemble_system", side_effect=_side_effect):
            result = integrate_campaign(root)

        assert not result.success
        outcome = result.ligands["A1"]
        assert outcome.systems["AA"].success is True
        assert outcome.systems["AG"].success is False
        assert "grompp failed for AG" in outcome.systems["AG"].error
        assert result.systems_succeeded == 1
        assert result.systems_total == 2

    def test_missing_ligpargen_files_reported(self, tmp_path):
        root = _prepared_campaign_root(tmp_path, {"AA": (30, 0, 0)})
        # Deliberately do not place ligpargen outputs.
        result = integrate_campaign(root)
        assert not result.success
        outcome = result.ligands["A1"]
        assert outcome.success is False
        assert "could not find" in outcome.error.lower()

    def test_invalid_ligpargen_topology_reported(self, tmp_path):
        root = _prepared_campaign_root(tmp_path, {"AA": (30, 0, 0)})
        ligpargen_dir = root / "simforge_campaign" / "ligands" / "A1" / "ligpargen"
        ligpargen_dir.mkdir(parents=True, exist_ok=True)
        (ligpargen_dir / "A1.itp").write_text("; empty, no moleculetype or atoms\n")
        (ligpargen_dir / "A1.gro").write_text("bad\n1\n garbage line\n0 0 0\n")

        result = integrate_campaign(root)
        assert not result.success
        assert result.ligands["A1"].error is not None

    def test_mismatched_heavy_atom_count_fails_pose_transfer(self, tmp_path):
        """A receptor pose with the wrong heavy-atom count must fail cleanly
        during pose reconstruction rather than silently truncating."""
        root = _prepared_campaign_root(tmp_path, {"AA": (30, 0, 0)})
        _place_fake_ligpargen_outputs(root)

        # Corrupt AA's ligand pose to drop one heavy atom.
        pose_path = root / "AA" / "simforge" / "ligand_pose.pdb"
        lines = pose_path.read_text().splitlines()
        het = [l for l in lines if l.startswith("HETATM")]
        kept = het[:-1]
        other = [l for l in lines if not l.startswith("HETATM")]
        pose_path.write_text("\n".join(kept + other) + "\n")

        result = integrate_campaign(root)
        assert not result.success
        outcome = result.ligands["A1"]
        assert outcome.systems["AA"].success is False
        assert "mismatch" in outcome.systems["AA"].error.lower()

    def test_manifest_missing_files_fails_clearly(self, tmp_path):
        root = _prepared_campaign_root(tmp_path, {"AA": (30, 0, 0)})
        _place_fake_ligpargen_outputs(root)
        (root / "AA" / "simforge" / "protein_only.pdb").unlink()

        result = integrate_campaign(root)
        assert not result.success
        assert "missing" in result.error.lower()
