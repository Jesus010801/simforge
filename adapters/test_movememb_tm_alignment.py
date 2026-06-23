# adapters/test_movememb_tm_alignment.py
"""
Phase 2 remediation: TM-aware bilayer alignment.

All tests use synthetic protein geometry chosen so the full-protein
Z-centre and the TM-residue Z-centre differ by a large, measurable
amount.  A protein with a long extracellular domain is the canonical
failure mode: using the full centre mis-places the bilayer core by
several nm.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from adapters.movememb_adapter import MoveMembAdapter, _tm_ca_center_z, _parse_gro


# ── GRO construction helpers ───────────────────────────────────────────────────

def _gro_line(resnum: int, resname: str, atomname: str, serial: int,
              x: float, y: float, z: float) -> str:
    """Format one GRO atom line (no velocity columns)."""
    return f"{resnum:5d}{resname:<5}{atomname:>5}{serial:5d}{x:8.3f}{y:8.3f}{z:8.3f}"


def _write_gro(path: Path, title: str, atoms: list[tuple], box=(20.0, 20.0, 20.0)) -> Path:
    """
    Write a minimal .gro file.
    atoms: list of (resnum, resname, atomname, x, y, z)
    """
    lines = [title, str(len(atoms))]
    for i, (rn, rname, aname, x, y, z) in enumerate(atoms, start=1):
        lines.append(_gro_line(rn, rname, aname, i, x, y, z))
    lines.append(f"{box[0]:.5f}  {box[1]:.5f}  {box[2]:.5f}")
    path.write_text("\n".join(lines) + "\n")
    return path


def _asymmetric_protein(tmp_path: Path) -> tuple[Path, set[int]]:
    """
    Protein with a large EC domain so full-centre ≠ TM-centre.

    Geometry (all CA, residues 1-15):
        Residues  1- 5  TM segment  Z = [3.0, 4.0, 5.0, 6.0, 7.0]  mean = 5.0 nm
        Residues  6-15  EC domain   Z = [8.0 … 17.0]                 mean = 12.5 nm

    TM CA centre  = 5.0 nm
    Full Z extent = [3.0, 17.0], centre = 10.0 nm   ← differs by 5 nm
    """
    atoms = []
    # TM residues 1-5
    for i, z in enumerate([3.0, 4.0, 5.0, 6.0, 7.0], start=1):
        atoms.append((i, "ALA", "CA", 5.0, 5.0, z))
    # EC residues 6-15
    for i, z in enumerate([8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0], start=6):
        atoms.append((i, "ALA", "CA", 5.0, 5.0, z))

    path = _write_gro(tmp_path / "protein.gro", "asymmetric protein", atoms)
    tm_residues = set(range(1, 6))   # residues 1-5
    return path, tm_residues


def _bilayer(tmp_path: Path, midplane_z: float = 10.0) -> Path:
    """Minimal DPPC bilayer with midplane at midplane_z."""
    half = 3.0
    atoms = [
        (1, "DPP", "O33", 5.0, 5.0, midplane_z - half),
        (2, "DPP", "O33", 5.0, 5.0, midplane_z + half),
    ]
    return _write_gro(tmp_path / "bilayer.gro", "DPPC bilayer", atoms)


# ── Unit tests for _tm_ca_center_z ────────────────────────────────────────────

def test_tm_ca_center_z_correct(tmp_path):
    """mean Z of CA atoms in TM residues only."""
    prot, tm_residues = _asymmetric_protein(tmp_path)
    _, lines, _ = _parse_gro(prot)
    z = _tm_ca_center_z(lines, tm_residues)
    assert z is not None
    assert abs(z - 5.0) < 1e-6, f"Expected TM centre 5.0, got {z}"


def test_tm_ca_center_z_excludes_ec_atoms(tmp_path):
    """EC CA atoms must not affect the TM centre."""
    prot, tm_residues = _asymmetric_protein(tmp_path)
    _, lines, _ = _parse_gro(prot)
    tm_z = _tm_ca_center_z(lines, tm_residues)
    # Full protein Z extent = [3, 17], centre = 10.0 — far from 5.0
    assert tm_z is not None
    full_centre = (3.0 + 17.0) / 2.0    # 10.0
    assert abs(tm_z - full_centre) > 3.0, "TM centre must differ significantly from full protein centre"


def test_tm_ca_center_z_returns_none_for_empty_set(tmp_path):
    prot, _ = _asymmetric_protein(tmp_path)
    _, lines, _ = _parse_gro(prot)
    assert _tm_ca_center_z(lines, set()) is None


def test_tm_ca_center_z_ignores_non_ca_atoms(tmp_path):
    """Only CA atoms count; CB / N / O in TM residues must be ignored."""
    atoms = [
        (1, "ALA", "CA",  5.0, 5.0, 4.0),   # CA — should be used
        (1, "ALA", "CB",  5.0, 5.0, 100.0),  # CB — must be ignored
        (1, "ALA", "N",   5.0, 5.0, 200.0),  # N  — must be ignored
    ]
    p = _write_gro(tmp_path / "p.gro", "test", atoms)
    _, lines, _ = _parse_gro(p)
    z = _tm_ca_center_z(lines, {1})
    assert z is not None
    assert abs(z - 4.0) < 1e-6


# ── Adapter: z_shift uses TM centre when annotation present ───────────────────

def test_tm_center_used_for_z_shift(tmp_path):
    """
    With TM annotation the shift must target TM-centre Z, not protein-centre Z.

    Geometry:
        TM CA centre  = 5.0 nm
        Bilayer midplane = 10.0 nm
        Expected shift  = 5.0 - 10.0 = -5.0 nm

    Without annotation:
        Full protein centre = 10.0 nm
        Expected shift      = 10.0 - 10.0 = 0.0 nm  (clearly different)
    """
    prot_gro, tm_residues = _asymmetric_protein(tmp_path)
    bil_gro  = _bilayer(tmp_path, midplane_z=10.0)
    out_gro  = tmp_path / "system.gro"

    adapter = MoveMembAdapter()
    result  = adapter.run(
        protein_gro=prot_gro,
        bilayer_gro=bil_gro,
        gro_out=out_gro,
        tm_residues=tm_residues,
    )

    assert result.success
    m = result.metadata
    assert m["alignment_method"] == "tm_center"
    assert m["tm_annotation_used"] is True
    assert m["tm_center_z"] is not None
    assert abs(m["tm_center_z"] - 5.0) < 1e-6
    assert abs(m["z_shift_nm"] - (-5.0)) < 1e-3


def test_large_ec_domain_does_not_shift_membrane_placement(tmp_path):
    """
    The difference in z_shift between TM-aware and fallback must be large
    (≥ 4 nm for this geometry), proving EC domain atoms are excluded.
    """
    prot_gro, tm_residues = _asymmetric_protein(tmp_path)
    bil_gro  = _bilayer(tmp_path, midplane_z=10.0)

    adapter = MoveMembAdapter()

    result_tm = adapter.run(
        protein_gro=prot_gro,
        bilayer_gro=bil_gro,
        gro_out=tmp_path / "sys_tm.gro",
        tm_residues=tm_residues,
    )
    result_fb = adapter.run(
        protein_gro=prot_gro,
        bilayer_gro=bil_gro,
        gro_out=tmp_path / "sys_fb.gro",
        tm_residues=None,
    )

    assert result_tm.success and result_fb.success

    shift_tm = result_tm.metadata["z_shift_nm"]   # −5.0
    shift_fb = result_fb.metadata["z_shift_nm"]   # 0.0 (protein centre = 10.0 = bilayer midplane)
    assert abs(shift_tm - shift_fb) >= 4.0, (
        f"TM-aware shift ({shift_tm:.3f}) and fallback shift ({shift_fb:.3f}) "
        f"should differ by ≥ 4 nm for asymmetric protein geometry"
    )


# ── Adapter: fallback when no TM annotation ───────────────────────────────────

def test_protein_center_fallback_when_no_annotation(tmp_path):
    """Without tm_residues the adapter uses full protein Z-centre."""
    prot_gro, _ = _asymmetric_protein(tmp_path)
    bil_gro = _bilayer(tmp_path, midplane_z=10.0)

    result = MoveMembAdapter().run(
        protein_gro=prot_gro,
        bilayer_gro=bil_gro,
        gro_out=tmp_path / "sys.gro",
        tm_residues=None,
    )

    assert result.success
    m = result.metadata
    assert m["alignment_method"] == "protein_center_fallback"
    assert m["tm_annotation_used"] is False
    assert m["tm_center_z"] is None
    # shift uses full protein centre = (3+17)/2 = 10.0; bilayer at 10.0 → shift=0.0
    assert abs(m["z_shift_nm"] - 0.0) < 1e-3


# ── Report fields: presence and types ─────────────────────────────────────────

def test_all_report_fields_present_with_tm_annotation(tmp_path):
    prot_gro, tm_residues = _asymmetric_protein(tmp_path)
    bil_gro = _bilayer(tmp_path)

    m = MoveMembAdapter().run(
        protein_gro=prot_gro,
        bilayer_gro=bil_gro,
        gro_out=tmp_path / "sys.gro",
        tm_residues=tm_residues,
    ).metadata

    assert isinstance(m["tm_center_z"],        float)
    assert isinstance(m["protein_center_z"],   float)
    assert isinstance(m["bilayer_midplane_z"],  float)
    assert m["alignment_method"]   == "tm_center"
    assert m["tm_annotation_used"] is True


def test_all_report_fields_present_without_tm_annotation(tmp_path):
    prot_gro, _ = _asymmetric_protein(tmp_path)
    bil_gro = _bilayer(tmp_path)

    m = MoveMembAdapter().run(
        protein_gro=prot_gro,
        bilayer_gro=bil_gro,
        gro_out=tmp_path / "sys.gro",
    ).metadata

    assert m["tm_center_z"]        is None
    assert isinstance(m["protein_center_z"],   float)
    assert isinstance(m["bilayer_midplane_z"],  float)
    assert m["alignment_method"]   == "protein_center_fallback"
    assert m["tm_annotation_used"] is False


# ── Generated script checks ───────────────────────────────────────────────────

def _find_step_dir(workspace: Path, step_id: str) -> Path | None:
    for d in sorted((workspace / "steps").iterdir()):
        if d.is_dir() and step_id in d.name:
            return d
    return None


@pytest.fixture(scope="module")
def membrane_workspace_with_tm(tmp_path_factory):
    """
    Workspace compiled from a YAML that includes TM annotation so
    TM_RESIDUES is baked into the generated run_embed.py.

    The augmented YAML is written to a temp dir alongside a copy of
    peptide_A.pdb so the workspace builder can stage it.
    """
    import shutil
    import textwrap
    from pathlib import Path as _P
    from core.compiler import SimulationCompiler
    from builders.workspace_builder import WorkspaceBuilder

    configs_dir = _P(__file__).parent.parent / "configs"
    yaml_src    = configs_dir / "membrane_test.yaml"
    base_yaml   = yaml_src.read_text()

    annotation_block = textwrap.dedent("""
    structural_annotation:
      membrane_topology:
        transmembrane_segments:
          - "5-25"
    """)
    augmented_yaml = base_yaml.rstrip() + "\n" + annotation_block

    tmp = tmp_path_factory.mktemp("ws_with_tm")
    aug_yaml_path = tmp / "membrane_with_tm.yaml"
    aug_yaml_path.write_text(augmented_yaml)

    # Copy the protein PDB so relative path 'peptide_A.pdb' resolves in tmp
    shutil.copy2(configs_dir / "peptide_A.pdb", tmp / "peptide_A.pdb")

    result = SimulationCompiler().compile(str(aug_yaml_path))
    ws = WorkspaceBuilder().build(result, workspace_path=str(tmp / "workspace"))
    return ws


def test_embed_script_has_tm_residues_set(membrane_workspace_with_tm):
    """run_embed.py must contain TM_RESIDUES as a concrete set, not None."""
    embed_dir = _find_step_dir(membrane_workspace_with_tm, "embed_in_bilayer")
    assert embed_dir is not None
    script = (embed_dir / "run_embed.py").read_text()
    assert "TM_RESIDUES = set(" in script, (
        "TM_RESIDUES should be a concrete set when tm_residues is in step params"
    )


def test_embed_script_passes_tm_residues_to_adapter(membrane_workspace_with_tm):
    """The adapter.run() call must include tm_residues=TM_RESIDUES."""
    embed_dir = _find_step_dir(membrane_workspace_with_tm, "embed_in_bilayer")
    script = (embed_dir / "run_embed.py").read_text()
    assert "tm_residues=TM_RESIDUES" in script


def test_embed_script_no_annotation_has_none(tmp_path_factory):
    """Without TM annotation, TM_RESIDUES must be None in the generated script."""
    from core.compiler import SimulationCompiler
    from builders.workspace_builder import WorkspaceBuilder
    from pathlib import Path as _P

    yaml_path = _P(__file__).parent.parent / "configs" / "membrane_test.yaml"
    tmp       = tmp_path_factory.mktemp("ws_no_tm")
    result    = SimulationCompiler().compile(str(yaml_path))
    ws        = WorkspaceBuilder().build(result, workspace_path=str(tmp))

    embed_dir = _find_step_dir(ws, "embed_in_bilayer")
    assert embed_dir is not None
    script = (embed_dir / "run_embed.py").read_text()
    assert "TM_RESIDUES = None" in script
