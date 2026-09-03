"""Tests for core.system_campaign — declarative multi-system orchestration."""

from __future__ import annotations

from textwrap import dedent

import pytest

from core.system_spec import SystemSpecError
from core.system_campaign import expand_campaign, run_campaign


_LIG = dedent("""\
    [ atomtypes ]
      opls_900  C900  12.011  0.000  A  3.5E-01  2.7E-01
    [ moleculetype ]
    {name}   3
    [ atoms ]
       1  opls_900  1  {name}  C1  1  0.0  12.011
""")
_COA = dedent("""\
    [ atomtypes ]
      opls_901  N901  14.007  0.000  A  3.2E-01  7.1E-01
    [ moleculetype ]
    COA   3
    [ atoms ]
       1  opls_901  1  COA  N1  1  0.0  14.007
""")
def _gro(name: str) -> str:
    atom = f"{1:5d}{name:<5s}{'C1':>5s}{1:5d}{2.0:8.3f}{2.0:8.3f}{2.0:8.3f}"
    return f"m\n{1:5d}\n{atom}\n   6.0   6.0   6.0\n"
_PROT_CHAIN = dedent("""\
    [ moleculetype ]
    Protein_chain_A   3
    [ atoms ]
       1  opls_238  1  ALA  N   1  -0.5  14.007
       2  opls_239  1  ALA  H   1   0.3   1.008
       3  opls_240  1  ALA  CA  1   0.1  12.011
""")
_PROT_GRO = dedent("""\
    protein
       3
        1ALA      N    1   1.000   1.000   1.000
        1ALA      H    2   1.050   1.080   1.000
        1ALA     CA    3   1.150   1.000   1.000
       6.0 6.0 6.0
""")


@pytest.fixture()
def campaign(tmp_path):
    (tmp_path / "protein_only.gro").write_text(_PROT_GRO)
    (tmp_path / "topol_Protein_chain_A.itp").write_text(_PROT_CHAIN)
    for lig in ("A1", "A3"):
        (tmp_path / f"{lig}.itp").write_text(_LIG.format(name=lig))
        (tmp_path / f"{lig}.gro").write_text(_gro(lig))
    (tmp_path / "COA.itp").write_text(_COA)
    (tmp_path / "COA.gro").write_text(_gro("COA"))
    p = tmp_path / "campaign.yaml"
    p.write_text(dedent("""
        campaign:
          name: series
          shared:
            protein:
              structure: protein_only.gro
              topology: auto
            forcefield: oplsaa
            water_model: spce
            components:
              - {id: COA, topology: COA.itp, coordinates: COA.gro, role: cofactor}
          systems:
            - name: A1
              components:
                - {id: A1, topology: A1.itp, coordinates: A1.gro, role: competitive_ligand}
            - name: A3
              components:
                - {id: A3, topology: A3.itp, coordinates: A3.gro, role: competitive_ligand}
    """))
    return p


def test_expand_campaign_merges_shared_protein_and_cofactor(campaign):
    name, specs, _ = expand_campaign(campaign)
    assert name == "series"
    assert [n for n, _ in specs] == ["A1", "A3"]
    a1 = dict(specs)["A1"]
    assert [c.id for c in a1.components] == ["COA", "A1"]  # shared first
    assert a1.protein_mode == "preparameterized"
    assert a1.forcefield == "oplsaa"


def test_run_campaign_builds_every_system(campaign, tmp_path):
    result = run_campaign(campaign, tmp_path / "out", run_grompp=False)
    assert result.total == 2
    assert result.succeeded == 2
    assert result.success
    assert (tmp_path / "out" / "A1" / "topol.top").exists()
    assert (tmp_path / "out" / "A3" / "provenance.json").exists()
    assert (tmp_path / "out" / "campaign_summary.yaml").exists()


def test_campaign_requires_shared_protein(tmp_path):
    p = tmp_path / "c.yaml"
    p.write_text("campaign:\n  shared: {}\n  systems:\n    - name: x\n")
    with pytest.raises(SystemSpecError, match="shared.protein is required"):
        expand_campaign(p)


def test_campaign_empty_systems_rejected(tmp_path):
    p = tmp_path / "c.yaml"
    p.write_text("campaign:\n  shared:\n    protein: {structure: p.gro}\n  systems: []\n")
    with pytest.raises(SystemSpecError, match="non-empty list"):
        expand_campaign(p)


def test_campaign_propagates_per_system_spec_errors(tmp_path):
    (tmp_path / "protein_only.gro").write_text(_PROT_GRO)
    (tmp_path / "topol_Protein_chain_A.itp").write_text(_PROT_CHAIN)
    p = tmp_path / "campaign.yaml"
    p.write_text(dedent("""
        campaign:
          shared:
            protein: {structure: protein_only.gro, topology: [topol_Protein_chain_A.itp]}
          systems:
            - name: A1
              components:
                - {id: A1, topology: missing.itp, coordinates: missing.gro}
    """))
    with pytest.raises(SystemSpecError, match="A1"):
        expand_campaign(p)
