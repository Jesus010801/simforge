"""Tests for the unified declarative system specification (core.system_spec)."""

from __future__ import annotations

from textwrap import dedent

import pytest

from core.system_spec import SystemSpecError, load_system_spec


_ITP = dedent("""\
    [ atomtypes ]
      opls_900  C900  12.011  0.000  A  3.5E-01  2.7E-01
    [ moleculetype ]
    LIG   3
    [ atoms ]
       1  opls_900  1  LIG  C1  1  0.0  12.011
""")

_GRO = dedent("""\
    lig
       1
        1LIG     C1    1   1.000   1.000   1.000
       5.0 5.0 5.0
""")

_PROT_CHAIN = dedent("""\
    [ moleculetype ]
    Protein_chain_A   3
    [ atoms ]
       1  opls_235  1  ALA  N  1  -0.5  14.007
       2  opls_236  1  ALA  CA 1   0.1  12.011
""")

_PROT_GRO = dedent("""\
    protein
       2
        1ALA      N    1   1.000   1.000   1.000
        1ALA     CA    2   1.100   1.000   1.000
       5.0 5.0 5.0
""")


def _write(tmp_path, **files):
    for name, content in files.items():
        (tmp_path / name).write_text(content)


def _spec(tmp_path, body: str):
    p = tmp_path / "system.yaml"
    p.write_text(dedent(body))
    return load_system_spec(p)


# ── raw protein mode ─────────────────────────────────────────────────────────

def test_valid_raw_protein_one_component(tmp_path):
    _write(tmp_path, **{"prot.pdb": "END\n", "A6.itp": _ITP, "A6.gro": _GRO})
    spec = _spec(tmp_path, """
        protein:
          structure: prot.pdb
        components:
          - id: A6
            topology: A6.itp
            coordinates: A6.gro
            role: inhibitor
    """)
    assert spec.protein_mode == "raw_pdb"
    assert spec.components[0].role == "inhibitor"
    assert spec.component_topology_paths["A6"].is_absolute()


def test_raw_protein_with_topology_is_rejected(tmp_path):
    _write(tmp_path, **{"prot.pdb": "END\n", "A6.itp": _ITP, "A6.gro": _GRO,
                        "topol_Protein_chain_A.itp": _PROT_CHAIN})
    with pytest.raises(SystemSpecError, match="must not be set for a raw"):
        _spec(tmp_path, """
            protein:
              structure: prot.pdb
              topology: [topol_Protein_chain_A.itp]
            components:
              - {id: A6, topology: A6.itp, coordinates: A6.gro}
        """)


# ── preparameterized protein mode ────────────────────────────────────────────

def test_valid_preparameterized_multi_component(tmp_path):
    _write(tmp_path, **{
        "protein_only.gro": _PROT_GRO,
        "topol_Protein_chain_A.itp": _PROT_CHAIN,
        "A6.itp": _ITP, "A6.gro": _GRO,
        "COA.itp": _ITP.replace("LIG", "COA").replace("opls_900", "opls_901"),
        "COA.gro": _GRO.replace("LIG", "COA"),
    })
    spec = _spec(tmp_path, """
        system:
          name: competitive
          protein:
            structure: protein_only.gro
            topology: [topol_Protein_chain_A.itp]
          components:
            - {id: A6, topology: A6.itp, coordinates: A6.gro, role: inhibitor}
            - {id: COA, topology: COA.itp, coordinates: COA.gro, role: cofactor}
          forcefield: oplsaa
          water_model: spce
    """)
    assert spec.protein_mode == "preparameterized"
    assert spec.name == "competitive"
    assert len(spec.components) == 2


def test_preparameterized_without_topology_is_rejected(tmp_path):
    _write(tmp_path, **{"protein_only.gro": _PROT_GRO, "A6.itp": _ITP, "A6.gro": _GRO})
    with pytest.raises(SystemSpecError, match="protein.topology is required"):
        _spec(tmp_path, """
            protein:
              structure: protein_only.gro
            components:
              - {id: A6, topology: A6.itp, coordinates: A6.gro}
        """)


# ── auto-discovery ──────────────────────────────────────────────────────────

def test_topology_auto_discovery_chains(tmp_path):
    _write(tmp_path, **{
        "protein_only.gro": _PROT_GRO,
        "topol_Protein_chain_A.itp": _PROT_CHAIN,
        "topol_Protein_chain_B.itp": _PROT_CHAIN.replace("chain_A", "chain_B"),
        "posre_Protein_chain_A.itp": "; posre\n",
        "posre_Protein_chain_B.itp": "; posre\n",
        "A6.itp": _ITP, "A6.gro": _GRO,
    })
    spec = _spec(tmp_path, """
        protein:
          structure: protein_only.gro
          topology: auto
          restraints: auto
        components:
          - {id: A6, topology: A6.itp, coordinates: A6.gro}
    """)
    assert [p.name for p in spec.protein_topology_paths] == [
        "topol_Protein_chain_A.itp", "topol_Protein_chain_B.itp"]
    assert len(spec.protein_restraint_paths) == 2
    assert any("auto-discovered" in n for n in spec.discovery_notes)


def test_ambiguous_auto_discovery_fails_with_candidates(tmp_path):
    _write(tmp_path, **{
        "protein_only.gro": _PROT_GRO,
        "topol_Protein_chain_A.itp": _PROT_CHAIN,
        "topol_Protein.itp": _PROT_CHAIN,
        "A6.itp": _ITP, "A6.gro": _GRO,
    })
    with pytest.raises(SystemSpecError, match="ambiguous"):
        _spec(tmp_path, """
            protein:
              structure: protein_only.gro
              topology: auto
            components:
              - {id: A6, topology: A6.itp, coordinates: A6.gro}
        """)


def test_auto_discovery_no_candidates_fails(tmp_path):
    _write(tmp_path, **{"protein_only.gro": _PROT_GRO, "A6.itp": _ITP, "A6.gro": _GRO})
    with pytest.raises(SystemSpecError, match="found no topol"):
        _spec(tmp_path, """
            protein:
              structure: protein_only.gro
              topology: auto
            components:
              - {id: A6, topology: A6.itp, coordinates: A6.gro}
        """)


# ── missing files / bad combos ──────────────────────────────────────────────

def test_missing_component_coordinates_rejected(tmp_path):
    _write(tmp_path, **{"prot.pdb": "END\n", "A6.itp": _ITP})
    with pytest.raises(SystemSpecError, match="coordinates not found"):
        _spec(tmp_path, """
            protein: {structure: prot.pdb}
            components:
              - {id: A6, topology: A6.itp, coordinates: A6.gro}
        """)


def test_no_components_rejected(tmp_path):
    _write(tmp_path, **{"prot.pdb": "END\n"})
    with pytest.raises(SystemSpecError, match="at least one component"):
        _spec(tmp_path, "protein: {structure: prot.pdb}\ncomponents: []\n")


def test_unsupported_structure_extension_rejected(tmp_path):
    _write(tmp_path, **{"prot.xyz": "x\n", "A6.itp": _ITP, "A6.gro": _GRO})
    with pytest.raises(SystemSpecError, match="unsupported extension"):
        _spec(tmp_path, """
            protein: {structure: prot.xyz}
            components:
              - {id: A6, topology: A6.itp, coordinates: A6.gro}
        """)


def test_duplicate_component_ids_rejected(tmp_path):
    _write(tmp_path, **{"prot.pdb": "END\n", "A6.itp": _ITP, "A6.gro": _GRO})
    with pytest.raises(SystemSpecError, match="duplicate component id"):
        _spec(tmp_path, """
            protein: {structure: prot.pdb}
            components:
              - {id: A6, topology: A6.itp, coordinates: A6.gro}
              - {id: A6, topology: A6.itp, coordinates: A6.gro}
        """)


def test_relative_paths_resolved_against_yaml_dir(tmp_path):
    sub = tmp_path / "inputs"
    sub.mkdir()
    _write(sub, **{"prot.pdb": "END\n", "A6.itp": _ITP, "A6.gro": _GRO})
    p = tmp_path / "system.yaml"
    p.write_text(dedent("""
        protein: {structure: inputs/prot.pdb}
        components:
          - {id: A6, topology: inputs/A6.itp, coordinates: inputs/A6.gro}
    """))
    spec = load_system_spec(p)
    assert spec.protein_structure_path == sub / "prot.pdb"
