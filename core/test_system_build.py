"""Tests for core.system_build — high-level SystemSpec -> assembled system."""

from __future__ import annotations

import json
from textwrap import dedent

import pytest

from core.system_spec import load_system_spec
from core.system_build import build_system, preflight_component_atomtypes
from ligand.integrate import ParameterizedMolecule


_LIG_ITP = dedent("""\
    [ atomtypes ]
      opls_900  C900  12.011  0.000  A  3.50000E-01  2.76144E-01
    [ moleculetype ]
    LIG   3
    [ atoms ]
    ;  nr  type  resnr  res  atom  cgnr  charge  mass
       1  opls_900  1  LIG  C1  1  0.000  12.011
""")

_COA_ITP = dedent("""\
    [ atomtypes ]
      opls_901  N901  14.007  0.000  A  3.25000E-01  7.11280E-01
    [ moleculetype ]
    COA   3
    [ atoms ]
       1  opls_901  1  COA  N1  1  0.000  14.007
""")

_LIG_GRO = dedent("""\
    lig
       1
        1LIG     C1    1   2.000   2.000   2.000
       6.0 6.0 6.0
""")

_COA_GRO = dedent("""\
    coa
       1
        1COA     N1    1   3.000   3.000   3.000
       6.0 6.0 6.0
""")

_PROT_CHAIN_A = dedent("""\
    [ moleculetype ]
    Protein_chain_A   3
    [ atoms ]
       1  opls_238  1  ALA  N   1  -0.50  14.007
       2  opls_239  1  ALA  H   1   0.30   1.008
       3  opls_240  1  ALA  CA  1   0.10  12.011
""")

_PROT_GRO = dedent("""\
    protein
       3
        1ALA      N    1   1.000   1.000   1.000
        1ALA      H    2   1.050   1.080   1.000
        1ALA     CA    3   1.150   1.000   1.000
       6.0 6.0 6.0
""")


def _make_preparam_system(tmp_path):
    (tmp_path / "protein_only.gro").write_text(_PROT_GRO)
    (tmp_path / "topol_Protein_chain_A.itp").write_text(_PROT_CHAIN_A)
    (tmp_path / "A6.itp").write_text(_LIG_ITP)
    (tmp_path / "A6.gro").write_text(_LIG_GRO)
    (tmp_path / "COA.itp").write_text(_COA_ITP)
    (tmp_path / "COA.gro").write_text(_COA_GRO)
    p = tmp_path / "system.yaml"
    p.write_text(dedent("""
        system:
          name: comp
          protein:
            structure: protein_only.gro
            topology: auto
          components:
            - {id: A6, topology: A6.itp, coordinates: A6.gro, role: inhibitor}
            - {id: COA, topology: COA.itp, coordinates: COA.gro, role: cofactor}
          forcefield: oplsaa
          water_model: spce
    """))
    return p


def test_build_records_dirty_source_tree_advisory(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "core.system_build.build_provenance",
        lambda **kw: {"environment": {"git": {"dirty": True, "commit": "abc123def456",
                                              "branch": "wip"}}, "warnings": []},
    )
    spec = load_system_spec(_make_preparam_system(tmp_path))
    outcome = build_system(spec, tmp_path / "built", run_grompp=False)
    assert outcome.provenance["reproducibility"]["clean_source_tree"] is False
    assert any("uncommitted changes" in w for w in outcome.warnings)


def test_build_records_clean_source_tree(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "core.system_build.build_provenance",
        lambda **kw: {"environment": {"git": {"dirty": False}}, "warnings": []},
    )
    spec = load_system_spec(_make_preparam_system(tmp_path))
    outcome = build_system(spec, tmp_path / "built", run_grompp=False)
    assert outcome.provenance["reproducibility"] == {"clean_source_tree": True}


def test_build_preparameterized_protein_plus_two_components(tmp_path):
    spec_path = _make_preparam_system(tmp_path)
    spec = load_system_spec(spec_path)
    out = tmp_path / "built"

    outcome = build_system(spec, out, run_grompp=False, command=["simforge", "build", str(spec_path)])

    assert outcome.success, outcome.errors
    assert (out / "topol.top").exists()
    assert (out / "complex.gro").exists()
    assert (out / "provenance.json").exists()
    assert (out / "validation_report.json").exists()

    # total atoms = 3 protein + 1 + 1
    assert outcome.assembly_result.total_atom_count == 5

    top = (out / "topol.top").read_text()
    # atomtypes included before component itps
    assert top.index("component_atomtypes.itp") < top.index("A6.itp")
    assert "Protein_chain_A" in top

    prov = json.loads((out / "provenance.json").read_text())
    assert prov["protein_mode"] == "preparameterized"
    assert prov["pdb2gmx_used"] is False
    assert prov["component_molecule_names"] == {"A6": "LIG", "COA": "COA"}
    assert prov["inputs"]["component_A6_topology"]["sha256"]
    assert any("auto-discovered" in n for n in prov["discovery_notes"])
    assert prov["environment"]["simforge_version"]


def test_build_records_atomtype_conflict_resolution(tmp_path):
    # A6 and COA both define opls_900 with *different* params -> namespaced
    (tmp_path / "protein_only.gro").write_text(_PROT_GRO)
    (tmp_path / "topol_Protein_chain_A.itp").write_text(_PROT_CHAIN_A)
    (tmp_path / "A6.itp").write_text(_LIG_ITP)
    (tmp_path / "A6.gro").write_text(_LIG_GRO)
    coa_conflict = _COA_ITP.replace("opls_901", "opls_900").replace("N901", "N900")
    (tmp_path / "COA.itp").write_text(coa_conflict)
    (tmp_path / "COA.gro").write_text(_COA_GRO)
    p = tmp_path / "system.yaml"
    p.write_text(dedent("""
        protein: {structure: protein_only.gro, topology: [topol_Protein_chain_A.itp]}
        components:
          - {id: A6, topology: A6.itp, coordinates: A6.gro}
          - {id: COA, topology: COA.itp, coordinates: COA.gro}
    """))
    spec = load_system_spec(p)
    outcome = build_system(spec, tmp_path / "built", run_grompp=False)
    assert outcome.success, outcome.errors
    renames = outcome.provenance["atomtype_renames"]
    assert renames  # at least one component had a rename recorded
    flat = {k: v for cid in renames for k, v in renames[cid].items()}
    assert "opls_900" in flat


# ── Phase 6 preflight ───────────────────────────────────────────────────────

def test_preflight_flags_component_with_no_atomtypes_block(tmp_path):
    bad = dedent("""\
        [ moleculetype ]
        BAD   3
        [ atoms ]
           1  lp_zzz_1  1  BAD  C1  1  0.0  12.011
    """)
    itp = tmp_path / "BAD.itp"
    itp.write_text(bad)
    gro = tmp_path / "BAD.gro"
    gro.write_text(_LIG_GRO.replace("LIG", "BAD"))
    comps = [ParameterizedMolecule("BAD", itp, gro)]
    errs = preflight_component_atomtypes(comps, "oplsaa")
    assert errs
    assert "lp_zzz_1" in errs[0]
    assert "atomtypes" in errs[0].lower()


def test_preflight_passes_for_selfcontained_component(tmp_path):
    itp = tmp_path / "A6.itp"
    itp.write_text(_LIG_ITP)
    gro = tmp_path / "A6.gro"
    gro.write_text(_LIG_GRO)
    errs = preflight_component_atomtypes([ParameterizedMolecule("A6", itp, gro)], "oplsaa")
    assert errs == []


def test_build_aborts_on_preflight_failure(tmp_path):
    (tmp_path / "protein_only.gro").write_text(_PROT_GRO)
    (tmp_path / "topol_Protein_chain_A.itp").write_text(_PROT_CHAIN_A)
    bad = dedent("""\
        [ moleculetype ]
        BAD   3
        [ atoms ]
           1  lp_missing_9  1  BAD  C1  1  0.0  12.011
    """)
    (tmp_path / "BAD.itp").write_text(bad)
    (tmp_path / "BAD.gro").write_text(_LIG_GRO.replace("LIG", "BAD"))
    p = tmp_path / "system.yaml"
    p.write_text(dedent("""
        protein: {structure: protein_only.gro, topology: [topol_Protein_chain_A.itp]}
        components:
          - {id: BAD, topology: BAD.itp, coordinates: BAD.gro}
    """))
    spec = load_system_spec(p)
    outcome = build_system(spec, tmp_path / "built", run_grompp=False)
    assert not outcome.success
    assert outcome.preflight_errors
    assert not (tmp_path / "built" / "topol.top").exists()


def test_build_raw_pdb_mode_maps_to_protein_pdb_kwarg(tmp_path, monkeypatch):
    (tmp_path / "prot.pdb").write_text("END\n")
    (tmp_path / "A6.itp").write_text(_LIG_ITP)
    (tmp_path / "A6.gro").write_text(_LIG_GRO)
    p = tmp_path / "system.yaml"
    p.write_text(dedent("""
        protein: {structure: prot.pdb}
        components:
          - {id: A6, topology: A6.itp, coordinates: A6.gro}
    """))
    spec = load_system_spec(p)

    seen = {}

    def _fake(**kw):
        seen.update(kw)
        raise RuntimeError("stop-after-capture")

    monkeypatch.setattr("core.system_build.assemble_system_multi", _fake)
    with pytest.raises(RuntimeError, match="stop-after-capture"):
        build_system(spec, tmp_path / "built", run_grompp=False)

    assert seen["protein_pdb"] == spec.protein_structure_path
    assert "protein_gro" not in seen
    assert seen["forcefield"] == "oplsaa"


@pytest.mark.skipif(
    __import__("shutil").which("gmx") is None, reason="gmx not in PATH",
)
def test_build_grompp_validation_when_gmx_present(tmp_path):
    spec_path = _make_preparam_system(tmp_path)
    spec = load_system_spec(spec_path)
    outcome = build_system(spec, tmp_path / "built", run_grompp=True)
    grompp = [c for c in outcome.validation["checks"] if c["check_id"] == "grompp_dry_run"]
    assert grompp and grompp[0]["status"] in ("PASS", "WARN", "FAIL")
