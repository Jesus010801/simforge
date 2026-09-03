"""
CLI tests for the platform-consolidation commands:
    simforge build / doctor / validate-system / inspect-run
"""

from __future__ import annotations

import json
from textwrap import dedent

import pytest
from typer.testing import CliRunner

from cli import cli

try:  # click >= 8.2 removed the mix_stderr kwarg (stderr is separate by default)
    runner = CliRunner(mix_stderr=False)
except TypeError:
    runner = CliRunner()


_LIG_ITP = dedent("""\
    [ atomtypes ]
      opls_900  C900  12.011  0.000  A  3.50000E-01  2.76144E-01
    [ moleculetype ]
    LIG   3
    [ atoms ]
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
_LIG_GRO = "lig\n   1\n    1LIG     C1    1   2.000   2.000   2.000\n   6.0 6.0 6.0\n"
_COA_GRO = "coa\n   1\n    1COA     N1    1   3.000   3.000   3.000\n   6.0 6.0 6.0\n"
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


@pytest.fixture()
def preparam_spec(tmp_path):
    (tmp_path / "protein_only.gro").write_text(_PROT_GRO)
    (tmp_path / "topol_Protein_chain_A.itp").write_text(_PROT_CHAIN_A)
    (tmp_path / "A6.itp").write_text(_LIG_ITP)
    (tmp_path / "A6.gro").write_text(_LIG_GRO)
    (tmp_path / "COA.itp").write_text(_COA_ITP)
    (tmp_path / "COA.gro").write_text(_COA_GRO)
    spec = tmp_path / "system.yaml"
    spec.write_text(dedent("""
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
    return spec


# ── doctor ──────────────────────────────────────────────────────────────────

def test_doctor_json_mode():
    result = runner.invoke(cli, ["doctor", "--json"])
    payload = json.loads(result.stdout)
    assert "overall" in payload and "checks" in payload
    names = {c["name"] for c in payload["checks"]}
    assert {"SimForge", "Python", "GROMACS", "RDKit"} <= names


def test_doctor_human_output_runs():
    result = runner.invoke(cli, ["doctor"])
    assert "SIMFORGE ENVIRONMENT" in result.stdout
    assert "Overall status" in result.stdout


# ── build ───────────────────────────────────────────────────────────────────

def test_build_preparameterized(preparam_spec, tmp_path):
    out = tmp_path / "built"
    result = runner.invoke(cli, ["build", str(preparam_spec), "--out", str(out), "--no-grompp"])
    assert result.exit_code == 0, result.stdout
    assert (out / "topol.top").exists()
    assert (out / "provenance.json").exists()
    assert "Build complete" in result.stdout
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["protein_mode"] == "preparameterized"
    assert prov["command"][:2] == ["simforge", "build"]


def test_build_invalid_spec_reports_clearly(tmp_path):
    spec = tmp_path / "bad.yaml"
    spec.write_text("protein: {structure: nope.gro}\ncomponents: []\n")
    result = runner.invoke(cli, ["build", str(spec), "--out", str(tmp_path / "o")])
    assert result.exit_code == 1
    assert "Invalid system specification" in result.stdout


def test_build_preflight_failure_blocks(tmp_path):
    (tmp_path / "protein_only.gro").write_text(_PROT_GRO)
    (tmp_path / "topol_Protein_chain_A.itp").write_text(_PROT_CHAIN_A)
    (tmp_path / "BAD.itp").write_text(dedent("""\
        [ moleculetype ]
        BAD   3
        [ atoms ]
           1  sf_unknown_1  1  BAD  C1  1  0.0  12.011
    """))
    (tmp_path / "BAD.gro").write_text(_LIG_GRO.replace("LIG", "BAD"))
    spec = tmp_path / "system.yaml"
    spec.write_text(dedent("""
        protein: {structure: protein_only.gro, topology: [topol_Protein_chain_A.itp]}
        components:
          - {id: BAD, topology: BAD.itp, coordinates: BAD.gro}
    """))
    result = runner.invoke(cli, ["build", str(spec), "--out", str(tmp_path / "o"), "--no-grompp"])
    assert result.exit_code == 1
    assert "preflight" in result.stdout.lower()


# ── validate-system ─────────────────────────────────────────────────────────

def test_validate_system_on_built_dir(preparam_spec, tmp_path):
    out = tmp_path / "built"
    runner.invoke(cli, ["build", str(preparam_spec), "--out", str(out), "--no-grompp"])
    result = runner.invoke(cli, ["validate-system", str(out), "--no-grompp"])
    assert result.exit_code == 0, result.stdout
    assert "SYSTEM VALIDATION" in result.stdout
    assert (out / "validation_report.json").exists()


def test_validate_system_missing_topology_exit_code(tmp_path):
    result = runner.invoke(cli, ["validate-system", str(tmp_path), "--no-grompp"])
    assert result.exit_code == 1


# ── inspect-run ─────────────────────────────────────────────────────────────

def test_inspect_run_summary(preparam_spec, tmp_path):
    out = tmp_path / "built"
    runner.invoke(cli, ["build", str(preparam_spec), "--out", str(out), "--no-grompp"])
    result = runner.invoke(cli, ["inspect-run", str(out)])
    assert result.exit_code == 0, result.stdout
    assert "Reproducibility Manifest" in result.stdout
    assert "Inputs" in result.stdout


def test_inspect_run_missing_provenance(tmp_path):
    result = runner.invoke(cli, ["inspect-run", str(tmp_path)])
    assert result.exit_code == 1


# ── backward compatibility ──────────────────────────────────────────────────

def test_existing_commands_still_registered():
    result = runner.invoke(cli, ["--help"])
    for cmd in ("compile", "run", "validate", "inspect", "ligand", "build", "doctor"):
        assert cmd in result.stdout
