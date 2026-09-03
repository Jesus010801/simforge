"""Tests for core.system_validate — structured validation of a built system dir."""

from __future__ import annotations

from textwrap import dedent

import pytest

from core.system_validate import BLOCKING, FAIL, PASS, SKIP, WARN, validate_system_dir


_ATOMTYPES = dedent("""\
    [ atomtypes ]
      sf_test_c  C  12.011  0.000  A  3.5E-01  2.7E-01
      sf_test_n  N  14.007  0.000  A  3.2E-01  7.1E-01
""")

_LIG_ITP = dedent("""\
    [ moleculetype ]
    LIG   3
    [ atoms ]
       1  sf_test_c  1  LIG  C1  1  0.0  12.011
       2  sf_test_n  1  LIG  N1  1  0.0  14.007
""")

_TOPOL = dedent("""\
    #include "oplsaa.ff/forcefield.itp"
    #include "component_atomtypes.itp"
    #include "LIG.itp"
    #include "oplsaa.ff/spce.itp"

    [ system ]
    test

    [ molecules ]
    LIG   1
    SOL   2
""")

_COMPLEX_GRO = dedent("""\
    complex
       2
        1LIG     C1    1   1.000   1.000   1.000
        1LIG     N1    2   1.200   1.000   1.000
       5.0 5.0 5.0
""")


def _build_dir(tmp_path, *, topol=_TOPOL, gro=_COMPLEX_GRO, lig=_LIG_ITP, atomtypes=_ATOMTYPES):
    (tmp_path / "topol.top").write_text(topol)
    (tmp_path / "complex.gro").write_text(gro)
    (tmp_path / "LIG.itp").write_text(lig)
    (tmp_path / "component_atomtypes.itp").write_text(atomtypes)
    return tmp_path


def _status(report, check_id):
    return next(c.status for c in report.checks if c.check_id == check_id)


def test_valid_system_passes_structural_checks(tmp_path):
    d = _build_dir(tmp_path)
    report = validate_system_dir(d, run_grompp=False)
    assert _status(report, "topology_present") == PASS
    assert _status(report, "coordinate_present") == PASS
    assert _status(report, "includes_resolve") == PASS
    assert _status(report, "molecule_table") == PASS
    assert _status(report, "atom_counts") in (PASS, SKIP)  # SOL unsized -> SKIP
    assert _status(report, "coordinate_validity") == PASS
    assert report.ok


def test_missing_include_is_blocking_failure(tmp_path):
    d = _build_dir(tmp_path)
    (d / "LIG.itp").unlink()
    report = validate_system_dir(d, run_grompp=False)
    assert _status(report, "includes_resolve") == FAIL
    assert not report.ok


def test_atom_count_mismatch_detected(tmp_path):
    # topology LIG has 2 atoms; give a 3-atom coordinate file and no SOL row
    gro = dedent("""\
        complex
           3
            1LIG     C1    1   1.000   1.000   1.000
            1LIG     N1    2   1.200   1.000   1.000
            1LIG     X1    3   1.400   1.000   1.000
           5.0 5.0 5.0
    """)
    topol = _TOPOL.replace("SOL   2\n", "")
    d = _build_dir(tmp_path, topol=topol, gro=gro)
    report = validate_system_dir(d, run_grompp=False)
    assert _status(report, "atom_counts") == FAIL
    assert not report.ok


def test_undefined_moleculetype_in_molecules_table(tmp_path):
    topol = _TOPOL.replace("LIG   1", "GHOST 1")
    d = _build_dir(tmp_path, topol=topol)
    report = validate_system_dir(d, run_grompp=False)
    assert _status(report, "molecule_table") == FAIL


def test_missing_atomtype_detected_against_forcefield(tmp_path):
    # LIG references sf_test_c/sf_test_n but atomtypes only defines sf_test_c
    d = _build_dir(tmp_path, atomtypes="[ atomtypes ]\n  sf_test_c  C  12.011  0.0  A  3.5E-01  2.7E-01\n")
    report = validate_system_dir(d, forcefield="oplsaa", run_grompp=False)
    st = _status(report, "atomtype_resolution")
    # sf_test_n exists neither locally nor in oplsaa -> FAIL when FF data present, else WARN
    assert st in (FAIL, WARN)


def test_nan_coordinate_is_blocking(tmp_path):
    gro = dedent("""\
        complex
           2
            1LIG     C1    1   1.000   1.000   1.000
            1LIG     N1    2     nan   1.000   1.000
           5.0 5.0 5.0
    """)
    d = _build_dir(tmp_path, gro=gro)
    report = validate_system_dir(d, run_grompp=False)
    assert _status(report, "coordinate_validity") == FAIL
    assert not report.ok


def test_exact_overlap_is_advisory_warning(tmp_path):
    gro = dedent("""\
        complex
           2
            1LIG     C1    1   1.000   1.000   1.000
            1LIG     N1    2   1.000   1.000   1.000
           5.0 5.0 5.0
    """)
    d = _build_dir(tmp_path, gro=gro)
    report = validate_system_dir(d, run_grompp=False)
    cv = next(c for c in report.checks if c.check_id == "coordinate_validity")
    assert cv.status == WARN
    assert report.ok  # advisory, not blocking


def test_report_is_json_serializable(tmp_path):
    import json
    d = _build_dir(tmp_path)
    report = validate_system_dir(d, run_grompp=False)
    json.dumps(report.to_dict())


def test_no_topology_is_immediate_failure(tmp_path):
    report = validate_system_dir(tmp_path, run_grompp=False)
    assert _status(report, "topology_present") == FAIL
    assert not report.ok


@pytest.mark.skipif(__import__("shutil").which("gmx") is None, reason="gmx not in PATH")
def test_grompp_runs_when_gmx_present(tmp_path):
    d = _build_dir(tmp_path)
    report = validate_system_dir(d, run_grompp=True)
    assert any(c.check_id == "grompp_dry_run" for c in report.checks)
