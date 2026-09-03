"""The shipped examples must stay runnable from a fresh clone.

- structural checks (no GROMACS needed): every referenced file exists and the
  spec loads + resolves.
- @pytest.mark.gmx: a real end-to-end build of the multicomponent example that
  doubles as the multi-component / atomtype-namespace-isolation regression —
  preparameterized protein + two independently parameterized components whose
  ITPs collide on `opls_800`, assembled and passed through
  `gmx grompp -maxwarn 0`.
"""
from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from core.system_spec import load_system_spec

_ROOT = Path(__file__).resolve().parent.parent
_EXAMPLES = _ROOT / "examples"
_CASES = ["protein_ligand_min", "multicomponent"]


@pytest.mark.parametrize("name", _CASES)
def test_example_spec_loads_and_all_inputs_present(name):
    spec_path = _EXAMPLES / name / "system.yaml"
    assert spec_path.is_file(), f"missing {spec_path}"

    spec = load_system_spec(spec_path)  # parses + resolves + validates paths

    assert Path(spec.protein_structure_path).is_file()
    for p in spec.protein_topology_paths:
        assert Path(p).is_file(), f"missing protein topology {p}"
    for p in spec.protein_restraint_paths:
        assert Path(p).is_file(), f"missing restraint {p}"
    assert spec.components, "example declares no components"
    for c in spec.components:
        assert Path(spec.component_topology_paths[c.id]).is_file()
        assert Path(spec.component_coordinate_paths[c.id]).is_file()


def test_multicomponent_example_has_the_atomtype_collision():
    """The whole point of the multicomponent example: both ITPs define opls_800."""
    mth = (_EXAMPLES / "multicomponent" / "MTH.itp").read_text()
    eol = (_EXAMPLES / "multicomponent" / "EOL.itp").read_text()
    assert "opls_800" in mth and "opls_800" in eol


@pytest.mark.gmx
@pytest.mark.skipif(shutil.which("gmx") is None, reason="gmx not on PATH")
def test_multicomponent_example_builds_and_grompps_clean(tmp_path):
    from core.system_build import build_system

    spec = load_system_spec(_EXAMPLES / "multicomponent" / "system.yaml")
    outcome = build_system(spec, tmp_path / "multi", run_grompp=True)

    assert outcome.success, outcome.errors
    assert (tmp_path / "multi" / "topol.top").is_file()
    assert (tmp_path / "multi" / "provenance.json").is_file()

    # atomtype namespace isolation actually happened and is recorded
    renames = outcome.provenance["atomtype_renames"]
    assert renames["MTH"]["opls_800"] == "MTH_opls_800"
    assert renames["EOL"]["opls_800"] == "EOL_opls_800"

    # the grompp gate really ran and really passed at -maxwarn 0
    checks = {c["check_id"]: c for c in outcome.validation["checks"]}
    assert checks["grompp_dry_run"]["status"] == "PASS", checks["grompp_dry_run"]
    assert outcome.validation["readiness"] in ("READY", "ready")
