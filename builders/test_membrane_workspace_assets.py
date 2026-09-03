# builders/test_membrane_workspace_assets.py
"""
Phase 1 remediation: membrane workspace self-containment.

Verifies that WorkspaceBuilder for protein-membrane systems:
  - Stages dppc512_whole.gro, inflategro-Jorge.pl, oplsaa_membrane.ff/ into
    workspace/membrane_assets/
  - Writes metadata/membrane_assets.json listing staged assets and sources
  - Generated embed script references only workspace-local paths (no Prot-Memb_FILES)
  - Generated shrink-loop script points INFLATEGRO at membrane_assets/
  - Generated topology script sets GMXLIB to membrane_assets/ for pdb2gmx
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from core.compiler import SimulationCompiler
from builders.workspace_builder import WorkspaceBuilder

MEMBRANE_YAML = Path(__file__).parent.parent / "configs" / "membrane_test.yaml"


@pytest.fixture(scope="module")
def membrane_workspace(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("membrane_ws")
    result = SimulationCompiler().compile(str(MEMBRANE_YAML))
    ws = WorkspaceBuilder().build(result, workspace_path=str(tmp))
    return ws


def _find_step_dir(workspace: Path, step_id: str) -> Path | None:
    for d in sorted((workspace / "steps").iterdir()):
        if d.is_dir() and step_id in d.name:
            return d
    return None


# ── Asset presence ─────────────────────────────────────────────────────────────

def test_membrane_assets_dir_exists(membrane_workspace):
    assert (membrane_workspace / "membrane_assets").is_dir()


def test_bilayer_gro_staged(membrane_workspace):
    assert (membrane_workspace / "membrane_assets" / "dppc512_whole.gro").is_file()


def test_inflategro_staged(membrane_workspace):
    assert (membrane_workspace / "membrane_assets" / "inflategro-Jorge.pl").is_file()


def test_ff_dir_staged(membrane_workspace):
    assert (membrane_workspace / "membrane_assets" / "oplsaa_membrane.ff").is_dir()


def test_ff_dir_has_core_files(membrane_workspace):
    ff = membrane_workspace / "membrane_assets" / "oplsaa_membrane.ff"
    assert (ff / "ffnonbonded.itp").is_file()
    assert (ff / "ffbonded.itp").is_file()
    assert (ff / "aminoacids.rtp").is_file()


# ── Metadata ───────────────────────────────────────────────────────────────────

def test_membrane_assets_metadata_exists(membrane_workspace):
    assert (membrane_workspace / "metadata" / "membrane_assets.json").is_file()


def test_membrane_assets_metadata_lists_all_assets(membrane_workspace):
    meta = json.loads(
        (membrane_workspace / "metadata" / "membrane_assets.json").read_text()
    )
    asset_names = {a["asset"] for a in meta["assets"]}
    assert "dppc512_whole.gro" in asset_names
    assert "inflategro-Jorge.pl" in asset_names
    assert "oplsaa_membrane.ff" in asset_names


def test_membrane_assets_metadata_has_source_paths(membrane_workspace):
    meta = json.loads(
        (membrane_workspace / "metadata" / "membrane_assets.json").read_text()
    )
    assert "source_dir" in meta
    assert "staged_at" in meta
    for asset in meta["assets"]:
        assert "source" in asset
        assert "destination" in asset


# ── Scripts: workspace-local paths only ───────────────────────────────────────

def test_embed_script_no_prot_memb_files_ref(membrane_workspace):
    """Generated embed script must not reference Prot-Memb_FILES."""
    embed_dir = _find_step_dir(membrane_workspace, "embed_in_bilayer")
    assert embed_dir is not None, "embed_in_bilayer step directory not found"
    script = (embed_dir / "run_embed.py").read_text()
    assert "Prot-Memb_FILES" not in script


def test_embed_script_uses_assets_dir(membrane_workspace):
    """Generated embed script must resolve bilayer via ASSETS_DIR from membrane_assets."""
    embed_dir = _find_step_dir(membrane_workspace, "embed_in_bilayer")
    script = (embed_dir / "run_embed.py").read_text()
    assert "ASSETS_DIR" in script
    assert "membrane_assets" in script


def test_embed_script_bilayer_path_from_assets(membrane_workspace):
    """bilayer_path must be derived from ASSETS_DIR, not Path(BILAYER_FILE) CWD lookup."""
    embed_dir = _find_step_dir(membrane_workspace, "embed_in_bilayer")
    script = (embed_dir / "run_embed.py").read_text()
    assert "bilayer_path = ASSETS_DIR / BILAYER_FILE" in script


def test_shrink_loop_inflategro_in_membrane_assets(membrane_workspace):
    """INFLATEGRO= in shrink loop must point into membrane_assets/, not bare filename."""
    embed_dir = _find_step_dir(membrane_workspace, "membrane_embedding")
    assert embed_dir is not None, "membrane_embedding step directory not found"
    script = (embed_dir / "run.sh").read_text()
    for line in script.splitlines():
        if line.startswith("INFLATEGRO="):
            assert "membrane_assets" in line, (
                f"INFLATEGRO must point into membrane_assets, got: {line}"
            )
            break
    else:
        pytest.fail("INFLATEGRO= assignment not found in shrink loop run.sh")


def test_topology_script_sets_gmxlib(membrane_workspace):
    """generate_protein_topology must set GMXLIB so pdb2gmx finds oplsaa_membrane.ff."""
    topo_dir = _find_step_dir(membrane_workspace, "generate_protein_topology")
    assert topo_dir is not None, "generate_protein_topology step directory not found"
    script = (topo_dir / "run_protein_topology.py").read_text()
    assert "GMXLIB" in script
    assert "membrane_assets" in script


def test_topology_script_passes_env_to_pdb2gmx(membrane_workspace):
    """pdb2gmx subprocess call must include env=_gmx_env."""
    topo_dir = _find_step_dir(membrane_workspace, "generate_protein_topology")
    assert topo_dir is not None, "generate_protein_topology step directory not found"
    script = (topo_dir / "run_protein_topology.py").read_text()
    assert "env=_gmx_env" in script
