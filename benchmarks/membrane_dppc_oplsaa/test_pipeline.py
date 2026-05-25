"""
benchmarks/membrane_dppc_oplsaa/test_pipeline.py

End-to-end compile + MDP verification for the DPPC + OPLS-AA membrane pipeline.

Covers:
  - Correct pipeline selection (MembraneWorkflowOPLSAA, not MDPipeline)
  - Exact step count and stages
  - Critical MDP parameters (semiisotropic, gen_vel, continuation, dt=0.001)
  - Fortran-free adapters available out of the box
  - Fallback: protein-only YAML does NOT select MembraneWorkflowOPLSAA
"""

from __future__ import annotations

import pytest
import tempfile
import textwrap
from pathlib import Path

from core.compiler import SimulationCompiler
from core.execution_models import AutomationLevel
from pipelines.membrane_pipeline import MembraneWorkflowOPLSAA
from pipelines.md_pipeline import MDPipeline
from adapters.water_deletor_adapter import WaterDeletorAdapter, _parse_gro, _delete_bilayer_waters
from adapters.movememb_adapter import MoveMembAdapter


# ── Fixtures ──────────────────────────────────────────────────────────────────

CONFIGS_DIR = Path(__file__).parent.parent.parent / "configs"
MEMBRANE_YAML = CONFIGS_DIR / "membrane_test.yaml"

@pytest.fixture(scope="session")
def membrane_result(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("membrane_run")
    compiler = SimulationCompiler()
    result = compiler.compile(str(MEMBRANE_YAML))
    return result


# ── Pipeline selection ────────────────────────────────────────────────────────

def test_pipeline_type_is_protein_membrane(membrane_result):
    assert membrane_result.state.inferred_system_type == "protein-membrane"


def test_correct_pipeline_selected(membrane_result):
    # pipeline type must trigger MembraneWorkflowOPLSAA, not the generic MDPipeline
    plan = membrane_result.plan
    assert plan.notes, "Pipeline should add notes"
    assert any("MembraneWorkflowOPLSAA" in n for n in plan.notes)


def test_step_count(membrane_result):
    # 11 base steps + 1 analysis_rmsd = 12
    assert len(membrane_result.plan.steps) == 12


def test_step_stages(membrane_result):
    stages = [s.stage.value for s in membrane_result.plan.steps]
    assert "membrane_embedding" in stages
    assert stages.count("equilibration") == 1
    assert stages.count("production") == 1


def test_manual_steps(membrane_result):
    # match_box_to_bilayer — promoted to AUTOMATED (MatchBoxBuilder).
    # embed_in_bilayer    — promoted to AUTOMATED (MoveMembAdapter Python).
    # orient_protein without structural_annotation = GUIDED (step_type=preparation, not manual).
    manual = [s for s in membrane_result.plan.steps if s.step_type.value == "manual"]
    step_ids = [s.step_id for s in manual]
    assert "match_box_to_bilayer" not in step_ids, (
        "match_box_to_bilayer is now AUTOMATED — should not appear in manual steps"
    )
    assert "embed_in_bilayer" not in step_ids, (
        "embed_in_bilayer is now AUTOMATED (MoveMembAdapter) — should not appear in manual steps"
    )


def test_embed_in_bilayer_is_automated(membrane_result):
    embed = next(s for s in membrane_result.plan.steps if s.step_id == "embed_in_bilayer")
    assert embed.step_type.value == "automatic"
    assert embed.automation_level == AutomationLevel.AUTOMATED


# ── MDP parameter verification ────────────────────────────────────────────────

@pytest.fixture(scope="session")
def workspace(membrane_result):
    runs_dir = Path("simforge_runs/protein-membrane")
    if not runs_dir.exists():
        pytest.skip("Workspace not materialized — run compile first")
    return runs_dir


def test_nvt_mdp_gen_vel(workspace):
    mdp = (workspace / "steps/10_equilibration/nvt.mdp").read_text()
    assert "gen_vel                 = yes" in mdp


def test_nvt_mdp_tc_grps_system(workspace):
    mdp = (workspace / "steps/10_equilibration/nvt.mdp").read_text()
    assert "tc-grps                 = system" in mdp


def test_nvt_mdp_rcoulomb_12(workspace):
    mdp = (workspace / "steps/10_equilibration/nvt.mdp").read_text()
    assert "rcoulomb                = 1.2" in mdp


def test_nvt_mdp_dispcorr(workspace):
    mdp = (workspace / "steps/10_equilibration/nvt.mdp").read_text()
    assert "DispCorr                = EnerPres" in mdp


def test_npt_mdp_continuation(workspace):
    mdp = (workspace / "steps/10_equilibration/npt.mdp").read_text()
    assert "continuation            = yes" in mdp


def test_npt_mdp_berendsen(workspace):
    mdp = (workspace / "steps/10_equilibration/npt.mdp").read_text()
    assert "pcoupl                  = Berendsen" in mdp


def test_npt_mdp_semiisotropic(workspace):
    mdp = (workspace / "steps/10_equilibration/npt.mdp").read_text()
    assert "pcoupltype              = semiisotropic" in mdp


def test_npt_mdp_ref_p_half(workspace):
    mdp = (workspace / "steps/10_equilibration/npt.mdp").read_text()
    assert "ref_p                   = 0.5  0.5" in mdp


def test_production_dt_001(workspace):
    mdp = (workspace / "steps/11_production_md/md.mdp").read_text()
    assert "dt                      = 0.001" in mdp, "OPLS-AA lipid stability requires dt=0.001"


def test_production_nose_hoover(workspace):
    mdp = (workspace / "steps/11_production_md/md.mdp").read_text()
    assert "tcoupl                  = Nose-Hoover" in mdp


def test_production_parrinello_rahman(workspace):
    mdp = (workspace / "steps/11_production_md/md.mdp").read_text()
    assert "pcoupl                  = Parrinello-Rahman" in mdp


def test_production_semiisotropic(workspace):
    mdp = (workspace / "steps/11_production_md/md.mdp").read_text()
    assert "pcoupltype              = semiisotropic" in mdp


def test_production_ref_p_1_bar(workspace):
    mdp = (workspace / "steps/11_production_md/md.mdp").read_text()
    assert "ref_p                   = 1.0  1.0" in mdp


def test_production_dispcorr(workspace):
    mdp = (workspace / "steps/11_production_md/md.mdp").read_text()
    assert "DispCorr                = EnerPres" in mdp


def test_generate_topology_uses_oplsaa_membrane(workspace):
    script = (workspace / "steps/04_generate_topology/run.sh").read_text()
    assert "oplsaa_membrane" in script


# ── Broken case: protein-only YAML should NOT use MembraneWorkflowOPLSAA ─────

PROTEIN_YAML = CONFIGS_DIR / "lysozyme_test.yaml"

def test_protein_only_uses_md_pipeline():
    if not PROTEIN_YAML.exists():
        pytest.skip("lysozyme_test.yaml not found")
    result = SimulationCompiler().compile(str(PROTEIN_YAML))
    assert result.state.inferred_system_type != "protein-membrane"
    notes_combined = " ".join(result.plan.notes)
    assert "MembraneWorkflowOPLSAA" not in notes_combined


# ── Adapter smoke tests ────────────────────────────────────────────────────────

def test_water_deletor_always_available():
    adapter = WaterDeletorAdapter()
    avail = adapter.check_availability()
    assert avail.available
    assert "python3" in avail.binary_path


def test_movememb_always_available():
    adapter = MoveMembAdapter()
    avail = adapter.check_availability()
    assert avail.available
    assert "python3" in avail.binary_path


def _make_gro(tmp_path: Path, atoms: list[tuple]) -> Path:
    """
    Build a minimal .gro file.
    atoms: list of (resnum, resname, atomname, x, y, z)
    """
    p = tmp_path / "test.gro"
    lines = ["test system", f"{len(atoms)}"]
    for i, (rn, rname, aname, x, y, z) in enumerate(atoms, start=1):
        lines.append(f"{rn:5d}{rname:<5}{aname:>5}{i:5d}{x:8.3f}{y:8.3f}{z:8.3f}")
    lines.append("10.00000  10.00000  10.00000")
    p.write_text("\n".join(lines) + "\n")
    return p


def test_water_deletor_removes_bilayer_water(tmp_path):
    """SOL OW inside bilayer [z_bot, z_top] must be removed."""
    atoms = [
        # Lipid top-leaflet headgroup (O33, z=4.0 > midplane=2.0)
        (1, "DPP",  "O33", 0.0, 0.0, 4.0),
        # Lipid bottom-leaflet headgroup (O33, z=0.0 < midplane=2.0)
        (2, "DPP",  "O33", 0.0, 0.0, 0.0),
        # Lipid middle atoms (C50, midplane z ≈ 2.0)
        (1, "DPP",  "C50", 0.0, 0.0, 2.0),
        (2, "DPP",  "C50", 0.0, 0.0, 2.0),
        # SOL inside bilayer: OW z=2.5 ∈ [0.0, 4.0] → MUST be deleted
        (3, "SOL",  "OW",  1.0, 0.0, 2.5),
        (3, "SOL",  "HW1", 1.0, 0.1, 2.5),
        (3, "SOL",  "HW2", 1.0, 0.0, 2.6),
        # SOL outside bilayer: OW z=6.0 > z_top=4.0 → kept
        (4, "SOL",  "OW",  1.0, 0.0, 6.0),
        (4, "SOL",  "HW1", 1.0, 0.1, 6.0),
        (4, "SOL",  "HW2", 1.0, 0.0, 6.1),
    ]
    gro_in  = _make_gro(tmp_path, atoms)
    gro_out = tmp_path / "out.gro"
    adapter = WaterDeletorAdapter()
    result  = adapter.run(gro_in=gro_in, gro_out=gro_out, ref_atom="O33", middle_atom="C50")

    assert result.success
    assert result.metadata["waters_removed"] == 1
    assert result.metadata["atoms_out"] == 7   # 4 lipid + 3 outside SOL
    assert gro_out.exists()


def test_water_deletor_keeps_exterior_water(tmp_path):
    """SOL outside bilayer Z range must not be removed."""
    atoms = [
        (1, "DPP", "O33", 0.0, 0.0, 4.0),
        (2, "DPP", "O33", 0.0, 0.0, 0.0),
        (1, "DPP", "C50", 0.0, 0.0, 2.0),
        (2, "DPP", "C50", 0.0, 0.0, 2.0),
        # OW at z=8.0 — well above z_top=4.0
        (3, "SOL", "OW",  0.0, 0.0, 8.0),
        (3, "SOL", "HW1", 0.0, 0.0, 8.1),
        (3, "SOL", "HW2", 0.0, 0.1, 8.0),
    ]
    gro_in  = _make_gro(tmp_path, atoms)
    gro_out = tmp_path / "out.gro"
    result  = WaterDeletorAdapter().run(
        gro_in=gro_in, gro_out=gro_out, ref_atom="O33", middle_atom="C50",
    )
    assert result.success
    assert result.metadata["waters_removed"] == 0
    assert result.metadata["atoms_out"] == 7


def test_movememb_z_shift_centres_bilayer(tmp_path):
    """Bilayer midplane must be shifted to protein Z-centre."""
    # Protein: Z from 3.0 to 5.0 → centre = 4.0
    prot_atoms = [
        (1, "ALA", "CA", 0.0, 0.0, 3.0),
        (1, "ALA", "C",  0.0, 0.0, 5.0),
    ]
    # Bilayer: Z from 0.0 to 4.0 → midplane = 2.0
    bil_atoms = [
        (1, "DPP", "O33", 0.0, 0.0, 0.0),
        (2, "DPP", "O33", 0.0, 0.0, 4.0),
    ]
    prot_gro = tmp_path / "prot.gro"
    prot_gro.write_text(
        "protein\n2\n"
        "    1ALA   CA    1   0.000   0.000   3.000\n"
        "    1ALA    C    2   0.000   0.000   5.000\n"
        "10.00000  10.00000  10.00000\n"
    )
    bil_gro = tmp_path / "bilayer.gro"
    bil_gro.write_text(
        "bilayer\n2\n"
        "    1DPP  O33    1   0.000   0.000   0.000\n"
        "    2DPP  O33    2   0.000   0.000   4.000\n"
        "10.00000  10.00000  10.00000\n"
    )
    gro_out = tmp_path / "system.gro"
    adapter = MoveMembAdapter()
    result  = adapter.run(protein_gro=prot_gro, bilayer_gro=bil_gro, gro_out=gro_out)

    assert result.success
    # Expected shift: prot_centre(4.0) - bil_midplane(2.0) = +2.0 nm
    assert abs(result.metadata["z_shift_nm"] - 2.0) < 1e-6
    assert gro_out.exists()
    # Combined file has 4 atoms (2 protein + 2 bilayer)
    _, atom_lines, _ = _parse_gro(gro_out)
    assert len(atom_lines) == 4
    # Bilayer atoms should now be at z=2.0 and z=6.0
    zs = sorted(float(l[36:44]) for l in atom_lines)
    assert abs(zs[0] - 2.0) < 1e-3   # protein CA at z=3.0 and shifted bil at z=2.0
