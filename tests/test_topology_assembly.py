# tests/test_topology_assembly.py
"""
Phase 11 — Topology assembly tests.

Tests:
  A. TopologyBuildError / is_mixed_system_gro
  B. check_pdb2gmx_input_safety guardrail
  C. generate_protein_topology builder produces correct workspace
  D. generate_topology (refactored) builder produces correct workspace
  E. assemble_system_topology builder produces correct workspace
  F. Pipeline step ordering (generate_protein_topology before orient_protein)
  G. Pipeline generate_topology engine changed from pdb2gmx to topology:assemble_embed
  H. solvate_membrane topology path updated to assemble_system_topology
  I. generate_topology run_topology.py has no _find_root / sys.path.insert
  J. embedding_builder still resolves system_processed.gro from generate_topology
"""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest


# ═══════════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════════

def _write_gro(path: Path, title: str, atom_lines: list[str], box: str = "10.0 10.0 10.0") -> None:
    lines = [title, f"{len(atom_lines):5d}"] + atom_lines + [box]
    path.write_text("\n".join(lines) + "\n")


def _make_mixed_gro(path: Path) -> None:
    """GRO with protein (ALA) + lipid (DPP) residues."""
    atoms = [
        "    1ALA      N    1   1.000   1.000   1.000",
        "    1ALA     CA    2   1.100   1.000   1.000",
        "    2DPP     P1    3   2.000   2.000   2.000",
        "    2DPP     O1    4   2.100   2.000   2.000",
    ]
    _write_gro(path, "Mixed system", atoms)


def _make_protein_only_gro(path: Path) -> None:
    """GRO with only protein residues."""
    atoms = [
        "    1ALA      N    1   1.000   1.000   1.000",
        "    1ALA     CA    2   1.100   1.000   1.000",
        "    2GLY      N    3   1.500   1.000   1.000",
    ]
    _write_gro(path, "Protein only", atoms)


def _make_lipid_only_gro(path: Path) -> None:
    """GRO with only lipid residues."""
    atoms = [
        "    1DPP     P1    1   2.000   2.000   2.000",
        "    1DPP     O1    2   2.100   2.000   2.000",
    ]
    _write_gro(path, "Lipid only", atoms)


def _build_pipeline_workspace(tmp_path: Path) -> Path:
    """Compile and build a membrane pipeline workspace."""
    from core.compiler import SimulationCompiler
    from builders.workspace_builder import WorkspaceBuilder
    import os

    yaml_path = Path("configs/membrane_orient_test.yaml").resolve()
    orig_cwd = os.getcwd()
    os.chdir(yaml_path.parent)
    try:
        result = SimulationCompiler().compile(str(yaml_path))
    finally:
        os.chdir(orig_cwd)

    return WorkspaceBuilder().build(result, workspace_path=tmp_path)


# ═══════════════════════════════════════════════════════════════════════════════
# A. TopologyBuildError / is_mixed_system_gro
# ═══════════════════════════════════════════════════════════════════════════════

class TestIsMixedSystemGro:

    def test_mixed_gro_detected(self, tmp_path):
        gro = tmp_path / "mixed.gro"
        _make_mixed_gro(gro)
        from core.topology_models import is_mixed_system_gro
        assert is_mixed_system_gro(gro) is True

    def test_protein_only_is_not_mixed(self, tmp_path):
        gro = tmp_path / "protein.gro"
        _make_protein_only_gro(gro)
        from core.topology_models import is_mixed_system_gro
        assert is_mixed_system_gro(gro) is False

    def test_lipid_only_is_not_mixed(self, tmp_path):
        gro = tmp_path / "lipid.gro"
        _make_lipid_only_gro(gro)
        from core.topology_models import is_mixed_system_gro
        assert is_mixed_system_gro(gro) is False

    def test_nonexistent_file_returns_false(self, tmp_path):
        from core.topology_models import is_mixed_system_gro
        assert is_mixed_system_gro(tmp_path / "nonexistent.gro") is False

    def test_solvent_triggers_mixed(self, tmp_path):
        gro = tmp_path / "prot_sol.gro"
        atoms = [
            "    1ALA      N    1   1.000   1.000   1.000",
            "    2SOL     OW    2   5.000   5.000   5.000",
        ]
        _write_gro(gro, "protein + water", atoms)
        from core.topology_models import is_mixed_system_gro
        assert is_mixed_system_gro(gro) is True


# ═══════════════════════════════════════════════════════════════════════════════
# B. check_pdb2gmx_input_safety guardrail
# ═══════════════════════════════════════════════════════════════════════════════

class TestTopologyGuardrail:

    def test_embed_in_bilayer_path_blocked(self, tmp_path):
        from validators.topology_guardrail import check_pdb2gmx_input_safety
        from core.topology_models import TopologyBuildError
        bad_path = tmp_path / "steps/03_embed_in_bilayer/system.gro"
        with pytest.raises(TopologyBuildError, match="embed_in_bilayer"):
            check_pdb2gmx_input_safety(bad_path)

    def test_membrane_embedding_path_blocked(self, tmp_path):
        from validators.topology_guardrail import check_pdb2gmx_input_safety
        from core.topology_models import TopologyBuildError
        bad_path = tmp_path / "steps/05_membrane_embedding/converged.gro"
        with pytest.raises(TopologyBuildError, match="membrane_embedding"):
            check_pdb2gmx_input_safety(bad_path)

    def test_converged_gro_path_blocked(self, tmp_path):
        from validators.topology_guardrail import check_pdb2gmx_input_safety
        from core.topology_models import TopologyBuildError
        bad_path = tmp_path / "some_dir/converged.gro"
        with pytest.raises(TopologyBuildError, match="converged"):
            check_pdb2gmx_input_safety(bad_path)

    def test_mixed_gro_blocked_by_content(self, tmp_path):
        from validators.topology_guardrail import check_pdb2gmx_input_safety
        from core.topology_models import TopologyBuildError
        gro = tmp_path / "mixed_system.gro"
        _make_mixed_gro(gro)
        with pytest.raises(TopologyBuildError, match="mixed"):
            check_pdb2gmx_input_safety(gro)

    def test_pdb_file_passes_guardrail(self, tmp_path):
        from validators.topology_guardrail import check_pdb2gmx_input_safety
        pdb = tmp_path / "protein_1.pdb"
        pdb.write_text("ATOM      1  N   ALA A   1       1.000   1.000   1.000\nEND\n")
        check_pdb2gmx_input_safety(pdb)  # must not raise

    def test_protein_only_gro_passes_guardrail(self, tmp_path):
        from validators.topology_guardrail import check_pdb2gmx_input_safety
        gro = tmp_path / "protein.gro"
        _make_protein_only_gro(gro)
        check_pdb2gmx_input_safety(gro)  # must not raise

    def test_lipid_only_gro_passes_guardrail(self, tmp_path):
        from validators.topology_guardrail import check_pdb2gmx_input_safety
        gro = tmp_path / "lipids_only.gro"
        _make_lipid_only_gro(gro)
        check_pdb2gmx_input_safety(gro)  # must not raise

    def test_nonexistent_gro_passes_guardrail(self, tmp_path):
        from validators.topology_guardrail import check_pdb2gmx_input_safety
        # Non-existent GRO not in embed path → passes (can't read content to check)
        check_pdb2gmx_input_safety(tmp_path / "will_be_created.gro")


# ═══════════════════════════════════════════════════════════════════════════════
# C. generate_protein_topology builder
# ═══════════════════════════════════════════════════════════════════════════════

class TestGenerateProteinTopologyBuilder:

    @pytest.fixture(scope="class")
    def step_dir(self, tmp_path_factory):
        from core.execution_models import SimulationStep, StepStage, StepType, AutomationLevel
        from builders.step_builders.preparation_builder import PreparationBuilder

        step = SimulationStep(
            step_id="generate_protein_topology",
            title="Generar topología proteína",
            stage=StepStage.PREPARATION,
            step_type=StepType.AUTOMATIC,
            automation_level=AutomationLevel.AUTOMATED,
            engine="gromacs:pdb2gmx_protein",
            target_components=["protein_1"],
            params={
                "source_file": "protein_1.pdb",
                "forcefield":  "opls-aa-membrane",
                "water_model": "none",
            },
        )
        d = tmp_path_factory.mktemp("gen_prot_top")
        PreparationBuilder().build(step, d, step_dir_map={})
        return d

    def test_run_protein_topology_py_created(self, step_dir):
        assert (step_dir / "run_protein_topology.py").exists()

    def test_run_sh_created(self, step_dir):
        assert (step_dir / "run.sh").exists()

    def test_run_sh_calls_run_protein_topology(self, step_dir):
        content = (step_dir / "run.sh").read_text()
        assert "run_protein_topology.py" in content

    def test_script_has_pdb2gmx(self, step_dir):
        content = (step_dir / "run_protein_topology.py").read_text()
        assert "pdb2gmx" in content

    def test_script_references_pdb_input(self, step_dir):
        content = (step_dir / "run_protein_topology.py").read_text()
        assert "protein_1.pdb" in content

    def test_script_no_find_root(self, step_dir):
        content = (step_dir / "run_protein_topology.py").read_text()
        assert "_find_root" not in content
        assert "sys.path.insert" not in content

    def test_script_uses_oplsaa_membrane_ff(self, step_dir):
        content = (step_dir / "run_protein_topology.py").read_text()
        assert "oplsaa_membrane" in content

    def test_metadata_correct(self, step_dir):
        meta = json.loads((step_dir / "metadata.json").read_text())
        assert meta["step_id"] == "generate_protein_topology"
        assert meta["engine"] == "gromacs:pdb2gmx_protein"
        assert "protein_processed.gro" in meta["expected_outputs"]
        assert "topol.top" in meta["expected_outputs"]
        assert "protein_topology_manifest.json" in meta["expected_outputs"]

    def test_script_writes_manifest_json(self, step_dir):
        content = (step_dir / "run_protein_topology.py").read_text()
        assert "protein_topology_manifest.json" in content

    def test_script_sets_gmxlib(self, step_dir):
        content = (step_dir / "run_protein_topology.py").read_text()
        assert "GMXLIB" in content


# ═══════════════════════════════════════════════════════════════════════════════
# D. generate_topology (refactored) builder
# ═══════════════════════════════════════════════════════════════════════════════

class TestGenerateTopologyRefactoredBuilder:

    @pytest.fixture(scope="class")
    def step_dir(self, tmp_path_factory):
        from core.execution_models import SimulationStep, StepStage, StepType, AutomationLevel
        from builders.step_builders.preparation_builder import PreparationBuilder

        step = SimulationStep(
            step_id="generate_topology",
            title="Ensamblar topología de embedding",
            stage=StepStage.PREPARATION,
            step_type=StepType.AUTOMATIC,
            automation_level=AutomationLevel.AUTOMATED,
            engine="topology:assemble_embed",
            depends_on=["embed_in_bilayer", "generate_protein_topology"],
            params={
                "source_step":  "embed_in_bilayer",
                "forcefield":   "opls-aa-membrane",
                "water_model":  "none",
            },
        )
        d = tmp_path_factory.mktemp("gen_top_refactored")
        PreparationBuilder().build(step, d, step_dir_map={})
        return d

    def test_run_topology_py_created(self, step_dir):
        assert (step_dir / "run_topology.py").exists()

    def test_run_sh_created(self, step_dir):
        assert (step_dir / "run.sh").exists()

    def test_run_sh_calls_run_topology(self, step_dir):
        content = (step_dir / "run.sh").read_text()
        assert "run_topology.py" in content

    def test_script_no_find_root(self, step_dir):
        content = (step_dir / "run_topology.py").read_text()
        assert "_find_root" not in content
        assert "sys.path.insert" not in content

    def test_script_does_not_run_pdb2gmx_on_system_gro(self, step_dir):
        """Core guardrail: the script must not invoke pdb2gmx on system.gro."""
        content = (step_dir / "run_topology.py").read_text()
        # Must not have pdb2gmx with embed_in_bilayer/system.gro as input
        assert '"embed_in_bilayer/system.gro"' not in content or "pdb2gmx" not in content.split('"embed_in_bilayer/system.gro"')[0].rsplit("pdb2gmx", 1)[-1]

    def test_script_reads_protein_topology(self, step_dir):
        content = (step_dir / "run_topology.py").read_text()
        assert "generate_protein_topology" in content
        assert "topol.top" in content

    def test_script_copies_embed_system_gro(self, step_dir):
        content = (step_dir / "run_topology.py").read_text()
        assert "system_processed.gro" in content
        assert "system.gro" in content

    def test_script_writes_topology_consistency_report(self, step_dir):
        content = (step_dir / "run_topology.py").read_text()
        assert "topology_consistency_report.json" in content

    def test_script_handles_strong_posre(self, step_dir):
        content = (step_dir / "run_topology.py").read_text()
        assert "strong_posre.itp" in content
        assert "STRONG_POSRES" in content

    def test_script_sets_gmxlib(self, step_dir):
        content = (step_dir / "run_topology.py").read_text()
        assert "GMXLIB" in content

    def test_metadata_correct(self, step_dir):
        meta = json.loads((step_dir / "metadata.json").read_text())
        assert meta["step_id"] == "generate_topology"
        assert meta["engine"] == "topology:assemble_embed"
        assert "system_processed.gro" in meta["expected_outputs"]
        assert "topol.top" in meta["expected_outputs"]
        assert "topology_consistency_report.json" in meta["expected_outputs"]
        assert meta["gate"]["type"] == "topology_consistency"


# ═══════════════════════════════════════════════════════════════════════════════
# E. assemble_system_topology builder
# ═══════════════════════════════════════════════════════════════════════════════

class TestAssembleSystemTopologyBuilder:

    @pytest.fixture(scope="class")
    def step_dir(self, tmp_path_factory):
        from core.execution_models import SimulationStep, StepStage, StepType, AutomationLevel
        from builders.step_builders.assembly_builder import AssemblyBuilder

        step = SimulationStep(
            step_id="assemble_system_topology",
            title="Ensamblar topología final",
            stage=StepStage.ASSEMBLY,
            step_type=StepType.AUTOMATIC,
            automation_level=AutomationLevel.AUTOMATED,
            engine="topology:assemble_system",
            depends_on=["membrane_embedding"],
            params={
                "forcefield":  "opls-aa-membrane",
                "water_model": "none",
            },
        )
        d = tmp_path_factory.mktemp("assemble_sys_top")
        AssemblyBuilder().build(step, d, step_dir_map={})
        return d

    def test_run_assemble_system_py_created(self, step_dir):
        assert (step_dir / "run_assemble_system.py").exists()

    def test_run_sh_created(self, step_dir):
        assert (step_dir / "run.sh").exists()

    def test_run_sh_calls_run_assemble_system(self, step_dir):
        content = (step_dir / "run.sh").read_text()
        assert "run_assemble_system.py" in content

    def test_script_reads_embed_in_bilayer_system_gro(self, step_dir):
        content = (step_dir / "run_assemble_system.py").read_text()
        assert "system.gro" in content
        assert "embed_in_bilayer" in content or "EMBED_DIR" in content

    def test_script_reads_protein_topology(self, step_dir):
        content = (step_dir / "run_assemble_system.py").read_text()
        assert "generate_protein_topology" in content
        assert "topol.top" in content

    def test_script_writes_topol_top(self, step_dir):
        content = (step_dir / "run_assemble_system.py").read_text()
        assert "topol.top" in content

    def test_script_writes_assembly_report(self, step_dir):
        content = (step_dir / "run_assemble_system.py").read_text()
        assert "topology_assembly_report.json" in content

    def test_script_no_find_root(self, step_dir):
        content = (step_dir / "run_assemble_system.py").read_text()
        assert "_find_root" not in content
        assert "sys.path.insert" not in content

    def test_metadata_correct(self, step_dir):
        meta = json.loads((step_dir / "metadata.json").read_text())
        assert meta["step_id"] == "assemble_system_topology"
        assert "topol.top" in meta["expected_outputs"]
        assert "topology_assembly_report.json" in meta["expected_outputs"]


# ═══════════════════════════════════════════════════════════════════════════════
# F–J. Pipeline integration tests (no GROMACS required)
# ═══════════════════════════════════════════════════════════════════════════════

class TestPhase11PipelineIntegration:

    @pytest.fixture(scope="class")
    def workspace(self, tmp_path_factory):
        ws = tmp_path_factory.mktemp("ws_phase11")
        return _build_pipeline_workspace(ws)

    def _step_dir(self, workspace: Path, name: str) -> Path:
        return next(
            d for d in (workspace / "steps").iterdir()
            if name in d.name
        )

    # F. Step ordering
    def test_generate_protein_topology_step_exists(self, workspace):
        step_names = [d.name for d in (workspace / "steps").iterdir()]
        assert any("generate_protein_topology" in n for n in step_names), (
            "generate_protein_topology step dir must exist"
        )

    def test_generate_protein_topology_before_orient_protein(self, workspace):
        """generate_protein_topology (step 0) must come before orient_protein."""
        step_dirs = sorted(d.name for d in (workspace / "steps").iterdir())
        gen_prot_idx = next(i for i, n in enumerate(step_dirs) if "generate_protein_topology" in n)
        orient_idx   = next(i for i, n in enumerate(step_dirs) if "orient_protein" in n)
        assert gen_prot_idx < orient_idx, (
            "generate_protein_topology must be ordered before orient_protein"
        )

    def test_assemble_system_topology_step_exists(self, workspace):
        step_names = [d.name for d in (workspace / "steps").iterdir()]
        assert any("assemble_system_topology" in n for n in step_names), (
            "assemble_system_topology step dir must exist"
        )

    def test_assemble_system_topology_after_membrane_embedding(self, workspace):
        """assemble_system_topology runs AFTER membrane_embedding (post-shrink DAG)."""
        step_dirs = sorted(d.name for d in (workspace / "steps").iterdir())
        assemble_idx = next(i for i, n in enumerate(step_dirs) if "assemble_system_topology" in n)
        memb_emb_idx = next(i for i, n in enumerate(step_dirs) if "membrane_embedding" in n)
        assert assemble_idx > memb_emb_idx, (
            "assemble_system_topology must come AFTER membrane_embedding "
            "(final topology from converged.gro)"
        )

    def test_no_generate_topology_step_in_dag(self, workspace):
        """generate_topology step is removed; assemble_system_topology is the single assembler."""
        step_dirs = [d.name for d in (workspace / "steps").iterdir()]
        has_gen_topo = any("generate_topology" in n and "generate_protein_topology" not in n for n in step_dirs)
        assert not has_gen_topo, (
            "generate_topology step must not exist; assemble_system_topology is the single assembler"
        )

    # G. assemble_system_topology engine is correct
    def test_assemble_system_topology_engine(self, workspace):
        topo_dir = self._step_dir(workspace, "assemble_system_topology")
        meta = json.loads((topo_dir / "metadata.json").read_text())
        assert meta["engine"] == "topology:assemble_system", (
            "assemble_system_topology must use engine='topology:assemble_system'"
        )

    def test_assemble_system_topology_no_pdb2gmx_on_system_gro(self, workspace):
        """run_assemble_system.py must not pass system.gro as -f to pdb2gmx (only lipids_only.gro is ok)."""
        topo_dir = self._step_dir(workspace, "assemble_system_topology")
        script = (topo_dir / "run_assemble_system.py").read_text()
        import re as _re
        f_args = _re.findall(r'"-f",\s+(?:str\()?"([^"]+)"', script)
        for f_arg in f_args:
            assert "system.gro" not in f_arg, (
                f"pdb2gmx -f argument must not be system.gro, got: {f_arg!r}"
            )

    # H. solvate_membrane topology path updated
    def test_solvate_membrane_uses_assemble_system_topology(self, workspace):
        solvate_dir = self._step_dir(workspace, "solvate_membrane")
        scripts = list(solvate_dir.glob("*.sh")) + list(solvate_dir.glob("*.py"))
        combined = " ".join(s.read_text() for s in scripts if s.exists())
        assert "assemble_system_topology" in combined, (
            "solvate_membrane must reference assemble_system_topology/topol.top"
        )

    # I. assemble_system_topology script quality
    def test_assemble_system_topology_run_script_no_find_root(self, workspace):
        topo_dir = self._step_dir(workspace, "assemble_system_topology")
        script = (topo_dir / "run_assemble_system.py").read_text()
        assert "_find_root" not in script
        assert "sys.path.insert" not in script

    def test_assemble_system_topology_produces_system_processed_gro(self, workspace):
        topo_dir = self._step_dir(workspace, "assemble_system_topology")
        script = (topo_dir / "run_assemble_system.py").read_text()
        assert "system_processed.gro" in script, (
            "run_assemble_system.py must write system_processed.gro for membrane_embedding"
        )

    def test_assemble_system_topology_references_protein_topology(self, workspace):
        topo_dir = self._step_dir(workspace, "assemble_system_topology")
        script = (topo_dir / "run_assemble_system.py").read_text()
        assert "generate_protein_topology" in script or "PROT_TOP_DIR" in script, (
            "run_assemble_system.py must read protein topology from generate_protein_topology"
        )

    # J. Embedding builder bootstrap-calls run_assemble_system.py (not system_processed.gro)
    def test_embedding_bootstrap_calls_run_assemble_system(self, workspace):
        """membrane_embedding/run.sh must call assemble_system_topology/run_assemble_system.py
        for bootstrap topology — NOT start from system_processed.gro."""
        embed_dir = self._step_dir(workspace, "membrane_embedding")
        script = (embed_dir / "run.sh").read_text()
        assert "assemble_system_topology" in script, (
            "membrane_embedding/run.sh must reference assemble_system_topology (for bootstrap topology)"
        )
        assert "run_assemble_system.py" in script, (
            "membrane_embedding/run.sh must call run_assemble_system.py for bootstrap topology"
        )
        assert "embed_in_bilayer" in script, (
            "membrane_embedding/run.sh must reference embed_in_bilayer/system.gro as input"
        )
        assert "system_processed.gro" not in script, (
            "membrane_embedding/run.sh must NOT use system_processed.gro as input "
            "(assemble_system_topology runs after embedding in the new DAG)"
        )

    def test_generate_protein_topology_script_exists_in_workspace(self, workspace):
        gen_prot_dir = self._step_dir(workspace, "generate_protein_topology")
        assert (gen_prot_dir / "run_protein_topology.py").exists()

    def test_generate_protein_topology_uses_inputs_pdb(self, workspace):
        import re as _re
        gen_prot_dir = self._step_dir(workspace, "generate_protein_topology")
        script = (gen_prot_dir / "run_protein_topology.py").read_text()
        assert ".pdb" in script, "run_protein_topology.py must use a PDB input file"
        f_args = _re.findall(r'"-f",\s+(?:str\()?"([^"]+)"', script)
        for arg in f_args:
            assert "system.gro" not in arg, (
                f"run_protein_topology.py -f argument must not be system.gro, got: {arg}"
            )


# ═══════════════════════════════════════════════════════════════════════════════
# K. Input materialization regression
#    Source PDB in a subdirectory (e.g. configs/Canal.pdb) must be materialized
#    into <workspace>/inputs/protein_1.pdb — and the generate_protein_topology
#    step must reference THAT artifact, not the original source path.
# ═══════════════════════════════════════════════════════════════════════════════

class TestInputMaterialization:
    """
    Regression for the bug where generate_protein_topology used os.path.relpath()
    from the step directory back to the *original* source PDB absolute path, producing
    paths like '../../inputs/configs/Canal.pdb' that don't exist in the workspace.

    The fix: always derive pdb_ref from inputs/{comp_id}{ext} — the canonical
    workspace-materialized artifact — regardless of where the source lives.
    """

    _YAML_TEMPLATE = """\
project:
  name: test_materialization
components:
  - id: protein_1
    role: protein
    file: {pdb_relpath}
structural_annotation:
  membrane_topology:
    extracellular_regions: ["1-10"]
    intracellular_regions: ["40-50"]
    transmembrane_segments:
      - residues: "11-39"
        label: "TM1"
environment:
  membrane:
    enabled: true
forcefields:
  protein: opls-aa
simulation_objectives:
  - membrane_insertion
"""

    def _build(self, tmp_path: Path, pdb_subdir: str = "configs"):
        """Build workspace from a YAML whose protein PDB is in a subdirectory."""
        import os
        from core.compiler import SimulationCompiler
        from builders.workspace_builder import WorkspaceBuilder

        # Put PDB in a subdirectory to simulate configs/Canal.pdb scenario
        pdb_dir = tmp_path / pdb_subdir
        pdb_dir.mkdir(parents=True, exist_ok=True)
        pdb_file = pdb_dir / "Canal.pdb"
        pdb_file.write_text(
            "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00\n"
            "TER\n"
            "END\n"
        )

        yaml_text = self._YAML_TEMPLATE.format(pdb_relpath=f"{pdb_subdir}/Canal.pdb")
        yaml_file = tmp_path / "test.yaml"
        yaml_file.write_text(yaml_text)

        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            state = SimulationCompiler().compile(str(yaml_file))
        finally:
            os.chdir(orig)

        ws = tmp_path / "workspace"
        WorkspaceBuilder().build(state, workspace_path=ws)
        return ws, pdb_file

    # ── materialization ────────────────────────────────────────────────────────

    def test_protein_pdb_materialized_into_inputs(self, tmp_path):
        ws, _ = self._build(tmp_path)
        materialized = ws / "inputs" / "protein_1.pdb"
        assert materialized.exists(), (
            "WorkspaceBuilder must copy the protein PDB to inputs/protein_1.pdb"
        )

    def test_materialized_pdb_has_ter_records(self, tmp_path):
        ws, _ = self._build(tmp_path)
        content = (ws / "inputs" / "protein_1.pdb").read_text()
        assert "TER" in content, "Materialized PDB must preserve TER records"

    def test_original_source_not_in_script(self, tmp_path):
        ws, pdb_file = self._build(tmp_path)
        gen_prot_dir = next((ws / "steps").glob("*_generate_protein_topology"))
        script = (gen_prot_dir / "run_protein_topology.py").read_text()
        # The step must not reference the original source path (configs/Canal.pdb)
        assert str(pdb_file) not in script, (
            "run_protein_topology.py must not contain the original absolute source path"
        )
        assert "configs/Canal.pdb" not in script, (
            "run_protein_topology.py must not reference configs/Canal.pdb"
        )

    def test_script_references_workspace_inputs(self, tmp_path):
        import re as _re
        ws, _ = self._build(tmp_path)
        gen_prot_dir = next((ws / "steps").glob("*_generate_protein_topology"))
        script = (gen_prot_dir / "run_protein_topology.py").read_text()
        # The -f argument must point inside inputs/
        f_args = _re.findall(r'"-f",\s+str\(PDB_INPUT\)', script)
        # PDB_INPUT is defined via a pdb_ref containing "inputs"
        assert "inputs" in script, (
            "run_protein_topology.py must reference the workspace inputs/ directory"
        )
        assert "protein_1.pdb" in script, (
            "run_protein_topology.py must reference the materialized protein_1.pdb"
        )

    def test_required_inputs_in_metadata_uses_workspace_path(self, tmp_path):
        ws, _ = self._build(tmp_path)
        gen_prot_dir = next((ws / "steps").glob("*_generate_protein_topology"))
        meta = json.loads((gen_prot_dir / "metadata.json").read_text())
        required = meta.get("required_inputs", [])
        assert required, "metadata.json must declare required_inputs"
        for p in required:
            assert "configs" not in p, (
                f"required_inputs must not reference the config subdirectory, got: {p}"
            )
            assert "Canal.pdb" not in p, (
                f"required_inputs must not reference the original filename Canal.pdb, got: {p}"
            )
            assert "protein_1.pdb" in p or p.endswith(".pdb"), (
                f"required_inputs must reference the materialized PDB artifact, got: {p}"
            )

    def test_preflight_passes_after_build(self, tmp_path):
        ws, _ = self._build(tmp_path)
        gen_prot_dir = next((ws / "steps").glob("*_generate_protein_topology"))
        meta = json.loads((gen_prot_dir / "metadata.json").read_text())
        required = meta.get("required_inputs", [])
        missing = [
            rel for rel in required
            if not (gen_prot_dir / rel).resolve().exists()
        ]
        assert not missing, (
            f"Preflight check failed — these required_inputs do not exist:\n"
            + "\n".join(f"  {r}" for r in missing)
        )

    def test_no_pdb2gmx_on_system_gro_in_generate_protein_topology(self, tmp_path):
        import re as _re
        ws, _ = self._build(tmp_path)
        gen_prot_dir = next((ws / "steps").glob("*_generate_protein_topology"))
        script = (gen_prot_dir / "run_protein_topology.py").read_text()
        f_args = _re.findall(r'"-f",\s+(?:str\([^)]+\)|"[^"]+")', script)
        for arg in f_args:
            assert "system.gro" not in arg, (
                f"generate_protein_topology must never call pdb2gmx on system.gro; got: {arg}"
            )
            assert "embed" not in arg.lower(), (
                f"generate_protein_topology -f argument must not reference an embed step; got: {arg}"
            )

    def test_generate_topology_no_pdb2gmx_on_system_gro(self, tmp_path):
        import re as _re
        ws, _ = self._build(tmp_path)
        gen_topo_dir = next((ws / "steps").glob("*_generate_topology"), None)
        if gen_topo_dir is None:
            return  # step not present in this pipeline variant
        run_py = gen_topo_dir / "run_topology.py"
        if not run_py.exists():
            return
        script = run_py.read_text()
        # pdb2gmx -f system.gro is the bug we're guarding against
        f_args = _re.findall(r'"-f",\s+(?:str\([^)]+\)|"[^"]+")', script)
        for arg in f_args:
            assert "system.gro" not in arg or "lipids_only" in arg or "lipid" in arg.lower(), (
                f"generate_topology must not call pdb2gmx on system.gro directly; got: {arg}"
            )
        # The direct system.gro -f pattern must not appear
        assert not _re.search(r'"-f",\s+.*system\.gro', script), (
            "generate_topology run_topology.py must never pass system.gro directly to pdb2gmx -f"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# L. Multichain pdb2gmx output — dynamic include graph, no posre.itp assumption
# ═══════════════════════════════════════════════════════════════════════════════

class TestMultichainTopology:
    """
    Regression for the bug where SimForge marked generate_protein_topology as
    failed when pdb2gmx produced chain-specific files (posre_Protein_chain_A.itp
    etc.) instead of a literal posre.itp.

    Core invariants:
    - posre.itp must NOT appear in expected_outputs for generate_protein_topology.
    - posre.itp must NOT appear in expected_outputs for generate_topology.
    - The _walk_includes function correctly discovers local include graph.
    - protein_topology_manifest.json records discovered files, not hardcoded names.
    - assemble_system_topology reads the manifest and uses relative paths.
    - ARTIFACT_VALIDATION_ERROR exists and is non-retryable.
    """

    # ── expected_outputs ────────────────────────────────────────────────────────

    def test_generate_protein_topology_no_posre_in_expected_outputs(self, tmp_path):
        from core.compiler import SimulationCompiler
        from builders.workspace_builder import WorkspaceBuilder
        import os
        yaml_text = """\
project:
  name: test_no_posre
components:
  - id: protein_1
    role: protein
    file: protein.pdb
environment:
  membrane:
    enabled: true
forcefields:
  protein: opls-aa
simulation_objectives:
  - membrane_insertion
"""
        (tmp_path / "protein.pdb").write_text(
            "ATOM      1  CA  ALA A   1       1.0   2.0   3.0  1.00  0.00\nTER\nEND\n"
        )
        (tmp_path / "test.yaml").write_text(yaml_text)
        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            state = SimulationCompiler().compile(str(tmp_path / "test.yaml"))
        finally:
            os.chdir(orig)
        ws = tmp_path / "ws"
        WorkspaceBuilder().build(state, workspace_path=ws)

        gen_prot_dir = next((ws / "steps").glob("*_generate_protein_topology"))
        meta = json.loads((gen_prot_dir / "metadata.json").read_text())
        assert "posre.itp" not in meta.get("expected_outputs", []), (
            "posre.itp must NOT be in generate_protein_topology expected_outputs — "
            "multichain pdb2gmx generates posre_Protein_chain_X.itp instead"
        )

    def test_generate_topology_no_posre_in_expected_outputs(self, tmp_path):
        from core.compiler import SimulationCompiler
        from builders.workspace_builder import WorkspaceBuilder
        import os
        yaml_text = """\
project:
  name: test_no_posre_gen_topo
components:
  - id: protein_1
    role: protein
    file: protein.pdb
environment:
  membrane:
    enabled: true
forcefields:
  protein: opls-aa
simulation_objectives:
  - membrane_insertion
"""
        (tmp_path / "protein.pdb").write_text(
            "ATOM      1  CA  ALA A   1       1.0   2.0   3.0  1.00  0.00\nTER\nEND\n"
        )
        (tmp_path / "test.yaml").write_text(yaml_text)
        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            state = SimulationCompiler().compile(str(tmp_path / "test.yaml"))
        finally:
            os.chdir(orig)
        ws = tmp_path / "ws"
        WorkspaceBuilder().build(state, workspace_path=ws)

        gen_topo_dir = next((ws / "steps").glob("*_generate_topology"), None)
        if gen_topo_dir is None:
            return
        meta = json.loads((gen_topo_dir / "metadata.json").read_text())
        assert "posre.itp" not in meta.get("expected_outputs", []), (
            "posre.itp must NOT be in generate_topology expected_outputs"
        )

    # ── _walk_includes ────────────────────────────────────────────────────────

    def test_walk_includes_discovers_chain_itps(self, tmp_path):
        """_walk_includes finds all local .itp files recursively."""
        # Simulate multichain pdb2gmx output
        top = tmp_path / "topol.top"
        top.write_text(
            '#include "oplsaa_membrane.ff/forcefield.itp"\n'
            '#include "topol_Protein_chain_A.itp"\n'
            '#include "topol_Protein_chain_B.itp"\n'
            '[ system ]\n'
            'Protein\n'
            '[ molecules ]\n'
            'Protein_chain_A 1\n'
            'Protein_chain_B 1\n'
        )
        (tmp_path / "topol_Protein_chain_A.itp").write_text(
            '[ moleculetype ]\nProtein_chain_A\n'
            '#ifdef POSRES\n#include "posre_Protein_chain_A.itp"\n#endif\n'
        )
        (tmp_path / "topol_Protein_chain_B.itp").write_text(
            '[ moleculetype ]\nProtein_chain_B\n'
            '#ifdef POSRES\n#include "posre_Protein_chain_B.itp"\n#endif\n'
        )
        (tmp_path / "posre_Protein_chain_A.itp").write_text("[ position_restraints ]\n")
        (tmp_path / "posre_Protein_chain_B.itp").write_text("[ position_restraints ]\n")

        # Import and run _walk_includes from the generated code (as module-level constant)
        from builders.step_builders.preparation_builder import _ITP_WALK_CODE
        ns = {}
        exec("from pathlib import Path\n" + _ITP_WALK_CODE, ns)
        _walk_includes = ns["_walk_includes"]

        local_itps, ff_includes = _walk_includes(top, tmp_path)

        assert "topol_Protein_chain_A.itp" in local_itps
        assert "topol_Protein_chain_B.itp" in local_itps
        assert "posre_Protein_chain_A.itp" in local_itps
        assert "posre_Protein_chain_B.itp" in local_itps
        assert "oplsaa_membrane.ff/forcefield.itp" in ff_includes

    def test_walk_includes_classifies_posre_vs_molecule(self, tmp_path):
        (tmp_path / "topol.top").write_text(
            '#include "ff/forcefield.itp"\n'
            '#include "topol_Protein_chain_A.itp"\n'
            '#include "topol_Protein_chain_B.itp"\n'
        )
        (tmp_path / "topol_Protein_chain_A.itp").write_text(
            '[ moleculetype ]\n#include "posre_Protein_chain_A.itp"\n'
        )
        (tmp_path / "topol_Protein_chain_B.itp").write_text(
            '[ moleculetype ]\n#include "posre_Protein_chain_B.itp"\n'
        )
        (tmp_path / "posre_Protein_chain_A.itp").write_text("")
        (tmp_path / "posre_Protein_chain_B.itp").write_text("")

        from builders.step_builders.preparation_builder import _ITP_WALK_CODE
        ns = {}
        exec("from pathlib import Path\n" + _ITP_WALK_CODE, ns)
        _walk = ns["_walk_includes"]
        local_itps, _ = _walk(tmp_path / "topol.top", tmp_path)

        molecule_itps = [f for f in local_itps if "posre" not in f.lower()]
        posre_itps    = [f for f in local_itps if "posre"     in f.lower()]
        assert len(molecule_itps) == 2
        assert len(posre_itps) == 2

    def test_walk_includes_single_chain_posre_itp(self, tmp_path):
        """Single-chain proteins have a bare posre.itp (still discovered)."""
        (tmp_path / "topol.top").write_text(
            '#include "ff/forcefield.itp"\n'
            '#include "topol_protein.itp"\n'
        )
        (tmp_path / "topol_protein.itp").write_text(
            '[ moleculetype ]\n#include "posre.itp"\n'
        )
        (tmp_path / "posre.itp").write_text("[ position_restraints ]\n")

        from builders.step_builders.preparation_builder import _ITP_WALK_CODE
        ns = {}
        exec("from pathlib import Path\n" + _ITP_WALK_CODE, ns)
        _walk = ns["_walk_includes"]
        local_itps, ff_includes = _walk(tmp_path / "topol.top", tmp_path)

        assert "topol_protein.itp" in local_itps
        assert "posre.itp" in local_itps
        posre_itps = [f for f in local_itps if "posre" in f.lower()]
        assert len(posre_itps) == 1

    def test_walk_includes_missing_file_goes_to_ff(self, tmp_path):
        """A referenced include that doesn't exist locally ends up in ff_includes."""
        (tmp_path / "topol.top").write_text(
            '#include "ff/forcefield.itp"\n'
            '#include "phantom.itp"\n'
        )
        from builders.step_builders.preparation_builder import _ITP_WALK_CODE
        ns = {}
        exec("from pathlib import Path\n" + _ITP_WALK_CODE, ns)
        _walk = ns["_walk_includes"]
        local_itps, ff_includes = _walk(tmp_path / "topol.top", tmp_path)

        assert "phantom.itp" not in local_itps
        assert "phantom.itp" in ff_includes

    # ── run_protein_topology.py script ────────────────────────────────────────

    def test_generated_script_has_walk_includes(self, tmp_path):
        """The generated run_protein_topology.py must embed _walk_includes."""
        from core.compiler import SimulationCompiler
        from builders.workspace_builder import WorkspaceBuilder
        import os
        yaml_text = """\
project:
  name: test_walk_in_script
components:
  - id: protein_1
    role: protein
    file: protein.pdb
environment:
  membrane:
    enabled: true
forcefields:
  protein: opls-aa
simulation_objectives:
  - membrane_insertion
"""
        (tmp_path / "protein.pdb").write_text("ATOM      1  CA  ALA A   1       1.0   2.0   3.0  1.00  0.00\n")
        (tmp_path / "test.yaml").write_text(yaml_text)
        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            state = SimulationCompiler().compile(str(tmp_path / "test.yaml"))
        finally:
            os.chdir(orig)
        ws = tmp_path / "ws"
        WorkspaceBuilder().build(state, workspace_path=ws)

        gen_prot_dir = next((ws / "steps").glob("*_generate_protein_topology"))
        script = (gen_prot_dir / "run_protein_topology.py").read_text()
        assert "_walk_includes" in script, (
            "run_protein_topology.py must embed the _walk_includes include graph walker"
        )
        assert "molecule_itps" in script, (
            "run_protein_topology.py must populate molecule_itps in manifest"
        )
        assert "position_restraint_itps" in script, (
            "run_protein_topology.py must populate position_restraint_itps in manifest"
        )

    # ── assemble_system_topology reads manifest ────────────────────────────────

    def test_assemble_system_topology_reads_manifest(self, tmp_path):
        from core.compiler import SimulationCompiler
        from builders.workspace_builder import WorkspaceBuilder
        import os
        yaml_text = """\
project:
  name: test_asm_manifest
components:
  - id: protein_1
    role: protein
    file: protein.pdb
environment:
  membrane:
    enabled: true
forcefields:
  protein: opls-aa
simulation_objectives:
  - membrane_insertion
"""
        (tmp_path / "protein.pdb").write_text("ATOM      1  CA  ALA A   1       1.0   2.0   3.0  1.00  0.00\n")
        (tmp_path / "test.yaml").write_text(yaml_text)
        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            state = SimulationCompiler().compile(str(tmp_path / "test.yaml"))
        finally:
            os.chdir(orig)
        ws = tmp_path / "ws"
        WorkspaceBuilder().build(state, workspace_path=ws)

        asm_dir = next((ws / "steps").glob("*_assemble_system_topology"), None)
        if asm_dir is None:
            return
        script = (asm_dir / "run_assemble_system.py").read_text()
        assert "protein_topology_manifest.json" in script, (
            "run_assemble_system.py must read protein_topology_manifest.json"
        )
        assert "_molecule_itps_p" in script, (
            "run_assemble_system.py must use manifest's molecule_itps to build includes"
        )
        assert "posre.itp" not in script or "_prot_dir_rel" in script, (
            "run_assemble_system.py must not hardcode posre.itp copy"
        )

    def test_assemble_system_topology_uses_relative_includes(self, tmp_path):
        from core.compiler import SimulationCompiler
        from builders.workspace_builder import WorkspaceBuilder
        import os
        yaml_text = """\
project:
  name: test_rel_includes
components:
  - id: protein_1
    role: protein
    file: protein.pdb
environment:
  membrane:
    enabled: true
forcefields:
  protein: opls-aa
simulation_objectives:
  - membrane_insertion
"""
        (tmp_path / "protein.pdb").write_text("ATOM      1  CA  ALA A   1       1.0   2.0   3.0  1.00  0.00\n")
        (tmp_path / "test.yaml").write_text(yaml_text)
        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            state = SimulationCompiler().compile(str(tmp_path / "test.yaml"))
        finally:
            os.chdir(orig)
        ws = tmp_path / "ws"
        WorkspaceBuilder().build(state, workspace_path=ws)

        asm_dir = next((ws / "steps").glob("*_assemble_system_topology"), None)
        if asm_dir is None:
            return
        script = (asm_dir / "run_assemble_system.py").read_text()
        # Must compute relative path from SCRIPT_DIR to PROT_TOP_DIR
        assert "_prot_dir_rel" in script, (
            "run_assemble_system.py must compute _prot_dir_rel to make includes portable"
        )
        assert 'path.relpath' in script, (
            "run_assemble_system.py must use os.path.relpath to compute protein include path"
        )

    # ── ARTIFACT_VALIDATION_ERROR ──────────────────────────────────────────────

    def test_artifact_validation_error_exists(self):
        from executors.remediation_models import ErrorCategory
        assert hasattr(ErrorCategory, "ARTIFACT_VALIDATION_ERROR")
        assert ErrorCategory.ARTIFACT_VALIDATION_ERROR.value == "artifact_validation_error"

    def test_artifact_validation_error_plan_is_non_retryable(self, tmp_path):
        from executors.remediation_models import ErrorCategory, ErrorSeverity, DiagnosisResult
        from executors.signal_detector import AdaptiveReasoner
        diag = DiagnosisResult(
            step_id        = "generate_protein_topology",
            step_dir       = str(tmp_path),
            category       = ErrorCategory.ARTIFACT_VALIDATION_ERROR,
            severity       = ErrorSeverity.FATAL,
            confidence     = 0.99,
            primary_signal = "posre.itp not found after exit=0",
            explanation    = "pdb2gmx succeeded but posre.itp is absent (multichain output)",
        )
        plan = AdaptiveReasoner().plan_remediation(diag, tmp_path)
        assert plan.is_applicable is False
        assert plan.max_retries == 0
        assert plan.requires_human is True


class TestTopologyDuplicateDefaults:
    """Regression tests for the 'Found a second defaults directive' bug.

    Root cause: pdb2gmx writes the FF include as an absolute path
    (e.g. #include "/abs/.../oplsaa_membrane.ff/forcefield.itp").
    _walk_includes resolves this via Path(base_d)/inc → absolute path exists → classified
    as LOCAL, not FF. The manifest stores forcefield.itp in molecule_itps.
    assemble_system_topology then adds its canonical FF include at the top AND
    a second include from molecule_itps → GROMACS sees two [ defaults ] sections.
    """

    # ── _walk_includes FF-path classification ──────────────────────────────────

    def test_walk_includes_skips_ff_dir_paths(self, tmp_path):
        """_walk_includes must classify .ff/ paths as FF includes, not local."""
        from builders.step_builders.preparation_builder import _ITP_WALK_CODE
        assert '.ff/' in _ITP_WALK_CODE, (
            "_ITP_WALK_CODE must check for .ff/ in _inc to classify FF includes"
        )
        assert "continue" in _ITP_WALK_CODE, (
            "_ITP_WALK_CODE must skip FF includes without recursing into them"
        )

    def test_walk_includes_absolute_ff_path_goes_to_ff_not_local(self, tmp_path):
        """Absolute-path FF include must land in _ff, not _local."""
        from builders.step_builders.preparation_builder import _ITP_WALK_CODE

        # Write a fake topol.top with an absolute-path forcefield include
        ff_dir = tmp_path / "membrane_assets" / "oplsaa_membrane.ff"
        ff_dir.mkdir(parents=True)
        ff_itp = ff_dir / "forcefield.itp"
        ff_itp.write_text("[ defaults ]\n1 3 yes 0.5 0.5\n#include \"ffnonbonded.itp\"\n")
        (ff_dir / "ffnonbonded.itp").write_text("; nonbonded params\n")

        step_dir = tmp_path / "step"
        step_dir.mkdir()
        topol = step_dir / "topol.top"
        topol.write_text(
            f'#include "{ff_itp}"\n'
            f'#include "topol_chain_A.itp"\n'
        )
        (step_dir / "topol_chain_A.itp").write_text("[ moleculetype ]\nProtein_chain_A 3\n")

        ns: dict = {}
        exec(compile(_ITP_WALK_CODE, "<_ITP_WALK_CODE>", "exec"), ns)
        _walk_includes = ns["_walk_includes"]

        local, ff = _walk_includes(topol, step_dir)
        assert "topol_chain_A.itp" in local, "chain itp must be local"
        assert str(ff_itp) not in local, (
            f"forcefield.itp absolute path must NOT be in local; got local={local}"
        )
        ff_flat = " ".join(ff)
        assert "forcefield" in ff_flat or "oplsaa_membrane.ff" in ff_flat, (
            f"forcefield.itp must be in ff list; got ff={ff}"
        )

    def test_molecule_itps_filter_excludes_ff_paths(self, tmp_path):
        """generate_protein_topology script must exclude .ff/ paths from molecule_itps."""
        from builders.workspace_builder import WorkspaceBuilder
        from core.compiler import SimulationCompiler
        import os

        (tmp_path / "protein.pdb").write_text(
            "ATOM      1  CA  ALA A   1       1.0   2.0   3.0  1.00  0.00\n"
        )
        (tmp_path / "test.yaml").write_text("""\
project:
  name: test_ff_filter
components:
  - id: protein_1
    role: protein
    file: protein.pdb
environment:
  membrane:
    enabled: true
forcefields:
  protein: opls-aa
simulation_objectives:
  - membrane_insertion
""")
        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            state = SimulationCompiler().compile(str(tmp_path / "test.yaml"))
        finally:
            os.chdir(orig)
        ws = tmp_path / "ws"
        WorkspaceBuilder().build(state, workspace_path=ws)
        prot_dir = next((ws / "steps").glob("*_generate_protein_topology"), None)
        if prot_dir is None:
            return
        script = (prot_dir / "run_protein_topology.py").read_text()
        assert "'.ff/' not in _f" in script or '".ff/" not in _f' in script, (
            "generate_protein_topology must filter out .ff/ paths from molecule_itps"
        )

    # ── assemble_system_topology FF guard ──────────────────────────────────────

    def test_assemble_topology_skips_ff_in_molecule_itps_loop(self, tmp_path):
        """assemble_system_topology must skip any FF file in the molecule_itps loop."""
        from builders.workspace_builder import WorkspaceBuilder
        from core.compiler import SimulationCompiler
        import os

        (tmp_path / "protein.pdb").write_text(
            "ATOM      1  CA  ALA A   1       1.0   2.0   3.0  1.00  0.00\n"
        )
        (tmp_path / "test.yaml").write_text("""\
project:
  name: test_asm_ff_guard
components:
  - id: protein_1
    role: protein
    file: protein.pdb
environment:
  membrane:
    enabled: true
forcefields:
  protein: opls-aa
simulation_objectives:
  - membrane_insertion
""")
        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            state = SimulationCompiler().compile(str(tmp_path / "test.yaml"))
        finally:
            os.chdir(orig)
        ws = tmp_path / "ws"
        WorkspaceBuilder().build(state, workspace_path=ws)
        asm_dir = next((ws / "steps").glob("*_assemble_system_topology"), None)
        if asm_dir is None:
            return
        script = (asm_dir / "run_assemble_system.py").read_text()
        assert "'.ff/' in _itp" in script or '".ff/" in _itp' in script, (
            "assemble_system_topology must guard the molecule_itps loop against .ff/ entries"
        )

    def test_assemble_topology_has_single_forcefield_include_line(self, tmp_path):
        """The assembled topology script must add forcefield.itp exactly once."""
        from builders.workspace_builder import WorkspaceBuilder
        from core.compiler import SimulationCompiler
        import os

        (tmp_path / "protein.pdb").write_text(
            "ATOM      1  CA  ALA A   1       1.0   2.0   3.0  1.00  0.00\n"
        )
        (tmp_path / "test.yaml").write_text("""\
project:
  name: test_single_ff
components:
  - id: protein_1
    role: protein
    file: protein.pdb
environment:
  membrane:
    enabled: true
forcefields:
  protein: opls-aa
simulation_objectives:
  - membrane_insertion
""")
        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            state = SimulationCompiler().compile(str(tmp_path / "test.yaml"))
        finally:
            os.chdir(orig)
        ws = tmp_path / "ws"
        WorkspaceBuilder().build(state, workspace_path=ws)
        asm_dir = next((ws / "steps").glob("*_assemble_system_topology"), None)
        if asm_dir is None:
            return
        script = (asm_dir / "run_assemble_system.py").read_text()
        ff_lines = [l for l in script.splitlines() if "forcefield.itp" in l and "#include" in l]
        assert len(ff_lines) == 1, (
            f"script must reference forcefield.itp exactly once; found:\n" +
            "\n".join(ff_lines)
        )

    # ── _clean_topology helper ─────────────────────────────────────────────────

    def test_lipid_itp_map_present_in_assemble_script(self, tmp_path):
        """assemble_system_topology must define _LIPID_ITP_MAP for static ITP lookup."""
        from builders.workspace_builder import WorkspaceBuilder
        from core.compiler import SimulationCompiler
        import os

        (tmp_path / "protein.pdb").write_text(
            "ATOM      1  CA  ALA A   1       1.0   2.0   3.0  1.00  0.00\n"
        )
        (tmp_path / "test.yaml").write_text("""\
project:
  name: test_lipid_itp_map
components:
  - id: protein_1
    role: protein
    file: protein.pdb
environment:
  membrane:
    enabled: true
forcefields:
  protein: opls-aa
simulation_objectives:
  - membrane_insertion
""")
        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            state = SimulationCompiler().compile(str(tmp_path / "test.yaml"))
        finally:
            os.chdir(orig)
        ws = tmp_path / "ws"
        WorkspaceBuilder().build(state, workspace_path=ws)
        asm_dir = next((ws / "steps").glob("*_assemble_system_topology"), None)
        if asm_dir is None:
            return
        script = (asm_dir / "run_assemble_system.py").read_text()
        assert "_LIPID_ITP_MAP" in script, (
            "assemble_system_topology must define _LIPID_ITP_MAP for static lipid ITP lookup "
            "(replaces pdb2gmx on lipids_only.gro — lipids are defined only in RTP, not as "
            "standalone ITPs that pdb2gmx can reliably process)"
        )
        assert "dpp.itp" in script, (
            "assemble_system_topology _LIPID_ITP_MAP must map DPP to dpp.itp"
        )

    def test_clean_topology_strips_ff_includes_from_standalone_itp(self, tmp_path):
        """_clean_topology must strip [ defaults ] and FF includes from standalone lipid ITPs."""
        import shutil, re as _re, os as _os
        from pathlib import Path

        ff_dir = tmp_path / "ff"
        ff_dir.mkdir()
        (ff_dir / "forcefield.itp").write_text("[ defaults ]\n1 3 yes 0.5 0.5\n")

        standalone_itp = tmp_path / "dppc.itp"
        standalone_itp.write_text(
            '#include "oplsaa_membrane.ff/forcefield.itp"\n'
            "[ defaults ]\n"
            "1 3 yes 0.5 0.5\n"
            "[ moleculetype ]\n"
            "DPP\n"
            "1\n"
            "[ atoms ]\n"
            ";   nr  type ...\n"
            "1   opls_800  1  DPP  C1  1  -0.3\n"
            "[ system ]\n"
            "DPP bilayer\n"
            "[ molecules ]\n"
            "DPP 100\n"
        )
        lipid_topol = tmp_path / "lipid_topol.top"
        lipid_topol.write_text(
            '#include "oplsaa_membrane.ff/forcefield.itp"\n'
            f'#include "{standalone_itp.name}"\n'
            "[ system ]\n"
            "Bilayer\n"
        )
        shutil.copy(standalone_itp, ff_dir / "dppc.itp")

        clean_code = r"""
import re, os
from pathlib import Path

def _clean_topology(text, ff_dir, sdir):
    lines, out, skip = text.splitlines(), [], False
    for ln in lines:
        _m = re.match(r'^\s*\[\s*(\w+)\s*\]', ln)
        if _m: skip = _m.group(1).lower() in ('defaults', 'system', 'molecules')
        if skip: continue
        _ls = ln.lstrip()
        if _ls.startswith('#include'):
            _pts = _ls.split('"')
            _inc = _pts[1] if len(_pts) >= 2 else ''
            if '.ff/' in _inc or _inc.endswith('.ff'): continue
            _cand = (ff_dir / Path(_inc).name)
            if not _cand.exists(): _cand = (Path(sdir) / _inc)
            if _cand.exists():
                _ct = _cand.read_text()
                if re.search(r'^\s*\[\s*defaults\s*\]', _ct, re.MULTILINE | re.IGNORECASE):
                    out.append(_clean_topology(_ct, ff_dir, sdir))
                else:
                    _rel = os.path.relpath(str(_cand), sdir)
                    out.append(f'#include "{_rel}"')
                continue
        out.append(ln)
    return '\n'.join(out)
"""
        ns: dict = {}
        exec(clean_code, ns)
        _clean_topology = ns["_clean_topology"]

        result = _clean_topology(lipid_topol.read_text(), ff_dir, str(tmp_path))
        assert "[ defaults ]" not in result
        assert "forcefield.itp" not in result
        assert "[ moleculetype ]" in result
        assert "[ atoms ]" in result
        assert "[ system ]" not in result
        assert "[ molecules ]" not in result

    # ── _count_defaults validator ──────────────────────────────────────────────

    def test_count_defaults_validator_present(self, tmp_path):
        """assemble_system_topology must validate exactly one [ defaults ] before finishing."""
        from builders.workspace_builder import WorkspaceBuilder
        from core.compiler import SimulationCompiler
        import os

        (tmp_path / "protein.pdb").write_text(
            "ATOM      1  CA  ALA A   1       1.0   2.0   3.0  1.00  0.00\n"
        )
        (tmp_path / "test.yaml").write_text("""\
project:
  name: test_count_defaults
components:
  - id: protein_1
    role: protein
    file: protein.pdb
environment:
  membrane:
    enabled: true
forcefields:
  protein: opls-aa
simulation_objectives:
  - membrane_insertion
""")
        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            state = SimulationCompiler().compile(str(tmp_path / "test.yaml"))
        finally:
            os.chdir(orig)
        ws = tmp_path / "ws"
        WorkspaceBuilder().build(state, workspace_path=ws)
        asm_dir = next((ws / "steps").glob("*_assemble_system_topology"), None)
        if asm_dir is None:
            return
        script = (asm_dir / "run_assemble_system.py").read_text()
        assert "_count_defaults" in script
        assert "TOPOLOGY_DUPLICATE_DEFAULTS_ERROR" in script

    def test_count_defaults_counts_like_gromacs(self, tmp_path):
        """_count_defaults must count siblings independently (GROMACS-style, no dedup)."""
        ff_dir = tmp_path / "ff"
        ff_dir.mkdir()
        ff_itp = ff_dir / "forcefield.itp"
        ff_itp.write_text("[ defaults ]\n1 3 yes 0.5 0.5\n")

        topol = tmp_path / "topol.top"
        topol.write_text(
            f'#include "{ff_itp}"\n'
            f'[ moleculetype ]\nProtein 3\n'
            f'#include "{ff_itp}"\n'
        )

        count_code = r"""
import re
from pathlib import Path

def _count_defaults(f, _anc=None):
    if _anc is None: _anc = set()
    key = str(Path(f).resolve())
    if key in _anc: return 0
    _child_anc = _anc | {key}
    n = 0
    try: txt = Path(f).read_text()
    except OSError: return 0
    for _l in txt.splitlines():
        if re.match(r'^\s*\[\s*defaults\s*\]', _l, re.IGNORECASE): n += 1
        _ls = _l.lstrip()
        if _ls.startswith('#include'):
            _pts = _ls.split('"')
            if len(_pts) >= 2:
                _cand = (Path(f).parent / _pts[1]).resolve()
                if _cand.exists(): n += _count_defaults(str(_cand), _child_anc)
    return n
"""
        ns: dict = {}
        exec(count_code, ns)
        _count_defaults = ns["_count_defaults"]
        assert _count_defaults(str(topol)) == 2, (
            "_count_defaults must count both includes of forcefield.itp (GROMACS-style)"
        )

    # ── Error category and signal detection ────────────────────────────────────

    def test_topology_duplicate_defaults_error_in_error_category(self):
        """TOPOLOGY_DUPLICATE_DEFAULTS_ERROR must be a valid ErrorCategory."""
        from executors.remediation_models import ErrorCategory
        assert hasattr(ErrorCategory, "TOPOLOGY_DUPLICATE_DEFAULTS_ERROR")
        assert ErrorCategory.TOPOLOGY_DUPLICATE_DEFAULTS_ERROR == "topology_duplicate_defaults_error"

    def test_second_defaults_directive_detected(self, tmp_path):
        """'Found a second defaults directive' must map to TOPOLOGY_DUPLICATE_DEFAULTS_ERROR."""
        from executors.signal_detector import _SIGNAL_PATTERNS
        from executors.remediation_models import ErrorCategory

        stderr = "Fatal error:\nFound a second defaults directive\n"
        matched = None
        for sig in _SIGNAL_PATTERNS:
            hit, _ = sig.match(stderr)
            if hit:
                matched = sig.category
                break
        assert matched == ErrorCategory.TOPOLOGY_DUPLICATE_DEFAULTS_ERROR, (
            f"'Found a second defaults directive' must map to TOPOLOGY_DUPLICATE_DEFAULTS_ERROR; "
            f"got {matched}"
        )

    def test_duplicate_defaults_plan_is_non_retryable(self, tmp_path):
        """TOPOLOGY_DUPLICATE_DEFAULTS_ERROR must produce a non-retryable plan."""
        from executors.remediation_models import ErrorCategory, ErrorSeverity, DiagnosisResult
        from executors.signal_detector import AdaptiveReasoner

        diag = DiagnosisResult(
            step_id        = "membrane_embedding",
            step_dir       = str(tmp_path),
            category       = ErrorCategory.TOPOLOGY_DUPLICATE_DEFAULTS_ERROR,
            severity       = ErrorSeverity.FATAL,
            confidence     = 0.99,
            primary_signal = "Found a second defaults directive",
            explanation    = "forcefield.itp included twice",
        )
        plan = AdaptiveReasoner().plan_remediation(diag, tmp_path)
        assert plan.is_applicable is False
        assert plan.max_retries == 0
        assert plan.requires_human is True

# ═══════════════════════════════════════════════════════════════════════════════
# K. Forcefield include path fix (Phase 11 — Task: FF include normalization)
# ═══════════════════════════════════════════════════════════════════════════════

_FF_TEST_YAML = """\
project:
  name: test_ff_include
components:
  - id: protein_1
    role: protein
    file: protein.pdb
environment:
  membrane:
    enabled: true
forcefields:
  protein: opls-aa
simulation_objectives:
  - membrane_insertion
"""


def _build_ff_workspace(tmp_path):
    from core.compiler import SimulationCompiler
    from builders.workspace_builder import WorkspaceBuilder
    import os
    (tmp_path / "protein.pdb").write_text(
        "ATOM      1  CA  ALA A   1       1.0   2.0   3.0  1.00  0.00\n"
    )
    (tmp_path / "test.yaml").write_text(_FF_TEST_YAML)
    orig = os.getcwd()
    os.chdir(tmp_path)
    try:
        state = SimulationCompiler().compile(str(tmp_path / "test.yaml"))
    finally:
        os.chdir(orig)
    ws = tmp_path / "ws"
    WorkspaceBuilder().build(state, workspace_path=ws)
    return ws


class TestForcefieldIncludeFix:
    """Regression tests for the forcefield include path normalization.

    Root cause of the bug: `_ff_incs_p` from the manifest contained whatever
    pdb2gmx wrote (GMXLIB-relative path like ``oplsaa_membrane.ff/forcefield.itp``).
    Writing that verbatim into the assembled topol.top caused grompp to fail
    when GMXLIB was not set in the embedding environment.

    Fix: always write a filesystem-relative path computed from ASSETS_DIR.
    """

    def test_generate_topology_script_uses_assets_ff_rel(self, tmp_path):
        ws = _build_ff_workspace(tmp_path)
        gen_top = next((ws / "steps").glob("*_generate_topology"), None)
        if gen_top is None:
            return
        script = (gen_top / "run_topology.py").read_text()
        assert "_assets_ff_rel" in script, (
            "run_topology.py must compute _assets_ff_rel for the forcefield include"
        )

    def test_generate_topology_no_verbatim_ff_incs_copy(self, tmp_path):
        ws = _build_ff_workspace(tmp_path)
        gen_top = next((ws / "steps").glob("*_generate_topology"), None)
        if gen_top is None:
            return
        script = (gen_top / "run_topology.py").read_text()
        assert "for _ffi in _ff_incs_p" not in script, (
            "run_topology.py must NOT loop over _ff_incs_p to copy verbatim FF includes"
        )

    def test_generate_topology_single_ff_include_line(self, tmp_path):
        ws = _build_ff_workspace(tmp_path)
        gen_top = next((ws / "steps").glob("*_generate_topology"), None)
        if gen_top is None:
            return
        script = (gen_top / "run_topology.py").read_text()
        ff_lines = [l for l in script.splitlines() if "forcefield.itp" in l]
        assert len(ff_lines) == 1, (
            f"Expected exactly 1 forcefield.itp reference in run_topology.py, got {len(ff_lines)}: {ff_lines}"
        )

    def test_generate_topology_no_ffnonbonded_in_script(self, tmp_path):
        ws = _build_ff_workspace(tmp_path)
        gen_top = next((ws / "steps").glob("*_generate_topology"), None)
        if gen_top is None:
            return
        script = (gen_top / "run_topology.py").read_text()
        assert "ffnonbonded" not in script, (
            "run_topology.py must never include ffnonbonded.itp directly"
        )

    def test_assemble_system_topology_uses_assets_ff_rel(self, tmp_path):
        ws = _build_ff_workspace(tmp_path)
        asm_dir = next((ws / "steps").glob("*_assemble_system_topology"), None)
        if asm_dir is None:
            return
        script = (asm_dir / "run_assemble_system.py").read_text()
        assert "_assets_ff_rel" in script, (
            "run_assemble_system.py must compute _assets_ff_rel for the forcefield include"
        )

    def test_assemble_system_topology_no_ffnonbonded(self, tmp_path):
        ws = _build_ff_workspace(tmp_path)
        asm_dir = next((ws / "steps").glob("*_assemble_system_topology"), None)
        if asm_dir is None:
            return
        script = (asm_dir / "run_assemble_system.py").read_text()
        assert "ffnonbonded" not in script, (
            "run_assemble_system.py must never include ffnonbonded.itp directly"
        )

    def test_generate_topology_has_include_graph_validator(self, tmp_path):
        ws = _build_ff_workspace(tmp_path)
        gen_top = next((ws / "steps").glob("*_generate_topology"), None)
        if gen_top is None:
            return
        script = (gen_top / "run_topology.py").read_text()
        assert "_chk_inc" in script, (
            "run_topology.py must contain the _chk_inc include-graph validator"
        )
        assert "_inc_errors" in script, (
            "run_topology.py must run the include-graph validator after writing topol.top"
        )

    def test_assemble_system_topology_has_include_graph_validator(self, tmp_path):
        ws = _build_ff_workspace(tmp_path)
        asm_dir = next((ws / "steps").glob("*_assemble_system_topology"), None)
        if asm_dir is None:
            return
        script = (asm_dir / "run_assemble_system.py").read_text()
        assert "_chk_inc" in script, (
            "run_assemble_system.py must contain the _chk_inc include-graph validator"
        )

    def test_embedding_shrink_loop_exports_gmxlib(self, tmp_path):
        ws = _build_ff_workspace(tmp_path)
        emb_dir = next((ws / "steps").glob("*_membrane_embedding"), None)
        if emb_dir is None:
            return
        run_sh = emb_dir / "run.sh"
        if not run_sh.exists():
            return
        content = run_sh.read_text()
        if "run_tm_aware.py" in content:
            return  # TM-aware backend — uses Python, different path
        assert "export GMXLIB" in content, (
            "run.sh (inflategro backend) must export GMXLIB so grompp finds the forcefield"
        )
        assert "membrane_assets" in content, (
            "GMXLIB in run.sh must point at the workspace membrane_assets directory"
        )

    def test_tm_aware_grompp_uses_ff_env(self, tmp_path):
        ws = _build_ff_workspace(tmp_path)
        emb_dir = next((ws / "steps").glob("*_membrane_embedding"), None)
        if emb_dir is None:
            return
        tm_script = emb_dir / "run_tm_aware.py"
        if not tm_script.exists():
            return
        script = tm_script.read_text()
        # The grompp call must pass env=ff_env, not run without env
        grompp_lines = [l for l in script.splitlines() if "grompp_cmd" in l and "run(" in l]
        for line in grompp_lines:
            assert "ff_env" in line, (
                f"grompp subprocess.run call must use env=ff_env so GMXLIB is set: {line!r}"
            )


# ═══════════════════════════════════════════════════════════════════════════════
# L. Include path normalization — absolute vs relative manifest paths
# ═══════════════════════════════════════════════════════════════════════════════

class TestIncludePathNormalization:
    """Regression tests for the manifest include path normalization bug.

    Root cause: pdb2gmx may write absolute paths in its generated topol.top
    (e.g. #include "/abs/path/to/topol_Protein_chain_A.itp").  _walk_includes
    then stores that absolute path in protein_topology_manifest.json. When the
    assembler naively prepends `../01_generate_protein_topology/`, the result is
    `../01_generate_protein_topology//abs/path/...` — an invalid path.

    Fix: always resolve each include via Path.resolve() relative to its source
    directory, then rewrite as os.path.relpath(abs, target_topology_dir).
    """

    def _build_workspace_with_manifest(self, tmp_path, molecule_itps):
        """Build a workspace and inject a custom manifest into generate_protein_topology."""
        from core.compiler import SimulationCompiler
        from builders.workspace_builder import WorkspaceBuilder
        import os, json as _json

        yaml_text = """\
project:
  name: test_path_norm
components:
  - id: protein_1
    role: protein
    file: protein.pdb
environment:
  membrane:
    enabled: true
forcefields:
  protein: opls-aa
simulation_objectives:
  - membrane_insertion
"""
        (tmp_path / "protein.pdb").write_text(
            "ATOM      1  CA  ALA A   1       1.0   2.0   3.0  1.00  0.00\n"
        )
        (tmp_path / "test.yaml").write_text(yaml_text)
        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            state = SimulationCompiler().compile(str(tmp_path / "test.yaml"))
        finally:
            os.chdir(orig)
        ws = tmp_path / "ws"
        WorkspaceBuilder().build(state, workspace_path=ws)

        # Inject a custom manifest with the given molecule_itps into generate_protein_topology
        prot_dir = next((ws / "steps").glob("*_generate_protein_topology"), None)
        if prot_dir:
            manifest = {
                "molecule_itps": molecule_itps,
                "forcefield_includes": ["oplsaa_membrane.ff/forcefield.itp"],
                "molecule_entries": [{"name": "Protein_chain_A", "count": 1}],
                "position_restraint_itps": [],
                "all_local_includes": molecule_itps,
                "passed": True,
                "validation_errors": [],
            }
            (prot_dir / "protein_topology_manifest.json").write_text(
                _json.dumps(manifest, indent=2)
            )
            # Create stub .itp files so the include validator can resolve them
            for itp in molecule_itps:
                p = Path(itp)
                if p.is_absolute():
                    p.parent.mkdir(parents=True, exist_ok=True)
                    p.write_text("; stub\n[ moleculetype ]\nProtein_chain_A 3\n")
                else:
                    (prot_dir / itp).write_text("; stub\n[ moleculetype ]\nProtein_chain_A 3\n")
        return ws

    def test_relative_manifest_path_writes_correct_include(self, tmp_path):
        ws = self._build_workspace_with_manifest(
            tmp_path, ["topol_Protein_chain_A.itp"]
        )
        asm_dir = next((ws / "steps").glob("*_assemble_system_topology"), None)
        if asm_dir is None:
            return
        script = (asm_dir / "run_assemble_system.py").read_text()
        # Must normalize relative to absolute then back to relative
        assert "_itp_abs" in script, (
            "assembler must normalize includes via _itp_abs = (PROT_TOP_DIR / _itp).resolve()"
        )
        assert "_itp_rel" in script, (
            "assembler must compute _itp_rel = relpath(_itp_abs, SCRIPT_DIR)"
        )

    def test_absolute_manifest_path_normalized(self, tmp_path):
        abs_itp = str(tmp_path / "abs_dir" / "topol_Protein_chain_A.itp")
        ws = self._build_workspace_with_manifest(tmp_path, [abs_itp])
        asm_dir = next((ws / "steps").glob("*_assemble_system_topology"), None)
        if asm_dir is None:
            return
        script = (asm_dir / "run_assemble_system.py").read_text()
        # The normalization code must handle absolute paths via Path(_itp).is_absolute()
        assert "is_absolute()" in script, (
            "assembler must check Path(_itp).is_absolute() to handle absolute manifest paths"
        )

    def test_no_prot_dir_rel_slash_itp_concatenation(self, tmp_path):
        """Generated scripts must never use `f'..{_prot_dir_rel}/{_itp}...'` to build includes."""
        ws = self._build_workspace_with_manifest(
            tmp_path, ["topol_Protein_chain_A.itp"]
        )
        for asm_dir in (ws / "steps").glob("*_assemble_system_topology"):
            script = (asm_dir / "run_assemble_system.py").read_text()
            # The exact bad pattern: f-string concat of _prot_dir_rel + "/" + _itp
            bad_lines = [l for l in script.splitlines()
                         if "_prot_dir_rel}" in l and "/_itp}" in l]
            assert not bad_lines, (
                f"assembler must not build includes as {{_prot_dir_rel}}/{{_itp}}: {bad_lines}"
            )

    def test_include_path_normalization_in_topology_assembler_script(self, tmp_path):
        """The assembler script must use the normalize pattern for every itp."""
        ws = self._build_workspace_with_manifest(
            tmp_path, ["topol_Protein_chain_A.itp"]
        )
        asm_dir = next((ws / "steps").glob("*_assemble_system_topology"), None)
        if asm_dir is None:
            return
        script = (asm_dir / "run_assemble_system.py").read_text()
        # Must resolve + relpath, not naive string concatenation
        assert "is_absolute()" in script
        assert "_itp_abs" in script
        assert "_itp_rel = _os.path.relpath" in script
        assert f'#include "{{}}"' not in script  # no old-style f-string concat

    def test_dag_has_no_generate_topology_step(self, tmp_path):
        """DAG must not contain generate_topology (the legacy assembler is removed)."""
        from core.compiler import SimulationCompiler
        import os, json as _json

        (tmp_path / "protein.pdb").write_text(
            "ATOM      1  CA  ALA A   1       1.0   2.0   3.0  1.00  0.00\n"
        )
        yaml_text = """\
project:
  name: test_dag_clean
components:
  - id: protein_1
    role: protein
    file: protein.pdb
environment:
  membrane:
    enabled: true
forcefields:
  protein: opls-aa
simulation_objectives:
  - membrane_insertion
"""
        (tmp_path / "test.yaml").write_text(yaml_text)
        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            result = SimulationCompiler().compile(str(tmp_path / "test.yaml"))
        finally:
            os.chdir(orig)
        step_ids = [s.step_id for s in result.plan.steps]
        assert "generate_topology" not in step_ids, (
            f"DAG must not contain generate_topology; found steps: {step_ids}"
        )
        assert "assemble_system_topology" in step_ids, (
            "DAG must contain assemble_system_topology as the single topology assembler"
        )

    def test_topology_include_resolution_error_is_non_retryable(self, tmp_path):
        from executors.remediation_models import ErrorCategory, ErrorSeverity, DiagnosisResult
        from executors.signal_detector import AdaptiveReasoner
        diag = DiagnosisResult(
            step_id        = "assemble_system_topology",
            step_dir       = str(tmp_path),
            category       = ErrorCategory.TOPOLOGY_INCLUDE_RESOLUTION_ERROR,
            severity       = ErrorSeverity.FATAL,
            confidence     = 0.99,
            primary_signal = "Cannot resolve '../01_generate_protein_topology//home/...' in topol.top",
            explanation    = "Absolute path concatenated with relative prefix",
        )
        plan = AdaptiveReasoner().plan_remediation(diag, tmp_path)
        assert plan.is_applicable is False
        assert plan.max_retries == 0
        assert plan.requires_human is True


# ═══════════════════════════════════════════════════════════════════════════════
# M. Undefined moleculetype — multichain protein + DPP 1 lipid regression
# ═══════════════════════════════════════════════════════════════════════════════

class TestUndefinedMoleculeType:
    """Regression tests for 'No such moleculetype Protein' (and 'DPP 478 vs DPP 1').

    Root cause: assemble_system_topology wrote 'Protein 1' (stale hardcoded generic
    name) when pdb2gmx generated 'Protein_chain_A 1' through 'Protein_chain_E 1' for
    a 5-chain GLP-1R. Also, lipid residue counting from GRO produced 'DPP 478' but the
    pdb2gmx-generated topology defines a single combined 'DPP' moleculetype (DPP 1).
    """

    # ── Script structure checks ────────────────────────────────────────────────

    def _get_asm_script(self, tmp_path: Path) -> str:
        from core.compiler import SimulationCompiler
        from builders.workspace_builder import WorkspaceBuilder
        import os
        yaml_text = """\
project:
  name: test_moltype
components:
  - id: protein_1
    role: protein
    file: protein.pdb
environment:
  membrane:
    enabled: true
forcefields:
  protein: opls-aa
simulation_objectives:
  - membrane_insertion
"""
        (tmp_path / "protein.pdb").write_text(
            "ATOM      1  CA  ALA A   1       1.0   2.0   3.0  1.00  0.00\n"
        )
        (tmp_path / "test.yaml").write_text(yaml_text)
        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            state = SimulationCompiler().compile(str(tmp_path / "test.yaml"))
        finally:
            os.chdir(orig)
        ws = tmp_path / "ws"
        WorkspaceBuilder().build(state, workspace_path=ws)
        asm_dir = next((ws / "steps").glob("*_assemble_system_topology"), None)
        if asm_dir is None:
            return ""
        return (asm_dir / "run_assemble_system.py").read_text()

    def test_script_uses_molecule_entries_primary(self, tmp_path):
        """run_assemble_system.py must use molecule_entries from manifest (priority 1)."""
        script = self._get_asm_script(tmp_path)
        assert "_mol_entries_p" in script, (
            "script must read molecule_entries (priority 1 protein mol source)"
        )
        assert "molecule_entries" in script

    def test_script_uses_molecule_names_fallback(self, tmp_path):
        """run_assemble_system.py must fall back to molecule_names (priority 2)."""
        script = self._get_asm_script(tmp_path)
        assert "molecule_names" in script, (
            "script must fall back to molecule_names when molecule_entries is absent"
        )
        assert "_mol_names_p" in script

    def test_script_parses_moleculetype_from_itp_tertiary(self, tmp_path):
        """run_assemble_system.py must parse [ moleculetype ] from ITPs as tertiary fallback."""
        script = self._get_asm_script(tmp_path)
        assert "moleculetype" in script.lower(), (
            "script must parse [ moleculetype ] from ITP files as tertiary fallback"
        )
        assert "_in_mt" in script or "in_mt" in script

    def test_script_uses_static_itp_for_lipid_topology(self, tmp_path):
        """run_assemble_system.py must use static ITP include, not pdb2gmx lipid_topol.top."""
        script = self._get_asm_script(tmp_path)
        assert "_LIPID_ITP_MAP" in script, (
            "script must use _LIPID_ITP_MAP for static ITP lookup — "
            "pdb2gmx on lipids_only.gro is removed (DPP is only in aminoacids.rtp)"
        )
        assert "_lip_itp_path" in script or "_lip_itp_parts" in script, (
            "script must resolve lipid ITP path from ASSETS_DIR/oplsaa_membrane.ff/"
        )

    def test_script_uses_gro_count_for_lipid_molecules(self, tmp_path):
        """Lipid [ molecules ] count must come from the GRO residue count, not pdb2gmx."""
        script = self._get_asm_script(tmp_path)
        # GRO residue count is always used (no pdb2gmx fallback / override needed)
        assert "lipid_counts" in script, (
            "script must use lipid_counts (from GRO) for [ molecules ] lipid entries"
        )
        assert "_lipid_mol_entries" not in script, (
            "script must NOT use _lipid_mol_entries — the pdb2gmx mol count approach "
            "is removed; GRO residue count is authoritative"
        )

    def test_script_has_collect_moleculetypes_validator(self, tmp_path):
        """run_assemble_system.py must define _collect_moleculetypes validator."""
        script = self._get_asm_script(tmp_path)
        assert "_collect_moleculetypes" in script, (
            "script must define _collect_moleculetypes to validate [ molecules ] consistency"
        )
        assert "TOPOLOGY_UNDEFINED_MOLECULETYPE_ERROR" in script

    # ── _collect_moleculetypes functional tests ───────────────────────────────

    def test_collect_moleculetypes_finds_chain_types(self, tmp_path):
        """_collect_moleculetypes must discover all [ moleculetype ] names recursively."""
        topol = tmp_path / "topol.top"
        topol.write_text(
            '#include "ff/forcefield.itp"\n'
            '#include "topol_Protein_chain_A.itp"\n'
            '#include "topol_Protein_chain_B.itp"\n'
            '#include "topol_DPP.itp"\n'
            '[ system ]\nSystem\n'
            '[ molecules ]\nProtein_chain_A 1\nProtein_chain_B 1\nDPP 1\n'
        )
        for chain in ("A", "B"):
            (tmp_path / f"topol_Protein_chain_{chain}.itp").write_text(
                f"[ moleculetype ]\nProtein_chain_{chain}    3\n"
                "[ atoms ]\n; ...\n"
            )
        (tmp_path / "topol_DPP.itp").write_text(
            "[ moleculetype ]\nDPP    1\n"
            "[ atoms ]\n; ...\n"
        )
        ff_dir = tmp_path / "ff"
        ff_dir.mkdir()
        (ff_dir / "forcefield.itp").write_text("[ defaults ]\n1 3 yes 0.5 0.5\n")

        code = r"""
import re
from pathlib import Path

def _collect_moleculetypes(f, visited=None):
    if visited is None: visited = set()
    key = str(Path(f).resolve())
    if key in visited: return []
    visited.add(key)
    result = []
    try: txt = Path(f).read_text()
    except OSError: return result
    _in_mt = False
    for _l in txt.splitlines():
        _ls = _l.lstrip()
        if re.match(r'^\[\s*moleculetype\s*\]', _ls, re.IGNORECASE):
            _in_mt = True; continue
        if _in_mt:
            if _ls and not _ls.startswith(';') and not _ls.startswith('['):
                result.append(_ls.split()[0]); _in_mt = False
            elif _ls.startswith('['):
                _in_mt = False
        if _ls.startswith('#include'):
            _pts = _ls.split('"')
            if len(_pts) >= 2:
                _cand = (Path(f).parent / _pts[1]).resolve()
                if _cand.exists(): result += _collect_moleculetypes(str(_cand), visited)
    return result
"""
        ns: dict = {}
        exec(code, ns)
        types = set(ns["_collect_moleculetypes"](str(topol)))
        assert "Protein_chain_A" in types
        assert "Protein_chain_B" in types
        assert "DPP" in types
        assert "Protein" not in types

    def test_undefined_type_detected(self, tmp_path):
        """Validator must catch 'Protein 1' when no [ moleculetype ] Protein is defined."""
        topol = tmp_path / "topol.top"
        topol.write_text(
            '#include "topol_Protein_chain_A.itp"\n'
            '[ system ]\nSystem\n'
            '[ molecules ]\nProtein 1\n'  # WRONG: should be Protein_chain_A
        )
        (tmp_path / "topol_Protein_chain_A.itp").write_text(
            "[ moleculetype ]\nProtein_chain_A    3\n"
        )

        code = r"""
import re
from pathlib import Path

def _collect_moleculetypes(f, visited=None):
    if visited is None: visited = set()
    key = str(Path(f).resolve())
    if key in visited: return []
    visited.add(key)
    result = []
    try: txt = Path(f).read_text()
    except OSError: return result
    _in_mt = False
    for _l in txt.splitlines():
        _ls = _l.lstrip()
        if re.match(r'^\[\s*moleculetype\s*\]', _ls, re.IGNORECASE):
            _in_mt = True; continue
        if _in_mt:
            if _ls and not _ls.startswith(';') and not _ls.startswith('['):
                result.append(_ls.split()[0]); _in_mt = False
            elif _ls.startswith('['):
                _in_mt = False
        if _ls.startswith('#include'):
            _pts = _ls.split('"')
            if len(_pts) >= 2:
                _cand = (Path(f).parent / _pts[1]).resolve()
                if _cand.exists(): result += _collect_moleculetypes(str(_cand), visited)
    return result
"""
        ns: dict = {}
        exec(code, ns)
        defined = set(ns["_collect_moleculetypes"](str(topol)))
        mol_text = topol.read_text()
        import re
        mol_m = re.search(r'^\s*\[\s*molecules\s*\](.*)', mol_text,
                          re.MULTILINE | re.IGNORECASE | re.DOTALL)
        mol_names = []
        if mol_m:
            for ml in mol_m.group(1).splitlines():
                ml = ml.strip()
                if ml and not ml.startswith(";"):
                    ps = ml.split()
                    if ps: mol_names.append(ps[0])
        undefined = [n for n in mol_names if n not in defined]
        assert "Protein" in undefined, (
            "Validator must detect 'Protein' as undefined when only 'Protein_chain_A' is defined"
        )
        assert "Protein_chain_A" not in undefined

    def test_known_working_topology_passes_validation(self, tmp_path):
        """The known working topology structure (5-chain + DPP 1 + SOL + ions) must validate."""
        ff_dir = tmp_path / "membrane_assets" / "oplsaa_membrane.ff"
        ff_dir.mkdir(parents=True)
        (ff_dir / "forcefield.itp").write_text("[ defaults ]\n1 3 yes 0.5 0.5\n")
        (ff_dir / "fbamem.itp").write_text("; water\n[ moleculetype ]\nSOL 1\n")
        (ff_dir / "ions.itp").write_text(
            "[ moleculetype ]\nNA    1\n[ moleculetype ]\nCL    1\n"
        )

        for chain in ("A", "B", "C", "D", "E"):
            (tmp_path / f"topol_Protein_chain_{chain}.itp").write_text(
                f"[ moleculetype ]\nProtein_chain_{chain}    3\n"
            )
        (tmp_path / "topol_DPP.itp").write_text("[ moleculetype ]\nDPP    1\n")

        topol = tmp_path / "topol.top"
        topol.write_text(
            '#include "membrane_assets/oplsaa_membrane.ff/forcefield.itp"\n'
            '#include "topol_Protein_chain_A.itp"\n'
            '#include "topol_Protein_chain_B.itp"\n'
            '#include "topol_Protein_chain_C.itp"\n'
            '#include "topol_Protein_chain_D.itp"\n'
            '#include "topol_Protein_chain_E.itp"\n'
            '#include "topol_DPP.itp"\n'
            '#include "membrane_assets/oplsaa_membrane.ff/fbamem.itp"\n'
            '#include "membrane_assets/oplsaa_membrane.ff/ions.itp"\n'
            '[ system ]\nIonic Channel in Water\n'
            '[ molecules ]\n'
            'Protein_chain_A     1\n'
            'Protein_chain_B     1\n'
            'Protein_chain_C     1\n'
            'Protein_chain_D     1\n'
            'Protein_chain_E     1\n'
            'DPP                 1\n'
            'SOL             75166\n'
            'NA                307\n'
            'CL                332\n'
        )

        code = r"""
import re
from pathlib import Path

def _collect_moleculetypes(f, visited=None):
    if visited is None: visited = set()
    key = str(Path(f).resolve())
    if key in visited: return []
    visited.add(key)
    result = []
    try: txt = Path(f).read_text()
    except OSError: return result
    _in_mt = False
    for _l in txt.splitlines():
        _ls = _l.lstrip()
        if re.match(r'^\[\s*moleculetype\s*\]', _ls, re.IGNORECASE):
            _in_mt = True; continue
        if _in_mt:
            if _ls and not _ls.startswith(';') and not _ls.startswith('['):
                result.append(_ls.split()[0]); _in_mt = False
            elif _ls.startswith('['):
                _in_mt = False
        if _ls.startswith('#include'):
            _pts = _ls.split('"')
            if len(_pts) >= 2:
                _cand = (Path(f).parent / _pts[1]).resolve()
                if _cand.exists(): result += _collect_moleculetypes(str(_cand), visited)
    return result
"""
        ns: dict = {}
        exec(code, ns)
        defined = set(ns["_collect_moleculetypes"](str(topol)))
        import re
        mol_text = topol.read_text()
        mol_m = re.search(r'^\s*\[\s*molecules\s*\](.*)', mol_text,
                          re.MULTILINE | re.IGNORECASE | re.DOTALL)
        mol_names = []
        if mol_m:
            for ml in mol_m.group(1).splitlines():
                ml = ml.strip()
                if ml and not ml.startswith(";"):
                    ps = ml.split()
                    if ps: mol_names.append(ps[0])
        undefined = [n for n in mol_names if n not in defined]
        assert undefined == [], (
            f"Known working topology must validate without undefined moleculetypes; "
            f"got undefined={undefined}, defined={defined}"
        )

    def test_dpp_one_preserved_from_pdb2gmx(self, tmp_path):
        """When lipid_topol.top has 'DPP 1', _lipid_mol_entries must be [('DPP', 1)]."""
        import re
        lipid_topol = tmp_path / "lipid_topol.top"
        lipid_topol.write_text(
            '#include "oplsaa_membrane.ff/forcefield.itp"\n'
            '#include "topol_DPP.itp"\n'
            '#include "oplsaa_membrane.ff/fbamem.itp"\n'
            '#include "oplsaa_membrane.ff/ions.itp"\n'
            '[ system ]\nDPP bilayer\n'
            '[ molecules ]\nDPP    1\n'
        )
        lt = lipid_topol.read_text()
        lip_mol_m = re.search(r'^\s*\[\s*molecules\s*\](.*)', lt,
                               re.MULTILINE | re.IGNORECASE | re.DOTALL)
        lipid_mol_entries = []
        if lip_mol_m:
            for lml in lip_mol_m.group(1).splitlines():
                lml = lml.strip()
                if lml and not lml.startswith(";"):
                    lps = lml.split()
                    if lps:
                        lipid_mol_entries.append((lps[0], int(lps[1]) if len(lps) > 1 else 1))
        assert lipid_mol_entries == [("DPP", 1)], (
            f"pdb2gmx 'DPP 1' must be preserved exactly; got {lipid_mol_entries}"
        )

    # ── Error category and signal ─────────────────────────────────────────────

    def test_topology_undefined_moleculetype_error_in_error_category(self):
        from executors.remediation_models import ErrorCategory
        assert hasattr(ErrorCategory, "TOPOLOGY_UNDEFINED_MOLECULETYPE_ERROR")
        assert ErrorCategory.TOPOLOGY_UNDEFINED_MOLECULETYPE_ERROR == "topology_undefined_moleculetype_error"

    def test_no_such_moleculetype_detected(self, tmp_path):
        """'No such moleculetype' must map to TOPOLOGY_UNDEFINED_MOLECULETYPE_ERROR."""
        from executors.signal_detector import _SIGNAL_PATTERNS
        from executors.remediation_models import ErrorCategory

        stderr = (
            "ERROR 1 [file topol.top, line 126734]:\n"
            "  No such moleculetype Protein\n\n"
            "There was 1 error in input file(s)\n"
        )
        matched = None
        for sig in _SIGNAL_PATTERNS:
            hit, _ = sig.match(stderr)
            if hit:
                matched = sig.category
                break
        assert matched == ErrorCategory.TOPOLOGY_UNDEFINED_MOLECULETYPE_ERROR, (
            f"'No such moleculetype' must map to TOPOLOGY_UNDEFINED_MOLECULETYPE_ERROR; "
            f"got {matched}"
        )

    def test_assembler_diagnostic_detected(self, tmp_path):
        """TOPOLOGY_UNDEFINED_MOLECULETYPE_ERROR prefix must be detectable in stderr."""
        from executors.signal_detector import _SIGNAL_PATTERNS
        from executors.remediation_models import ErrorCategory

        stderr = (
            "TOPOLOGY_UNDEFINED_MOLECULETYPE_ERROR: [ molecules ] contains 'Protein', "
            "but no included [ moleculetype ] named 'Protein' exists.\n"
        )
        matched = None
        for sig in _SIGNAL_PATTERNS:
            hit, _ = sig.match(stderr)
            if hit:
                matched = sig.category
                break
        assert matched == ErrorCategory.TOPOLOGY_UNDEFINED_MOLECULETYPE_ERROR

    def test_undefined_moleculetype_plan_is_non_retryable(self, tmp_path):
        from executors.remediation_models import ErrorCategory, ErrorSeverity, DiagnosisResult
        from executors.signal_detector import AdaptiveReasoner

        diag = DiagnosisResult(
            step_id        = "assemble_system_topology",
            step_dir       = str(tmp_path),
            category       = ErrorCategory.TOPOLOGY_UNDEFINED_MOLECULETYPE_ERROR,
            severity       = ErrorSeverity.FATAL,
            confidence     = 0.99,
            primary_signal = "No such moleculetype Protein",
            explanation    = "Protein 1 in [ molecules ] but only Protein_chain_A defined",
        )
        plan = AdaptiveReasoner().plan_remediation(diag, tmp_path)
        assert plan.is_applicable is False
        assert plan.max_retries == 0
        assert plan.requires_human is True
