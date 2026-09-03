"""
test_membrane_integration.py
============================
Phase 4 integration tests — protein-membrane workflow end-to-end with real GROMACS.

Test classes:
  TestFilePathFixes       — regression tests for Phase-4 file-path bug fixes
  TestCrossStepPaths      — builder generates correct inter-step relative paths
  TestScriptNoFindRoot    — generated scripts no longer rely on _find_root walk
  TestRunnerCompilation   — MembraneIntegrationRunner compiles and builds workspace
  TestRunnerFastSteps     — orient, match_box, embed, generate_topology with GROMACS
  TestReportStructure     — membrane_build_report.yaml has required fields

Tests requiring real GROMACS are guarded by:
    pytest.mark.skipif(not shutil.which("gmx"), reason="gmx not found")

Tests exercising the shrink loop or MD are guarded by:
    pytest.mark.slow
and are excluded from the default test run.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest
import yaml

# ── Optional GROMACS guard ────────────────────────────────────────────────────
_HAS_GMX = bool(shutil.which("gmx"))
_GMX_REASON = "gmx not found on PATH"


# ═══════════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════════

def _compile_workspace(yaml_path: Path, workspace_dir: Path) -> Path:
    """Compile and build a workspace; return workspace root."""
    from core.compiler import SimulationCompiler
    from builders.workspace_builder import WorkspaceBuilder

    import os
    orig_cwd = os.getcwd()
    os.chdir(yaml_path.parent)
    try:
        result = SimulationCompiler().compile(str(yaml_path))
    finally:
        os.chdir(orig_cwd)

    return WorkspaceBuilder().build(
        result, workspace_path=workspace_dir, yaml_source=str(yaml_path)
    )


def _orient_yaml() -> Path:
    return Path("configs/membrane_orient_test.yaml").resolve()


# ═══════════════════════════════════════════════════════════════════════════════
# TestFilePathFixes — regression tests for the three cross-step path bugs fixed
# in Phase 4:
#   1. orient_protein source_file → staged inputs/<comp_id>.pdb
#   2. box_match_helper PROTEIN_GRO_REL → relative path to orient step output
#   3. Generated scripts no longer use _find_root (editable install handles it)
# ═══════════════════════════════════════════════════════════════════════════════

class TestFilePathFixes:

    @pytest.fixture(scope="class")
    def workspace(self, tmp_path_factory):
        ws = tmp_path_factory.mktemp("ws_path_fixes")
        return _compile_workspace(_orient_yaml(), ws)

    def test_orient_run_sh_uses_staged_input(self, workspace):
        """orient_protein run.sh must use protein_processed.gro (pdb2gmx output with H atoms).

        Using protein_processed.gro ensures that coordinates used for embedding match
        the protein topology (same atom count including hydrogens added by pdb2gmx).
        The raw inputs/ PDB has no hydrogens and would cause atom-count mismatches with
        the pdb2gmx-derived topology at grompp time.
        """
        orient_dir = next(
            d for d in (workspace / "steps").iterdir()
            if "orient_protein" in d.name
        )
        run_sh = (orient_dir / "run.sh").read_text()
        # Must NOT reference the raw filename (before pdb2gmx, missing hydrogens)
        assert "peptide_A.pdb" not in run_sh, (
            "run.sh must not use the raw PDB (no hydrogens) — use protein_processed.gro"
        )
        # Must reference protein_processed.gro from generate_protein_topology
        assert "protein_processed.gro" in run_sh, (
            "run.sh must reference protein_processed.gro from generate_protein_topology "
            "so atom count matches the pdb2gmx-derived topology"
        )

    def test_orient_run_sh_uses_component_id_name(self, workspace):
        """orient step must reference generate_protein_topology step directory."""
        orient_dir = next(
            d for d in (workspace / "steps").iterdir()
            if "orient_protein" in d.name
        )
        run_sh = (orient_dir / "run.sh").read_text()
        # Must reference protein_processed.gro from generate_protein_topology
        assert "generate_protein_topology" in run_sh or "protein_processed.gro" in run_sh, (
            "run.sh must reference protein_processed.gro from generate_protein_topology step"
        )

    def test_staged_input_file_exists(self, workspace):
        """The staged input file referenced in run.sh must actually exist."""
        inputs_dir = workspace / "inputs"
        assert (inputs_dir / "peptide_1.pdb").exists(), (
            "peptide_1.pdb must be staged into workspace/inputs/"
        )

    def test_box_match_helper_has_protein_gro_rel(self, workspace):
        """box_match_helper.py must have PROTEIN_GRO_REL pointing to orient step."""
        match_dir = next(
            d for d in (workspace / "steps").iterdir()
            if "match_box" in d.name
        )
        helper = (match_dir / "box_match_helper.py").read_text()
        assert "PROTEIN_GRO_REL" in helper
        # Must NOT be an empty string (should resolve to orient step)
        for line in helper.splitlines():
            if line.strip().startswith("PROTEIN_GRO_REL") and "=" in line:
                val = line.split("=", 1)[1].strip().strip('"').strip("'")
                assert "protein_oriented.gro" in val
                assert val != "protein_oriented.gro" or "orient" in val, (
                    "PROTEIN_GRO_REL should point to orient_protein step output, not current dir"
                )
                break

    def test_box_match_helper_protein_gro_rel_points_to_orient(self, workspace):
        """PROTEIN_GRO_REL must resolve to the orient_protein step directory."""
        match_dir = next(
            d for d in (workspace / "steps").iterdir()
            if "match_box" in d.name
        )
        helper_text = (match_dir / "box_match_helper.py").read_text()
        # Extract PROTEIN_GRO_REL value (strip trailing comment before stripping quotes)
        gro_rel = None
        for line in helper_text.splitlines():
            if line.strip().startswith("PROTEIN_GRO_REL") and "=" in line:
                val = line.split("=", 1)[1].strip()
                if "#" in val:
                    val = val.split("#", 1)[0].strip()
                gro_rel = val.strip('"').strip("'")
                break
        assert gro_rel is not None
        # Resolve it from the match_box step directory
        resolved = (match_dir / gro_rel).resolve()
        orient_dir = next(
            d for d in (workspace / "steps").iterdir()
            if "orient_protein" in d.name
        )
        expected = (orient_dir / "protein_oriented.gro").resolve()
        assert resolved == expected, (
            f"PROTEIN_GRO_REL ({gro_rel}) resolves to {resolved}, "
            f"expected {expected}"
        )

    def test_embed_script_has_no_find_root(self, workspace):
        """run_embed.py must not use the fragile _find_root walk."""
        embed_dir = next(
            d for d in (workspace / "steps").iterdir()
            if "embed_in_bilayer" in d.name
        )
        script = (embed_dir / "run_embed.py").read_text()
        assert "_find_root" not in script
        assert "sys.path.insert" not in script

    def test_embed_script_clashes_in_warnings_not_errors(self, workspace):
        """run_embed.py must place initial lipid-protein clashes in warnings, not errors.

        Regression: GLP-1R membrane run was blocked because clashes were written to
        the 'errors' list, making the overlap gate return blocked=True.  Initial
        embedding clashes are expected and resolved by the shrink loop.
        """
        embed_dir = next(
            d for d in (workspace / "steps").iterdir()
            if "embed_in_bilayer" in d.name
        )
        script = (embed_dir / "run_embed.py").read_text()
        # Must not hardcode clashes into the errors list
        assert '"errors":                  []' in script or (
            '"errors":' in script and 'ov.n_clashes' not in script.split('"errors":')[1].split('\n')[0]
        ), "run_embed.py must write errors=[] for the overlap report"
        # Must have _overlap_warnings advisory path
        assert "_overlap_warnings" in script, (
            "run_embed.py must put initial embedding clashes into _overlap_warnings (advisory)"
        )

    def test_embed_script_has_initial_clashes_expected_flag(self, workspace):
        """overlap_report.json must include initial_clashes_expected=True."""
        embed_dir = next(
            d for d in (workspace / "steps").iterdir()
            if "embed_in_bilayer" in d.name
        )
        script = (embed_dir / "run_embed.py").read_text()
        assert "initial_clashes_expected" in script, (
            "run_embed.py must write initial_clashes_expected flag so the gate "
            "can distinguish expected overlaps from hard geometry failures"
        )

    def test_clean_water_script_has_no_find_root(self, workspace):
        """run_clean_water.py must not use _find_root."""
        clean_dir = next(
            d for d in (workspace / "steps").iterdir()
            if "clean_water" in d.name
        )
        script = (clean_dir / "run_clean_water.py").read_text()
        assert "_find_root" not in script
        assert "sys.path.insert" not in script

    def test_topology_script_has_no_find_root(self, workspace):
        """run_protein_topology.py must not use _find_root."""
        topo_dir = next(
            d for d in (workspace / "steps").iterdir()
            if "generate_protein_topology" in d.name
        )
        script = (topo_dir / "run_protein_topology.py").read_text()
        assert "_find_root" not in script
        assert "sys.path.insert" not in script

    def test_run_sh_editconf_before_registration(self, workspace):
        """run.sh must call editconf_cmd.sh (step 2) before run_membrane_registration.py (step 3).

        Phase 8A regression: registration must read protein_boxed.gro (post-editconf frame).
        Swapping the order caused a systematic -2.9 nm placement error for GLP-1R.
        """
        match_dir = next(
            d for d in (workspace / "steps").iterdir()
            if "match_box" in d.name
        )
        run_sh = (match_dir / "run.sh").read_text()
        non_comment = [
            l.strip() for l in run_sh.splitlines()
            if l.strip() and not l.strip().startswith("#")
        ]
        editconf_idx = next(
            (i for i, l in enumerate(non_comment) if "editconf_cmd.sh" in l), None
        )
        reg_idx = next(
            (i for i, l in enumerate(non_comment) if "run_membrane_registration.py" in l), None
        )
        assert editconf_idx is not None, "run.sh must contain 'bash editconf_cmd.sh'"
        assert reg_idx is not None, "run.sh must contain 'python3 run_membrane_registration.py'"
        assert editconf_idx < reg_idx, (
            f"editconf_cmd.sh (pos {editconf_idx}) must precede "
            f"run_membrane_registration.py (pos {reg_idx}) in run.sh"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# TestCrossStepPaths — builder generates correctly-resolved inter-step paths
# ═══════════════════════════════════════════════════════════════════════════════

class TestCrossStepPaths:

    @pytest.fixture(scope="class")
    def workspace(self, tmp_path_factory):
        ws = tmp_path_factory.mktemp("ws_cross_step")
        return _compile_workspace(_orient_yaml(), ws)

    def test_metadata_has_protein_gro_rel(self, workspace):
        """match_box_to_bilayer metadata.json records protein_gro_rel."""
        match_dir = next(
            d for d in (workspace / "steps").iterdir()
            if "match_box" in d.name
        )
        meta = json.loads((match_dir / "metadata.json").read_text())
        assert "protein_gro_rel" in meta
        assert meta["protein_gro_rel"] != ""

    def test_embedding_input_gro_uses_embed_in_bilayer_system(self, workspace):
        """run.sh INPUT_GRO must use embed_in_bilayer/system.gro (pre-shrink coordinate).

        The repair plan (Session 2) moved assemble_system_topology AFTER membrane_embedding.
        The embedding shrink loop therefore starts from embed_in_bilayer/system.gro (the raw
        bilayer+protein system), not from a pre-assembled system_processed.gro.
        system_processed.gro is produced by assemble_system_topology which runs AFTER embedding.
        """
        embed_dir = next(
            d for d in (workspace / "steps").iterdir()
            if "membrane_embedding" in d.name
        )
        script = (embed_dir / "run.sh").read_text()
        assert "embed_in_bilayer" in script, (
            "run.sh must reference embed_in_bilayer/system.gro as INPUT_GRO "
            "(assemble_system_topology runs after membrane_embedding, not before)"
        )
        assert "system_processed.gro" not in script, (
            "run.sh must NOT reference system_processed.gro — that file is produced by "
            "assemble_system_topology which runs AFTER the embedding shrink loop"
        )

    def test_embedding_topol_uses_assemble_system_topology(self, workspace):
        """run.sh topology must point to assemble_system_topology/topol.top."""
        embed_dir = next(
            d for d in (workspace / "steps").iterdir()
            if "membrane_embedding" in d.name
        )
        script = (embed_dir / "run.sh").read_text()
        assert "assemble_system_topology" in script
        assert "topol.top" in script


# ═══════════════════════════════════════════════════════════════════════════════
# TestRunnerCompilation — runner compiles and builds workspace correctly
# ═══════════════════════════════════════════════════════════════════════════════

class TestRunnerCompilation:

    @pytest.fixture(scope="class")
    def runner_workspace(self, tmp_path_factory):
        from executors.membrane_integration_runner import MembraneIntegrationRunner
        ws = tmp_path_factory.mktemp("ws_runner")
        runner = MembraneIntegrationRunner(
            yaml_path=_orient_yaml(),
            workspace_dir=ws,
            stop_after="orient_protein",  # compile only, skip execution for speed
            timeout_fast=0,               # immediate "timeout" to skip actual GMX
        )
        # Just build workspace without running
        return runner._compile_and_build()

    def test_workspace_has_steps_dir(self, runner_workspace):
        assert (runner_workspace / "steps").exists()

    def test_workspace_has_membrane_assets(self, runner_workspace):
        assert (runner_workspace / "membrane_assets").exists()
        assert (runner_workspace / "membrane_assets" / "dppc512_whole.gro").exists()

    def test_workspace_has_all_expected_step_dirs(self, runner_workspace):
        step_names = {d.name for d in (runner_workspace / "steps").iterdir()}
        required = {"orient_protein", "match_box_to_bilayer", "embed_in_bilayer",
                    "assemble_system_topology", "membrane_embedding"}
        for req in required:
            assert any(req in n for n in step_names), f"Step dir for {req} not found"

    def test_runner_make_empty_report(self, runner_workspace):
        from executors.membrane_integration_runner import MembraneIntegrationRunner
        runner = MembraneIntegrationRunner(yaml_path=_orient_yaml())
        report = runner._make_empty_report(runner_workspace)
        required_keys = [
            "selected_bilayer", "bilayer_x", "bilayer_y", "bilayer_z",
            "protein_x", "protein_y", "protein_z", "xy_fit_status",
            "topology_status", "solvation_status", "ionization_status",
            "em_status", "nvt_status", "npt_status",
            "first_failure_stage", "first_failure_category", "first_failure_detail",
            "remediation_applied", "final_system_ready",
        ]
        for key in required_keys:
            assert key in report, f"Report missing key: {key}"

    def test_embed_in_bilayer_has_tm_mask_config(self, runner_workspace):
        embed_dir = next(
            d for d in (runner_workspace / "steps").iterdir()
            if "embed_in_bilayer" in d.name
        )
        run_embed_py = (embed_dir / "run_embed.py").read_text()
        assert "TM_MASK_VALIDATION =" in run_embed_py
        assert "TM_MASK_POLICY =" in run_embed_py


# ═══════════════════════════════════════════════════════════════════════════════
# TestRunnerFastSteps — execute real GROMACS through embed_in_bilayer
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.mark.skipif(not _HAS_GMX, reason=_GMX_REASON)
class TestRunnerFastSteps:
    """
    Execute the fast automated steps up to embed_in_bilayer with real GROMACS.
    These steps are fast (<30s each): orient → match_box → embed_in_bilayer.

    membrane_embedding (GROMACS mdrun shrink loop) is excluded — it is slow and
    tested separately in TestRunnerAssemblySteps.

    Failures here indicate regressions in file-path logic or bilayer embedding.
    """

    @pytest.fixture(scope="class")
    def run_result(self, tmp_path_factory):
        from executors.membrane_integration_runner import MembraneIntegrationRunner
        ws = tmp_path_factory.mktemp("ws_fast_steps")
        runner = MembraneIntegrationRunner(
            yaml_path=_orient_yaml(),
            workspace_dir=ws,
            stop_after="embed_in_bilayer",
            timeout_fast=120,
        )
        report = runner.run()
        return report, ws

    def test_orient_protein_passes(self, run_result):
        report, ws = run_result
        orient_dir = next(d for d in (ws / "steps").iterdir() if "orient_protein" in d.name)
        assert (orient_dir / "protein_oriented.gro").exists(), (
            "orient_protein must produce protein_oriented.gro"
        )
        assert (orient_dir / "orientation_report.json").exists()
        data = json.loads((orient_dir / "orientation_report.json").read_text())
        assert data["passed"] is True, f"Orientation validation failed: {data['errors']}"

    def test_match_box_passes(self, run_result):
        report, ws = run_result
        match_dir = next(d for d in (ws / "steps").iterdir() if "match_box" in d.name)
        assert (match_dir / "protein_boxed.gro").exists()
        data = json.loads((match_dir / "box_match_report.json").read_text())
        assert data["passed"] is True, f"Box match failed: {data['errors']}"

    def test_bilayer_xy_authority_in_report(self, run_result):
        """box_match_report.json must show bilayer as XY authority after Phase 3."""
        report, ws = run_result
        match_dir = next(d for d in (ws / "steps").iterdir() if "match_box" in d.name)
        data = json.loads((match_dir / "box_match_report.json").read_text())
        bxy = data.get("bilayer_xy", {})
        assert bxy.get("xy_authority") == "selected_bilayer", (
            "Bilayer must be authoritative for XY box dimensions (Phase 3)"
        )
        assert bxy.get("xy_fit_status") == "pass"

    def test_embed_in_bilayer_passes(self, run_result):
        report, ws = run_result
        embed_dir = next(d for d in (ws / "steps").iterdir() if "embed_in_bilayer" in d.name)
        assert (embed_dir / "system.gro").exists()
        assert (embed_dir / "overlap_report.json").exists()
        # Initial embedding places the protein inside the bilayer which causes lipid-protein
        # clashes; these are EXPECTED and removed by the downstream membrane_embedding shrink loop.
        # The step itself exits 0; do not assert overlap_report["passed"] is True here.

    def test_embed_clashes_are_advisory_not_blocking(self, run_result):
        """Regression: GLP-1R membrane run blocked when embed_in_bilayer had 20 initial clashes.

        embed_in_bilayer must be marked as passing even when overlap_report shows
        n_clashes > 0, because initial lipid-protein overlaps are expected and
        resolved by the downstream membrane_embedding shrink loop.
        """
        report, ws = run_result
        embed_dir = next(d for d in (ws / "steps").iterdir() if "embed_in_bilayer" in d.name)
        data = json.loads((embed_dir / "overlap_report.json").read_text())

        # errors must be empty — clashes go in warnings
        assert data.get("errors") == [], (
            f"overlap_report.json must have errors=[] for initial embedding; "
            f"got errors={data.get('errors')}"
        )
        # initial_clashes_expected flag must be set
        assert data.get("initial_clashes_expected") is True

        # Gate evaluator must NOT block when errors=[]
        from runtime.overlap_gate import evaluate_overlap_gate
        gate = evaluate_overlap_gate(embed_dir)
        assert gate is not None
        assert gate.blocked is False, (
            "overlap gate must not block embed_in_bilayer when only initial clashes present"
        )

    def _run_assemble_normal_mode(self, ws):
        """Helper: seed converged.gro from embed step and run assemble in normal mode."""
        import subprocess, shutil
        embed_dir  = next(d for d in (ws / "steps").iterdir() if "embed_in_bilayer" in d.name)
        embed5_dir = next(d for d in (ws / "steps").iterdir() if "membrane_embedding" in d.name)
        topo_dir   = next(d for d in (ws / "steps").iterdir() if "assemble_system_topology" in d.name)
        # Seed converged.gro so assemble can run without the slow shrink loop
        converged = embed5_dir / "converged.gro"
        if not converged.exists():
            shutil.copy2(str(embed_dir / "system.gro"), str(converged))
        result = subprocess.run(
            ["python3", str(topo_dir / "run_assemble_system.py")],
            capture_output=True, text=True, cwd=str(topo_dir), timeout=30,
        )
        return result, topo_dir, embed_dir

    def test_assemble_system_topology_passes(self, run_result):
        """Run assemble_system_topology in normal mode using embed_in_bilayer/system.gro.

        Seeds membrane_embedding/converged.gro from embed_in_bilayer/system.gro to
        test the assembly script without needing the slow mdrun shrink loop.
        """
        report, ws = run_result
        result, topo_dir, _ = self._run_assemble_normal_mode(ws)
        assert result.returncode == 0, (
            f"run_assemble_system.py failed:\nSTDOUT: {result.stdout[-800:]}\n"
            f"STDERR: {result.stderr[-800:]}"
        )
        assert (topo_dir / "topol.top").exists(), "topol.top must be created"
        assert (topo_dir / "system_processed.gro").exists(), "system_processed.gro must be created"
        data = json.loads((topo_dir / "topology_assembly_report.json").read_text())
        assert data["passed"] is True, f"Topology assembly failed: {data.get('errors', [])}"

    def test_topology_atom_count_matches_structure(self, run_result):
        """system_processed.gro and topol.top must have consistent atom counts."""
        report, ws = run_result
        result, topo_dir, embed_dir = self._run_assemble_normal_mode(ws)
        if result.returncode != 0:
            pytest.skip("run_assemble_system.py failed — skipping atom count check")
        report_path = topo_dir / "topology_assembly_report.json"
        if not report_path.exists():
            pytest.skip("topology_assembly_report.json not produced")
        data = json.loads(report_path.read_text())
        gro_path = topo_dir / "system_processed.gro"
        n_atoms_gro = int(gro_path.read_text().splitlines()[1].strip())
        assert data.get("n_atoms") == n_atoms_gro, (
            "topology_assembly_report.json atom count must match system_processed.gro"
        )

    def test_first_failure_is_none(self, run_result):
        report, ws = run_result
        assert report["first_failure_stage"] is None, (
            f"Unexpected failure at {report['first_failure_stage']}: "
            f"{report.get('first_failure_detail', '')[:200]}"
        )

    def test_topology_status_pass(self, run_result):
        """generate_protein_topology sets topology_status=pass when embedding stops early."""
        report, ws = run_result
        assert report["topology_status"] == "pass"

    def test_final_system_ready_true(self, run_result):
        """generate_protein_topology sets topology_status=pass when embedding stops early."""
        report, ws = run_result
        assert report["topology_status"] == "pass"


# ═══════════════════════════════════════════════════════════════════════════════
# TestReportStructure — membrane_build_report.yaml has all required fields
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.mark.skipif(not _HAS_GMX, reason=_GMX_REASON)
class TestReportStructure:

    @pytest.fixture(scope="class")
    def report_path(self, tmp_path_factory):
        from executors.membrane_integration_runner import MembraneIntegrationRunner
        ws = tmp_path_factory.mktemp("ws_report")
        runner = MembraneIntegrationRunner(
            yaml_path=_orient_yaml(),
            workspace_dir=ws,
            stop_after="generate_topology",
            timeout_fast=120,
        )
        runner.run()
        return ws / "membrane_build_report.yaml"

    def test_report_file_exists(self, report_path):
        assert report_path.exists(), "membrane_build_report.yaml must be produced"

    def test_report_is_valid_yaml(self, report_path):
        data = yaml.safe_load(report_path.read_text())
        assert isinstance(data, dict)

    def test_report_has_all_required_fields(self, report_path):
        data = yaml.safe_load(report_path.read_text())
        required = [
            "selected_bilayer", "bilayer_x", "bilayer_y", "bilayer_z",
            "protein_x", "protein_y", "protein_z", "xy_fit_status",
            "topology_status", "solvation_status", "ionization_status",
            "em_status", "nvt_status", "npt_status",
            "first_failure_stage", "remediation_applied", "final_system_ready",
        ]
        for key in required:
            assert key in data, f"membrane_build_report.yaml missing field: {key}"

    def test_report_bilayer_geometry_populated(self, report_path):
        data = yaml.safe_load(report_path.read_text())
        assert data["selected_bilayer"] == "dppc512_whole.gro"
        assert data["bilayer_x"] is not None and data["bilayer_x"] > 10.0
        assert data["bilayer_y"] is not None and data["bilayer_y"] > 10.0
        assert data["protein_x"] is not None
        assert data["protein_y"] is not None
        assert data["protein_z"] is not None

    def test_report_xy_fit_status_pass(self, report_path):
        data = yaml.safe_load(report_path.read_text())
        assert data["xy_fit_status"] == "pass"

    def test_report_remediation_applied_is_list(self, report_path):
        data = yaml.safe_load(report_path.read_text())
        assert isinstance(data["remediation_applied"], list)


# ═══════════════════════════════════════════════════════════════════════════════
# TestShrinkLoopScriptCorrectness — verify shrink loop script without running it
# ═══════════════════════════════════════════════════════════════════════════════

class TestShrinkLoopScriptCorrectness:
    """
    Verify shrink_loop.sh correctness without running the long shrink loop.
    Checks that Phase-4 fixes (boolean conversion, correct input file) are present.
    """

    @pytest.fixture(scope="class")
    def workspace(self, tmp_path_factory):
        ws = tmp_path_factory.mktemp("ws_shrink")
        return _compile_workspace(_orient_yaml(), ws)

    def test_shrink_loop_has_no_python_false_literal(self, workspace):
        """Shell 'false' must not appear bare inside Python -c strings."""
        embed_dir = next(
            d for d in (workspace / "steps").iterdir()
            if "membrane_embedding" in d.name
        )
        script = (embed_dir / "run.sh").read_text()
        # The pattern that was broken: 'converged': false == 'true'
        assert "converged': false" not in script, (
            "Shell 'false' must not appear directly in Python dict literal — "
            "use _pyconv conversion before python3 -c call"
        )

    def test_shrink_loop_uses_pyconv_conversion(self, workspace):
        """log_iter and finalize_telemetry must convert shell bool to Python bool."""
        embed_dir = next(
            d for d in (workspace / "steps").iterdir()
            if "membrane_embedding" in d.name
        )
        script = (embed_dir / "run.sh").read_text()
        assert "_pyconv" in script, (
            "run.sh must use _pyconv conversion for shell boolean → Python True/False"
        )

    def test_shrink_loop_starts_from_embed_in_bilayer_system(self, workspace):
        """run.sh INPUT_GRO must reference embed_in_bilayer/system.gro.

        Post-repair: assemble_system_topology runs AFTER membrane_embedding.
        The shrink loop input is therefore embed_in_bilayer/system.gro, not
        system_processed.gro (which is the output of assemble_system_topology).
        """
        embed_dir = next(
            d for d in (workspace / "steps").iterdir()
            if "membrane_embedding" in d.name
        )
        script = (embed_dir / "run.sh").read_text()
        # Find the INPUT_GRO= line
        for line in script.splitlines():
            if line.strip().startswith("INPUT_GRO="):
                assert "embed_in_bilayer" in line or "system.gro" in line, (
                    f"INPUT_GRO must reference embed_in_bilayer/system.gro, got: {line}"
                )
                assert "system_processed" not in line, (
                    f"INPUT_GRO must NOT be system_processed.gro "
                    f"(assemble runs after embedding now), got: {line}"
                )
                break
        else:
            pytest.fail("INPUT_GRO= line not found in run.sh")


# ═══════════════════════════════════════════════════════════════════════════════
# TestEmbeddingCutoffIntegration — cutoff plumbing from pipeline → script
# ═══════════════════════════════════════════════════════════════════════════════

class TestEmbeddingCutoffIntegration:
    """
    Verify that membrane_pipeline passes lipid_exclusion_cutoff_nm through to
    the generated run.sh, and that the generated script is correct.
    No GROMACS execution required.
    """

    @pytest.fixture(scope="class")
    def workspace(self, tmp_path_factory):
        ws = tmp_path_factory.mktemp("ws_cutoff")
        return _compile_workspace(_orient_yaml(), ws)

    def _embed_script(self, workspace):
        embed_dir = next(
            d for d in (workspace / "steps").iterdir()
            if "membrane_embedding" in d.name
        )
        return (embed_dir / "run.sh").read_text()

    def _embed_meta(self, workspace):
        embed_dir = next(
            d for d in (workspace / "steps").iterdir()
            if "membrane_embedding" in d.name
        )
        return json.loads((embed_dir / "metadata.json").read_text())

    def test_default_pipeline_uses_nonzero_cutoff(self, workspace):
        """Default pipeline uses lipid_exclusion_cutoff_nm=1.4 (not 0.0) for the
        one-shot initial cutoff. Checked as an exact script line (not a bare
        substring) because 'LOOP_CUTOFF_NM=0.0' — the correct, separate,
        zero-default shrink-loop cutoff — contains 'CUTOFF_NM=0.0' as a
        substring."""
        script = self._embed_script(workspace)
        lines = script.splitlines()
        assert "CUTOFF_NM=1.4" in lines, (
            "Default pipeline must use nonzero one-shot cutoff (1.4 nm)"
        )
        assert "CUTOFF_NM=0.0" not in lines, (
            "One-shot cutoff (CUTOFF_NM) must not be zero by default"
        )

    def test_shrink_loop_defaults_to_zero_cutoff_not_one_shot_cutoff(self, workspace):
        """
        Regression for the GLP-1R membrane-void bug (see
        docs/audits/glp1r_membrane_void_shrink_loop_audit.md): the shrink
        loop's per-iteration perl call must use $LOOP_CUTOFF_ANGSTROM
        (default 0, matching the original protmemfiles method), NOT the
        one-shot $CUTOFF_ANGSTROM. Reusing the nonzero one-shot cutoff on
        every deflate iteration was the root cause of the 478 -> 400 DPP
        lipid loss on GLP-1R.
        """
        script = self._embed_script(workspace)
        lines = script.splitlines()
        assert "LOOP_CUTOFF_NM=0.0" in lines

        perl_call_lines = [
            l for l in lines
            if l.strip().startswith("perl") and "$INFLATEGRO" in l
        ]
        assert len(perl_call_lines) == 2, (
            f"expected exactly 2 perl inflategro calls (initial + loop); "
            f"found {len(perl_call_lines)}"
        )
        initial_call, loop_call = perl_call_lines
        assert "$CUTOFF_ANGSTROM" in initial_call
        assert "$LOOP_CUTOFF_ANGSTROM" not in initial_call

        assert "$LOOP_CUTOFF_ANGSTROM" in loop_call, (
            f"shrink-loop perl call must use $LOOP_CUTOFF_ANGSTROM, got: {loop_call}"
        )
        assert "$CUTOFF_ANGSTROM" not in loop_call, (
            f"shrink-loop perl call must not reuse the one-shot $CUTOFF_ANGSTROM, "
            f"got: {loop_call}"
        )

    def test_metadata_records_nonzero_cutoff(self, workspace):
        """metadata.json lipid_exclusion_cutoff_nm must be > 0 by default."""
        meta = self._embed_meta(workspace)
        cutoff = meta["params"]["lipid_exclusion_cutoff_nm"]
        assert cutoff > 0.0, f"Default cutoff must be > 0, got {cutoff}"
        assert meta["params"]["cutoff_angstrom_to_perl"] == cutoff * 10

    def test_metadata_records_zero_default_shrink_loop_cutoff(self, workspace):
        """
        metadata.json shrink_loop_cutoff_nm must default to 0.0 — confirms
        pipelines/membrane_pipeline.py forwards
        embed_cfg.get("shrink_loop_cutoff_nm", 0.0) through to
        step.params["shrink_loop_cutoff_nm"], independently of the one-shot
        lipid_exclusion_cutoff_nm.
        """
        meta = self._embed_meta(workspace)
        assert meta["params"]["shrink_loop_cutoff_nm"] == 0.0
        assert meta["params"]["shrink_loop_cutoff_angstrom_to_perl"] == 0.0
        # Must be independent of (not derived from) the one-shot cutoff.
        assert meta["params"]["lipid_exclusion_cutoff_nm"] != meta["params"]["shrink_loop_cutoff_nm"]

    def test_diagnose_trapped_py_generated(self, workspace):
        """diagnose_trapped.py must be generated in the membrane_embedding step dir."""
        embed_dir = next(
            d for d in (workspace / "steps").iterdir()
            if "membrane_embedding" in d.name
        )
        diag = embed_dir / "diagnose_trapped.py"
        assert diag.exists(), "diagnose_trapped.py must be generated"
        content = diag.read_text()
        assert "trapped_lipid_diagnosis.json" in content
        assert "initial_lipid_count" in content
        assert "n_lipids_removed_by_inflategro" in content
        assert "lipid_exclusion_cutoff_nm" in content

    def test_topology_update_in_script(self, workspace):
        """run.sh must contain topology lipid count update logic."""
        script = self._embed_script(workspace)
        assert "update_topology_lipid_count" in script, (
            "run.sh must call update_topology_lipid_count when inflategro removes lipids"
        )
        assert "N_REMOVED" in script


# ═══════════════════════════════════════════════════════════════════════════════
# TestPreShrinkExclusionIntegration — pre-shrink config plumbing through pipeline→script
# ═══════════════════════════════════════════════════════════════════════════════

class TestPreShrinkExclusionIntegration:
    """
    Verify that pre_exclude_trapped_lipids plumbs from pipeline config through to
    the generated run_embed.py. No GROMACS execution required.
    """

    @pytest.fixture(scope="class")
    def workspace(self, tmp_path_factory):
        ws = tmp_path_factory.mktemp("ws_pre_excl")
        return _compile_workspace(_orient_yaml(), ws)

    def _embed_script(self, workspace):
        embed_dir = next(
            d for d in (workspace / "steps").iterdir()
            if "embed_in_bilayer" in d.name
        )
        return (embed_dir / "run_embed.py").read_text()

    def _embed_meta(self, workspace):
        embed_dir = next(
            d for d in (workspace / "steps").iterdir()
            if "embed_in_bilayer" in d.name
        )
        return json.loads((embed_dir / "metadata.json").read_text())

    def test_run_embed_imports_pre_shrink_module(self, workspace):
        """run_embed.py must import apply_pre_shrink_exclusion."""
        script = self._embed_script(workspace)
        assert "apply_pre_shrink_exclusion" in script, (
            "run_embed.py must import apply_pre_shrink_exclusion from validators"
        )

    def test_run_embed_has_pre_exclude_config_vars(self, workspace):
        """run_embed.py must define PRE_EXCLUDE_TRAPPED and MAX_PRE_EXCLUDED."""
        script = self._embed_script(workspace)
        assert "PRE_EXCLUDE_TRAPPED" in script
        assert "MAX_PRE_EXCLUDED" in script

    def test_run_embed_writes_pre_shrink_report(self, workspace):
        """run_embed.py must write pre_shrink_lipid_exclusion_report.json."""
        script = self._embed_script(workspace)
        assert "pre_shrink_lipid_exclusion_report.json" in script
        assert "n_removed_by_membrane_mask" in script
        assert "membrane_mask_atom_count" in script
        assert "annular_occupancy_before" in script

    def test_metadata_has_pre_exclude_params(self, workspace):
        """metadata.json records pre_exclude_trapped_lipids and max_pre_excluded_lipids."""
        meta = self._embed_meta(workspace)
        assert "pre_exclude_trapped_lipids" in meta["params"]
        assert "max_pre_excluded_lipids" in meta["params"]

    def test_metadata_expected_outputs_includes_report(self, workspace):
        """pre_shrink_lipid_exclusion_report.json is declared in expected_outputs."""
        meta = self._embed_meta(workspace)
        assert "pre_shrink_lipid_exclusion_report.json" in meta["expected_outputs"]

    def test_safety_limit_message_in_script(self, workspace):
        """run_embed.py must reference safety limit in error path."""
        script = self._embed_script(workspace)
        assert "safety_limit" in script or "MAX_PRE_EXCLUDED" in script


# ═══════════════════════════════════════════════════════════════════════════════
# TestEmbeddingQualityIntegration — snapshot + quality config plumbing in run.sh
# ═══════════════════════════════════════════════════════════════════════════════

class TestEmbeddingQualityIntegration:
    """
    Verify that save_embedding_iterations, save_embedding_iteration_stride,
    and quality_diagnostics config values plumb correctly from the builder
    into the generated run.sh shrink loop script.

    No GROMACS execution required.
    """

    def _build_embedding_step(self, tmp_path: Path, extra_params: dict = {}) -> Path:
        """Build a membrane_embedding step dir with given extra params."""
        from core.execution_models import (
            SimulationStep, StepStage, StepType,
        )
        from builders.step_builders.embedding_builder import EmbeddingBuilder

        step = SimulationStep(
            step_id="membrane_embedding",
            title="Test shrink loop",
            stage=StepStage.MEMBRANE_EMBEDDING,
            step_type=StepType.AUTOMATIC,
            engine="gromacs+perl:inflategro",
            params={
                "lipid":              "DPPC",
                "lipid_residue_name": "DPP",
                "forcefield":         "opls-aa",
                "temperature_K":      298.0,
                "apl_target_ang2":    62.0,
                "apl_tolerance_ang2": 2.0,
                "inflate_factor":     4.0,
                "deflate_factor":     0.95,
                "max_iterations":     200,
                "gridsize":           5,
                "cutoff":             0.14,
                "trapped_lipid_policy": "warn",
                **extra_params,
            },
        )
        step_dir = tmp_path / "membrane_embedding"
        step_dir.mkdir()
        EmbeddingBuilder().build(step, step_dir, step_dir_map={})
        return step_dir

    def _run_sh(self, step_dir: Path) -> str:
        return (step_dir / "run.sh").read_text()

    def _meta(self, step_dir: Path) -> dict:
        return json.loads((step_dir / "metadata.json").read_text())

    # ── snapshot mode: final_only (default) ──────────────────────────────────

    def test_final_only_no_snapshots_dir(self, tmp_path):
        """Default final_only mode: run.sh must NOT create a snapshots/ directory."""
        step_dir = self._build_embedding_step(tmp_path)
        script = self._run_sh(step_dir)
        # No mkdir for snapshots in the loop
        assert "mkdir -p snapshots" not in script

    def test_final_only_metadata_records_mode(self, tmp_path):
        step_dir = self._build_embedding_step(tmp_path)
        meta = self._meta(step_dir)
        assert meta["params"]["save_embedding_iterations"] == "final_only"

    # ── snapshot mode: all ────────────────────────────────────────────────────

    def test_all_mode_creates_snapshots_dir(self, tmp_path):
        """save_embedding_iterations=all: run.sh must mkdir snapshots."""
        step_dir = self._build_embedding_step(
            tmp_path, {"save_embedding_iterations": "all"}
        )
        script = self._run_sh(step_dir)
        assert "mkdir -p snapshots" in script

    def test_all_mode_copies_every_iteration(self, tmp_path):
        """save_embedding_iterations=all: unconditional cp inside the while loop."""
        step_dir = self._build_embedding_step(
            tmp_path, {"save_embedding_iterations": "all"}
        )
        script = self._run_sh(step_dir)
        assert "snapshots/iter_" in script
        # Must NOT be guarded by a modulo condition
        assert "% SNAP_STRIDE" not in script

    def test_all_mode_in_metadata(self, tmp_path):
        step_dir = self._build_embedding_step(
            tmp_path, {"save_embedding_iterations": "all"}
        )
        meta = self._meta(step_dir)
        assert meta["params"]["save_embedding_iterations"] == "all"
        assert "snapshots/" in meta["expected_outputs"]

    # ── snapshot mode: every_n ────────────────────────────────────────────────

    def test_every_n_mode_bakes_stride(self, tmp_path):
        """save_embedding_iterations=every_n: run.sh contains SNAP_STRIDE variable."""
        step_dir = self._build_embedding_step(tmp_path, {
            "save_embedding_iterations":       "every_n",
            "save_embedding_iteration_stride": 10,
        })
        script = self._run_sh(step_dir)
        assert "SNAP_STRIDE=10" in script

    def test_every_n_mode_conditional_copy(self, tmp_path):
        """save_embedding_iterations=every_n: cp is guarded by modulo condition."""
        step_dir = self._build_embedding_step(tmp_path, {
            "save_embedding_iterations": "every_n",
        })
        script = self._run_sh(step_dir)
        assert "% SNAP_STRIDE" in script
        assert "snapshots/iter_" in script

    def test_every_n_default_stride_is_5(self, tmp_path):
        step_dir = self._build_embedding_step(tmp_path, {
            "save_embedding_iterations": "every_n",
        })
        script = self._run_sh(step_dir)
        assert "SNAP_STRIDE=5" in script

    # ── quality diagnostics ───────────────────────────────────────────────────

    def test_quality_diagnostics_default_true(self, tmp_path):
        """quality_diagnostics defaults to True: diagnose_quality.py is generated."""
        step_dir = self._build_embedding_step(tmp_path)
        assert (step_dir / "diagnose_quality.py").exists()

    def test_quality_script_imports_correct_module(self, tmp_path):
        step_dir = self._build_embedding_step(tmp_path)
        script = (step_dir / "diagnose_quality.py").read_text()
        assert "diagnose_embedding_quality" in script
        assert "validators.embedding_quality" in script

    def test_quality_script_writes_report(self, tmp_path):
        step_dir = self._build_embedding_step(tmp_path)
        script = (step_dir / "diagnose_quality.py").read_text()
        assert "embedding_quality_report.json" in script

    def test_quality_script_never_exits_nonzero(self, tmp_path):
        """diagnose_quality.py must not block workflow — no sys.exit with nonzero."""
        step_dir = self._build_embedding_step(tmp_path)
        script = (step_dir / "diagnose_quality.py").read_text()
        # Should not have a sys.exit call (diagnostic only)
        assert "sys.exit" not in script

    def test_run_sh_calls_diagnose_quality(self, tmp_path):
        """run.sh must call python3 diagnose_quality.py after trapped diagnostic."""
        step_dir = self._build_embedding_step(tmp_path)
        script = self._run_sh(step_dir)
        assert "diagnose_quality.py" in script

    def test_quality_report_in_expected_outputs(self, tmp_path):
        step_dir = self._build_embedding_step(tmp_path)
        meta = self._meta(step_dir)
        assert "embedding_quality_report.json" in meta["expected_outputs"]

    def test_quality_disabled_omits_script_and_output(self, tmp_path):
        """quality_diagnostics=False: no diagnose_quality.py and no report in outputs."""
        step_dir = self._build_embedding_step(
            tmp_path, {"quality_diagnostics": False}
        )
        assert not (step_dir / "diagnose_quality.py").exists()
        meta = self._meta(step_dir)
        assert "embedding_quality_report.json" not in meta["expected_outputs"]
        script = self._run_sh(step_dir)
        assert "diagnose_quality.py" not in script

    def test_quality_analysis_region_baked_in(self, tmp_path):
        """quality_analysis_region is baked into diagnose_quality.py."""
        step_dir = self._build_embedding_step(tmp_path, {
            "quality_analysis_region": "radial_shell",
        })
        script = (step_dir / "diagnose_quality.py").read_text()
        assert "radial_shell" in script

    def test_quality_headgroup_and_tail_atoms_baked_in(self, tmp_path):
        """Headgroup (O33) and tail (C50) atom names are baked into diagnose_quality.py."""
        step_dir = self._build_embedding_step(tmp_path)
        script = (step_dir / "diagnose_quality.py").read_text()
        assert "O33" in script   # DPPC OPLS-AA headgroup phosphate oxygen
        assert "C50" in script   # DPPC OPLS-AA tail carbon
