import pytest
import shutil
import json
import subprocess
import sys
from pathlib import Path
from core.compiler import SimulationCompiler
from builders.workspace_builder import WorkspaceBuilder

_HAS_GMX = bool(shutil.which("gmx"))

class TestEmbeddingBackend:

    def test_invalid_backend_fails_validation(self, tmp_path):
        yaml_content = """
project:
  name: test_invalid
components:
  - id: peptide_1
    role: protein
    file: peptide_A.pdb
structural_annotation:
  membrane_topology:
    extracellular_regions: ["1-10"]
    intracellular_regions: ["40-50"]
    transmembrane_segments:
      - residues: "11-39"
        label: "TM1"
        helix_type: "alpha"
environment:
  membrane:
    enabled: true
membrane:
  embedding:
    backend: "invalid_backend"
forcefields:
  protein: opls-aa
simulation_objectives:
  - membrane_insertion
"""
        yaml_file = tmp_path / "invalid.yaml"
        yaml_file.write_text(yaml_content)
        (tmp_path / "peptide_A.pdb").write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00")

        with pytest.raises(ValueError, match="Invalid membrane.embedding.backend"):
            import os
            orig = os.getcwd()
            os.chdir(tmp_path)
            try:
                SimulationCompiler().compile(str(yaml_file))
            finally:
                os.chdir(orig)

    def test_inflategro_backend_generation(self, tmp_path):
        yaml_content = """
project:
  name: test_inflategro
components:
  - id: peptide_1
    role: protein
    file: peptide_A.pdb
structural_annotation:
  membrane_topology:
    extracellular_regions: ["1-10"]
    intracellular_regions: ["40-50"]
    transmembrane_segments:
      - residues: "11-39"
        label: "TM1"
        helix_type: "alpha"
environment:
  membrane:
    enabled: true
membrane:
  embedding:
    backend: "inflategro"
forcefields:
  protein: opls-aa
simulation_objectives:
  - membrane_insertion
"""
        yaml_file = tmp_path / "inflategro.yaml"
        yaml_file.write_text(yaml_content)
        (tmp_path / "peptide_A.pdb").write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00")

        import os
        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            state = SimulationCompiler().compile(str(yaml_file))
        finally:
            os.chdir(orig)

        ws_dir = tmp_path / "ws"
        WorkspaceBuilder().build(state, workspace_path=ws_dir, yaml_source=str(yaml_file))

        embed_step_dir = next((ws_dir / "steps").glob("*_membrane_embedding"), None)
        assert embed_step_dir is not None and embed_step_dir.exists()
        metadata = json.loads((embed_step_dir / "metadata.json").read_text())
        assert metadata["engine"] == "gromacs+perl:inflategro"
        assert metadata["gate"]["type"] == "apl_report"

        run_sh = (embed_step_dir / "run.sh").read_text()
        assert "inflategro-Jorge.pl" in run_sh
        assert not (embed_step_dir / "run_tm_aware.py").exists()

    def test_default_inflategro_cutoff_matches_protmemfiles_tutorial(self, tmp_path):
        """
        Regression for docs/audits/simforge_vs_protmemfiles_audit.md Fix 1.

        The default lipid-exclusion cutoff passed to inflategro-Jorge.pl must
        match the original protmemfiles tutorial command exactly:

            docs/Prot-Memb_FILES/tutorial_membrana.txt:36
            perl inflategro-Jorge.pl system.gro 4 DPP 14 system_inflated.gro 5 area.dat

        The "14" argument becomes 1.4 nm internally ($cutoff = $ARGV[3]*0.1 in
        inflategro.pl). A previous default of cutoff_nm=0.14 was 10x too small
        and left protein-clashing lipids essentially unremoved.
        """
        yaml_content = """
project:
  name: test_default_cutoff
components:
  - id: peptide_1
    role: protein
    file: peptide_A.pdb
structural_annotation:
  membrane_topology:
    extracellular_regions: ["1-10"]
    intracellular_regions: ["40-50"]
    transmembrane_segments:
      - residues: "11-39"
        label: "TM1"
        helix_type: "alpha"
environment:
  membrane:
    enabled: true
membrane:
  embedding:
    backend: "inflategro"
forcefields:
  protein: opls-aa
simulation_objectives:
  - membrane_insertion
"""
        yaml_file = tmp_path / "default_cutoff.yaml"
        yaml_file.write_text(yaml_content)
        (tmp_path / "peptide_A.pdb").write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00")

        import os
        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            state = SimulationCompiler().compile(str(yaml_file))
        finally:
            os.chdir(orig)

        ws_dir = tmp_path / "ws"
        WorkspaceBuilder().build(state, workspace_path=ws_dir, yaml_source=str(yaml_file))

        embed_step_dir = next((ws_dir / "steps").glob("*_membrane_embedding"), None)
        assert embed_step_dir is not None and embed_step_dir.exists()

        run_sh = (embed_step_dir / "run.sh").read_text()
        assert "CUTOFF_NM=1.4" in run_sh, (
            "Default inflategro cutoff must be 1.4 nm (pinned to "
            "docs/Prot-Memb_FILES/tutorial_membrana.txt:36), not 0.14 nm"
        )
        assert "CUTOFF_ANGSTROM=14" in run_sh, (
            "Default cutoff in Å units passed to inflategro-Jorge.pl must be 14 "
            "(-> 1.4 nm internally), matching the original tutorial command"
        )

        # Also verify the pipeline-level default (the value actually baked into
        # step.params before EmbeddingBuilder ever sees it — this is the real
        # source of truth for compiled workflows).
        embed_step = next(
            s for s in state.plan.steps if s.step_id == "membrane_embedding"
        )
        assert embed_step.params.get("cutoff") == 1.4

    def test_shrink_loop_uses_zero_cutoff_not_one_shot_cutoff(self, tmp_path):
        """
        Regression for the GLP-1R membrane-void bug (see
        docs/audits/glp1r_membrane_void_shrink_loop_audit.md).

        The original protmemfiles method applies the overlap-removal cutoff
        exactly ONCE (docs/Prot-Memb_FILES/tutorial_membrana.txt:36), then
        every shrink-loop deflate call uses a literal cutoff of 0
        (docs/Prot-Memb_FILES/ScriptCamilo-Jorge.sh:18,38;
        run_inflategro.sh:41,48) — the loop only relaxes/compresses lipid
        packing, it never deletes lipids again.

        SimForge previously reused the same nonzero one-shot cutoff on every
        shrink-loop iteration, which cumulatively deleted lipids from the
        annulus around the protein as the box shrank (observed: 478 -> 400
        DPP over 29 iterations on GLP-1R, even though the global APL
        correctly converged to 63 Å²). This test asserts the generated
        run.sh keeps the one-shot cutoff for the initial call but uses a
        separate, zero-default loop cutoff for every shrink-loop iteration.
        """
        yaml_content = """
project:
  name: test_shrink_loop_cutoff
components:
  - id: peptide_1
    role: protein
    file: peptide_A.pdb
structural_annotation:
  membrane_topology:
    extracellular_regions: ["1-10"]
    intracellular_regions: ["40-50"]
    transmembrane_segments:
      - residues: "11-39"
        label: "TM1"
        helix_type: "alpha"
environment:
  membrane:
    enabled: true
membrane:
  embedding:
    backend: "inflategro"
forcefields:
  protein: opls-aa
simulation_objectives:
  - membrane_insertion
"""
        yaml_file = tmp_path / "shrink_loop_cutoff.yaml"
        yaml_file.write_text(yaml_content)
        (tmp_path / "peptide_A.pdb").write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00")

        import os
        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            state = SimulationCompiler().compile(str(yaml_file))
        finally:
            os.chdir(orig)

        ws_dir = tmp_path / "ws"
        WorkspaceBuilder().build(state, workspace_path=ws_dir, yaml_source=str(yaml_file))

        embed_step_dir = next((ws_dir / "steps").glob("*_membrane_embedding"), None)
        assert embed_step_dir is not None and embed_step_dir.exists()

        run_sh = (embed_step_dir / "run.sh").read_text()

        # One-shot cutoff (initial call, before the loop) is unchanged.
        assert "CUTOFF_NM=1.4" in run_sh
        assert "CUTOFF_ANGSTROM=14" in run_sh

        # New: separate, zero-default loop cutoff.
        assert "LOOP_CUTOFF_NM=0.0" in run_sh, (
            "shrink-loop cutoff must default to 0.0 nm, matching the original "
            "protmemfiles method's cutoff=0 during every deflate iteration"
        )
        assert "LOOP_CUTOFF_ANGSTROM=0" in run_sh

        lines = run_sh.splitlines()
        while_idx = next(i for i, l in enumerate(lines) if l.strip().startswith("while ["))
        done_idx = next(i for i, l in enumerate(lines[while_idx:], start=while_idx) if l.strip() == "done")

        # The initial (Step 1, before the loop) perl call must use the one-shot cutoff.
        initial_perl_idx = next(
            i for i, l in enumerate(lines[:while_idx]) if 'perl "$INFLATEGRO"' in l
        )
        assert "$CUTOFF_ANGSTROM" in lines[initial_perl_idx]
        assert "$LOOP_CUTOFF_ANGSTROM" not in lines[initial_perl_idx]

        # The shrink loop's perl call must use the loop cutoff, NOT the one-shot cutoff.
        loop_perl_idx = next(
            i for i in range(while_idx, done_idx) if 'perl "$INFLATEGRO"' in lines[i]
        )
        assert "$LOOP_CUTOFF_ANGSTROM" in lines[loop_perl_idx], (
            "shrink-loop inflategro call must use $LOOP_CUTOFF_ANGSTROM"
        )
        assert "$CUTOFF_ANGSTROM" not in lines[loop_perl_idx], (
            "shrink-loop inflategro call must not reuse the one-shot $CUTOFF_ANGSTROM "
            f"(note: '$CUTOFF_ANGSTROM' is not a substring of '$LOOP_CUTOFF_ANGSTROM'); "
            f"line was: {lines[loop_perl_idx]!r}"
        )

    def test_generated_script_syncs_topology_before_every_grompp(self, tmp_path):
        """
        Regression: after inflategro removes lipids from work.gro, the
        topology's [ molecules ] lipid count must be synchronized with the
        actual coordinate file before every subsequent grompp call. Without
        this, grompp fails with a hard atom-count mismatch (see
        validators/topology_sync.py docstring for the real-run evidence).

        This asserts the generated run.sh (inflategro backend) contains the
        sync_topology.py helper and calls it before each of the three grompp
        invocations: the post-inflation minimization (Step 2), every
        shrink-loop iteration (Step 3), and the final minimization (Step 5).
        """
        yaml_content = """
project:
  name: test_topology_sync
components:
  - id: peptide_1
    role: protein
    file: peptide_A.pdb
structural_annotation:
  membrane_topology:
    extracellular_regions: ["1-10"]
    intracellular_regions: ["40-50"]
    transmembrane_segments:
      - residues: "11-39"
        label: "TM1"
        helix_type: "alpha"
environment:
  membrane:
    enabled: true
membrane:
  embedding:
    backend: "inflategro"
forcefields:
  protein: opls-aa
simulation_objectives:
  - membrane_insertion
"""
        yaml_file = tmp_path / "topology_sync.yaml"
        yaml_file.write_text(yaml_content)
        (tmp_path / "peptide_A.pdb").write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00")

        import os
        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            state = SimulationCompiler().compile(str(yaml_file))
        finally:
            os.chdir(orig)

        ws_dir = tmp_path / "ws"
        WorkspaceBuilder().build(state, workspace_path=ws_dir, yaml_source=str(yaml_file))

        embed_step_dir = next((ws_dir / "steps").glob("*_membrane_embedding"), None)
        assert embed_step_dir is not None and embed_step_dir.exists()

        # The sync helper script itself must be generated.
        sync_script = embed_step_dir / "sync_topology.py"
        assert sync_script.exists(), "EmbeddingBuilder must generate sync_topology.py"
        sync_text = sync_script.read_text()
        assert "sync_topology_lipid_count" in sync_text
        assert "validators.topology_sync" in sync_text

        run_sh = (embed_step_dir / "run.sh").read_text()
        lines = run_sh.splitlines()

        grompp_indices = [i for i, l in enumerate(lines) if l.strip().startswith("gmx grompp")]
        sync_indices = [i for i, l in enumerate(lines) if "sync_topology.py" in l]

        assert len(grompp_indices) == 3, (
            f"expected 3 grompp calls (Step 2, shrink-loop, Step 5); found {len(grompp_indices)}"
        )
        assert len(sync_indices) >= 3, (
            "expected at least 3 sync_topology.py calls (one before each grompp call); "
            f"found {len(sync_indices)}"
        )

        # Every grompp call must be preceded (not necessarily immediately, but
        # somewhere earlier in the script) by at least one sync call.
        for gi in grompp_indices:
            assert any(si < gi for si in sync_indices), (
                f"grompp call at line {gi} ('{lines[gi].strip()}') has no preceding "
                "sync_topology.py call"
            )

        # The shrink-loop's sync call must be inside the while loop (indented),
        # i.e. it re-syncs on every iteration, not just once before the loop.
        while_idx = next(i for i, l in enumerate(lines) if l.strip().startswith("while ["))
        done_idx = next(i for i, l in enumerate(lines[while_idx:], start=while_idx) if l.strip() == "done")
        loop_sync_indices = [si for si in sync_indices if while_idx < si < done_idx]
        assert loop_sync_indices, (
            "shrink loop must call sync_topology.py on every iteration, not just once"
        )
        loop_grompp_indices = [gi for gi in grompp_indices if while_idx < gi < done_idx]
        assert loop_grompp_indices, "shrink loop must contain its own grompp call"
        assert loop_sync_indices[0] < loop_grompp_indices[0], (
            "sync_topology.py must run before grompp inside the shrink loop"
        )

    def test_tm_aware_backend_generation(self, tmp_path):
        yaml_content = """
project:
  name: test_tm_aware
components:
  - id: peptide_1
    role: protein
    file: peptide_A.pdb
structural_annotation:
  membrane_topology:
    extracellular_regions: ["1-10"]
    intracellular_regions: ["40-50"]
    transmembrane_segments:
      - residues: "11-39"
        label: "TM1"
        helix_type: "alpha"
environment:
  membrane:
    enabled: true
membrane:
  embedding:
    backend: "tm_aware"
forcefields:
  protein: opls-aa
simulation_objectives:
  - membrane_insertion
"""
        yaml_file = tmp_path / "tm_aware.yaml"
        yaml_file.write_text(yaml_content)
        (tmp_path / "peptide_A.pdb").write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00")

        import os
        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            state = SimulationCompiler().compile(str(yaml_file))
        finally:
            os.chdir(orig)

        ws_dir = tmp_path / "ws"
        WorkspaceBuilder().build(state, workspace_path=ws_dir, yaml_source=str(yaml_file))

        embed_step_dir = next((ws_dir / "steps").glob("*_membrane_embedding"), None)
        assert embed_step_dir is not None and embed_step_dir.exists()
        metadata = json.loads((embed_step_dir / "metadata.json").read_text())
        assert metadata["engine"] == "gromacs:tm_aware_embedding"
        assert metadata["gate"]["type"] == "tm_aware_report"

        assert (embed_step_dir / "run_tm_aware.py").exists()
        run_sh = (embed_step_dir / "run.sh").read_text()
        assert "run_tm_aware.py" in run_sh

    def test_tm_aware_missing_annotation_blocks(self, tmp_path):
        yaml_content = """
project:
  name: test_tm_aware_no_tm
components:
  - id: peptide_1
    role: protein
    file: peptide_A.pdb
environment:
  membrane:
    enabled: true
membrane:
  embedding:
    backend: "tm_aware"
forcefields:
  protein: opls-aa
simulation_objectives:
  - membrane_insertion
"""
        yaml_file = tmp_path / "tm_aware_no_tm.yaml"
        yaml_file.write_text(yaml_content)
        (tmp_path / "peptide_A.pdb").write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00")

        import os
        orig = os.getcwd()
        os.chdir(tmp_path)
        try:
            state = SimulationCompiler().compile(str(yaml_file))
        finally:
            os.chdir(orig)

        ws_dir = tmp_path / "ws"
        WorkspaceBuilder().build(state, workspace_path=ws_dir, yaml_source=str(yaml_file))

        embed_step_dir = next((ws_dir / "steps").glob("*_membrane_embedding"), None)
        assert embed_step_dir is not None and embed_step_dir.exists()
        ret = subprocess.run([sys.executable, "run_tm_aware.py"], cwd=embed_step_dir, capture_output=True, text=True)
        assert ret.returncode == 1
        assert "TM-aware backend requires" in ret.stderr

        backend_report = json.loads((embed_step_dir / "tm_aware_backend_report.json").read_text())
        assert backend_report["quality_passed"] is False
        assert "TM-aware backend requires" in backend_report["blocking_reason"]

    @pytest.mark.skipif(not _HAS_GMX, reason="gmx not found on PATH")
    @pytest.mark.skipif(
        not Path("simforge_runs/protein-membrane/runs/2026-06-24_01-38-59/steps/03_embed_in_bilayer/system.gro").exists(),
        reason="GLP-1R test data not found"
    )
    def test_glp1r_tm_aware_backend_execution(self, tmp_path):
        from core.execution_models import SimulationStep, StepStage, StepType
        from builders.step_builders.embedding_builder import EmbeddingBuilder

        step_dir = tmp_path / "embedding"
        step_dir.mkdir()

        glp1r_embed_dir = Path("simforge_runs/protein-membrane/runs/2026-06-24_01-38-59/steps/03_embed_in_bilayer")
        glp1r_topol_dir = Path("simforge_runs/protein-membrane/runs/2026-06-24_01-38-59/steps/04_generate_topology")
        ff_src = Path("simforge_runs/protein-membrane/runs/2026-06-24_01-38-59/membrane_assets/oplsaa_membrane.ff")

        # Skip if topology is in pre-repair format (pdb2gmx on system.gro → only Protein 1
        # in [ molecules ]; no DPP). The new tm_aware backend uses update_topology_lipid_count_py
        # which requires DPP to already be present in [ molecules ] — added by bootstrap topology
        # (run_assemble_system.py). Regenerate test data from a new-workflow run to re-enable.
        _old_topol = glp1r_topol_dir / "topol.top"
        if _old_topol.exists():
            _tail = _old_topol.read_text()[-500:]
            if "\nDPP" not in _tail and "DPP " not in _tail:
                pytest.skip(
                    "Pre-repair topology (DPP not in [ molecules ]): "
                    "regenerate test data from a new-workflow run to re-enable."
                )

        # Set up coordinates & topology
        input_gro = step_dir / "system_processed.gro"
        shutil.copy(glp1r_topol_dir / "system_processed.gro", input_gro)

        topol_top = step_dir / "topol.top"
        shutil.copy(glp1r_topol_dir / "topol.top", topol_top)

        for itp in glp1r_topol_dir.glob("*.itp"):
            shutil.copy(itp, step_dir / itp.name)

        if ff_src.exists():
            shutil.copytree(ff_src, step_dir / "oplsaa_membrane.ff")

        step = SimulationStep(
            step_id="membrane_embedding",
            title="TM-aware embedding",
            stage=StepStage.MEMBRANE_EMBEDDING,
            step_type=StepType.AUTOMATIC,
            engine="gromacs:tm_aware_embedding",
            depends_on=["embed_in_bilayer", "generate_topology"],
            params={
                "backend": "tm_aware",
                "lipid": "DPPC",
                "lipid_residue_name": "DPP",
                "forcefield": "opls-aa",
                "temperature_K": 310.0,
                "input_gro": str(input_gro),
                "topol_top": str(topol_top),
                "tm_residues": "143-172,175-210,214-242,249-272,277-296,307-328,342-365",
                # Thresholds are permissive for this pre-shrink test dataset:
                # the system is post-exclusion before bilayer re-equilibration,
                # so voids are large and burial score reflects the disturbed bilayer.
                "max_void_fraction_local": 0.60,
                "max_void_area_local_nm2": 5.0,
                "min_tm_burial_score": 0.40,
                # Some residual detector disagreement on test data
                "max_trapped_lipids": 10,
                # GPCR inter-TM loop atoms naturally lie in the bilayer Z range
                "max_soluble_domain_core_atoms": 2000,
                "tm_aware_max_removed_lipids": 100,
                "save_embedding_iterations": "final_only",
                "save_embedding_iteration_stride": 5,
                "quality_diagnostics": True,
            }
        )

        # Build step with correct step directory mappings
        step_dir_map = {
            "generate_topology": glp1r_topol_dir,
            "embed_in_bilayer": glp1r_embed_dir
        }
        EmbeddingBuilder().build(step, step_dir, step_dir_map)

        # Run script
        ret = subprocess.run([sys.executable, "run_tm_aware.py"], cwd=step_dir, capture_output=True, text=True)
        assert ret.returncode == 0, f"run_tm_aware.py failed: {ret.stderr}\nSTDOUT:\n{ret.stdout}"

        # Verify output coordinates & reports
        assert (step_dir / "converged.gro").exists()
        assert (step_dir / "tm_mask_report.json").exists()
        assert (step_dir / "tm_lipid_classification_report.json").exists()
        assert (step_dir / "tm_aware_lipid_exclusion_report.json").exists()
        assert (step_dir / "embedding_quality_report.json").exists()
        assert (step_dir / "tm_aware_backend_report.json").exists()

        backend_report = json.loads((step_dir / "tm_aware_backend_report.json").read_text())
        assert backend_report["quality_passed"] is True
        assert backend_report["final_trapped_lipid_count"] <= 10
        assert backend_report["n_lipids_removed"] > 0
        # Pre-shrink system; burial score reflects disturbed bilayer (not converged)
        assert backend_report["final_tm_burial_score"] > 0.0

        # Verify that quality gates block when thresholds are violated
        # Set min_tm_burial_score very high to force a failure
        step_fail = SimulationStep(
            step_id="membrane_embedding",
            title="TM-aware embedding",
            stage=StepStage.MEMBRANE_EMBEDDING,
            step_type=StepType.AUTOMATIC,
            engine="gromacs:tm_aware_embedding",
            depends_on=["embed_in_bilayer", "generate_topology"],
            params={
                "backend": "tm_aware",
                "lipid": "DPPC",
                "lipid_residue_name": "DPP",
                "forcefield": "opls-aa",
                "temperature_K": 310.0,
                "input_gro": str(input_gro),
                "topol_top": str(topol_top),
                "tm_residues": "143-172,175-210,214-242,249-272,277-296,307-328,342-365",
                # All other gates are permissive so only tm_burial_score triggers failure
                "max_void_fraction_local": 0.60,
                "max_void_area_local_nm2": 5.0,
                "min_tm_burial_score": 0.99,  # Intentionally strict — forces failure
                "max_trapped_lipids": 10,
                "max_soluble_domain_core_atoms": 2000,
                "tm_aware_max_removed_lipids": 100,
                "save_embedding_iterations": "final_only",
                "save_embedding_iteration_stride": 5,
                "quality_diagnostics": True,
            }
        )
        step_dir_fail = tmp_path / "embedding_fail"
        step_dir_fail.mkdir()
        shutil.copy(input_gro, step_dir_fail / "system_processed.gro")
        shutil.copy(topol_top, step_dir_fail / "topol.top")
        for itp in glp1r_topol_dir.glob("*.itp"):
            shutil.copy(itp, step_dir_fail / itp.name)
        if ff_src.exists():
            shutil.copytree(ff_src, step_dir_fail / "oplsaa_membrane.ff")

        EmbeddingBuilder().build(step_fail, step_dir_fail, step_dir_map)
        ret_fail = subprocess.run([sys.executable, "run_tm_aware.py"], cwd=step_dir_fail, capture_output=True, text=True)
        assert ret_fail.returncode == 1

        backend_report_fail = json.loads((step_dir_fail / "tm_aware_backend_report.json").read_text())
        assert backend_report_fail["quality_passed"] is False
        assert "tm_burial_score" in backend_report_fail["blocking_reason"]


@pytest.mark.parametrize(
    ("embedding_config", "expected_guard"),
    [
        ("", "if false ; then"),
        ("    annular_repair_enabled: true\n    annular_repair_allow_insertion: false", "if true ; then"),
    ],
    ids=["defaults-disabled", "explicitly-enabled"],
)
def test_annular_repair_compile_and_literal_rendering(
    tmp_path, embedding_config, expected_guard, monkeypatch
):
    """Compile succeeds and annular repair flags render as shell-safe literals."""
    yaml_content = f"""
project:
  name: test_annular_repair
components:
  - id: peptide_1
    role: protein
    file: peptide_A.pdb
structural_annotation:
  membrane_topology:
    extracellular_regions: ["1-10"]
    intracellular_regions: ["40-50"]
    transmembrane_segments:
      - residues: "11-39"
        label: "TM1"
        helix_type: "alpha"
environment:
  membrane:
    enabled: true
membrane:
  embedding:
    backend: inflategro
{embedding_config}
forcefields:
  protein: opls-aa
simulation_objectives:
  - membrane_insertion
"""
    yaml_file = tmp_path / "annular_repair.yaml"
    yaml_file.write_text(yaml_content)
    (tmp_path / "peptide_A.pdb").write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00\n"
    )

    import os
    orig = os.getcwd()
    os.chdir(tmp_path)
    try:
        state = SimulationCompiler().compile(str(yaml_file))
    finally:
        os.chdir(orig)

    ws_dir = tmp_path / "ws"
    WorkspaceBuilder().build(state, workspace_path=ws_dir, yaml_source=str(yaml_file))
    embed_step_dir = next((ws_dir / "steps").glob("*_membrane_embedding"))
    run_sh = (embed_step_dir / "run.sh").read_text()

    assert expected_guard in run_sh
    assert "allow_insertion=False" in run_sh
    assert "set(TM_RESIDUES)" not in run_sh
    assert "[11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39]" in run_sh
    assert "annular_repair_enabled" not in run_sh
    assert "annular_repair_allow_insertion" not in run_sh

    phase = run_sh.split("# ── Phase M2B:", 1)[1].split("# ── Step 7:", 1)[0]
    if expected_guard == "if false ; then":
        completed = subprocess.run(
            ["bash", "-c", phase], cwd=embed_step_dir, capture_output=True, text=True
        )
        assert completed.returncode == 0, completed.stderr
        report = json.loads((embed_step_dir / "annular_repair_report.json").read_text())
        assert report["enabled"] is False
        assert report["status"] == "disabled"
        assert report["coordinate_file_modified"] is False
        assert report["topology_modified"] is False
    else:
        import types
        calls = []
        fake_module = types.ModuleType("validators.annular_repair")
        fake_module.repair_annular_packing = lambda *args, **kwargs: calls.append((args, kwargs))
        monkeypatch.setitem(sys.modules, "validators.annular_repair", fake_module)
        enabled_python = phase.split("python3 - <<'PY'\n", 1)[1].split("\nPY", 1)[0]
        exec(enabled_python, {})
        assert calls
        args, kwargs = calls[0]
        assert args[0] == "converged.gro"
        assert args[1] == list(range(11, 40))
        assert kwargs == {"enabled": True, "allow_insertion": False}
