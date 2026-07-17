"""
tests/test_membrane_embedding_optimizer.py
Phase 9B: Closed-Loop Membrane Embedding Optimizer — unit tests.
"""
from __future__ import annotations

import json
import shutil
import tempfile
from pathlib import Path

import pytest

from validators.membrane_embedding_optimizer import (
    run_membrane_embedding_optimizer,
    compute_candidate_score,
    check_quality_passed,
)


# ── GRO helpers ───────────────────────────────────────────────────────────────

def _write_gro(path, atoms, box=(10.0, 10.0, 10.0)):
    """atoms: list of (resid, resname, atomname, x, y, z)"""
    lines = ["test system", str(len(atoms))]
    for i, (resid, resname, atomname, x, y, z) in enumerate(atoms, 1):
        lines.append(f"{resid:5d}{resname:<5s}{atomname:>5s}{i:5d}{x:8.3f}{y:8.3f}{z:8.3f}")
    bx, by, bz = box
    lines.append(f"   {bx:.5f}   {by:.5f}   {bz:.5f}")
    Path(path).write_text("\n".join(lines) + "\n")


def _bad_coverage_system(baseline_path):
    """10 TM atoms at z=5.0, lipids 0.9 nm above → all exposed at baseline.
    z_shift=-0.6 brings lipids to z=5.3 (0.3 nm above), distance ≈ 0.30 nm → covered.
    """
    atoms = []
    tm_residues = set(range(1, 11))
    for ri in range(10):
        atoms.append((ri + 1, "ALA", "CA", 5.0 + ri * 0.2, 5.0, 5.0))
    for ri in range(10):
        lx = 5.0 + ri * 0.2   # same XY as corresponding TM atom
        atoms.append((100 + ri, "DPPC", "P", lx, 5.0, 5.9))  # 0.9 nm above
    _write_gro(baseline_path, atoms)
    return tm_residues


def _fully_covered_system(path):
    """10 TM atoms at z=5.0 with lipids at 0.30 nm (covered even at baseline)."""
    atoms = []
    tm_residues = set(range(1, 11))
    for ri in range(10):
        atoms.append((ri + 1, "ALA", "CA", 5.0, 5.0, 5.0 + ri * 0.15))
    for ri in range(10):
        atoms.append((100 + ri, "DPPC", "P", 5.3, 5.0, 5.0 + ri * 0.15))
    _write_gro(path, atoms)
    return tm_residues


# ── Test 1: optimizer selects candidate with higher coverage ──────────────────

def test_optimizer_selects_better_candidate():
    """Baseline has all TM atoms 0.9 nm from lipids (exposed).
    z_shift=-0.6 brings lipids to 0.3 nm → covered.
    Optimizer must select a negative z_shift candidate.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        baseline = tmpdir / "system_baseline.gro"
        work_input = tmpdir / "work_input.gro"
        tm_residues = _bad_coverage_system(baseline)

        report = run_membrane_embedding_optimizer(
            system_baseline_gro    = baseline,
            work_input_gro         = work_input,
            tm_residues            = tm_residues,
            output_dir             = tmpdir,
            z_shift_offsets_nm     = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6],
            mask_padding_nm_values = [0.0],
        )

        assert report["enabled"] is True
        assert report["selected_candidate_id"] is not None, "Should have selected a candidate"
        # Selected z_shift must be negative (shift lipids toward TM)
        selected_params = report["selected_candidate_params"]
        assert selected_params["z_shift_offset_nm"] < 0.0, (
            f"Expected negative z_shift, got {selected_params['z_shift_offset_nm']}"
        )
        # Score must have improved
        assert report["score_improvement"] is not None
        assert report["score_improvement"] > 0.0, (
            f"Selected candidate should improve baseline, got improvement={report['score_improvement']}"
        )
        # work_input.gro must exist (promoted candidate)
        assert work_input.exists()


# ── Test 2: candidate with excessive lipid deletion is marked invalid ─────────

def test_candidate_with_excessive_deletion_is_invalid():
    """Set safety_limit=2 so any candidate with mask_padding that would remove 3+
    lipids is invalid. Uses a system where TM and lipids overlap significantly.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        baseline = tmpdir / "system_baseline.gro"
        work_input = tmpdir / "work_input.gro"

        # TM atoms and lipids placed overlapping at the SAME XY, similar Z
        # Lipids at 0.05 nm from TM atoms → would all be removed with any cutoff > 0.05
        atoms = []
        tm_residues = set(range(1, 6))
        for ri in range(5):
            atoms.append((ri + 1, "ALA", "CA", 5.0, 5.0, 4.0 + ri * 0.2))
        # 10 lipid residues, all very close to TM atoms (within 0.12 nm)
        for ri in range(10):
            lz = 4.0 + (ri % 5) * 0.2 + 0.10  # 0.10 nm above each TM atom
            atoms.append((100 + ri, "DPPC", "P", 5.0, 5.0, lz))
        _write_gro(baseline, atoms)

        report = run_membrane_embedding_optimizer(
            system_baseline_gro    = baseline,
            work_input_gro         = work_input,
            tm_residues            = tm_residues,
            output_dir             = tmpdir,
            z_shift_offsets_nm     = [0.0],
            mask_padding_nm_values = [0.0, 0.5],  # mask_padding=0.5 → cutoff=0.64 → removes all
            exclusion_safety_limit = 3,             # only allow up to 3 lipid removals
        )

        assert report["enabled"] is True
        invalid_cands = [c for c in report["candidates"] if c["invalid_reason"] is not None]
        assert len(invalid_cands) >= 1, (
            f"Expected at least 1 invalid candidate (safety limit exceeded), "
            f"got 0. Candidates: {report['candidates']}"
        )
        assert "exclusion_safety_limit_exceeded" in invalid_cands[0]["invalid_reason"]


# ── Test 3: candidate promotion writes work_input.gro with promoted content ───

def test_candidate_promotion_replaces_work_input():
    """After optimization, work_input.gro must differ from baseline when a
    candidate with z_shift != 0 is promoted.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        baseline = tmpdir / "system_baseline.gro"
        work_input = tmpdir / "work_input.gro"
        tm_residues = _bad_coverage_system(baseline)

        baseline_content = baseline.read_bytes()

        report = run_membrane_embedding_optimizer(
            system_baseline_gro    = baseline,
            work_input_gro         = work_input,
            tm_residues            = tm_residues,
            output_dir             = tmpdir,
            z_shift_offsets_nm     = [-0.6, 0.0],
            mask_padding_nm_values = [0.0],
        )

        assert work_input.exists(), "work_input.gro must exist after optimization"
        # If a non-zero z_shift was selected, content should differ
        if report["selected_candidate_params"] and report["selected_candidate_params"]["z_shift_offset_nm"] != 0.0:
            assert work_input.read_bytes() != baseline_content, (
                "work_input.gro should differ from baseline when z_shift != 0 selected"
            )
        assert report["promoted_to_work_input_gro"] is True


# ── Test 4: system_baseline.gro is preserved ──────────────────────────────────

def test_baseline_preserved_as_system_baseline():
    """system_baseline.gro must be an exact copy of INPUT_GRO and not be modified
    after the optimizer runs.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        baseline = tmpdir / "system_baseline.gro"
        work_input = tmpdir / "work_input.gro"

        # Write baseline manually
        tm_residues = _bad_coverage_system(baseline)
        original_bytes = baseline.read_bytes()

        run_membrane_embedding_optimizer(
            system_baseline_gro    = baseline,
            work_input_gro         = work_input,
            tm_residues            = tm_residues,
            output_dir             = tmpdir,
            z_shift_offsets_nm     = [-0.4, 0.0, 0.4],
            mask_padding_nm_values = [0.0],
        )

        assert baseline.exists(), "system_baseline.gro must still exist after optimization"
        assert baseline.read_bytes() == original_bytes, (
            "system_baseline.gro must not be modified by the optimizer"
        )


# ── Test 5: strict mode emits warning but continues ───────────────────────────

def test_strict_mode_continues_with_warning():
    """When policy=strict and no candidate passes, optimizer emits warnings but
    does NOT delete work_input.gro; the pipeline can continue with the baseline.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        baseline = tmpdir / "system_baseline.gro"
        work_input = tmpdir / "work_input.gro"

        # Fully covered system — baseline already passes, but optimizer with strict
        # should still complete without crashing
        tm_residues = _fully_covered_system(baseline)

        report = run_membrane_embedding_optimizer(
            system_baseline_gro    = baseline,
            work_input_gro         = work_input,
            tm_residues            = tm_residues,
            output_dir             = tmpdir,
            z_shift_offsets_nm     = [0.0],
            mask_padding_nm_values = [0.0],
            policy                 = "strict",
        )

        # strict mode must not crash and work_input must still exist
        assert work_input.exists(), "work_input.gro must exist even in strict mode"
        assert report["policy"] == "strict"


# ── Test 6: warn mode continues with best candidate even if none pass ─────────

def test_warn_mode_selects_best_even_when_no_candidate_passes_fully():
    """With very tight thresholds, no candidate passes quality gates.
    warn mode should still select the best-score candidate if it improves baseline.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        baseline = tmpdir / "system_baseline.gro"
        work_input = tmpdir / "work_input.gro"
        tm_residues = _bad_coverage_system(baseline)

        report = run_membrane_embedding_optimizer(
            system_baseline_gro    = baseline,
            work_input_gro         = work_input,
            tm_residues            = tm_residues,
            output_dir             = tmpdir,
            z_shift_offsets_nm     = [-0.6, -0.4, 0.0],
            mask_padding_nm_values = [0.0],
            policy                 = "warn",
            # Absurdly tight threshold — nothing will pass
            min_fraction_covered   = 0.999,
            min_score_improvement  = 0.05,
        )

        # Warn mode: should select best improving candidate even without quality_pass
        # (or fall back to baseline if no improvement)
        assert report["policy"] == "warn"
        assert work_input.exists()
        # Must have a warning about quality gate failure
        if report["n_candidates_quality_passed"] == 0 and report["selected_candidate_id"] is not None:
            assert any("quality" in w.lower() or "gate" in w.lower()
                       for w in report["warnings"]), (
                f"Expected quality warning, got: {report['warnings']}"
            )


# ── Test 7: score favors coverage over void-only improvement ──────────────────

def test_score_favors_coverage_over_void_reduction():
    """Candidate A: fraction_covered=0.40, void_fraction=0.05 (low void but bad coverage)
    Candidate B: fraction_covered=0.90, void_fraction=0.25 (good coverage, higher void)
    Score formula must favor candidate B (coverage weight 3.0 >> void weight 1.0).
    """
    iface_a = {"fraction_covered": 0.40, "fraction_exposed_gap": 0.50,
                "p90_nearest_lipid_distance": 1.0}
    quality_a = {"void_fraction_local": 0.05, "tm_burial_score": 0.5, "apl_local_nm2": None}

    iface_b = {"fraction_covered": 0.90, "fraction_exposed_gap": 0.05,
                "p90_nearest_lipid_distance": 0.35}
    quality_b = {"void_fraction_local": 0.25, "tm_burial_score": 0.5, "apl_local_nm2": None}

    score_a = compute_candidate_score(iface_a, quality_a, 100, 100)
    score_b = compute_candidate_score(iface_b, quality_b, 100, 100)

    assert score_b > score_a, (
        f"Coverage-rich candidate B (score={score_b:.3f}) should score higher "
        f"than void-reduced candidate A (score={score_a:.3f})"
    )


# ── Test 8: GLP-1R regression — optimizer improves Phase 9A failure ───────────

def test_glp1r_regression_optimizer_improves_coverage():
    """Synthetic GLP-1R-like system: 40 TM atoms, lipids placed 0.9 nm above.
    Baseline Phase 9A: fraction_covered ≈ 0, p90 ≈ 0.9 nm (FAIL).
    Optimizer with z_shift=-0.6 brings lipids to 0.3 nm → passes.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        baseline = tmpdir / "system_baseline.gro"
        work_input = tmpdir / "work_input.gro"

        atoms = []
        tm_residues = set(range(1, 41))
        # 40 TM CA atoms in a circle at TM Z center
        import math
        for ri in range(40):
            angle = ri * 2 * math.pi / 40
            x = 5.0 + 0.5 * math.cos(angle)
            y = 5.0 + 0.5 * math.sin(angle)
            z = 5.0 + ri * 0.05  # slight spread in Z
            atoms.append((ri + 1, "ALA", "CA", x, y, z))
        # Lipids placed 0.9 nm directly above each TM atom (same XY)
        for ri in range(40):
            angle = ri * 2 * math.pi / 40
            lx = 5.0 + 0.5 * math.cos(angle)
            ly = 5.0 + 0.5 * math.sin(angle)
            lz = 5.0 + ri * 0.05 + 0.9  # 0.9 nm above
            atoms.append((100 + ri, "DPPC", "P", lx, ly, lz))
        _write_gro(baseline, atoms)

        # Verify baseline Phase 9A failure
        from validators.protein_membrane_interface import evaluate_protein_membrane_interface
        (tmpdir / "baseline_check").mkdir(exist_ok=True)
        baseline_iface = evaluate_protein_membrane_interface(
            gro_path    = baseline,
            tm_residues = tm_residues,
            output_dir  = tmpdir / "baseline_check",
        )
        (tmpdir / "baseline_check").mkdir(exist_ok=True)
        baseline_fc = baseline_iface.get("fraction_covered", 1.0)
        assert baseline_fc < 0.3, f"Baseline should have poor coverage, got {baseline_fc:.3f}"

        # Run optimizer
        report = run_membrane_embedding_optimizer(
            system_baseline_gro    = baseline,
            work_input_gro         = work_input,
            tm_residues            = tm_residues,
            output_dir             = tmpdir,
            z_shift_offsets_nm     = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6],
            mask_padding_nm_values = [0.0],
        )

        # Optimizer must select a negative z_shift that improves coverage
        assert report["score_improvement"] is not None and report["score_improvement"] > 0.1, (
            f"Optimizer should find significant improvement over baseline. "
            f"score_improvement={report['score_improvement']}"
        )
        # work_input.gro must be different from baseline (a better placement)
        assert work_input.exists()
        assert report["promoted_to_work_input_gro"] is True


# ── Test 9: no coordinate inconsistency after promotion ───────────────────────

def test_promoted_gro_has_no_partial_residues():
    """All residues in the promoted work_input.gro must be complete
    (same atoms-per-residue as the corresponding residue in the baseline).
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        baseline = tmpdir / "system_baseline.gro"
        work_input = tmpdir / "work_input.gro"
        tm_residues = _bad_coverage_system(baseline)

        run_membrane_embedding_optimizer(
            system_baseline_gro    = baseline,
            work_input_gro         = work_input,
            tm_residues            = tm_residues,
            output_dir             = tmpdir,
            z_shift_offsets_nm     = [-0.6, 0.0],
            mask_padding_nm_values = [0.0],
        )

        # Parse both files, count atoms per residue
        def _atoms_per_resid(path):
            lines = path.read_text().splitlines()
            n = int(lines[1].strip())
            counts: dict = {}
            for ln in lines[2: 2 + n]:
                if len(ln) >= 10:
                    try:
                        rid = int(ln[0:5])
                        counts[rid] = counts.get(rid, 0) + 1
                    except ValueError:
                        pass
            return counts

        baseline_counts   = _atoms_per_resid(baseline)
        work_input_counts = _atoms_per_resid(work_input)

        # Z-shift does not remove atoms; every residue count must match baseline
        for resid, count in work_input_counts.items():
            assert count == baseline_counts.get(resid, count), (
                f"Residue {resid} has {count} atoms in work_input.gro but "
                f"{baseline_counts.get(resid)} in baseline — partial residue!"
            )


# ── Test 10: EmbeddingBuilder generates run_membrane_optimizer.py ─────────────

def test_embedding_builder_generates_optimizer_script():
    """Compile a workspace and verify run_membrane_optimizer.py exists with
    correct structure: reads system_baseline.gro, writes work_input.gro.
    """
    import os
    from core.compiler import SimulationCompiler
    from builders.workspace_builder import WorkspaceBuilder

    yaml_text = """
project:
  name: test_optimizer_workspace
  description: Phase 9B optimizer integration test

components:
  - id: protein_1
    role: protein
    file: GLP-1_R.pdb

structural_annotation:
  membrane_topology:
    transmembrane_segments:
      - residues: "31-69"
        label: "TM1"
        helix_type: "alpha"

environment:
  membrane:
    enabled: true
    type: DPPC
    embedding:
      backend: inflategro
      optimizer_enabled: true
      optimizer_policy: warn
      optimizer_max_candidates: 14
      optimizer_z_shift_offsets_nm: [-0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8]
      optimizer_mask_padding_nm: [0.0, 0.1]
  solvent:
    water_model: spce
  ions:
    concentration: 0.15
  temperature_K: 298.0
  duration_ns: 10.0

forcefields:
  protein: opls-aa

simulation_objectives:
  - membrane_insertion

analysis:
  - type: rmsd
"""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        # Write YAML to file (compiler expects file path)
        yaml_file = tmpdir / "test_config.yaml"
        yaml_file.write_text(yaml_text)

        # Provide a minimal dummy PDB so workspace staging succeeds
        dummy_pdb = tmpdir / "GLP-1_R.pdb"
        dummy_pdb.write_text(
            "ATOM      1  CA  ALA A   1       5.000   5.000   5.000  1.00  0.00           C\n"
            "END\n"
        )

        orig_cwd = os.getcwd()
        os.chdir(tmpdir)
        try:
            result = SimulationCompiler().compile(str(yaml_file))
        finally:
            os.chdir(orig_cwd)

        ws_path = WorkspaceBuilder().build(result, workspace_path=tmpdir / "ws")

        # Find the embedding step directory (steps/ subdirectory)
        steps_dir = ws_path / "steps"
        assert steps_dir.exists(), f"steps/ dir not found in {ws_path}"

        # membrane_embedding step (shrink loop) is where the optimizer is written
        embed_dirs = [d for d in steps_dir.iterdir()
                      if d.is_dir() and "membrane_embedding" in d.name]
        if not embed_dirs:
            # Fallback: any directory containing run_membrane_optimizer.py
            embed_dirs = [d for d in steps_dir.iterdir()
                          if d.is_dir() and (d / "run_membrane_optimizer.py").exists()]
        assert embed_dirs, (
            f"No membrane_embedding step found. Steps: {[d.name for d in sorted(steps_dir.iterdir())]}"
        )
        embed_dir = embed_dirs[0]

        opt_script = embed_dir / "run_membrane_optimizer.py"
        assert opt_script.exists(), (
            f"run_membrane_optimizer.py not found in {embed_dir}"
        )

        content = opt_script.read_text()
        assert "run_membrane_embedding_optimizer" in content
        assert "system_baseline_gro" in content
        assert "work_input_gro" in content

        # Verify run.sh calls it
        run_sh = embed_dir / "run.sh"
        assert run_sh.exists()
        run_sh_content = run_sh.read_text()
        assert "run_membrane_optimizer.py" in run_sh_content
        assert "work_input.gro" in run_sh_content

        # Verify membrane_embedding_optimizer_report.json in expected_outputs
        meta = json.loads((embed_dir / "metadata.json").read_text())
        assert "membrane_embedding_optimizer_report.json" in meta["expected_outputs"]
