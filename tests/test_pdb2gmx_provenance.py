"""
tests/test_pdb2gmx_provenance.py
─────────────────────────────────
Prove that the repair plan (docs/audits/protein_membrane_pipeline_architecture_audit.md)
is enforced in every generated workspace:

1. No generated script calls pdb2gmx on a mixed protein-lipid GRO (work.gro /
   embed_in_bilayer/system.gro / membrane_embedding/system.gro).
2. assemble_system_topology runs AFTER membrane_embedding in the step ordering.
3. solvate_membrane consumes the same converged.gro evaluated by Phase 9A.
4. assemble_system_topology reads converged.gro (post-shrink), not system.gro (pre-shrink).
5. membrane_embedding/run.sh calls run_assemble_system.py for bootstrap topology
   rather than using a pre-assembled topology from before the shrink loop.
"""

import re
import json
import pytest
from pathlib import Path

# ── Shared fixture ──────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def workspace(tmp_path_factory):
    """Build a membrane pipeline workspace and return the root path."""
    import os
    from core.compiler import SimulationCompiler
    from builders.workspace_builder import WorkspaceBuilder

    cfg_path = Path(__file__).resolve().parent.parent / "configs" / "membrane_orient_test.yaml"
    if not cfg_path.exists():
        pytest.skip(f"Membrane test config not found: {cfg_path}")

    root = tmp_path_factory.mktemp("pdb2gmx_provenance")
    orig_cwd = os.getcwd()
    os.chdir(cfg_path.parent)
    try:
        result = SimulationCompiler().compile(str(cfg_path))
    finally:
        os.chdir(orig_cwd)
    WorkspaceBuilder().build(result, workspace_path=root)
    return root


def _step_dir(workspace: Path, step_id: str) -> Path:
    steps = workspace / "steps"
    matches = [d for d in steps.iterdir() if step_id in d.name]
    assert matches, f"Step directory for '{step_id}' not found in {steps}"
    return matches[0]


def _all_generated_scripts(workspace: Path):
    """Yield (path, text) for every generated .sh and .py in the workspace steps/."""
    for p in (workspace / "steps").rglob("*.sh"):
        yield p, p.read_text()
    for p in (workspace / "steps").rglob("*.py"):
        yield p, p.read_text()


# ═══════════════════════════════════════════════════════════════════════════
# 1. pdb2gmx must never run on a mixed protein-lipid GRO
# ═══════════════════════════════════════════════════════════════════════════

class TestPdb2gmxSafety:
    """pdb2gmx is only allowed on: original protein PDB and lipids_only.gro."""

    FORBIDDEN_PATTERNS = [
        # Explicit mixed-system GRO names passed to -f
        r'"-f",\s+["\']?(?:[^"\']*/)?(work\.gro)["\']?',
        r'"-f",\s+["\']?(?:[^"\']*/)?(system\.gro)["\']?',
        r'"-f",\s+"?work\.gro"?',
        r'"-f",\s+"?system\.gro"?',
        # Bash -f flags
        r'-f\s+work\.gro\b',
        r'-f\s+"[^"]*system\.gro"',
        r'-f\s+\$[A-Z_]*GRO\s',      # -f $VARIABLE_GRO (may be a mixed GRO)
    ]

    ALLOWED_F_ARGS = {
        "lipids_only.gro",           # safe: lipid-only input
        "protein_1.pdb",             # safe: original protein PDB
        "inputs/protein_1.pdb",
        ".pdb",                       # any .pdb is allowed
    }

    def test_no_pdb2gmx_on_work_gro(self, workspace):
        """No generated script may call pdb2gmx with -f work.gro."""
        violations = []
        for path, text in _all_generated_scripts(workspace):
            if "pdb2gmx" not in text:
                continue
            if re.search(r'-f["\s]+(?:["\'])?(?:[^"\']*\/)?work\.gro', text):
                violations.append(str(path))
            # Python list form: "-f", "work.gro"
            f_args = re.findall(r'"-f",\s+"([^"]+)"', text)
            for arg in f_args:
                if "work.gro" in arg and "lipids_only" not in arg:
                    violations.append(f"{path}: -f {arg!r}")
        assert not violations, (
            "pdb2gmx must NEVER run on work.gro (mixed protein-lipid GRO).\n"
            "Violations:\n" + "\n".join(violations)
        )

    def test_no_pdb2gmx_on_embed_or_embed_system_gro(self, workspace):
        """No generated script may call pdb2gmx with -f pointing at embed_in_bilayer/system.gro
        or membrane_embedding/system.gro or work.gro from any subdirectory."""
        violations = []
        forbidden_re = re.compile(
            r'pdb2gmx.*?-f["\s,]+(?:["\'])?'
            r'(?:\.\./[^"\']*(?:embed_in_bilayer|membrane_embedding)/[^"\']*system\.gro'
            r'|work\.gro'
            r')',
            re.DOTALL,
        )
        for path, text in _all_generated_scripts(workspace):
            if "pdb2gmx" not in text:
                continue
            if forbidden_re.search(text):
                violations.append(str(path))
        assert not violations, (
            "pdb2gmx must not run on embed_in_bilayer/system.gro or "
            "membrane_embedding system.gro.\nViolations:\n" + "\n".join(violations)
        )

    def test_pdb2gmx_in_assemble_only_on_lipids_only(self, workspace):
        """In assemble_system_topology, pdb2gmx must only be called with lipids_only.gro."""
        assemble_dir = _step_dir(workspace, "assemble_system_topology")
        script = (assemble_dir / "run_assemble_system.py").read_text()
        f_args = re.findall(r'"-f",\s+(?:str\()?"([^"]+)"', script)
        for arg in f_args:
            assert "lipids_only" in arg or ".pdb" in arg, (
                f"assemble_system_topology pdb2gmx -f argument must be lipids_only.gro "
                f"or a .pdb file, got: {arg!r}"
            )

    def test_tm_aware_script_no_pdb2gmx(self, workspace):
        """run_tm_aware.py must not call pdb2gmx on work.gro or system.gro."""
        embed_dir = _step_dir(workspace, "membrane_embedding")
        tm_script = embed_dir / "run_tm_aware.py"
        if not tm_script.exists():
            pytest.skip("run_tm_aware.py not present (inflategro backend)")
        text = tm_script.read_text()
        # pdb2gmx must not appear at all in the TM-aware backend
        assert "pdb2gmx" not in text, (
            "run_tm_aware.py must not call pdb2gmx — use update_topology_lipid_count_py instead"
        )

    def test_generate_protein_topology_only_uses_pdb(self, workspace):
        """generate_protein_topology must only pass a .pdb to pdb2gmx."""
        gen_dir = _step_dir(workspace, "generate_protein_topology")
        script = (gen_dir / "run_protein_topology.py").read_text()
        f_args = re.findall(r'"-f",\s+(?:str\()?"([^"]+)"', script)
        for arg in f_args:
            assert arg.endswith(".pdb") or "protein" in arg, (
                f"generate_protein_topology pdb2gmx -f must be a .pdb file, got: {arg!r}"
            )


# ═══════════════════════════════════════════════════════════════════════════
# 2. DAG ordering: assemble_system_topology AFTER membrane_embedding
# ═══════════════════════════════════════════════════════════════════════════

class TestDagOrdering:
    def test_assemble_after_membrane_embedding(self, workspace):
        step_dirs = sorted(d.name for d in (workspace / "steps").iterdir())
        assemble_idx  = next(i for i, n in enumerate(step_dirs) if "assemble_system_topology" in n)
        memb_emb_idx  = next(i for i, n in enumerate(step_dirs) if "membrane_embedding" in n)
        assert assemble_idx > memb_emb_idx, (
            "assemble_system_topology must come AFTER membrane_embedding "
            "(final topology from converged.gro, not pre-shrink system.gro)"
        )

    def test_solvate_after_assemble(self, workspace):
        step_dirs = sorted(d.name for d in (workspace / "steps").iterdir())
        assemble_idx = next(i for i, n in enumerate(step_dirs) if "assemble_system_topology" in n)
        solvate_idx  = next(i for i, n in enumerate(step_dirs) if "solvate_membrane" in n)
        assert solvate_idx > assemble_idx, (
            "solvate_membrane must come AFTER assemble_system_topology"
        )

    def test_membrane_embedding_before_assemble(self, workspace):
        step_dirs = sorted(d.name for d in (workspace / "steps").iterdir())
        memb_idx     = next(i for i, n in enumerate(step_dirs) if "membrane_embedding" in n)
        assemble_idx = next(i for i, n in enumerate(step_dirs) if "assemble_system_topology" in n)
        assert memb_idx < assemble_idx, (
            "membrane_embedding must run before assemble_system_topology in the step ordering"
        )


# ═══════════════════════════════════════════════════════════════════════════
# 3. solvate_membrane and Phase 9A evaluate the same file (converged.gro)
# ═══════════════════════════════════════════════════════════════════════════

class TestCoordinateLineage:
    def test_solvate_uses_converged_gro(self, workspace):
        """solvate_membrane/run.sh must reference membrane_embedding/converged.gro."""
        solvate_dir = _step_dir(workspace, "solvate_membrane")
        script = (solvate_dir / "run.sh").read_text()
        assert "converged.gro" in script, (
            "solvate_membrane/run.sh must reference converged.gro from membrane_embedding"
        )
        assert "membrane_embedding" in script or "EMBED_DIR" in script, (
            "solvate_membrane must reference the membrane_embedding directory"
        )

    def test_phase9a_evaluates_converged_gro(self, workspace):
        """run_protein_membrane_interface.py must reference converged.gro."""
        embed_dir = _step_dir(workspace, "membrane_embedding")
        iface_script = embed_dir / "run_protein_membrane_interface.py"
        if not iface_script.exists():
            pytest.skip("run_protein_membrane_interface.py not present")
        text = iface_script.read_text()
        assert "converged.gro" in text, (
            "Phase 9A (run_protein_membrane_interface.py) must evaluate converged.gro"
        )

    def test_solvate_topology_from_assemble(self, workspace):
        """solvate_membrane must use assemble_system_topology/topol.top (post-shrink counts)."""
        solvate_dir = _step_dir(workspace, "solvate_membrane")
        script = (solvate_dir / "run.sh").read_text()
        assert "assemble_system_topology" in script, (
            "solvate_membrane must use assemble_system_topology/topol.top "
            "(which now has final counts from converged.gro)"
        )


# ═══════════════════════════════════════════════════════════════════════════
# 4. assemble_system_topology reads converged.gro (not pre-shrink system.gro)
# ═══════════════════════════════════════════════════════════════════════════

class TestTopologySourceCoordinate:
    def test_assemble_script_tries_converged_gro_first(self, workspace):
        """run_assemble_system.py must try converged.gro before falling back to system.gro."""
        assemble_dir = _step_dir(workspace, "assemble_system_topology")
        script = (assemble_dir / "run_assemble_system.py").read_text()
        assert "converged.gro" in script, (
            "run_assemble_system.py must try EMBED_DIR/converged.gro as primary input "
            "(final post-shrink coordinate file)"
        )

    def test_assemble_embed_dir_is_membrane_embedding(self, workspace):
        """EMBED_DIR in run_assemble_system.py should point to membrane_embedding/."""
        assemble_dir = _step_dir(workspace, "assemble_system_topology")
        script = (assemble_dir / "run_assemble_system.py").read_text()
        assert "membrane_embedding" in script, (
            "run_assemble_system.py EMBED_DIR must point to membrane_embedding/ "
            "(not embed_in_bilayer/)"
        )

    def test_assemble_argv_bootstrap_supported(self, workspace):
        """run_assemble_system.py must support argv[1]=GRO argv[2]=OUT_TOP for bootstrap."""
        assemble_dir = _step_dir(workspace, "assemble_system_topology")
        script = (assemble_dir / "run_assemble_system.py").read_text()
        assert "_gro_override" in script or "sys.argv[1]" in script, (
            "run_assemble_system.py must support argv[1] GRO override for bootstrap mode"
        )
        assert "_top_override" in script or "sys.argv[2]" in script, (
            "run_assemble_system.py must support argv[2] output topology override for bootstrap mode"
        )


# ═══════════════════════════════════════════════════════════════════════════
# 5. membrane_embedding bootstrap-calls run_assemble_system.py
# ═══════════════════════════════════════════════════════════════════════════

class TestEmbeddingBootstrap:
    def test_shrink_loop_bootstrap_call(self, workspace):
        """membrane_embedding/run.sh must call run_assemble_system.py for bootstrap topology."""
        embed_dir = _step_dir(workspace, "membrane_embedding")
        run_sh = embed_dir / "run.sh"
        if not run_sh.exists():
            pytest.skip("run.sh not present (tm_aware backend only)")
        text = run_sh.read_text()
        assert "run_assemble_system.py" in text, (
            "membrane_embedding/run.sh must call run_assemble_system.py "
            "to generate bootstrap_topol.top before the shrink loop"
        )
        assert "bootstrap_topol.top" in text, (
            "membrane_embedding/run.sh must reference bootstrap_topol.top"
        )

    def test_shrink_loop_input_from_embed_in_bilayer(self, workspace):
        """membrane_embedding/run.sh must use embed_in_bilayer/system.gro as input GRO."""
        embed_dir = _step_dir(workspace, "membrane_embedding")
        run_sh = embed_dir / "run.sh"
        if not run_sh.exists():
            pytest.skip("run.sh not present")
        text = run_sh.read_text()
        assert "embed_in_bilayer" in text, (
            "membrane_embedding/run.sh must reference embed_in_bilayer/system.gro as INPUT_GRO"
        )
        assert "system_processed.gro" not in text, (
            "membrane_embedding/run.sh must NOT reference system_processed.gro "
            "(assemble_system_topology runs after embedding, not before)"
        )

    def test_tm_aware_bootstrap_call(self, workspace):
        """run_tm_aware.py must call run_assemble_system.py if TOPOL is bootstrap_topol.top."""
        embed_dir = _step_dir(workspace, "membrane_embedding")
        tm_script = embed_dir / "run_tm_aware.py"
        if not tm_script.exists():
            pytest.skip("run_tm_aware.py not present (inflategro backend)")
        text = tm_script.read_text()
        assert "run_assemble_system.py" in text, (
            "run_tm_aware.py must call run_assemble_system.py for bootstrap topology "
            "(assemble_system_topology runs after embedding)"
        )
        assert "bootstrap_topol.top" in text, (
            "run_tm_aware.py must reference bootstrap_topol.top"
        )
