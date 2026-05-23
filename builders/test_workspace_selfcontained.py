"""
Tests for workspace self-containment and relocatability.

A workspace must be fully self-contained:
  - All external input files are copied into workspace/inputs/ at compile time
  - Scripts reference files via workspace-relative paths (never absolute, never CWD)
  - Moving the workspace to a different directory does not break execution
  - Missing source files fail at compile time, not runtime

Error that motivated these tests:
  gmx pdb2gmx fails because the executor runs from
  simforge_runs/protein/steps/01_prepare_protein_1/
  and the script references 'protein_1.pdb' without copying it to the workspace.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from core.compiler import SimulationCompiler
from builders.workspace_builder import WorkspaceBuilder


PROTEIN_YAML  = "configs/lysozyme_test.yaml"
MEMBRANE_YAML = "configs/membrane_test.yaml"


# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def protein_workspace(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("selfcontained_protein")
    result = SimulationCompiler().compile(PROTEIN_YAML)
    return WorkspaceBuilder().build(result, output_dir=str(tmp))


@pytest.fixture(scope="module")
def relocated_workspace(tmp_path_factory, protein_workspace):
    """Workspace moved to a completely different directory."""
    dest = tmp_path_factory.mktemp("relocated") / "moved_workspace"
    shutil.copytree(protein_workspace, dest)
    return dest


# ── inputs/ directory ─────────────────────────────────────────────────────────

class TestInputsDirectory:
    def test_inputs_dir_created(self, protein_workspace):
        assert (protein_workspace / "inputs").is_dir(), (
            "workspace/inputs/ must exist after compilation"
        )

    def test_inputs_contains_protein_pdb(self, protein_workspace):
        inputs = list((protein_workspace / "inputs").iterdir())
        pdbs = [f for f in inputs if f.suffix in (".pdb", ".gro", ".cif")]
        assert len(pdbs) >= 1, (
            f"workspace/inputs/ must contain at least one structure file; "
            f"found: {[f.name for f in inputs]}"
        )

    def test_inputs_files_are_non_empty(self, protein_workspace):
        for f in (protein_workspace / "inputs").iterdir():
            assert f.stat().st_size > 0, (
                f"Staged input file is empty: {f.name}"
            )

    def test_inputs_named_by_component_id(self, protein_workspace):
        """Files must be named {component_id}.{ext} for unambiguous identification."""
        inputs_dir = protein_workspace / "inputs"
        # For the protein workflow there's one component: protein_1
        stems = {f.stem for f in inputs_dir.iterdir() if f.is_file()}
        assert "protein_1" in stems, (
            f"inputs/ must contain protein_1.pdb; found stems: {stems}"
        )


# ── Protein preparation script ────────────────────────────────────────────────

class TestProteinPrepScript:
    def _get_prep_dir(self, workspace):
        return next(
            (p for p in (workspace / "steps").iterdir()
             if "prepare_protein" in p.name),
            None,
        )

    def test_script_references_inputs_dir(self, protein_workspace):
        prep_dir = self._get_prep_dir(protein_workspace)
        assert prep_dir is not None
        script = (prep_dir / "run.sh").read_text()
        assert "INPUTS_DIR" in script, (
            "prepare_protein/run.sh must set INPUTS_DIR pointing to workspace/inputs/"
        )

    def test_script_uses_inputs_dir_for_pdb(self, protein_workspace):
        prep_dir = self._get_prep_dir(protein_workspace)
        script = (prep_dir / "run.sh").read_text()
        # gmx pdb2gmx -f "$INPUTS_DIR/..." not just -f protein_1.pdb
        assert '-f "$INPUTS_DIR/' in script, (
            "gmx pdb2gmx must reference the PDB via $INPUTS_DIR, "
            "not a bare filename (which requires CWD = project root)"
        )

    def test_inputs_dir_is_relative_not_absolute(self, protein_workspace):
        """INPUTS_DIR must be a relative path so the workspace is relocatable."""
        import re
        prep_dir = self._get_prep_dir(protein_workspace)
        for sh in prep_dir.glob("*.sh"):
            for line in sh.read_text().splitlines():
                m = re.match(r'INPUTS_DIR="([^"]+)"', line.strip())
                if m:
                    inputs_val = m.group(1)
                    assert not inputs_val.startswith("/"), (
                        f"{sh.name}: INPUTS_DIR must be a relative path "
                        f"(got absolute: {inputs_val!r})"
                    )

    def test_no_hardcoded_absolute_variable_assignment(self, protein_workspace):
        """No DIR or PATH variable may be assigned an absolute path outside workspace."""
        import re
        workspace_str = str(protein_workspace.resolve())
        prep_dir = self._get_prep_dir(protein_workspace)
        assign_pat = re.compile(r'^[A-Z_]+=["\'](/[^"\']+)["\']')
        for sh in prep_dir.glob("*.sh"):
            for line in sh.read_text().splitlines():
                m = assign_pat.match(line.strip())
                if m:
                    abs_path = m.group(1)
                    assert abs_path.startswith(workspace_str), (
                        f"{sh.name}: variable assigned absolute external path: "
                        f"{line.strip()!r}"
                    )

    def test_required_inputs_lists_inputs_dir_pdb(self, protein_workspace):
        prep_dir = self._get_prep_dir(protein_workspace)
        meta = json.loads((prep_dir / "metadata.json").read_text())
        required = meta.get("required_inputs", [])
        assert any("inputs" in r and ".pdb" in r for r in required), (
            f"prepare_protein metadata.json must list inputs/protein_1.pdb "
            f"in required_inputs; got {required}"
        )

    def test_required_inputs_resolve_within_workspace(self, protein_workspace):
        """All required_inputs paths must exist relative to their step_dir."""
        for step_dir in (protein_workspace / "steps").iterdir():
            if not step_dir.is_dir():
                continue
            meta_path = step_dir / "metadata.json"
            if not meta_path.exists():
                continue
            meta = json.loads(meta_path.read_text())
            for rel_path in meta.get("required_inputs", []):
                resolved = (step_dir / rel_path).resolve()
                # Skip paths that go into other step dirs (checked at runtime)
                # Only check inputs/ paths here
                if "inputs" in rel_path:
                    assert resolved.exists(), (
                        f"{step_dir.name}: required_input '{rel_path}' "
                        f"does not exist at {resolved}"
                    )


# ── Relocatability ────────────────────────────────────────────────────────────

class TestRelocatability:
    def test_scripts_work_after_relocation(self, relocated_workspace):
        """All workspace-internal relative paths must resolve after moving."""
        for step_dir in (relocated_workspace / "steps").iterdir():
            if not step_dir.is_dir():
                continue
            meta_path = step_dir / "metadata.json"
            if not meta_path.exists():
                continue
            meta = json.loads(meta_path.read_text())
            for rel_path in meta.get("required_inputs", []):
                # Only check inputs/ refs (inter-step refs aren't staged yet)
                if "inputs" not in rel_path:
                    continue
                resolved = (step_dir / rel_path).resolve()
                assert resolved.exists(), (
                    f"After relocation, {step_dir.name}: "
                    f"required_input '{rel_path}' broken at {resolved}"
                )

    def test_inputs_dir_preserved_after_relocation(self, relocated_workspace):
        inputs_dir = relocated_workspace / "inputs"
        assert inputs_dir.is_dir()
        assert any(inputs_dir.iterdir()), "inputs/ must not be empty after relocation"

    def test_no_original_paths_in_scripts(self, protein_workspace, relocated_workspace):
        """Scripts in the relocated workspace must not reference the original path."""
        original_path = str(protein_workspace.resolve())
        for step_dir in (relocated_workspace / "steps").iterdir():
            if not step_dir.is_dir():
                continue
            for sh in step_dir.glob("*.sh"):
                content = sh.read_text()
                assert original_path not in content, (
                    f"{step_dir.name}/{sh.name} contains hardcoded original path: "
                    f"{original_path}"
                )


# ── Compile-time failure for missing source files ─────────────────────────────

class TestMissingSourceFileFails:
    def test_missing_pdb_fails_at_compile_time(self, tmp_path):
        """
        If a component's source file doesn't exist, WorkspaceBuilder must raise
        FileNotFoundError at compile time — before any scripts are generated.
        """
        import yaml

        yaml_content = {
            "project": {"name": "test_missing"},
            "components": [
                {"id": "protein_1", "role": "protein", "file": str(tmp_path / "nonexistent.pdb")}
            ],
            "environment": {
                "solvent": {"water_model": "spce"},
                "ions": {"concentration": 0.15},
                "temperature_K": 300.0,
                "duration_ns": 0.1,
            },
            "forcefields": {"protein": "opls-aa"},
            "simulation_objectives": ["stability"],
        }

        yaml_path = tmp_path / "test.yaml"
        yaml_path.write_text(yaml.dump(yaml_content))

        result = SimulationCompiler().compile(str(yaml_path))

        with pytest.raises(FileNotFoundError, match="nonexistent.pdb"):
            WorkspaceBuilder().build(result, output_dir=str(tmp_path / "runs"))
