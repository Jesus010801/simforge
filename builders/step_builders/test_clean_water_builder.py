"""
Tests for the clean_water assembly step.

Covers:
  - Builder generates run.sh + run_clean_water.py
  - metadata.json marks automation_level=automated
  - clean_water_report.json declared in expected_outputs
  - SOL count update logic (regex contract)
  - No protein/lipid atoms deleted by WaterDeletorAdapter
  - Dry-run executor treats clean_water as executable (not skipped)
  - Water gate reads clean_water_report.json (new) and water_report.json (compat)
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from core.compiler import SimulationCompiler
from builders.workspace_builder import WorkspaceBuilder
from executors.shell_executor import ShellExecutor
from executors.execution_state import StepStatus

MEMBRANE_YAML      = "configs/membrane_test.yaml"
MEMBRANE_ORIENT_YAML = "configs/membrane_orient_test.yaml"


# ── Shared fixtures ───────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def membrane_workspace(tmp_path_factory):
    result = SimulationCompiler().compile(MEMBRANE_YAML)
    tmp = tmp_path_factory.mktemp("clean_water_ws")
    return WorkspaceBuilder().build(result, output_dir=str(tmp))


@pytest.fixture(scope="module")
def orient_workspace(tmp_path_factory):
    """Full-automation workspace (orient_protein=AUTOMATED) so DAG is not blocked."""
    result = SimulationCompiler().compile(MEMBRANE_ORIENT_YAML)
    tmp = tmp_path_factory.mktemp("clean_water_orient_ws")
    return WorkspaceBuilder().build(result, output_dir=str(tmp))


@pytest.fixture(scope="module")
def clean_water_dir(membrane_workspace):
    steps = membrane_workspace / "steps"
    d = next((p for p in sorted(steps.iterdir()) if "clean_water" in p.name), None)
    assert d is not None, "clean_water step directory not found in compiled workspace"
    return d


@pytest.fixture(scope="module")
def clean_water_meta(clean_water_dir):
    return json.loads((clean_water_dir / "metadata.json").read_text())


@pytest.fixture(scope="module")
def dry_run_records(orient_workspace):
    # Uses the orient workspace so orient_protein is AUTOMATED and the DAG
    # flows through without blocking any downstream step.
    execution = ShellExecutor(orient_workspace, dry_run=True).run()
    return {r.step_id: r for r in execution.steps}


# ── GRO helpers for unit tests ─────────────────────────────────────────────────

def _make_gro(path: Path, atoms: list[tuple[str, str, str, float, float, float]]) -> Path:
    """
    Write a minimal GRO file to *path* (treated as a file, not a directory).
    atoms: list of (resnum_str, resname, atomname, x, y, z) in nm.
    """
    lines = ["test system", str(len(atoms))]
    for i, (resnum, resname, atomname, x, y, z) in enumerate(atoms, start=1):
        lines.append(
            f"{resnum:>5}{resname:<5}{atomname:>5}{i:5d}"
            f"{x:8.3f}{y:8.3f}{z:8.3f}"
        )
    lines.append("10.0  10.0  10.0")
    path.write_text("\n".join(lines) + "\n")
    return path


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Builder outputs
# ═══════════════════════════════════════════════════════════════════════════════

class TestCleanWaterBuilderOutputs:
    def test_run_sh_exists(self, clean_water_dir):
        assert (clean_water_dir / "run.sh").exists(), (
            "AssemblyBuilder must generate run.sh for clean_water"
        )

    def test_run_sh_delegates_to_python_helper(self, clean_water_dir):
        content = (clean_water_dir / "run.sh").read_text()
        assert "run_clean_water.py" in content, (
            "run.sh must delegate execution to run_clean_water.py"
        )

    def test_python_helper_exists(self, clean_water_dir):
        assert (clean_water_dir / "run_clean_water.py").exists(), (
            "AssemblyBuilder must generate run_clean_water.py for clean_water"
        )

    def test_python_helper_imports_adapter(self, clean_water_dir):
        content = (clean_water_dir / "run_clean_water.py").read_text()
        assert "WaterDeletorAdapter" in content, (
            "run_clean_water.py must import and use WaterDeletorAdapter"
        )

    def test_python_helper_writes_clean_water_report(self, clean_water_dir):
        content = (clean_water_dir / "run_clean_water.py").read_text()
        assert "clean_water_report.json" in content, (
            "run_clean_water.py must write clean_water_report.json"
        )

    def test_python_helper_counts_input_sol(self, clean_water_dir):
        content = (clean_water_dir / "run_clean_water.py").read_text()
        assert "input_water_count" in content, (
            "run_clean_water.py must count input SOL molecules for the report"
        )

    def test_python_helper_tracks_topology_updated(self, clean_water_dir):
        content = (clean_water_dir / "run_clean_water.py").read_text()
        assert "topology_updated" in content, (
            "run_clean_water.py must track whether topology update succeeded"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Metadata contract
# ═══════════════════════════════════════════════════════════════════════════════

class TestCleanWaterMetadata:
    def test_automation_level_is_automated(self, clean_water_meta):
        assert clean_water_meta.get("automation_level") == "automated", (
            "clean_water metadata.json must have automation_level=automated"
        )

    def test_step_type_present_for_backward_compat(self, clean_water_meta):
        # Backward compat: old runtimes that read step_type must not treat it as manual
        step_type = clean_water_meta.get("step_type", "automatic")
        assert step_type not in ("manual", "external", "validation"), (
            f"step_type='{step_type}' would cause old runtimes to skip clean_water"
        )

    def test_clean_water_report_in_expected_outputs(self, clean_water_meta):
        outputs = clean_water_meta.get("expected_outputs", [])
        assert "clean_water_report.json" in outputs, (
            "metadata expected_outputs must declare clean_water_report.json"
        )

    def test_system_clean_gro_in_expected_outputs(self, clean_water_meta):
        outputs = clean_water_meta.get("expected_outputs", [])
        assert "system_clean.gro" in outputs, (
            "metadata expected_outputs must declare system_clean.gro"
        )

    def test_topol_in_expected_outputs(self, clean_water_meta):
        outputs = clean_water_meta.get("expected_outputs", [])
        assert "topol.top" in outputs, (
            "metadata expected_outputs must declare topol.top "
            "(SOL count is updated here)"
        )

    def test_water_report_in_expected_outputs_for_gate_compat(self, clean_water_meta):
        outputs = clean_water_meta.get("expected_outputs", [])
        assert "water_report.json" in outputs, (
            "metadata expected_outputs must still declare water_report.json "
            "for gate backward compatibility"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# 3. Dry-run: clean_water must NOT be skipped
# ═══════════════════════════════════════════════════════════════════════════════

class TestCleanWaterDryRun:
    def test_dry_run_does_not_skip_clean_water(self, dry_run_records):
        record = dry_run_records.get("clean_water")
        assert record is not None, "clean_water step not found in dry-run records"
        assert record.status != StepStatus.SKIPPED, (
            f"clean_water must not be SKIPPED in dry-run — got status={record.status}. "
            "Check that automation_level=automated is set in metadata.json."
        )

    def test_dry_run_clean_water_completes(self, dry_run_records):
        record = dry_run_records.get("clean_water")
        assert record is not None
        assert record.status == StepStatus.DONE, (
            f"clean_water dry-run status should be DONE, got {record.status}"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# 4. SOL count update logic
# ═══════════════════════════════════════════════════════════════════════════════

class TestSolCountUpdateLogic:
    """Unit tests for the regex-based SOL count update embedded in run_clean_water.py."""

    def _update_sol_count(self, text: str, n_removed: int) -> str:
        lines = text.splitlines()
        out = []
        for line in lines:
            m = re.match(r'^(SOL)\s+(\d+)', line)
            if m:
                old = int(m.group(2))
                new = old - n_removed
                line = f"SOL              {new}"
            out.append(line)
        return "\n".join(out)

    def test_sol_count_decremented_by_removed(self):
        topol = "DPPC             512\nSOL              1500\nNA               5\n"
        result = self._update_sol_count(topol, 45)
        assert "SOL              1455" in result

    def test_non_sol_lines_unchanged(self):
        topol = "DPPC             512\nSOL              1000\nNA               10\n"
        result = self._update_sol_count(topol, 100)
        assert "DPPC             512" in result
        assert "NA               10" in result

    def test_zero_removed_leaves_count_unchanged(self):
        topol = "SOL              800\n"
        result = self._update_sol_count(topol, 0)
        assert "SOL              800" in result

    def test_multiple_sol_lines_all_updated(self):
        # Edge case: two SOL entries (shouldn't happen in practice, but must not crash)
        topol = "SOL              500\nSOL              300\n"
        result = self._update_sol_count(topol, 10)
        assert "SOL              490" in result
        assert "SOL              290" in result


# ═══════════════════════════════════════════════════════════════════════════════
# 5. Adapter: protein and lipid atoms are not deleted
# ═══════════════════════════════════════════════════════════════════════════════

class TestAdapterPreservesNonWaterAtoms:
    """
    WaterDeletorAdapter must only remove SOL residues.
    Protein and lipid atoms must survive deletion regardless of their Z position.
    """

    def _write_bilayer_gro(self, path: Path) -> Path:
        """
        Minimal bilayer GRO:
          - 2 DPPC lipids: ref_atom O33 at z=3.0 (top) and z=1.0 (bot)
          - 2 DPPC lipids: middle_atom C50 at z=2.0 (midplane)
          - 1 LYS protein atom at z=2.0 (inside bilayer zone — must NOT be deleted)
          - 3 SOL molecules: 1 outside (z=0.1), 2 inside bilayer zone (z=1.5, z=2.5)
        """
        atoms = [
            # DPPC lipids — headgroup ref (O33)
            ("1",  "DPPC", "O33",  5.0, 5.0, 3.0),   # top leaflet
            ("2",  "DPPC", "O33",  5.0, 5.0, 1.0),   # bot leaflet
            # DPPC lipids — tail middle (C50) at midplane
            ("3",  "DPPC", "C50",  5.0, 5.0, 2.0),
            ("4",  "DPPC", "C50",  5.0, 5.0, 2.0),
            # Protein atom inside bilayer zone (must survive)
            ("5",  "LYS",  "CA",   5.0, 5.0, 2.0),
            # SOL outside bilayer (z < z_bot=1.0) — must survive
            ("6",  "SOL",  "OW",   5.0, 5.0, 0.1),
            ("6",  "SOL",  "HW1",  5.1, 5.0, 0.1),
            ("6",  "SOL",  "HW2",  5.0, 5.1, 0.1),
            # SOL inside bilayer (z_bot <= z <= z_top) — must be deleted
            ("7",  "SOL",  "OW",   5.0, 5.0, 1.5),
            ("7",  "SOL",  "HW1",  5.1, 5.0, 1.5),
            ("7",  "SOL",  "HW2",  5.0, 5.1, 1.5),
            ("8",  "SOL",  "OW",   5.0, 5.0, 2.5),
            ("8",  "SOL",  "HW1",  5.1, 5.0, 2.5),
            ("8",  "SOL",  "HW2",  5.0, 5.1, 2.5),
        ]
        return _make_gro(path, atoms)

    def test_protein_atom_survives_deletion(self, tmp_path):
        from adapters.water_deletor_adapter import WaterDeletorAdapter
        gro_in  = self._write_bilayer_gro(tmp_path / "in.gro")
        gro_out = tmp_path / "out.gro"
        result  = WaterDeletorAdapter().run(
            gro_in=gro_in, gro_out=gro_out,
            ref_atom="O33", middle_atom="C50",
        )
        assert result.success
        out_text = gro_out.read_text()
        assert "LYS" in out_text, (
            "WaterDeletorAdapter deleted a protein (LYS) atom — only SOL must be removed"
        )

    def test_lipid_atoms_survive_deletion(self, tmp_path):
        from adapters.water_deletor_adapter import WaterDeletorAdapter
        gro_in  = self._write_bilayer_gro(tmp_path / "in.gro")
        gro_out = tmp_path / "out.gro"
        result  = WaterDeletorAdapter().run(
            gro_in=gro_in, gro_out=gro_out,
            ref_atom="O33", middle_atom="C50",
        )
        assert result.success
        out_text = gro_out.read_text()
        assert "DPPC" in out_text, (
            "WaterDeletorAdapter deleted DPPC (lipid) atoms — only SOL must be removed"
        )

    def test_only_bilayer_sol_is_removed(self, tmp_path):
        from adapters.water_deletor_adapter import WaterDeletorAdapter
        gro_in  = self._write_bilayer_gro(tmp_path / "in.gro")
        gro_out = tmp_path / "out.gro"
        result  = WaterDeletorAdapter().run(
            gro_in=gro_in, gro_out=gro_out,
            ref_atom="O33", middle_atom="C50",
        )
        assert result.success
        assert result.metadata["waters_removed"] == 2, (
            "Expected exactly 2 SOL molecules removed (those inside bilayer zone), "
            f"got {result.metadata['waters_removed']}"
        )

    def test_outside_sol_survives(self, tmp_path):
        from adapters.water_deletor_adapter import WaterDeletorAdapter
        gro_in  = self._write_bilayer_gro(tmp_path / "in.gro")
        gro_out = tmp_path / "out.gro"
        result  = WaterDeletorAdapter().run(
            gro_in=gro_in, gro_out=gro_out,
            ref_atom="O33", middle_atom="C50",
        )
        assert result.success
        lines = gro_out.read_text().splitlines()
        sol_lines = [l for l in lines if "SOL" in l and "OW" in l]
        assert len(sol_lines) == 1, (
            f"SOL outside bilayer must survive deletion; got {len(sol_lines)} OW lines"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# 6. Water gate: report priority
# ═══════════════════════════════════════════════════════════════════════════════

class TestWaterGateReportPriority:
    """Gate reads clean_water_report.json first; falls back to water_report.json."""

    def _write(self, path: Path, data: dict) -> None:
        path.write_text(json.dumps(data))

    def test_gate_absent_when_no_report(self, tmp_path):
        from runtime.water_gate import evaluate_water_gate
        assert evaluate_water_gate(tmp_path) is None

    def test_gate_reads_clean_water_report(self, tmp_path):
        from runtime.water_gate import evaluate_water_gate
        self._write(tmp_path / "clean_water_report.json", {
            "input_water_count":   100,
            "removed_water_count": 10,
            "final_water_count":   0,
            "cutoff_used":         {"z_bot_nm": 1.0, "z_top_nm": 3.0},
            "output_gro_path":     "/tmp/system_clean.gro",
            "topology_updated":    True,
        })
        result = evaluate_water_gate(tmp_path)
        assert result is not None
        assert result.passed

    def test_gate_passes_when_final_count_zero(self, tmp_path):
        from runtime.water_gate import evaluate_water_gate
        self._write(tmp_path / "clean_water_report.json", {
            "final_water_count": 0,
            "topology_updated":  True,
        })
        result = evaluate_water_gate(tmp_path)
        assert result.passed
        assert not result.blocked

    def test_gate_warns_when_few_remain(self, tmp_path):
        from runtime.water_gate import evaluate_water_gate
        self._write(tmp_path / "clean_water_report.json", {
            "final_water_count": 3,
            "topology_updated":  True,
        })
        result = evaluate_water_gate(tmp_path)
        assert not result.passed
        assert not result.blocked
        assert len(result.warnings) > 0

    def test_gate_blocks_when_many_remain(self, tmp_path):
        from runtime.water_gate import evaluate_water_gate
        self._write(tmp_path / "clean_water_report.json", {
            "final_water_count": 20,
            "topology_updated":  True,
        })
        result = evaluate_water_gate(tmp_path)
        assert result.blocked
        assert len(result.errors) > 0

    def test_gate_falls_back_to_water_report_when_primary_absent(self, tmp_path):
        from runtime.water_gate import evaluate_water_gate
        self._write(tmp_path / "water_report.json", {
            "passed":             True,
            "n_waters_remaining": 0,
            "errors":             [],
            "warnings":           [],
            "confidence":         1.0,
        })
        result = evaluate_water_gate(tmp_path)
        assert result is not None
        assert result.passed

    def test_gate_prefers_clean_water_report_over_water_report(self, tmp_path):
        from runtime.water_gate import evaluate_water_gate
        # Primary says 0 waters remaining (pass), legacy says 20 (block)
        self._write(tmp_path / "clean_water_report.json", {
            "final_water_count": 0,
            "topology_updated":  True,
        })
        self._write(tmp_path / "water_report.json", {
            "passed":             False,
            "n_waters_remaining": 20,
            "errors":             ["many waters remain"],
            "warnings":           [],
            "confidence":         1.0,
        })
        result = evaluate_water_gate(tmp_path)
        assert result.passed, (
            "Gate must use clean_water_report.json (pass) rather than "
            "water_report.json (block) when both are present"
        )

    def test_gate_handles_corrupted_primary_gracefully(self, tmp_path):
        from runtime.water_gate import evaluate_water_gate
        (tmp_path / "clean_water_report.json").write_text("NOT JSON {{{")
        result = evaluate_water_gate(tmp_path)
        assert result is None
