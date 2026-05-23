"""Integration tests for runtime/scientific_summary.py — no GROMACS required."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from runtime.scientific_summary import generate_summary, ScientificSummary


# ─────────────────────────────────────────────────────────────────────────────
# Workspace fixture helpers
# ─────────────────────────────────────────────────────────────────────────────

def _minimal_workspace(tmp_path: Path) -> Path:
    """Workspace with no XVG files → empty best-effort summary."""
    (tmp_path / "steps").mkdir()
    (tmp_path / "metadata").mkdir()
    return tmp_path


def _workspace_with_xvg(tmp_path: Path, content: str, name: str = "rmsd.xvg",
                         subdir: str = "steps/01_analysis") -> Path:
    """Workspace with one XVG file in a step subdirectory."""
    (tmp_path / "metadata").mkdir(exist_ok=True)
    step_dir = tmp_path / subdir
    step_dir.mkdir(parents=True, exist_ok=True)
    (step_dir / name).write_text(content)
    return tmp_path


_RMSD_XVG = """\
# RMSD
@ title "RMSD"
@ xaxis label "Time (ps)"
@ yaxis label "RMSD (nm)"
@ s0 legend "Backbone"
    0.0  0.000
  100.0  0.100
  200.0  0.150
  300.0  0.155
  400.0  0.152
  500.0  0.153
  600.0  0.151
  700.0  0.154
  800.0  0.152
  900.0  0.153
 1000.0  0.154
"""

_ENERGY_XVG = """\
# Energy
@ title "Potential Energy"
@ xaxis label "Time (ps)"
@ yaxis label "Energy (kJ/mol)"
@ s0 legend "Potential"
    0.0  -50000.0
  100.0  -50010.0
  200.0  -50005.0
  300.0  -50008.0
  400.0  -50003.0
"""

_MALFORMED_XVG = """\
@ title "Bad"
this is not a number line
also bad
0.0 0.1
1.0 0.2
"""


# ─────────────────────────────────────────────────────────────────────────────
# Tests
# ─────────────────────────────────────────────────────────────────────────────

class TestGenerateSummaryEmpty:
    def test_empty_workspace_no_crash(self, tmp_path):
        ws  = _minimal_workspace(tmp_path)
        res = generate_summary(ws)
        assert isinstance(res, ScientificSummary)

    def test_empty_workspace_no_analyses(self, tmp_path):
        ws  = _minimal_workspace(tmp_path)
        res = generate_summary(ws)
        assert res.analyses == []

    def test_empty_workspace_not_converged(self, tmp_path):
        """No XVG → best-effort, converged=False."""
        ws  = _minimal_workspace(tmp_path)
        res = generate_summary(ws)
        assert res.converged is False

    def test_empty_workspace_no_runtime_ns(self, tmp_path):
        ws  = _minimal_workspace(tmp_path)
        res = generate_summary(ws)
        assert res.runtime_ns is None

    def test_nonexistent_workspace_no_crash(self, tmp_path):
        """Workspace that does not exist → empty summary, no exception."""
        ws  = tmp_path / "nonexistent"
        res = generate_summary(ws)
        assert isinstance(res, ScientificSummary)
        assert res.analyses == []


class TestGenerateSummaryWithRMSD:
    def test_rmsd_xvg_detected(self, tmp_path):
        ws  = _workspace_with_xvg(tmp_path, _RMSD_XVG, "rmsd.xvg")
        res = generate_summary(ws)
        assert len(res.analyses) == 1
        assert res.analyses[0]["kind"] == "rmsd"

    def test_rmsd_verdict_populated(self, tmp_path):
        ws  = _workspace_with_xvg(tmp_path, _RMSD_XVG, "rmsd.xvg")
        res = generate_summary(ws)
        assert res.rmsd_verdict != ""

    def test_runtime_ns_estimated_from_xvg(self, tmp_path):
        ws  = _workspace_with_xvg(tmp_path, _RMSD_XVG, "rmsd.xvg")
        res = generate_summary(ws)
        # XVG time axis ends at 1000 ps = 1 ns
        assert res.runtime_ns == pytest.approx(1.0)

    def test_workspace_stored_on_result(self, tmp_path):
        ws  = _workspace_with_xvg(tmp_path, _RMSD_XVG, "rmsd.xvg")
        res = generate_summary(ws)
        assert res.workspace == ws


class TestGenerateSummaryWithEnergy:
    def test_energy_xvg_detected(self, tmp_path):
        ws  = _workspace_with_xvg(tmp_path, _ENERGY_XVG, "energy.xvg")
        res = generate_summary(ws)
        assert len(res.analyses) == 1
        assert res.analyses[0]["kind"] == "energy"

    def test_energy_verdict_populated(self, tmp_path):
        ws  = _workspace_with_xvg(tmp_path, _ENERGY_XVG, "energy.xvg")
        res = generate_summary(ws)
        assert res.energy_verdict != ""


class TestGenerateSummaryMalformed:
    def test_malformed_xvg_skipped_or_warned(self, tmp_path):
        ws  = _workspace_with_xvg(tmp_path, _MALFORMED_XVG, "bad.xvg")
        # Should not crash; may or may not add a warning depending on parse outcome
        res = generate_summary(ws)
        assert isinstance(res, ScientificSummary)

    def test_malformed_plus_valid_both_processed(self, tmp_path):
        """A malformed file doesn't block processing of a valid file."""
        step_dir = tmp_path / "steps" / "01_analysis"
        step_dir.mkdir(parents=True)
        (tmp_path / "metadata").mkdir(exist_ok=True)
        (step_dir / "bad.xvg").write_text(_MALFORMED_XVG)
        (step_dir / "rmsd.xvg").write_text(_RMSD_XVG)
        res = generate_summary(tmp_path)
        # At least the valid RMSD file should produce an analysis entry
        kinds = [a["kind"] for a in res.analyses]
        assert "rmsd" in kinds


class TestSummaryOutputFormats:
    def test_as_dict_serializable(self, tmp_path):
        ws  = _workspace_with_xvg(tmp_path, _RMSD_XVG)
        res = generate_summary(ws)
        d   = res.as_dict()
        # Must be JSON-serializable
        serialized = json.dumps(d, default=str)
        parsed     = json.loads(serialized)
        assert "converged" in parsed
        assert "analyses" in parsed

    def test_as_markdown_contains_header(self, tmp_path):
        ws  = _workspace_with_xvg(tmp_path, _RMSD_XVG)
        res = generate_summary(ws)
        md  = res.as_markdown()
        assert "# Scientific Summary" in md

    def test_as_markdown_contains_rmsd_section(self, tmp_path):
        ws  = _workspace_with_xvg(tmp_path, _RMSD_XVG)
        res = generate_summary(ws)
        md  = res.as_markdown()
        assert "RMSD" in md

    def test_as_dict_has_expected_keys(self, tmp_path):
        ws  = _minimal_workspace(tmp_path)
        res = generate_summary(ws)
        d   = res.as_dict()
        for key in ("workspace", "converged", "rmsd_verdict", "energy_verdict",
                    "runtime_ns", "n_analyses", "analyses", "warnings"):
            assert key in d, f"Missing key: {key}"

    def test_write_summary_json_to_workspace(self, tmp_path):
        """generate_summary can have its dict written to JSON without error."""
        ws  = _workspace_with_xvg(tmp_path, _RMSD_XVG)
        res = generate_summary(ws)
        out = tmp_path / "metadata" / "scientific_summary.json"
        out.write_text(json.dumps(res.as_dict(), default=str, indent=2))
        assert out.exists()
        parsed = json.loads(out.read_text())
        assert "converged" in parsed

    def test_write_summary_md_to_workspace(self, tmp_path):
        ws  = _workspace_with_xvg(tmp_path, _RMSD_XVG)
        res = generate_summary(ws)
        out = tmp_path / "metadata" / "scientific_summary.md"
        out.write_text(res.as_markdown())
        assert out.exists()
        assert "Scientific Summary" in out.read_text()


class TestExecutionBackend:
    """Sprint 1b — ExecutionBackend abstraction."""

    def test_local_backend_importable(self):
        from runtime.execution_backend import LocalSubprocessBackend, ExecutionBackend
        backend = LocalSubprocessBackend()
        assert isinstance(backend, ExecutionBackend)

    def test_remediation_executor_accepts_backend(self, tmp_path):
        from runtime.execution_backend import LocalSubprocessBackend
        from executors.remediation_executor import RemediationExecutor
        backend = LocalSubprocessBackend()
        rem = RemediationExecutor(tmp_path, dry_run=True, backend=backend)
        assert rem._backend is backend

    def test_remediation_executor_no_backend_is_none(self, tmp_path):
        from executors.remediation_executor import RemediationExecutor
        rem = RemediationExecutor(tmp_path, dry_run=True)
        assert rem._backend is None

    def test_runtime_executor_has_backend(self, tmp_path):
        """RuntimeExecutor creates a LocalSubprocessBackend on init."""
        from runtime.executor import RuntimeExecutor
        from runtime.execution_backend import LocalSubprocessBackend
        # create a minimal workspace so __init__ doesn't fail
        (tmp_path / "metadata").mkdir()
        (tmp_path / "metadata" / "execution_manifest.json").write_text(
            json.dumps({"system_type": "test", "steps": []})
        )
        executor = RuntimeExecutor(tmp_path, dry_run=True)
        assert isinstance(executor._backend, LocalSubprocessBackend)
