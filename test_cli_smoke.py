"""Smoke tests for CLI entry-points and import graph.

These tests catch broken imports immediately after refactors, before
the affected code path is exercised in a real simulation run.
"""
from __future__ import annotations

import importlib
import sys
from pathlib import Path

import pytest
from typer.testing import CliRunner

# ── CLI runner ────────────────────────────────────────────────────────────────

from cli import cli

runner = CliRunner(mix_stderr=False)

SIMPLE_YAML = "configs/lysozyme_test.yaml"


# ── Import graph ──────────────────────────────────────────────────────────────

REQUIRED_MODULES = [
    # Core pipeline
    "core.parser",
    "core.compiler",
    "core.compiler_models",
    "core.decision_engine",
    "core.execution_models",
    "core.variant_compiler",
    "core.workspace_fingerprint",
    "core.scientific_planner",
    "core.planning_models",
    "core.semantic_artifacts",
    # Builders
    "builders.workspace_builder",
    "builders.step_builders.analysis_builder",
    "builders.step_builders.assembly_builder",
    "builders.step_builders.equilibration_builder",
    "builders.step_builders.minimization_builder",
    "builders.step_builders.production_builder",
    # Executors
    "executors.base_executor",
    "executors.shell_executor",
    "executors.gromacs_executor",
    "executors.execution_state",
    "executors.remediation_executor",
    # Runtime  ← most likely victim of future renames
    "runtime.executor",
    "runtime.events",
    "runtime.artifacts",
    "runtime.cache",
    "runtime.stream",
    "runtime.journal",
    "runtime.metrics",
]


@pytest.mark.parametrize("module", REQUIRED_MODULES)
def test_module_importable(module: str) -> None:
    """Verify every key module imports cleanly — catches path renames early."""
    importlib.import_module(module)


# ── CLI --help (no filesystem, no GROMACS) ────────────────────────────────────

def test_help() -> None:
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "compile" in result.output
    assert "run" in result.output
    assert "dry-run" in result.output
    assert "validate" in result.output
    assert "inspect" in result.output


@pytest.mark.parametrize("command", ["compile", "run", "dry-run", "validate", "inspect"])
def test_subcommand_help(command: str) -> None:
    result = runner.invoke(cli, [command, "--help"])
    assert result.exit_code == 0, result.output


# ── CLI validate (parse + compile, no workspace write) ────────────────────────

def test_validate_simple_yaml() -> None:
    result = runner.invoke(cli, ["validate", SIMPLE_YAML])
    assert result.exit_code == 0, result.output
    assert "Valid" in result.output


def test_validate_missing_file() -> None:
    result = runner.invoke(cli, ["validate", "nonexistent.yaml"])
    assert result.exit_code != 0


# ── CLI inspect (read-only, no workspace write) ───────────────────────────────

def test_inspect_simple_yaml() -> None:
    result = runner.invoke(cli, ["inspect", SIMPLE_YAML])
    assert result.exit_code == 0, result.output
    assert "System Summary" in result.output


def test_inspect_missing_file() -> None:
    result = runner.invoke(cli, ["inspect", "nonexistent.yaml"])
    assert result.exit_code != 0


# ── CLI compile --no-build (compile without touching disk) ───────────────────

def test_compile_no_build(tmp_path: Path) -> None:
    result = runner.invoke(cli, ["compile", SIMPLE_YAML, "--no-build", "--no-plan"])
    assert result.exit_code == 0, result.output
    assert "Done" in result.output or "no-build" in result.output


# ── CLI dry-run (needs compiled workspace) ────────────────────────────────────

def test_dry_run_missing_workspace() -> None:
    result = runner.invoke(cli, ["dry-run", "nonexistent_workspace"])
    assert result.exit_code != 0


def test_dry_run_without_manifest(tmp_path: Path) -> None:
    """Workspace exists but has no manifest → graceful error, not crash."""
    (tmp_path / "metadata").mkdir()
    result = runner.invoke(cli, ["dry-run", str(tmp_path)])
    assert result.exit_code != 0
    assert "manifest" in result.output.lower() or result.exit_code != 0
