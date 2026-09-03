# conftest.py — shared fixtures for the full SimForge test suite
from __future__ import annotations

import os
from pathlib import Path

import pytest

# Deterministic CLI help/output rendering across environments. GitHub Actions
# sets FORCE_COLOR, which makes rich/typer split option flags across ANSI style
# spans (e.g. `--legacy` renders as `-` + `-legacy`), so a plain-substring check
# like `"--legacy" in help_text` fails there but passes locally. Force no-color,
# wide output for the whole test session (also propagates to subprocess-based
# CLI tests via the environment).
os.environ["NO_COLOR"] = "1"
os.environ.pop("FORCE_COLOR", None)
os.environ.setdefault("COLUMNS", "200")

from core.parser import parse_yaml
from core.decision_engine import build_simulation_plan
from core.compiler import SimulationCompiler
from builders.workspace_builder import WorkspaceBuilder
from executors.shell_executor import ShellExecutor


YAML_CONFIG = "configs/hmg_competition.yaml"


@pytest.fixture(scope="session")
def state():
    return parse_yaml(YAML_CONFIG)


@pytest.fixture(scope="session")
def plan(state):
    return build_simulation_plan(state)


@pytest.fixture(scope="session")
def compilation_result():
    return SimulationCompiler().compile(YAML_CONFIG)


@pytest.fixture(scope="session")
def workspace(compilation_result, tmp_path_factory):
    tmp = tmp_path_factory.mktemp("workspace")
    return WorkspaceBuilder().build(compilation_result, output_dir=str(tmp))


@pytest.fixture(scope="session")
def execution_state(workspace):
    return ShellExecutor(workspace, dry_run=True).run()
