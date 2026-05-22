# conftest.py — shared fixtures for the full SimForge test suite
from __future__ import annotations

from pathlib import Path

import pytest

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
