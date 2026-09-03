# conftest.py — shared fixtures for the full SimForge test suite
from __future__ import annotations

# ── Deterministic CLI help/output rendering across environments ─────────────
# MUST run before anything imports typer / rich. On CI, typer.rich_utils sets
# FORCE_TERMINAL=True whenever GITHUB_ACTIONS / FORCE_COLOR / PY_COLORS is set,
# and Rich treats the runner as a terminal, so it renders styled spans and
# wraps to width — `--legacy` becomes `-` + `-legacy`, numbers get bolded —
# and plain-substring assertions on help / CLI output pass locally but fail in
# CI. Neutralise both: typer's own `_TYPER_FORCE_DISABLE_TERMINAL` opt-out and
# Rich's `TTY_COMPATIBLE=0` ("not a terminal"). These also reach the
# subprocess-based CLI tests via the inherited environment.
import os

os.environ["_TYPER_FORCE_DISABLE_TERMINAL"] = "1"
os.environ["TTY_COMPATIBLE"] = "0"
os.environ["NO_COLOR"] = "1"
os.environ.pop("FORCE_COLOR", None)
os.environ.setdefault("COLUMNS", "200")

from pathlib import Path  # noqa: E402

import pytest  # noqa: E402

from core.parser import parse_yaml  # noqa: E402
from core.decision_engine import build_simulation_plan  # noqa: E402
from core.compiler import SimulationCompiler  # noqa: E402
from builders.workspace_builder import WorkspaceBuilder  # noqa: E402
from executors.shell_executor import ShellExecutor  # noqa: E402


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
