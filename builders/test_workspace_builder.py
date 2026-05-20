from rich import print
from rich.console import Console

from core.compiler import (
    SimulationCompiler,
)

from builders.workspace_builder import (
    WorkspaceBuilder,
)

console = Console()


# ═══════════════════════════════════════════════════════════════════════════════
# Compile workflow
# ═══════════════════════════════════════════════════════════════════════════════

compiler = SimulationCompiler()

result = compiler.compile(
    "configs/hmg_competition.yaml"
)


# ═══════════════════════════════════════════════════════════════════════════════
# Build workspace
# ═══════════════════════════════════════════════════════════════════════════════

builder = WorkspaceBuilder()

workspace = builder.build(result)

console.rule("[bold green]Workspace Created[/bold green]")

print()

print(
    f"  Workspace path: [bold]{workspace}[/bold]"
)

print()

console.rule(
    "[bold cyan]Workspace Build Complete[/bold cyan]"
)