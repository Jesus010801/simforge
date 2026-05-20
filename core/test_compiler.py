from rich import print
from rich.console import Console

from core.compiler import (
    SimulationCompiler,
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
# Summary
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold green]SimForge Compiler[/bold green]")

print()

for item in result.summary:

    print(f"  • {item}")


# ═══════════════════════════════════════════════════════════════════════════════
# User workflow
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold cyan]User Workflow[/bold cyan]")

for i, step in enumerate(result.user_view, start=1):

    print(
        f"\n  {i:02d}. {step}"
    )


# ═══════════════════════════════════════════════════════════════════════════════
# Mermaid
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold magenta]Mermaid Graph[/bold magenta]")

print()

print(
    result.mermaid_graph
)

print()

console.rule(
    "[bold green]Compilation Complete[/bold green]"
)