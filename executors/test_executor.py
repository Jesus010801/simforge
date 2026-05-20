# executors/test_executor.py
from rich import print
from rich.console import Console

from core.compiler import SimulationCompiler
from builders.workspace_builder import WorkspaceBuilder
from executors.shell_executor import ShellExecutor
from executors.execution_state import StepStatus

console = Console()


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Compilar
# ═══════════════════════════════════════════════════════════════════════════════

compiler = SimulationCompiler()
result   = compiler.compile("configs/hmg_competition.yaml")


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Construir workspace
# ═══════════════════════════════════════════════════════════════════════════════

builder   = WorkspaceBuilder()
workspace = builder.build(result)


# ═══════════════════════════════════════════════════════════════════════════════
# 3. Ejecutar (dry-run)
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold green]SimForge Executor — Dry Run[/bold green]")

executor = ShellExecutor(
    workspace_path = workspace,
    dry_run        = True,
)

state = executor.run()


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Reporte final
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold cyan]Execution Report[/bold cyan]")

status_colors = {
    StepStatus.DONE:    "green",
    StepStatus.FAILED:  "red",
    StepStatus.SKIPPED: "yellow",
    StepStatus.BLOCKED: "red",
    StepStatus.PENDING: "dim",
    StepStatus.RUNNING: "blue",
}

print()

for record in state.steps:

    color = status_colors.get(record.status, "white")
    elapsed = f"({record.elapsed_s:.1f}s)" if record.elapsed_s else ""

    print(
        f"  [{color}]{record.status.value:8}[/{color}]  "
        f"{record.step_id}  [dim]{elapsed}[/dim]"
    )

    if record.error_message:
        print(f"           [red]{record.error_message}[/red]")

    if record.outputs_missing:
        print(f"           [yellow]missing: {record.outputs_missing}[/yellow]")

print()

console.rule("[bold white]Summary[/bold white]")

print(f"\n  Done    : [green]{state.n_done()}[/green]")
print(f"  Failed  : [red]{state.n_failed()}[/red]")
print(f"  Pending : [dim]{state.n_pending()}[/dim]")
print(f"  Complete: {'[green]yes[/green]' if state.is_complete else '[red]no[/red]'}")
print(f"\n  State saved → {workspace}/execution_state.json")
print()

console.rule("[bold green]Executor Test Complete[/bold green]")
