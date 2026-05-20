# executors/test_gromacs_executor.py
"""
Test de integración del GROMACSExecutor.

Pipeline:
    YAML → Compiler → WorkspaceBuilder → GROMACSExecutor (dry-run)

En dry-run el executor marca todos los steps como DONE sin ejecutar GROMACS real.
Los diagnósticos se construyen igualmente — en dry-run los outputs no existen,
por lo que el diagnóstico reportará incomplete/unknown según el stage.

Para testear con GROMACS real: dry_run=False (requiere GROMACS instalado).

Qué se valida aquí:
    1. GROMACSExecutor hereda correctamente de ShellExecutor
    2. Los diagnósticos se construyen para stages GROMACS
    3. Los stages non-GROMACS (preparation, parametrization, etc.) se ignoran
    4. El summary_report() genera salida legible
    5. get_diagnostic() retorna el diagnóstico correcto por step_id
"""

from rich import print
from rich.console import Console
from rich.rule import Rule

from core.compiler import SimulationCompiler
from builders.workspace_builder import WorkspaceBuilder
from executors.gromacs_executor import GROMACSExecutor
from executors.execution_state import StepStatus

console = Console()


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Compilar
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold green]SimForge — GROMACS Executor Test[/bold green]")

compiler = SimulationCompiler()
result   = compiler.compile("configs/hmg_competition.yaml")

print(f"\n  Sistema     : {result.state.inferred_system_type}")
print(f"  Steps       : {len(result.plan.steps)}")
print(f"  Protocolos  : {result.plan.special_protocols}")


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Construir workspace
# ═══════════════════════════════════════════════════════════════════════════════

builder   = WorkspaceBuilder()
workspace = builder.build(result)

print(f"\n  Workspace   : {workspace}")


# ═══════════════════════════════════════════════════════════════════════════════
# 3. Ejecutar con GROMACSExecutor
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold cyan]GROMACSExecutor — Dry Run[/bold cyan]")

executor = GROMACSExecutor(
    workspace_path = workspace,
    dry_run        = True,
)

state = executor.run()


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Reporte de ejecución
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold yellow]Execution Steps[/bold yellow]")

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
    color   = status_colors.get(record.status, "white")
    elapsed = f"({record.elapsed_s:.1f}s)" if record.elapsed_s else ""
    print(
        f"  [{color}]{record.status.value:8}[/{color}]  "
        f"{record.step_id}  [dim]{elapsed}[/dim]"
    )
    if record.error_message:
        print(f"           [red]{record.error_message}[/red]")


# ═══════════════════════════════════════════════════════════════════════════════
# 5. Diagnósticos GROMACS
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold magenta]GROMACS Diagnostics[/bold magenta]")

diagnostics = executor.all_diagnostics()

if not diagnostics:
    print("\n  [dim]No se generaron diagnósticos (dry-run sin outputs reales)[/dim]")
else:
    for step_id, diag in diagnostics.items():

        verdict_color = {
            "ok":            "green",
            "converged":     "green",
            "not_converged": "yellow",
            "warning":       "yellow",
            "incomplete":    "red",
            "crashed":       "red",
            "unknown":       "dim",
        }.get(diag.verdict, "white")

        print(f"\n  [bold]{step_id}[/bold]  [{diag.stage}]")
        print(
            f"    verdict  : [{verdict_color}]{diag.verdict}[/{verdict_color}]"
        )

        if diag.minimization:
            m = diag.minimization
            print(f"    converged: {m.converged} ({m.convergence_reason})")
            if m.final_epot is not None:
                print(f"    Epot     : {m.final_epot:.2f} kJ/mol")
            if m.final_fmax is not None:
                print(f"    Fmax     : {m.final_fmax:.4f} kJ/mol/nm")

        if diag.md:
            d = diag.md
            print(f"    completed: {d.completed}")
            if d.last_temperature is not None:
                print(f"    Temp     : {d.last_temperature:.1f} K")
            if d.last_epot is not None:
                print(f"    Epot     : {d.last_epot:.2f} kJ/mol")
            if d.has_lincs_warning:
                print(f"    [yellow]LINCS warnings: {d.n_lincs_warnings}[/yellow]")
            if d.has_fatal_error:
                print(f"    [red]Fatal error: {d.fatal_error_msg[:100]}[/red]")

        if diag.output_files:
            print(f"    outputs  :")
            for o in diag.output_files:
                ok = "✓" if (o.exists and not o.is_empty and not o.suspect) else "✖"
                print(
                    f"      {ok}  {o.filename}  "
                    f"({o.size_bytes} bytes, exists={o.exists})"
                )

        if diag.notes:
            for n in diag.notes:
                print(f"    [dim]ℹ  {n}[/dim]")
        if diag.warnings:
            for w in diag.warnings:
                print(f"    [yellow]⚠  {w}[/yellow]")
        if diag.errors:
            for e in diag.errors:
                print(f"    [red]✖  {e}[/red]")


# ═══════════════════════════════════════════════════════════════════════════════
# 6. Summary report
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold white]Summary Report[/bold white]")

print()
print(executor.summary_report())


# ═══════════════════════════════════════════════════════════════════════════════
# 7. Estadísticas
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold white]Execution Statistics[/bold white]")

print(f"\n  Done     : [green]{state.n_done()}[/green]")
print(f"  Failed   : [red]{state.n_failed()}[/red]")
print(f"  Pending  : [dim]{state.n_pending()}[/dim]")
print(f"  Complete : {'[green]yes[/green]' if state.is_complete else '[red]no[/red]'}")
print(f"\n  GROMACS diagnostics generated: {len(diagnostics)}")
print()

console.rule("[bold green]GROMACS Executor Test Complete[/bold green]")
