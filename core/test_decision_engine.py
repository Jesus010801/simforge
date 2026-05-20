# core/test_decision_engine.py
"""
Test de integración completo del decision engine.

Pipeline:

    YAML
      ↓
    SystemState enriquecido
      ↓
    Decision Engine
      ↓
    SimulationPlan

Principio:
    Este archivo NO ejecuta lógica manual.
    Todo debe leerse desde:

        - state.*
        - plan.*

Objetivo:
    Verificar que el planner transforma correctamente
    reasoning → acciones ejecutables.
"""

from rich import print
from rich.console import Console
from rich.rule import Rule

from core.parser import parse_yaml
from core.decision_engine import build_simulation_plan

from core.execution_models import (
    PlanStatus,
)

console = Console()


# ═══════════════════════════════════════════════════════════════════════════════
# Pipeline completo
# ═══════════════════════════════════════════════════════════════════════════════

state = parse_yaml("configs/hmg_competition.yaml")

plan = build_simulation_plan(state)


# ═══════════════════════════════════════════════════════════════════════════════
# Encabezado
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold green]SimForge — Decision Engine[/bold green]")

print(f"\n  [bold]Proyecto[/bold]        : {state.project.name}")
print(f"  [bold]Tipo inferido[/bold]   : {state.inferred_system_type}")
print(f"  [bold]Componentes[/bold]     : {state.component_ids()}")
print(f"  [bold]Objetivos[/bold]       : {state.simulation_objectives}")


# ═══════════════════════════════════════════════════════════════════════════════
# Estado global del plan
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold cyan]SimulationPlan[/bold cyan]")

status_colors = {
    PlanStatus.READY:         "green",
    PlanStatus.NEEDS_REVIEW: "yellow",
    PlanStatus.BLOCKED:      "red",
    PlanStatus.EXPERIMENTAL: "magenta",
}

status_color = status_colors.get(plan.status, "white")

print(f"\n  [bold]Estado[/bold]                : [{status_color}]{plan.status.value}[/{status_color}]")
print(f"  [bold]Tipo sistema[/bold]          : {plan.inferred_system_type}")
print(f"  [bold]Steps[/bold]                 : {len(plan.steps)}")
print(f"  [bold]Issues bloqueantes[/bold]   : {len(plan.blocking_issues)}")
print(f"  [bold]Protocols especiales[/bold] : {len(plan.special_protocols)}")
print(f"  [bold]Checklist[/bold]            : {len(plan.checklist)}")


# ═══════════════════════════════════════════════════════════════════════════════
# Notes
# ═══════════════════════════════════════════════════════════════════════════════

if plan.notes:

    print("\n  [bold]── Notes[/bold]")

    for note in plan.notes:
        print(f"    ℹ  {note}")


# ═══════════════════════════════════════════════════════════════════════════════
# Blocking Issues
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold red]Blocking Issues[/bold red]")

if not plan.blocking_issues:

    print("\n  [green]✓ No se detectaron issues bloqueantes[/green]")

else:

    for issue in plan.blocking_issues:

        sev_color = {
            "low": "cyan",
            "medium": "yellow",
            "high": "red",
        }.get(issue.severity.value, "white")

        print(
            f"\n  [{sev_color}]✖ [{issue.severity.value.upper()}][/{sev_color}] "
            f"{issue.source}"
        )

        print(f"      {issue.message}")

        if issue.resolution:
            print(f"      [dim]{issue.resolution}[/dim]")


# ═══════════════════════════════════════════════════════════════════════════════
# Special Protocols
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold magenta]Special Protocols[/bold magenta]")

if not plan.special_protocols:

    print("\n  [green]✓ No se requieren protocolos especiales[/green]")

else:

    for protocol in plan.special_protocols:
        print(f"\n  →  {protocol}")


# ═══════════════════════════════════════════════════════════════════════════════
# Workflow completo
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold yellow]Simulation Workflow[/bold yellow]")

for i, step in enumerate(plan.steps, start=1):

    stage_color = {
        
    "preparation":       "cyan",

    "validation":        "bright_yellow",

    "parametrization":   "yellow",

    "assembly":          "magenta",

    "minimization":      "green",

    "equilibration":     "blue",

    "enhanced_sampling": "bright_magenta",

    "production":        "red",

    "analysis":          "white",

    }.get(step.stage.value, "white")

    print(
        f"\n  [{stage_color}]{i:02d}. "
        f"{step.title}[/{stage_color}]"
    )

    print(f"      step_id     : {step.step_id}")
    print(f"      stage       : {step.stage.value}")
    print(f"      engine      : {step.engine}")

    if step.target_components:
        print(f"      targets     : {step.target_components}")

    print(f"      required    : {step.required}")
    print(f"      blocking    : {step.blocking}")

    if step.depends_on:
        print(f"      depends_on  : {step.depends_on}")

    if step.notes:
        print(f"      notes:")
        for note in step.notes:
            print(f"        - {note}")


# ═══════════════════════════════════════════════════════════════════════════════
# Checklist previo a ejecución
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold blue]Execution Checklist[/bold blue]")

if not plan.checklist:

    print("\n  [green]✓ Checklist vacío[/green]")

else:

    for item in plan.checklist:

        icon = "✓" if not item.required else "•"

        print(f"\n  {icon}  {item.message}")

        if item.component:
            print(f"      componente: {item.component}")


# ═══════════════════════════════════════════════════════════════════════════════
# Estadísticas del workflow
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold white]Workflow Statistics[/bold white]")

stage_counts = {}

for step in plan.steps:

    stage = step.stage.value

    if stage not in stage_counts:
        stage_counts[stage] = 0

    stage_counts[stage] += 1

print()

for stage, count in stage_counts.items():
    print(f"  {stage:<18} : {count}")

blocking_steps = sum(
    1 for s in plan.steps if s.blocking
)

required_steps = sum(
    1 for s in plan.steps if s.required
)

print(f"\n  blocking steps    : {blocking_steps}")
print(f"  required steps    : {required_steps}")
print(f"  optional steps    : {len(plan.steps) - required_steps}")


# ═══════════════════════════════════════════════════════════════════════════════
# DAG simplificado
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold green]Dependency Graph[/bold green]")

for step in plan.steps:

    if step.depends_on:

        deps = ", ".join(step.depends_on)

        print(
            f"\n  {step.step_id}"
        )

        print(
            f"      ← depends on: {deps}"
        )

    else:

        print(
            f"\n  {step.step_id}"
        )

        print(
            f"      ← root step"
        )


console.rule("[bold green]Decision Engine Test Complete[/bold green]")