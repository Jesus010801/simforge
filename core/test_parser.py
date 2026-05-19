# core/test_parser.py

from core.parser import parse_yaml
from rich import print
from rich.console import Console
from rich.table import Table

console = Console()
state = parse_yaml("configs/hmg_competition.yaml")

print("\n[bold green]SystemState cargado correctamente[/bold green]\n")
print(f"  Proyecto     : {state.project.name}")
print(f"  Descripción  : {state.project.description}")
print(f"  Tipo inferido: {state.inferred_system_type}")
print(f"  Componentes  : {state.component_ids()}")
print(f"  Objetivos    : {state.simulation_objectives}")
print(f"  Membrana     : {state.has_membrane()}")
print(f"  FF proteína  : {state.forcefields.protein}")

# ─── Warnings ────────────────────────────────────────────────────────────────
print(f"\n[bold yellow]Warnings ({len(state.warnings)})[/bold yellow]")
for w in state.warnings:
    print(f"  ⚠  [{w.severity.value}] {w.message}"
          + (f" → {w.target}" if w.target else ""))

# ─── Risks ───────────────────────────────────────────────────────────────────
print(f"\n[bold red]Risks ({len(state.risks)})[/bold red]")
for r in state.risks:
    print(f"  ✖  [{r.severity.value}] {r.message}"
          + (f" → {r.target}" if r.target else ""))

# ─── Recommendations ─────────────────────────────────────────────────────────
print(f"\n[bold cyan]Recommendations ({len(state.recommendations)})[/bold cyan]")
for rec in state.recommendations:
    print(f"  →  {rec.message}")
    if rec.action:
        print(f"     [dim]{rec.action}[/dim]")