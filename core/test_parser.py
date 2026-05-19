# core/test_parser.py

from core.parser import parse_yaml
from rich import print

state = parse_yaml("configs/hmg_competition.yaml")

print("\n[bold green]SystemState cargado correctamente[/bold green]\n")
print(f"  Proyecto     : {state.project.name}")
print(f"  Descripción  : {state.project.description}")
print(f"  Tipo inferido: {state.inferred_system_type}")
print(f"  Componentes  : {state.component_ids()}")
print(f"  Objetivos    : {state.simulation_objectives}")
print(f"  Membrana     : {state.has_membrane()}")
print(f"  Bio context  : membrane_associated → {state.has_biological_context('membrane_associated')}")
print(f"  Restraints   : {[r.type for r in state.restraints]}")
print(f"  Análisis     : {[a.type for a in state.analysis]}")
print(f"  FF proteína  : {state.forcefields.protein}")
print(f"  Warnings     : {state.warnings}")
print(f"\n[bold yellow]Warnings ({len(state.warnings)})[/bold yellow]")
for w in state.warnings:
    print(f"  ⚠  {w}")
