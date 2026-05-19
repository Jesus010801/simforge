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

# ─── Protein Validator ────────────────────────────────────────────────────────
from validators.protein_validator import validate_protein
from pathlib import Path

print("\n[bold green]── Protein Validator ──[/bold green]\n")

proteins = state.get_components_by_role("protein")

if proteins:
    for protein in proteins:
        pdb_path = Path(protein.file)
        if pdb_path.exists():
            pv = validate_protein(pdb_path)
            print(f"  Componente   : {protein.id}")
            print(f"  Archivo      : {pv.source_file}")
            print(f"  Cadenas      : {pv.chains}")
            print(f"  Residuos     : {pv.total_residues}")
            print(f"  Hidrógenos   : {pv.has_hydrogens}")
            print(f"  HETATM       : {pv.has_hetatm}")
            print(f"  Faltantes    : {len(pv.missing_residues)}")
            print(f"  Terminales   : {pv.exposed_termini}")
            print(f"  Oligómero    : {pv.likely_oligomer} → {pv.oligomeric_state}")
            print(f"  Tamaño cadenas: {pv.chain_sizes}")

            print(f"\n[yellow]Warnings ({len(pv.warnings)})[/yellow]")
            for w in pv.warnings:
                print(f"  ⚠  [{w.severity.value}] {w.message}")

            print(f"\n[red]Risks ({len(pv.risks)})[/red]")
            for r in pv.risks:
                print(f"  ✖  [{r.severity.value}] {r.message}")

            print(f"\n[cyan]Recommendations ({len(pv.recommendations)})[/cyan]")
            for rec in pv.recommendations:
                print(f"  →  {rec.message}")
                if rec.action:
                    print(f"     [dim]{rec.action}[/dim]")
        else:
            print(f"  [yellow]Archivo no encontrado: {pdb_path}[/yellow]")
else:
    print("  [yellow]No hay componentes con role: protein[/yellow]")