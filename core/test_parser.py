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

# ─── Ligand Validator ─────────────────────────────────────────────────────────

from pathlib import Path
from validators.ligand_validator import validate_ligand

BASE_DIR = Path(__file__).resolve().parent.parent

print("\n[bold green]── Ligand Validator ──[/bold green]\n")

ligand_roles = [
    "substrate",
    "competitive_ligand",
]

for role in ligand_roles:

    components = state.get_components_by_role(role)

    for comp in components:

        # ─── Resolver ruta absoluta ──────────────────────────────────────────
        lig_path = BASE_DIR / comp.file

        if lig_path.exists():

            lv = validate_ligand(
                lig_path,
                role = comp.role,
            )

            print(
                f"  Componente   : "
                f"{comp.id} [{lv.role}]"
            )

            print(
                f"  Archivo      : "
                f"{lv.source_file}"
            )

            print(
                f"  Parser       : "
                f"{lv.parser_used}"
            )

            print(
                f"  Completo     : "
                f"{lv.is_complete}"
            )

            print(
                f"  Átomos       : "
                f"{lv.n_atoms}"
                f"  |  Bonds: {lv.n_bonds}"
            )

            print(
                f"  Elementos    : "
                f"{sorted(lv.atom_elements)}"
            )

            print(
                f"  Carga neta   : "
                f"{lv.net_charge:+d}"
            )

            print(
                f"  Aromático    : "
                f"{lv.has_aromatic} "
                f"({lv.aromatic_atoms} átomos)"
            )

            print(
                f"  Flexibilidad : "
                f"{lv.estimated_flexibility} "
                f"({lv.n_rotatable_bonds} rot. bonds)"
            )

            print(
                f"  Polaridad    : "
                f"{lv.estimated_polarity}"
            )

            print(
                f"  Param. dif.  : "
                f"{lv.parametrization_difficulty}"
            )
            

            # ─── Geometry / Polarity descriptors ───────────────────────────

            if hasattr(lv, "geometry") and lv.geometry:

                geom = lv.geometry

                print(
                    f"  Shape          : "
                    f"{geom.shape.shape_class.value}"
                )

                print(
                    f"  Planario       : "
                    f"{geom.planarity.is_planar} "
                    f"({geom.planarity.rmsd_from_plane:.3f}Å)"
                )

                print(
                    f"  Rg             : "
                    f"{geom.shape.radius_of_gyration:.2f}Å"
                )

            if hasattr(lv, "polarity") and lv.polarity:

                pol = lv.polarity

                print(
                    f"  logP           : "
                    f"{pol.logp_estimate:.1f} "
                    f"({pol.logp_class})"
                )

                print(
                    f"  HBD/HBA        : "
                    f"{pol.hbond_donors} / "
                    f"{pol.hbond_acceptors}"
                )

                print(
                    f"  Lipinski       : "
                    f"{pol.lipinski_compliant}"
                )

                if pol.functional_groups:

                    fg_str = ", ".join(
                        sorted({
                            fg.name
                            for fg in pol.functional_groups
                        })
                    )

                    print(
                        f"  Grupos         : "
                        f"{fg_str}"
                    )

                print(
                    f"  Solubilidad    : "
                    f"{pol.solubility_class}"
                )

            # ─── Warnings ────────────────────────────────────────────────────

            print(
                f"\n  [yellow]"
                f"Warnings ({len(lv.warnings)})"
                f"[/yellow]"
            )

            for w in lv.warnings:

                print(
                    f"    ⚠  "
                    f"[{w.severity.value}] "
                    f"{w.message}"
                    + (
                        f" → {w.target}"
                        if w.target
                        else ""
                    )
                )

            # ─── Risks ───────────────────────────────────────────────────────

            print(
                f"  [red]"
                f"Risks ({len(lv.risks)})"
                f"[/red]"
            )

            for r in lv.risks:

                print(
                    f"    ✖  "
                    f"[{r.severity.value}] "
                    f"{r.message}"
                    + (
                        f" → {r.target}"
                        if r.target
                        else ""
                    )
                )

            # ─── Recommendations ─────────────────────────────────────────────

            print(
                f"  [cyan]"
                f"Recommendations "
                f"({len(lv.recommendations)})"
                f"[/cyan]"
            )

            for rec in lv.recommendations:

                print(
                    f"    →  "
                    f"{rec.message}"
                )

                if rec.action:

                    print(
                        f"       [dim]"
                        f"{rec.action}"
                        f"[/dim]"
                    )

            print()

        else:

            print(
                f"  [yellow]"
                f"Archivo no encontrado: "
                f"{lig_path}"
                f"[/yellow]\n"
            )