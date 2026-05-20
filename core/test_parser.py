# core/test_parser.py
"""
Test de integración del pipeline completo de SimForge.

Lee el SystemState enriquecido y muestra todos los campos
por componente: config → validation → descriptors → reasoning.

Principio: este archivo NUNCA llama directamente a validators
ni a descriptores. Todo se lee desde state.components[i].* y
state.global_reasoning.
"""

from core.parser import parse_yaml
from rich import print
from rich.console import Console
from rich.rule import Rule

console = Console()

# ═══════════════════════════════════════════════════════════════════════════════
# Carga del pipeline completo
# ═══════════════════════════════════════════════════════════════════════════════

state = parse_yaml("configs/hmg_competition.yaml")


# ═══════════════════════════════════════════════════════════════════════════════
# Resumen global del sistema
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold green]SimForge — SystemState[/bold green]")

print(f"\n  [bold]Proyecto[/bold]      : {state.project.name}")
print(f"  [bold]Descripción[/bold]   : {state.project.description}")
print(f"  [bold]Tipo inferido[/bold] : {state.inferred_system_type}")
print(f"  [bold]Componentes[/bold]   : {state.component_ids()}")
print(f"  [bold]Objetivos[/bold]     : {state.simulation_objectives}")
print(f"  [bold]Membrana[/bold]      : {state.has_membrane()}")
print(f"  [bold]FF proteína[/bold]   : {state.forcefields.protein}")
print(f"  [bold]FF ligandos[/bold]   : {state.forcefields.ligands}")


# ═══════════════════════════════════════════════════════════════════════════════
# Global reasoning
# ═══════════════════════════════════════════════════════════════════════════════

gr = state.global_reasoning
console.rule("[bold cyan]Global Reasoning[/bold cyan]")

ready_color = "green" if gr.system_is_ready else "red"
print(f"\n  Sistema listo     : [{ready_color}]{gr.system_is_ready}[/{ready_color}]")
print(f"  Errores bloquean  : {gr.has_blocking_errors}")
print(f"  Sampling especial : {gr.needs_special_sampling}")
print(f"  Validados / Total : {gr.n_components_validated} / {len(state.components)}")
print(f"  Con errores       : {gr.n_components_with_errors}")

for note in gr.notes:
    print(f"\n  [dim]ℹ  {note}[/dim]")

# Warnings del inference pipeline (compatibilidad hacia atrás)
if state.warnings:
    print(f"\n  [yellow]Warnings del sistema ({len(state.warnings)})[/yellow]")
    for w in state.warnings:
        print(f"    ⚠  [{w.severity.value}] {w.message}"
              + (f"  → {w.target}" if w.target else ""))

if state.risks:
    print(f"\n  [red]Risks del sistema ({len(state.risks)})[/red]")
    for r in state.risks:
        print(f"    ✖  [{r.severity.value}] {r.message}"
              + (f"  → {r.target}" if r.target else ""))

if state.recommendations:
    print(f"\n  [cyan]Recommendations del sistema ({len(state.recommendations)})[/cyan]")
    for rec in state.recommendations:
        print(f"    →  {rec.message}")
        if rec.action:
            print(f"       [dim]{rec.action}[/dim]")


# ═══════════════════════════════════════════════════════════════════════════════
# Por componente: config → validation → descriptors → reasoning
# ═══════════════════════════════════════════════════════════════════════════════

for comp in state.components:
    console.rule(f"[bold white]Componente: {comp.id}  [{comp.role}][/bold white]")

    # ── Config ────────────────────────────────────────────────────────────────
    print(f"\n  [bold]── Config[/bold]")
    print(f"    id               : {comp.id}")
    print(f"    role             : {comp.role}")
    print(f"    file             : {comp.file}")
    print(f"    bio_context      : {comp.biological_context}")
    print(f"    validado         : {comp.is_validated}")
    print(f"    tiene descriptors: {comp.has_descriptors}")

    # ── Validation ────────────────────────────────────────────────────────────
    val = comp.validation
    if val:
        print(f"\n  [bold]── Validation[/bold]  [{val.validator_used}]")
        v_color = "green" if val.is_valid else "red"
        print(f"    is_valid         : [{v_color}]{val.is_valid}[/{v_color}]")
        if val.validation_error:
            print(f"    error            : [red]{val.validation_error}[/red]")

        # Datos específicos según tipo de validator
        d = val.data
        if val.validator_used == "protein_validator" and d:
            print(f"    cadenas          : {d.get('chains', [])}")
            print(f"    residuos         : {d.get('total_residues', 0)}")
            print(f"    hidrógenos       : {d.get('has_hydrogens', False)}")
            print(f"    HETATM           : {d.get('has_hetatm', [])}")
            print(f"    faltantes        : {len(d.get('missing_residues', []))}")
            print(f"    terminales exp.  : {d.get('exposed_termini', [])}")
            print(f"    oligómero        : {d.get('likely_oligomer', False)} → {d.get('oligomeric_state', None)}")
            print(f"    tamaño cadenas   : {d.get('chain_sizes', {})}")

        elif val.validator_used == "ligand_validator" and d:
            print(f"    parser           : {d.get('parser_used', '?')}")
            print(f"    átomos           : {d.get('n_atoms', 0)}")
            print(f"    bonds            : {d.get('n_bonds', 0)}")
            print(f"    elementos        : {sorted(d.get('atom_elements', []))}")
            print(f"    carga neta       : {d.get('net_charge', 0):+d}")
            print(f"    aromático        : {d.get('has_aromatic', False)} ({d.get('aromatic_atoms', 0)} átomos)")
            print(f"    flexibilidad     : {d.get('estimated_flexibility', '?')} ({d.get('n_rotatable_bonds', 0)} rot.bonds)")
            print(f"    polaridad        : {d.get('estimated_polarity', '?')}")
            print(f"    param. dif.      : {d.get('parametrization_difficulty', '?')}")
            print(f"    Lipinski         : {d.get('lipinski_compliant', True)}")

        if val.warnings:
            print(f"\n    [yellow]Warnings validation ({len(val.warnings)})[/yellow]")
            for w in val.warnings:
                print(f"      ⚠  [{w.severity.value}] {w.message}")
        if val.risks:
            print(f"    [red]Risks validation ({len(val.risks)})[/red]")
            for r in val.risks:
                print(f"      ✖  [{r.severity.value}] {r.message}")
        if val.recommendations:
            print(f"    [cyan]Recommendations validation ({len(val.recommendations)})[/cyan]")
            for rec in val.recommendations:
                print(f"      →  {rec.message}")
                if rec.action:
                    print(f"         [dim]{rec.action}[/dim]")
    else:
        print(f"\n  [dim]── Validation: no ejecutada para este rol[/dim]")

    # ── Descriptors ───────────────────────────────────────────────────────────
    desc = comp.descriptors
    if desc:
        print(f"\n  [bold]── Descriptors[/bold]")
        print(f"    n_heavy_atoms    : {desc.n_heavy_atoms}")
        print(f"    rings            : {desc.n_aromatic_rings} aromáticos, {desc.n_fused_aromatic} fusionados")
        print(f"    flexibilidad     : {desc.flexibility_class}  (eff={desc.effective_rot_bonds:.1f}, score={0.0})")
        print(f"    scaffold         : {desc.scaffold_rigidity}")
        print(f"    shape            : {desc.shape_class}")
        print(f"    planitud global  : {desc.is_globally_planar} (RMS={desc.global_planarity_rms:.3f}Å)")
        print(f"    Rg               : {desc.radius_of_gyration:.2f}Å")
        print(f"    polaridad        : {desc.polarity_class} (score={desc.polarity_score:.2f})")
        print(f"    HBD / HBA        : {desc.hb_donors} / {desc.hb_acceptors}")
        print(f"    carga neta       : {desc.net_charge:+d}")
        print(f"    anfipático       : {desc.amphipathic_class}")
        print(f"    grupos func.     : {desc.n_functional_groups}")
        print(f"    Lipinski         : HBD={desc.lipinski_hbd} HBA={desc.lipinski_hba} → {'✓' if desc.passes_lipinski else '✗'}")
        if desc.sampling_recommendation:
            print(f"    muestreo         : [dim]{desc.sampling_recommendation[:80]}...[/dim]")
    else:
        print(f"\n  [dim]── Descriptors: no calculados para este rol[/dim]")

    # ── Reasoning ─────────────────────────────────────────────────────────────
    reas = comp.reasoning
    if reas:
        print(f"\n  [bold]── Reasoning[/bold]")
        print(f"    sampling especial  : {reas.needs_special_sampling}")
        print(f"    revisar param.     : {reas.needs_parametrization_review}")
        print(f"    verificar prot.    : {reas.needs_protonation_check}")
        print(f"    validar pose       : {reas.needs_pose_validation}")
        print(f"    dificultad param.  : {reas.parametrization_difficulty}")

        if reas.notes:
            for note in reas.notes:
                print(f"    [dim]ℹ  {note}[/dim]")

        if reas.warnings:
            print(f"\n    [yellow]Warnings reasoning ({len(reas.warnings)})[/yellow]")
            for w in reas.warnings:
                print(f"      ⚠  [{w.severity.value}] {w.message}")
        if reas.risks:
            print(f"    [red]Risks reasoning ({len(reas.risks)})[/red]")
            for r in reas.risks:
                print(f"      ✖  [{r.severity.value}] {r.message}")
        if reas.recommendations:
            print(f"    [cyan]Recommendations reasoning ({len(reas.recommendations)})[/cyan]")
            for rec in reas.recommendations:
                print(f"      →  {rec.message}")
                if rec.action:
                    print(f"         [dim]{rec.action}[/dim]")
    else:
        print(f"\n  [dim]── Reasoning: no ejecutado para este rol[/dim]")


# ═══════════════════════════════════════════════════════════════════════════════
# Vista agregada de todos los issues del sistema
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold red]Issues agregados del sistema completo[/bold red]")

all_warnings = state.collect_all_warnings()
all_risks    = state.collect_all_risks()
all_recs     = state.collect_all_recommendations()

print(f"\n  [yellow]Total warnings : {len(all_warnings)}[/yellow]")
for src, w in all_warnings:
    print(f"    [{src}] ⚠  [{w.severity.value}] {w.message}")

print(f"\n  [red]Total risks    : {len(all_risks)}[/red]")
for src, r in all_risks:
    print(f"    [{src}] ✖  [{r.severity.value}] {r.message}")

print(f"\n  [cyan]Total recs     : {len(all_recs)}[/cyan]")
for src, rec in all_recs:
    print(f"    [{src}] →  {rec.message}")

console.rule()
