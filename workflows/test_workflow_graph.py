# workflows/test_workflow_graph.py
"""
Test de integración del WorkflowGraph.

Pipeline completo:
    YAML → parse_yaml() → build_simulation_plan() → WorkflowGraph

Verifica:
    - DAG válido (sin ciclos, sin dependencias rotas)
    - Orden topológico correcto
    - Ondas de ejecución paralela
    - Camino crítico
    - Vista por stage
    - Vista de usuario
    - Estadísticas del grafo
    - Export Mermaid agrupado por stage
"""

from rich import print
from rich.console import Console
from rich.table import Table

from core.parser import parse_yaml
from core.decision_engine import build_simulation_plan
from workflows.workflow_graph import WorkflowGraph
from core.execution_models import StepType

console = Console()


# ═══════════════════════════════════════════════════════════════════════════════
# Pipeline
# ═══════════════════════════════════════════════════════════════════════════════

state = parse_yaml("configs/hmg_competition.yaml")
plan  = build_simulation_plan(state)
graph = WorkflowGraph(plan)


# ═══════════════════════════════════════════════════════════════════════════════
# Validación del DAG
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold green]Workflow Validation[/bold green]")

try:
    graph.validate()
    print(f"\n  [green]✓ DAG válido[/green] — {graph}")

except Exception as e:
    print(f"\n  [red]✖ Validation failed:[/red] {e}")
    raise


# ═══════════════════════════════════════════════════════════════════════════════
# Estadísticas del grafo
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold white]Graph Statistics[/bold white]")

s = graph.stats()

print(f"\n  Steps totales     : {s['n_steps']}")
print(f"  Edges (deps)      : {s['n_edges']}")
print(f"  Root steps        : {s['n_root_steps']}")
print(f"  Leaf steps        : {s['n_leaf_steps']}")
print(f"  Blocking steps    : {s['n_blocking_steps']}")
print(f"  Parallel waves    : {s['n_parallel_waves']}")
print(f"  Critical path len : {s['critical_path_len']}")
print(f"  Plan status       : {s['plan_status']}")
print(f"  Blocking issues   : {s['n_blocking_issues']}")

print(f"\n  [bold]Steps por stage:[/bold]")
for stage, count in s["stage_counts"].items():
    print(f"    {stage:<20} : {count}")

print(f"\n  [bold]Steps por tipo:[/bold]")
for stype, count in s["type_counts"].items():
    print(f"    {stype:<20} : {count}")


# ═══════════════════════════════════════════════════════════════════════════════
# Nodos raíz y hoja
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold cyan]Roots & Leaves[/bold cyan]")

roots = graph.root_steps()
leaves = graph.leaf_steps()

print(f"\n  [bold]Root steps[/bold] (sin dependencias):")
for r in roots:
    print(f"    → {r.step_id}  [{r.stage.value}]")

print(f"\n  [bold]Leaf steps[/bold] (sin sucesores):")
for l in leaves:
    print(f"    → {l.step_id}  [{l.stage.value}]")


# ═══════════════════════════════════════════════════════════════════════════════
# Camino crítico
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold red]Critical Path[/bold red]")

cp = graph.critical_path()

print(f"\n  Longitud: {len(cp)} steps\n")

for i, step in enumerate(cp):
    connector = "└─" if i == len(cp) - 1 else "├─"
    blocking_tag = " [red][BLOCKING][/red]" if step.blocking else ""
    print(f"  {connector} {i+1:02d}. {step.step_id}{blocking_tag}")
    print(f"       stage: {step.stage.value}")


# ═══════════════════════════════════════════════════════════════════════════════
# Ondas de ejecución paralela
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold magenta]Parallel Execution Waves[/bold magenta]")

waves = graph.to_parallel_waves()

for i, wave in enumerate(waves, start=1):
    parallel_tag = " [dim](parallelizable)[/dim]" if len(wave) > 1 else ""
    print(f"\n  Wave {i}{parallel_tag}")
    for step in wave:
        icon = "⚙" if step.step_type.value == "automatic" else "✋"
        blocking_tag = " [red]BLOCKING[/red]" if step.blocking else ""
        print(f"    {icon}  {step.step_id}  [{step.stage.value}]{blocking_tag}")


# ═══════════════════════════════════════════════════════════════════════════════
# Vista por stage
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold yellow]Stage View[/bold yellow]")

stage_view = graph.to_stage_view()

for stage, steps in stage_view.items():
    print(f"\n  [bold]{stage.value.upper()}[/bold]  ({len(steps)} steps)")
    for step in steps:
        deps = f"  ← {step.depends_on}" if step.depends_on else ""
        print(f"    • {step.step_id}{deps}")


# ═══════════════════════════════════════════════════════════════════════════════
# Orden de ejecución completo
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold cyan]Topological Execution Order[/bold cyan]")

ordered = graph.to_execution_view()

for i, step in enumerate(ordered, start=1):
    type_color = {
        "automatic":  "green",
        "manual":     "yellow",
        "external":   "magenta",
        "validation": "cyan",
    }.get(step.step_type.value, "white")

    blocking_tag = "  [red][BLOCKING][/red]" if step.blocking else ""

    print(f"\n  {i:02d}. [bold]{step.title}[/bold]{blocking_tag}")
    print(f"       step_id   : {step.step_id}")
    print(f"       stage     : {step.stage.value}")
    print(f"       type      : [{type_color}]{step.step_type.value}[/{type_color}]")
    print(f"       engine    : {step.engine}")

    if step.depends_on:
        print(f"       depends_on: {step.depends_on}")

    if step.target_components:
        print(f"       targets   : {step.target_components}")

    if step.notes:
        for note in step.notes:
            print(f"       [dim]ℹ  {note}[/dim]")


# ═══════════════════════════════════════════════════════════════════════════════
# Helpers de consulta — ejemplo con un step del medio
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold white]Graph Queries (ejemplo)[/bold white]")

# Usar el primer step de assembly como ejemplo de consulta
example_id = "assemble_system"
if example_id in graph:

    print(f"\n  Consulta sobre: [bold]{example_id}[/bold]")

    preds = graph.predecessors(example_id)
    succs = graph.successors(example_id)
    ancs  = graph.ancestors(example_id)
    descs = graph.descendants(example_id)

    print(f"\n  Predecesores directos ({len(preds)}):")
    for p in preds:
        print(f"    ← {p.step_id}")

    print(f"\n  Sucesores directos ({len(succs)}):")
    for s in succs:
        print(f"    → {s.step_id}")

    print(f"\n  Ancestros transitivos ({len(ancs)}):")
    for a in ancs:
        print(f"    ↑ {a.step_id}")

    print(f"\n  Descendientes transitivos ({len(descs)}):")
    for d in descs:
        print(f"    ↓ {d.step_id}")

else:
    print(f"\n  [dim]Step '{example_id}' no encontrado en este plan[/dim]")


# ═══════════════════════════════════════════════════════════════════════════════
# Vista de usuario
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold yellow]User Workflow View[/bold yellow]")

user_steps = graph.to_user_view()

for i, step_text in enumerate(user_steps, start=1):
    print(f"\n  {i:02d}. {step_text}")


# ═══════════════════════════════════════════════════════════════════════════════
# Mermaid export
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold magenta]Mermaid (agrupado por stage)[/bold magenta]")

print()
print(graph.render_mermaid(group_by_stage=True))
print()

console.rule("[bold magenta]Mermaid (plano)[/bold magenta]")

print()
print(graph.render_mermaid(group_by_stage=False))
print()

console.rule("[bold green]Workflow Graph Test Complete[/bold green]")