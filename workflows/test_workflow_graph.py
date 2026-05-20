from rich import print
from rich.console import Console

from core.parser import parse_yaml
from core.decision_engine import (
    build_simulation_plan,
)

from workflows.workflow_graph import (
    WorkflowGraph,
)

console = Console()


# ═══════════════════════════════════════════════════════════════════════════════
# Build pipeline
# ═══════════════════════════════════════════════════════════════════════════════

state = parse_yaml(
    "configs/hmg_competition.yaml"
)

plan = build_simulation_plan(state)

graph = WorkflowGraph(plan)


# ═══════════════════════════════════════════════════════════════════════════════
# Validate graph
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold green]Workflow Validation[/bold green]")

try:

    graph.validate()

    print(
        "\n  [green]✓ Workflow DAG válido[/green]"
    )

except Exception as e:

    print(
        f"\n  [red]✖ Validation failed:[/red] {e}"
    )

    raise


# ═══════════════════════════════════════════════════════════════════════════════
# Execution order
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold cyan]Topological Execution Order[/bold cyan]")

ordered = graph.to_execution_view()

for i, step in enumerate(ordered, start=1):

    print(
        f"\n  {i:02d}. "
        f"[bold]{step.title}[/bold]"
    )

    print(
        f"      stage      : {step.stage.value}"
    )

    print(
        f"      step_type  : {step.step_type.value}"
    )

    if step.depends_on:

        print(
            f"      depends_on : {step.depends_on}"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# User View
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold yellow]User Workflow View[/bold yellow]")

user_steps = graph.to_user_view()

for i, step in enumerate(user_steps, start=1):

    print(
        f"\n  {i:02d}. {step}"
    )


# ═══════════════════════════════════════════════════════════════════════════════
# Mermaid export
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold magenta]Mermaid Export[/bold magenta]")

print()

print(
    graph.render_mermaid()
)

print()

console.rule(
    "[bold green]Workflow Graph Test Complete[/bold green]"
)