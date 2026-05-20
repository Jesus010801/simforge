# workflows/workflow_graph.py

from __future__ import annotations

from collections import defaultdict, deque

from core.execution_models import (
    SimulationPlan,
    SimulationStep,
)


# ═══════════════════════════════════════════════════════════════════════════════
# Workflow Graph
# ═══════════════════════════════════════════════════════════════════════════════

class WorkflowGraph:
    """
    Representación DAG del workflow de simulación.

    Responsabilidades:
        - construir grafo
        - validar estructura
        - ordenar topológicamente
        - generar vistas simplificadas
    """

    def __init__(self, plan: SimulationPlan):

        self.plan = plan

        self.nodes: dict[str, SimulationStep] = {}

        self.edges: dict[str, list[str]] = defaultdict(list)

        self.reverse_edges: dict[str, list[str]] = defaultdict(list)

        self.sorted_steps: list[SimulationStep] = []

        self._build_graph()


    # ═══════════════════════════════════════════════════════════════════════════
    # Build graph
    # ═══════════════════════════════════════════════════════════════════════════

    def _build_graph(self) -> None:

        # ── registrar nodos ───────────────────────────────────────────────────

        for step in self.plan.steps:

            if step.step_id in self.nodes:

                raise ValueError(
                    f"Duplicated step_id detected: {step.step_id}"
                )

            self.nodes[step.step_id] = step

        # ── registrar dependencias ───────────────────────────────────────────

        for step in self.plan.steps:

            for dep in step.depends_on:

                if dep not in self.nodes:

                    raise ValueError(
                        f"Missing dependency '{dep}' "
                        f"for step '{step.step_id}'"
                    )

                self.edges[dep].append(step.step_id)

                self.reverse_edges[step.step_id].append(dep)


    # ═══════════════════════════════════════════════════════════════════════════
    # Validation
    # ═══════════════════════════════════════════════════════════════════════════

    def validate(self) -> None:
        """
        Verifica que el workflow sea un DAG válido.
        """

        self.topological_sort()

        if len(self.sorted_steps) != len(self.nodes):

            raise ValueError(
                "Workflow contains cycles."
            )


    # ═══════════════════════════════════════════════════════════════════════════
    # Topological sorting
    # ═══════════════════════════════════════════════════════════════════════════

    def topological_sort(self) -> list[SimulationStep]:
        """
        Ordenamiento topológico usando Kahn's Algorithm.
        """

        in_degree = {
            node: len(self.reverse_edges[node])
            for node in self.nodes
        }

        queue = deque(
            [
                node
                for node, degree in in_degree.items()
                if degree == 0
            ]
        )

        ordered = []

        while queue:

            current = queue.popleft()

            ordered.append(
                self.nodes[current]
            )

            for neighbor in self.edges[current]:

                in_degree[neighbor] -= 1

                if in_degree[neighbor] == 0:

                    queue.append(neighbor)

        self.sorted_steps = ordered

        return ordered


    # ═══════════════════════════════════════════════════════════════════════════
    # Views
    # ═══════════════════════════════════════════════════════════════════════════

    def to_execution_view(self) -> list[SimulationStep]:
        """
        Vista completa y ordenada del workflow.
        """

        if not self.sorted_steps:
            self.topological_sort()

        return self.sorted_steps


    def to_user_view(self) -> list[str]:
        """
        Vista simplificada para usuarios.
        """

        if not self.sorted_steps:
            self.topological_sort()

        user_steps = []

        for step in self.sorted_steps:

            # ── ocultar detalles internos ────────────────────────────────────

            if step.step_id.startswith("prepare_"):
                continue

            if step.step_id.startswith("validate_pose_"):
                continue

            # ── simplificación manual ───────────────────────────────────────

            if step.step_id.startswith(
                "review_parametrization"
            ):

                user_steps.append(
                    "Manual ligand parametrization review required"
                )

                continue

            user_steps.append(step.title)

        return user_steps


    # ═══════════════════════════════════════════════════════════════════════════
    # Mermaid export
    # ═══════════════════════════════════════════════════════════════════════════

    def render_mermaid(self) -> str:
        """
        Exporta el workflow como Mermaid graph.
        """

        lines = [
            "graph TD"
        ]

        for step in self.to_execution_view():

            lines.append(
                f'    {step.step_id}["{step.title}"]'
            )

        for step in self.to_execution_view():

            for dep in step.depends_on:

                lines.append(
                    f"    {dep} --> {step.step_id}"
                )

        return "\n".join(lines)