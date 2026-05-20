# core/compiler.py

from __future__ import annotations

from core.models import SystemState

from core.parser import parse_yaml

from core.decision_engine import (
    build_simulation_plan,
)

from core.compiler_models import (
    CompilationResult,
)

from workflows.workflow_graph import (
    WorkflowGraph,
)

from pipelines.md_pipeline import (
    MDPipeline,
)

from pipelines.inhibition_pipeline import (
    InhibitionPipeline,
)



# ═══════════════════════════════════════════════════════════════════════════════
# Simulation Compiler
# ═══════════════════════════════════════════════════════════════════════════════

class SimulationCompiler:
    """
    API pública principal de SimForge.

    Pipeline:

        YAML
          ↓
        SystemState
          ↓
        SimulationPlan
          ↓
        WorkflowGraph
          ↓
        CompilationResult
    """

    def compile(
        self,
        yaml_path: str,
    ) -> CompilationResult:

        # ────────────────────────────────────────────────────────────────────
        # Parse YAML
        # ────────────────────────────────────────────────────────────────────

        state = parse_yaml(yaml_path)

        # ────────────────────────────────────────────────────────────────────
        # Build semantic plan
        # ────────────────────────────────────────────────────────────────────

        pipeline = self._select_pipeline(
            state
        )

        plan = pipeline.build_plan(
            state
        )

        # ────────────────────────────────────────────────────────────────────
        # Build workflow graph
        # ────────────────────────────────────────────────────────────────────

        graph = WorkflowGraph(plan)

        graph.validate()
        execution_order = (
            graph.to_execution_view()
        )

        # ────────────────────────────────────────────────────────────────────
        # Build result
        # ────────────────────────────────────────────────────────────────────

        summary = []

        summary.append(
            f"System type: {state.inferred_system_type}"
        )

        summary.append(
            f"Workflow steps: {len(plan.steps)}"
        )

        summary.append(
            f"Blocking issues: {len(plan.blocking_issues)}"
        )

        summary.append(
            f"Special protocols: {len(plan.special_protocols)}"
        )

        return CompilationResult(
            state=state,

            plan=plan,
            
            execution_order=execution_order,

            user_view=graph.to_user_view(),

            mermaid_graph=graph.render_mermaid(),

            workflow_valid=True,

            summary=summary,
        )
    
    def _select_pipeline(
        self,
        state: SystemState,
    ):

    # ────────────────────────────────────────────────────────────────────
    # Inhibition systems
    # ────────────────────────────────────────────────────────────────────

        if state.inferred_system_type == (
        "competitive-inhibition"
        ):

         return InhibitionPipeline()

    # ────────────────────────────────────────────────────────────────────
    # Default MD
    # ────────────────────────────────────────────────────────────────────

        return MDPipeline()