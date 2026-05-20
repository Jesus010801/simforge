# core/compiler_models.py

from __future__ import annotations

from pydantic import BaseModel

from core.models import SystemState

from core.execution_models import (
    SimulationPlan,
)

from workflows.workflow_graph import (
    WorkflowGraph,
)

from core.execution_models import (
    SimulationStep,
)

# ═══════════════════════════════════════════════════════════════════════════════
# Compilation Result
# ═══════════════════════════════════════════════════════════════════════════════

class CompilationResult(BaseModel):
    """
    Resultado completo de compilación de SimForge.

    Contiene todas las capas derivadas desde YAML.
    """

    state: SystemState

    plan: SimulationPlan
    
    execution_order: list[SimulationStep]

    user_view: list[str]

    mermaid_graph: str

    workflow_valid: bool = True

    summary: list[str] = []