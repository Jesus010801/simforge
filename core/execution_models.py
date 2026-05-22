# core/execution_models.py

from __future__ import annotations

from enum import Enum
from typing import Any, Optional
from pydantic import BaseModel

from core.models import Severity


# ═══════════════════════════════════════════════════════════════════════════════
# Estado global del plan
# ═══════════════════════════════════════════════════════════════════════════════

class PlanStatus(str, Enum):

    READY = "ready"

    BLOCKED = "blocked"

    NEEDS_REVIEW = "needs_review"

    EXPERIMENTAL = "experimental"


# ═══════════════════════════════════════════════════════════════════════════════
# Etapas del pipeline de simulación
# ═══════════════════════════════════════════════════════════════════════════════

class StepType(str, Enum):

    AUTOMATIC = "automatic"

    MANUAL = "manual"

    EXTERNAL = "external"

    VALIDATION = "validation"

class StepStage(str, Enum):

    PREPARATION = "preparation"

    VALIDATION = "validation"

    PARAMETRIZATION = "parametrization"

    ASSEMBLY = "assembly"

    MINIMIZATION = "minimization"

    EQUILIBRATION = "equilibration"

    ENHANCED_SAMPLING = "enhanced_sampling"

    PRODUCTION = "production"

    ANALYSIS = "analysis"

# ═══════════════════════════════════════════════════════════════════════════════
# Step individual del workflow
# ═══════════════════════════════════════════════════════════════════════════════

class SimulationStep(BaseModel):

    step_id: str

    title: str

    stage: StepStage

    step_type: StepType = StepType.AUTOMATIC

    engine: str

    target_components: list[str] = []

    required: bool = True

    blocking: bool = False

    depends_on: list[str] = []

    notes: list[str] = []

    # IR payload: parámetros engine-specific generados por decision_engine.
    # Los builders leen desde aquí en lugar de hardcodear valores.
    params: dict[str, Any] = {}

    # Condición semántica que causó la inclusión de este step.
    # Metadata de auditoría — no afecta ejecución.
    condition: Optional[str] = None


# ═══════════════════════════════════════════════════════════════════════════════
# Problemas bloqueantes
# ═══════════════════════════════════════════════════════════════════════════════

class BlockingIssue(BaseModel):

    source: str

    message: str

    severity: Severity

    resolution: str | None = None


# ═══════════════════════════════════════════════════════════════════════════════
# Checklist previo a ejecución
# ═══════════════════════════════════════════════════════════════════════════════

class CheckItem(BaseModel):

    message: str

    required: bool = True

    component: str | None = None


# ═══════════════════════════════════════════════════════════════════════════════
# Plan global de simulación
# ═══════════════════════════════════════════════════════════════════════════════

class WorkflowPolicy(BaseModel):
    """
    Decisiones científicas a nivel del pipeline completo.

    Generada por decision_engine desde SystemState.
    Fuente de verdad para duración, temperatura, presión y
    estrategia de sampling — ningún builder hardcodea estos valores.
    """

    production_time_ns:    float = 10.0
    equilibration_time_ns: float = 0.1
    minimization_steps:    int   = 50_000
    temperature_K:         float = 300.0
    pressure_bar:          float = 1.0
    timestep_ps:           float = 0.002
    enhanced_sampling:     bool  = False
    sampling_method:       str   = "standard"   # "standard" | "REST2" | "metadynamics"


class SimulationPlan(BaseModel):

    status: PlanStatus

    inferred_system_type: str | None = None

    # Política de workflow: decisiones científicas globales del pipeline.
    workflow_policy: WorkflowPolicy = WorkflowPolicy()

    blocking_issues: list[BlockingIssue] = []

    steps: list[SimulationStep] = []

    special_protocols: list[str] = []

    checklist: list[CheckItem] = []

    notes: list[str] = []