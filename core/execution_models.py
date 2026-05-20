# core/execution_models.py

from __future__ import annotations

from enum import Enum
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

class SimulationPlan(BaseModel):

    status: PlanStatus

    inferred_system_type: str | None = None

    blocking_issues: list[BlockingIssue] = []

    steps: list[SimulationStep] = []

    special_protocols: list[str] = []

    checklist: list[CheckItem] = []

    notes: list[str] = []