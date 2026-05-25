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


class AutomationLevel(str, Enum):
    """
    Grado de automatización de un step de simulación.

    MANUAL        — El usuario hace todo a mano; SimForge no genera scripts.
    GUIDED        — SimForge genera instrucciones (README.md), el usuario las
                    sigue manualmente. Sin script ejecutable.
    SEMI_AUTOMATED — SimForge genera y ejecuta el script, pero puede requerir
                    confirmación o input del usuario en runtime.
    AUTOMATED     — Totalmente automático; no requiere intervención del usuario.

    El runtime decide si ejecutar o saltar un step basándose en este campo.
    Cuando no está presente en metadata.json, cae de vuelta al campo step_type
    legado para compatibilidad con workspaces antiguos.
    """
    MANUAL         = "manual"
    GUIDED         = "guided"
    SEMI_AUTOMATED = "semi_automated"
    AUTOMATED      = "automated"

    # Mapeado de AutomationLevel → step_type legado para serialización cruzada
    @classmethod
    def from_step_type(cls, step_type: "StepType") -> "AutomationLevel":
        """Derive automation level from legacy StepType."""
        return {
            StepType.AUTOMATIC:  cls.AUTOMATED,
            StepType.MANUAL:     cls.GUIDED,
            StepType.EXTERNAL:   cls.MANUAL,
            StepType.VALIDATION: cls.GUIDED,
        }.get(step_type, cls.AUTOMATED)

    @property
    def needs_user(self) -> bool:
        """True when the step cannot run unattended."""
        return self in (self.MANUAL, self.GUIDED)

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

    MEMBRANE_EMBEDDING = "membrane_embedding"

# ═══════════════════════════════════════════════════════════════════════════════
# Step individual del workflow
# ═══════════════════════════════════════════════════════════════════════════════

class SimulationStep(BaseModel):

    step_id: str

    title: str

    stage: StepStage

    step_type: StepType = StepType.AUTOMATIC

    # automation_level is the authoritative field for runtime skip decisions.
    # When set explicitly it overrides the legacy step_type heuristic.
    # Defaults to None so pipelines that haven't been updated yet continue to
    # work via the step_type fallback in the executor.
    automation_level: Optional[AutomationLevel] = None

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

    def effective_automation_level(self) -> AutomationLevel:
        """Return automation_level, falling back to step_type derivation."""
        if self.automation_level is not None:
            return self.automation_level
        return AutomationLevel.from_step_type(self.step_type)


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
    hardware:              str   = "auto"        # "auto" | "gpu" | "cpu"
    # ── Derived from WorkflowHints ─────────────────────────────────────────────
    semiisotropic_coupling: bool = False   # semiisotropic NPT (membrane systems)
    membrane_required:      bool = False   # system embeds a bilayer
    extended_equilibration: bool = False   # longer NVT/NPT phase
    extended_production:    bool = False   # longer production run


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