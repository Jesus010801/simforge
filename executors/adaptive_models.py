# executors/adaptive_models.py
"""
Modelos de output del adaptive reasoning engine.

Estos modelos representan las conclusiones del reasoning post-ejecución:
qué salió mal, qué se puede remediar automáticamente, qué requiere
intervención humana, y si el pipeline puede continuar.

Árbol de datos:

    AdaptiveReasoningResult
    │
    ├── verdict              ReasoningVerdict   — continuar / revisar / abortar
    ├── step_analyses        list[StepAnalysis] — análisis por step GROMACS
    ├── remediation_plan     RemediationPlan    — qué hacer para corregir
    └── notes / warnings / errors               — resumen textual

Principio arquitectónico:
    - Los modelos son Pydantic puros — sin lógica
    - La lógica vive en adaptive_reasoning.py
    - El decision engine original NO se modifica — este es un segundo
      nivel de reasoning que ocurre POST-ejecución
    - Compatible hacia atrás: si no hay diagnósticos GROMACS, el
      reasoning produce un resultado vacío sin errores
"""

from __future__ import annotations

from enum import Enum
from typing import Optional
from pydantic import BaseModel

from core.models import Severity


# ═══════════════════════════════════════════════════════════════════════════════
# Veredicto global
# ═══════════════════════════════════════════════════════════════════════════════

class ReasoningVerdict(str, Enum):

    CONTINUE    = "continue"
    # Todos los steps completaron sin issues críticos.
    # El pipeline puede continuar al siguiente stage.

    REVIEW      = "review"
    # Hay advertencias que el usuario debería revisar antes de continuar.
    # El pipeline puede continuar pero con supervisión.

    REMEDIATE   = "remediate"
    # Se detectaron problemas corregibles. El reasoning generó un
    # RemediationPlan con acciones concretas. Ejecutar el plan
    # y re-run antes de continuar.

    ABORT       = "abort"
    # Error crítico que no puede remediarse automáticamente.
    # Requiere intervención manual antes de cualquier re-ejecución.


# ═══════════════════════════════════════════════════════════════════════════════
# Análisis por step
# ═══════════════════════════════════════════════════════════════════════════════

class StepAnalysisVerdict(str, Enum):

    OK              = "ok"
    NOT_CONVERGED   = "not_converged"
    NEEDS_REVIEW    = "needs_review"
    REMEDIABLE      = "remediable"
    FATAL           = "fatal"


class StepAnalysis(BaseModel):
    """
    Análisis de un step individual post-ejecución.

    Consume GROMACSStepDiagnostic y produce una interpretación
    contextual — el diagnóstico dice qué pasó, el análisis dice
    qué significa y qué hacer.
    """

    step_id:        str
    stage:          str
    verdict:        StepAnalysisVerdict = StepAnalysisVerdict.OK
    severity:       Severity            = Severity.LOW

    # Interpretación científica del diagnóstico
    interpretation: str   = ""

    # Qué acción se recomienda (texto libre para el usuario)
    recommended_action: str = ""

    # Si existe una remediación automática posible
    has_remediation:    bool = False
    remediation_id:     Optional[str] = None   # referencia a RemediationStep

    # Métricas clave extraídas del diagnóstico (para display)
    key_metrics:        dict = {}

    # Notas adicionales
    notes:              list[str] = []


# ═══════════════════════════════════════════════════════════════════════════════
# Plan de remediación
# ═══════════════════════════════════════════════════════════════════════════════

class RemediationTarget(str, Enum):
    """Qué tipo de artefacto modifica esta remediación."""

    MDP_FILE        = "mdp_file"        # parámetro en .mdp
    TOPOLOGY        = "topology"        # topol.top o .itp
    STRUCTURE       = "structure"       # .gro o .pdb
    SCRIPT          = "script"          # run.sh
    MANUAL_ACTION   = "manual_action"   # requiere intervención humana
    CHECKPOINT      = "checkpoint"      # reanudar desde .cpt


class RemediationStep(BaseModel):
    """
    Acción concreta de remediación para un problema detectado.

    Una remediación puede ser:
        - automática: el sistema puede aplicarla sin intervención humana
        - semi-automática: genera el cambio pero el usuario debe confirmar
        - manual: solo describe qué hacer, el usuario ejecuta

    El adaptive reasoning NO aplica remediaciones — las propone.
    Un futuro RemediationExecutor podría aplicarlas automáticamente.
    """

    remediation_id:     str
    step_id:            str        # step al que aplica
    target:             RemediationTarget
    priority:           Severity   # HIGH = aplicar primero

    # Descripción del problema
    problem:            str
    root_cause:         str = ""

    # Acción propuesta
    action:             str        # descripción legible
    automatic:          bool = False   # True si puede aplicarse sin confirmación

    # Para remediaciones de MDP: parámetro específico a cambiar
    mdp_parameter:      Optional[str]   = None
    mdp_current_value:  Optional[str]   = None
    mdp_suggested_value: Optional[str]  = None

    # Para remediaciones de script
    script_change:      Optional[str]   = None

    # Referencia al archivo afectado
    target_file:        Optional[str]   = None

    # Notas adicionales
    notes:              list[str] = []


class RemediationPlan(BaseModel):
    """
    Plan completo de remediación post-ejecución.

    Generado por adaptive_reasoning.py cuando hay problemas corregibles.
    Contiene pasos ordenados por prioridad.
    """

    steps:              list[RemediationStep] = []

    # Resumen ejecutivo
    n_automatic:        int  = 0   # pasos que pueden aplicarse solos
    n_manual:           int  = 0   # pasos que requieren intervención
    estimated_effort:   str  = "unknown"   # "minutes" / "hours" / "days"

    # ¿Tiene sentido re-ejecutar después de aplicar este plan?
    rerun_recommended:  bool = True
    rerun_from_step:    Optional[str] = None   # step_id desde donde reiniciar

    @property
    def is_empty(self) -> bool:
        return len(self.steps) == 0

    def sorted_steps(self) -> list[RemediationStep]:
        """Pasos ordenados: HIGH primero, luego MEDIUM, luego LOW."""
        order = {Severity.HIGH: 0, Severity.MEDIUM: 1, Severity.LOW: 2}
        return sorted(self.steps, key=lambda s: order.get(s.priority, 3))


# ═══════════════════════════════════════════════════════════════════════════════
# Resultado completo del adaptive reasoning
# ═══════════════════════════════════════════════════════════════════════════════

class AdaptiveReasoningResult(BaseModel):
    """
    Output completo del adaptive reasoning engine.

    Contrato entre adaptive_reasoning.py y cualquier capa que lo consuma
    (UI, CLI, futuro orchestrator).

    Se construye desde:
        - WorkspaceExecutionState   (qué pasó en la ejecución)
        - dict[GROMACSStepDiagnostic] (qué detectó el executor GROMACS)
        - SystemState opcional      (contexto del sistema original)
    """

    # Veredicto global
    verdict:            ReasoningVerdict = ReasoningVerdict.CONTINUE

    # Análisis por step
    step_analyses:      list[StepAnalysis]  = []

    # Plan de remediación (vacío si verdict == CONTINUE o ABORT)
    remediation_plan:   RemediationPlan     = RemediationPlan()

    # Resumen textual
    summary:            str   = ""
    notes:              list[str] = []
    warnings:           list[str] = []
    errors:             list[str] = []

    # Contexto de la ejecución analizada
    workspace_path:     str   = ""
    system_type:        Optional[str] = None
    n_steps_analyzed:   int   = 0
    n_steps_ok:         int   = 0
    n_steps_failed:     int   = 0

    @property
    def has_critical_errors(self) -> bool:
        return self.verdict == ReasoningVerdict.ABORT

    @property
    def can_continue(self) -> bool:
        return self.verdict in (
            ReasoningVerdict.CONTINUE,
            ReasoningVerdict.REVIEW,
        )
