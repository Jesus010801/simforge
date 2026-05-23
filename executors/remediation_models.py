# executors/remediation_models.py
"""
Modelos de datos para el loop adaptativo de remediación.

Separación de responsabilidades:
    DiagnosisResult     → qué salió mal y por qué (percepción)
    RemediationAction   → qué hacer para corregirlo (decisión)
    RemediationPlan     → conjunto de acciones para un step fallido
    RemediationRecord   → historial de un ciclo execute→diagnose→remediate→retry

Estos modelos son el contrato entre:
    signal_detector.py    → produce DiagnosisResult + RemediationPlan
    remediation_executor.py → consume RemediationPlan y escribe RemediationRecord

Invariante de diseño:
    - Los modelos NO contienen lógica
    - Los modelos son serializables (Pydantic) para persistencia en disco
    - Un RemediationRecord se añade a StepExecutionRecord.remediations
      (campo nuevo que agregaremos a execution_state.py)
"""

from __future__ import annotations

from enum import Enum
from datetime import datetime
from typing import Optional
from pydantic import BaseModel, Field


# ═══════════════════════════════════════════════════════════════════════════════
# Categorías de error conocidas (percepción)
# ═══════════════════════════════════════════════════════════════════════════════

class ErrorCategory(str, Enum):
    # GROMACS — errores de integración
    LINCS_WARNING       = "lincs_warning"       # LINCS bond angle warnings
    LINCS_FATAL         = "lincs_fatal"         # LINCS failure → crash
    NAN_ENERGY          = "nan_energy"          # NaN/Inf en energías
    EXPLODING_SYSTEM    = "exploding_system"    # velocidades/coordenadas extremas

    # GROMACS — convergencia
    FMAX_NOT_CONVERGED  = "fmax_not_converged"  # minimización no convergió (Fmax > tol)
    POOR_EQUILIBRATION  = "poor_equilibration"  # temperatura/presión inestable

    # GROMACS — topología / parámetros
    MISSING_PARAMETER   = "missing_parameter"   # parámetro FF no encontrado
    ATOM_TYPE_MISMATCH  = "atom_type_mismatch"  # tipo de átomo no reconocido
    CHARGE_IMBALANCE    = "charge_imbalance"    # carga total ≠ 0 en sistema

    # GROMACS — I/O
    MISSING_INPUT_FILE  = "missing_input_file"  # archivo de entrada no encontrado
    CORRUPT_CHECKPOINT  = "corrupt_checkpoint"  # .cpt corrupto

    # Ejecución general
    INTERACTIVE_BLOCK   = "interactive_block"   # proceso bloqueado esperando input TTY
    TIMEOUT             = "timeout"             # step excedió límite de tiempo
    NONZERO_EXIT        = "nonzero_exit"        # exit code ≠ 0 sin categoría específica
    MISSING_OUTPUT      = "missing_output"      # outputs esperados no generados

    # Desconocido
    UNKNOWN             = "unknown"


class ErrorSeverity(str, Enum):
    RECOVERABLE   = "recoverable"   # se puede corregir automáticamente
    NEEDS_REVIEW  = "needs_review"  # remediación posible pero requiere confirmación
    FATAL         = "fatal"         # no se puede remediar sin intervención humana


# ═══════════════════════════════════════════════════════════════════════════════
# Diagnóstico
# ═══════════════════════════════════════════════════════════════════════════════

class DiagnosisResult(BaseModel):
    """
    Resultado del análisis de un step fallido.

    Producido por AdaptiveReasoner.diagnose().
    Describe exactamente qué ocurrió y con qué certeza.
    """
    step_id:          str
    step_dir:         str

    # Clasificación
    category:         ErrorCategory
    severity:         ErrorSeverity
    confidence:       float           # 0.0 – 1.0

    # Evidencia
    primary_signal:   str             # línea o patrón que disparó el diagnóstico
    evidence_lines:   list[str] = []  # líneas relevantes del log
    exit_code:        Optional[int] = None

    # Contexto
    stage:            str = ""        # stage del step (minimization, equilibration, ...)
    engine:           str = ""        # engine del step (gromacs, ...)

    # Texto libre
    explanation:      str = ""        # por qué falló (para el usuario)
    reasoning:        str = ""        # trazabilidad interna del diagnóstico

    diagnosed_at:     datetime = Field(default_factory=datetime.now)


# ═══════════════════════════════════════════════════════════════════════════════
# Acciones de remediación
# ═══════════════════════════════════════════════════════════════════════════════

class ActionType(str, Enum):
    PATCH_MDP          = "patch_mdp"          # modificar parámetro en archivo .mdp
    PATCH_SCRIPT       = "patch_script"       # modificar línea en script .sh
    WRITE_FILE         = "write_file"         # escribir archivo nuevo (override completo)
    DELETE_FILE        = "delete_file"        # eliminar archivo corrupto
    COPY_FILE          = "copy_file"          # copiar archivo desde otra ubicación
    RESET_STEP         = "reset_step"         # limpiar outputs del step para retry
    INJECT_RESTRAINTS  = "inject_restraints"  # agregar position restraints al MDP
    SCALE_TIMESTEP     = "scale_timestep"     # reducir dt en MDP
    REDUCE_TEMPERATURE = "reduce_temperature" # bajar ref_t inicial en MDP
    LOG_ONLY           = "log_only"           # solo registrar, sin modificar archivos


class RemediationAction(BaseModel):
    """
    Acción atómica de corrección.

    Cada acción modifica exactamente un aspecto del workspace.
    El RemediationExecutor las aplica en orden secuencial.
    """
    action_type:    ActionType
    description:    str           # qué hace esta acción (para el usuario)

    # Parámetros por tipo de acción
    target_file:    Optional[str] = None    # path relativo al step_dir
    patch_key:      Optional[str] = None    # clave MDP a modificar (ej: "dt")
    patch_value:    Optional[str] = None    # nuevo valor (ej: "0.001")
    patch_old:      Optional[str] = None    # valor anterior (para log)
    file_content:   Optional[str] = None    # contenido completo (WRITE_FILE)
    source_file:    Optional[str] = None    # origen (COPY_FILE)

    # Metadatos
    rationale:      str = ""      # por qué esta acción es la correcta
    confidence:     float = 1.0   # confianza en que esta acción es apropiada
    is_reversible:  bool  = True   # ¿se puede deshacer?


# ═══════════════════════════════════════════════════════════════════════════════
# Plan de remediación
# ═══════════════════════════════════════════════════════════════════════════════

class RemediationPlan(BaseModel):
    """
    Conjunto de acciones para remediar un step fallido.

    Producido por AdaptiveReasoner.plan_remediation().
    Consumido por RemediationExecutor.apply().
    """
    step_id:          str
    diagnosis:        DiagnosisResult

    actions:          list[RemediationAction] = []
    is_applicable:    bool  = True    # False si el error es fatal
    requires_human:   bool  = False   # True si necesita confirmación
    max_retries:      int   = 2       # cuántas veces intentar este plan

    # Texto libre
    strategy:         str = ""        # descripción de la estrategia
    expected_outcome: str = ""        # qué debería pasar si funciona

    planned_at:       datetime = Field(default_factory=datetime.now)


# ═══════════════════════════════════════════════════════════════════════════════
# Registro histórico de un ciclo de remediación
# ═══════════════════════════════════════════════════════════════════════════════

class RemediationStatus(str, Enum):
    PENDING   = "pending"    # plan creado, no aplicado
    APPLIED   = "applied"    # acciones aplicadas, retry en curso
    SUCCEEDED = "succeeded"  # retry exitoso después de remediación
    FAILED    = "failed"     # retry falló incluso después de remediación
    SKIPPED   = "skipped"    # no aplicable o fatal, sin retry
    ABORTED   = "aborted"    # cancelado por el usuario


class RemediationRecord(BaseModel):
    """
    Registro persistente de un ciclo completo de remediación.

    Se añade a StepExecutionRecord.remediations[].
    Permite reconstruir exactamente qué se intentó y qué pasó.
    """
    attempt_number:   int             # 1-based: primer intento = 1
    plan:             RemediationPlan
    status:           RemediationStatus = RemediationStatus.PENDING

    # Timing
    started_at:       Optional[datetime] = None
    finished_at:      Optional[datetime] = None
    elapsed_s:        Optional[float]    = None

    # Resultado del retry
    retry_exit_code:  Optional[int]  = None
    retry_stdout:     str            = ""
    retry_stderr:     str            = ""

    # Qué archivos se modificaron
    files_modified:   list[str]      = []
    files_backed_up:  list[str]      = []

    # Notas del executor
    executor_notes:   list[str]      = []