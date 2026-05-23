# executors/execution_state.py
"""
Modelos de estado de ejecución en tiempo real.

Estos modelos representan el estado de un workspace mientras se ejecuta,
no el plan (que vive en execution_models.py).

La distinción es importante:
    SimulationPlan     → qué hay que hacer (intención)
    ExecutionState     → qué pasó realmente (realidad)

El executor escribe en ExecutionState conforme avanza.
El adaptive reasoning lee desde ExecutionState para decidir si continuar.

Campos añadidos en StepExecutionRecord:

    gromacs_diagnostic : Optional[dict]
        Serialización de GROMACSStepDiagnostic producida por GROMACSExecutor
        después de cada step. Se guarda como dict para que execution_state.json
        sea completamente serializable sin importar GROMACSStepDiagnostic aquí
        (evita dependencia circular executors → executors).
        El RemediationExecutor lo reconstruye cuando necesita:
            GROMACSStepDiagnostic(**record.gromacs_diagnostic)

    remediations : list[RemediationRecord]
        Historial completo de ciclos diagnose→plan→apply→retry para este step.
        Reemplaza el dict global _REMEDIATION_HISTORY que existía en
        remediation_executor.py. El estado es ahora completamente local
        al record y se persiste con él en execution_state.json.
"""

from __future__ import annotations

from enum import Enum
from datetime import datetime
from typing import Any, Optional
from pydantic import BaseModel, Field


# ═══════════════════════════════════════════════════════════════════════════════
# Estado de un step individual
# ═══════════════════════════════════════════════════════════════════════════════

class StepStatus(str, Enum):

    PENDING     = "pending"      # no ha corrido
    RUNNING     = "running"      # corriendo ahora
    DONE        = "done"         # terminó exitosamente
    FAILED      = "failed"       # terminó con error
    SKIPPED     = "skipped"      # saltado por decisión del executor
    BLOCKED     = "blocked"      # no puede correr por dependencia fallida
    RECOVERABLE = "recoverable"  # interrumpido con checkpoint disponible


class StepExecutionRecord(BaseModel):
    """
    Registro completo de la ejecución de un step individual.

    Es la unidad de estado que el executor actualiza en tiempo real
    y que el reasoner lee para decidir si remediar.

    Invariante de serialización:
        Todos los campos deben ser serializables a JSON sin pérdida.
        gromacs_diagnostic se almacena como dict (no como modelo tipado)
        para evitar dependencias circulares en el módulo de estado.
        remediations se serializa inline como lista de dicts Pydantic.
    """

    step_id:     str
    step_dir:    str                     # path absoluto al directorio del step
    depends_on:  list[str] = []          # step_ids que deben completar antes
    status:      StepStatus = StepStatus.PENDING

    # ── Timing ────────────────────────────────────────────────────────────────
    started_at:  Optional[datetime] = None
    finished_at: Optional[datetime] = None
    elapsed_s:   Optional[float]    = None

    # ── Output de ejecución ───────────────────────────────────────────────────
    stdout:      str            = ""
    stderr:      str            = ""
    exit_code:   Optional[int]  = None

    # ── Diagnóstico de ejecución ──────────────────────────────────────────────
    error_message:   Optional[str] = None
    retry_count:     int           = 0

    # ── Archivos generados ────────────────────────────────────────────────────
    outputs_found:   list[str] = []   # archivos que existen post-ejecución
    outputs_missing: list[str] = []   # archivos esperados que no aparecieron

    # ── Percepción rica de GROMACS ────────────────────────────────────────────
    # Serialización de GROMACSStepDiagnostic producida por GROMACSExecutor.
    # Se almacena como dict para mantener execution_state.py libre de imports
    # del executor GROMACS (evita dependencia circular).
    #
    # Escritura (en GROMACSExecutor._run_step):
    #     record.gromacs_diagnostic = diag.model_dump()
    #
    # Lectura (en RemediationExecutor.diagnose_and_plan):
    #     if record.gromacs_diagnostic:
    #         rich_diag = GROMACSStepDiagnostic(**record.gromacs_diagnostic)
    gromacs_diagnostic: Optional[dict[str, Any]] = None

    # ── Historial de remediación ──────────────────────────────────────────────
    # Registro de cada ciclo diagnose→plan→apply→retry para este step.
    # Reemplaza el dict global _REMEDIATION_HISTORY de remediation_executor.py.
    #
    # Invariante: se añade un RemediationRecord por cada intento,
    # independientemente de si tuvo éxito o no.
    # El campo se serializa como lista de dicts para JSON completo.
    remediations: list[dict[str, Any]] = Field(default_factory=list)

    # ── Helpers de remediación ────────────────────────────────────────────────

    def n_remediations(self) -> int:
        """Número de ciclos de remediación intentados."""
        return len(self.remediations)

    def last_remediation(self) -> Optional[dict[str, Any]]:
        """Último RemediationRecord serializado, o None si no hay ninguno."""
        return self.remediations[-1] if self.remediations else None

    def add_remediation(self, record: Any) -> None:
        """
        Añade un RemediationRecord al historial.

        Acepta tanto el modelo Pydantic como el dict serializado
        para no forzar el import de remediation_models aquí.
        """
        if hasattr(record, "model_dump"):
            self.remediations.append(record.model_dump())
        else:
            self.remediations.append(record)


# ═══════════════════════════════════════════════════════════════════════════════
# Estado global del workspace en ejecución
# ═══════════════════════════════════════════════════════════════════════════════

class WorkspaceExecutionState(BaseModel):
    """
    Estado completo de la ejecución de un workspace.

    Se serializa a execution_state.json dentro del workspace
    para poder reanudar ejecuciones interrumpidas.
    """

    workspace_path:  str
    system_type:     Optional[str] = None

    # ── Estado global ─────────────────────────────────────────────────────────
    started_at:      Optional[datetime] = None
    finished_at:     Optional[datetime] = None
    is_complete:     bool               = False
    was_interrupted: bool               = False

    # ── Steps ─────────────────────────────────────────────────────────────────
    steps:           list[StepExecutionRecord] = []

    # ── Modo de ejecución ─────────────────────────────────────────────────────
    dry_run:         bool = True   # True = no ejecuta comandos reales

    # ── Log global ────────────────────────────────────────────────────────────
    log_lines:       list[str] = Field(default_factory=list)

    # ── Helpers de consulta ───────────────────────────────────────────────────

    def get_step(self, step_id: str) -> Optional[StepExecutionRecord]:
        for s in self.steps:
            if s.step_id == step_id:
                return s
        return None

    def n_done(self) -> int:
        return sum(1 for s in self.steps if s.status == StepStatus.DONE)

    def n_failed(self) -> int:
        return sum(1 for s in self.steps if s.status == StepStatus.FAILED)

    def n_pending(self) -> int:
        return sum(1 for s in self.steps if s.status == StepStatus.PENDING)

    def n_remediated(self) -> int:
        """Steps que tuvieron al menos un ciclo de remediación."""
        return sum(1 for s in self.steps if s.remediations)

    def has_failures(self) -> bool:
        return self.n_failed() > 0

    def all_done(self) -> bool:
        return all(
            s.status in (StepStatus.DONE, StepStatus.SKIPPED)
            for s in self.steps
        )

    def steps_with_gromacs_diagnostic(self) -> list[StepExecutionRecord]:
        """Steps que tienen diagnóstico rico de GROMACS disponible."""
        return [s for s in self.steps if s.gromacs_diagnostic is not None]