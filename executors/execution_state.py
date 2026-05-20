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
"""

from __future__ import annotations

from enum import Enum
from datetime import datetime
from pathlib import Path
from pydantic import BaseModel, Field


# ═══════════════════════════════════════════════════════════════════════════════
# Estado de un step individual
# ═══════════════════════════════════════════════════════════════════════════════

class StepStatus(str, Enum):

    PENDING   = "pending"    # no ha corrido
    RUNNING   = "running"    # corriendo ahora
    DONE      = "done"       # terminó exitosamente
    FAILED    = "failed"     # terminó con error
    SKIPPED   = "skipped"    # saltado por decisión del executor
    BLOCKED   = "blocked"    # no puede correr por dependencia fallida


class StepExecutionRecord(BaseModel):
    """
    Registro completo de la ejecución de un step individual.
    """

    step_id:     str
    step_dir:    str                    # path absoluto al directorio del step
    status:      StepStatus = StepStatus.PENDING

    # Timing
    started_at:  datetime | None = None
    finished_at: datetime | None = None
    elapsed_s:   float | None    = None

    # Output
    stdout:      str  = ""
    stderr:      str  = ""
    exit_code:   int | None = None

    # Diagnóstico
    error_message:  str | None = None
    retry_count:    int        = 0

    # Archivos generados confirmados
    outputs_found:  list[str] = []     # archivos que existen post-ejecución
    outputs_missing: list[str] = []    # archivos esperados que no aparecieron


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
    system_type:     str | None = None

    # Estado global
    started_at:      datetime | None = None
    finished_at:     datetime | None = None
    is_complete:     bool            = False
    was_interrupted: bool            = False

    # Steps
    steps:           list[StepExecutionRecord] = []

    # Modo de ejecución
    dry_run:         bool = True       # True = no ejecuta comandos reales

    # Log global
    log_lines:       list[str] = Field(default_factory=list)

    # ── Helpers ──────────────────────────────────────────────────────────────

    def get_step(self, step_id: str) -> StepExecutionRecord | None:
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

    def has_failures(self) -> bool:
        return self.n_failed() > 0

    def all_done(self) -> bool:
        return all(
            s.status in (StepStatus.DONE, StepStatus.SKIPPED)
            for s in self.steps
        )
