# executors/base_executor.py
"""
Interfaz abstracta del executor.

Cada executor concreto implementa _run_step() con su propia estrategia:
    ShellExecutor    → subprocess real
    DryRunExecutor   → simula ejecución sin tocar disco ni procesos
    GROMACSExecutor  → futuro: orquestación específica de GROMACS

El executor nunca sabe qué contiene un step — solo sabe:
    - dónde está el directorio del step
    - qué script ejecutar
    - qué archivos espera encontrar después

La lógica científica vive en los builders y en el decision engine.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from datetime import datetime
import json

from executors.execution_state import (
    WorkspaceExecutionState,
    StepExecutionRecord,
    StepStatus,
)


# ═══════════════════════════════════════════════════════════════════════════════
# Base Executor
# ═══════════════════════════════════════════════════════════════════════════════

class BaseExecutor(ABC):
    """
    Executor abstracto de workspaces SimForge.

    Responsabilidades:
        - iterar steps en orden
        - gestionar estado de ejecución
        - detectar fallos
        - serializar estado a disco
        - NO interpretar contenido científico
    """

    def __init__(
        self,
        workspace_path: str | Path,
        dry_run: bool = True,
    ):
        self.workspace_path = Path(workspace_path)
        self.dry_run        = dry_run
        self.state: WorkspaceExecutionState | None = None

    # ── API pública ───────────────────────────────────────────────────────────

    def run(self) -> WorkspaceExecutionState:
        """
        Ejecuta el workspace completo en orden topológico.
        Retorna el estado final de ejecución.
        """
        self.state = self._initialize_state()
        self._log(f"Executor iniciado — dry_run={self.dry_run}")
        self._log(f"Workspace: {self.workspace_path}")

        self.state.started_at = datetime.now()

        try:
            for record in self.state.steps:

                # Verificar que las dependencias completaron
                if self._is_blocked(record):
                    record.status = StepStatus.BLOCKED
                    self._log(f"[BLOCKED]  {record.step_id}")
                    continue

                # Ejecutar
                self._log(f"[START]    {record.step_id}")
                record.status     = StepStatus.RUNNING
                record.started_at = datetime.now()
                self._save_state()

                self._run_step(record)

                record.finished_at = datetime.now()
                record.elapsed_s   = (
                    record.finished_at - record.started_at
                ).total_seconds()

                if record.status == StepStatus.RUNNING:
                    record.status = StepStatus.DONE

                self._log(
                    f"[{record.status.value.upper():8}] "
                    f"{record.step_id} "
                    f"({record.elapsed_s:.1f}s)"
                )

                self._save_state()

                # Si un step bloqueante falló → detener
                if record.status == StepStatus.FAILED:
                    step_dir = Path(record.step_dir)
                    meta     = self._read_metadata(step_dir)
                    if meta.get("blocking", False):
                        self._log(
                            f"[ABORT] Step bloqueante falló: {record.step_id}"
                        )
                        break

        except KeyboardInterrupt:
            self.state.was_interrupted = True
            self._log("[INTERRUPTED] Ejecución interrumpida por el usuario")

        finally:
            self.state.finished_at = datetime.now()
            self.state.is_complete = self.state.all_done()
            self._save_state()

        return self.state

    # ── Método abstracto — cada executor lo implementa ───────────────────────

    @abstractmethod
    def _run_step(self, record: StepExecutionRecord) -> None:
        """
        Ejecuta un step individual y escribe en record.
        Debe actualizar: record.status, record.stdout, record.stderr,
        record.exit_code, record.outputs_found, record.outputs_missing.
        """
        raise NotImplementedError

    # ── Inicialización ────────────────────────────────────────────────────────

    def _initialize_state(self) -> WorkspaceExecutionState:
        """
        Construye el WorkspaceExecutionState desde el workspace en disco.

        Orden de preferencia:
            1. execution_state.json   → reanudar ejecución previa (solo non-dry-run)
            2. execution_manifest.json → manifest generado por WorkspaceBuilder (fuente de verdad del DAG)
            3. filesystem scan         → fallback para workspaces sin manifest (backward compat)
        """
        steps_dir  = self.workspace_path / "steps"
        state_file = self.workspace_path / "execution_state.json"
        manifest_file = self.workspace_path / "metadata" / "execution_manifest.json"

        # ── 1. Reanudar si existe estado previo (solo modo real) ──────────────
        if state_file.exists() and not self.dry_run:
            raw = json.loads(state_file.read_text())
            return WorkspaceExecutionState(**raw)

        # Limpiar estado anterior si existe
        if state_file.exists():
            state_file.unlink()

        # ── 2. Manifest-driven execution ──────────────────────────────────────
        if manifest_file.exists():
            manifest   = json.loads(manifest_file.read_text())
            system_type = manifest.get("system_type")
            records: list[StepExecutionRecord] = []

            for entry in manifest["steps"]:
                step_dir = steps_dir / entry["dir_name"]
                if not step_dir.exists():
                    raise RuntimeError(
                        f"[Executor] step_dir no encontrado para '{entry['step_id']}': {step_dir}\n"
                        "El workspace puede estar corrupto. Reconstruye con WorkspaceBuilder."
                    )
                records.append(
                    StepExecutionRecord(
                        step_id    = entry["step_id"],
                        step_dir   = str(step_dir.resolve()),
                        depends_on = entry.get("depends_on", []),
                    )
                )

            self._log(
                f"Manifest cargado → {len(records)} steps "
                f"desde {manifest_file.relative_to(self.workspace_path)}"
            )
            return WorkspaceExecutionState(
                workspace_path = str(self.workspace_path.resolve()),
                system_type    = system_type,
                dry_run        = self.dry_run,
                steps          = records,
            )

        # ── 3. Fallback: filesystem scan (backward compat) ────────────────────
        self._log("[WARN] execution_manifest.json no encontrado — usando filesystem scan")

        meta_file = self.workspace_path / "metadata" / "summary.json"
        system_type = None
        if meta_file.exists():
            meta = json.loads(meta_file.read_text())
            system_type = meta.get("system_type")

        records = []
        if steps_dir.exists():
            for step_dir in sorted(steps_dir.iterdir()):
                if not step_dir.is_dir():
                    continue
                parts   = step_dir.name.split("_", 1)
                step_id = parts[1] if len(parts) == 2 else step_dir.name
                records.append(
                    StepExecutionRecord(
                        step_id  = step_id,
                        step_dir = str(step_dir.resolve()),
                    )
                )

        return WorkspaceExecutionState(
            workspace_path = str(self.workspace_path.resolve()),
            system_type    = system_type,
            dry_run        = self.dry_run,
            steps          = records,
        )

    # ── Dependencias ──────────────────────────────────────────────────────────

    def _is_blocked(self, record: StepExecutionRecord) -> bool:
        """
        Un step está bloqueado si alguno de sus predecesores directos falló.

        Si record.depends_on está poblado (manifest-driven): usa dependencias reales del DAG.
        Si está vacío (filesystem scan / workspace antiguo): fallback conservador secuencial.
        """
        if record.depends_on:
            step_status = {r.step_id: r.status for r in self.state.steps}
            return any(
                step_status.get(dep) == StepStatus.FAILED
                for dep in record.depends_on
            )

        # Fallback conservador: cualquier fallo previo bloquea los siguientes
        for other in self.state.steps:
            if other.step_id == record.step_id:
                break
            if other.status == StepStatus.FAILED:
                return True
        return False

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _log(self, message: str) -> None:
        ts  = datetime.now().strftime("%H:%M:%S")
        line = f"[{ts}] {message}"
        print(line)
        if self.state:
            self.state.log_lines.append(line)

    def _save_state(self) -> None:
        if self.state is None:
            return
        state_file = self.workspace_path / "execution_state.json"
        state_file.write_text(
            self.state.model_dump_json(indent=2)
        )

    def _read_metadata(self, step_dir: Path) -> dict:
        meta_file = step_dir / "metadata.json"
        if meta_file.exists():
            return json.loads(meta_file.read_text())
        return {}

    def _check_expected_outputs(
        self,
        step_dir: Path,
        expected: list[str],
    ) -> tuple[list[str], list[str]]:
        """
        Verifica qué archivos esperados existen post-ejecución.
        Retorna (found, missing).
        """
        found   = []
        missing = []
        for name in expected:
            path = step_dir / name
            if path.exists():
                found.append(name)
            else:
                missing.append(name)
        return found, missing
