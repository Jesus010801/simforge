# executors/remediation_executor.py
"""
RemediationExecutor — ejecutor del loop adaptativo de SimForge.

Loop completo:
    execute → diagnose → plan → apply → retry → [repeat | escalate]

Responsabilidades:
    apply(plan)         → aplica las acciones del RemediationPlan al workspace
    retry(record)       → relanza el step afectado via ShellExecutor
    run_adaptive(state) → orquesta el loop completo sobre un WorkspaceExecutionState

Separación de responsabilidades:
    AdaptiveReasoner    → qué falló y qué hacer (NO toca disco)
    RemediationExecutor → aplica los cambios y relanza (SÍ toca disco)
    ShellExecutor       → ejecuta el script (sin saber de remediación)

Detalles de aplicación por ActionType:

    PATCH_MDP
        Lee el .mdp, encuentra la línea con `patch_key = valor_anterior`,
        la reemplaza con `patch_key = patch_value`.
        Guarda un backup .mdp.bak antes de modificar.

    PATCH_SCRIPT
        Igual que PATCH_MDP pero sobre archivos .sh.

    WRITE_FILE
        Escribe file_content completo en target_file.
        Guarda backup si el archivo ya existe.

    DELETE_FILE
        Elimina target_file. Guarda backup antes.

    COPY_FILE
        Copia source_file → target_file dentro del step_dir.

    RESET_STEP
        Elimina outputs conocidos del step (*.xtc, *.edr, *.log, *.cpt,
        *.gro de output, *.tpr) pero preserva inputs y scripts.

    INJECT_RESTRAINTS
        Caso especial de PATCH_MDP: agrega o modifica el campo `define`.

    SCALE_TIMESTEP
        Caso especial de PATCH_MDP sobre el campo `dt`.

    REDUCE_TEMPERATURE
        Caso especial de PATCH_MDP sobre `ref_t`.

    LOG_ONLY
        Solo registra la acción, no toca archivos.

Persistencia:
    Cada RemediationRecord se añade a record.remediations (campo nuevo
    en StepExecutionRecord via execution_state.py actualizado).
    El estado completo se serializa en execution_state.json.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Optional

from executors.execution_state import (
    WorkspaceExecutionState,
    StepExecutionRecord,
    StepStatus,
)
from executors.remediation_models import (
    RemediationPlan,
    RemediationAction,
    RemediationRecord,
    RemediationStatus,
    ActionType,
)
from executors.adaptive_reasoner import AdaptiveReasoner


# ═══════════════════════════════════════════════════════════════════════════════
# Aplicación de acciones individuales
# ═══════════════════════════════════════════════════════════════════════════════

def _backup(path: Path) -> Path:
    """Crea backup .bak del archivo. Retorna el path del backup."""
    backup_path = path.with_suffix(path.suffix + ".bak")
    shutil.copy2(path, backup_path)
    return backup_path


def _patch_text_file(
    path:      Path,
    key:       str,
    new_value: str,
    old_value: Optional[str] = None,
) -> bool:
    """
    Parcha un archivo de texto plano (MDP o script) reemplazando
    la primera ocurrencia de `key = valor` por `key = new_value`.

    Formato soportado:
        MDP:    key = value          (con espacios opcionales)
        script: KEY=value            (sin espacios, para variables de shell)

    Retorna True si el patch tuvo éxito.
    """
    if not path.exists():
        return False

    content = path.read_text()
    original = content

    # Patrón MDP: `key = valor_anterior` (con cualquier whitespace)
    mdp_pattern = re.compile(
        rf"^(\s*{re.escape(key)}\s*=\s*)(.+?)(\s*(?:;.*)?)$",
        re.MULTILINE,
    )

    if mdp_pattern.search(content):
        # Preservar comentarios inline al final de la línea
        content = mdp_pattern.sub(
            lambda m: f"{m.group(1)}{new_value}{m.group(3)}",
            content,
            count=1,
        )
    else:
        # Clave no existe en el archivo → agregar al final de la sección
        content += f"\n{key}                    = {new_value}\n"

    if content == original:
        return False   # sin cambios

    path.write_text(content)
    return True


def _reset_step_outputs(step_dir: Path) -> list[str]:
    """
    Elimina outputs del intento fallido preservando inputs y scripts.

    Elimina: *.xtc, *.edr, *.log, *.cpt, *.trr, *_out.gro, *.tpr
    Preserva: *.mdp, *.sh, *.itp, *.top, metadata.json, README.md,
              *_in.gro, *.pdb
    """
    _REMOVE_EXTENSIONS  = {".xtc", ".edr", ".log", ".cpt", ".trr", ".tpr"}
    _REMOVE_PATTERNS    = re.compile(r"(^md\.|^nvt\.|^npt\.|^em\.)")
    _PRESERVE_EXTENSIONS = {".mdp", ".sh", ".itp", ".top", ".pdb", ".bak", ".dat"}
    _PRESERVE_NAMES      = {"metadata.json", "README.md", "ions.mdp",
                             "analysis_config.json", "plumed_template.dat"}

    removed: list[str] = []

    for f in step_dir.iterdir():
        if not f.is_file():
            continue
        if f.name in _PRESERVE_NAMES:
            continue
        if f.suffix in _PRESERVE_EXTENSIONS:
            continue

        # Eliminar outputs típicos de GROMACS
        if f.suffix in _REMOVE_EXTENSIONS:
            f.unlink()
            removed.append(f.name)
        elif _REMOVE_PATTERNS.match(f.name) and f.suffix == ".gro":
            f.unlink()
            removed.append(f.name)

    return removed


def _apply_action(
    action:   RemediationAction,
    step_dir: Path,
    log:      list[str],
) -> tuple[bool, list[str]]:
    """
    Aplica una RemediationAction al step_dir.

    Retorna (success, files_modified).
    """
    modified: list[str] = []

    # ── LOG_ONLY ──────────────────────────────────────────────────────────────
    if action.action_type == ActionType.LOG_ONLY:
        log.append(f"  [LOG]  {action.description}")
        return True, []

    # ── RESET_STEP ────────────────────────────────────────────────────────────
    if action.action_type == ActionType.RESET_STEP:
        removed = _reset_step_outputs(step_dir)
        log.append(f"  [RESET] Eliminados {len(removed)} archivos de output")
        for f in removed:
            log.append(f"    - {f}")
        return True, removed

    # ── DELETE_FILE ───────────────────────────────────────────────────────────
    if action.action_type == ActionType.DELETE_FILE:
        if not action.target_file:
            log.append("  [ERROR] DELETE_FILE sin target_file")
            return False, []
        target = step_dir / action.target_file
        if target.exists():
            backup = _backup(target)
            target.unlink()
            log.append(f"  [DELETE] {action.target_file} (backup: {backup.name})")
            modified.append(action.target_file)
            return True, modified
        log.append(f"  [SKIP] {action.target_file} no existe — nada que eliminar")
        return True, []

    # ── WRITE_FILE ────────────────────────────────────────────────────────────
    if action.action_type == ActionType.WRITE_FILE:
        if not action.target_file or action.file_content is None:
            log.append("  [ERROR] WRITE_FILE sin target_file o file_content")
            return False, []
        target = step_dir / action.target_file
        if target.exists():
            _backup(target)
        target.write_text(action.file_content)
        log.append(f"  [WRITE] {action.target_file} ({len(action.file_content)} chars)")
        modified.append(action.target_file)
        return True, modified

    # ── COPY_FILE ─────────────────────────────────────────────────────────────
    if action.action_type == ActionType.COPY_FILE:
        if not action.source_file or not action.target_file:
            log.append("  [ERROR] COPY_FILE sin source_file o target_file")
            return False, []
        src = Path(action.source_file)
        dst = step_dir / action.target_file
        if not src.exists():
            log.append(f"  [ERROR] Fuente no existe: {src}")
            return False, []
        if dst.exists():
            _backup(dst)
        shutil.copy2(src, dst)
        log.append(f"  [COPY]  {src.name} → {dst.name}")
        modified.append(action.target_file)
        return True, modified

    # ── PATCH_MDP / INJECT_RESTRAINTS / SCALE_TIMESTEP / REDUCE_TEMPERATURE ──
    if action.action_type in (
        ActionType.PATCH_MDP,
        ActionType.INJECT_RESTRAINTS,
        ActionType.SCALE_TIMESTEP,
        ActionType.REDUCE_TEMPERATURE,
    ):
        if not action.target_file or not action.patch_key or action.patch_value is None:
            log.append("  [ERROR] PATCH_MDP sin target_file, patch_key o patch_value")
            return False, []

        target = step_dir / action.target_file
        if not target.exists():
            # Buscar cualquier .mdp en el directorio
            mdp_candidates = list(step_dir.glob("*.mdp"))
            if not mdp_candidates:
                log.append(f"  [ERROR] No se encontró {action.target_file} ni ningún .mdp")
                return False, []
            target = mdp_candidates[0]
            log.append(f"  [INFO] Usando {target.name} en lugar de {action.target_file}")

        _backup(target)
        success = _patch_text_file(target, action.patch_key, action.patch_value)
        if success:
            old = f" (era: {action.patch_old})" if action.patch_old else ""
            log.append(
                f"  [PATCH] {target.name}: {action.patch_key} = {action.patch_value}{old}"
            )
            modified.append(target.name)
        else:
            log.append(f"  [WARN]  Patch sin efecto en {target.name} — sin cambios")
        return True, modified

    # ── PATCH_SCRIPT ──────────────────────────────────────────────────────────
    if action.action_type == ActionType.PATCH_SCRIPT:
        if not action.target_file or not action.patch_key or action.patch_value is None:
            log.append("  [ERROR] PATCH_SCRIPT incompleto")
            return False, []
        target = step_dir / action.target_file
        if not target.exists():
            log.append(f"  [ERROR] Script no encontrado: {action.target_file}")
            return False, []
        _backup(target)
        success = _patch_text_file(target, action.patch_key, action.patch_value)
        if success:
            log.append(f"  [PATCH_SH] {action.target_file}: {action.patch_key} = {action.patch_value}")
            modified.append(action.target_file)
        return True, modified

    log.append(f"  [WARN] ActionType no manejado: {action.action_type}")
    return False, []


# ═══════════════════════════════════════════════════════════════════════════════
# RemediationExecutor
# ═══════════════════════════════════════════════════════════════════════════════

class RemediationExecutor:
    """
    Ejecutor del loop adaptativo execute → diagnose → remediate → retry.

    Uso básico:
        # Después de que ShellExecutor produce un estado con fallos:
        rem = RemediationExecutor(workspace_path, dry_run=True)
        final_state = rem.run_adaptive(execution_state)

    Uso avanzado (control manual):
        plan = rem.diagnose_and_plan(record)
        rem.apply(plan, record)
        rem.retry(record, execution_state)
    """

    def __init__(
        self,
        workspace_path: str | Path,
        dry_run:        bool = True,
        max_global_retries: int = 3,
    ):
        self.workspace_path      = Path(workspace_path)
        self.dry_run             = dry_run
        self.max_global_retries  = max_global_retries
        self.reasoner            = AdaptiveReasoner()
        self._log_lines: list[str] = []

    # ─── API pública ──────────────────────────────────────────────────────────

    def run_adaptive(
        self,
        state: WorkspaceExecutionState,
    ) -> WorkspaceExecutionState:
        """
        Loop adaptativo completo sobre un estado de ejecución.

        Para cada step fallido:
            1. Diagnosticar
            2. Planificar remediación
            3. Aplicar cambios
            4. Retry
            5. Repetir hasta max_retries o éxito

        Después de remediar todos los steps posibles, verifica si
        los steps bloqueados pueden desbloquearse y los relanza.
        """
        self._log("=" * 60)
        self._log("RemediationExecutor — Loop adaptativo iniciado")
        self._log(f"Workspace: {self.workspace_path}")
        self._log(f"dry_run  : {self.dry_run}")
        self._log("=" * 60)

        failed_records = [
            r for r in state.steps
            if r.status == StepStatus.FAILED
        ]

        if not failed_records:
            self._log("No hay steps fallidos. Nada que remediar.")
            return state

        self._log(f"Steps fallidos detectados: {len(failed_records)}")

        for record in failed_records:
            self._remediate_record(record, state)

        # ── Desbloquear steps que dependían de los ahora corregidos ──────────
        self._unblock_and_retry_dependents(state)

        # ── Persistir estado final ─────────────────────────────────────────────
        self._save_state(state)
        self._log("=" * 60)
        self._log(
            f"Loop adaptativo completado — "
            f"done={state.n_done()} failed={state.n_failed()} "
            f"pending={state.n_pending()}"
        )

        return state

    def diagnose_and_plan(
        self,
        record: StepExecutionRecord,
    ) -> RemediationPlan:
        """
        Diagnostica un step fallido y genera el plan de remediación.
        No aplica ningún cambio al disco.
        """
        step_dir = Path(record.step_dir)
        diag     = self.reasoner.diagnose(record, step_dir)
        plan     = self.reasoner.plan_remediation(diag, step_dir)
        return plan

    def apply(
        self,
        plan:   RemediationPlan,
        record: StepExecutionRecord,
    ) -> RemediationRecord:
        """
        Aplica las acciones del plan al disco.

        Retorna un RemediationRecord con el resultado de la aplicación.
        """
        step_dir    = Path(record.step_dir)
        attempt_num = record.n_remediations() + 1

        rem_record = RemediationRecord(
            attempt_number = attempt_num,
            plan           = plan,
            status         = RemediationStatus.PENDING,
            started_at     = datetime.now(),
        )

        if not plan.is_applicable:
            rem_record.status = RemediationStatus.SKIPPED
            rem_record.executor_notes.append(
                "Plan marcado como no aplicable (error fatal). "
                "Requiere intervención humana."
            )
            self._log(f"  [SKIP] {record.step_id}: plan no aplicable — {plan.strategy}")
            record.add_remediation(rem_record)
            return rem_record

        if self.dry_run:
            return self._dry_run_apply(plan, record, rem_record, step_dir)

        # ── Aplicar acciones ───────────────────────────────────────────────────
        all_modified: list[str] = []
        action_log:   list[str] = []
        all_success = True

        for i, action in enumerate(plan.actions, start=1):
            action_log.append(
                f"  Acción {i}/{len(plan.actions)}: [{action.action_type.value}] "
                f"{action.description}"
            )
            success, modified = _apply_action(action, step_dir, action_log)
            all_modified.extend(modified)
            if not success:
                all_success = False
                action_log.append(f"    ✖ Acción {i} falló")

        rem_record.files_modified  = all_modified
        rem_record.executor_notes  = action_log
        rem_record.status          = RemediationStatus.APPLIED if all_success else RemediationStatus.FAILED
        rem_record.finished_at     = datetime.now()
        rem_record.elapsed_s       = (
            rem_record.finished_at - rem_record.started_at
        ).total_seconds()

        for line in action_log:
            self._log(line)

        record.add_remediation(rem_record)
        return rem_record

    def retry(
        self,
        record: StepExecutionRecord,
        state:  WorkspaceExecutionState,
    ) -> bool:
        """
        Relanza un step individual via ShellExecutor.

        Actualiza record.status y record.retry_count.
        Retorna True si el retry tuvo éxito.
        """
        from executors.shell_executor import ShellExecutor

        self._log(f"  [RETRY] {record.step_id} (intento #{record.retry_count + 1})")

        # Reset del estado del record para el retry
        record.status      = StepStatus.RUNNING
        record.started_at  = datetime.now()
        record.stdout      = ""
        record.stderr      = ""
        record.exit_code   = None
        record.error_message = None
        record.outputs_found   = []
        record.outputs_missing = []

        if self.dry_run:
            # Simular retry exitoso en dry-run
            record.status      = StepStatus.DONE
            record.finished_at = datetime.now()
            record.elapsed_s   = 0.1
            record.retry_count += 1
            self._log(f"  [DRY]   retry simulado exitoso para {record.step_id}")
            return True

        # ── Retry real vía ShellExecutor ─────────────────────────────────────
        # Creamos un executor temporal con el workspace_path y
        # lo usamos solo para ejecutar este step puntual.
        executor = ShellExecutor(
            workspace_path = self.workspace_path,
            dry_run        = False,
        )
        executor.state = state

        try:
            executor._run_step(record)
        except Exception as e:
            record.status        = StepStatus.FAILED
            record.error_message = str(e)

        record.finished_at = datetime.now()
        record.elapsed_s   = (record.finished_at - record.started_at).total_seconds()
        record.retry_count += 1

        success = record.status == StepStatus.DONE
        icon    = "✓" if success else "✖"
        self._log(
            f"  [{icon}] Retry {record.step_id}: {record.status.value} "
            f"({record.elapsed_s:.1f}s)"
        )
        return success

    # ─── Internos ─────────────────────────────────────────────────────────────

    def _remediate_record(
        self,
        record: StepExecutionRecord,
        state:  WorkspaceExecutionState,
    ) -> None:
        """
        Ciclo completo diagnose → plan → apply → retry para un record.
        Respeta max_retries del plan y max_global_retries.
        """
        self._log(f"\n{'─'*50}")
        self._log(f"Remediando: {record.step_id}")

        attempt = 0

        while record.status == StepStatus.FAILED:
            attempt += 1

            if attempt > self.max_global_retries:
                self._log(
                    f"  [ABORT] {record.step_id}: "
                    f"máximo de intentos globales alcanzado ({self.max_global_retries})"
                )
                break

            # ── Diagnosticar ──────────────────────────────────────────────────
            plan = self.diagnose_and_plan(record)
            self._log(
                f"  [DIAG]  {record.step_id}: "
                f"{plan.diagnosis.category.value} "
                f"(conf={plan.diagnosis.confidence:.0%})"
            )
            self._log(f"  [PLAN]  {plan.strategy}")

            # ── Verificar si ya se superó max_retries del plan ────────────────
            if record.n_remediations() >= plan.max_retries:
                self._log(
                    f"  [ABORT] {record.step_id}: "
                    f"max_retries del plan alcanzado ({plan.max_retries})"
                )
                break

            if not plan.is_applicable:
                self._log(
                    f"  [FATAL] {record.step_id}: "
                    "error no remediable automáticamente"
                )
                self._log(f"  → {plan.expected_outcome}")
                break

            # ── Aplicar remediación ───────────────────────────────────────────
            rem_record = self.apply(plan, record)

            # Persistir inmediatamente: si retry falla o se interrumpe,
            # la remediación aplicada queda registrada en execution_state.json.
            self._save_state(state)

            if rem_record.status == RemediationStatus.SKIPPED:
                break

            # ── Retry ─────────────────────────────────────────────────────────
            success = self.retry(record, state)

            # ── Actualizar RemediationRecord con resultado del retry ───────────
            # Reconstruir desde record.remediations para actualizar el dict
            # persistido, ya que rem_record es el objeto Python en memoria.
            rem_record.status         = (
                RemediationStatus.SUCCEEDED if success
                else RemediationStatus.FAILED
            )
            rem_record.retry_exit_code = record.exit_code
            rem_record.retry_stdout    = record.stdout[:2000]
            rem_record.retry_stderr    = record.stderr[:2000]

            # Reemplazar el último dict serializado con el estado actualizado
            if record.remediations:
                record.remediations[-1] = rem_record.model_dump()

            self._save_state(state)

            if success:
                self._log(
                    f"  [OK]    {record.step_id}: remediado exitosamente "
                    f"en intento {attempt}"
                )
                break

    def _unblock_and_retry_dependents(
        self,
        state: WorkspaceExecutionState,
    ) -> None:
        """
        Después de aplicar remediaciones, algunos steps bloqueados
        pueden haberse desbloqueado si sus dependencias ahora están DONE.

        Lee depends_on desde el metadata.json de cada step bloqueado.
        Si todos los predecesores están DONE → retry el step.
        """
        blocked = [r for r in state.steps if r.status == StepStatus.BLOCKED]
        if not blocked:
            return

        step_status = {r.step_id: r.status for r in state.steps}

        for record in blocked:
            step_dir  = Path(record.step_dir)
            meta_file = step_dir / "metadata.json"

            depends_on: list[str] = []
            if meta_file.exists():
                try:
                    meta       = json.loads(meta_file.read_text())
                    depends_on = meta.get("depends_on", [])
                except Exception:
                    pass

            if not depends_on:
                # Sin información de dependencias → intentar desbloquear
                can_unblock = True
            else:
                can_unblock = all(
                    step_status.get(dep) == StepStatus.DONE
                    for dep in depends_on
                )

            if can_unblock:
                self._log(f"\n  [UNBLOCK] {record.step_id}: predecesores OK → retry")
                record.status = StepStatus.FAILED   # forzar para entrar en retry
                self._remediate_record(record, state)

    def _dry_run_apply(
        self,
        plan:       RemediationPlan,
        record:     StepExecutionRecord,
        rem_record: RemediationRecord,
        step_dir:   Path,
    ) -> RemediationRecord:
        """Simula la aplicación del plan sin tocar el disco."""
        self._log(f"  [DRY] Remediación de {record.step_id}:")
        self._log(f"        Estrategia: {plan.strategy}")
        self._log(f"        Categoría : {plan.diagnosis.category.value}")
        self._log(f"        Severidad : {plan.diagnosis.severity.value}")

        simulated_files: list[str] = []
        for i, action in enumerate(plan.actions, start=1):
            self._log(
                f"    {i}. [{action.action_type.value}] {action.description}"
            )
            if action.target_file:
                self._log(f"       → {action.target_file}")
            if action.patch_key:
                old = f" (era: {action.patch_old})" if action.patch_old else ""
                self._log(
                    f"       → {action.patch_key} = {action.patch_value}{old}"
                )
            if action.target_file and action.action_type != ActionType.LOG_ONLY:
                simulated_files.append(action.target_file)

        rem_record.status         = RemediationStatus.APPLIED
        rem_record.files_modified = simulated_files
        rem_record.executor_notes = [
            "[dry-run] Acciones simuladas — ningún archivo fue modificado"
        ]
        rem_record.finished_at = datetime.now()
        rem_record.elapsed_s   = 0.0

        record.add_remediation(rem_record)
        return rem_record

    def _save_state(self, state: WorkspaceExecutionState) -> None:
        state_file = self.workspace_path / "execution_state.json"
        try:
            state_file.write_text(state.model_dump_json(indent=2))
        except Exception as e:
            self._log(f"  [WARN] No se pudo persistir estado: {e}")

    def _log(self, message: str) -> None:
        ts   = datetime.now().strftime("%H:%M:%S")
        line = f"[{ts}] {message}"
        print(line)
        self._log_lines.append(line)

    # ─── Acceso al historial ──────────────────────────────────────────────────

    def remediation_history(self, record: StepExecutionRecord) -> list[RemediationRecord]:
        """Retorna el historial de remediaciones de un step reconstruido desde el record."""
        return [RemediationRecord(**r) for r in record.remediations]

    def full_report(self, state: WorkspaceExecutionState) -> dict:
        """
        Resumen del loop adaptativo completo para logging y auditoría.
        Lee desde state.steps — la fuente de verdad persistida.
        """
        all_records = [
            RemediationRecord(**r)
            for step in state.steps
            for r in step.remediations
        ]

        by_category: dict[str, int] = {}
        for r in all_records:
            cat = r.plan.diagnosis.category.value
            by_category[cat] = by_category.get(cat, 0) + 1

        steps_with_remediations = sum(
            1 for step in state.steps if step.remediations
        )

        return {
            "total_steps_remediated": steps_with_remediations,
            "total_attempts":         len(all_records),
            "total_succeeded":        sum(1 for r in all_records if r.status == RemediationStatus.SUCCEEDED),
            "total_failed":           sum(1 for r in all_records if r.status == RemediationStatus.FAILED),
            "total_skipped":          sum(1 for r in all_records if r.status == RemediationStatus.SKIPPED),
            "by_error_category":      by_category,
            "log_lines":              len(self._log_lines),
        }