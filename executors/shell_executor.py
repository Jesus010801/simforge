# executors/shell_executor.py
"""
Shell Executor — ejecutor concreto de SimForge.

Dos modos controlados por dry_run:

    dry_run=True  (default):
        Simula la ejecución sin lanzar ningún proceso.
        Verifica que los scripts existen, lee metadata,
        imprime lo que haría. Útil para validar un workspace
        antes de enviarlo a un cluster.

    dry_run=False:
        Ejecuta los scripts reales via subprocess.
        Captura stdout/stderr en tiempo real.
        Verifica outputs esperados post-ejecución.

El ShellExecutor no sabe nada de GROMACS ni de química.
Solo sabe ejecutar scripts de shell y verificar archivos.
"""

from __future__ import annotations

from pathlib import Path
import subprocess
import json

from executors.base_executor import BaseExecutor
from executors.execution_state import (
    StepExecutionRecord,
    StepStatus,
)


# ═══════════════════════════════════════════════════════════════════════════════
# Shell Executor
# ═══════════════════════════════════════════════════════════════════════════════

class ShellExecutor(BaseExecutor):
    """
    Ejecutor de scripts de shell con soporte dry-run.
    """

    # Nombre del script principal por stage
    # El executor busca estos en orden — usa el primero que encuentre
    _SCRIPT_CANDIDATES: list[str] = [
        "run.sh",
        "run_md.sh",
        "run_nvt.sh",
        "run_analysis.sh",
        "commands.sh",
    ]

    def _run_step(self, record: StepExecutionRecord) -> None:

        step_dir = Path(record.step_dir)
        meta     = self._read_metadata(step_dir)

        # ── Detectar script ──────────────────────────────────────────────────
        script = self._find_script(step_dir)

        if script is None:
            # Verificar si es un step manual o externo — skip esperado
            step_type = meta.get("step_type", "automatic")
            if step_type in ("manual", "external", "validation"):
                record.status        = StepStatus.SKIPPED
                record.error_message = f"Step tipo '{step_type}' — requiere acción manual"
                self._log(f"  [manual] {step_dir.name} → ver README.md")
            else:
                record.status        = StepStatus.SKIPPED
                record.error_message = "No se encontró script ejecutable en el directorio"
                self._log(f"  [skip]   sin script en {step_dir.name}")
            return

        # ── Leer outputs esperados ────────────────────────────────────────────
        expected_outputs: list[str] = meta.get("expected_outputs", [])

        # ── Dry run ──────────────────────────────────────────────────────────
        if self.dry_run:
            self._dry_run_step(record, step_dir, script, expected_outputs, meta)
            return

        # ── Ejecución real ───────────────────────────────────────────────────
        self._real_run_step(record, step_dir, script, expected_outputs)

    # ── Dry run ───────────────────────────────────────────────────────────────

    def _dry_run_step(
        self,
        record:           StepExecutionRecord,
        step_dir:         Path,
        script:           Path,
        expected_outputs: list[str],
        meta:             dict,
    ) -> None:

        self._log(f"  [dry]  script   : {script.name}")
        self._log(f"  [dry]  engine   : {meta.get('engine', '?')}")
        self._log(f"  [dry]  stage    : {meta.get('stage', '?')}")

        if expected_outputs:
            self._log(
                f"  [dry]  outputs  : {', '.join(expected_outputs)}"
            )

        # Mostrar primeras líneas del script
        try:
            lines = script.read_text().splitlines()
            preview = [l for l in lines if l.strip() and not l.startswith("#")][:4]
            for line in preview:
                self._log(f"  [dry]  > {line.strip()}")
        except Exception:
            pass

        record.status        = StepStatus.DONE
        record.stdout        = "[dry-run: no output]"
        record.exit_code     = 0
        record.outputs_found = expected_outputs   # en dry-run asumimos OK

    # ── Ejecución real ────────────────────────────────────────────────────────

    def _real_run_step(
        self,
        record:           StepExecutionRecord,
        step_dir:         Path,
        script:           Path,
        expected_outputs: list[str],
    ) -> None:

        try:
            result = subprocess.run(
                ["bash", str(script)],
                cwd     = str(step_dir),
                capture_output = True,
                text    = True,
                timeout = 86400,   # 24h máximo por step
            )

            record.stdout    = result.stdout
            record.stderr    = result.stderr
            record.exit_code = result.returncode

            if result.returncode != 0:
                record.status        = StepStatus.FAILED
                record.error_message = (
                    f"Exit code {result.returncode}. "
                    f"Stderr: {result.stderr[:500]}"
                )
                return

            # Verificar outputs esperados
            found, missing = self._check_expected_outputs(
                step_dir, expected_outputs
            )
            record.outputs_found   = found
            record.outputs_missing = missing

            if missing:
                record.status        = StepStatus.FAILED
                record.error_message = (
                    f"Outputs esperados no encontrados: {missing}"
                )
            else:
                record.status = StepStatus.DONE

        except subprocess.TimeoutExpired:
            record.status        = StepStatus.FAILED
            record.error_message = "Timeout: el step excedió 24 horas"

        except Exception as e:
            record.status        = StepStatus.FAILED
            record.error_message = str(e)

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _find_script(self, step_dir: Path) -> Path | None:
        for name in self._SCRIPT_CANDIDATES:
            candidate = step_dir / name
            if candidate.exists():
                return candidate
        return None
