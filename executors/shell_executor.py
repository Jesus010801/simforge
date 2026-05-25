# executors/shell_executor.py
"""
Shell Executor — ejecutor concreto de SimForge.

Dos modos controlados por dry_run:

    dry_run=True:
        Simula la ejecución sin lanzar ningún proceso.
        Verifica que los scripts existen, lee metadata,
        imprime lo que haría. Útil para validar un workspace
        antes de enviarlo a un cluster.

    dry_run=False (default para `simforge run`):
        Ejecuta los scripts reales via subprocess con streaming.
        Muestra stdout/stderr en tiempo real con prefijo [gmx].
        Detecta hangs y procesos interactivos.
        Verifica outputs esperados post-ejecución.

El ShellExecutor no sabe nada de GROMACS ni de química.
Solo sabe ejecutar scripts de shell y verificar archivos.
"""

from __future__ import annotations

from pathlib import Path
import json

from rich.console import Console

from executors.base_executor import BaseExecutor
from executors.execution_state import StepExecutionRecord, StepStatus
from executors.stream import run_streaming

_console = Console(highlight=False)

# ═══════════════════════════════════════════════════════════════════════════════
# Shell Executor
# ═══════════════════════════════════════════════════════════════════════════════

class ShellExecutor(BaseExecutor):
    """
    Ejecutor de scripts de shell con soporte dry-run y streaming real.
    """

    _SCRIPT_CANDIDATES: list[str] = [
        "run.sh",
        "run_md.sh",
        "run_nvt.sh",
        "run_analysis.sh",
        "commands.sh",
    ]

    def _run_step(self, record: StepExecutionRecord) -> None:
        step_dir  = Path(record.step_dir)
        meta      = self._read_metadata(step_dir)
        step_type = meta.get("step_type", "automatic")

        automation_level = meta.get("automation_level")
        if automation_level is not None:
            needs_user = automation_level in ("manual", "guided")
        else:
            needs_user = step_type in ("manual", "external", "validation")

        if needs_user:
            label = automation_level or step_type
            record.status        = StepStatus.SKIPPED
            record.error_message = f"Step '{label}' — requiere acción del usuario"
            self._log(f"  [manual] {step_dir.name} → ver README.md")
            return

        script = self._find_script(step_dir)
        if script is None:
            record.status        = StepStatus.SKIPPED
            record.error_message = "No se encontró script ejecutable en el directorio"
            self._log(f"  [skip]   sin script en {step_dir.name}")
            return

        expected_outputs: list[str] = meta.get("expected_outputs", [])

        if self.dry_run:
            self._dry_run_step(record, step_dir, script, expected_outputs, meta)
        else:
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
            self._log(f"  [dry]  outputs  : {', '.join(expected_outputs)}")

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
        record.outputs_found = expected_outputs

    # ── Ejecución real con streaming ──────────────────────────────────────────

    def _real_run_step(
        self,
        record:           StepExecutionRecord,
        step_dir:         Path,
        script:           Path,
        expected_outputs: list[str],
    ) -> None:
        step_name = step_dir.name

        def _on_stdout(line: str) -> None:
            if line.strip():
                _console.print(f"  [dim][gmx][/dim] {line}")

        def _on_stderr(line: str) -> None:
            # GROMACS writes most output to stderr; display but de-emphasize
            if line.strip():
                _console.print(f"  [dim][gmx][/dim] {line}")

        def _on_heartbeat(elapsed_s: float) -> None:
            mins = int(elapsed_s // 60)
            secs = int(elapsed_s % 60)
            elapsed_str = f"{mins}m {secs}s" if mins else f"{secs}s"
            _console.print(
                f"  [yellow]⏳[/yellow] [dim]{step_name} still running "
                f"(no output for {elapsed_str})...[/dim]"
            )

        try:
            result = run_streaming(
                cmd=["bash", str(script)],
                cwd=step_dir,
                on_stdout=_on_stdout,
                on_stderr=_on_stderr,
                on_heartbeat=_on_heartbeat,
                timeout_s=86400,
                heartbeat_interval_s=30,
                hang_timeout_s=600,
            )
        except Exception as exc:
            record.status        = StepStatus.FAILED
            record.error_message = str(exc)
            return

        # Merge stderr into stdout for single-field storage
        # (GROMACS logs extensively to stderr)
        combined = "\n".join(filter(None, [result.stdout, result.stderr]))
        record.stdout    = combined
        record.stderr    = result.stderr
        record.exit_code = result.returncode

        if result.timed_out:
            record.status        = StepStatus.FAILED
            record.error_message = "Timeout: el step excedió 24 horas"
            return

        if result.hang_killed:
            _console.print(
                f"\n  [red]✗[/red]  [bold]{step_name}[/bold] fue terminado por hang: "
                "no emitió output por 10 minutos.\n"
                "  [dim]Causa probable: comando GROMACS esperando input interactivo "
                "(ej: gmx editconf -princ sin pipe de selección).[/dim]\n"
            )
            record.status        = StepStatus.FAILED
            record.error_message = (
                "Proceso terminado por hang (10 min sin output). "
                "Probablemente esperando input interactivo del terminal."
            )
            return

        if result.likely_interactive:
            _console.print(
                f"\n  [yellow]⚠[/yellow]  [bold]{step_name}[/bold] appears to be waiting "
                "for interactive input.\n"
                "  [dim]GROMACS commands should not require TTY. "
                "Check the script for missing -noconfirm / -yes flags.[/dim]\n"
            )

        if result.returncode != 0:
            record.status        = StepStatus.FAILED
            record.error_message = (
                f"Exit code {result.returncode}. "
                f"Last stderr: {result.stderr[-500:] if result.stderr else '(none)'}"
            )
            return

        found, missing = self._check_expected_outputs(step_dir, expected_outputs)
        record.outputs_found   = found
        record.outputs_missing = missing

        if missing:
            record.status        = StepStatus.FAILED
            record.error_message = f"Outputs esperados no encontrados: {missing}"
        else:
            record.status = StepStatus.DONE

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _find_script(self, step_dir: Path) -> Path | None:
        for name in self._SCRIPT_CANDIDATES:
            candidate = step_dir / name
            if candidate.exists():
                return candidate
        return None
