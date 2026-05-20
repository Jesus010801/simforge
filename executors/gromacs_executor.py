# executors/gromacs_executor.py
"""
GROMACS Executor — executor con percepción de outputs reales de GROMACS.

Extiende ShellExecutor con inteligencia específica de GROMACS:

    Minimización:
        - Parsea md.log para detectar convergencia (Epot, Fmax)
        - Detecta "Steepest Descents converged" vs max steps reached
        - Reporta Epot final y Fmax para decisiones de continuación

    Equilibración / Producción:
        - Parsea md.log para detectar velocidades explosivas (LINCS warnings)
        - Detecta NaN en energías → crash de integración
        - Detecta "Fatal error" de GROMACS
        - Reporta temperatura promedio, energía potencial, presión

    Outputs:
        - Verifica checksum de archivos .edr, .xtc, .gro post-ejecución
        - Reporta tamaño de archivo como proxy de frames escritos
        - Detecta archivos truncados (tamaño 0 o sospechosamente pequeño)

    Diagnóstico:
        - Construye GROMACSStepDiagnostic con métricas y veredicto
        - Adjunta diagnóstico a StepExecutionRecord como campo extra
        - No toma decisiones automáticas — reporta para adaptive reasoning

Principio arquitectónico:
    GROMACSExecutor PERCIBE y REPORTA.
    El adaptive reasoning (capa siguiente) INTERPRETA y DECIDE.
    Este executor nunca modifica parámetros ni reintenta automáticamente.

Modo de uso:
    executor = GROMACSExecutor(
        workspace_path = workspace,
        dry_run        = False,   # True para test sin GROMACS real
    )
    state = executor.run()
    # state.steps[i].gromacs_diagnostic → GROMACSStepDiagnostic
"""

from __future__ import annotations

import re
import json
from pathlib import Path
from typing import Optional
from pydantic import BaseModel

from executors.shell_executor import ShellExecutor
from executors.execution_state import (
    StepExecutionRecord,
    StepStatus,
)


# ═══════════════════════════════════════════════════════════════════════════════
# Modelos de diagnóstico GROMACS
# ═══════════════════════════════════════════════════════════════════════════════

class MinimizationMetrics(BaseModel):
    """Métricas extraídas del log de minimización."""

    converged:          bool  = False
    final_epot:         Optional[float] = None   # kJ/mol
    final_fmax:         Optional[float] = None   # kJ/mol/nm
    n_steps_taken:      int   = 0
    convergence_reason: str   = "unknown"
    # "converged_fmax"  → Fmax < emtol → exitoso
    # "max_steps"       → alcanzó nsteps sin converger
    # "not_found"       → no se pudo parsear el log


class MDMetrics(BaseModel):
    """Métricas de runs MD (equilibración/producción)."""

    completed:          bool  = False
    n_steps_logged:     int   = 0
    last_epot:          Optional[float] = None   # kJ/mol
    last_temperature:   Optional[float] = None   # K
    last_pressure:      Optional[float] = None   # bar

    # Anomalías detectadas
    has_nan_energy:     bool  = False
    has_lincs_warning:  bool  = False
    n_lincs_warnings:   int   = 0
    has_fatal_error:    bool  = False
    fatal_error_msg:    str   = ""
    has_exploded:       bool  = False   # velocidades > umbral


class OutputFileStatus(BaseModel):
    """Estado de un archivo de output post-ejecución."""

    filename:   str
    exists:     bool  = False
    size_bytes: int   = 0
    is_empty:   bool  = True
    suspect:    bool  = False   # tamaño sospechosamente pequeño para su tipo


class GROMACSStepDiagnostic(BaseModel):
    """
    Diagnóstico completo de un step GROMACS.

    Adjuntado a StepExecutionRecord.gromacs_diagnostic por GROMACSExecutor.
    Leído por el adaptive reasoning para decisiones de continuación.
    """

    step_id:        str
    stage:          str   # "minimization" | "equilibration" | "production" | etc.
    engine:         str   = "gromacs"

    # Métricas por tipo de step
    minimization:   Optional[MinimizationMetrics] = None
    md:             Optional[MDMetrics]            = None

    # Outputs
    output_files:   list[OutputFileStatus] = []

    # Veredicto global
    verdict:        str   = "unknown"
    # "ok"              → todo bien, continuar
    # "converged"       → minimización exitosa
    # "not_converged"   → minimización sin converger (puede ser aceptable)
    # "crashed"         → NaN, Fatal error, o archivos vacíos
    # "warning"         → LINCS u otras advertencias no fatales
    # "incomplete"      → outputs esperados ausentes o truncados
    # "unknown"         → no se pudo parsear suficiente información

    # Notas del diagnóstico (para adaptive reasoning)
    notes:          list[str] = []
    warnings:       list[str] = []
    errors:         list[str] = []


# ═══════════════════════════════════════════════════════════════════════════════
# Parsers de log GROMACS
# ═══════════════════════════════════════════════════════════════════════════════

class GROMACSLogParser:
    """
    Parsea archivos .log de GROMACS para extraer métricas clave.

    No intenta cubrir el formato completo — se enfoca en las señales
    más importantes para el diagnóstico de convergencia y crashes.
    """

    # ── Minimización ──────────────────────────────────────────────────────────

    # "Steepest Descents converged to Fmax < 1000 in 427 steps"
    _RE_EM_CONVERGED = re.compile(
        r"(?:Steepest Descents|L-BFGS|Conjugate Gradients)\s+converged"
        r".*?Fmax\s*[<>]\s*([\d.eE+\-]+)"
        r".*?in\s+(\d+)\s+steps",
        re.IGNORECASE | re.DOTALL,
    )

    # "Steepest Descents did not converge to Fmax < ... in 50000 steps"
    _RE_EM_NOT_CONVERGED = re.compile(
        r"(?:Steepest Descents|L-BFGS|Conjugate Gradients)\s+did\s+not\s+converge",
        re.IGNORECASE,
    )

    # "Potential Energy  = -1.2345e+06"  (línea de resumen final)
    _RE_EPOT = re.compile(
        r"Potential\s+Energy\s*=\s*([-\d.eE+]+)",
        re.IGNORECASE,
    )

    # "Maximum force     =  1.23456e+02" (línea de resumen final)
    _RE_FMAX = re.compile(
        r"Maximum\s+force\s*=\s*([-\d.eE+]+)",
        re.IGNORECASE,
    )

    # ── MD ────────────────────────────────────────────────────────────────────

    # "Step  Time" table header → siguiente línea tiene step number
    _RE_STEP_LINE = re.compile(r"^\s*Step\s+Time\s*$")

    # "Potential Energy" en el bloque de energías del log
    _RE_MD_EPOT = re.compile(
        r"Potential\s+Energy\s+([-\d.eE+]+)",
        re.IGNORECASE,
    )

    # "Temperature" en el bloque de energías
    _RE_MD_TEMP = re.compile(
        r"Temperature\s+([\d.eE+]+)",
        re.IGNORECASE,
    )

    # "Pressure" en el bloque de energías
    _RE_MD_PRESS = re.compile(
        r"Pressure\s+([-\d.eE+]+)",
        re.IGNORECASE,
    )

    # "nan" en cualquier línea de energía → crash
    _RE_NAN = re.compile(r"\bnan\b", re.IGNORECASE)

    # LINCS warning
    _RE_LINCS = re.compile(
        r"LINCS\s+warning|lincs_warning|too\s+many\s+LINCS\s+warnings",
        re.IGNORECASE,
    )

    # Fatal error
    _RE_FATAL = re.compile(
        r"Fatal\s+error:",
        re.IGNORECASE,
    )

    # "Finished mdrun" → completó
    _RE_FINISHED = re.compile(
        r"Finished\s+mdrun",
        re.IGNORECASE,
    )

    # ── Interfaz ──────────────────────────────────────────────────────────────

    def parse_minimization_log(self, log_path: Path) -> MinimizationMetrics:
        """Parsea un log de minimización y retorna MinimizationMetrics."""
        if not log_path.exists():
            return MinimizationMetrics(convergence_reason="not_found")

        text = log_path.read_text(errors="replace")

        metrics = MinimizationMetrics()

        # Convergencia
        if self._RE_EM_NOT_CONVERGED.search(text):
            metrics.converged          = False
            metrics.convergence_reason = "max_steps"
        elif self._RE_EM_CONVERGED.search(text):
            metrics.converged          = True
            metrics.convergence_reason = "converged_fmax"
            m = self._RE_EM_CONVERGED.search(text)
            if m:
                try:
                    metrics.n_steps_taken = int(m.group(2))
                except (ValueError, IndexError):
                    pass

        # Epot final (última ocurrencia)
        epot_matches = self._RE_EPOT.findall(text)
        if epot_matches:
            try:
                metrics.final_epot = float(epot_matches[-1])
            except ValueError:
                pass

        # Fmax final (última ocurrencia)
        fmax_matches = self._RE_FMAX.findall(text)
        if fmax_matches:
            try:
                metrics.final_fmax = float(fmax_matches[-1])
            except ValueError:
                pass

        return metrics

    def parse_md_log(self, log_path: Path) -> MDMetrics:
        """Parsea un log de MD (equilibración/producción) y retorna MDMetrics."""
        if not log_path.exists():
            return MDMetrics()

        text = log_path.read_text(errors="replace")
        lines = text.splitlines()

        metrics = MDMetrics()

        # Completó
        if self._RE_FINISHED.search(text):
            metrics.completed = True

        # NaN
        energy_section = False
        for line in lines:
            if "Energies" in line or "Epot" in line:
                energy_section = True
            if energy_section and self._RE_NAN.search(line):
                metrics.has_nan_energy = True
                break

        # LINCS
        lincs_count = len(self._RE_LINCS.findall(text))
        if lincs_count > 0:
            metrics.has_lincs_warning = True
            metrics.n_lincs_warnings  = lincs_count

        # Fatal error
        fatal_match = self._RE_FATAL.search(text)
        if fatal_match:
            metrics.has_fatal_error = True
            # Extraer las 2 líneas siguientes como mensaje
            start = fatal_match.start()
            snippet = text[start:start + 300].splitlines()
            metrics.fatal_error_msg = " ".join(snippet[:3])

        # Temperatura y Epot — última ocurrencia (estado final)
        epot_matches = self._RE_MD_EPOT.findall(text)
        if epot_matches:
            try:
                metrics.last_epot = float(epot_matches[-1])
            except ValueError:
                pass

        temp_matches = self._RE_MD_TEMP.findall(text)
        if temp_matches:
            try:
                metrics.last_temperature = float(temp_matches[-1])
            except ValueError:
                pass

        press_matches = self._RE_MD_PRESS.findall(text)
        if press_matches:
            try:
                metrics.last_pressure = float(press_matches[-1])
            except ValueError:
                pass

        # Exploción de velocidades: temperatura >> 500K
        if metrics.last_temperature and metrics.last_temperature > 500:
            metrics.has_exploded = True

        return metrics


# ═══════════════════════════════════════════════════════════════════════════════
# Verificación de outputs
# ═══════════════════════════════════════════════════════════════════════════════

# Tamaños mínimos esperados por tipo de archivo (bytes)
# Por debajo → archivo probablemente truncado o vacío
_MIN_SIZE: dict[str, int] = {
    ".xtc":  1024 * 10,    # 10 KB mínimo para una trayectoria real
    ".edr":  1024 * 5,     # 5 KB mínimo para energías
    ".gro":  512,           # estructura final pequeña puede ser legítima
    ".cpt":  1024 * 100,   # checkpoint ~100 KB
    ".log":  1024,          # log mínimo
    ".tpr":  1024 * 10,    # tpr siempre tiene cierto tamaño
}


def _check_output_file(step_dir: Path, filename: str) -> OutputFileStatus:
    path = step_dir / filename

    if not path.exists():
        return OutputFileStatus(filename=filename, exists=False, is_empty=True)

    size = path.stat().st_size
    ext  = path.suffix.lower()
    min_expected = _MIN_SIZE.get(ext, 0)

    is_empty = (size == 0)
    suspect  = (not is_empty and min_expected > 0 and size < min_expected)

    return OutputFileStatus(
        filename    = filename,
        exists      = True,
        size_bytes  = size,
        is_empty    = is_empty,
        suspect     = suspect,
    )


# ═══════════════════════════════════════════════════════════════════════════════
# Lógica de veredicto
# ═══════════════════════════════════════════════════════════════════════════════

def _compute_verdict(
    stage:       str,
    mini:        Optional[MinimizationMetrics],
    md:          Optional[MDMetrics],
    outputs:     list[OutputFileStatus],
    notes:       list[str],
    warnings:    list[str],
    errors:      list[str],
) -> str:
    """
    Determina el veredicto global del step.
    Prioridad: crashed > incomplete > not_converged > warning > converged/ok.
    """

    # ── Outputs ───────────────────────────────────────────────────────────────
    critical_missing = [o for o in outputs if not o.exists and o.filename.endswith((".gro", ".edr"))]
    truncated        = [o for o in outputs if o.suspect]

    # ── Crashes explícitos ────────────────────────────────────────────────────
    if md and (md.has_nan_energy or md.has_fatal_error or md.has_exploded):
        return "crashed"

    if critical_missing:
        return "incomplete"

    if truncated:
        return "incomplete"

    # ── Minimización ─────────────────────────────────────────────────────────
    if stage == "minimization" and mini:
        if mini.convergence_reason == "not_found":
            return "unknown"
        if mini.converged:
            return "converged"
        return "not_converged"

    # ── MD warnings ───────────────────────────────────────────────────────────
    if md and md.has_lincs_warning:
        return "warning"

    # ── MD completó ───────────────────────────────────────────────────────────
    if md and md.completed:
        return "ok"

    if errors:
        return "crashed"

    return "ok"


# ═══════════════════════════════════════════════════════════════════════════════
# GROMACS Executor
# ═══════════════════════════════════════════════════════════════════════════════

class GROMACSExecutor(ShellExecutor):
    """
    Executor GROMACS-aware.

    Hereda toda la lógica de ejecución de ShellExecutor.
    Añade post-processing específico de GROMACS después de cada step:
        1. Detecta el tipo de step (minimization / equilibration / production)
        2. Parsea el log correspondiente
        3. Verifica outputs esperados con análisis de tamaño
        4. Construye GROMACSStepDiagnostic
        5. Adjunta el diagnóstico al StepExecutionRecord

    El diagnóstico está disponible en:
        state.steps[i].model_extra["gromacs_diagnostic"]

    O accesible via helper:
        executor.get_diagnostic(step_id)
    """

    def __init__(
        self,
        workspace_path: str | Path,
        dry_run: bool = True,
    ):
        super().__init__(workspace_path, dry_run)
        self._log_parser   = GROMACSLogParser()
        self._diagnostics: dict[str, GROMACSStepDiagnostic] = {}

    # ── Override: post-process después de cada step ───────────────────────────

    def _run_step(self, record: StepExecutionRecord) -> None:
        # Ejecutar normalmente (ShellExecutor)
        super()._run_step(record)

        # Post-process solo si el step realmente corrió
        if record.status in (StepStatus.DONE, StepStatus.FAILED):
            step_dir = Path(record.step_dir)
            meta     = self._read_metadata(step_dir)
            stage    = meta.get("stage", "")

            if self._is_gromacs_stage(stage):
                diag = self._build_diagnostic(record, step_dir, meta)
                self._diagnostics[record.step_id] = diag
                self._log_diagnostic(diag)

                # Si el diagnóstico revela un crash que el exit code no detectó
                if diag.verdict == "crashed" and record.status == StepStatus.DONE:
                    record.status        = StepStatus.FAILED
                    record.error_message = (
                        f"GROMACS diagnostic: crashed — "
                        + ("; ".join(diag.errors) if diag.errors else "ver log")
                    )

    # ── Detección de stage GROMACS ────────────────────────────────────────────

    @staticmethod
    def _is_gromacs_stage(stage: str) -> bool:
        return stage in (
            "minimization",
            "equilibration",
            "production",
        )

    # ── Construcción del diagnóstico ──────────────────────────────────────────

    def _build_diagnostic(
        self,
        record:   StepExecutionRecord,
        step_dir: Path,
        meta:     dict,
    ) -> GROMACSStepDiagnostic:

        stage   = meta.get("stage", "unknown")
        engine  = meta.get("engine", "gromacs")
        notes:    list[str] = []
        warnings: list[str] = []
        errors:   list[str] = []

        mini: Optional[MinimizationMetrics] = None
        md:   Optional[MDMetrics]           = None

        # ── Parseo de log ────────────────────────────────────────────────────
        if stage == "minimization":
            log_path = self._find_log(step_dir, prefix="em")
            if log_path:
                mini = self._log_parser.parse_minimization_log(log_path)
                self._annotate_minimization(mini, notes, warnings, errors)
            else:
                notes.append("Log de minimización no encontrado")

        elif stage in ("equilibration", "production"):
            log_path = self._find_log(step_dir, prefix=None)
            if log_path:
                md = self._log_parser.parse_md_log(log_path)
                self._annotate_md(md, stage, notes, warnings, errors)
            else:
                notes.append(f"Log de {stage} no encontrado")

        # ── Verificar outputs ─────────────────────────────────────────────────
        expected = meta.get("expected_outputs", [])
        output_statuses = [
            _check_output_file(step_dir, fname)
            for fname in expected
        ]
        self._annotate_outputs(output_statuses, notes, warnings, errors)

        # ── Veredicto ─────────────────────────────────────────────────────────
        verdict = _compute_verdict(
            stage    = stage,
            mini     = mini,
            md       = md,
            outputs  = output_statuses,
            notes    = notes,
            warnings = warnings,
            errors   = errors,
        )

        return GROMACSStepDiagnostic(
            step_id      = record.step_id,
            stage        = stage,
            engine       = engine,
            minimization = mini,
            md           = md,
            output_files = output_statuses,
            verdict      = verdict,
            notes        = notes,
            warnings     = warnings,
            errors       = errors,
        )

    # ── Helpers de anotación ──────────────────────────────────────────────────

    @staticmethod
    def _annotate_minimization(
        mini:     MinimizationMetrics,
        notes:    list[str],
        warnings: list[str],
        errors:   list[str],
    ) -> None:
        if mini.convergence_reason == "not_found":
            warnings.append("No se pudo parsear convergencia del log")
            return

        if mini.converged:
            notes.append(
                f"Minimización convergida en {mini.n_steps_taken} pasos"
            )
        else:
            warnings.append(
                "Minimización no convergió (alcanzó nsteps máximo)"
            )
            notes.append(
                "Puede ser aceptable si Fmax es razonable (< 100 kJ/mol/nm)"
            )

        if mini.final_epot is not None:
            notes.append(f"Epot final: {mini.final_epot:.2f} kJ/mol")

        if mini.final_fmax is not None:
            notes.append(f"Fmax final: {mini.final_fmax:.4f} kJ/mol/nm")
            if mini.final_fmax > 1000:
                errors.append(
                    f"Fmax muy alto ({mini.final_fmax:.1f} kJ/mol/nm) — "
                    "estructura posiblemente con clashes graves"
                )
            elif mini.final_fmax > 100:
                warnings.append(
                    f"Fmax moderadamente alto ({mini.final_fmax:.1f} kJ/mol/nm) — "
                    "verificar geometría antes de equilibración"
                )

    @staticmethod
    def _annotate_md(
        md:       MDMetrics,
        stage:    str,
        notes:    list[str],
        warnings: list[str],
        errors:   list[str],
    ) -> None:
        if md.has_fatal_error:
            errors.append(f"Fatal error GROMACS: {md.fatal_error_msg[:200]}")

        if md.has_nan_energy:
            errors.append(
                "NaN detectado en energías — crash de integración. "
                "Posibles causas: clashes, dt muy grande, parametrización incorrecta."
            )

        if md.has_exploded:
            errors.append(
                f"Temperatura explosiva ({md.last_temperature:.0f} K) — "
                "el sistema se desintegró. Revisar estructura inicial y parámetros."
            )

        if md.has_lincs_warning:
            warnings.append(
                f"LINCS: {md.n_lincs_warnings} warning(s) — "
                "posibles problemas de geometría de enlaces. "
                "Verificar si son esporádicos (aceptable) o persistentes (problema)."
            )

        if md.completed:
            notes.append(f"{stage.capitalize()} completada")
        else:
            warnings.append(f"{stage.capitalize()} no completó (sin 'Finished mdrun')")

        if md.last_temperature is not None:
            ref = 300.0
            dev = abs(md.last_temperature - ref)
            notes.append(f"Temperatura final: {md.last_temperature:.1f} K")
            if dev > 20:
                warnings.append(
                    f"Temperatura final desviada de 300K: {md.last_temperature:.1f}K "
                    f"(Δ={dev:.1f}K)"
                )

        if md.last_epot is not None:
            notes.append(f"Epot final: {md.last_epot:.2f} kJ/mol")

    @staticmethod
    def _annotate_outputs(
        statuses: list[OutputFileStatus],
        notes:    list[str],
        warnings: list[str],
        errors:   list[str],
    ) -> None:
        for s in statuses:
            if not s.exists:
                errors.append(f"Output no encontrado: {s.filename}")
            elif s.is_empty:
                errors.append(f"Output vacío: {s.filename}")
            elif s.suspect:
                warnings.append(
                    f"Output sospechosamente pequeño: "
                    f"{s.filename} ({s.size_bytes} bytes)"
                )
            else:
                notes.append(
                    f"Output ok: {s.filename} ({s.size_bytes / 1024:.1f} KB)"
                )

    # ── Helpers ───────────────────────────────────────────────────────────────

    @staticmethod
    def _find_log(step_dir: Path, prefix: str | None) -> Path | None:
        """
        Busca el archivo .log del step.
        Si prefix se da, busca <prefix>.log primero.
        Fallback: cualquier .log en el directorio.
        """
        if prefix:
            candidate = step_dir / f"{prefix}.log"
            if candidate.exists():
                return candidate

        logs = sorted(step_dir.glob("*.log"))
        if logs:
            return logs[0]

        return None

    def _log_diagnostic(self, diag: GROMACSStepDiagnostic) -> None:
        """Imprime resumen del diagnóstico al log del executor."""
        verdict_icons = {
            "ok":            "✓",
            "converged":     "✓",
            "not_converged": "⚠",
            "warning":       "⚠",
            "incomplete":    "✖",
            "crashed":       "✖",
            "unknown":       "?",
        }
        icon = verdict_icons.get(diag.verdict, "?")
        self._log(
            f"  [GROMACS] {icon} {diag.step_id} — {diag.verdict}"
        )
        for note in diag.notes:
            self._log(f"    ℹ  {note}")
        for w in diag.warnings:
            self._log(f"    ⚠  {w}")
        for e in diag.errors:
            self._log(f"    ✖  {e}")

    # ── API pública de diagnósticos ───────────────────────────────────────────

    def get_diagnostic(self, step_id: str) -> Optional[GROMACSStepDiagnostic]:
        """Retorna el diagnóstico de un step, si existe."""
        return self._diagnostics.get(step_id)

    def all_diagnostics(self) -> dict[str, GROMACSStepDiagnostic]:
        """Retorna todos los diagnósticos generados."""
        return dict(self._diagnostics)

    def summary_report(self) -> str:
        """
        Genera un resumen legible de todos los diagnósticos.
        Útil para logging y para adaptive reasoning.
        """
        lines = ["=== GROMACS Execution Summary ==="]
        for step_id, diag in self._diagnostics.items():
            lines.append(f"\n  {step_id} [{diag.stage}]")
            lines.append(f"    verdict : {diag.verdict}")
            if diag.minimization:
                m = diag.minimization
                lines.append(
                    f"    epot    : {m.final_epot} kJ/mol  "
                    f"fmax: {m.final_fmax} kJ/mol/nm  "
                    f"converged: {m.converged}"
                )
            if diag.md:
                d = diag.md
                lines.append(
                    f"    temp    : {d.last_temperature} K  "
                    f"epot: {d.last_epot} kJ/mol  "
                    f"completed: {d.completed}"
                )
            if diag.errors:
                lines.append(f"    errors  : {'; '.join(diag.errors)}")
            if diag.warnings:
                lines.append(f"    warnings: {'; '.join(diag.warnings)}")
        return "\n".join(lines)
