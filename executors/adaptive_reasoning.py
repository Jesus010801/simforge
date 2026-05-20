# executors/adaptive_reasoning.py
"""
Adaptive Reasoning Engine — razonamiento post-ejecución de SimForge.

Lee el estado real de ejecución (WorkspaceExecutionState) y los
diagnósticos GROMACS (GROMACSStepDiagnostic) para producir un
AdaptiveReasoningResult con veredicto, análisis por step, y
plan de remediación accionable.

Principios arquitectónicos:
    - NO ejecuta comandos
    - NO modifica archivos
    - NO reutiliza lógica del decision engine original
    - Solo RAZONA sobre lo que ya ocurrió
    - Produce recomendaciones — un RemediationExecutor futuro las aplica

Pipeline:
    WorkspaceExecutionState
    + dict[str, GROMACSStepDiagnostic]
    + SystemState (opcional, para contexto)
        ↓
    run_adaptive_reasoning()
        ↓
    AdaptiveReasoningResult

Reglas de razonamiento:

    Minimización:
        converged + Fmax < 100          → OK, continuar
        not_converged + Fmax < 100      → REVIEW, probablemente aceptable
        not_converged + Fmax > 100      → REMEDIATE: bajar emtol o más steps
        not_converged + Fmax > 1000     → ABORT: clashes estructurales graves

    Equilibración / Producción:
        completed + sin warnings        → OK
        LINCS esporádico (< 5)          → REVIEW
        LINCS persistente (≥ 5)         → REMEDIATE: revisar parametrización
        NaN en energías                 → ABORT: crash de integración
        Fatal error                     → ABORT
        temperatura > 500K              → ABORT: explosión

    Outputs:
        missing critical (.gro, .edr)   → ABORT si es step bloqueante
        truncated (.xtc < 10KB)         → REVIEW si producción

Interfaz pública:
    run_adaptive_reasoning(
        exec_state:   WorkspaceExecutionState,
        diagnostics:  dict[str, GROMACSStepDiagnostic],
        system_state: SystemState | None = None,
    ) -> AdaptiveReasoningResult
"""

from __future__ import annotations

from typing import Optional

from executors.execution_state import (
    WorkspaceExecutionState,
    StepStatus,
)

from executors.gromacs_executor import (
    GROMACSStepDiagnostic,
    MinimizationMetrics,
    MDMetrics,
)

from executors.adaptive_models import (
    AdaptiveReasoningResult,
    ReasoningVerdict,
    StepAnalysis,
    StepAnalysisVerdict,
    RemediationPlan,
    RemediationStep,
    RemediationTarget,
)

from core.models import Severity


# ═══════════════════════════════════════════════════════════════════════════════
# API pública
# ═══════════════════════════════════════════════════════════════════════════════

def run_adaptive_reasoning(
    exec_state:   WorkspaceExecutionState,
    diagnostics:  dict[str, GROMACSStepDiagnostic],
    system_state=None,   # SystemState | None — opcional para contexto futuro
) -> AdaptiveReasoningResult:
    """
    Punto de entrada del adaptive reasoning.

    Recibe el estado post-ejecución y los diagnósticos GROMACS
    y produce un AdaptiveReasoningResult completo.
    """

    result = AdaptiveReasoningResult(
        workspace_path  = exec_state.workspace_path,
        system_type     = exec_state.system_type,
        n_steps_analyzed = len(exec_state.steps),
    )

    # ── Analizar cada step ────────────────────────────────────────────────────

    remediation_steps: list[RemediationStep] = []

    for record in exec_state.steps:

        diag = diagnostics.get(record.step_id)

        if diag is None:
            # Step sin diagnóstico GROMACS (preparation, assembly, etc.)
            # Solo verificar si falló
            if record.status == StepStatus.FAILED:
                analysis = _analyze_non_gromacs_failure(record)
                result.step_analyses.append(analysis)
                result.n_steps_failed += 1
                if analysis.verdict == StepAnalysisVerdict.FATAL:
                    result.errors.append(
                        f"{record.step_id}: {analysis.interpretation}"
                    )
            else:
                result.n_steps_ok += 1
            continue

        # Step con diagnóstico GROMACS
        analysis, new_remediations = _analyze_gromacs_step(record, diag)
        result.step_analyses.append(analysis)
        remediation_steps.extend(new_remediations)

        if analysis.verdict in (
            StepAnalysisVerdict.OK,
            StepAnalysisVerdict.NOT_CONVERGED,
        ) and analysis.severity != Severity.HIGH:
            result.n_steps_ok += 1
        else:
            result.n_steps_failed += 1

    # ── Construir plan de remediación ─────────────────────────────────────────

    result.remediation_plan = _build_remediation_plan(
        remediation_steps,
        exec_state,
    )

    # ── Veredicto global ──────────────────────────────────────────────────────

    result.verdict = _compute_global_verdict(result)

    # ── Resumen textual ───────────────────────────────────────────────────────

    _build_summary(result)

    return result


# ═══════════════════════════════════════════════════════════════════════════════
# Análisis por step
# ═══════════════════════════════════════════════════════════════════════════════

def _analyze_non_gromacs_failure(
    record,
) -> StepAnalysis:
    """
    Analiza un step que falló sin diagnóstico GROMACS
    (preparation, parametrization, assembly, etc.)
    """
    return StepAnalysis(
        step_id         = record.step_id,
        stage           = "unknown",
        verdict         = StepAnalysisVerdict.FATAL,
        severity        = Severity.HIGH,
        interpretation  = (
            f"Step falló sin diagnóstico disponible. "
            f"Error: {record.error_message or 'desconocido'}"
        ),
        recommended_action = (
            "Revisar stderr del step y corregir manualmente antes de re-ejecutar."
        ),
    )


def _analyze_gromacs_step(
    record,
    diag: GROMACSStepDiagnostic,
) -> tuple[StepAnalysis, list[RemediationStep]]:
    """
    Analiza un step GROMACS con diagnóstico completo.
    Retorna (StepAnalysis, lista de RemediationStep propuestos).
    """
    stage = diag.stage

    if stage == "minimization":
        return _analyze_minimization(record, diag)
    elif stage in ("equilibration", "production"):
        return _analyze_md(record, diag)
    else:
        return _analyze_generic_gromacs(record, diag), []


# ── Minimización ──────────────────────────────────────────────────────────────

def _analyze_minimization(
    record,
    diag: GROMACSStepDiagnostic,
) -> tuple[StepAnalysis, list[RemediationStep]]:

    mini = diag.minimization
    remediations: list[RemediationStep] = []

    # Sin datos de minimización
    if mini is None or mini.convergence_reason == "not_found":
        return StepAnalysis(
            step_id        = record.step_id,
            stage          = "minimization",
            verdict        = StepAnalysisVerdict.NEEDS_REVIEW,
            severity       = Severity.MEDIUM,
            interpretation = "No se pudo leer el log de minimización.",
            recommended_action = "Verificar que em.log existe y es legible.",
        ), []

    fmax = mini.final_fmax
    epot = mini.final_epot

    key_metrics = {}
    if fmax is not None:
        key_metrics["fmax_kJ_mol_nm"] = fmax
    if epot is not None:
        key_metrics["epot_kJ_mol"] = epot
    key_metrics["converged"] = mini.converged
    key_metrics["steps_taken"] = mini.n_steps_taken

    # ── Caso 1: convergió limpiamente ────────────────────────────────────────
    if mini.converged and (fmax is None or fmax < 100):
        return StepAnalysis(
            step_id        = record.step_id,
            stage          = "minimization",
            verdict        = StepAnalysisVerdict.OK,
            severity       = Severity.LOW,
            interpretation = (
                f"Minimización convergida exitosamente "
                f"(Fmax={fmax:.2f} kJ/mol/nm, Epot={epot:.1f} kJ/mol)."
            ),
            recommended_action = "Continuar a equilibración.",
            key_metrics    = key_metrics,
        ), []

    # ── Caso 2: no convergió pero Fmax razonable (< 100) ────────────────────
    if not mini.converged and fmax is not None and fmax < 100:
        return StepAnalysis(
            step_id        = record.step_id,
            stage          = "minimization",
            verdict        = StepAnalysisVerdict.NOT_CONVERGED,
            severity       = Severity.LOW,
            interpretation = (
                f"Minimización no convergió formalmente pero Fmax={fmax:.2f} "
                f"kJ/mol/nm — aceptable para continuar a equilibración."
            ),
            recommended_action = (
                "Continuar a equilibración. Si hay problemas en NVT, "
                "volver a minimización con más nsteps."
            ),
            key_metrics    = key_metrics,
            notes = [
                "Fmax < 100 kJ/mol/nm es generalmente suficiente para iniciar equilibración.",
            ],
        ), []

    # ── Caso 3: no convergió, Fmax moderado (100-1000) → remediar ────────────
    if fmax is not None and 100 <= fmax <= 1000:
        rem = RemediationStep(
            remediation_id      = f"rem_mini_nsteps_{record.step_id}",
            step_id             = record.step_id,
            target              = RemediationTarget.MDP_FILE,
            priority            = Severity.MEDIUM,
            problem             = f"Fmax={fmax:.1f} kJ/mol/nm — demasiado alto para equilibración estable.",
            root_cause          = "nsteps insuficiente o emtol demasiado restrictivo.",
            action              = "Aumentar nsteps en em.mdp (de 50000 a 100000) y re-minimizar.",
            automatic           = True,
            mdp_parameter       = "nsteps",
            mdp_current_value   = "50000",
            mdp_suggested_value = "100000",
            target_file         = "em.mdp",
            notes = [
                "Alternativamente: reducir emtol de 1000 a 5000 kJ/mol/nm "
                "si la convergencia formal no es necesaria.",
            ],
        )
        remediations.append(rem)

        return StepAnalysis(
            step_id             = record.step_id,
            stage               = "minimization",
            verdict             = StepAnalysisVerdict.REMEDIABLE,
            severity            = Severity.MEDIUM,
            interpretation      = (
                f"Fmax={fmax:.1f} kJ/mol/nm — el sistema no está suficientemente "
                f"minimizado para equilibración segura."
            ),
            recommended_action  = "Aplicar remediación: aumentar nsteps en em.mdp.",
            has_remediation     = True,
            remediation_id      = rem.remediation_id,
            key_metrics         = key_metrics,
        ), remediations

    # ── Caso 4: Fmax > 1000 → clashes estructurales → abortar ────────────────
    if fmax is not None and fmax > 1000:
        return StepAnalysis(
            step_id        = record.step_id,
            stage          = "minimization",
            verdict        = StepAnalysisVerdict.FATAL,
            severity       = Severity.HIGH,
            interpretation = (
                f"Fmax={fmax:.1f} kJ/mol/nm — clashes estructurales graves. "
                f"La minimización no puede resolver estos conflictos."
            ),
            recommended_action = (
                "Revisar estructura inicial: buscar átomos solapados, "
                "ligandos mal posicionados, o errores en el assembly. "
                "Corregir con herramientas de preparación (pdb2gmx, VMD, PyMOL) "
                "antes de re-minimizar."
            ),
            key_metrics    = key_metrics,
            notes = [
                "Fmax > 1000 kJ/mol/nm generalmente indica clashes < 1Å entre átomos.",
                "Herramientas: 'gmx check -f system.gro' para detectar overlaps.",
            ],
        ), []

    # ── Fallback ──────────────────────────────────────────────────────────────
    return StepAnalysis(
        step_id        = record.step_id,
        stage          = "minimization",
        verdict        = StepAnalysisVerdict.NEEDS_REVIEW,
        severity       = Severity.MEDIUM,
        interpretation = "Minimización completó pero no se pudo evaluar convergencia.",
        recommended_action = "Revisar em.log manualmente.",
        key_metrics    = key_metrics,
    ), []


# ── Equilibración / Producción ────────────────────────────────────────────────

def _analyze_md(
    record,
    diag: GROMACSStepDiagnostic,
) -> tuple[StepAnalysis, list[RemediationStep]]:

    md = diag.md
    stage = diag.stage
    remediations: list[RemediationStep] = []

    if md is None:
        return StepAnalysis(
            step_id        = record.step_id,
            stage          = stage,
            verdict        = StepAnalysisVerdict.NEEDS_REVIEW,
            severity       = Severity.MEDIUM,
            interpretation = f"No se pudo leer el log de {stage}.",
            recommended_action = f"Verificar que el log de {stage} existe y es legible.",
        ), []

    key_metrics: dict = {}
    if md.last_temperature is not None:
        key_metrics["temperature_K"]  = md.last_temperature
    if md.last_epot is not None:
        key_metrics["epot_kJ_mol"]     = md.last_epot
    if md.last_pressure is not None:
        key_metrics["pressure_bar"]    = md.last_pressure
    key_metrics["completed"]           = md.completed
    key_metrics["lincs_warnings"]      = md.n_lincs_warnings

    # ── Crash: NaN ───────────────────────────────────────────────────────────
    if md.has_nan_energy:
        return StepAnalysis(
            step_id        = record.step_id,
            stage          = stage,
            verdict        = StepAnalysisVerdict.FATAL,
            severity       = Severity.HIGH,
            interpretation = (
                f"NaN en energías durante {stage} — crash de integración. "
                "El sistema divergió numéricamente."
            ),
            recommended_action = (
                "Causas comunes: (1) parámetros de ligando incorrectos, "
                "(2) clashes no resueltos en minimización, "
                "(3) dt demasiado grande. "
                "Verificar parametrización y reducir dt de 0.002 a 0.001 ps como prueba."
            ),
            key_metrics    = key_metrics,
        ), []

    # ── Crash: Fatal error ───────────────────────────────────────────────────
    if md.has_fatal_error:
        return StepAnalysis(
            step_id        = record.step_id,
            stage          = stage,
            verdict        = StepAnalysisVerdict.FATAL,
            severity       = Severity.HIGH,
            interpretation = (
                f"Fatal error de GROMACS en {stage}: {md.fatal_error_msg[:150]}"
            ),
            recommended_action = (
                "Leer el mensaje completo del Fatal error en el log. "
                "Los Fatal errors de GROMACS suelen ser autodescriptivos."
            ),
            key_metrics    = key_metrics,
        ), []

    # ── Crash: explosión de temperatura ─────────────────────────────────────
    if md.has_exploded:
        return StepAnalysis(
            step_id        = record.step_id,
            stage          = stage,
            verdict        = StepAnalysisVerdict.FATAL,
            severity       = Severity.HIGH,
            interpretation = (
                f"Temperatura explosiva ({md.last_temperature:.0f} K) en {stage}. "
                "El sistema se desintegró."
            ),
            recommended_action = (
                "Verificar: (1) equilibración previa completó correctamente, "
                "(2) restraints de posición activos en NVT inicial, "
                "(3) parametrización de ligandos (cargas, constantes de fuerza)."
            ),
            key_metrics    = key_metrics,
        ), []

    # ── LINCS persistente (≥ 5 warnings) ────────────────────────────────────
    if md.has_lincs_warning and md.n_lincs_warnings >= 5:
        rem = RemediationStep(
            remediation_id      = f"rem_lincs_{record.step_id}",
            step_id             = record.step_id,
            target              = RemediationTarget.MDP_FILE,
            priority            = Severity.MEDIUM,
            problem             = (
                f"LINCS persistente: {md.n_lincs_warnings} warnings — "
                "posibles problemas de geometría de enlaces."
            ),
            root_cause          = (
                "Ligando flexible con enlaces difíciles de constrainear, "
                "o parametrización con constantes de fuerza incorrectas."
            ),
            action              = (
                "Opción A: cambiar constraints=h-bonds a constraints=all-bonds en el MDP. "
                "Opción B: reducir dt de 0.002 a 0.001 ps."
            ),
            automatic           = True,
            mdp_parameter       = "dt",
            mdp_current_value   = "0.002",
            mdp_suggested_value = "0.001",
            target_file         = f"{'nvt' if stage == 'equilibration' else 'md'}.mdp",
            notes = [
                "LINCS esporádico (1-4 warnings) es generalmente inofensivo.",
                "LINCS persistente puede indicar parametrización incorrecta del ligando.",
            ],
        )
        remediations.append(rem)

        return StepAnalysis(
            step_id             = record.step_id,
            stage               = stage,
            verdict             = StepAnalysisVerdict.REMEDIABLE,
            severity            = Severity.MEDIUM,
            interpretation      = (
                f"LINCS persistente en {stage} ({md.n_lincs_warnings} warnings): "
                "problemas de geometría de enlaces. Corregible."
            ),
            recommended_action  = "Aplicar remediación: reducir dt o cambiar constraints.",
            has_remediation     = True,
            remediation_id      = rem.remediation_id,
            key_metrics         = key_metrics,
        ), remediations

    # ── LINCS esporádico (< 5) ────────────────────────────────────────────────
    if md.has_lincs_warning and md.n_lincs_warnings < 5:
        return StepAnalysis(
            step_id        = record.step_id,
            stage          = stage,
            verdict        = StepAnalysisVerdict.NEEDS_REVIEW,
            severity       = Severity.LOW,
            interpretation = (
                f"LINCS esporádico en {stage} ({md.n_lincs_warnings} warnings): "
                "generalmente inofensivo pero merece atención."
            ),
            recommended_action = (
                "Revisar si los warnings están concentrados al inicio (aceptable) "
                "o distribuidos durante toda la simulación (problema)."
            ),
            key_metrics    = key_metrics,
        ), []

    # ── Temperatura desviada ─────────────────────────────────────────────────
    if md.last_temperature is not None:
        dev = abs(md.last_temperature - 300.0)
        if dev > 20:
            return StepAnalysis(
                step_id        = record.step_id,
                stage          = stage,
                verdict        = StepAnalysisVerdict.NEEDS_REVIEW,
                severity       = Severity.MEDIUM,
                interpretation = (
                    f"Temperatura final desviada: {md.last_temperature:.1f} K "
                    f"(Δ={dev:.1f} K de 300 K objetivo)."
                ),
                recommended_action = (
                    "Verificar tau_t y grupos de temperatura en el MDP. "
                    "Una desviación persistente puede indicar acoplamiento incorrecto."
                ),
                key_metrics    = key_metrics,
            ), []

    # ── Completó sin issues ───────────────────────────────────────────────────
    if md.completed:
        interpretation = f"{stage.capitalize()} completada sin issues detectados."
        if md.last_temperature:
            interpretation += f" Temperatura final: {md.last_temperature:.1f} K."

        return StepAnalysis(
            step_id        = record.step_id,
            stage          = stage,
            verdict        = StepAnalysisVerdict.OK,
            severity       = Severity.LOW,
            interpretation = interpretation,
            recommended_action = (
                "Continuar al siguiente step."
                if stage == "equilibration"
                else "Proceder al análisis."
            ),
            key_metrics    = key_metrics,
        ), []

    # ── No completó (sin "Finished mdrun") ───────────────────────────────────
    return StepAnalysis(
        step_id        = record.step_id,
        stage          = stage,
        verdict        = StepAnalysisVerdict.NEEDS_REVIEW,
        severity       = Severity.MEDIUM,
        interpretation = (
            f"{stage.capitalize()} no tiene 'Finished mdrun' en el log — "
            "posiblemente interrumpida o incompleta."
        ),
        recommended_action = (
            "Verificar si el job fue interrumpido (timeout en cluster, Ctrl+C). "
            "Si hay checkpoint .cpt, reanudar con '-cpi npt.cpt'."
        ),
        key_metrics    = key_metrics,
    ), []


# ── Genérico ──────────────────────────────────────────────────────────────────

def _analyze_generic_gromacs(
    record,
    diag: GROMACSStepDiagnostic,
) -> StepAnalysis:
    verdict_map = {
        "ok":           StepAnalysisVerdict.OK,
        "converged":    StepAnalysisVerdict.OK,
        "warning":      StepAnalysisVerdict.NEEDS_REVIEW,
        "not_converged": StepAnalysisVerdict.NOT_CONVERGED,
        "incomplete":   StepAnalysisVerdict.REMEDIABLE,
        "crashed":      StepAnalysisVerdict.FATAL,
        "unknown":      StepAnalysisVerdict.NEEDS_REVIEW,
    }
    return StepAnalysis(
        step_id        = record.step_id,
        stage          = diag.stage,
        verdict        = verdict_map.get(diag.verdict, StepAnalysisVerdict.NEEDS_REVIEW),
        severity       = Severity.HIGH if diag.verdict == "crashed" else Severity.LOW,
        interpretation = f"Step GROMACS ({diag.stage}): veredicto = {diag.verdict}",
        notes          = diag.notes,
    )


# ═══════════════════════════════════════════════════════════════════════════════
# Plan de remediación
# ═══════════════════════════════════════════════════════════════════════════════

def _build_remediation_plan(
    steps:      list[RemediationStep],
    exec_state: WorkspaceExecutionState,
) -> RemediationPlan:

    if not steps:
        return RemediationPlan()

    n_auto   = sum(1 for s in steps if s.automatic)
    n_manual = len(steps) - n_auto

    # Estimar esfuerzo según cantidad y tipo
    if n_manual == 0 and n_auto <= 3:
        effort = "minutes"
    elif n_manual <= 2:
        effort = "hours"
    else:
        effort = "days"

    # Desde qué step reiniciar: el primero con remediación HIGH
    rerun_from = None
    high_priority = [s for s in steps if s.priority == Severity.HIGH]
    if high_priority:
        rerun_from = high_priority[0].step_id
    elif steps:
        rerun_from = steps[0].step_id

    return RemediationPlan(
        steps              = steps,
        n_automatic        = n_auto,
        n_manual           = n_manual,
        estimated_effort   = effort,
        rerun_recommended  = True,
        rerun_from_step    = rerun_from,
    )


# ═══════════════════════════════════════════════════════════════════════════════
# Veredicto global
# ═══════════════════════════════════════════════════════════════════════════════

def _compute_global_verdict(
    result: AdaptiveReasoningResult,
) -> ReasoningVerdict:
    """
    Agrega los análisis individuales en un veredicto global.
    Prioridad: ABORT > REMEDIATE > REVIEW > CONTINUE.
    """
    has_fatal      = any(
        a.verdict == StepAnalysisVerdict.FATAL
        for a in result.step_analyses
    )
    has_remediable = any(
        a.verdict == StepAnalysisVerdict.REMEDIABLE
        for a in result.step_analyses
    )
    has_review     = any(
        a.verdict in (StepAnalysisVerdict.NEEDS_REVIEW, StepAnalysisVerdict.NOT_CONVERGED)
        for a in result.step_analyses
    )

    if has_fatal:
        return ReasoningVerdict.ABORT

    if has_remediable or not result.remediation_plan.is_empty:
        return ReasoningVerdict.REMEDIATE

    if has_review:
        return ReasoningVerdict.REVIEW

    return ReasoningVerdict.CONTINUE


# ═══════════════════════════════════════════════════════════════════════════════
# Resumen textual
# ═══════════════════════════════════════════════════════════════════════════════

def _build_summary(result: AdaptiveReasoningResult) -> None:
    """Construye el campo summary y los notes/warnings/errors globales."""

    verdict_messages = {
        ReasoningVerdict.CONTINUE:  "Sistema ejecutó sin problemas. Listo para continuar.",
        ReasoningVerdict.REVIEW:    "Ejecución completó con advertencias menores. Revisar antes de continuar.",
        ReasoningVerdict.REMEDIATE: "Se detectaron problemas corregibles. Aplicar plan de remediación.",
        ReasoningVerdict.ABORT:     "Error crítico detectado. Requiere corrección manual antes de re-ejecutar.",
    }

    result.summary = verdict_messages.get(result.verdict, "")

    # Agregar notas de steps con problemas
    for analysis in result.step_analyses:
        if analysis.verdict == StepAnalysisVerdict.FATAL:
            result.errors.append(
                f"[{analysis.step_id}] {analysis.interpretation}"
            )
        elif analysis.verdict == StepAnalysisVerdict.REMEDIABLE:
            result.warnings.append(
                f"[{analysis.step_id}] {analysis.interpretation}"
            )
        elif analysis.verdict == StepAnalysisVerdict.NEEDS_REVIEW:
            result.notes.append(
                f"[{analysis.step_id}] {analysis.interpretation}"
            )

    # Nota de remediación
    plan = result.remediation_plan
    if not plan.is_empty:
        result.notes.append(
            f"Plan de remediación: {len(plan.steps)} paso(s), "
            f"{plan.n_automatic} automático(s), "
            f"{plan.n_manual} manual(es). "
            f"Esfuerzo estimado: {plan.estimated_effort}."
        )
        if plan.rerun_from_step:
            result.notes.append(
                f"Re-ejecutar desde: {plan.rerun_from_step}"
            )
