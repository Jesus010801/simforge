# executors/test_adaptive_reasoning.py
"""
Test de integración del Adaptive Reasoning Engine.

Pipeline completo:
    YAML → Compiler → WorkspaceBuilder → GROMACSExecutor (dry-run)
        → run_adaptive_reasoning()
        → AdaptiveReasoningResult

En dry-run los diagnósticos reportarán outputs incompletos (no existen
archivos reales), por lo que el veredicto será REVIEW o REMEDIATE.
Esto valida la lógica de reasoning sin necesitar GROMACS instalado.

Para simular diferentes escenarios de post-ejecución, el test también
construye diagnósticos sintéticos para verificar las reglas de reasoning:
    - minimización convergida → CONTINUE
    - minimización con Fmax alto → REMEDIATE
    - crash NaN → ABORT
    - LINCS persistente → REMEDIATE
"""

from rich import print
from rich.console import Console

from core.compiler import SimulationCompiler
from builders.workspace_builder import WorkspaceBuilder
from executors.gromacs_executor import (
    GROMACSExecutor,
    GROMACSStepDiagnostic,
    MinimizationMetrics,
    MDMetrics,
    OutputFileStatus,
)
from executors.adaptive_reasoning import run_adaptive_reasoning
from executors.execution_reasoning_models import (
    ReasoningVerdict,
    StepAnalysisVerdict,
)
from executors.execution_state import (
    WorkspaceExecutionState,
    StepExecutionRecord,
    StepStatus,
)

console = Console()


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Pipeline real (dry-run)
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold green]Adaptive Reasoning — Pipeline Real (dry-run)[/bold green]")

compiler  = SimulationCompiler()
result    = compiler.compile("configs/hmg_competition.yaml")
builder   = WorkspaceBuilder()
workspace = builder.build(result)

executor = GROMACSExecutor(workspace_path=workspace, dry_run=True)
exec_state = executor.run()

reasoning_result = run_adaptive_reasoning(
    exec_state  = exec_state,
    diagnostics = executor.all_diagnostics(),
)

verdict_color = {
    ReasoningVerdict.CONTINUE:  "green",
    ReasoningVerdict.REVIEW:    "yellow",
    ReasoningVerdict.REMEDIATE: "bright_yellow",
    ReasoningVerdict.ABORT:     "red",
}.get(reasoning_result.verdict, "white")

print(f"\n  Veredicto global : [{verdict_color}]{reasoning_result.verdict.value}[/{verdict_color}]")
print(f"  Resumen          : {reasoning_result.summary}")
print(f"  Steps analizados : {reasoning_result.n_steps_analyzed}")
print(f"  Steps ok         : [green]{reasoning_result.n_steps_ok}[/green]")
print(f"  Steps con issues : [red]{reasoning_result.n_steps_failed}[/red]")

if reasoning_result.step_analyses:
    console.rule("[bold cyan]Step Analyses[/bold cyan]")
    for a in reasoning_result.step_analyses:
        v_color = {
            StepAnalysisVerdict.OK:           "green",
            StepAnalysisVerdict.NOT_CONVERGED: "yellow",
            StepAnalysisVerdict.NEEDS_REVIEW: "yellow",
            StepAnalysisVerdict.REMEDIABLE:   "bright_yellow",
            StepAnalysisVerdict.FATAL:        "red",
        }.get(a.verdict, "white")
        print(f"\n  [{v_color}]{a.verdict.value:15}[/{v_color}]  {a.step_id}")
        print(f"    {a.interpretation}")
        if a.recommended_action:
            print(f"    → {a.recommended_action}")

if not reasoning_result.remediation_plan.is_empty:
    console.rule("[bold magenta]Remediation Plan[/bold magenta]")
    plan = reasoning_result.remediation_plan
    print(f"\n  Pasos     : {len(plan.steps)}")
    print(f"  Automático: {plan.n_automatic}")
    print(f"  Manual    : {plan.n_manual}")
    print(f"  Esfuerzo  : {plan.estimated_effort}")
    print(f"  Re-run desde: {plan.rerun_from_step}")
    for step in plan.sorted_steps():
        print(f"\n  [{step.priority.value}] {step.remediation_id}")
        print(f"    Problema : {step.problem}")
        print(f"    Acción   : {step.action}")
        if step.mdp_parameter:
            print(
                f"    MDP      : {step.mdp_parameter} "
                f"{step.mdp_current_value} → {step.mdp_suggested_value}"
            )


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Escenarios sintéticos
# ═══════════════════════════════════════════════════════════════════════════════

console.rule("[bold white]Escenarios Sintéticos[/bold white]")


def _make_exec_state(step_ids: list[str]) -> WorkspaceExecutionState:
    return WorkspaceExecutionState(
        workspace_path = "/tmp/test_workspace",
        system_type    = "competitive-inhibition",
        dry_run        = True,
        steps = [
            StepExecutionRecord(
                step_id  = sid,
                step_dir = f"/tmp/test_workspace/steps/{sid}",
                status   = StepStatus.DONE,
            )
            for sid in step_ids
        ],
    )


# ── Escenario A: minimización perfecta ───────────────────────────────────────

state_a = _make_exec_state(["energy_minimization"])
diags_a = {
    "energy_minimization": GROMACSStepDiagnostic(
        step_id      = "energy_minimization",
        stage        = "minimization",
        engine       = "gromacs",
        minimization = MinimizationMetrics(
            converged          = True,
            final_epot         = -850000.0,
            final_fmax         = 42.5,
            n_steps_taken      = 427,
            convergence_reason = "converged_fmax",
        ),
        verdict = "converged",
    )
}
r_a = run_adaptive_reasoning(state_a, diags_a)
print(f"\n  [A] Minimización perfecta     → [{('green' if r_a.verdict == ReasoningVerdict.CONTINUE else 'red')}]{r_a.verdict.value}[/{'green' if r_a.verdict == ReasoningVerdict.CONTINUE else 'red'}]  ✓" if r_a.verdict == ReasoningVerdict.CONTINUE else f"\n  [A] Minimización perfecta     → {r_a.verdict.value}  ✗")

# ── Escenario B: Fmax alto → remediar ────────────────────────────────────────

state_b = _make_exec_state(["energy_minimization"])
diags_b = {
    "energy_minimization": GROMACSStepDiagnostic(
        step_id      = "energy_minimization",
        stage        = "minimization",
        engine       = "gromacs",
        minimization = MinimizationMetrics(
            converged          = False,
            final_epot         = -200000.0,
            final_fmax         = 450.0,
            n_steps_taken      = 50000,
            convergence_reason = "max_steps",
        ),
        verdict = "not_converged",
    )
}
r_b = run_adaptive_reasoning(state_b, diags_b)
ok_b = r_b.verdict == ReasoningVerdict.REMEDIATE
print(f"\n  [B] Fmax=450 kJ/mol/nm        → [{'green' if ok_b else 'red'}]{r_b.verdict.value}[/{'green' if ok_b else 'red'}]  {'✓' if ok_b else '✗'}")

# ── Escenario C: crash NaN → abortar ─────────────────────────────────────────

state_c = _make_exec_state(["equilibration"])
diags_c = {
    "equilibration": GROMACSStepDiagnostic(
        step_id = "equilibration",
        stage   = "equilibration",
        engine  = "gromacs",
        md = MDMetrics(
            completed      = False,
            has_nan_energy = True,
            has_fatal_error = False,
        ),
        verdict = "crashed",
    )
}
r_c = run_adaptive_reasoning(state_c, diags_c)
ok_c = r_c.verdict == ReasoningVerdict.ABORT
print(f"\n  [C] NaN en equilibración      → [{'green' if ok_c else 'red'}]{r_c.verdict.value}[/{'green' if ok_c else 'red'}]  {'✓' if ok_c else '✗'}")

# ── Escenario D: LINCS persistente → remediar ────────────────────────────────

state_d = _make_exec_state(["equilibration"])
diags_d = {
    "equilibration": GROMACSStepDiagnostic(
        step_id = "equilibration",
        stage   = "equilibration",
        engine  = "gromacs",
        md = MDMetrics(
            completed          = True,
            has_lincs_warning  = True,
            n_lincs_warnings   = 12,
            last_temperature   = 301.5,
            last_epot          = -750000.0,
        ),
        verdict = "warning",
    )
}
r_d = run_adaptive_reasoning(state_d, diags_d)
ok_d = r_d.verdict == ReasoningVerdict.REMEDIATE
print(f"\n  [D] LINCS x12 en equilibración → [{'green' if ok_d else 'red'}]{r_d.verdict.value}[/{'green' if ok_d else 'red'}]  {'✓' if ok_d else '✗'}")

# ── Escenario E: temperatura explosiva → abortar ─────────────────────────────

state_e = _make_exec_state(["equilibration"])
diags_e = {
    "equilibration": GROMACSStepDiagnostic(
        step_id = "equilibration",
        stage   = "equilibration",
        engine  = "gromacs",
        md = MDMetrics(
            completed        = False,
            has_exploded     = True,
            last_temperature = 12500.0,
        ),
        verdict = "crashed",
    )
}
r_e = run_adaptive_reasoning(state_e, diags_e)
ok_e = r_e.verdict == ReasoningVerdict.ABORT
print(f"\n  [E] Temperatura 12500K         → [{'green' if ok_e else 'red'}]{r_e.verdict.value}[/{'green' if ok_e else 'red'}]  {'✓' if ok_e else '✗'}")

# ── Escenario F: LINCS esporádico → review ───────────────────────────────────

state_f = _make_exec_state(["production_md"])
diags_f = {
    "production_md": GROMACSStepDiagnostic(
        step_id = "production_md",
        stage   = "production",
        engine  = "gromacs",
        md = MDMetrics(
            completed         = True,
            has_lincs_warning = True,
            n_lincs_warnings  = 2,
            last_temperature  = 299.8,
            last_epot         = -820000.0,
        ),
        verdict = "warning",
    )
}
r_f = run_adaptive_reasoning(state_f, diags_f)
ok_f = r_f.verdict == ReasoningVerdict.REVIEW
print(f"\n  [F] LINCS x2 en producción    → [{'green' if ok_f else 'red'}]{r_f.verdict.value}[/{'green' if ok_f else 'red'}]  {'✓' if ok_f else '✗'}")


# ── Resumen de escenarios ─────────────────────────────────────────────────────

all_ok = all([
    r_a.verdict == ReasoningVerdict.CONTINUE,
    r_b.verdict == ReasoningVerdict.REMEDIATE,
    r_c.verdict == ReasoningVerdict.ABORT,
    r_d.verdict == ReasoningVerdict.REMEDIATE,
    r_e.verdict == ReasoningVerdict.ABORT,
    r_f.verdict == ReasoningVerdict.REVIEW,
])

print()
console.rule(
    f"[bold {'green' if all_ok else 'red'}]"
    f"{'✓ Todos los escenarios correctos' if all_ok else '✗ Algún escenario falló'}"
    f"[/bold {'green' if all_ok else 'red'}]"
)
print()
console.rule("[bold green]Adaptive Reasoning Test Complete[/bold green]")
