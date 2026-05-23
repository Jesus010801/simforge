# core/report_generator.py
"""
Compile Report Generator.

Generates metadata/compile_report.md — a professional scientific summary
of what SimForge compiled, why, what risks were detected, and what to expect.

Called by WorkspaceBuilder after materialization.
"""

from __future__ import annotations

import subprocess
from datetime import datetime
from pathlib import Path
from typing import Optional

from core.compiler_models import CompilationResult
from core.execution_models import StepStage


# ─── Software version helpers ─────────────────────────────────────────────────

def _gmx_version() -> str:
    try:
        out = subprocess.run(
            ["gmx", "--version"],
            capture_output=True, text=True, timeout=5,
        )
        for line in out.stdout.splitlines():
            if "GROMACS version" in line:
                return line.strip().split()[-1]
    except Exception:
        pass
    return "not detected"


def _python_version() -> str:
    import sys
    return f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"


# ─── Runtime estimator ───────────────────────────────────────────────────────

_NS_PER_DAY_CPU  = 5.0   # rough estimate: CPU-only GROMACS
_NS_PER_DAY_GPU  = 80.0  # rough estimate: GPU-accelerated GROMACS


def _estimate_runtime(result: CompilationResult) -> str:
    """
    Very rough estimate — meant to set expectations, not be precise.
    Returns a human-readable string.
    """
    policy     = result.plan.workflow_policy
    prod_ns    = policy.production_time_ns or 10.0
    hardware   = getattr(policy, "hardware", "auto")

    ns_per_day = _NS_PER_DAY_GPU if hardware == "gpu" else _NS_PER_DAY_CPU
    prod_days  = prod_ns / ns_per_day

    # Non-production overhead (prep, param, minim, equil, analysis) ≈ 3–10 h
    overhead_h = 6.0
    total_h    = prod_days * 24 + overhead_h

    if total_h < 2:
        return f"~{total_h*60:.0f} minutes"
    if total_h < 48:
        return f"~{total_h:.0f} hours"
    return f"~{total_h/24:.1f} days"


# ─── Public API ───────────────────────────────────────────────────────────────

def generate_compile_report(
    result: CompilationResult,
    workspace: Path,
    planning_session: Optional[dict] = None,
) -> Path:
    """
    Write metadata/compile_report.md and return its path.
    """
    report_path = workspace / "metadata" / "compile_report.md"
    content     = _build_report(result, workspace, planning_session)
    report_path.write_text(content, encoding="utf-8")
    return report_path


def _build_report(
    result: CompilationResult,
    workspace: Path,
    planning_session: Optional[dict],
) -> str:
    state  = result.state
    plan   = result.plan
    policy = plan.workflow_policy

    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    lines: list[str] = []
    a = lines.append

    a(f"# SimForge Compile Report")
    a(f"")
    a(f"Generated: {now}  ")
    a(f"Workspace: `{workspace}`")
    a(f"")

    # ── 1. System overview ────────────────────────────────────────────────────
    a(f"## System Overview")
    a(f"")
    a(f"| Field | Value |")
    a(f"|---|---|")
    a(f"| System type | `{state.inferred_system_type or 'unknown'}` |")
    a(f"| Project name | `{state.project.name if state.project else 'unnamed'}` |")
    a(f"| Components | {len(state.components)} |")
    a(f"| Steps | {len(plan.steps)} |")
    a(f"| Forcefield (protein) | `{state.forcefields.protein or '?'}` |")
    if state.forcefields.ligands:
        a(f"| Forcefield (ligands) | `{state.forcefields.ligands}` |")
    a(f"| Water model | `{policy.water_model or '?'}` |")
    a(f"| Temperature | {policy.temperature_K} K |")
    a(f"| Pressure | {policy.pressure_bar} bar |")
    a(f"| Production time | {policy.production_time_ns} ns |")
    a(f"| Hardware | `{getattr(policy, 'hardware', 'auto')}` |")
    a(f"")

    # ── 2. Components ─────────────────────────────────────────────────────────
    a(f"## Components")
    a(f"")
    for comp in state.components:
        val_status = "✓ valid" if (comp.validation and comp.validation.is_valid) else "✗ issues"
        a(f"- **{comp.id}** (`{comp.role}`) — `{comp.file}`  [{val_status}]")
    a(f"")

    # ── 3. Inferred biology ───────────────────────────────────────────────────
    a(f"## Inferred Biology")
    a(f"")
    gr = state.global_reasoning
    if gr.notes:
        for note in gr.notes:
            a(f"- {note}")
    else:
        a(f"- No global notes.")
    if plan.special_protocols:
        a(f"- Special protocols: {', '.join(plan.special_protocols)}")
    a(f"")

    # ── 4. Workflow DAG ───────────────────────────────────────────────────────
    a(f"## Workflow DAG")
    a(f"")
    a(f"```")
    current_stage = None
    for i, step in enumerate(result.execution_order, 1):
        stage = step.stage.value if step.stage else "unknown"
        if stage != current_stage:
            if current_stage is not None:
                a(f"")
            a(f"── {stage.upper()} ──")
            current_stage = stage
        deps = f"  ← {', '.join(step.depends_on)}" if step.depends_on else ""
        a(f"  [{i:02d}] {step.step_id}  ({step.engine}){deps}")
    a(f"```")
    a(f"")

    # ── 5. Scientific warnings ────────────────────────────────────────────────
    all_warnings = list(gr.warnings)
    for comp in state.components:
        if comp.reasoning:
            all_warnings.extend(comp.reasoning.warnings)

    all_risks = list(gr.risks)
    for comp in state.components:
        all_risks.extend(comp.all_risks)

    if all_warnings or all_risks:
        a(f"## Scientific Warnings & Risks")
        a(f"")
        for w in all_warnings:
            icon = "⚠" if w.severity.value in ("high", "medium") else "ℹ"
            a(f"- {icon} `[{w.severity.value.upper()}]` {w.message}")
        for r in all_risks:
            a(f"- ✗ `[RISK]` {r.message}")
        a(f"")

    # ── 6. Blocking issues ────────────────────────────────────────────────────
    if plan.blocking_issues:
        a(f"## Blocking Issues")
        a(f"")
        for issue in plan.blocking_issues:
            a(f"- **{issue.severity.value.upper()}** [{issue.source}]: {issue.message}")
        a(f"")

    # ── 7. Expected outputs ───────────────────────────────────────────────────
    a(f"## Expected Outputs")
    a(f"")
    a(f"```")
    a(f"{workspace.name}/")
    a(f"  steps/")
    for step in result.execution_order:
        stage = step.stage.value if step.stage else "?"
        a(f"    {step.step_id}/")
    a(f"  metadata/")
    a(f"    execution_manifest.json")
    a(f"    compile_report.md")
    a(f"    planning_session.json")
    a(f"  analysis/                   ← populated post-run")
    a(f"  execution_state.json        ← written during run")
    a(f"```")
    a(f"")

    # ── 8. Reproducibility ───────────────────────────────────────────────────
    a(f"## Reproducibility")
    a(f"")
    a(f"| Software | Version |")
    a(f"|---|---|")
    a(f"| GROMACS | {_gmx_version()} |")
    a(f"| Python | {_python_version()} |")
    a(f"| SimForge | 0.1.0-dev |")
    a(f"")
    a(f"To reproduce this run:")
    a(f"```bash")
    yaml_hint = state.project.name if state.project else "config"
    a(f"simforge compile configs/{yaml_hint}.yaml")
    a(f"simforge run {workspace}")
    a(f"```")
    a(f"")

    # ── 9. Estimated runtime ─────────────────────────────────────────────────
    a(f"## Runtime Estimate")
    a(f"")
    rt = _estimate_runtime(result)
    a(f"> **{rt}** (rough estimate; depends on system size, hardware, and convergence)")
    a(f"")
    a(f"- Production: {policy.production_time_ns} ns")
    a(f"- Hardware: `{getattr(policy, 'hardware', 'auto')}`")
    a(f"- Timestep: {getattr(policy, 'dt', 0.002)} ps")
    a(f"")

    # ── 10. Planning decisions ────────────────────────────────────────────────
    if planning_session and planning_session.get("answers"):
        a(f"## Planning Decisions")
        a(f"")
        for ans in planning_session["answers"]:
            a(f"- **{ans.get('question_id', '?')}**: {ans.get('selected_label', '?')}")
        a(f"")

    a(f"---")
    a(f"*Generated by SimForge — molecular simulation workflow compiler*")
    a(f"")

    return "\n".join(lines)
