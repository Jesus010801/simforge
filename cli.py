#!/usr/bin/env python3
# cli.py — SimForge command-line interface

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.table import Table
from rich import box
from rich.panel import Panel
from rich.text import Text
from rich.tree import Tree

app     = Console()
cli     = typer.Typer(
    name="simforge",
    help="Molecular simulation workflow compiler.",
    add_completion=False,
)


# ═══════════════════════════════════════════════════════════════════════════════
# Compile helpers
# ═══════════════════════════════════════════════════════════════════════════════

_LIGAND_ROLES = {"ligand", "competitive_ligand", "substrate", "cofactor", "inhibitor"}

_STAGE_COLORS = {
    "preparation":       "cyan",
    "parametrization":   "blue",
    "assembly":          "blue",
    "minimization":      "yellow",
    "equilibration":     "yellow",
    "production":        "green",
    "analysis":          "magenta",
    "enhanced_sampling": "red",
    "membrane_embedding":"red",
    "validation":        "dim",
    "review":            "dim",
}


def _stage_ok(n: int, total: int, label: str, elapsed: float, note: str = "") -> None:
    note_str = f"  [dim]{note}[/dim]" if note else ""
    app.print(
        f"  [green]✓[/green] [dim][{n}/{total}][/dim] "
        f"{label:<38} [dim]{elapsed:.2f}s[/dim]{note_str}"
    )


def _collect_scientific_warnings(state, plan) -> list[tuple[str, str]]:
    """Returns list of (message, color) tuples."""
    warns: list[tuple[str, str]] = []

    # Membrane-associated without membrane enabled
    for comp in state.components:
        ctx = getattr(comp, "biological_context", None) or []
        if "membrane_associated" in ctx:
            mem = getattr(state.environment, "membrane", None)
            if not (mem and getattr(mem, "enabled", False)):
                warns.append((
                    f"membrane_associated protein ({comp.id}) without membrane in config",
                    "yellow",
                ))

    # High-flexibility ligands
    for comp in state.components:
        if comp.role in _LIGAND_ROLES:
            desc = getattr(comp, "descriptors", None)
            if desc and desc.flexibility_class in ("flexible", "very_flexible"):
                warns.append((
                    f"High-flexibility ligand {comp.id} "
                    f"({desc.n_rotatable_bonds} rotatable bonds) — consider longer equilibration",
                    "yellow",
                ))
            if desc and abs(desc.net_charge) > 2:
                warns.append((
                    f"Ligand {comp.id} has net charge {desc.net_charge:+d} — verify ion neutralization",
                    "yellow",
                ))

    # Blocking issues from plan
    for issue in plan.blocking_issues:
        color = "red" if issue.severity.value in ("critical", "error") else "yellow"
        warns.append((issue.message, color))

    # Enhanced sampling
    if plan.workflow_policy.enhanced_sampling:
        warns.append((
            f"Enhanced sampling: {plan.workflow_policy.sampling_method} — "
            f"verify replica count and temperature range",
            "cyan",
        ))

    # Long run
    if plan.workflow_policy.production_time_ns >= 100:
        warns.append((
            f"Long production ({plan.workflow_policy.production_time_ns} ns) — "
            f"verify disk space and GPU availability",
            "yellow",
        ))

    return warns


def _build_dag_tree(steps: list) -> Tree:
    """Build a Rich Tree from execution order using primary-parent algorithm."""
    step_ids   = [s.step_id for s in steps]
    id_to_step = {s.step_id: s for s in steps}

    # Each node's primary parent = its latest dependency in execution order
    primary_parent: dict[str, str | None] = {}
    for step in steps:
        deps = [d for d in (step.depends_on or []) if d in id_to_step]
        primary_parent[step.step_id] = (
            None if not deps
            else max(deps, key=lambda d: step_ids.index(d))
        )

    children: dict[str, list[str]] = {sid: [] for sid in step_ids}
    roots: list[str] = []
    for sid, parent in primary_parent.items():
        if parent is None:
            roots.append(sid)
        else:
            children[parent].append(sid)

    def _render(node, sid: str) -> None:
        step  = id_to_step.get(sid)
        color = _STAGE_COLORS.get(step.stage.value if step else "", "white")
        stage = f"  [dim]{step.stage.value}[/dim]" if step else ""
        child_node = node.add(f"[{color}]{sid}[/{color}]{stage}")
        for cid in children.get(sid, []):
            _render(child_node, cid)

    tree = Tree("[bold dim]Workflow DAG[/bold dim]")
    for rid in roots:
        _render(tree, rid)
    return tree


def _show_system_summary(result) -> None:
    state  = result.state
    plan   = result.plan
    policy = plan.workflow_policy

    comp_lines = []
    for c in state.components:
        desc = getattr(c, "descriptors", None)
        extra = ""
        if desc:
            if desc.n_heavy_atoms:
                extra += f"  {desc.n_heavy_atoms} heavy atoms"
            if desc.flexibility_class != "unknown":
                extra += f"  [{desc.flexibility_class}]"
            if desc.net_charge:
                extra += f"  charge {desc.net_charge:+d}"
        comp_lines.append(f"  [dim]•[/dim] {c.id}  [dim]{c.role}[/dim]{extra}")

    dur   = policy.production_time_ns
    gpu_h = max(0.1, dur * 0.2)
    cpu_h = dur * 4.0
    physio = " [dim](37 °C — physiological)[/dim]" if abs(policy.temperature_K - 309.65) < 0.1 else ""

    lines = [
        f"  System type     [bold]{state.inferred_system_type}[/bold]",
        f"  Temperature     {policy.temperature_K} K{physio}",
        f"  Pressure        {policy.pressure_bar} bar",
        f"  Duration        {dur} ns  ({int(dur / policy.timestep_ps * 1000):,} steps)",
        f"  Workflow        {len(plan.steps)} steps",
        "",
        "  [bold]Components[/bold]",
    ] + comp_lines + [
        "",
        f"  [bold]Runtime estimate[/bold] [dim](rough — depends on hardware & system size)[/dim]",
        f"  GPU ≈ {gpu_h:.1f} h     CPU ≈ {cpu_h:.0f} h",
    ]

    if policy.enhanced_sampling:
        lines.append(f"\n  [red]Enhanced sampling:[/red] {policy.sampling_method}")

    app.print(Panel("\n".join(lines), title="System Summary", border_style="cyan", padding=(0, 2)))


def _show_scientific_warnings(warnings: list[tuple[str, str]]) -> None:
    if not warnings:
        return
    table = Table(box=box.SIMPLE, show_header=False, padding=(0, 1))
    table.add_column("icon", width=2, no_wrap=True)
    table.add_column("message")
    for msg, color in warnings:
        icon = "⚠" if color == "yellow" else ("✗" if color == "red" else "◆")
        table.add_row(f"[{color}]{icon}[/{color}]", f"[{color}]{msg}[/{color}]")
    app.print(Panel(table, title="Scientific Warnings", border_style="yellow", padding=(0, 1)))


def _inspect_variants(yaml_path: Path) -> None:
    from core.variant_compiler import parse_variant_yaml

    app.print(f"\n[bold cyan]SimForge Inspect[/bold cyan]  [magenta]variants[/magenta]  "
              f"[dim]{yaml_path}[/dim]\n")

    try:
        manifest = parse_variant_yaml(yaml_path)
    except Exception as e:
        app.print(f"[red]Parse failed:[/red] {e}")
        raise typer.Exit(1)

    app.print(f"  Project:   [bold]{manifest.project_name}[/bold]")
    app.print(f"  Variants:  {len(manifest.variants)}\n")

    for spec in manifest.variants:
        exists = "[green]✓[/green]" if Path(spec.file).exists() else "[red]✗[/red]"
        label  = spec.label or spec.variant_id
        app.print(f"  {exists}  [bold]{spec.variant_id}[/bold]  [dim]({label})[/dim]")
        app.print(f"       [dim]{spec.file}[/dim]")

    sw = manifest.shared_workflow
    app.print("\n  [bold]Shared workflow:[/bold]")
    for k, v in sw.items():
        if k != "analysis":
            app.print(f"    {k}: {v}")
    if "analysis" in sw:
        app.print(f"    analysis: {sw['analysis']}")
    app.print()


def _write_compile_report(result, workspace: Path, yaml_path: Path, warnings: list) -> None:
    """Write workspace/metadata/compile_report.md."""
    from datetime import datetime

    state  = result.state
    plan   = result.plan
    policy = plan.workflow_policy

    comps_md = "\n".join(
        f"- `{c.id}` ({c.role})" + (f" — `{c.file}`" if c.file else "")
        for c in state.components
    )

    steps_md = "| # | Step | Stage | Engine |\n|---|------|-------|--------|\n"
    for i, step in enumerate(result.execution_order, 1):
        steps_md += f"| {i} | `{step.step_id}` | {step.stage.value} | {step.engine} |\n"

    warns_md = "\n".join(f"- {msg}" for msg, _ in warnings) or "None"

    env   = state.environment
    water = getattr(getattr(env, "solvent", None), "water_model", "?")
    ions  = getattr(getattr(env, "ions", None), "concentration", "?")

    report = f"""# SimForge Compile Report

**Generated:** {datetime.now().isoformat(timespec="seconds")}
**Config:** `{yaml_path}`
**SimForge version:** 0.1.0

---

## System

| Field | Value |
|-------|-------|
| Type | {state.inferred_system_type} |
| Temperature | {policy.temperature_K} K |
| Pressure | {policy.pressure_bar} bar |
| Duration | {policy.production_time_ns} ns |
| Water model | {water} |
| Ion concentration | {ions} M |
| Hardware | {getattr(policy, "hardware", "auto")} |
| Enhanced sampling | {"yes — " + policy.sampling_method if policy.enhanced_sampling else "no"} |

## Components

{comps_md}

## Workflow ({len(result.execution_order)} steps)

{steps_md}
## Scientific Warnings

{warns_md}

## Reproducibility

- Config snapshot: `{yaml_path}`
- Execution manifest: `metadata/execution_manifest.json`
- Planning session: `metadata/planning_session.json`
"""

    (workspace / "metadata" / "compile_report.md").write_text(report)


# ═══════════════════════════════════════════════════════════════════════════════
# compile
# ═══════════════════════════════════════════════════════════════════════════════

@cli.command()
def compile(
    yaml_path:  Path = typer.Argument(..., help="YAML config file to compile."),
    output_dir: str  = typer.Option("simforge_runs", "--output-dir", "-o", help="Root output directory."),
    no_build:   bool = typer.Option(False, "--no-build", help="Only compile — skip workspace materialization."),
    no_plan:    bool = typer.Option(False, "--no-plan",  help="Skip scientific planning dialogue."),
):
    """Compile a YAML config into an executable workspace."""
    from core.compiler import SimulationCompiler
    from core.parser import parse_yaml
    from builders.workspace_builder import WorkspaceBuilder

    if not yaml_path.exists():
        app.print(Panel(
            f"[red]File not found:[/red] {yaml_path}\n\n"
            f"Suggested fix: check the path or run [dim]simforge init[/dim]",
            title="❌ Compile Failed", border_style="red",
        ))
        raise typer.Exit(1)

    # ── Variants dispatch ─────────────────────────────────────────────────────
    try:
        from core.variant_compiler import is_variant_yaml
        if is_variant_yaml(yaml_path):
            _compile_variants(yaml_path, output_dir, no_build)
            return
    except Exception:
        pass  # not a variants YAML — continue normally

    t_start = time.perf_counter()
    app.print(Panel(
        f"[bold cyan]SimForge Compiler[/bold cyan]  [dim]{yaml_path}[/dim]",
        border_style="cyan", padding=(0, 2),
    ))

    # ── [1/6] Parse YAML ──────────────────────────────────────────────────────
    t = time.perf_counter()
    try:
        with app.status("  [dim][1/6][/dim] Parsing YAML..."):
            state = parse_yaml(str(yaml_path))
    except Exception as e:
        app.print(Panel(
            f"[red]Parse error:[/red] {e}\n\nCheck YAML syntax and component file paths.",
            title="❌ Compile Failed", border_style="red",
        ))
        raise typer.Exit(1)
    _stage_ok(1, 6, "Parsing YAML", time.perf_counter() - t,
              note=f"{len(state.components)} component(s)  system: {state.inferred_system_type}")

    # ── [2/6] Scientific planning ──────────────────────────────────────────────
    t = time.perf_counter()
    session_record = _run_planning_dialogue(state, yaml_path, skip=no_plan)
    if session_record is None:
        raise typer.Exit(1)
    _stage_ok(2, 6, "Scientific planning", time.perf_counter() - t,
              note="skipped" if (no_plan or session_record.skipped) else
                   f"{len(session_record.answers or [])} answer(s) applied")

    # ── [3/6] Build execution plan ────────────────────────────────────────────
    t = time.perf_counter()
    try:
        with app.status("  [dim][3/6][/dim] Building execution plan..."):
            result = SimulationCompiler().compile_from_state(state)
    except Exception as e:
        app.print(Panel(
            f"[red]Compilation error:[/red] {e}",
            title="❌ Compile Failed", border_style="red",
        ))
        raise typer.Exit(1)
    _stage_ok(3, 6, "Building execution plan", time.perf_counter() - t,
              note=f"{len(result.plan.steps)} steps planned")

    # ── [4/6] Workflow DAG + warnings ─────────────────────────────────────────
    t = time.perf_counter()
    warnings = _collect_scientific_warnings(result.state, result.plan)
    _stage_ok(4, 6, "Computing workflow DAG", time.perf_counter() - t,
              note=f"{len(result.execution_order)} nodes  "
                   f"{len(warnings)} warning(s)  "
                   f"{len(result.plan.blocking_issues)} blocking issue(s)")

    # ── [5/6] Materialize workspace ───────────────────────────────────────────
    workspace = None
    if not no_build:
        t = time.perf_counter()
        try:
            with app.status("  [dim][5/6][/dim] Materializing workspace..."):
                workspace = WorkspaceBuilder().build(result, output_dir=output_dir)
        except Exception as e:
            app.print(Panel(
                f"[red]Workspace error:[/red] {e}",
                title="❌ Compile Failed", border_style="red",
            ))
            raise typer.Exit(1)
        _stage_ok(5, 6, "Materializing workspace", time.perf_counter() - t,
                  note=str(workspace))

        # ── [6/6] Write reports ───────────────────────────────────────────────
        t = time.perf_counter()
        with app.status("  [dim][6/6][/dim] Writing reports..."):
            _save_planning_session(session_record, workspace / "metadata" / "planning_session.json")
            _write_compile_report(result, workspace, yaml_path, warnings)
        _stage_ok(6, 6, "Writing reports", time.perf_counter() - t,
                  note="compile_report.md  execution_manifest.json")
    else:
        app.print("  [dim]--no-build: steps 5–6 skipped[/dim]")

    # ── System summary ────────────────────────────────────────────────────────
    app.print()
    _show_system_summary(result)

    # ── Scientific warnings ───────────────────────────────────────────────────
    _show_scientific_warnings(warnings)

    # ── DAG preview ───────────────────────────────────────────────────────────
    app.print(Panel(
        _build_dag_tree(result.execution_order),
        title="Workflow DAG", border_style="dim", padding=(0, 2),
    ))

    # ── Blocking issues ───────────────────────────────────────────────────────
    if result.plan.blocking_issues:
        app.print("[red]Blocking issues:[/red]")
        for issue in result.plan.blocking_issues:
            app.print(f"  [red]✗[/red] [{issue.severity.value}] {issue.source}: {issue.message}")
        app.print()

    # ── Final panel ───────────────────────────────────────────────────────────
    total = time.perf_counter() - t_start
    if workspace:
        app.print(Panel(
            f"[bold green]✓[/bold green] Workspace  → [cyan]{workspace}[/cyan]\n"
            f"  Report     → [dim]{workspace}/metadata/compile_report.md[/dim]\n"
            f"  Manifest   → [dim]{workspace}/metadata/execution_manifest.json[/dim]\n\n"
            f"Compiled in [bold]{total:.2f}s[/bold]\n\n"
            f"Next:  [dim]simforge run {workspace}[/dim]",
            title="✓ Done", border_style="green", padding=(0, 2),
        ))
    else:
        app.print(f"\n[dim]--no-build: workspace not materialized  ({total:.2f}s)[/dim]\n")


# ═══════════════════════════════════════════════════════════════════════════════
# variants
# ═══════════════════════════════════════════════════════════════════════════════

def _compile_variants(yaml_path: Path, output_dir: str, no_build: bool) -> None:
    """Handle a variants YAML: compile each variant + build comparative workspace."""
    from core.variant_compiler import (
        parse_variant_yaml, compile_variants, build_comparative_workspace,
    )

    app.print(Panel(
        f"[bold cyan]SimForge Variants[/bold cyan]  [dim]{yaml_path}[/dim]",
        border_style="magenta", padding=(0, 2),
    ))

    try:
        manifest = parse_variant_yaml(yaml_path)
    except Exception as e:
        app.print(f"[red]Variants parse failed:[/red] {e}")
        raise typer.Exit(1)

    app.print(f"  Project:  [bold]{manifest.project_name}[/bold]")
    app.print(f"  Variants: {len(manifest.variants)}")
    for spec in manifest.variants:
        label   = spec.label or spec.variant_id
        exists  = "[green]✓[/green]" if Path(spec.file).exists() else "[red]✗[/red]"
        app.print(f"    {exists}  [cyan]{spec.variant_id}[/cyan]  ({label})  [dim]{spec.file}[/dim]")
    app.print()

    with app.status("  Compiling variants..."):
        manifest = compile_variants(manifest, output_dir=output_dir, no_build=no_build)

    app.print()
    for spec in manifest.variants:
        if spec.variant_id in manifest.errors:
            app.print(f"  [red]✗[/red]  {spec.variant_id}  — {manifest.errors[spec.variant_id]}")
        elif spec.variant_id in manifest.workspaces:
            app.print(
                f"  [green]✓[/green]  [bold]{spec.variant_id}[/bold]  →  "
                f"[cyan]{manifest.workspaces[spec.variant_id]}[/cyan]"
            )

    if not no_build and manifest.workspaces:
        comp_dir = build_comparative_workspace(manifest, output_dir=output_dir)
        app.print(Panel(
            f"  Comparative workspace:  [cyan]{comp_dir}[/cyan]\n\n"
            f"  Run each variant with:  [dim]simforge run <workspace>[/dim]\n"
            f"  Then compare results:   [dim]see {comp_dir}/README.md[/dim]",
            title="✓ Variants compiled",
            border_style="magenta",
        ))


# ═══════════════════════════════════════════════════════════════════════════════
# run
# ═══════════════════════════════════════════════════════════════════════════════

def _execute_workspace(
    workspace:     Path,
    dry_run:       bool,
    executor_type: str,
    no_confirm:    bool = False,
) -> None:
    """Shared execution logic for `run` and `dry-run` commands."""
    if not workspace.exists():
        app.print(f"[red]Error:[/red] Workspace not found: {workspace}")
        raise typer.Exit(1)

    manifest_file = workspace / "metadata" / "execution_manifest.json"
    if not manifest_file.exists():
        app.print(f"[red]Error:[/red] No execution manifest found. Run [cyan]simforge compile[/cyan] first.")
        raise typer.Exit(1)

    if dry_run:
        app.print(f"\n[bold cyan]SimForge[/bold cyan]  [yellow]dry-run[/yellow]  "
                  f"[green]{workspace}[/green]\n")
    else:
        app.print(f"\n[bold cyan]SimForge[/bold cyan]  [bold]run[/bold]  "
                  f"[green]{workspace}[/green]\n")

    # Pre-run summary for real execution
    if not dry_run:
        _print_prerun_summary(workspace, manifest_file)
        if not no_confirm:
            confirmed = typer.confirm("\nProceed with real execution?")
            if not confirmed:
                app.print("[dim]Aborted.[/dim]")
                raise typer.Exit(0)

    if executor_type == "gromacs":
        from executors.gromacs_executor import GROMACSExecutor
        executor = GROMACSExecutor(workspace, dry_run=dry_run)
    else:
        from executors.shell_executor import ShellExecutor
        executor = ShellExecutor(workspace, dry_run=dry_run)

    app.print()

    try:
        state = executor.run()
    except Exception as e:
        app.print(f"[red]Execution error:[/red] {e}")
        raise typer.Exit(1)

    done    = sum(1 for s in state.steps if s.status.value == "done")
    failed  = sum(1 for s in state.steps if s.status.value == "failed")
    skipped = sum(1 for s in state.steps if s.status.value in ("skipped", "blocked"))
    total   = len(state.steps)

    app.print()
    if failed == 0:
        app.print(Panel(
            f"  [green]✓[/green]  {done}/{total} steps completed  "
            f"[dim](skipped: {skipped})[/dim]",
            border_style="green",
            title="Run complete",
        ))
    else:
        failed_lines = "\n".join(
            f"  [red]·[/red]  {s.step_id}  exit={s.exit_code}"
            for s in state.steps if s.status.value == "failed"
        )
        app.print(Panel(
            f"  [red]✗[/red]  {failed} step(s) failed  |  {done} done  |  {skipped} skipped\n"
            + failed_lines,
            border_style="red",
            title="Run failed",
        ))

    if not dry_run and done > 0:
        _print_output_paths(workspace)

    if failed:
        raise typer.Exit(1)


def _print_prerun_summary(workspace: Path, manifest_file: Path) -> None:
    """Show estimated cost/time before a real run."""
    manifest    = json.loads(manifest_file.read_text())
    n_steps     = manifest.get("n_steps", "?")
    system_type = manifest.get("system_type", "unknown")

    runtime_est = None
    report_path = workspace / "metadata" / "compile_report.md"
    if report_path.exists():
        for line in report_path.read_text().splitlines():
            if line.startswith(">") and "estimate" in line.lower():
                runtime_est = line.lstrip("> ").strip()
                break

    lines = [
        f"  Workspace     [cyan]{workspace}[/cyan]",
        f"  System type   {system_type}",
        f"  Steps         {n_steps}",
    ]
    if runtime_est:
        lines.append(f"  Runtime est.  {runtime_est}")
    lines += [
        "",
        "  [dim]This will execute real GROMACS commands and may run for hours.[/dim]",
    ]
    app.print(Panel("\n".join(lines), title="Pre-run summary", border_style="yellow"))


def _print_output_paths(workspace: Path) -> None:
    """Show where outputs landed after a run."""
    app.print("[bold]Results:[/bold]")
    app.print(f"  Workspace:  [cyan]{workspace}[/cyan]")

    steps_dir = workspace / "steps"
    if steps_dir.exists():
        for step_dir in sorted(steps_dir.iterdir(), reverse=True):
            if not step_dir.is_dir():
                continue
            xtc = step_dir / "md.xtc"
            gro = step_dir / "md.gro"
            if xtc.exists():
                app.print(f"  Trajectory: [cyan]{xtc}[/cyan]")
            if gro.exists():
                app.print(f"  Structure:  [cyan]{gro}[/cyan]")
            if xtc.exists() or gro.exists():
                break

    analysis_dir = workspace / "analysis"
    if analysis_dir.exists():
        app.print(f"  Analysis:   [cyan]{analysis_dir}[/cyan]")

    report = workspace / "metadata" / "compile_report.md"
    if report.exists():
        app.print(f"  Report:     [cyan]{report}[/cyan]")
    app.print()


@cli.command()
def run(
    workspace:     Path = typer.Argument(..., help="Workspace directory to execute."),
    executor_type: str  = typer.Option("shell", "--executor", "-e", help="Executor type: shell | gromacs."),
    no_confirm:    bool = typer.Option(False, "--no-confirm", help="Skip confirmation prompt (for automation)."),
):
    """Execute a compiled workspace with real GROMACS commands."""
    _execute_workspace(workspace, dry_run=False, executor_type=executor_type, no_confirm=no_confirm)


@cli.command(name="dry-run")
def dry_run_cmd(
    workspace:     Path = typer.Argument(..., help="Workspace directory to simulate."),
    executor_type: str  = typer.Option("shell", "--executor", "-e", help="Executor type: shell | gromacs."),
):
    """Simulate execution without running real GROMACS commands."""
    _execute_workspace(workspace, dry_run=True, executor_type=executor_type)


# ═══════════════════════════════════════════════════════════════════════════════
# validate
# ═══════════════════════════════════════════════════════════════════════════════

@cli.command()
def validate(
    yaml_path: Path = typer.Argument(..., help="YAML config to validate."),
):
    """Parse and validate a config without generating a workspace."""
    from core.parser import parse_yaml
    from core.compiler import SimulationCompiler

    if not yaml_path.exists():
        app.print(f"[red]Error:[/red] File not found: {yaml_path}")
        raise typer.Exit(1)

    app.print(f"\n[bold cyan]SimForge Validate[/bold cyan]  [dim]{yaml_path}[/dim]\n")

    try:
        with app.status("  Parsing..."):
            state = parse_yaml(str(yaml_path))
        app.print(f"  [green]✓[/green] Parse        {len(state.components)} component(s)  system: {state.inferred_system_type}")
    except Exception as e:
        app.print(Panel(
            f"[red]Parse failed:[/red] {e}",
            title="❌ Validation Failed", border_style="red",
        ))
        raise typer.Exit(1)

    try:
        with app.status("  Compiling plan..."):
            result = SimulationCompiler().compile_from_state(state)
        app.print(f"  [green]✓[/green] Plan         {len(result.plan.steps)} steps  {len(result.execution_order)} in DAG")
    except Exception as e:
        app.print(Panel(
            f"[red]Compilation failed:[/red] {e}",
            title="❌ Validation Failed", border_style="red",
        ))
        raise typer.Exit(1)

    warnings = _collect_scientific_warnings(result.state, result.plan)
    _show_scientific_warnings(warnings)

    if result.plan.blocking_issues:
        app.print(f"\n  [red]✗[/red] {len(result.plan.blocking_issues)} blocking issue(s):")
        for issue in result.plan.blocking_issues:
            app.print(f"    [red]•[/red] [{issue.severity.value}] {issue.message}")
        raise typer.Exit(1)

    app.print(
        f"\n  [bold green]✓ Valid[/bold green]  "
        f"{len(warnings)} warning(s)  0 blocking issues\n"
    )


# ═══════════════════════════════════════════════════════════════════════════════
# inspect
# ═══════════════════════════════════════════════════════════════════════════════

@cli.command()
def inspect(
    yaml_path: Path = typer.Argument(..., help="YAML config to inspect."),
):
    """Show parsed IR, descriptors, and DAG without writing any files."""
    from core.parser import parse_yaml
    from core.compiler import SimulationCompiler

    if not yaml_path.exists():
        app.print(f"[red]Error:[/red] File not found: {yaml_path}")
        raise typer.Exit(1)

    # Variants dispatch
    try:
        from core.variant_compiler import is_variant_yaml, parse_variant_yaml
        if is_variant_yaml(yaml_path):
            _inspect_variants(yaml_path)
            return
    except Exception:
        pass

    app.print(f"\n[bold cyan]SimForge Inspect[/bold cyan]  [dim]{yaml_path}[/dim]\n")

    try:
        with app.status("  Parsing + compiling..."):
            state  = parse_yaml(str(yaml_path))
            result = SimulationCompiler().compile_from_state(state)
    except Exception as e:
        app.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(1)

    # System summary
    _show_system_summary(result)

    # Component descriptors table
    desc_table = Table(box=box.SIMPLE, show_header=True, header_style="bold dim")
    desc_table.add_column("Component")
    desc_table.add_column("Role")
    desc_table.add_column("Heavy atoms", justify="right")
    desc_table.add_column("Flexibility")
    desc_table.add_column("Charge", justify="right")
    desc_table.add_column("H-bond D/A", justify="right")
    for comp in state.components:
        d = getattr(comp, "descriptors", None)
        desc_table.add_row(
            comp.id,
            comp.role,
            str(d.n_heavy_atoms) if d and d.n_heavy_atoms else "—",
            d.flexibility_class if d else "—",
            f"{d.net_charge:+d}" if d else "—",
            f"{d.hb_donors}/{d.hb_acceptors}" if d else "—",
        )
    app.print(Panel(desc_table, title="Component Descriptors", border_style="dim"))

    # DAG preview
    app.print(Panel(
        _build_dag_tree(result.execution_order),
        title="Workflow DAG", border_style="dim", padding=(0, 2),
    ))

    # Key step params (IR snapshot)
    params_table = Table(box=box.SIMPLE, show_header=True, header_style="bold dim")
    params_table.add_column("Step", style="white")
    params_table.add_column("Key IR params", overflow="fold", style="dim")
    _KEY_PARAMS = {"emtol", "nsteps", "temperature", "box_type", "box_distance", "hardware", "concentration"}
    for step in result.execution_order:
        shown = {k: v for k, v in (step.params or {}).items() if k in _KEY_PARAMS}
        if shown:
            params_table.add_row(step.step_id, "  ".join(f"{k}={v}" for k, v in shown.items()))
    app.print(Panel(params_table, title="Step Parameters (IR snapshot)", border_style="dim"))

    # Warnings
    _show_scientific_warnings(_collect_scientific_warnings(state, result.plan))
    app.print()


# ═══════════════════════════════════════════════════════════════════════════════
# status
# ═══════════════════════════════════════════════════════════════════════════════

_STATUS_STYLE = {
    "done":    ("green",  "✓"),
    "failed":  ("red",    "✗"),
    "running": ("yellow", "▶"),
    "skipped": ("dim",    "–"),
    "blocked": ("dim",    "⊘"),
    "pending": ("dim",    "·"),
}

@cli.command()
def status(
    workspace: Path = typer.Argument(..., help="Workspace directory to inspect."),
    verbose: bool   = typer.Option(False, "--verbose", "-v", help="Show timing and error details."),
):
    """Show execution status of a workspace."""

    state_file = workspace / "execution_state.json"
    manifest_file = workspace / "metadata" / "execution_manifest.json"

    if not state_file.exists():
        if manifest_file.exists():
            manifest = json.loads(manifest_file.read_text())
            app.print(f"\n[cyan]{workspace}[/cyan]  [dim]compiled, not yet executed[/dim]")
            app.print(f"  {manifest['n_steps']} steps  system: {manifest.get('system_type', '?')}")
            app.print(f"\nRun:  [dim]simforge run {workspace}[/dim]\n")
        else:
            app.print(f"[red]Error:[/red] No state or manifest found in {workspace}")
            raise typer.Exit(1)
        return

    raw   = json.loads(state_file.read_text())
    steps = raw.get("steps", [])

    # Enrich steps with stage info from manifest (not stored in execution_state)
    if manifest_file.exists():
        manifest_steps = {e["step_id"]: e for e in json.loads(manifest_file.read_text()).get("steps", [])}
        for s in steps:
            if s["step_id"] in manifest_steps:
                s.setdefault("stage", manifest_steps[s["step_id"]].get("stage", ""))

    done    = sum(1 for s in steps if s["status"] == "done")
    failed  = sum(1 for s in steps if s["status"] == "failed")
    total   = len(steps)
    dry_run = raw.get("dry_run", False)

    header = f"[bold cyan]{workspace}[/bold cyan]"
    if dry_run:
        header += "  [dim](dry-run)[/dim]"
    if raw.get("is_complete"):
        header += "  [green]complete[/green]"
    elif raw.get("was_interrupted"):
        header += "  [yellow]interrupted[/yellow]"

    app.print(f"\n{header}\n")

    table = Table(box=box.SIMPLE, show_header=True, header_style="bold")
    table.add_column("",        width=2)
    table.add_column("Step",    style="white")
    table.add_column("Stage",   style="cyan", no_wrap=True)
    table.add_column("Status",  no_wrap=True)
    if verbose:
        table.add_column("Time", style="dim", justify="right")
        table.add_column("Note", style="dim")

    for s in steps:
        status_val = s["status"]
        color, icon = _STATUS_STYLE.get(status_val, ("dim", "?"))
        status_text = Text(f"{icon} {status_val}", style=color)

        if verbose:
            elapsed = f"{s.get('elapsed_s', 0):.1f}s" if s.get("elapsed_s") else ""
            note    = s.get("error_message", "")[:60] if status_val == "failed" else ""
            table.add_row("", s["step_id"], s.get("stage", ""), status_text, elapsed, note)
        else:
            table.add_row("", s["step_id"], s.get("stage", ""), status_text)

    app.print(table)

    # Summary line
    pct = int(done / total * 100) if total else 0
    bar_filled = int(pct / 5)
    bar = "█" * bar_filled + "░" * (20 - bar_filled)
    app.print(f"  [{bar}] {done}/{total} steps  ({pct}%)")
    if failed:
        app.print(f"\n  [red]{failed} step(s) failed.[/red]  Run [dim]simforge status {workspace} --verbose[/dim] for details.\n")
    else:
        app.print()


# ═══════════════════════════════════════════════════════════════════════════════
# init
# ═══════════════════════════════════════════════════════════════════════════════

def _ask(prompt: str, options: list[tuple[str, str]], default: int = 1) -> str:
    """Numbered menu prompt. Returns the key of the selected option."""
    app.print(f"\n  [bold]{prompt}[/bold]")
    for i, (key, label) in enumerate(options, 1):
        marker = "[green]●[/green]" if i == default else " "
        app.print(f"    [{i}] {marker} {label}")
    while True:
        raw = input(f"\n  Select [1-{len(options)}] (Enter = {default}): ").strip()
        if raw == "":
            return options[default - 1][0]
        if raw.isdigit() and 1 <= int(raw) <= len(options):
            return options[int(raw) - 1][0]
        app.print("  [red]Opción inválida.[/red]")


def _ask_multi(prompt: str, options: list[tuple[str, str]], defaults: list[int]) -> list[str]:
    """Multi-select numbered menu. Returns list of selected keys."""
    app.print(f"\n  [bold]{prompt}[/bold]")
    app.print("  [dim](separados por coma, ej: 1,2,3)[/dim]")
    for i, (key, label) in enumerate(options, 1):
        marker = "[green]●[/green]" if i in defaults else " "
        app.print(f"    [{i}] {marker} {label}")
    default_str = ",".join(str(d) for d in defaults)
    while True:
        raw = input(f"\n  Select (Enter = {default_str}): ").strip()
        if raw == "":
            return [options[d - 1][0] for d in defaults]
        parts = [p.strip() for p in raw.split(",")]
        if all(p.isdigit() and 1 <= int(p) <= len(options) for p in parts):
            return [options[int(p) - 1][0] for p in parts]
        app.print("  [red]Entrada inválida.[/red]")


def _ask_float(prompt: str, default: float) -> float:
    """Free-text float prompt with default."""
    raw = input(f"\n  {prompt} (Enter = {default}): ").strip()
    if raw == "":
        return default
    try:
        return float(raw)
    except ValueError:
        app.print("  [yellow]Valor inválido, usando default.[/yellow]")
        return default


def _ask_str(prompt: str, default: str = "") -> str:
    """Free-text string prompt."""
    suffix = f" (Enter = {default})" if default else ""
    raw = input(f"\n  {prompt}{suffix}: ").strip()
    return raw if raw else default


@cli.command()
def init(
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Output YAML path (default: configs/<name>.yaml)."),
):
    """Interactive wizard — genera un YAML de configuración para SimForge."""

    app.print(Panel(
        "[bold cyan]SimForge Init[/bold cyan]  [dim]generador de configuración[/dim]",
        border_style="cyan",
        padding=(0, 2),
    ))

    # ── 1. Nombre del proyecto ────────────────────────────────────────────────
    name = _ask_str("Nombre del proyecto")
    if not name:
        app.print("[red]El nombre es obligatorio.[/red]")
        raise typer.Exit(1)
    name = name.replace(" ", "_")

    # ── 2. Tipo de sistema ────────────────────────────────────────────────────
    system_type = _ask("Tipo de sistema", [
        ("protein_in_water",  "Proteína / péptido en agua"),
        ("protein_ligand",    "Proteína + ligando (inhibición, binding)"),
        ("protein_membrane",  "Proteína + membrana"),
    ])

    # ── 3. Componentes ────────────────────────────────────────────────────────
    app.print("\n  [bold]Archivos PDB[/bold]")
    protein_file = _ask_str("  Archivo PDB de la proteína/péptido")
    if not protein_file:
        app.print("[red]El PDB de proteína es obligatorio.[/red]")
        raise typer.Exit(1)

    extra_components: list[dict] = []

    if system_type == "protein_ligand":
        n_extra = int(_ask("¿Cuántos ligandos/sustratos?", [
            ("1", "1"), ("2", "2"), ("3", "3"),
        ]))
        for i in range(1, n_extra + 1):
            role = _ask(f"Rol del componente {i}", [
                ("competitive_ligand", "Ligando competitivo"),
                ("substrate",          "Sustrato"),
                ("ligand",             "Ligando (genérico)"),
                ("cofactor",           "Cofactor"),
            ])
            pdb = _ask_str(f"  Archivo PDB del {role}_{i}")
            extra_components.append({"role": role, "file": pdb, "index": i})

    # ── 4. Forcefield ─────────────────────────────────────────────────────────
    ff_protein = _ask("Forcefield (proteína)", [
        ("opls-aa",   "OPLS-AA  (péptidos, proteínas pequeñas)"),
        ("charmm36",  "CHARMM36 (proteínas, membrana)"),
        ("amber99sb", "AMBER99SB-ILDN"),
    ])

    ff_ligands = None
    if system_type == "protein_ligand":
        ff_ligands = _ask("Forcefield (ligandos)", [
            ("gaff2",  "GAFF2   (AMBER ecosystem)"),
            ("cgenff", "CGenFF  (CHARMM ecosystem)"),
            ("opls",   "OPLS-AA (ligandos simples)"),
        ])

    # ── 5. Modelo de agua ─────────────────────────────────────────────────────
    water_model = _ask("Modelo de agua", [
        ("spce",   "SPC/E   (recomendado para OPLS-AA)"),
        ("tip3p",  "TIP3P   (recomendado para CHARMM/AMBER)"),
        ("tip4p",  "TIP4P-Ew"),
    ])

    # ── 6. Temperatura ────────────────────────────────────────────────────────
    temp_choice = _ask("Temperatura", [
        ("300",    "300 K  (estándar)"),
        ("309.65", "309.65 K  (37 °C — fisiológico humano)"),
        ("custom", "Otra..."),
    ])
    if temp_choice == "custom":
        temperature = _ask_float("Temperatura (K)", 300.0)
    else:
        temperature = float(temp_choice)

    # ── 7. Duración ───────────────────────────────────────────────────────────
    dur_choice = _ask("Duración de la simulación de producción", [
        ("1",    " 1 ns  (exploración rápida)"),
        ("10",   "10 ns  (estabilidad)"),
        ("50",   "50 ns  (binding, inhibición)"),
        ("100",  "100 ns"),
        ("custom", "Otra..."),
    ])
    if dur_choice == "custom":
        duration = _ask_float("Duración (ns)", 10.0)
    else:
        duration = float(dur_choice)

    # ── 8. Análisis ───────────────────────────────────────────────────────────
    all_analyses = [
        ("rmsd",             "RMSD  (estabilidad global)"),
        ("rmsf",             "RMSF  (flexibilidad por residuo)"),
        ("energy",           "Energía  (Ep, Ek, T, P)"),
        ("hydrogen_bonds",   "Puentes de hidrógeno"),
        ("radius_of_gyration", "Radio de giro"),
        ("distance_analysis", "Distancias entre grupos"),
    ]
    default_analyses = [1, 2, 3]
    if system_type == "protein_ligand":
        default_analyses = [1, 2, 3, 4, 6]

    selected_analyses = _ask_multi("Análisis a incluir", all_analyses, default_analyses)

    # ── 9. Hardware ───────────────────────────────────────────────────────────
    hardware = _ask("Hardware disponible para mdrun", [
        ("auto", "Auto-detectar en runtime (recomendado)"),
        ("gpu",  "GPU  (NVIDIA, siempre usar flags GPU)"),
        ("cpu",  "CPU  (sin GPU)"),
    ])

    # ── 10. Objetivos de simulación ───────────────────────────────────────────
    objectives_map = {
        "protein_in_water": ["stability"],
        "protein_ligand":   ["competitive_binding", "active_site_stability"],
        "protein_membrane": ["stability", "membrane_protein_dynamics"],
    }
    objectives = objectives_map[system_type]

    # ── Construir YAML ────────────────────────────────────────────────────────
    components_yaml = f"  - id: protein_1\n    role: protein\n    file: {protein_file}\n"
    for ec in extra_components:
        role_id = f"{ec['role']}_{ec['index']}"
        components_yaml += (
            f"\n  - id: {role_id}\n"
            f"    role: {ec['role']}\n"
            f"    file: {ec['file']}\n"
        )

    forcefields_yaml = f"  protein: {ff_protein}\n"
    if ff_ligands:
        forcefields_yaml += f"  ligands: {ff_ligands}\n"

    objectives_yaml = "\n".join(f"  - {o}" for o in objectives)

    analyses_yaml = ""
    for a in selected_analyses:
        analyses_yaml += f"\n  - type: {a}"
        if a == "distance_analysis" and system_type == "protein_ligand" and extra_components:
            first = extra_components[0]
            second = extra_components[1] if len(extra_components) > 1 else extra_components[0]
            g1 = f"{first['role']}_{first['index']}"
            g2 = f"{second['role']}_{second['index']}"
            analyses_yaml += f"\n    selection:\n      group1: {g1}\n      group2: {g2}"

    hardware_line = f"\nhardware: {hardware}\n" if hardware != "auto" else ""

    yaml_content = f"""project:
  name: {name}

components:
{components_yaml}
environment:
  solvent:
    water_model: {water_model}
  ions:
    concentration: 0.154
  temperature_K: {temperature}
  duration_ns: {duration}

forcefields:
{forcefields_yaml}
simulation_objectives:
{objectives_yaml}

analysis:{analyses_yaml}
{hardware_line}"""

    # ── Escribir archivo ──────────────────────────────────────────────────────
    if output is None:
        configs_dir = Path("configs")
        configs_dir.mkdir(exist_ok=True)
        output = configs_dir / f"{name}.yaml"

    output.write_text(yaml_content)

    app.print(f"\n[bold green]✓[/bold green] Config generado → [cyan]{output}[/cyan]")
    app.print(f"\nSiguiente:  [dim]simforge compile {output}[/dim]\n")


# ═══════════════════════════════════════════════════════════════════════════════
# Scientific planning dialogue
# ═══════════════════════════════════════════════════════════════════════════════

_KIND_LABELS = {
    "missing_parameter":   ("?",  "dim"),
    "policy_selection":    ("◆",  "cyan"),
    "risk_confirmation":   ("⚠",  "yellow"),
    "scientific_tradeoff": ("⟁",  "magenta"),
    "protocol_selection":  ("◎",  "blue"),
}

_KIND_TITLES = {
    "missing_parameter":   "Missing Parameter",
    "policy_selection":    "Policy Selection",
    "risk_confirmation":   "Risk Confirmation",
    "scientific_tradeoff": "Scientific Tradeoff",
    "protocol_selection":  "Protocol Selection",
}


def _run_planning_dialogue(
    state,
    yaml_path: Path,
    skip: bool = False,
) -> Optional[object]:
    """
    Detecta preguntas de planning, muestra el diálogo y aplica patches al state.

    Retorna PlanningSessionRecord, o None si el usuario eligió abortar.
    """
    from core.scientific_planner import detect_questions, apply_patches
    from core.planning_models import PlanningSessionRecord, PlanningAnswer
    from datetime import datetime

    questions = detect_questions(state)

    record = PlanningSessionRecord(
        created_at  = datetime.now().isoformat(timespec="seconds"),
        config_path = str(yaml_path),
        skipped     = skip or len(questions) == 0,
    )

    if not questions or skip:
        if questions and skip:
            app.print(f"[dim]--no-plan: skipping {len(questions)} planning question(s)[/dim]\n")
        return record

    app.print(f"[bold]Scientific Planning[/bold]  [dim]{len(questions)} question(s) detected[/dim]\n")

    answers: list[PlanningAnswer] = []
    total = len(questions)

    for idx, q in enumerate(questions, 1):
        icon, color = _KIND_LABELS.get(q.kind.value, ("?", "dim"))
        title = _KIND_TITLES.get(q.kind.value, q.kind.value)

        # ── Header ───────────────────────────────────────────────────────────
        app.rule(f"[{color}]{icon} [{idx}/{total}] {title}[/{color}]")
        app.print()

        # ── Context ───────────────────────────────────────────────────────────
        for line in q.context.splitlines():
            app.print(f"  {line}")
        app.print()

        # ── Question ──────────────────────────────────────────────────────────
        app.print(f"  [bold]{q.question}[/bold]")
        app.print()

        # ── Options ───────────────────────────────────────────────────────────
        for i, opt in enumerate(q.options, 1):
            marker = "[bold green]●[/bold green]" if opt.is_default else " "
            app.print(f"  [{i}] {marker} [white]{opt.label}[/white]")
            app.print(f"       [dim]{opt.description}[/dim]")
        app.print()

        # ── Prompt ────────────────────────────────────────────────────────────
        default_idx = next(
            (i for i, o in enumerate(q.options, 1) if o.is_default), 1
        )

        while True:
            raw = input(f"  Select [1-{len(q.options)}] (Enter = {default_idx}): ").strip()
            if raw == "":
                choice = default_idx
                break
            if raw.isdigit() and 1 <= int(raw) <= len(q.options):
                choice = int(raw)
                break
            app.print(f"  [red]Invalid choice.[/red] Enter a number between 1 and {len(q.options)}.")

        selected = q.options[choice - 1]
        app.print(f"\n  [green]→[/green] {selected.label}\n")

        # ── Abort check ───────────────────────────────────────────────────────
        if selected.state_patch.get("abort"):
            app.print("[yellow]Compile aborted by planning decision.[/yellow]")
            app.print(f"  Reason: {selected.description}\n")
            return None

        answer = PlanningAnswer(
            question_id    = q.id,
            selected_key   = selected.key,
            selected_label = selected.label,
            state_patch    = {k: v for k, v in selected.state_patch.items() if k != "abort"},
        )
        answers.append(answer)

    # ── Apply patches ─────────────────────────────────────────────────────────
    all_patches: list[dict] = [a.state_patch for a in answers if a.state_patch]
    apply_patches(state, answers)

    # ── Show summary ──────────────────────────────────────────────────────────
    if all_patches:
        app.rule("[dim]Planning complete[/dim]")
        app.print()
        for a in answers:
            if a.state_patch:
                app.print(f"  [green]✓[/green] {a.selected_label}")
        app.print()

    # ── Build record ──────────────────────────────────────────────────────────
    record.questions_shown = [q.model_dump(mode="json") for q in questions]
    record.answers         = answers
    record.patches_applied = all_patches

    return record


def _save_planning_session(record, path: Path) -> None:
    """Persiste la sesión de planning en metadata/ del workspace."""
    if record is None:
        return
    try:
        path.write_text(
            json.dumps(record.model_dump(mode="json"), indent=4)
        )
    except Exception:
        pass   # no bloquear compile por fallo de persistencia


# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    cli()
