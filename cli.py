#!/usr/bin/env python3
# cli.py — SimForge command-line interface

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.table import Table
from rich import box
from rich.panel import Panel
from rich.text import Text

app     = Console()
cli     = typer.Typer(
    name="simforge",
    help="Molecular simulation workflow compiler.",
    add_completion=False,
)


# ═══════════════════════════════════════════════════════════════════════════════
# compile
# ═══════════════════════════════════════════════════════════════════════════════

@cli.command()
def compile(
    yaml_path: Path = typer.Argument(..., help="YAML config file to compile."),
    output_dir: str = typer.Option("simforge_runs", "--output-dir", "-o", help="Root output directory."),
    no_build: bool  = typer.Option(False, "--no-build", help="Only compile — skip workspace materialization."),
):
    """Compile a YAML config into an executable workspace."""

    from core.compiler import SimulationCompiler
    from builders.workspace_builder import WorkspaceBuilder

    if not yaml_path.exists():
        app.print(f"[red]Error:[/red] File not found: {yaml_path}")
        raise typer.Exit(1)

    app.print(f"\n[bold cyan]SimForge[/bold cyan]  compiling [green]{yaml_path}[/green]\n")

    # ── Compile ───────────────────────────────────────────────────────────────
    try:
        result = SimulationCompiler().compile(str(yaml_path))
    except Exception as e:
        app.print(f"[red]Compilation failed:[/red] {e}")
        raise typer.Exit(1)

    plan  = result.plan
    state = result.state

    # ── Summary panel ─────────────────────────────────────────────────────────
    lines = [
        f"  System type    {state.inferred_system_type}",
        f"  Steps          {len(plan.steps)}",
        f"  Policy         {plan.workflow_policy.production_time_ns}ns production  "
        f"T={plan.workflow_policy.temperature_K}K  "
        f"P={plan.workflow_policy.pressure_bar}bar",
    ]
    if plan.workflow_policy.enhanced_sampling:
        lines.append(f"  Sampling       {plan.workflow_policy.sampling_method}")
    if plan.blocking_issues:
        lines.append(f"  [yellow]Blocking issues  {len(plan.blocking_issues)}[/yellow]")
    if plan.special_protocols:
        lines.append(f"  Protocols      {', '.join(plan.special_protocols)}")

    app.print(Panel("\n".join(lines), title="Compilation result", border_style="cyan"))

    # ── Workflow table ────────────────────────────────────────────────────────
    table = Table(box=box.SIMPLE, show_header=True, header_style="bold")
    table.add_column("#",       style="dim",   width=3)
    table.add_column("Step",    style="white")
    table.add_column("Stage",   style="cyan")
    table.add_column("Engine",  style="dim")
    table.add_column("Type",    style="dim")

    for i, step in enumerate(result.execution_order, 1):
        table.add_row(
            str(i),
            step.step_id,
            step.stage.value,
            step.engine,
            step.step_type.value,
        )

    app.print(table)

    # ── Blocking issues ───────────────────────────────────────────────────────
    if plan.blocking_issues:
        app.print("[yellow]Blocking issues:[/yellow]")
        for issue in plan.blocking_issues:
            app.print(f"  [{issue.severity.value}] {issue.source}: {issue.message}")
        app.print()

    # ── Materialize ───────────────────────────────────────────────────────────
    if not no_build:
        app.print("[dim]Materializing workspace...[/dim]")
        try:
            workspace = WorkspaceBuilder().build(result, output_dir=output_dir)
        except Exception as e:
            app.print(f"[red]Workspace build failed:[/red] {e}")
            raise typer.Exit(1)

        app.print(f"\n[bold green]✓[/bold green] Workspace ready → [cyan]{workspace}[/cyan]")
        app.print(f"  Manifest → {workspace}/metadata/execution_manifest.json")
        app.print(f"\nNext:  [dim]simforge run {workspace}[/dim]\n")
    else:
        app.print("[dim](--no-build: workspace not materialized)[/dim]\n")


# ═══════════════════════════════════════════════════════════════════════════════
# run
# ═══════════════════════════════════════════════════════════════════════════════

@cli.command()
def run(
    workspace: Path = typer.Argument(..., help="Workspace directory to execute."),
    dry_run: bool   = typer.Option(True,  "--dry-run/--real", help="Dry-run (default) or real execution."),
    executor_type: str = typer.Option("shell", "--executor", "-e", help="Executor type: shell | gromacs."),
):
    """Execute a compiled workspace."""

    if not workspace.exists():
        app.print(f"[red]Error:[/red] Workspace not found: {workspace}")
        raise typer.Exit(1)

    manifest_file = workspace / "metadata" / "execution_manifest.json"
    if not manifest_file.exists():
        app.print(f"[red]Error:[/red] No execution manifest found. Run [cyan]simforge compile[/cyan] first.")
        raise typer.Exit(1)

    mode_label = "[yellow]DRY-RUN[/yellow]" if dry_run else "[bold red]REAL[/bold red]"
    app.print(f"\n[bold cyan]SimForge[/bold cyan]  {mode_label}  [green]{workspace}[/green]\n")

    if not dry_run:
        confirmed = typer.confirm("Real execution will run GROMACS commands. Continue?")
        if not confirmed:
            raise typer.Exit(0)

    # ── Select executor ───────────────────────────────────────────────────────
    if executor_type == "gromacs":
        from executors.gromacs_executor import GROMACSExecutor
        executor = GROMACSExecutor(workspace, dry_run=dry_run)
    else:
        from executors.shell_executor import ShellExecutor
        executor = ShellExecutor(workspace, dry_run=dry_run)

    # ── Run ───────────────────────────────────────────────────────────────────
    try:
        state = executor.run()
    except Exception as e:
        app.print(f"[red]Execution error:[/red] {e}")
        raise typer.Exit(1)

    # ── Result summary ────────────────────────────────────────────────────────
    done    = sum(1 for s in state.steps if s.status.value == "done")
    failed  = sum(1 for s in state.steps if s.status.value == "failed")
    skipped = sum(1 for s in state.steps if s.status.value in ("skipped", "blocked"))

    color = "green" if failed == 0 else "red"
    app.print(f"\n[{color}]{'✓' if failed == 0 else '✗'}[/{color}]  "
              f"Done: {done}  Failed: {failed}  Skipped: {skipped}  "
              f"Complete: {'yes' if state.is_complete else 'no'}\n")

    if failed:
        app.print("[red]Failed steps:[/red]")
        for s in state.steps:
            if s.status.value == "failed":
                app.print(f"  {s.step_id}  exit={s.exit_code}")
        raise typer.Exit(1)


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

if __name__ == "__main__":
    cli()
