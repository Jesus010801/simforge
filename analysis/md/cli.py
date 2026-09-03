"""CLI implementation for `simforge analyze-md`.

Exposes `analyze_md_fn` — a plain function with Typer annotations — for
registration in the root cli.py:

    from analysis.md.cli import analyze_md_fn
    cli.command(name="analyze-md")(analyze_md_fn)

This keeps cli.py changes to two lines while placing all implementation here.

Naming trade-off: the existing `simforge analyze` command (XVG quality
classification) is registered as @cli.command(). Typer does not support the
same name as both a plain command and a sub-app callback, so we use the name
"analyze-md" to avoid breaking backward compatibility. A future unified
`simforge analyze xvg / analyze md` hierarchy can be introduced in a separate
refactor that retires the old flat command.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated, Optional

import typer
from rich.console import Console

from analysis.md.discovery import discover_md_run
from analysis.md.observables import plan_observables
from analysis.md.provenance import build_provenance
from analysis.md.report import write_report

_console = Console()


def analyze_md_fn(
    run_dir: Annotated[str, typer.Argument(help="Path to the MD run directory")],
    out: Annotated[str, typer.Option("--out", "-o", help="Output directory")] = "analysis",
    dry_run: Annotated[bool, typer.Option("--dry-run", help="Plan only — no GROMACS commands executed")] = False,
    trajectory: Annotated[Optional[str], typer.Option("--trajectory", help="Override primary trajectory (.xtc/.trr)")] = None,
    topology: Annotated[Optional[str], typer.Option("--topology", help="Override topology input (.tpr)")] = None,
    structure: Annotated[Optional[str], typer.Option("--structure", help="Override reference structure (.gro/.pdb)")] = None,
    energy: Annotated[Optional[str], typer.Option("--energy", help="Override energy file (.edr)")] = None,
    index: Annotated[Optional[str], typer.Option("--index", help="Override index file (.ndx)")] = None,
    log: Annotated[Optional[str], typer.Option("--log", help="Override primary log file")] = None,
) -> None:
    """Discover GROMACS artifacts in RUN_DIR and plan analysis observables."""
    run_path = Path(run_dir)
    if not run_path.exists():
        _console.print(f"[red]Error:[/red] Run directory not found: {run_path}")
        raise typer.Exit(1)

    out_path = Path(out)

    overrides: dict[str, str] = {}
    if trajectory is not None:
        overrides["trajectory"] = trajectory
    if topology is not None:
        overrides["tpr"] = topology
    if structure is not None:
        overrides["structure"] = structure
    if energy is not None:
        overrides["edr"] = energy
    if index is not None:
        overrides["index"] = index
    if log is not None:
        overrides["log"] = log

    with _console.status(f"  Discovering MD run in [cyan]{run_path}[/cyan] ..."):
        run = discover_md_run(run_path, overrides=overrides or None)

    _console.print(f"  [green]✓[/green] Discovered {len(run.all_files)} GROMACS artifact(s)")

    for w in run.warnings:
        icon = "[yellow]⚠[/yellow]" if w.severity == "warn" else "[cyan]ℹ[/cyan]"
        _console.print(f"  {icon} {w.message}")

    plan       = plan_observables(run)
    provenance = build_provenance(run, out_path, dry_run)

    if dry_run:
        meta_dir = out_path / "metadata"
        meta_dir.mkdir(parents=True, exist_ok=True)

        (meta_dir / "run_discovery.json").write_text(
            json.dumps(run.to_dict(), indent=2, default=str)
        )
        (meta_dir / "analysis_plan.json").write_text(
            json.dumps(plan.to_dict(), indent=2)
        )
        (meta_dir / "provenance.json").write_text(
            json.dumps(provenance.to_dict(), indent=2)
        )
        write_report(run, plan, provenance, out_path)

        planned = sum(1 for o in plan.observables if o.status in ("planned", "needs_selection"))
        avail   = sum(1 for o in plan.observables if o.status == "already_available")
        unavail = sum(1 for o in plan.observables if o.status == "unavailable")

        _console.print(f"\n  [bold]Plan written to[/bold] [cyan]{out_path}/[/cyan]")
        _console.print(
            f"  Observables: {planned} planned, {avail} already available, {unavail} unavailable"
        )
        _console.print("\n  [dim]Dry-run complete — no GROMACS commands were executed.[/dim]")
    else:
        _console.print(
            "\n  [yellow]Live execution not implemented yet.[/yellow] "
            "Re-run with [bold]--dry-run[/bold] to generate a plan."
        )
