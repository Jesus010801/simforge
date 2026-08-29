"""Thin Typer entry points for ``simforge study inspect|analyze|analyses``.

All scientific logic lives in ``analysis.campaign`` — these functions only parse
input, call the orchestration layer, and render results.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Annotated, Optional

import typer
from rich.console import Console

from analysis.campaign.models import AnalysisStatus, Severity, ValidationState

_console = Console()


# ═══════════════════════════════════════════════════════════════════════════════
# simforge study analyses
# ═══════════════════════════════════════════════════════════════════════════════

def study_analyses_fn(
    as_json: Annotated[bool, typer.Option("--json", help="Emit the catalog as JSON.")] = False,
) -> None:
    """List every registered analysis (id, name, category, required components)."""
    from analysis.campaign.observables import registry
    registry.ensure_loaded()
    catalog = registry.catalog()
    if as_json:
        print(json.dumps(catalog, indent=2))
        raise typer.Exit(0)
    _console.print("[bold]Registered study analyses[/bold]\n")
    for spec in catalog:
        req = ", ".join(spec["required_components"]) or "—"
        _console.print(f"  [cyan]{spec['id']}[/cyan]  [dim]({spec['category']})[/dim]")
        _console.print(f"      {spec['display_name']}")
        _console.print(f"      {spec['description']}")
        _console.print(f"      [dim]requires: {req}[/dim]\n")
    _console.print(
        "[dim]Run:  simforge study analyze <dir> --analysis <id> [--analysis <id> ...][/dim]"
    )


# ═══════════════════════════════════════════════════════════════════════════════
# simforge study inspect
# ═══════════════════════════════════════════════════════════════════════════════

def study_inspect_fn(
    path: Annotated[Path, typer.Argument(help="Study directory to inspect.")] = Path("."),
    output: Annotated[Optional[str], typer.Option("--output", "-o",
        help="Output directory (default: <study>/simforge_analysis).")] = None,
    manifest: Annotated[Optional[str], typer.Option("--manifest",
        help="Use an existing study_manifest.yaml instead of rediscovering.")] = None,
    no_trajectory_inspection: Annotated[bool, typer.Option("--no-trajectory-inspection",
        help="Skip `gmx check` (faster; no frame/box metadata).")] = False,
    strong_fingerprint: Annotated[bool, typer.Option("--strong-fingerprint",
        help="Full-content hash of every source file (slow for large trajectories).")] = False,
    gmx: Annotated[str, typer.Option("--gmx", help="GROMACS binary.")] = "gmx",
    as_json: Annotated[bool, typer.Option("--json", help="Emit machine-readable JSON.")] = False,
) -> None:
    """Discover a study, build its manifest, infer components, inspect trajectories
    and validate — without running any scientific analysis."""
    from analysis.campaign.orchestration.study_analyzer import run_inspect

    if not path.is_dir():
        _console.print(f"[red]Error:[/red] not a directory: {path}")
        raise typer.Exit(1)

    result = run_inspect(
        path, output_dir=output, manifest_path=manifest, gmx=gmx,
        strong_fingerprint=strong_fingerprint,
        inspect_trajectories=not no_trajectory_inspection,
    )

    if as_json:
        print(json.dumps(result.to_dict(), indent=2))
        raise typer.Exit(0 if result.validation.study_state != ValidationState.INVALID else 1)

    _render_inspect(result)
    raise typer.Exit(0 if result.validation.study_state != ValidationState.INVALID else 1)


# ═══════════════════════════════════════════════════════════════════════════════
# simforge study analyze
# ═══════════════════════════════════════════════════════════════════════════════

def study_analyze_fn(
    path: Annotated[Path, typer.Argument(help="Study directory to analyse.")] = Path("."),
    analysis: Annotated[Optional[list[str]], typer.Option("--analysis",
        help="Analysis id to run (repeatable). NOTHING runs unless you pass this "
             "or select interactively.")] = None,
    output: Annotated[Optional[str], typer.Option("--output", "-o",
        help="Output directory (default: <study>/simforge_analysis).")] = None,
    manifest: Annotated[Optional[str], typer.Option("--manifest",
        help="Use an existing (optionally hand-corrected) study_manifest.yaml.")] = None,
    non_interactive: Annotated[bool, typer.Option("--non-interactive",
        help="Never prompt; require --analysis.")] = False,
    dry_run: Annotated[bool, typer.Option("--dry-run",
        help="Validate + construct commands but do not execute GROMACS.")] = False,
    force: Annotated[bool, typer.Option("--force",
        help="Ignore cached trajectory views / results.")] = False,
    gmx: Annotated[str, typer.Option("--gmx", help="GROMACS binary.")] = "gmx",
    as_json: Annotated[bool, typer.Option("--json", help="Emit machine-readable JSON.")] = False,
) -> None:
    """Discover + validate a study, then run ONLY the analyses you explicitly select.

    No analysis is ever selected by default.
    """
    from analysis.campaign.observables import registry
    from analysis.campaign.orchestration.study_analyzer import run_analyze
    registry.ensure_loaded()

    if not path.is_dir():
        _console.print(f"[red]Error:[/red] not a directory: {path}")
        raise typer.Exit(1)

    requested = list(analysis or [])
    interactive = sys.stdin.isatty() and not non_interactive and not as_json

    if not requested:
        if interactive:
            requested = _interactive_select()
            if requested is None:
                _console.print("[yellow]Cancelled — no analysis run.[/yellow]")
                raise typer.Exit(0)
        else:
            _console.print(
                "[yellow]No analysis selected.[/yellow] Nothing was run.\n"
                "  Select analyses explicitly, e.g.:\n"
                "    simforge study analyze . --analysis rmsd-receptor --analysis rmsd-complex\n"
                "  List available analyses:\n"
                "    simforge study analyses"
            )
            raise typer.Exit(2)

    result = run_analyze(
        path, requested, output_dir=output, manifest_path=manifest, gmx=gmx,
        dry_run=dry_run, force=force,
    )

    if as_json:
        print(json.dumps(result.to_dict(), indent=2))
    else:
        _render_inspect(result)
        _render_analysis_results(result)

    errs = result.error_warnings()
    if errs and not result.results:
        for w in errs:
            _console.print(f"[red]Error[/red] ({w.code}): {w.message}", highlight=False)
        raise typer.Exit(1)
    any_success = any(r.status in (AnalysisStatus.SUCCESS, AnalysisStatus.PLANNED,
                                   AnalysisStatus.CACHED) for r in result.results)
    raise typer.Exit(0 if any_success or not result.results else 1)


# ═══════════════════════════════════════════════════════════════════════════════
# rendering helpers
# ═══════════════════════════════════════════════════════════════════════════════

def _render_inspect(result) -> None:
    m = result.manifest
    _console.print(f"\n[bold cyan]SimForge Study[/bold cyan]  [dim]{result.study_root}[/dim]")
    _console.print(f"[dim]output: {result.output_dir}[/dim]\n")
    for line in (result.validation.lines if result.validation else []):
        _console.print(line)
    if m and m.ambiguities:
        _console.print("\n[yellow]Ambiguities (edit study_manifest.yaml → 'resolution'):[/yellow]")
        for a in m.ambiguities:
            _console.print(f"  [yellow]•[/yellow] {a.system_id} — {a.message}")
            if a.options:
                _console.print(f"      options: {a.options}")
    warns = [w for w in (result.warnings or []) if w.severity in (Severity.WARN, Severity.ERROR)]
    if warns:
        _console.print("\n[yellow]Warnings:[/yellow]")
        for w in warns[:25]:
            colour = "red" if w.severity == Severity.ERROR else "yellow"
            _console.print(f"  [{colour}]{w.severity}[/{colour}] [{w.code}] {w.message}")


def _render_analysis_results(result) -> None:
    if not result.results:
        return
    _console.print("\n[bold]Analysis results[/bold]")
    _STYLE = {
        AnalysisStatus.SUCCESS: "green", AnalysisStatus.CACHED: "green",
        AnalysisStatus.PLANNED: "cyan", AnalysisStatus.SKIPPED: "dim",
        AnalysisStatus.REVIEW_REQUIRED: "yellow", AnalysisStatus.FAILED: "red",
        AnalysisStatus.ERROR: "red",
    }
    for r in result.results:
        colour = _STYLE.get(r.status, "white")
        _console.print(
            f"  [{colour}]{r.status:<16}[/{colour}] {r.system_id}  [dim]{r.analysis_id}[/dim]")
        if r.status == AnalysisStatus.SUCCESS and r.data_summary:
            ds = r.data_summary
            if "mean_nm" in ds:
                _console.print(
                    f"      mean={ds['mean_nm']} nm  final={ds.get('final_nm')} nm  "
                    f"plateau={ds.get('plateau_mean_nm')} nm  ({ds.get('n_frames')} frames)")
        if r.status in (AnalysisStatus.SKIPPED, AnalysisStatus.REVIEW_REQUIRED,
                        AnalysisStatus.FAILED) and r.message:
            _console.print(f"      [dim]{r.message}[/dim]")


def _interactive_select() -> Optional[list[str]]:
    from analysis.campaign.observables import registry
    specs = registry.catalog()
    _console.print("[bold]Select analyses[/bold] (none selected by default)\n")
    for i, s in enumerate(specs, 1):
        _console.print(f"  [{i}] {s['id']}  [dim]{s['display_name']}[/dim]")
    _console.print("\n  Enter numbers separated by commas (e.g. 1,3), 'all', or blank to cancel.")
    try:
        raw = input("> ").strip()
    except (EOFError, KeyboardInterrupt):
        return None
    if not raw:
        return None
    if raw.lower() == "all":
        return [s["id"] for s in specs]
    chosen: list[str] = []
    for tok in raw.replace(" ", "").split(","):
        if not tok:
            continue
        if not tok.isdigit() or not (1 <= int(tok) <= len(specs)):
            _console.print(f"[red]Ignoring invalid selection '{tok}'[/red]")
            continue
        chosen.append(specs[int(tok) - 1]["id"])
    return chosen or None
