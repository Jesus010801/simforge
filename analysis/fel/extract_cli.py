"""CLI handler for `simforge fel extract-minima`."""
from __future__ import annotations

from pathlib import Path
from typing import Annotated, Optional

import typer

from analysis.fel.models import ExtractConfig
from analysis.fel.extract import run_extract_minima


def fel_extract_minima_fn(
    fel_dir: Annotated[
        str,
        typer.Argument(help="Path to an existing FEL output directory."),
    ],
    trajectory: Annotated[
        str,
        typer.Option("--trajectory", "-f", help="GROMACS trajectory file (.xtc/.trr)."),
    ],
    topology: Annotated[
        str,
        typer.Option("--topology", "-s", help="GROMACS topology/run input (.tpr)."),
    ],
    structure: Annotated[
        Optional[str],
        typer.Option("--structure", help="Reference structure file (.gro/.pdb)."),
    ] = None,
    index: Annotated[
        Optional[str],
        typer.Option("--index", "-n", help="GROMACS index file (.ndx)."),
    ] = None,
    selection: Annotated[
        str,
        typer.Option("--selection", help="Group name passed to gmx trjconv (e.g. Protein)."),
    ] = "Protein",
    n_minima: Annotated[
        int,
        typer.Option("--n-minima", help="Maximum number of minima to extract."),
    ] = 5,
    min_distance_bins: Annotated[
        int,
        typer.Option("--min-distance-bins", help="Minimum bin-space separation between selected minima."),
    ] = 2,
    representative: Annotated[
        str,
        typer.Option("--representative", help="Strategy for representative frame selection."),
    ] = "nearest-frame",
    output_format: Annotated[
        str,
        typer.Option("--format", help="Output structure format: pdb, gro, or both."),
    ] = "pdb",
    gmx: Annotated[
        str,
        typer.Option("--gmx", help="GROMACS binary name or path."),
    ] = "gmx",
    out: Annotated[
        str,
        typer.Option("--out", "-o", help="Output directory for representative states."),
    ] = "representative_states",
    dry_run: Annotated[
        bool,
        typer.Option("--dry-run", help="Plan extraction without calling gmx."),
    ] = False,
    execute: Annotated[
        bool,
        typer.Option("--execute", help="Run gmx trjconv to extract structures."),
    ] = False,
    force: Annotated[
        bool,
        typer.Option("--force", help="Overwrite existing output directory."),
    ] = False,
    max_delta_g_kj: Annotated[
        str,
        typer.Option(
            "--max-delta-g-kj",
            help=(
                "Keep only minima with ΔG ≤ this value (kJ/mol). "
                "Use 'none' to disable the filter. Default: 2.5"
            ),
        ),
    ] = "2.5",
    min_count: Annotated[
        int,
        typer.Option(
            "--min-count",
            help="Exclude bins with fewer than this many trajectory frames. Default: 1",
        ),
    ] = 1,
    system_name: Annotated[
        Optional[str],
        typer.Option(
            "--system-name",
            help=(
                "System label used in output filenames, e.g. 'HMGCoA_R'. "
                "Inferred from FEL dir / topology / trajectory stem when omitted."
            ),
        ),
    ] = None,
    compat_receptor_name: Annotated[
        bool,
        typer.Option(
            "--compat-receptor-name",
            help="Also create receptor.{ext} symlink alongside the descriptive filename.",
        ),
    ] = False,
) -> None:
    """Detect free-energy minima and map them to representative trajectory frames.

    \b
    Dry-run (plan only):
        simforge fel extract-minima FEL_DIR --trajectory md.xtc --topology md.tpr --dry-run

    \b
    Execute (extract structures):
        simforge fel extract-minima FEL_DIR --trajectory md.xtc --topology md.tpr --execute

    \b
    Disable energy filter:
        simforge fel extract-minima FEL_DIR ... --max-delta-g-kj none --dry-run
    """
    if not dry_run and not execute:
        typer.echo(
            "Error: specify --dry-run (plan only) or --execute (run gmx).",
            err=True,
        )
        raise typer.Exit(1)

    if output_format not in ("pdb", "gro", "both"):
        typer.echo(
            f"Error: --format must be 'pdb', 'gro', or 'both' (got '{output_format}').",
            err=True,
        )
        raise typer.Exit(1)

    # Parse --max-delta-g-kj
    try:
        max_dg: Optional[float] = (
            None if max_delta_g_kj.strip().lower() == "none"
            else float(max_delta_g_kj)
        )
    except ValueError:
        typer.echo(
            f"Error: --max-delta-g-kj must be a number or 'none' (got '{max_delta_g_kj}').",
            err=True,
        )
        raise typer.Exit(1)

    fel_path = Path(fel_dir)
    if not fel_path.is_dir():
        typer.echo(f"Error: FEL directory not found: {fel_dir}", err=True)
        raise typer.Exit(1)

    out_path = Path(out)
    if out_path.exists() and not force:
        typer.echo(
            f"Output directory already exists: {out}. "
            "Use --force to overwrite or choose a new --out directory.",
            err=True,
        )
        raise typer.Exit(1)

    config = ExtractConfig(
        fel_dir=fel_path,
        trajectory=Path(trajectory),
        topology=Path(topology),
        structure=Path(structure) if structure else None,
        index=Path(index) if index else None,
        selection=selection,
        n_minima=n_minima,
        min_distance_bins=min_distance_bins,
        representative=representative,
        format=output_format,
        gmx=gmx,
        out_dir=out_path,
        dry_run=dry_run,
        execute=execute,
        force=force,
        max_delta_g_kj=max_dg,
        min_count=min_count,
        system_name=system_name,
        compat_receptor_name=compat_receptor_name,
    )

    result = run_extract_minima(config)

    # ── Report outcome ────────────────────────────────────────────────────────
    error_warns = [w for w in result.warnings if w.severity == "error"]
    if error_warns:
        for w in error_warns:
            typer.echo(f"Error [{w.code}]: {w.message}", err=True)
        raise typer.Exit(1)

    mode_label = "DRY RUN" if dry_run else "completed"
    typer.echo(f"extract-minima {mode_label} → {out}/")
    typer.echo(f"  Minima found: {len(result.minima)}")
    if result.rejected_minima:
        typer.echo(f"  Rejected (energy filter): {len(result.rejected_minima)}")
    typer.echo(f"  Report: {out}/report.md")

    other_warns = [w for w in result.warnings if w.severity != "error"]
    for w in other_warns:
        icon = "⚠" if w.severity == "warn" else "ℹ"
        typer.echo(f"  {icon} [{w.code}] {w.message}")
