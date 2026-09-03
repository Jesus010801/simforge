"""CLI entry point for `simforge fel run`."""
from __future__ import annotations

from pathlib import Path
from typing import Annotated, Optional

import typer

from analysis.fel.models import FELConfig
from analysis.fel.workflow import run_fel
from analysis.fel.xvg import discover_feature_xvg, feature_label


def fel_run_fn(
    run_dir: Annotated[
        Optional[str],
        typer.Argument(help="MD run directory for feature auto-discovery mode."),
    ] = None,
    xvg_x: Annotated[
        Optional[str],
        typer.Option("--xvg-x", help="XVG file for the X-axis feature."),
    ] = None,
    xvg_y: Annotated[
        Optional[str],
        typer.Option("--xvg-y", help="XVG file for the Y-axis feature."),
    ] = None,
    x_column: Annotated[
        int,
        typer.Option("--x-column", help="1-based data column index for X axis."),
    ] = 1,
    y_column: Annotated[
        int,
        typer.Option("--y-column", help="1-based data column index for Y axis."),
    ] = 1,
    x_label: Annotated[
        Optional[str],
        typer.Option("--x-label", help="Display label for X axis (default: column name)."),
    ] = None,
    y_label: Annotated[
        Optional[str],
        typer.Option("--y-label", help="Display label for Y axis (default: column name)."),
    ] = None,
    temperature: Annotated[
        float,
        typer.Option("--temperature", help="Simulation temperature in Kelvin."),
    ] = 300.0,
    bins: Annotated[
        str,
        typer.Option("--bins", help='Number of histogram bins as "NX,NY".'),
    ] = "50,50",
    out: Annotated[
        str,
        typer.Option("--out", "-o", help="Output directory."),
    ] = "fel_analysis",
    time_unit_x: Annotated[
        str,
        typer.Option("--time-unit-x", help='Force X time unit: "ps", "ns", or "auto".'),
    ] = "auto",
    time_unit_y: Annotated[
        str,
        typer.Option("--time-unit-y", help='Force Y time unit: "ps", "ns", or "auto".'),
    ] = "auto",
    time_tolerance: Annotated[
        float,
        typer.Option("--time-tolerance", help="Max allowed time-axis diff (ns)."),
    ] = 1e-3,
    empty_bin_policy: Annotated[
        str,
        typer.Option("--empty-bins", help='"nan" (default) or "cap".'),
    ] = "nan",
    energy_cap_kj: Annotated[
        float,
        typer.Option("--energy-cap-kj", help="Energy cap for empty bins in kJ/mol."),
    ] = 50.0,
    features: Annotated[
        Optional[str],
        typer.Option("--features", help='Comma-separated feature list for discovery mode, e.g. "rmsd,rg".'),
    ] = None,
    dry_run: Annotated[
        bool,
        typer.Option("--dry-run", help="Validate inputs and write provenance only."),
    ] = False,
) -> None:
    """Compute a 2D Free Energy Landscape from two GROMACS XVG observables.

    \b
    Explicit mode:
        simforge fel run --xvg-x rmsd.xvg --xvg-y gyrate.xvg --temperature 300

    \b
    Discovery mode (auto-find XVG files in a run directory):
        simforge fel run /path/to/run --features rmsd,rg --temperature 300
    """
    # ── Resolve XVG paths ────────────────────────────────────────────────────
    resolved_x: Optional[Path] = None
    resolved_y: Optional[Path] = None
    resolved_xl: Optional[str] = x_label
    resolved_yl: Optional[str] = y_label

    if run_dir is not None and features is not None:
        # ── Discovery mode ──────────────────────────────────────────────────
        run_path = Path(run_dir)
        if not run_path.is_dir():
            typer.echo(f"[red]Error:[/red] '{run_dir}' is not a directory.", err=True)
            raise typer.Exit(1)

        feature_list = [f.strip() for f in features.split(",") if f.strip()]
        if len(feature_list) != 2:
            typer.echo(
                f"Error: --features must list exactly 2 features separated by comma "
                f"(got {len(feature_list)}: {feature_list}).",
                err=True,
            )
            raise typer.Exit(1)

        feat_x, feat_y = feature_list
        resolved_x = discover_feature_xvg(run_path, feat_x)
        resolved_y = discover_feature_xvg(run_path, feat_y)

        missing: list[str] = []
        if resolved_x is None:
            missing.append(f"'{feat_x}' (expected one of: rmsd.xvg, gyrate.xvg, …)")
        if resolved_y is None:
            missing.append(f"'{feat_y}' (expected one of: rmsd.xvg, gyrate.xvg, …)")
        if missing:
            typer.echo(
                f"Error: Could not auto-discover XVG files for: {', '.join(missing)}. "
                f"Use --xvg-x / --xvg-y to specify paths explicitly.",
                err=True,
            )
            raise typer.Exit(1)

        if resolved_xl is None:
            resolved_xl = feature_label(feat_x)
        if resolved_yl is None:
            resolved_yl = feature_label(feat_y)

    elif xvg_x is not None and xvg_y is not None:
        # ── Explicit mode ───────────────────────────────────────────────────
        resolved_x = Path(xvg_x)
        resolved_y = Path(xvg_y)

    else:
        typer.echo(
            "Error: provide either:\n"
            "  --xvg-x FILE --xvg-y FILE  (explicit mode)\n"
            "  RUN_DIR --features FEAT1,FEAT2  (discovery mode)",
            err=True,
        )
        raise typer.Exit(1)

    # ── Parse bins ───────────────────────────────────────────────────────────
    try:
        parts = bins.split(",")
        if len(parts) != 2:
            raise ValueError
        n_bins_x, n_bins_y = int(parts[0].strip()), int(parts[1].strip())
        if n_bins_x < 1 or n_bins_y < 1:
            raise ValueError
    except (ValueError, IndexError):
        typer.echo(
            f"Error: --bins must be 'NX,NY' with positive integers (got '{bins}').",
            err=True,
        )
        raise typer.Exit(1)

    # ── Validate time-unit options ───────────────────────────────────────────
    _VALID_UNITS = {"auto", "ps", "ns", "us"}
    for opt, val in (("--time-unit-x", time_unit_x), ("--time-unit-y", time_unit_y)):
        if val not in _VALID_UNITS:
            typer.echo(
                f"Error: {opt} must be one of {sorted(_VALID_UNITS)} (got '{val}').",
                err=True,
            )
            raise typer.Exit(1)

    if empty_bin_policy not in ("nan", "cap"):
        typer.echo(
            f"Error: --empty-bins must be 'nan' or 'cap' (got '{empty_bin_policy}').",
            err=True,
        )
        raise typer.Exit(1)

    # ── Build config and run ─────────────────────────────────────────────────
    config = FELConfig(
        xvg_x           = resolved_x,
        xvg_y           = resolved_y,
        x_column        = x_column,
        y_column        = y_column,
        x_label         = resolved_xl or "X",
        y_label         = resolved_yl or "Y",
        temperature_K   = temperature,
        n_bins_x        = n_bins_x,
        n_bins_y        = n_bins_y,
        out_dir         = Path(out),
        time_unit_x     = time_unit_x,
        time_unit_y     = time_unit_y,
        time_tolerance  = time_tolerance,
        empty_bin_policy = empty_bin_policy,
        energy_cap_kj   = energy_cap_kj,
        dry_run         = dry_run,
    )

    result = run_fel(config)

    # ── Report outcome ───────────────────────────────────────────────────────
    error_warns = [w for w in result.warnings if w.severity == "error"]
    if error_warns:
        for w in error_warns:
            typer.echo(f"Error [{w.code}]: {w.message}", err=True)
        raise typer.Exit(1)

    mode_label = "DRY RUN" if dry_run else "completed"
    typer.echo(f"FEL run {mode_label} → {out}/")
    typer.echo(f"  Frames: {result.features.n_frames if result.features else 'N/A'}")
    if result.surface and not dry_run:
        typer.echo(f"  Bins:   {n_bins_x} × {n_bins_y}")
        typer.echo(f"  Output: {out}/report.md")
