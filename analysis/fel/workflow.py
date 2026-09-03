"""Orchestrate a full FEL computation from config to output files."""
from __future__ import annotations

import json
from pathlib import Path

from analysis.fel.models import FELConfig, FELRunResult, FELWarning
from analysis.fel.xvg import load_xvg_series
from analysis.fel.features import build_feature_table
from analysis.fel.surface import compute_fel_surface
from analysis.fel.provenance import build_provenance
from analysis.fel.report import write_report, write_minima_readme
from analysis.fel.plot import render_plots


def run_fel(config: FELConfig) -> FELRunResult:
    """Execute a full FEL run and write all output files.

    In dry-run mode the surface is not computed; only the XVG headers and
    feature-table alignment are validated, and the output directory tree plus
    provenance are written so the user can inspect what would happen.
    """
    result = FELRunResult(config=config)
    out_dir = config.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "metadata").mkdir(exist_ok=True)
    (out_dir / "data").mkdir(exist_ok=True)

    # ── 1. Parse XVGs ────────────────────────────────────────────────────────
    try:
        result.series_x = load_xvg_series(
            config.xvg_x, config.x_column,
            override_unit=None if config.time_unit_x == "auto" else config.time_unit_x,
        )
    except (FileNotFoundError, ValueError) as exc:
        result.warnings.append(FELWarning("xvg_x_error", str(exc), "error"))
        _write_provenance(result)
        return result

    try:
        result.series_y = load_xvg_series(
            config.xvg_y, config.y_column,
            override_unit=None if config.time_unit_y == "auto" else config.time_unit_y,
        )
    except (FileNotFoundError, ValueError) as exc:
        result.warnings.append(FELWarning("xvg_y_error", str(exc), "error"))
        _write_provenance(result)
        return result

    # ── 2. Build feature table ───────────────────────────────────────────────
    try:
        result.features = build_feature_table(
            result.series_x,
            result.series_y,
            x_label=config.x_label,
            y_label=config.y_label,
            time_tolerance=config.time_tolerance,
        )
    except ValueError as exc:
        result.warnings.append(FELWarning("feature_alignment_error", str(exc), "error"))
        _write_provenance(result)
        return result

    # ── 3. Write features.csv; propagate unit metadata to config ────────────
    feat_csv = out_dir / "data" / "features.csv"
    feat_csv.write_text("\n".join(result.features.to_csv_lines()) + "\n")
    result.output_files.append(feat_csv)

    # Propagate detected feature units back into config for config.json
    config.x_unit = result.features.x_unit
    config.y_unit = result.features.y_unit

    if config.dry_run:
        result.warnings.append(FELWarning(
            "dry_run", "Dry run — FEL surface not computed.", "info"
        ))
        _finalize(result, out_dir)
        return result

    # ── 4. Compute FEL surface ───────────────────────────────────────────────
    try:
        result.surface = compute_fel_surface(
            result.features.x,
            result.features.y,
            n_bins_x=config.n_bins_x,
            n_bins_y=config.n_bins_y,
            temperature_K=config.temperature_K,
            empty_bin_policy=config.empty_bin_policy,
            energy_cap_kj=config.energy_cap_kj,
            x_label=config.x_label,
            y_label=config.y_label,
            x_unit=result.features.x_unit,
            y_unit=result.features.y_unit,
        )
    except ValueError as exc:
        result.warnings.append(FELWarning("surface_error", str(exc), "error"))
        _finalize(result, out_dir)
        return result

    # ── 5. Write 2D CSV outputs ──────────────────────────────────────────────
    surf = result.surface
    counts_csv = out_dir / "data" / "histogram_counts.csv"
    prob_csv   = out_dir / "data" / "probability.csv"
    G_csv      = out_dir / "data" / "free_energy_surface.csv"

    counts_csv.write_text("\n".join(surf.counts_csv()) + "\n")
    prob_csv.write_text("\n".join(surf.probability_csv()) + "\n")
    G_csv.write_text("\n".join(surf.free_energy_csv()) + "\n")

    result.output_files += [counts_csv, prob_csv, G_csv]

    # ── 6. Optional plots ────────────────────────────────────────────────────
    result.output_files.extend(render_plots(result, out_dir))

    _finalize(result, out_dir)
    return result


# ── Internal helpers ──────────────────────────────────────────────────────────

def _finalize(result: FELRunResult, out_dir: Path) -> None:
    """Write config, provenance, warnings, minima placeholder, and report."""
    config_path = out_dir / "config.json"
    config_path.write_text(json.dumps(result.config.to_dict(), indent=2) + "\n")
    result.output_files.append(config_path)

    prov = build_provenance(result)
    prov_path = out_dir / "metadata" / "provenance.json"
    prov_path.write_text(json.dumps(prov.to_dict(), indent=2) + "\n")
    result.output_files.append(prov_path)

    all_warnings = list(result.warnings)
    if result.features:
        all_warnings += result.features.warnings
    if result.surface:
        all_warnings += result.surface.warnings
    warn_path = out_dir / "metadata" / "warnings.json"
    warn_path.write_text(json.dumps([w.to_dict() for w in all_warnings], indent=2) + "\n")
    result.output_files.append(warn_path)

    minima_readme = write_minima_readme(out_dir)
    result.output_files.append(minima_readme)

    report_path = write_report(result, prov, out_dir, result.output_files)
    result.output_files.append(report_path)


def _write_provenance(result: FELRunResult) -> None:
    """Write minimal provenance when aborting early due to errors."""
    out_dir = result.config.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "metadata").mkdir(exist_ok=True)
    prov = build_provenance(result)
    (out_dir / "metadata" / "provenance.json").write_text(
        json.dumps(prov.to_dict(), indent=2) + "\n"
    )


