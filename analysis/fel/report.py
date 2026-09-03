"""Markdown report writer for FEL runs."""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Optional

from analysis.fel.models import FELRunResult
from analysis.fel.provenance import FELProvenance

_MINIMA_README = """\
# Minima Extraction — Phase 3

This directory is reserved for free-energy minima extraction.
Minima extraction is **not yet implemented** and will be added in Phase 3.

When Phase 3 is complete, this directory will contain:
- `minima.json`     — coordinates and energies of detected free-energy minima
- `minima.csv`      — tabular format
- Frame indices and trajectory snapshots closest to each minimum
"""


def write_minima_readme(out_dir: Path) -> Path:
    minima_dir = out_dir / "minima"
    minima_dir.mkdir(parents=True, exist_ok=True)
    p = minima_dir / "README.md"
    p.write_text(_MINIMA_README)
    return p


def write_report(
    result: FELRunResult,
    prov: FELProvenance,
    out_dir: Path,
    output_files: list[Path],
) -> Path:
    lines: list[str] = [
        "# SimForge FEL Run Report",
        "",
        f"**Generated:** {prov.timestamp}  ",
        f"**SimForge version:** {prov.simforge_version}",
        "",
        "---",
        "",
        "## Inputs",
        "",
        f"- **X axis XVG:** `{result.config.xvg_x}`",
        f"- **Y axis XVG:** `{result.config.xvg_y}`",
        f"- **X column (1-based):** {result.config.x_column}",
        f"- **Y column (1-based):** {result.config.y_column}",
        f"- **X label:** {result.config.x_label}",
        f"- **Y label:** {result.config.y_label}",
        "",
    ]

    if result.features:
        ft = result.features
        lines += [
            "## Feature Table",
            "",
            f"- **Frames used:** {ft.n_frames}",
            f"- **Time conversion:** {ft.time_conversion_note}",
            f"- **Time unit (normalized):** {ft.time_unit_normalized}",
            f"- **X time unit (original):** {ft.x_unit_original or 'unknown (assumed ps)'}",
            f"- **Y time unit (original):** {ft.y_unit_original or 'unknown (assumed ps)'}",
            f"- **X feature unit:** {ft.x_unit or 'unknown'}",
            f"- **Y feature unit:** {ft.y_unit or 'unknown'}",
            "",
        ]

    lines += [
        "## FEL Configuration",
        "",
        f"- **Temperature:** {result.config.temperature_K} K",
        f"- **Bins (X × Y):** {result.config.n_bins_x} × {result.config.n_bins_y}",
        f"- **Empty bin policy:** {result.config.empty_bin_policy}",
    ]
    if result.config.empty_bin_policy == "cap":
        lines.append(f"- **Energy cap:** {result.config.energy_cap_kj} kJ/mol")
    lines.append("")

    if result.surface:
        surf = result.surface
        lines += [
            "## FEL Computation",
            "",
            f"- **Formula:** `{surf.formula}`",
            f"- **Frames used:** {surf.n_frames_used}",
            f"- **X range:** [{surf.x_edges[0]:.4g}, {surf.x_edges[-1]:.4g}]",
            f"- **Y range:** [{surf.y_edges[0]:.4g}, {surf.y_edges[-1]:.4g}]",
            "",
        ]

    # Collect all warnings (needed for Visualizations section and Warnings section)
    all_warnings = list(result.warnings)
    if result.features:
        all_warnings += result.features.warnings
    if result.surface:
        all_warnings += result.surface.warnings

    # Visualizations section
    _PLOT_DESC = {
        "fel_contour.png":    "2D contour plot",
        "fel_surface_3d.png": "3D surface plot (static)",
        "fel_surface_3d.html": "Interactive 3D surface",
    }
    plots_dir = out_dir / "plots"
    plot_files = sorted(
        [p for p in output_files if p.parent == plots_dir],
        key=lambda p: p.name,
    )
    lines += ["## Visualizations", ""]
    if plot_files:
        for pf in plot_files:
            rel = pf.relative_to(out_dir)
            desc = _PLOT_DESC.get(pf.name, pf.name)
            lines.append(f"- {desc}: `{rel}`")
    else:
        lines.append("No plots generated (matplotlib not available).")
    no_plotly = any(w.code == "no_plotly" for w in all_warnings)
    if no_plotly:
        lines.append("")
        lines.append("> Interactive 3D plot was skipped because Plotly is not installed.")
    lines.append("")

    if all_warnings:
        lines += ["## Warnings", ""]
        for w in all_warnings:
            icon = "⚠" if w.severity == "warn" else ("✗" if w.severity == "error" else "ℹ")
            lines.append(f"- {icon} `{w.code}`: {w.message}")
        lines.append("")

    lines += [
        "## Output Files",
        "",
    ]
    for p in sorted(output_files, key=str):
        rel = p.relative_to(out_dir) if p.is_relative_to(out_dir) else p
        lines.append(f"- `{rel}`")
    lines.append("")

    # Minima section: written by extract-minima after the fact
    minima_csv = out_dir / "minima" / "frame_mapping.csv"
    if minima_csv.exists():
        lines += _minima_table_section(minima_csv, out_dir)

    lines += [
        "---",
        "",
        "## Phase Status",
        "",
        "- ✓ **Phase 1** — MD run discovery and observable planning",
        "- ✓ **Phase 2** — FEL surface computation (this run)",
        "- ○ **Phase 3** — Free-energy minima extraction (`simforge fel extract-minima`)",
        "- ○ **Phase 4** — Docking ensemble generation",
        "",
    ]

    if result.config.dry_run:
        lines.insert(3, "**Mode:** DRY RUN — surface was not computed.\n")

    report_path = out_dir / "report.md"
    report_path.write_text("\n".join(lines) + "\n")
    return report_path


def _minima_table_section(frame_mapping_csv: Path, out_dir: Path) -> list[str]:
    """Build a Markdown table from frame_mapping.csv (written by extract-minima)."""
    try:
        with frame_mapping_csv.open(newline="") as fh:
            rows = list(csv.DictReader(fh))
    except Exception:
        return []
    if not rows:
        return []

    x_unit = rows[0].get("x_unit", "")
    y_unit = rows[0].get("y_unit", "")
    x_hdr = f"X ({x_unit})" if x_unit else "X"
    y_hdr = f"Y ({y_unit})" if y_unit else "Y"

    lines = [
        "## Free-energy Minima",
        "",
        f"| Minimum | {x_hdr} | {y_hdr} | ΔG (kJ/mol) | Time (ns) | Frame |",
        f"|---------|{'-' * (len(x_hdr)+2)}|{'-' * (len(y_hdr)+2)}|-------------|-----------|-------|",
    ]
    for r in rows:
        lines.append(
            f"| M{r['minimum_id']} "
            f"| {float(r['target_x']):.4g} "
            f"| {float(r['target_y']):.4g} "
            f"| {float(r['free_energy_kJ_mol']):.4g} "
            f"| {float(r['selected_time_ns']):.4g} "
            f"| {r['selected_frame_index']} |"
        )
    lines.append("")
    return lines
