"""Report writer for the representative_states/ output directory."""
from __future__ import annotations

from pathlib import Path
from typing import Optional

from analysis.fel.models import ExtractResult
from analysis.fel.xvg import default_unit_for_label as _default_unit_for_label


def write_extract_report(
    result: ExtractResult,
    out_dir: Path,
    fel_meta: dict,
) -> Path:
    x_label = fel_meta.get("x_label", "x")
    y_label = fel_meta.get("y_label", "y")
    x_unit: Optional[str] = fel_meta.get("x_unit") or _default_unit_for_label(x_label)
    y_unit: Optional[str] = fel_meta.get("y_unit") or _default_unit_for_label(y_label)
    xl = f"{x_label} ({x_unit})" if x_unit else x_label
    yl = f"{y_label} ({y_unit})" if y_unit else y_label

    cfg = result.config
    mode = "DRY RUN" if cfg.dry_run else "EXECUTE"
    n_found = len(result.minima)
    n_rejected = len(result.rejected_minima)

    filter_note = ""
    if cfg.max_delta_g_kj is not None:
        filter_note = f"{cfg.max_delta_g_kj:.2g} kJ/mol"
        if n_rejected:
            filter_note += f" ({n_rejected} rejected)"
    else:
        filter_note = "disabled"

    lines: list[str] = [
        "# SimForge FEL Extraction Report",
        "",
        f"**Mode:** {mode}  ",
        f"**Source FEL directory:** `{cfg.fel_dir}`  ",
        f"**Output directory:** `{out_dir}`",
        "",
        "---",
        "",
        "## Configuration",
        "",
        f"- **Topology:** `{cfg.topology}`",
        f"- **Trajectory:** `{cfg.trajectory}`",
        *(([f"- **Structure:** `{cfg.structure}`"]) if cfg.structure else []),
        *(([f"- **Index:** `{cfg.index}`"]) if cfg.index else []),
        f"- **Selection:** `{cfg.selection}`",
        f"- **N minima requested:** {cfg.n_minima}",
        f"- **N minima extracted:** {n_found}",
        f"- **Min distance bins:** {cfg.min_distance_bins}",
        f"- **ΔG filter:** {filter_note}",
        f"- **Min bin count:** {cfg.min_count}",
        f"- **Output format:** {cfg.format}",
        f"- **GROMACS binary:** `{cfg.gmx}`",
        "",
        "## Feature Units",
        "",
        f"- **X axis:** {xl}",
        f"- **Y axis:** {yl}",
        "- **Time axis:** ns (normalized)",
        "",
        "## Minima and Mapped Frames",
        "",
    ]

    if result.frame_map:
        from analysis.fel.extract import resolve_system_name, build_output_filename
        system_name = resolve_system_name(result.config)

        count_col = any(m.count is not None for m in result.minima)
        count_hdr = " | Count" if count_col else ""
        count_sep = "-------|" if count_col else ""
        lines += [
            f"| Display Name | {xl} | {yl} | ΔG (kJ/mol) | Time (ns) | Output File | Status |{count_hdr}",
            f"|--------------|{'-'*(len(xl)+2)}|{'-'*(len(yl)+2)}|-------------|"
            f"-----------|-------------|--------|{count_sep}",
        ]
        fm_by_id = {fm.minimum_id: fm for fm in result.frame_map}
        sr_by_id = {sr.minimum_id: sr for sr in result.state_results}
        fmt = result.config.format if result.config.format != "both" else "pdb"
        for m in result.minima:
            fm = fm_by_id.get(m.minimum_id)
            sr = sr_by_id.get(m.minimum_id)
            status = sr.status if sr else "—"
            if fm:
                display_name = (
                    f"{system_name} FEL-M{m.minimum_id:02d} "
                    f"at {fm.selected_time_ns:.1f} ns"
                )
                out_fname = build_output_filename(
                    system_name, m.minimum_id, fm.selected_time_ns, fmt
                )
                count_cell = f" | {m.count}" if count_col else ""
                lines.append(
                    f"| {display_name} "
                    f"| {fm.target_x:.4g} "
                    f"| {fm.target_y:.4g} "
                    f"| {fm.free_energy_kJ_mol:.4g} "
                    f"| {fm.selected_time_ns:.1f} "
                    f"| `{out_fname}` "
                    f"| {status} |{count_cell}"
                )
    else:
        lines.append("*No frame mapping available.*")

    if result.rejected_minima:
        lines += [
            "",
            f"### Rejected by ΔG filter (>{cfg.max_delta_g_kj:.2g} kJ/mol)",
            "",
            f"| Minimum | {xl} | {yl} | ΔG (kJ/mol) |",
            f"|---------|{'-'*(len(xl)+2)}|{'-'*(len(yl)+2)}|-------------|",
        ]
        for i, m in enumerate(result.rejected_minima, start=n_found + 1):
            lines.append(
                f"| M{i} (rejected) | {m.x:.4g} | {m.y:.4g} | {m.free_energy_kJ_mol:.4g} |"
            )

    lines += ["", "## Extraction Commands", ""]
    for sr in result.state_results:
        for cmd in sr.planned_commands:
            lines += ["```bash", cmd, "```", ""]

    if result.warnings:
        lines += ["## Warnings", ""]
        for w in result.warnings:
            icon = "⚠" if w.severity == "warn" else ("✗" if w.severity == "error" else "ℹ")
            lines.append(f"- {icon} `{w.code}`: {w.message}")
        lines.append("")

    lines += [
        "---",
        "",
        "## Phase Status",
        "",
        "- ✓ **Phase 2** — FEL surface computation",
        "- ✓ **Phase 3** — Minima detection and representative extraction (this run)",
        "- ○ **Phase 4** — Docking ensemble generation (not yet implemented)",
        "",
    ]

    report_path = out_dir / "report.md"
    report_path.write_text("\n".join(lines) + "\n")
    return report_path
