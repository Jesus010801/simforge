"""Orchestration for `simforge fel extract-minima`.

Reads an existing FEL directory, detects local minima, maps them to frames,
plans (and optionally executes) gmx trjconv extraction commands.
"""
from __future__ import annotations

import csv
import datetime
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional

from analysis.fel.models import (
    ExtractConfig, ExtractResult, FELMinimum, FELWarning, FrameMapping, StateResult,
)
from analysis.fel.minima import detect_local_minima, write_minima
from analysis.fel.frame_map import map_frames_to_minima, write_frame_mapping
from analysis.fel.plot import render_minima_plots
from analysis.fel.extract_report import write_extract_report
from analysis.fel.xvg import default_unit_for_label

# Public alias so extract_report can import without creating a circular dep.
_default_unit_for_label = default_unit_for_label

# ── Filename utilities ────────────────────────────────────────────────────────

def sanitize_system_name(name: str) -> str:
    """Produce a filesystem-safe system name component.

    Rules: spaces/dots → underscore; strip everything except [A-Za-z0-9_-];
    collapse repeated underscores; strip leading/trailing punctuation.
    """
    s = name.replace(" ", "_").replace(".", "_")
    s = re.sub(r"[^\w-]", "", s)          # keep letters, digits, _, -
    s = re.sub(r"_+", "_", s)             # collapse repeated underscores
    s = s.strip("_-")
    return s or "system"


def format_time_ns(t: float, decimals: int = 1) -> str:
    """Format time in ns to a safe filename token (decimal point → 'p').

    Example: 44.7 → '44p7', 100.0 → '100p0'.
    """
    return f"{t:.{decimals}f}".replace(".", "p")


def build_output_filename(system_name: str, minimum_id: int, time_ns: float, fmt: str) -> str:
    """Return the descriptive output filename for one representative state.

    Pattern: ``{system_name}__FEL-M{minimum_id:02d}__t{time_safe}ns.{fmt}``
    """
    return f"{system_name}__FEL-M{minimum_id:02d}__t{format_time_ns(time_ns)}ns.{fmt}"


def resolve_system_name(config: ExtractConfig) -> str:
    """Determine the system name, sanitized, from config or path stems."""
    if config.system_name:
        return sanitize_system_name(config.system_name)

    # Prefer the FEL directory name (most descriptive)
    candidate = config.fel_dir.name
    if candidate and candidate not in (".", ".."):
        s = sanitize_system_name(candidate)
        if s:
            return s

    # Topology stem (e.g. "md" from "md.tpr")
    s = sanitize_system_name(config.topology.stem)
    if s:
        return s

    # Trajectory stem
    s = sanitize_system_name(config.trajectory.stem)
    if s:
        return s

    return "system"


# ── Main workflow ─────────────────────────────────────────────────────────────

def run_extract_minima(config: ExtractConfig) -> ExtractResult:
    """Full extract-minima workflow: detect → filter → map → plan → [execute] → report."""
    result = ExtractResult(config=config)
    out_dir = config.out_dir

    # ── --force: wipe stale output before writing new results ────────────────
    if config.force and out_dir.exists():
        shutil.rmtree(out_dir)

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "metadata").mkdir(exist_ok=True)

    # ── Resolve system name once ──────────────────────────────────────────────
    system_name = resolve_system_name(config)

    # ── Load FEL metadata ────────────────────────────────────────────────────
    try:
        fel_meta = _load_fel_meta(config.fel_dir)
    except FileNotFoundError as exc:
        result.warnings.append(FELWarning("fel_meta_missing", str(exc), "error"))
        _write_warnings(result, out_dir)
        return result

    x_label = fel_meta.get("x_label", "x")
    y_label = fel_meta.get("y_label", "y")
    # Apply label-based unit defaults when the XVG ylabel had no unit annotation
    x_unit: Optional[str] = fel_meta.get("x_unit") or default_unit_for_label(x_label)
    y_unit: Optional[str] = fel_meta.get("y_unit") or default_unit_for_label(y_label)

    # ── 1. Detect local minima ────────────────────────────────────────────────
    try:
        minima, min_warns = detect_local_minima(
            config.fel_dir,
            n_minima=config.n_minima,
            min_distance_bins=config.min_distance_bins,
        )
    except (FileNotFoundError, ValueError) as exc:
        result.warnings.append(FELWarning("minima_error", str(exc), "error"))
        _write_warnings(result, out_dir)
        return result

    result.warnings.extend(min_warns)

    # ── 2. Energy filter ─────────────────────────────────────────────────────
    minima, rejected = _apply_energy_filter(minima, config.max_delta_g_kj)
    result.rejected_minima = rejected
    if rejected:
        n_all = len(minima) + len(rejected)
        result.warnings.append(FELWarning(
            "energy_filter",
            f"Only {len(minima)} of {n_all} local minima passed the "
            f"ΔG cutoff of {config.max_delta_g_kj:.2g} kJ/mol.",
            "info",
        ))

    # ── 3. Count enrichment + min_count filter ────────────────────────────────
    counts_csv = config.fel_dir / "data" / "histogram_counts.csv"
    minima = _enrich_with_counts(minima, counts_csv, config.min_count, result.warnings)

    result.minima = minima

    minima_dir = config.fel_dir / "minima"
    minima_dir.mkdir(exist_ok=True)
    result.output_files.extend(
        write_minima(minima, minima_dir, x_label, y_label, x_unit, y_unit)
    )

    # ── 4. Map frames to minima ───────────────────────────────────────────────
    features_csv = config.fel_dir / "data" / "features.csv"
    try:
        frame_map, map_warns = map_frames_to_minima(
            minima, features_csv, x_label, y_label, x_unit, y_unit,
        )
    except (FileNotFoundError, ValueError) as exc:
        result.warnings.append(FELWarning("frame_map_error", str(exc), "error"))
        _write_warnings(result, out_dir)
        return result

    result.frame_map = frame_map
    result.warnings.extend(map_warns)
    result.output_files.extend(write_frame_mapping(frame_map, minima_dir))

    # ── 5. Plan/execute extraction for each state ─────────────────────────────
    for minimum in minima:
        mapping = next(m for m in frame_map if m.minimum_id == minimum.minimum_id)
        state_dir = out_dir / f"state_{minimum.minimum_id:03d}"
        state_dir.mkdir(exist_ok=True)

        state_res = _process_state(
            minimum, mapping, config, state_dir,
            x_label, y_label, system_name,
        )
        result.state_results.append(state_res)

    # ── 6. Write planned_commands.sh ─────────────────────────────────────────
    sh_path = _write_planned_commands(result, out_dir, x_label, y_label, x_unit, y_unit)
    result.output_files.append(sh_path)

    # ── 7. Minima overlay plots ───────────────────────────────────────────────
    surface_csv = config.fel_dir / "data" / "free_energy_surface.csv"
    plots_dir = config.fel_dir / "plots"
    plots_dir.mkdir(exist_ok=True)
    try:
        plot_paths = render_minima_plots(
            surface_csv, minima, frame_map,
            x_label, y_label, x_unit, y_unit,
            plots_dir, result.warnings,
        )
        result.output_files.extend(plot_paths)
    except Exception as exc:
        result.warnings.append(FELWarning("overlay_plot_error", str(exc), "warn"))

    # ── 8. Metadata files ─────────────────────────────────────────────────────
    _write_metadata(result, out_dir, fel_meta)

    # ── 9. Report ─────────────────────────────────────────────────────────────
    report_path = write_extract_report(result, out_dir, fel_meta)
    result.output_files.append(report_path)

    # ── 10. Copy minima outputs to out_dir for self-contained access ──────────
    for fname, src in {
        "minima_summary.csv": minima_dir / "minima.csv",
        "frame_mapping.csv":  minima_dir / "frame_mapping.csv",
    }.items():
        if src.exists():
            dst = out_dir / fname
            shutil.copy2(str(src), str(dst))
            result.output_files.append(dst)

    return result


# ── Internal helpers ──────────────────────────────────────────────────────────

def _apply_energy_filter(
    minima: list[FELMinimum],
    max_delta_g_kj: Optional[float],
) -> tuple[list[FELMinimum], list[FELMinimum]]:
    """Split minima into (kept, rejected) based on max ΔG cutoff.

    When all minima exceed the cutoff, the global minimum is kept so the caller
    always receives at least one state.
    """
    if max_delta_g_kj is None or not minima:
        return minima, []
    kept = [m for m in minima if m.free_energy_kJ_mol <= max_delta_g_kj]
    rejected = [m for m in minima if m.free_energy_kJ_mol > max_delta_g_kj]
    if not kept:
        kept = [minima[0]]
        rejected = minima[1:]
    return kept, rejected


def _enrich_with_counts(
    minima: list[FELMinimum],
    counts_csv: Path,
    min_count: int,
    warnings: list[FELWarning],
) -> list[FELMinimum]:
    """Attach frame counts/probability from histogram_counts.csv; apply min_count filter."""
    if not counts_csv.exists():
        return minima

    count_map: dict[tuple[str, str], int] = {}
    try:
        with counts_csv.open(newline="") as fh:
            reader = csv.reader(fh)
            next(reader)  # header
            for row in reader:
                if len(row) < 3:
                    continue
                try:
                    xk = f"{float(row[0]):.6g}"
                    yk = f"{float(row[1]):.6g}"
                    count_map[(xk, yk)] = int(float(row[2]))
                except (ValueError, IndexError):
                    continue
    except Exception:
        return minima

    if not count_map:
        return minima

    total = sum(count_map.values()) or 1

    for m in minima:
        key = (f"{m.x:.6g}", f"{m.y:.6g}")
        c = count_map.get(key)
        if c is not None:
            m.count = c
            m.probability = c / total

    if min_count > 1:
        kept = [m for m in minima if m.count is None or m.count >= min_count]
        filtered = [m for m in minima if m.count is not None and m.count < min_count]
        if filtered:
            warnings.append(FELWarning(
                "min_count_filter",
                f"{len(filtered)} minima removed: bin count < {min_count}.",
                "info",
            ))
        return kept

    return minima


def _load_fel_meta(fel_dir: Path) -> dict:
    config_path = fel_dir / "config.json"
    if not config_path.exists():
        raise FileNotFoundError(
            f"FEL config not found: {config_path}. "
            "Run `simforge fel run` first."
        )
    return json.loads(config_path.read_text())


def _process_state(
    minimum: FELMinimum,
    mapping: FrameMapping,
    config: ExtractConfig,
    state_dir: Path,
    x_label: str,
    y_label: str,
    system_name: str,
) -> StateResult:
    dump_time_ps = mapping.selected_time_ns * 1000.0
    formats = ["pdb", "gro"] if config.format == "both" else [config.format]

    shell_cmds: list[str] = []
    cmd_lists: list[list[str]] = []
    output_paths: list[Path] = []
    output_filenames: list[str] = []

    # Build descriptive filename pattern
    filename_pattern = "{system_name}__FEL-M{mid:02d}__t{time_safe}ns.{ext}"
    display_name = (
        f"{system_name} FEL-M{minimum.minimum_id:02d} "
        f"at {mapping.selected_time_ns:.1f} ns"
    )

    for fmt in formats:
        fname = build_output_filename(
            system_name, minimum.minimum_id, mapping.selected_time_ns, fmt
        )
        out_file = state_dir / fname
        output_paths.append(out_file)
        output_filenames.append(fname)

        cmd_args = [
            config.gmx, "trjconv",
            "-s", str(config.topology),
            "-f", str(config.trajectory),
            "-o", str(out_file),
            "-dump", f"{dump_time_ps:.3f}",
        ]
        if config.index:
            cmd_args += ["-n", str(config.index)]
        cmd_lists.append(cmd_args)

        shell_cmds.append(f'echo "{config.selection}" | {" ".join(cmd_args)}')

    meta: dict = {
        "minimum_id": minimum.minimum_id,
        "display_name": display_name,
        "system_name": system_name,
        "filename_pattern": filename_pattern,
        "output_filenames": output_filenames,
        "output_paths": [str(p) for p in output_paths],
        "target_x": mapping.target_x,
        "target_y": mapping.target_y,
        "x_label": x_label,
        "y_label": y_label,
        "x_unit": mapping.x_unit,
        "y_unit": mapping.y_unit,
        "free_energy_kJ_mol": mapping.free_energy_kJ_mol,
        "selected_time_ns": mapping.selected_time_ns,
        "dump_time_ps": dump_time_ps,
        "selected_frame_index": mapping.selected_frame_index,
        "frame_x": mapping.frame_x,
        "frame_y": mapping.frame_y,
        "distance_in_feature_space": mapping.distance_in_feature_space,
        "planned_commands": shell_cmds,
        "status": "planned",
    }
    if minimum.count is not None:
        meta["count"] = minimum.count
    if minimum.probability is not None:
        meta["probability"] = minimum.probability

    if config.compat_receptor_name:
        # Write shell symlink commands for backwards compatibility
        for fmt, fname in zip(formats, output_filenames):
            compat = state_dir / f"receptor.{fmt}"
            meta.setdefault("compat_paths", []).append(str(compat))

    if config.dry_run:
        (state_dir / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
        return StateResult(
            minimum_id=minimum.minimum_id,
            planned_commands=shell_cmds,
            planned_cmd_lists=cmd_lists,
            output_paths=output_paths,
            status="planned",
        )

    # Execute mode
    return_codes: list[int] = []
    stdout_parts: list[str] = []
    stderr_parts: list[str] = []
    status = "completed"

    for cmd_args, out_file in zip(cmd_lists, output_paths):
        try:
            proc = subprocess.run(
                cmd_args,
                input=f"{config.selection}\n",
                text=True,
                capture_output=True,
            )
            return_codes.append(proc.returncode)
            stdout_parts.append(proc.stdout or "")
            stderr_parts.append(proc.stderr or "")
            if proc.returncode != 0:
                status = "failed"
            elif config.compat_receptor_name and out_file.exists():
                # Create symlink alias receptor.ext → descriptive name
                fmt = out_file.suffix.lstrip(".")
                compat = state_dir / f"receptor.{fmt}"
                if not compat.exists():
                    try:
                        compat.symlink_to(out_file.name)
                    except OSError:
                        pass  # non-critical
        except FileNotFoundError as exc:
            return_codes.append(-1)
            stderr_parts.append(str(exc))
            status = "failed"

    meta["status"] = status
    meta["return_codes"] = return_codes
    (state_dir / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")

    return StateResult(
        minimum_id=minimum.minimum_id,
        planned_commands=shell_cmds,
        planned_cmd_lists=cmd_lists,
        output_paths=output_paths,
        status=status,
        return_codes=return_codes,
        stdout_lines=stdout_parts,
        stderr_lines=stderr_parts,
    )


def _write_planned_commands(
    result: ExtractResult,
    out_dir: Path,
    x_label: str,
    y_label: str,
    x_unit: Optional[str],
    y_unit: Optional[str],
) -> Path:
    xl = f"{x_label} ({x_unit})" if x_unit else x_label
    yl = f"{y_label} ({y_unit})" if y_unit else y_label

    lines = ["#!/bin/bash", "# SimForge FEL extraction commands", ""]
    for sr, m, fm in zip(result.state_results, result.minima, result.frame_map):
        # Derive display_name from the first output filename if available
        if sr.output_paths:
            fname = sr.output_paths[0].name
        else:
            fname = f"FEL-M{m.minimum_id:02d}"
        lines += [
            f"# {fname}",
            f"# {xl}: {fm.target_x:.4g}  {yl}: {fm.target_y:.4g}"
            f"  ΔG: {fm.free_energy_kJ_mol:.4g} kJ/mol",
            f"# Frame index: {fm.selected_frame_index},"
            f" Time: {fm.selected_time_ns:.4g} ns"
            f" (dump_time: {fm.selected_time_ns * 1000:.3f} ps)",
            "",
        ]
        for cmd in sr.planned_commands:
            lines.append(cmd)
        lines.append("")

    sh_path = out_dir / "planned_commands.sh"
    sh_path.write_text("\n".join(lines) + "\n")
    sh_path.chmod(0o755)
    return sh_path


def _write_metadata(result: ExtractResult, out_dir: Path, fel_meta: dict) -> None:
    try:
        from simforge import __version__ as _v
        version = _v
    except Exception:
        version = "unknown"

    prov = {
        "timestamp": datetime.datetime.utcnow().isoformat() + "Z",
        "simforge_version": version,
        "argv": sys.argv[:],
        "config": result.config.to_dict(),
        "n_minima_requested": result.config.n_minima,
        "n_minima_found": len(result.minima),
        "n_minima_rejected": len(result.rejected_minima),
        "max_delta_g_kj": result.config.max_delta_g_kj,
        "fel_x_label": fel_meta.get("x_label"),
        "fel_y_label": fel_meta.get("y_label"),
        "fel_x_unit": fel_meta.get("x_unit"),
        "fel_y_unit": fel_meta.get("y_unit"),
    }
    (out_dir / "metadata" / "extraction_provenance.json").write_text(
        json.dumps(prov, indent=2) + "\n"
    )

    cmds_meta = [sr.to_dict() for sr in result.state_results]
    (out_dir / "metadata" / "extraction_commands.json").write_text(
        json.dumps(cmds_meta, indent=2) + "\n"
    )

    if result.rejected_minima:
        rej_data = [m.to_dict() for m in result.rejected_minima]
        (out_dir / "metadata" / "rejected_minima.json").write_text(
            json.dumps(rej_data, indent=2) + "\n"
        )

    _write_warnings(result, out_dir)


def _write_warnings(result: ExtractResult, out_dir: Path) -> None:
    (out_dir / "metadata").mkdir(exist_ok=True)
    (out_dir / "metadata" / "warnings.json").write_text(
        json.dumps([w.to_dict() for w in result.warnings], indent=2) + "\n"
    )
