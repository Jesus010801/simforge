"""Map FEL local minima to nearest trajectory frames via features.csv."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Optional

from analysis.fel.models import FELMinimum, FrameMapping, FELWarning


# ── Public API ────────────────────────────────────────────────────────────────

def map_frames_to_minima(
    minima: list[FELMinimum],
    features_csv: Path,
    x_label: str,
    y_label: str,
    x_unit: Optional[str] = None,
    y_unit: Optional[str] = None,
) -> tuple[list[FrameMapping], list[FELWarning]]:
    """For each minimum, find the nearest trajectory frame in feature space.

    Distance metric: Euclidean in (x, y) feature space.

    Deterministic tie-breaking order:
    1. lower Euclidean distance
    2. lower time_ns
    3. lower row index
    """
    if not features_csv.exists():
        raise FileNotFoundError(
            f"Features CSV not found: {features_csv}. "
            "Run `simforge fel run` first."
        )

    frames = _read_features_csv(features_csv, x_label, y_label)
    if not frames:
        raise ValueError(f"No usable frames in {features_csv}.")

    warnings: list[FELWarning] = []
    mappings: list[FrameMapping] = []

    for minimum in minima:
        best_idx, best_dist = _nearest_frame(frames, minimum.x, minimum.y)
        t_ns, fx, fy = frames[best_idx]
        mappings.append(FrameMapping(
            minimum_id=minimum.minimum_id,
            target_x=minimum.x,
            target_y=minimum.y,
            x_unit=x_unit,
            y_unit=y_unit,
            free_energy_kJ_mol=minimum.free_energy_kJ_mol,
            selected_time_ns=t_ns,
            selected_frame_index=best_idx,
            frame_x=fx,
            frame_y=fy,
            distance_in_feature_space=best_dist,
        ))

    return mappings, warnings


def write_frame_mapping(
    frame_map: list[FrameMapping],
    out_dir: Path,
) -> list[Path]:
    """Write ``frame_mapping.csv`` and ``frame_mapping.json`` to *out_dir*."""
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = out_dir / "frame_mapping.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "minimum_id", "target_x", "target_y", "x_unit", "y_unit",
            "free_energy_kJ_mol", "selected_time_ns",
            "selected_frame_index", "frame_x", "frame_y",
            "distance_in_feature_space",
        ])
        for m in frame_map:
            w.writerow([
                m.minimum_id,
                f"{m.target_x:.6g}", f"{m.target_y:.6g}",
                m.x_unit or "", m.y_unit or "",
                f"{m.free_energy_kJ_mol:.6g}",
                f"{m.selected_time_ns:.6g}",
                m.selected_frame_index,
                f"{m.frame_x:.6g}", f"{m.frame_y:.6g}",
                f"{m.distance_in_feature_space:.6g}",
            ])

    json_path = out_dir / "frame_mapping.json"
    json_path.write_text(json.dumps([m.to_dict() for m in frame_map], indent=2) + "\n")

    return [csv_path, json_path]


# ── Internal helpers ──────────────────────────────────────────────────────────

def _read_features_csv(
    path: Path,
    x_label: str,
    y_label: str,
) -> list[tuple[float, float, float]]:
    """Read features.csv and return list of (time_ns, x, y)."""
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        rows: list[tuple[float, float, float]] = []
        for row in reader:
            try:
                t = float(row["time_ns"])
                x = float(row[x_label])
                y = float(row[y_label])
                rows.append((t, x, y))
            except (KeyError, ValueError):
                continue
    return rows


def _nearest_frame(
    frames: list[tuple[float, float, float]],
    tx: float,
    ty: float,
) -> tuple[int, float]:
    """Return (index, distance) of the frame nearest to (tx, ty).

    Tie-break: lower distance → lower time_ns → lower row index.
    """
    best_idx = 0
    best_dist = float("inf")
    best_t = float("inf")

    for i, (t, fx, fy) in enumerate(frames):
        dist = math.sqrt((fx - tx) ** 2 + (fy - ty) ** 2)
        if dist < best_dist - 1e-12:
            best_dist, best_idx, best_t = dist, i, t
        elif abs(dist - best_dist) < 1e-12:
            if t < best_t - 1e-12:
                best_idx, best_t = i, t
            elif abs(t - best_t) < 1e-12 and i < best_idx:
                best_idx = i

    return best_idx, best_dist
