"""Local free-energy minima detection from an existing FEL directory."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Optional

from analysis.fel.models import FELMinimum, FELWarning


# ── Public API ────────────────────────────────────────────────────────────────

def detect_local_minima(
    fel_dir: Path,
    n_minima: int = 5,
    min_distance_bins: int = 2,
) -> tuple[list[FELMinimum], list[FELWarning]]:
    """Detect local minima in the FEL surface from an existing FEL directory.

    Reads ``data/free_energy_surface.csv``.  NaN bins are ignored.  Uses
    8-neighborhood comparison; deduplicates with Euclidean bin-space distance;
    falls back to the global finite minimum when no local minima exist.
    """
    surface_csv = fel_dir / "data" / "free_energy_surface.csv"
    if not surface_csv.exists():
        raise FileNotFoundError(
            f"FEL surface CSV not found: {surface_csv}. "
            "Run `simforge fel run` first."
        )

    rows, x_label, y_label = _read_surface_csv(surface_csv)
    if not rows:
        raise ValueError(f"No finite energy rows in {surface_csv}.")

    xs = sorted(set(r[0] for r in rows))
    ys = sorted(set(r[1] for r in rows))
    nx, ny = len(xs), len(ys)

    xi = {v: i for i, v in enumerate(xs)}
    yi = {v: i for i, v in enumerate(ys)}

    G: list[list[float]] = [[float("nan")] * ny for _ in range(nx)]
    for x, y, g in rows:
        G[xi[x]][yi[y]] = g

    candidates = _find_local_minima(G, xs, ys, nx, ny)
    candidates.sort(key=lambda m: m.free_energy_kJ_mol)

    warnings: list[FELWarning] = []

    if not candidates:
        warnings.append(FELWarning(
            "no_local_minima",
            "No local minima found in the FEL surface; falling back to global minimum.",
            "warn",
        ))
        gmin = _global_minimum(G, xs, ys, nx, ny)
        if gmin is None:
            raise ValueError("No finite energy bins — cannot detect minima.")
        candidates = [gmin]

    selected = _deduplicate(candidates, min_distance_bins, n_minima)
    for i, m in enumerate(selected):
        m.minimum_id = i + 1

    return selected, warnings


def write_minima(
    minima: list[FELMinimum],
    out_dir: Path,
    x_label: str,
    y_label: str,
    x_unit: Optional[str] = None,
    y_unit: Optional[str] = None,
) -> list[Path]:
    """Write ``minima.csv`` and ``minima.json`` to *out_dir*."""
    out_dir.mkdir(parents=True, exist_ok=True)

    xl = f"{x_label} ({x_unit})" if x_unit else x_label
    yl = f"{y_label} ({y_unit})" if y_unit else y_label

    csv_path = out_dir / "minima.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["minimum_id", xl, yl, "free_energy_kJ_mol", "x_bin", "y_bin"])
        for m in minima:
            w.writerow([
                m.minimum_id,
                f"{m.x:.6g}", f"{m.y:.6g}",
                f"{m.free_energy_kJ_mol:.6g}",
                m.x_bin, m.y_bin,
            ])

    json_path = out_dir / "minima.json"
    json_path.write_text(json.dumps([m.to_dict() for m in minima], indent=2) + "\n")

    return [csv_path, json_path]


# ── Internal helpers ──────────────────────────────────────────────────────────

def _read_surface_csv(
    path: Path,
) -> tuple[list[tuple[float, float, float]], str, str]:
    """Read the flat FEL surface CSV. Returns (rows, x_label, y_label).

    NaN values are excluded from rows.
    """
    with path.open(newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        x_col = header[0][2:] if header[0].startswith("x_") else header[0]
        y_col = header[1][2:] if header[1].startswith("y_") else header[1]
        rows: list[tuple[float, float, float]] = []
        for row in reader:
            if len(row) < 3:
                continue
            try:
                g = float(row[2])
                if g != g:  # NaN
                    continue
                rows.append((float(row[0]), float(row[1]), g))
            except (ValueError, IndexError):
                continue
    return rows, x_col, y_col


def _find_local_minima(
    G: list[list[float]],
    xs: list[float],
    ys: list[float],
    nx: int,
    ny: int,
) -> list[FELMinimum]:
    """Return all finite bins strictly lower than all finite 8-neighbors."""
    NEIGHBORS = [(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]
    result: list[FELMinimum] = []
    for ix in range(nx):
        for iy in range(ny):
            g = G[ix][iy]
            if g != g:  # NaN
                continue
            is_local_min = True
            for dx, dy in NEIGHBORS:
                nx2, ny2 = ix + dx, iy + dy
                if 0 <= nx2 < nx and 0 <= ny2 < ny:
                    g2 = G[nx2][ny2]
                    if g2 == g2 and g2 <= g:
                        is_local_min = False
                        break
            if is_local_min:
                result.append(FELMinimum(
                    minimum_id=0,
                    x_bin=ix, y_bin=iy,
                    x=xs[ix], y=ys[iy],
                    free_energy_kJ_mol=g,
                ))
    return result


def _global_minimum(
    G: list[list[float]],
    xs: list[float],
    ys: list[float],
    nx: int,
    ny: int,
) -> Optional[FELMinimum]:
    best_g = float("inf")
    best: Optional[tuple[int, int]] = None
    for ix in range(nx):
        for iy in range(ny):
            g = G[ix][iy]
            if g == g and g < best_g:
                best_g = g
                best = (ix, iy)
    if best is None:
        return None
    ix, iy = best
    return FELMinimum(minimum_id=0, x_bin=ix, y_bin=iy,
                      x=xs[ix], y=ys[iy], free_energy_kJ_mol=best_g)


def _deduplicate(
    candidates: list[FELMinimum],
    min_distance_bins: float,
    n_minima: int,
) -> list[FELMinimum]:
    """Greedily select minima separated by at least min_distance_bins."""
    selected: list[FELMinimum] = []
    for cand in candidates:
        if len(selected) >= n_minima:
            break
        too_close = any(
            math.sqrt((cand.x_bin - s.x_bin) ** 2 + (cand.y_bin - s.y_bin) ** 2)
            < min_distance_bins
            for s in selected
        )
        if not too_close:
            selected.append(cand)
    return selected
