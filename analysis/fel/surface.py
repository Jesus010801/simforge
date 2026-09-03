"""2D FEL surface computation.

Uses NumPy when available; falls back to a pure-Python implementation
for small datasets and tests.

Formula
-------
    P(x, y)  = N(x, y) / Σ N
    G(x, y)  = -k_B · T · ln(P(x, y))   [kJ/mol]
    G_shifted = G - min(G_finite)          shifts minimum to 0

Constants
---------
    k_B = 0.00831446261815324 kJ mol⁻¹ K⁻¹   (NIST 2018 CODATA)
"""
from __future__ import annotations

import math
from typing import Optional

from analysis.fel.models import FELSurface, FELWarning

KB_KJ_MOL_K: float = 0.00831446261815324
_NAN = float("nan")

_FORMULA = (
    "P(x,y) = N(x,y) / sum(N); "
    "G(x,y) = -k_B * T * ln(P(x,y)) [kJ/mol]; "
    "G shifted so min(G_finite) = 0; "
    f"k_B = {KB_KJ_MOL_K} kJ mol^-1 K^-1"
)


# ── Public entry point ────────────────────────────────────────────────────────

def compute_fel_surface(
    x: list[float],
    y: list[float],
    *,
    n_bins_x: int = 50,
    n_bins_y: int = 50,
    temperature_K: float = 300.0,
    empty_bin_policy: str = "nan",  # "nan" | "cap"
    energy_cap_kj: float = 50.0,
    x_label: str = "x",
    y_label: str = "y",
    x_unit: Optional[str] = None,
    y_unit: Optional[str] = None,
) -> FELSurface:
    """Compute a 2D free-energy surface from two equal-length coordinate arrays.

    Parameters
    ----------
    x, y:
        Coordinate arrays (equal length, ≥ 2 values each).
    n_bins_x, n_bins_y:
        Number of histogram bins along each axis.
    temperature_K:
        Simulation temperature in Kelvin.
    empty_bin_policy:
        ``"nan"``  → empty bins are NaN (default, recommended).
        ``"cap"``  → empty bins receive energy_cap_kj.
    energy_cap_kj:
        Energy value assigned to empty bins when ``empty_bin_policy="cap"``.
    """
    if len(x) != len(y):
        raise ValueError(f"x and y must have equal length ({len(x)} vs {len(y)})")
    if len(x) < 2:
        raise ValueError("At least 2 data points are required to compute a FEL surface")

    warnings: list[FELWarning] = []

    try:
        import numpy as np
        counts_arr, x_edges, y_edges = np.histogram2d(
            x, y, bins=[n_bins_x, n_bins_y]
        )
        x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
        y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])

        total = counts_arr.sum()
        prob_arr = counts_arr / total

        empty_mask = prob_arr == 0
        n_empty = int(empty_mask.sum())

        # G = -kB * T * ln(P)
        with np.errstate(divide="ignore", invalid="ignore"):
            G_arr = -KB_KJ_MOL_K * temperature_K * np.log(
                np.where(empty_mask, np.nan, prob_arr)
            )

        # Shift so minimum finite G is 0
        G_min = float(np.nanmin(G_arr))
        G_arr = G_arr - G_min

        if empty_bin_policy == "cap":
            G_arr = np.where(np.isnan(G_arr), energy_cap_kj, G_arr)

        counts_list  = counts_arr.tolist()
        prob_list    = prob_arr.tolist()
        G_list       = [[G_arr[ix, iy] for iy in range(n_bins_y)] for ix in range(n_bins_x)]
        x_edges_list = x_edges.tolist()
        y_edges_list = y_edges.tolist()
        x_cent_list  = x_centers.tolist()
        y_cent_list  = y_centers.tolist()

    except ImportError:
        (counts_list, x_edges_list, y_edges_list,
         x_cent_list, y_cent_list, n_empty) = _histogram2d_pure(x, y, n_bins_x, n_bins_y)
        prob_list, G_list = _compute_prob_and_G(
            counts_list, temperature_K, empty_bin_policy, energy_cap_kj
        )

    if n_empty > 0:
        pct = 100.0 * n_empty / (n_bins_x * n_bins_y)
        warnings.append(FELWarning(
            "empty_bins",
            f"{n_empty} empty bin(s) ({pct:.1f}%) — "
            f"policy: {empty_bin_policy}"
            + (f", cap at {energy_cap_kj} kJ/mol" if empty_bin_policy == "cap" else ""),
            severity="info",
        ))

    return FELSurface(
        x_edges      = x_edges_list,
        y_edges      = y_edges_list,
        x_centers    = x_cent_list,
        y_centers    = y_cent_list,
        counts       = counts_list,
        probability  = prob_list,
        free_energy  = G_list,
        n_frames_used = len(x),
        temperature_K = temperature_K,
        kB_kj_mol_K  = KB_KJ_MOL_K,
        x_label      = x_label,
        y_label      = y_label,
        formula      = _FORMULA,
        warnings     = warnings,
        x_unit       = x_unit,
        y_unit       = y_unit,
    )


# ── Pure-Python fallback ──────────────────────────────────────────────────────

def _histogram2d_pure(
    x: list[float],
    y: list[float],
    nx: int,
    ny: int,
) -> tuple[list[list[float]], list[float], list[float], list[float], list[float], int]:
    """Minimal pure-Python 2D histogram."""
    x_min, x_max = min(x), max(x)
    y_min, y_max = min(y), max(y)

    # Avoid zero-width bins when all values are equal
    x_range = x_max - x_min or 1.0
    y_range = y_max - y_min or 1.0

    x_step = x_range / nx
    y_step = y_range / ny

    x_edges = [x_min + i * x_step for i in range(nx + 1)]
    y_edges = [y_min + i * y_step for i in range(ny + 1)]
    x_cents = [x_edges[i] + x_step / 2 for i in range(nx)]
    y_cents = [y_edges[i] + y_step / 2 for i in range(ny)]

    counts: list[list[float]] = [[0.0] * ny for _ in range(nx)]
    for xi, yi in zip(x, y):
        ix = min(int((xi - x_min) / x_step), nx - 1)
        iy = min(int((yi - y_min) / y_step), ny - 1)
        counts[ix][iy] += 1.0

    n_empty = sum(1 for row in counts for v in row if v == 0)
    return counts, x_edges, y_edges, x_cents, y_cents, n_empty


def _compute_prob_and_G(
    counts: list[list[float]],
    T: float,
    policy: str,
    cap: float,
) -> tuple[list[list[float]], list[list[float]]]:
    total = sum(v for row in counts for v in row)
    if total == 0:
        raise ValueError("All histogram bins are empty — no data to compute FEL")

    prob: list[list[float]] = []
    G_raw: list[list[float]] = []

    for row in counts:
        p_row: list[float] = []
        g_row: list[float] = []
        for c in row:
            if c == 0:
                p_row.append(0.0)
                g_row.append(_NAN)
            else:
                p = c / total
                p_row.append(p)
                g_row.append(-KB_KJ_MOL_K * T * math.log(p))
        prob.append(p_row)
        G_raw.append(g_row)

    # Shift so min finite G = 0
    G_min = min(
        v for row in G_raw for v in row
        if v == v  # nan check: nan != nan is True
    )
    G: list[list[float]] = []
    for row in G_raw:
        new_row: list[float] = []
        for v in row:
            if v != v:  # NaN
                new_row.append(cap if policy == "cap" else _NAN)
            else:
                new_row.append(v - G_min)
        G.append(new_row)

    return prob, G
