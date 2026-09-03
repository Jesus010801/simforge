"""FEL visualization: 2D contour, static 3D surface, optional interactive 3D HTML,
and minima overlay plots."""
from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Optional, TYPE_CHECKING

from analysis.fel.models import FELRunResult, FELSurface, FELWarning

if TYPE_CHECKING:
    from analysis.fel.models import FELMinimum, FrameMapping


def _axis_label(label: str, unit: Optional[str]) -> str:
    """Combine a feature label with its unit, e.g. 'RMSD' + 'nm' → 'RMSD (nm)'."""
    return f"{label} ({unit})" if unit else label


def render_plots(result: FELRunResult, out_dir: Path) -> list[Path]:
    """Render FEL plots into <out_dir>/plots/. Returns list of created paths."""
    if result.surface is None:
        return []

    plots_dir = out_dir / "plots"
    plots_dir.mkdir(exist_ok=True)
    generated: list[Path] = []

    try:
        import matplotlib  # noqa: F401
    except ImportError:
        result.warnings.append(FELWarning(
            "no_matplotlib",
            "matplotlib not available — 2D/3D plots skipped.",
            "info",
        ))
    else:
        for fn in (_plot_2d_contour, _plot_3d_surface):
            p = fn(result.surface, plots_dir)
            if p:
                generated.append(p)

    p = _plot_interactive(result, plots_dir)
    if p:
        generated.append(p)

    return generated


# ── Internal helpers ──────────────────────────────────────────────────────────

def _build_grid(surf: FELSurface):
    """Return (X, Y, Z) numpy arrays, shape (ny, nx), for contourf/surface.

    free_energy is stored [ix][iy]; here we transpose so Z[iy, ix] maps
    correctly onto X = x_centers (columns) and Y = y_centers (rows).
    """
    import numpy as np
    nx = len(surf.x_centers)
    ny = len(surf.y_centers)
    X, Y = np.meshgrid(surf.x_centers, surf.y_centers)
    Z = np.array(
        [[surf.free_energy[ix][iy] for ix in range(nx)] for iy in range(ny)],
        dtype=float,
    )
    return X, Y, Z


def _plot_2d_contour(surf: FELSurface, plots_dir: Path) -> Optional[Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    X, Y, Z = _build_grid(surf)

    fig, ax = plt.subplots(figsize=(7, 6))
    cf = ax.contourf(X, Y, Z, levels=20, cmap="RdYlBu_r")
    cbar = fig.colorbar(cf, ax=ax)
    cbar.set_label("ΔG (kJ/mol)")
    ax.set_xlabel(_axis_label(surf.x_label, surf.x_unit))
    ax.set_ylabel(_axis_label(surf.y_label, surf.y_unit))
    ax.set_title("Free Energy Landscape")

    path = plots_dir / "fel_contour.png"
    fig.savefig(str(path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def _plot_3d_surface(surf: FELSurface, plots_dir: Path) -> Optional[Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    X, Y, Z = _build_grid(surf)
    Z_masked = np.ma.masked_invalid(Z)

    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(X, Y, Z_masked, cmap="RdYlBu_r", edgecolor="none", alpha=0.9)
    ax.set_xlabel(_axis_label(surf.x_label, surf.x_unit))
    ax.set_ylabel(_axis_label(surf.y_label, surf.y_unit))
    ax.set_zlabel("ΔG (kJ/mol)")
    ax.set_title("Free Energy Landscape")

    path = plots_dir / "fel_surface_3d.png"
    fig.savefig(str(path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def _plot_interactive(result: FELRunResult, plots_dir: Path) -> Optional[Path]:
    try:
        import plotly.graph_objects as go
    except ImportError:
        result.warnings.append(FELWarning(
            "no_plotly",
            "Interactive 3D plot was skipped because Plotly is not installed.",
            "info",
        ))
        return None

    surf = result.surface
    _, _, Z = _build_grid(surf)

    fig = go.Figure(data=[go.Surface(
        x=surf.x_centers,
        y=surf.y_centers,
        z=Z.tolist(),
        colorscale="RdYlBu",
        colorbar=dict(title="ΔG (kJ/mol)"),
        hovertemplate=(
            f"{_axis_label(surf.x_label, surf.x_unit)}: %{{x:.3f}}<br>"
            f"{_axis_label(surf.y_label, surf.y_unit)}: %{{y:.3f}}<br>"
            "ΔG: %{z:.2f} kJ/mol<extra></extra>"
        ),
    )])
    fig.update_layout(
        title="Free Energy Landscape",
        scene=dict(
            xaxis_title=_axis_label(surf.x_label, surf.x_unit),
            yaxis_title=_axis_label(surf.y_label, surf.y_unit),
            zaxis_title="ΔG (kJ/mol)",
        ),
        margin=dict(l=0, r=0, t=40, b=0),
    )

    path = plots_dir / "fel_surface_3d.html"
    fig.write_html(str(path))
    return path


# ── Minima overlay plots ──────────────────────────────────────────────────────

def render_minima_plots(
    surface_csv: Path,
    minima: list,
    frame_map: list,
    x_label: str,
    y_label: str,
    x_unit: Optional[str],
    y_unit: Optional[str],
    plots_dir: Path,
    warnings_out: list,
) -> list[Path]:
    """Overlay minima markers on the FEL contour and generate 3D interactive plot.

    Reads the flat surface CSV; produces:
      - ``fel_minima_overlay.png``  (always if matplotlib available)
      - ``fel_minima_3d.html``      (only if Plotly available)
    """
    generated: list[Path] = []

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return generated

    grid_data = _read_surface_for_overlay(surface_csv)
    if grid_data is None:
        return generated

    xs, ys, Z = grid_data
    X, Y = np.meshgrid(xs, ys)

    xl = _axis_label(x_label, x_unit)
    yl = _axis_label(y_label, y_unit)

    fig, ax = plt.subplots(figsize=(7, 6))
    cf = ax.contourf(X, Y, Z, levels=20, cmap="RdYlBu_r")
    cbar = fig.colorbar(cf, ax=ax)
    cbar.set_label("ΔG (kJ/mol)")
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)
    ax.set_title("Free Energy Landscape — Minima")

    _MARKER_COLORS = ["white", "yellow", "cyan", "lime", "orange"]
    for m in minima:
        color = _MARKER_COLORS[(m.minimum_id - 1) % len(_MARKER_COLORS)]
        ax.scatter(m.x, m.y, color=color, edgecolors="black", s=120, zorder=5)
        ax.annotate(
            f"M{m.minimum_id}",
            (m.x, m.y),
            textcoords="offset points", xytext=(6, 4),
            fontsize=9, fontweight="bold", color="white",
            path_effects=_text_outline(),
        )

    path = plots_dir / "fel_minima_overlay.png"
    fig.savefig(str(path), dpi=150, bbox_inches="tight")
    plt.close(fig)
    generated.append(path)

    # Interactive Plotly version with minima
    try:
        import plotly.graph_objects as go
    except ImportError:
        warnings_out.append(FELWarning(
            "no_plotly_minima",
            "Interactive minima 3D plot skipped — Plotly not installed.",
            "info",
        ))
        return generated

    traces = [go.Surface(
        x=xs, y=ys, z=Z.tolist(),
        colorscale="RdYlBu",
        colorbar=dict(title="ΔG (kJ/mol)"),
        opacity=0.85,
        showscale=True,
    )]
    for m in minima:
        traces.append(go.Scatter3d(
            x=[m.x], y=[m.y], z=[m.free_energy_kJ_mol],
            mode="markers+text",
            text=[f"M{m.minimum_id}"],
            textposition="top center",
            marker=dict(size=8, color="white", line=dict(color="black", width=2)),
            hovertemplate=(
                f"M{m.minimum_id}<br>"
                f"{xl}: %{{x:.3f}}<br>"
                f"{yl}: %{{y:.3f}}<br>"
                "ΔG: %{z:.2f} kJ/mol<extra></extra>"
            ),
        ))

    pfig = go.Figure(data=traces)
    pfig.update_layout(
        title="Free Energy Landscape — Minima",
        scene=dict(
            xaxis_title=xl,
            yaxis_title=yl,
            zaxis_title="ΔG (kJ/mol)",
        ),
        margin=dict(l=0, r=0, t=40, b=0),
    )
    html_path = plots_dir / "fel_minima_3d.html"
    pfig.write_html(str(html_path))
    generated.append(html_path)

    return generated


def _read_surface_for_overlay(
    surface_csv: Path,
) -> Optional[tuple]:
    """Read flat FEL CSV and return (x_centers, y_centers, Z_array).

    Z has shape (ny, nx) for use with meshgrid(xs, ys).
    """
    try:
        import numpy as np
    except ImportError:
        return None

    if not surface_csv.exists():
        return None

    xs_set: set[float] = set()
    ys_set: set[float] = set()
    data: dict[tuple[float, float], float] = {}

    with surface_csv.open(newline="") as fh:
        reader = csv.reader(fh)
        next(reader)  # header
        for row in reader:
            if len(row) < 3:
                continue
            try:
                x, y, g = float(row[0]), float(row[1]), float(row[2])
            except (ValueError, IndexError):
                continue
            if g != g:  # skip NaN
                continue
            xs_set.add(x)
            ys_set.add(y)
            data[(x, y)] = g

    if not data:
        return None

    xs = sorted(xs_set)
    ys = sorted(ys_set)
    nx, ny = len(xs), len(ys)
    Z = np.full((ny, nx), float("nan"))
    for i, y in enumerate(ys):
        for j, x in enumerate(xs):
            if (x, y) in data:
                Z[i, j] = data[(x, y)]

    return xs, ys, Z


def _text_outline():
    """Return matplotlib path effects for readable text on busy backgrounds."""
    try:
        import matplotlib.patheffects as pe
        return [pe.withStroke(linewidth=2, foreground="black")]
    except ImportError:
        return []
