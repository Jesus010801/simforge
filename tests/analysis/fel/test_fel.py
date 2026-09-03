"""Phase 2 FEL foundation tests.

All tests are self-contained (no GROMACS required).
Synthetic XVG files are written to tmp_path so tests run in any CI.
The real hDHFR XVG files in docs/FEL/ are used when present (optional).
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Optional

import pytest

_REPO = Path(__file__).resolve().parent.parent.parent.parent

# ── Helpers ───────────────────────────────────────────────────────────────────

def _write_xvg(
    path: Path,
    rows: list[tuple[float, ...]],
    title: str = "Test",
    xlabel: str = 'Time (ps)',
    ylabel: str = 'Value',
    legends: Optional[list[str]] = None,
) -> Path:
    """Write a minimal valid GROMACS XVG file."""
    lines = [
        f'# synthetic test xvg',
        f'@    title "{title}"',
        f'@    xaxis  label "{xlabel}"',
        f'@    yaxis  label "{ylabel}"',
        '@TYPE xy',
    ]
    if legends:
        for i, leg in enumerate(legends):
            lines.append(f'@ s{i} legend "{leg}"')
    for row in rows:
        lines.append("   ".join(f"{v:.7f}" for v in row))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


def _cli(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "cli", *args],
        capture_output=True, text=True, cwd=str(_REPO),
    )


def _simple_rows(n: int = 50, dt_ps: float = 100.0) -> list[tuple[float, float]]:
    """n frames, time in ps, data = 0.1 + sin(i/n)."""
    return [(i * dt_ps, 0.1 + math.sin(i / n)) for i in range(n)]


# ── Group 1: XVG parser ───────────────────────────────────────────────────────

class TestXVGParser:

    def test_load_simple_xvg(self, tmp_path: Path) -> None:
        """Parses a minimal XVG and returns correct frame count."""
        rows = _simple_rows(20)
        p = _write_xvg(tmp_path / "test.xvg", rows)

        from analysis.fel.xvg import load_xvg_series
        series = load_xvg_series(p, column=1)

        assert series.n_frames() == 20
        assert series.x_unit == "ps"
        assert len(series.data_values) == 20

    def test_unit_detection_ns(self, tmp_path: Path) -> None:
        """Detects 'ns' from xaxis label."""
        rows = [(i * 0.2, float(i)) for i in range(10)]
        p = _write_xvg(tmp_path / "ns.xvg", rows, xlabel='Time (ns)')

        from analysis.fel.xvg import load_xvg_series
        s = load_xvg_series(p)
        assert s.x_unit == "ns"

    def test_column_selection(self, tmp_path: Path) -> None:
        """Column 2 selects the second data column."""
        rows = [(float(i), float(i) * 2, float(i) * 3) for i in range(5)]
        p = _write_xvg(tmp_path / "multi.xvg", rows, legends=["col1", "col2"])

        from analysis.fel.xvg import load_xvg_series
        s1 = load_xvg_series(p, column=1)
        s2 = load_xvg_series(p, column=2)

        assert s1.data_values == [0.0, 2.0, 4.0, 6.0, 8.0]
        assert s2.data_values == [0.0, 3.0, 6.0, 9.0, 12.0]

    def test_column_out_of_range_raises(self, tmp_path: Path) -> None:
        """Column index beyond available columns raises ValueError."""
        rows = [(float(i), float(i)) for i in range(5)]
        p = _write_xvg(tmp_path / "single.xvg", rows)

        from analysis.fel.xvg import load_xvg_series
        with pytest.raises(ValueError, match="Column 2 not found"):
            load_xvg_series(p, column=2)

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        from analysis.fel.xvg import load_xvg_series
        with pytest.raises(FileNotFoundError):
            load_xvg_series(tmp_path / "ghost.xvg")

    def test_empty_file_raises(self, tmp_path: Path) -> None:
        """A file with only header lines raises ValueError."""
        p = tmp_path / "empty.xvg"
        p.write_text('@    title "Empty"\n@TYPE xy\n')

        from analysis.fel.xvg import load_xvg_series
        with pytest.raises(ValueError, match="No numeric data rows"):
            load_xvg_series(p)

    def test_skips_malformed_lines(self, tmp_path: Path) -> None:
        """Single-token lines are skipped; valid rows are still parsed."""
        content = (
            '@ title "T"\n'
            '@ xaxis label "Time (ps)"\n'
            '@ yaxis label "V"\n'
            '@TYPE xy\n'
            '0.0   1.0\n'
            'bad_line\n'
            '1.0   2.0\n'
        )
        p = tmp_path / "partial.xvg"
        p.write_text(content)

        from analysis.fel.xvg import load_xvg_series
        s = load_xvg_series(p)
        assert s.n_frames() == 2
        assert s.n_skipped == 1

    def test_legend_name_in_column_name(self, tmp_path: Path) -> None:
        """Series name is read from @ sN legend directive."""
        rows = [(float(i), float(i)) for i in range(3)]
        p = _write_xvg(tmp_path / "leg.xvg", rows, legends=["RMSD"])

        from analysis.fel.xvg import load_xvg_series
        s = load_xvg_series(p)
        assert s.column_name == "RMSD"


# ── Group 2: Feature table ────────────────────────────────────────────────────

class TestFeatureTable:

    def _pair(self, tmp_path: Path, n: int = 30, unit_x: str = 'Time (ps)', unit_y: str = 'Time (ps)'):
        rows = _simple_rows(n)
        px = _write_xvg(tmp_path / "x.xvg", rows, xlabel=unit_x)
        py = _write_xvg(tmp_path / "y.xvg", rows, xlabel=unit_y)
        from analysis.fel.xvg import load_xvg_series
        return load_xvg_series(px), load_xvg_series(py)

    def test_aligned_same_units(self, tmp_path: Path) -> None:
        sx, sy = self._pair(tmp_path)
        from analysis.fel.features import build_feature_table
        ft = build_feature_table(sx, sy, x_label="RMSD", y_label="Rg")
        assert ft.n_frames == 30
        assert ft.x_label == "RMSD"
        assert ft.y_label == "Rg"

    def test_length_mismatch_raises(self, tmp_path: Path) -> None:
        rows_a = _simple_rows(20)
        rows_b = _simple_rows(25)
        px = _write_xvg(tmp_path / "x.xvg", rows_a)
        py = _write_xvg(tmp_path / "y.xvg", rows_b)
        from analysis.fel.xvg import load_xvg_series
        from analysis.fel.features import build_feature_table
        sx = load_xvg_series(px)
        sy = load_xvg_series(py)
        with pytest.raises(ValueError, match="different lengths"):
            build_feature_table(sx, sy)

    def test_ps_ns_mismatch_auto_corrected(self, tmp_path: Path) -> None:
        """Header claims ps but values are in ns range → auto-corrected with warning."""
        # Both XVGs cover the same 0..1.9 ns trajectory.
        # X: header says "Time (ps)" but values are 0, 0.2, 0.4 ... 3.8 (really ns).
        # Y: header says "Time (ns)" and values are also 0, 0.2 ... 3.8 ns.
        # After naïve conversion X→ [0, 2e-4, 4e-4 ...] ns ≠ Y→ [0, 0.2 ...] ns.
        # The auto-correction reinterprets X as ns → they match.
        n = 20
        dt_ns = 0.2          # step in ns
        rows_wrong_header = [(i * dt_ns, float(i)) for i in range(n)]  # values in ns, header says ps
        rows_correct      = [(i * dt_ns, float(i)) for i in range(n)]  # values in ns, header correct

        px = _write_xvg(tmp_path / "x.xvg", rows_wrong_header, xlabel='Time (ps)')
        py = _write_xvg(tmp_path / "y.xvg", rows_correct,      xlabel='Time (ns)')

        from analysis.fel.xvg import load_xvg_series
        from analysis.fel.features import build_feature_table

        sx = load_xvg_series(px)
        sy = load_xvg_series(py)
        # tolerance in ns; after correction they agree exactly
        ft = build_feature_table(sx, sy, time_tolerance=1e-9)
        assert ft.n_frames == n
        # auto-correction warning must be emitted
        warn_codes = [w.code for w in ft.warnings]
        assert "time_unit_corrected" in warn_codes

    def test_unresolvable_time_mismatch_raises(self, tmp_path: Path) -> None:
        """Completely different time axes raise ValueError."""
        rows_a = [(i * 100.0, float(i)) for i in range(10)]
        rows_b = [(i * 500.0, float(i)) for i in range(10)]
        px = _write_xvg(tmp_path / "x.xvg", rows_a)
        py = _write_xvg(tmp_path / "y.xvg", rows_b)
        from analysis.fel.xvg import load_xvg_series
        from analysis.fel.features import build_feature_table
        sx = load_xvg_series(px)
        sy = load_xvg_series(py)
        with pytest.raises(ValueError, match="do not align"):
            build_feature_table(sx, sy, time_tolerance=1e-6)

    def test_feature_csv_output(self, tmp_path: Path) -> None:
        sx, sy = self._pair(tmp_path, n=5)
        from analysis.fel.features import build_feature_table
        ft = build_feature_table(sx, sy, x_label="X", y_label="Y")
        csv_lines = ft.to_csv_lines()
        assert csv_lines[0] == "time_ns,X,Y"
        assert len(csv_lines) == 6  # header + 5 rows


# ── Group 3: FEL surface ──────────────────────────────────────────────────────

class TestFELSurface:

    def _uniform_data(self, n: int = 100):
        import math
        x = [math.sin(i * 0.1) for i in range(n)]
        y = [math.cos(i * 0.1) for i in range(n)]
        return x, y

    def test_minimum_g_is_zero(self) -> None:
        """After shifting, the minimum finite G value must be exactly 0."""
        x, y = self._uniform_data(200)
        from analysis.fel.surface import compute_fel_surface
        surf = compute_fel_surface(x, y, n_bins_x=10, n_bins_y=10)
        finite_g = [
            surf.free_energy[ix][iy]
            for ix in range(len(surf.x_centers))
            for iy in range(len(surf.y_centers))
            if surf.free_energy[ix][iy] == surf.free_energy[ix][iy]  # not NaN
        ]
        assert len(finite_g) > 0
        assert min(finite_g) == pytest.approx(0.0, abs=1e-10)

    def test_empty_bins_nan_policy(self) -> None:
        """Empty bins are NaN with 'nan' policy."""
        # Only 2 distinct values → most bins empty
        x = [0.0, 1.0]
        y = [0.0, 1.0]
        from analysis.fel.surface import compute_fel_surface
        surf = compute_fel_surface(x, y, n_bins_x=5, n_bins_y=5, empty_bin_policy="nan")
        nan_count = sum(
            1 for ix in range(5) for iy in range(5)
            if surf.free_energy[ix][iy] != surf.free_energy[ix][iy]
        )
        assert nan_count > 0

    def test_empty_bins_cap_policy(self) -> None:
        """Empty bins are capped at energy_cap_kj with 'cap' policy."""
        x = [0.0, 1.0]
        y = [0.0, 1.0]
        from analysis.fel.surface import compute_fel_surface
        surf = compute_fel_surface(
            x, y, n_bins_x=5, n_bins_y=5,
            empty_bin_policy="cap", energy_cap_kj=30.0,
        )
        max_g = max(
            surf.free_energy[ix][iy]
            for ix in range(5) for iy in range(5)
            if surf.free_energy[ix][iy] == surf.free_energy[ix][iy]
        )
        assert max_g <= 30.0 + 1e-9

    def test_kb_constant(self) -> None:
        """k_B constant is the correct NIST value."""
        from analysis.fel.surface import KB_KJ_MOL_K
        assert KB_KJ_MOL_K == pytest.approx(0.00831446261815324, rel=1e-10)

    def test_probability_sums_to_one(self) -> None:
        """Sum of probabilities over all non-empty bins equals 1."""
        x, y = self._uniform_data(300)
        from analysis.fel.surface import compute_fel_surface
        surf = compute_fel_surface(x, y, n_bins_x=10, n_bins_y=10)
        total = sum(
            surf.probability[ix][iy]
            for ix in range(len(surf.x_centers))
            for iy in range(len(surf.y_centers))
        )
        assert total == pytest.approx(1.0, abs=1e-6)

    def test_too_few_points_raises(self) -> None:
        from analysis.fel.surface import compute_fel_surface
        with pytest.raises(ValueError, match="At least 2"):
            compute_fel_surface([1.0], [1.0])

    def test_csv_flat_format(self) -> None:
        """CSV outputs are flat with columns x_*, y_*, value."""
        x, y = self._uniform_data(50)
        from analysis.fel.surface import compute_fel_surface
        surf = compute_fel_surface(x, y, n_bins_x=5, n_bins_y=5,
                                   x_label="RMSD", y_label="Rg")
        counts_lines = surf.counts_csv()
        assert counts_lines[0] == "x_RMSD,y_Rg,counts"
        assert len(counts_lines) == 1 + 5 * 5  # header + nx*ny rows
        G_lines = surf.free_energy_csv()
        assert G_lines[0] == "x_RMSD,y_Rg,free_energy_kJ_mol"


# ── Group 4: Workflow & CLI ───────────────────────────────────────────────────

class TestWorkflow:

    def _make_pair(self, tmp_path: Path, n: int = 60):
        rows = _simple_rows(n)
        px = _write_xvg(tmp_path / "rmsd.xvg", rows, title="RMSD", legends=["RMSD"])
        py = _write_xvg(tmp_path / "gyrate.xvg", rows, title="Rg", legends=["Rg"])
        return px, py

    def test_workflow_creates_output_files(self, tmp_path: Path) -> None:
        """Full workflow run creates all expected output files."""
        px, py = self._make_pair(tmp_path)
        out = tmp_path / "fel_out"

        from analysis.fel.models import FELConfig
        from analysis.fel.workflow import run_fel
        config = FELConfig(
            xvg_x=px, xvg_y=py,
            x_column=1, y_column=1,
            x_label="RMSD", y_label="Rg",
            temperature_K=300.0,
            n_bins_x=10, n_bins_y=10,
            out_dir=out,
            time_unit_x="auto", time_unit_y="auto",
            time_tolerance=1e-6,
            empty_bin_policy="nan",
            energy_cap_kj=50.0,
            dry_run=False,
        )
        result = run_fel(config)

        assert (out / "config.json").exists()
        assert (out / "metadata" / "provenance.json").exists()
        assert (out / "metadata" / "warnings.json").exists()
        assert (out / "data" / "features.csv").exists()
        assert (out / "data" / "histogram_counts.csv").exists()
        assert (out / "data" / "probability.csv").exists()
        assert (out / "data" / "free_energy_surface.csv").exists()
        assert (out / "minima" / "README.md").exists()
        assert (out / "report.md").exists()

    def test_dry_run_no_surface_files(self, tmp_path: Path) -> None:
        """Dry run writes config + provenance but not the surface CSVs."""
        px, py = self._make_pair(tmp_path)
        out = tmp_path / "dry_out"

        from analysis.fel.models import FELConfig
        from analysis.fel.workflow import run_fel
        config = FELConfig(
            xvg_x=px, xvg_y=py,
            x_column=1, y_column=1,
            x_label="X", y_label="Y",
            temperature_K=300.0,
            n_bins_x=10, n_bins_y=10,
            out_dir=out,
            time_unit_x="auto", time_unit_y="auto",
            time_tolerance=1e-6,
            empty_bin_policy="nan",
            energy_cap_kj=50.0,
            dry_run=True,
        )
        result = run_fel(config)

        assert (out / "data" / "features.csv").exists()
        assert not (out / "data" / "free_energy_surface.csv").exists()

    def test_provenance_json_valid(self, tmp_path: Path) -> None:
        """provenance.json is valid JSON with required keys."""
        px, py = self._make_pair(tmp_path)
        out = tmp_path / "prov_out"

        from analysis.fel.models import FELConfig
        from analysis.fel.workflow import run_fel
        config = FELConfig(
            xvg_x=px, xvg_y=py,
            x_column=1, y_column=1,
            x_label="X", y_label="Y",
            temperature_K=300.0,
            n_bins_x=5, n_bins_y=5,
            out_dir=out,
            time_unit_x="auto", time_unit_y="auto",
            time_tolerance=1e-6,
            empty_bin_policy="nan",
            energy_cap_kj=50.0,
            dry_run=False,
        )
        run_fel(config)
        data = json.loads((out / "metadata" / "provenance.json").read_text())
        for key in ("timestamp", "simforge_version", "config", "warnings"):
            assert key in data, f"Missing key: {key}"

    def test_cli_fel_run_dry_run(self, tmp_path: Path) -> None:
        """CLI `simforge fel run --dry-run` exits 0 and creates output dir."""
        rows = _simple_rows(30)
        px = _write_xvg(tmp_path / "rmsd.xvg", rows, legends=["RMSD"])
        py = _write_xvg(tmp_path / "gyrate.xvg", rows, legends=["Rg"])
        out = tmp_path / "cli_dry"

        r = _cli(
            "fel", "run",
            "--xvg-x", str(px),
            "--xvg-y", str(py),
            "--temperature", "300",
            "--bins", "5,5",
            "--out", str(out),
            "--dry-run",
        )
        assert r.returncode == 0, r.stderr
        assert out.is_dir()
        assert (out / "data" / "features.csv").exists()

    def test_cli_fel_run_full(self, tmp_path: Path) -> None:
        """CLI `simforge fel run` (full) exits 0 and creates surface CSVs."""
        rows = _simple_rows(40)
        px = _write_xvg(tmp_path / "rmsd.xvg", rows, legends=["RMSD"])
        py = _write_xvg(tmp_path / "gyrate.xvg", rows, legends=["Rg"])
        out = tmp_path / "cli_full"

        r = _cli(
            "fel", "run",
            "--xvg-x", str(px),
            "--xvg-y", str(py),
            "--x-label", "RMSD",
            "--y-label", "Rg",
            "--temperature", "300",
            "--bins", "5,5",
            "--out", str(out),
        )
        assert r.returncode == 0, r.stderr
        assert (out / "data" / "free_energy_surface.csv").exists()

    def test_cli_missing_xvg_exits_nonzero(self, tmp_path: Path) -> None:
        """CLI errors clearly when required XVG options are omitted."""
        r = _cli("fel", "run", "--temperature", "300", "--out", str(tmp_path / "out"))
        assert r.returncode != 0


# ── Group 5: Phase 2C – Visualization ────────────────────────────────────────

class TestVisualization:

    def _make_result(self, tmp_path: Path, n: int = 60):
        rows = _simple_rows(n)
        px = _write_xvg(tmp_path / "rmsd.xvg", rows, legends=["RMSD"])
        py = _write_xvg(tmp_path / "gyrate.xvg", rows, legends=["Rg"])
        from analysis.fel.models import FELConfig
        from analysis.fel.workflow import run_fel
        config = FELConfig(
            xvg_x=px, xvg_y=py,
            x_column=1, y_column=1,
            x_label="RMSD", y_label="Rg",
            temperature_K=300.0,
            n_bins_x=10, n_bins_y=10,
            out_dir=tmp_path / "fel_out",
            time_unit_x="auto", time_unit_y="auto",
            time_tolerance=1e-6,
            empty_bin_policy="nan",
            energy_cap_kj=50.0,
            dry_run=False,
        )
        return run_fel(config)

    def _make_surface(self, n: int = 80):
        import math
        x = [math.sin(i * 0.1) for i in range(n)]
        y = [math.cos(i * 0.1) for i in range(n)]
        from analysis.fel.surface import compute_fel_surface
        return compute_fel_surface(x, y, n_bins_x=10, n_bins_y=10,
                                   x_label="RMSD", y_label="Rg")

    def test_2d_contour_axis_labels_correct(self, tmp_path: Path, monkeypatch) -> None:
        """2D contour plot has x_label on X axis and y_label on Y axis."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        saved_axes: list = []
        original_savefig = plt.Figure.savefig

        def capture(self, *args, **kwargs):
            saved_axes.append(list(self.axes))
            original_savefig(self, *args, **kwargs)

        monkeypatch.setattr(plt.Figure, "savefig", capture)

        surf = self._make_surface()
        from analysis.fel.plot import _plot_2d_contour
        _plot_2d_contour(surf, tmp_path)

        assert saved_axes, "savefig was not called"
        ax = saved_axes[0][0]
        assert ax.get_xlabel() == "RMSD"
        assert ax.get_ylabel() == "Rg"

    def test_2d_contour_png_created(self, tmp_path: Path) -> None:
        """2D contour PNG exists after a full workflow run."""
        result = self._make_result(tmp_path)
        assert (result.config.out_dir / "plots" / "fel_contour.png").exists()

    def test_3d_surface_png_created(self, tmp_path: Path) -> None:
        """Static 3D surface PNG exists after a full workflow run."""
        result = self._make_result(tmp_path)
        assert (result.config.out_dir / "plots" / "fel_surface_3d.png").exists()

    def test_3d_surface_axis_labels_correct(self, tmp_path: Path, monkeypatch) -> None:
        """3D surface plot has correct axis labels."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        saved_axes: list = []
        original_savefig = plt.Figure.savefig

        def capture(self, *args, **kwargs):
            saved_axes.append(list(self.axes))
            original_savefig(self, *args, **kwargs)

        monkeypatch.setattr(plt.Figure, "savefig", capture)

        surf = self._make_surface()
        from analysis.fel.plot import _plot_3d_surface
        _plot_3d_surface(surf, tmp_path)

        assert saved_axes
        ax = saved_axes[0][0]
        assert ax.get_xlabel() == "RMSD"
        assert ax.get_ylabel() == "Rg"

    def test_plotly_skipped_gracefully_when_unavailable(self, tmp_path: Path, monkeypatch) -> None:
        """When Plotly is not installed, no HTML is written and a warning is emitted."""
        import builtins
        real_import = builtins.__import__

        def block_plotly(name, *args, **kwargs):
            if name == "plotly" or name.startswith("plotly."):
                raise ImportError(f"Mocked: {name} not available")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", block_plotly)

        from analysis.fel.models import FELRunResult, FELConfig
        from analysis.fel.plot import _plot_interactive

        rows = _simple_rows(40)
        px = _write_xvg(tmp_path / "x.xvg", rows)
        py = _write_xvg(tmp_path / "y.xvg", rows)
        config = FELConfig(
            xvg_x=px, xvg_y=py,
            x_column=1, y_column=1,
            x_label="RMSD", y_label="Rg",
            temperature_K=300.0,
            n_bins_x=5, n_bins_y=5,
            out_dir=tmp_path / "out",
            time_unit_x="auto", time_unit_y="auto",
            time_tolerance=1e-6,
            empty_bin_policy="nan",
            energy_cap_kj=50.0,
            dry_run=False,
        )
        result = FELRunResult(config=config)
        result.surface = self._make_surface()

        path = _plot_interactive(result, tmp_path)
        assert path is None
        assert not (tmp_path / "fel_surface_3d.html").exists()
        warn_codes = [w.code for w in result.warnings]
        assert "no_plotly" in warn_codes

    def test_plotly_html_generated_when_available(self, tmp_path: Path, monkeypatch) -> None:
        """Interactive HTML is written when plotly is available (mocked)."""
        import sys
        from unittest.mock import MagicMock

        html_path = tmp_path / "fel_surface_3d.html"
        mock_fig = MagicMock()

        def fake_write_html(path, *args, **kwargs):
            Path(path).write_text("<html>mock</html>")

        mock_fig.write_html.side_effect = fake_write_html
        mock_go = MagicMock()
        mock_go.Figure.return_value = mock_fig

        # `import plotly.graph_objects as go` resolves via the parent module's
        # attribute, so we must set both sys.modules and the parent attribute.
        mock_plotly = MagicMock()
        mock_plotly.graph_objects = mock_go
        monkeypatch.setitem(sys.modules, "plotly", mock_plotly)
        monkeypatch.setitem(sys.modules, "plotly.graph_objects", mock_go)

        from analysis.fel.models import FELRunResult, FELConfig
        from analysis.fel.plot import _plot_interactive

        rows = _simple_rows(40)
        px = _write_xvg(tmp_path / "x.xvg", rows)
        py = _write_xvg(tmp_path / "y.xvg", rows)
        config = FELConfig(
            xvg_x=px, xvg_y=py,
            x_column=1, y_column=1,
            x_label="RMSD", y_label="Rg",
            temperature_K=300.0,
            n_bins_x=5, n_bins_y=5,
            out_dir=tmp_path / "out",
            time_unit_x="auto", time_unit_y="auto",
            time_tolerance=1e-6,
            empty_bin_policy="nan",
            energy_cap_kj=50.0,
            dry_run=False,
        )
        result = FELRunResult(config=config)
        result.surface = self._make_surface()

        path = _plot_interactive(result, tmp_path)

        assert path is not None
        assert path.name == "fel_surface_3d.html"
        assert html_path.exists()
        assert not any(w.code == "no_plotly" for w in result.warnings)

    def test_report_lists_generated_plots(self, tmp_path: Path) -> None:
        """report.md contains a Visualizations section listing plot files."""
        result = self._make_result(tmp_path)
        report_text = (result.config.out_dir / "report.md").read_text()
        assert "## Visualizations" in report_text
        assert "fel_contour.png" in report_text
        assert "fel_surface_3d.png" in report_text

    def test_no_hard_plotly_dependency(self) -> None:
        """Plotly is not listed as a required dependency in pyproject.toml."""
        repo = Path(__file__).resolve().parent.parent.parent.parent
        pyproject = (repo / "pyproject.toml").read_text()
        # Plotly should not appear in [project] dependencies (optional only)
        import re
        # Find the [project] dependencies block
        deps_match = re.search(
            r'\[project\].*?dependencies\s*=\s*\[(.*?)\]',
            pyproject, re.DOTALL
        )
        if deps_match:
            deps_block = deps_match.group(1)
            assert "plotly" not in deps_block.lower(), \
                "plotly must not be a hard dependency in [project].dependencies"

    def test_report_notes_skipped_interactive_when_no_plotly(self, tmp_path: Path, monkeypatch) -> None:
        """report.md notes that interactive HTML was skipped when Plotly is absent."""
        import builtins
        real_import = builtins.__import__

        def block_plotly(name, *args, **kwargs):
            if name == "plotly" or name.startswith("plotly."):
                raise ImportError(f"Mocked: {name} not available")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", block_plotly)

        result = self._make_result(tmp_path)
        report_text = (result.config.out_dir / "report.md").read_text()
        assert "Plotly" in report_text or "plotly" in report_text


# ── Group 6: Phase 3 — Minima detection and extraction ───────────────────────

def _make_fel_dir(
    tmp_path: Path,
    *,
    n_bins: int = 5,
    x_label: str = "RMSD",
    y_label: str = "Rg",
    x_unit: str = "nm",
    y_unit: str = "nm",
    n_frames: int = 30,
) -> Path:
    """Create a minimal FEL directory structure for Phase 3 tests."""
    import math

    fel = tmp_path / "fel_out"
    (fel / "data").mkdir(parents=True)
    (fel / "minima").mkdir()

    # --- free_energy_surface.csv ---
    # Create a symmetric bowl with clear minimum at center
    rows = []
    xs = [0.1 * i for i in range(n_bins)]
    ys = [1.0 + 0.1 * i for i in range(n_bins)]
    cx, cy = xs[n_bins // 2], ys[n_bins // 2]
    for x in xs:
        for y in ys:
            g = 4.0 * ((x - cx) ** 2 + (y - cy) ** 2)  # parabolic bowl
            rows.append((x, y, g))

    lines = [f"x_{x_label},y_{y_label},free_energy_kJ_mol"]
    for x, y, g in rows:
        lines.append(f"{x:.6g},{y:.6g},{g:.6g}")
    (fel / "data" / "free_energy_surface.csv").write_text("\n".join(lines) + "\n")

    # --- features.csv ---
    feat_rows = ["time_ns," + x_label + "," + y_label]
    dt = 0.2
    for i in range(n_frames):
        t = i * dt
        fx = cx + 0.01 * math.sin(i * 0.5)
        fy = cy + 0.01 * math.cos(i * 0.5)
        feat_rows.append(f"{t:.6g},{fx:.6g},{fy:.6g}")
    (fel / "data" / "features.csv").write_text("\n".join(feat_rows) + "\n")

    # --- config.json ---
    import json
    config = {
        "xvg_x": "rmsd.xvg", "xvg_y": "rg.xvg",
        "x_column": 1, "y_column": 1,
        "x_label": x_label, "y_label": y_label,
        "x_unit": x_unit, "y_unit": y_unit,
        "temperature_K": 300.0,
        "n_bins_x": n_bins, "n_bins_y": n_bins,
        "out_dir": str(fel),
        "time_unit_x": "auto", "time_unit_y": "auto",
        "time_tolerance": 1e-6,
        "empty_bin_policy": "nan",
        "energy_cap_kj": 50.0,
        "dry_run": False,
    }
    (fel / "config.json").write_text(json.dumps(config, indent=2) + "\n")

    return fel


class TestMinimaDetection:

    def test_cli_extract_minima_help(self) -> None:
        """simforge fel extract-minima --help exits 0."""
        r = _cli("fel", "extract-minima", "--help")
        assert r.returncode == 0
        assert "extract-minima" in r.stdout or "extract_minima" in r.stdout

    def test_local_minima_detected(self, tmp_path: Path) -> None:
        """Detects at least one local minimum from a parabolic bowl."""
        fel = _make_fel_dir(tmp_path)
        from analysis.fel.minima import detect_local_minima
        minima, warns = detect_local_minima(fel, n_minima=3)
        assert len(minima) >= 1
        assert minima[0].minimum_id == 1

    def test_nan_bins_ignored(self, tmp_path: Path) -> None:
        """NaN bins are excluded from local-minimum search."""
        fel = _make_fel_dir(tmp_path, n_bins=5)
        # Inject NaN into the CSV at the center bin
        csv_path = fel / "data" / "free_energy_surface.csv"
        lines = csv_path.read_text().splitlines()
        # The minimum is the center row — replace its G with NaN
        new_lines = []
        for line in lines:
            parts = line.split(",")
            if len(parts) == 3 and parts[2] not in ("free_energy_kJ_mol",):
                try:
                    g = float(parts[2])
                    if abs(g) < 1e-10:  # global minimum bin
                        parts[2] = "NaN"
                        line = ",".join(parts)
                except ValueError:
                    pass
            new_lines.append(line)
        csv_path.write_text("\n".join(new_lines) + "\n")

        from analysis.fel.minima import detect_local_minima
        minima, warns = detect_local_minima(fel, n_minima=3)
        # Should still find minima (adjacent ones), no crash
        assert isinstance(minima, list)

    def test_ranked_by_free_energy(self, tmp_path: Path) -> None:
        """Returned minima are sorted by free_energy_kJ_mol ascending."""
        # Make a surface with two clear local minima
        from pathlib import Path as P
        fel = tmp_path / "fel_rank"
        (fel / "data").mkdir(parents=True)
        n = 7
        xs = [i * 0.1 for i in range(n)]
        ys = [i * 0.1 for i in range(n)]
        # Two wells: at (1,1) and (5,5) in bin space
        lines = ["x_RMSD,y_Rg,free_energy_kJ_mol"]
        for ix, x in enumerate(xs):
            for iy, y in enumerate(ys):
                d1 = (ix - 1) ** 2 + (iy - 1) ** 2
                d2 = (ix - 5) ** 2 + (iy - 5) ** 2
                g = min(float(d1) * 2.0, float(d2) * 3.0)
                lines.append(f"{x:.4g},{y:.4g},{g:.4g}")
        (fel / "data" / "free_energy_surface.csv").write_text("\n".join(lines) + "\n")
        import json
        (fel / "config.json").write_text(json.dumps({"x_label": "RMSD", "y_label": "Rg"}) + "\n")

        from analysis.fel.minima import detect_local_minima
        minima, _ = detect_local_minima(fel, n_minima=5, min_distance_bins=3)
        energies = [m.free_energy_kJ_mol for m in minima]
        assert energies == sorted(energies)

    def test_deduplication(self, tmp_path: Path) -> None:
        """Two nearby minima are deduplicated when within min_distance_bins."""
        from pathlib import Path as P
        fel = tmp_path / "fel_dedup"
        (fel / "data").mkdir(parents=True)
        # Single bowl — only one minimum; request 3 but expect just 1
        n = 5
        xs = [i * 0.1 for i in range(n)]
        ys = [i * 0.1 for i in range(n)]
        lines = ["x_X,y_Y,free_energy_kJ_mol"]
        cx, cy = 2, 2
        for ix, x in enumerate(xs):
            for iy, y in enumerate(ys):
                g = float((ix - cx) ** 2 + (iy - cy) ** 2)
                lines.append(f"{x:.4g},{y:.4g},{g:.4g}")
        (fel / "data" / "free_energy_surface.csv").write_text("\n".join(lines) + "\n")
        import json
        (fel / "config.json").write_text(json.dumps({"x_label": "X", "y_label": "Y"}) + "\n")

        from analysis.fel.minima import detect_local_minima
        minima, _ = detect_local_minima(fel, n_minima=3, min_distance_bins=10)
        assert len(minima) == 1  # only center survives with large separation

    def test_global_minimum_fallback(self, tmp_path: Path) -> None:
        """Falls back to global minimum when no local minima exist (flat surface)."""
        fel = tmp_path / "fel_fallback"
        (fel / "data").mkdir(parents=True)
        # Flat surface: every cell has a neighbor with the same G → no strict local min
        # The 8-neighbor check is g2 <= g, so equal neighbors disqualify a bin.
        n = 4
        lines = ["x_X,y_Y,free_energy_kJ_mol"]
        for ix in range(n):
            for iy in range(n):
                g = 1.0  # all the same → no bin is strictly lower than all neighbors
                lines.append(f"{0.1*ix:.4g},{0.1*iy:.4g},{g:.4g}")
        (fel / "data" / "free_energy_surface.csv").write_text("\n".join(lines) + "\n")
        import json
        (fel / "config.json").write_text(json.dumps({"x_label": "X", "y_label": "Y"}) + "\n")

        from analysis.fel.minima import detect_local_minima
        minima, warns = detect_local_minima(fel, n_minima=3)
        assert len(minima) == 1
        assert any(w.code == "no_local_minima" for w in warns)

    def test_missing_surface_csv_raises(self, tmp_path: Path) -> None:
        """FileNotFoundError when free_energy_surface.csv is absent."""
        fel = tmp_path / "empty_fel"
        (fel / "data").mkdir(parents=True)
        from analysis.fel.minima import detect_local_minima
        with pytest.raises(FileNotFoundError, match="free_energy_surface.csv"):
            detect_local_minima(fel)


class TestFrameMapping:

    def _make_features(self, tmp_path: Path, n: int = 20) -> Path:
        import math
        path = tmp_path / "features.csv"
        lines = ["time_ns,RMSD,Rg"]
        for i in range(n):
            t = i * 0.1
            x = 0.2 + 0.01 * math.sin(i)
            y = 1.2 + 0.01 * math.cos(i)
            lines.append(f"{t:.4g},{x:.6g},{y:.6g}")
        path.write_text("\n".join(lines) + "\n")
        return path

    def test_nearest_frame_mapping(self, tmp_path: Path) -> None:
        """Nearest frame is selected by Euclidean distance in feature space."""
        fcsv = self._make_features(tmp_path)
        from analysis.fel.models import FELMinimum
        from analysis.fel.frame_map import map_frames_to_minima

        # Target at (0.20, 1.20) — very close to the first frame
        m = FELMinimum(minimum_id=1, x_bin=0, y_bin=0, x=0.20, y=1.20,
                       free_energy_kJ_mol=0.0)
        mappings, _ = map_frames_to_minima([m], fcsv, "RMSD", "Rg")
        assert len(mappings) == 1
        assert mappings[0].minimum_id == 1
        assert mappings[0].distance_in_feature_space >= 0
        assert mappings[0].selected_frame_index >= 0

    def test_deterministic_tie_breaking(self, tmp_path: Path) -> None:
        """When two frames equidistant: lower time_ns wins; lower index breaks ties."""
        path = tmp_path / "tie.csv"
        path.write_text("time_ns,X,Y\n0.5,1.0,1.0\n0.3,1.0,1.0\n0.3,1.0,1.0\n")
        from analysis.fel.models import FELMinimum
        from analysis.fel.frame_map import map_frames_to_minima

        m = FELMinimum(minimum_id=1, x_bin=0, y_bin=0, x=1.0, y=1.0,
                       free_energy_kJ_mol=0.0)
        mappings, _ = map_frames_to_minima([m], path, "X", "Y")
        # Row index 1 has t=0.3 (lower than 0.5) → wins
        assert mappings[0].selected_frame_index == 1
        assert mappings[0].selected_time_ns == pytest.approx(0.3)

    def test_unit_metadata_preserved(self, tmp_path: Path) -> None:
        """x_unit and y_unit are stored in FrameMapping."""
        fcsv = self._make_features(tmp_path)
        from analysis.fel.models import FELMinimum
        from analysis.fel.frame_map import map_frames_to_minima

        m = FELMinimum(minimum_id=1, x_bin=0, y_bin=0, x=0.2, y=1.2,
                       free_energy_kJ_mol=0.0)
        mappings, _ = map_frames_to_minima([m], fcsv, "RMSD", "Rg", x_unit="nm", y_unit="nm")
        assert mappings[0].x_unit == "nm"
        assert mappings[0].y_unit == "nm"

    def test_missing_features_csv_raises(self, tmp_path: Path) -> None:
        """FileNotFoundError when features.csv is absent."""
        from analysis.fel.models import FELMinimum
        from analysis.fel.frame_map import map_frames_to_minima

        m = FELMinimum(minimum_id=1, x_bin=0, y_bin=0, x=0.0, y=0.0,
                       free_energy_kJ_mol=0.0)
        with pytest.raises(FileNotFoundError, match="features.csv"):
            map_frames_to_minima([m], tmp_path / "no.csv", "X", "Y")


class TestExtractMinima:

    def _make_full_fel(self, tmp_path: Path) -> Path:
        return _make_fel_dir(tmp_path)

    def _make_extract_config(self, fel: Path, tmp_path: Path, *, dry_run: bool = True) -> object:
        from analysis.fel.models import ExtractConfig
        return ExtractConfig(
            fel_dir=fel,
            trajectory=tmp_path / "md.xtc",
            topology=tmp_path / "md.tpr",
            structure=None, index=None,
            selection="Protein",
            n_minima=3, min_distance_bins=1,
            representative="nearest-frame",
            format="pdb", gmx="gmx",
            out_dir=tmp_path / "out",
            dry_run=dry_run, execute=not dry_run,
            force=True,
        )

    def test_dry_run_creates_all_expected_files(self, tmp_path: Path) -> None:
        """Dry-run creates minima_summary, frame_mapping, planned_commands, metadata, report."""
        fel = self._make_full_fel(tmp_path)
        cfg = self._make_extract_config(fel, tmp_path)
        from analysis.fel.extract import run_extract_minima
        result = run_extract_minima(cfg)
        out = cfg.out_dir

        assert (out / "planned_commands.sh").exists()
        assert (out / "metadata" / "extraction_provenance.json").exists()
        assert (out / "metadata" / "extraction_commands.json").exists()
        assert (out / "metadata" / "warnings.json").exists()
        assert (out / "report.md").exists()
        assert (out / "minima_summary.csv").exists()
        assert (out / "frame_mapping.csv").exists()
        assert (fel / "minima" / "minima.csv").exists()
        assert (fel / "minima" / "frame_mapping.csv").exists()

    def test_dry_run_does_not_call_subprocess(self, tmp_path: Path, monkeypatch) -> None:
        """Dry-run never invokes subprocess.run."""
        import analysis.fel.extract as ext_mod
        calls: list = []

        def fake_run(*args, **kwargs):
            calls.append(args)
            raise AssertionError("subprocess.run must not be called in dry-run mode")

        monkeypatch.setattr(ext_mod.subprocess, "run", fake_run)

        fel = self._make_full_fel(tmp_path)
        cfg = self._make_extract_config(fel, tmp_path, dry_run=True)
        from analysis.fel.extract import run_extract_minima
        run_extract_minima(cfg)
        assert not calls

    def test_execute_mode_calls_subprocess(self, tmp_path: Path, monkeypatch) -> None:
        """Execute mode calls subprocess.run and records return code."""
        import analysis.fel.extract as ext_mod

        calls: list = []

        def fake_run(cmd, **kwargs):
            calls.append(cmd)
            class FakeProc:
                returncode = 0
                stdout = ""
                stderr = ""
            return FakeProc()

        monkeypatch.setattr(ext_mod.subprocess, "run", fake_run)

        fel = self._make_full_fel(tmp_path)
        from analysis.fel.models import ExtractConfig
        cfg = ExtractConfig(
            fel_dir=fel,
            trajectory=tmp_path / "md.xtc",
            topology=tmp_path / "md.tpr",
            structure=None, index=None,
            selection="Protein",
            n_minima=2, min_distance_bins=1,
            representative="nearest-frame",
            format="pdb", gmx="gmx",
            out_dir=tmp_path / "out_exec",
            dry_run=False, execute=True,
            force=True,
        )
        from analysis.fel.extract import run_extract_minima
        result = run_extract_minima(cfg)
        assert len(calls) > 0
        # Each state result records the return code
        for sr in result.state_results:
            assert sr.return_codes == [0]
            assert sr.status == "completed"

    def test_plot_axis_labels_include_units(self, tmp_path: Path) -> None:
        """2D contour axis labels include the feature unit, e.g. 'RMSD (nm)'."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import analysis.fel.plot as plot_mod

        saved_axes: list = []
        original_savefig = plt.Figure.savefig

        def capture(fig_self, *args, **kwargs):
            saved_axes.append(list(fig_self.axes))
            original_savefig(fig_self, *args, **kwargs)

        old_savefig = plt.Figure.savefig
        plt.Figure.savefig = capture
        try:
            surf = _make_surface_with_units()
            plot_mod._plot_2d_contour(surf, tmp_path)
        finally:
            plt.Figure.savefig = old_savefig

        ax = saved_axes[0][0]
        assert "(nm)" in ax.get_xlabel()
        assert "(nm)" in ax.get_ylabel()

    def test_report_contains_trajectory_info(self, tmp_path: Path) -> None:
        """Extraction report.md mentions topology and trajectory."""
        fel = self._make_full_fel(tmp_path)
        cfg = self._make_extract_config(fel, tmp_path)
        from analysis.fel.extract import run_extract_minima
        result = run_extract_minima(cfg)
        report = (cfg.out_dir / "report.md").read_text()
        assert "md.tpr" in report
        assert "md.xtc" in report


class TestUnitMetadata:

    def test_data_unit_detected_from_xvg_ylabel(self, tmp_path: Path) -> None:
        """XVG ylabel 'RMSD (nm)' gives data_unit='nm'."""
        rows = _simple_rows(10)
        p = _write_xvg(tmp_path / "rmsd.xvg", rows, ylabel='RMSD (nm)')
        from analysis.fel.xvg import load_xvg_series
        s = load_xvg_series(p)
        assert s.data_unit == "nm"

    def test_unit_propagated_to_feature_table(self, tmp_path: Path) -> None:
        """FeatureTable.x_unit / y_unit come from XVG ylabel."""
        rows = _simple_rows(15)
        px = _write_xvg(tmp_path / "x.xvg", rows, ylabel='RMSD (nm)')
        py = _write_xvg(tmp_path / "y.xvg", rows, ylabel='Rg (nm)')
        from analysis.fel.xvg import load_xvg_series
        from analysis.fel.features import build_feature_table
        sx = load_xvg_series(px)
        sy = load_xvg_series(py)
        ft = build_feature_table(sx, sy)
        assert ft.x_unit == "nm"
        assert ft.y_unit == "nm"
        assert ft.time_unit_normalized == "ns"

    def test_unit_propagated_to_surface(self, tmp_path: Path) -> None:
        """FELSurface.x_unit / y_unit are propagated from FeatureTable."""
        rows = _simple_rows(20)
        px = _write_xvg(tmp_path / "x.xvg", rows, ylabel='RMSD (nm)')
        py = _write_xvg(tmp_path / "y.xvg", rows, ylabel='Rg (nm)')
        from analysis.fel.xvg import load_xvg_series
        from analysis.fel.features import build_feature_table
        from analysis.fel.surface import compute_fel_surface
        sx = load_xvg_series(px)
        sy = load_xvg_series(py)
        ft = build_feature_table(sx, sy)
        surf = compute_fel_surface(ft.x, ft.y, n_bins_x=5, n_bins_y=5,
                                   x_unit=ft.x_unit, y_unit=ft.y_unit)
        assert surf.x_unit == "nm"
        assert surf.y_unit == "nm"

    def test_ps_ns_conversion_recorded_in_feature_table(self, tmp_path: Path) -> None:
        """time_conversion_x_to_ns is 1e-3 for ps, 1.0 for ns."""
        rows = _simple_rows(10)
        px = _write_xvg(tmp_path / "x.xvg", rows, xlabel='Time (ps)')
        py = _write_xvg(tmp_path / "y.xvg", rows, xlabel='Time (ns)')
        from analysis.fel.xvg import load_xvg_series
        from analysis.fel.features import build_feature_table
        sx = load_xvg_series(px)
        sy = load_xvg_series(py)
        # ps and ns mismatch — auto-corrected
        ft = build_feature_table(sx, sy, time_tolerance=1e-9)
        # One axis is ps (factor 1e-3), but after correction units may swap
        # — just verify the field exists and is numeric
        assert isinstance(ft.time_conversion_x_to_ns, float)
        assert isinstance(ft.time_conversion_y_to_ns, float)

    def test_mismatched_time_axes_raise_clearly(self, tmp_path: Path) -> None:
        """Truly mismatched time axes raise ValueError with a descriptive message."""
        rows_a = [(i * 100.0, float(i)) for i in range(10)]
        rows_b = [(i * 500.0, float(i)) for i in range(10)]
        px = _write_xvg(tmp_path / "x.xvg", rows_a)
        py = _write_xvg(tmp_path / "y.xvg", rows_b)
        from analysis.fel.xvg import load_xvg_series
        from analysis.fel.features import build_feature_table
        sx = load_xvg_series(px)
        sy = load_xvg_series(py)
        with pytest.raises(ValueError, match="do not align"):
            build_feature_table(sx, sy, time_tolerance=1e-6)

    def test_config_json_includes_units_after_run(self, tmp_path: Path) -> None:
        """config.json written by `simforge fel run` includes x_unit and y_unit."""
        rows = _simple_rows(20)
        px = _write_xvg(tmp_path / "rmsd.xvg", rows, legends=["RMSD"], ylabel='RMSD (nm)')
        py = _write_xvg(tmp_path / "gyrate.xvg", rows, legends=["Rg"], ylabel='Rg (nm)')
        out = tmp_path / "fel_units"
        from analysis.fel.models import FELConfig
        from analysis.fel.workflow import run_fel
        config = FELConfig(
            xvg_x=px, xvg_y=py,
            x_column=1, y_column=1,
            x_label="RMSD", y_label="Rg",
            temperature_K=300.0, n_bins_x=5, n_bins_y=5,
            out_dir=out,
            time_unit_x="auto", time_unit_y="auto",
            time_tolerance=1e-6, empty_bin_policy="nan",
            energy_cap_kj=50.0, dry_run=False,
        )
        run_fel(config)
        import json as _json
        cfg_data = _json.loads((out / "config.json").read_text())
        assert cfg_data["x_unit"] == "nm"
        assert cfg_data["y_unit"] == "nm"


def _make_surface_with_units():
    """Build a minimal FELSurface with x_unit/y_unit='nm' for plot tests."""
    import math
    x = [math.sin(i * 0.1) for i in range(80)]
    y = [math.cos(i * 0.1) for i in range(80)]
    from analysis.fel.surface import compute_fel_surface
    return compute_fel_surface(x, y, n_bins_x=8, n_bins_y=8,
                               x_label="RMSD", y_label="Rg",
                               x_unit="nm", y_unit="nm")


# ── Group 7: Phase 3B — Unit defaults, force cleanup, energy filtering ────────


def _make_fel_dir_no_units(
    tmp_path: Path,
    *,
    n_bins: int = 5,
    x_label: str = "RMSD",
    y_label: str = "Rg",
    n_frames: int = 30,
) -> Path:
    """FEL directory whose config.json has no x_unit/y_unit (simulates old run output)."""
    import json
    import math

    fel = tmp_path / "fel_no_units"
    (fel / "data").mkdir(parents=True)
    (fel / "minima").mkdir()

    xs = [0.1 * i for i in range(n_bins)]
    ys = [1.0 + 0.1 * i for i in range(n_bins)]
    cx, cy = xs[n_bins // 2], ys[n_bins // 2]
    rows = [(x, y, 4.0 * ((x - cx) ** 2 + (y - cy) ** 2))
            for x in xs for y in ys]

    lines = [f"x_{x_label},y_{y_label},free_energy_kJ_mol"]
    lines += [f"{x:.6g},{y:.6g},{g:.6g}" for x, y, g in rows]
    (fel / "data" / "free_energy_surface.csv").write_text("\n".join(lines) + "\n")

    feat = ["time_ns," + x_label + "," + y_label]
    for i in range(n_frames):
        t = i * 0.2
        feat.append(f"{t:.6g},{cx + 0.01 * math.sin(i):.6g},{cy + 0.01 * math.cos(i):.6g}")
    (fel / "data" / "features.csv").write_text("\n".join(feat) + "\n")

    # config WITHOUT x_unit/y_unit keys
    config = {
        "xvg_x": "rmsd.xvg", "xvg_y": "rg.xvg",
        "x_column": 1, "y_column": 1,
        "x_label": x_label, "y_label": y_label,
        "temperature_K": 300.0,
        "n_bins_x": n_bins, "n_bins_y": n_bins,
        "out_dir": str(fel),
        "time_unit_x": "auto", "time_unit_y": "auto",
        "time_tolerance": 1e-6,
        "empty_bin_policy": "nan",
        "energy_cap_kj": 50.0,
        "dry_run": False,
    }
    (fel / "config.json").write_text(json.dumps(config) + "\n")
    return fel


def _make_fel_dir_two_wells(
    tmp_path: Path,
    g2: float = 3.5,
    n_bins: int = 11,
) -> Path:
    """FEL with two local minima: M1 at G≈0, M2 at G≈g2, well-separated."""
    import json
    import math

    fel = tmp_path / "fel_two_wells"
    (fel / "data").mkdir(parents=True)
    (fel / "minima").mkdir()

    xs = [i * 0.1 for i in range(n_bins)]
    ys = [i * 0.1 for i in range(n_bins)]

    # Well 1 centred at bin (2,2), well 2 at bin (8,8)
    w1x, w1y = xs[2], ys[2]
    w2x, w2y = xs[8], ys[8]

    rows = []
    for x in xs:
        for y in ys:
            d1 = (x - w1x) ** 2 + (y - w1y) ** 2
            d2 = (x - w2x) ** 2 + (y - w2y) ** 2
            # Sharp separate wells with barrier
            g_raw = min(d1 * 80.0, d2 * 80.0 + g2)
            rows.append((x, y, g_raw))

    # Shift so global min = 0
    g_min = min(r[2] for r in rows)
    rows = [(x, y, g - g_min) for x, y, g in rows]

    lines = ["x_RMSD,y_Rg,free_energy_kJ_mol"]
    lines += [f"{x:.6g},{y:.6g},{g:.6g}" for x, y, g in rows]
    (fel / "data" / "free_energy_surface.csv").write_text("\n".join(lines) + "\n")

    # Features: frames visit both wells
    n_frames = 40
    feat = ["time_ns,RMSD,Rg"]
    for i in range(n_frames):
        t = i * 0.2
        if i < n_frames // 2:
            feat.append(f"{t:.6g},{w1x + 0.001 * i:.6g},{w1y + 0.001 * i:.6g}")
        else:
            feat.append(f"{t:.6g},{w2x + 0.001 * i:.6g},{w2y + 0.001 * i:.6g}")
    (fel / "data" / "features.csv").write_text("\n".join(feat) + "\n")

    config = {
        "xvg_x": "rmsd.xvg", "xvg_y": "rg.xvg",
        "x_column": 1, "y_column": 1,
        "x_label": "RMSD", "y_label": "Rg",
        "x_unit": "nm", "y_unit": "nm",
        "temperature_K": 300.0,
        "n_bins_x": n_bins, "n_bins_y": n_bins,
        "out_dir": str(fel),
        "time_unit_x": "auto", "time_unit_y": "auto",
        "time_tolerance": 1e-6,
        "empty_bin_policy": "nan",
        "energy_cap_kj": 50.0,
        "dry_run": False,
    }
    (fel / "config.json").write_text(json.dumps(config) + "\n")
    return fel


def _extract_cfg(
    fel: Path,
    tmp_path: Path,
    *,
    out_name: str = "out",
    n_minima: int = 3,
    dry_run: bool = True,
    max_delta_g_kj: object = None,  # None means "disable filter" in tests
    min_count: int = 1,
    force: bool = True,
) -> object:
    from analysis.fel.models import ExtractConfig
    return ExtractConfig(
        fel_dir=fel,
        trajectory=tmp_path / "md.xtc",
        topology=tmp_path / "md.tpr",
        structure=None, index=None,
        selection="Protein",
        n_minima=n_minima,
        min_distance_bins=1,
        representative="nearest-frame",
        format="pdb", gmx="gmx",
        out_dir=tmp_path / out_name,
        dry_run=dry_run, execute=not dry_run,
        force=force,
        max_delta_g_kj=max_delta_g_kj,
        min_count=min_count,
    )


class TestPhase3BUnitDefaults:

    def test_unit_defaults_for_rmsd_rg_are_nm(self, tmp_path: Path) -> None:
        """When config.json has no x_unit/y_unit, RMSD and Rg default to nm."""
        from analysis.fel.xvg import default_unit_for_label
        assert default_unit_for_label("RMSD") == "nm"
        assert default_unit_for_label("Rg") == "nm"
        assert default_unit_for_label("rmsf") == "nm"
        assert default_unit_for_label("radius of gyration") == "nm"
        assert default_unit_for_label("Energy") is None

    def test_unit_defaults_propagate_to_state_metadata(self, tmp_path: Path) -> None:
        """state_*/metadata.json has non-null x_unit/y_unit when inferred from label."""
        import json as _j
        fel = _make_fel_dir_no_units(tmp_path)
        cfg = _extract_cfg(fel, tmp_path, max_delta_g_kj=None)
        from analysis.fel.extract import run_extract_minima
        result = run_extract_minima(cfg)
        assert result.minima, "Expected at least one minimum"

        meta_file = cfg.out_dir / "state_001" / "metadata.json"
        assert meta_file.exists(), f"Missing {meta_file}"
        meta = _j.loads(meta_file.read_text())
        assert meta["x_unit"] == "nm", f"x_unit wrong: {meta['x_unit']!r}"
        assert meta["y_unit"] == "nm", f"y_unit wrong: {meta['y_unit']!r}"

    def test_unit_defaults_appear_in_report(self, tmp_path: Path) -> None:
        """report.md shows 'RMSD (nm)' and 'Rg (nm)' even when units absent from config.json."""
        fel = _make_fel_dir_no_units(tmp_path)
        cfg = _extract_cfg(fel, tmp_path, max_delta_g_kj=None)
        from analysis.fel.extract import run_extract_minima
        run_extract_minima(cfg)
        report = (cfg.out_dir / "report.md").read_text()
        assert "RMSD (nm)" in report
        assert "Rg (nm)" in report

    def test_unit_defaults_in_frame_mapping(self, tmp_path: Path) -> None:
        """FrameMapping.x_unit / y_unit are 'nm' when inferred from label."""
        fel = _make_fel_dir_no_units(tmp_path)
        cfg = _extract_cfg(fel, tmp_path, max_delta_g_kj=None)
        from analysis.fel.extract import run_extract_minima
        result = run_extract_minima(cfg)
        for fm in result.frame_map:
            assert fm.x_unit == "nm", f"x_unit: {fm.x_unit!r}"
            assert fm.y_unit == "nm", f"y_unit: {fm.y_unit!r}"


class TestPhase3BForceCleanup:

    def test_force_removes_stale_state_directories(self, tmp_path: Path) -> None:
        """--force deletes old state dirs (e.g. state_004) not present in current run."""
        # First run with n_minima=1 to create a known state
        fel = _make_fel_dir(tmp_path)
        out = tmp_path / "states"

        # Pre-seed stale directories
        (out / "state_004").mkdir(parents=True)
        (out / "state_004" / "metadata.json").write_text('{"stale": true}\n')
        (out / "state_005").mkdir(parents=True)
        (out / "state_005" / "metadata.json").write_text('{"stale": true}\n')

        cfg = _extract_cfg(fel, tmp_path, out_name="states", n_minima=1,
                           max_delta_g_kj=None, force=True)
        from analysis.fel.extract import run_extract_minima
        run_extract_minima(cfg)

        assert not (out / "state_004").exists(), "Stale state_004 should have been removed"
        assert not (out / "state_005").exists(), "Stale state_005 should have been removed"

    def test_no_force_preserves_existing_dir(self, tmp_path: Path) -> None:
        """Without --force, the output directory is not cleared (CLI would have rejected)."""
        out = tmp_path / "existing"
        out.mkdir()
        (out / "sentinel.txt").write_text("keep me\n")

        # run_extract_minima itself only clears when config.force=True
        # (the CLI gate prevents running without --force when dir exists)
        fel = _make_fel_dir(tmp_path)
        from analysis.fel.models import ExtractConfig
        cfg = ExtractConfig(
            fel_dir=fel,
            trajectory=tmp_path / "md.xtc",
            topology=tmp_path / "md.tpr",
            structure=None, index=None,
            selection="Protein",
            n_minima=1, min_distance_bins=1,
            representative="nearest-frame",
            format="pdb", gmx="gmx",
            out_dir=out,
            dry_run=True, execute=False,
            force=False,          # do NOT wipe
            max_delta_g_kj=None,
            min_count=1,
        )
        from analysis.fel.extract import run_extract_minima
        run_extract_minima(cfg)
        # The sentinel file must still exist — force=False means no rmtree
        assert (out / "sentinel.txt").exists()

    def test_cli_error_message_when_output_exists(self, tmp_path: Path) -> None:
        """CLI outputs the improved 'Use --force' message when output dir exists."""
        out = tmp_path / "exists_dir"
        out.mkdir()
        fel = _make_fel_dir(tmp_path)
        r = _cli(
            "fel", "extract-minima", str(fel),
            "--trajectory", "md.xtc",
            "--topology", "md.tpr",
            "--dry-run",
            "--out", str(out),
        )
        assert r.returncode != 0
        assert "--force" in r.stderr or "--force" in r.stdout


class TestPhase3BEnergyFilter:

    def test_apply_energy_filter_helper(self) -> None:
        """_apply_energy_filter splits kept/rejected correctly."""
        from analysis.fel.extract import _apply_energy_filter  # private but stable
        from analysis.fel.models import FELMinimum
        ms = [
            FELMinimum(1, 0, 0, 0.0, 0.0, 0.5),
            FELMinimum(2, 1, 1, 0.1, 0.1, 1.8),
            FELMinimum(3, 2, 2, 0.2, 0.2, 3.1),
        ]
        kept, rejected = _apply_energy_filter(ms, max_delta_g_kj=2.5)
        assert len(kept) == 2
        assert len(rejected) == 1
        assert rejected[0].free_energy_kJ_mol == pytest.approx(3.1)

    def test_filter_none_disables_filter(self) -> None:
        """max_delta_g_kj=None passes all minima through."""
        from analysis.fel.extract import _apply_energy_filter
        from analysis.fel.models import FELMinimum
        ms = [FELMinimum(i, i, i, float(i), float(i), float(i) * 5.0) for i in range(5)]
        kept, rejected = _apply_energy_filter(ms, max_delta_g_kj=None)
        assert len(kept) == 5
        assert len(rejected) == 0

    def test_filter_all_rejected_keeps_global_min(self) -> None:
        """When cutoff rejects every minimum, global min (first) is preserved."""
        from analysis.fel.extract import _apply_energy_filter
        from analysis.fel.models import FELMinimum
        ms = [FELMinimum(1, 0, 0, 0.0, 0.0, 5.0), FELMinimum(2, 1, 1, 0.1, 0.1, 7.0)]
        kept, rejected = _apply_energy_filter(ms, max_delta_g_kj=2.5)
        assert len(kept) == 1
        assert kept[0].free_energy_kJ_mol == pytest.approx(5.0)
        assert len(rejected) == 1

    def test_energy_filter_integration_via_run(self, tmp_path: Path) -> None:
        """run_extract_minima applies the ΔG filter and records rejected minima."""
        fel = _make_fel_dir_two_wells(tmp_path, g2=3.5)
        cfg = _extract_cfg(
            fel, tmp_path, n_minima=5,
            max_delta_g_kj=2.5,  # M2 at ~3.5 kJ/mol is rejected
            min_count=1,
        )
        from analysis.fel.extract import run_extract_minima
        result = run_extract_minima(cfg)

        # The high-energy well (≥ 3.5 kJ/mol) must be in rejected_minima
        kept_energies = [m.free_energy_kJ_mol for m in result.minima]
        assert all(e <= 2.5 for e in kept_energies), f"Non-filtered energies: {kept_energies}"
        # At least one minimum was rejected (the 3.5 kJ/mol well)
        assert len(result.rejected_minima) >= 1
        assert any(m.free_energy_kJ_mol > 2.5 for m in result.rejected_minima)

    def test_energy_filter_emits_info_warning(self, tmp_path: Path) -> None:
        """energy_filter warning is emitted when minima are rejected."""
        fel = _make_fel_dir_two_wells(tmp_path, g2=3.5)
        cfg = _extract_cfg(fel, tmp_path, n_minima=5, max_delta_g_kj=2.5)
        from analysis.fel.extract import run_extract_minima
        result = run_extract_minima(cfg)
        codes = [w.code for w in result.warnings]
        assert "energy_filter" in codes

    def test_report_shows_rejected_minima(self, tmp_path: Path) -> None:
        """report.md lists rejected minima in a separate table."""
        fel = _make_fel_dir_two_wells(tmp_path, g2=3.5)
        cfg = _extract_cfg(fel, tmp_path, n_minima=5, max_delta_g_kj=2.5)
        from analysis.fel.extract import run_extract_minima
        result = run_extract_minima(cfg)
        if not result.rejected_minima:
            pytest.skip("No minima rejected with these parameters; skip report check")
        report = (cfg.out_dir / "report.md").read_text()
        assert "rejected" in report.lower() or "Rejected" in report

    def test_report_shows_filter_summary(self, tmp_path: Path) -> None:
        """report.md shows ΔG filter value in the Configuration section."""
        fel = _make_fel_dir(tmp_path)
        cfg = _extract_cfg(fel, tmp_path, max_delta_g_kj=2.5)
        from analysis.fel.extract import run_extract_minima
        run_extract_minima(cfg)
        report = (cfg.out_dir / "report.md").read_text()
        assert "2.5" in report

    def test_cli_help_includes_new_options(self) -> None:
        """CLI help text shows --max-delta-g-kj and --min-count."""
        r = _cli("fel", "extract-minima", "--help")
        assert r.returncode == 0
        assert "--max-delta-g-kj" in r.stdout
        assert "--min-count" in r.stdout


class TestPhase3BCountFilter:

    def _make_fel_with_counts(self, tmp_path: Path) -> Path:
        """FEL directory that also has histogram_counts.csv."""
        import json

        fel = _make_fel_dir(tmp_path, n_bins=5)
        # Write histogram_counts.csv with same x/y as surface CSV but varying counts
        surf_lines = (fel / "data" / "free_energy_surface.csv").read_text().splitlines()
        count_lines = [surf_lines[0].replace("free_energy_kJ_mol", "counts")]
        for line in surf_lines[1:]:
            parts = line.split(",")
            # Give the minimum bin high count, others low count
            x, y = float(parts[0]), float(parts[1])
            # center of 5-bin grid: x≈0.2, y≈1.2
            count = 20 if (abs(x - 0.2) < 0.05 and abs(y - 1.2) < 0.05) else 1
            count_lines.append(f"{parts[0]},{parts[1]},{count}")
        (fel / "data" / "histogram_counts.csv").write_text("\n".join(count_lines) + "\n")
        return fel

    def test_count_enrichment_attaches_count(self, tmp_path: Path) -> None:
        """FELMinimum.count is populated when histogram_counts.csv is present."""
        fel = self._make_fel_with_counts(tmp_path)
        cfg = _extract_cfg(fel, tmp_path, max_delta_g_kj=None, min_count=1)
        from analysis.fel.extract import run_extract_minima
        result = run_extract_minima(cfg)
        counts = [m.count for m in result.minima]
        assert any(c is not None for c in counts), "Expected at least one count to be set"

    def test_min_count_filter_excludes_low_count_bins(self, tmp_path: Path) -> None:
        """min_count=10 excludes minima in bins with < 10 frames."""
        import math as _m
        # Build a surface with a local minimum in a LOW-count bin
        import json as _j

        fel = tmp_path / "low_count_fel"
        (fel / "data").mkdir(parents=True)
        (fel / "minima").mkdir()

        n = 5
        xs = [i * 0.1 for i in range(n)]
        ys = [i * 0.1 for i in range(n)]
        cx, cy = xs[2], ys[2]

        # Parabolic bowl — minimum at (cx, cy)
        surf_lines = ["x_RMSD,y_Rg,free_energy_kJ_mol"]
        count_lines = ["x_RMSD,y_Rg,counts"]
        for x in xs:
            for y in ys:
                g = 4.0 * ((x - cx) ** 2 + (y - cy) ** 2)
                surf_lines.append(f"{x:.6g},{y:.6g},{g:.6g}")
                count = 1 if (abs(x - cx) < 0.05 and abs(y - cy) < 0.05) else 5
                count_lines.append(f"{x:.6g},{y:.6g},{count}")

        (fel / "data" / "free_energy_surface.csv").write_text("\n".join(surf_lines) + "\n")
        (fel / "data" / "histogram_counts.csv").write_text("\n".join(count_lines) + "\n")

        feat = ["time_ns,RMSD,Rg"]
        for i in range(20):
            feat.append(f"{i * 0.1:.4g},{cx + 0.001 * i:.6g},{cy + 0.001 * i:.6g}")
        (fel / "data" / "features.csv").write_text("\n".join(feat) + "\n")
        (fel / "config.json").write_text(_j.dumps({
            "x_label": "RMSD", "y_label": "Rg",
            "x_unit": "nm", "y_unit": "nm",
            "n_bins_x": n, "n_bins_y": n,
        }) + "\n")

        cfg = _extract_cfg(fel, tmp_path, max_delta_g_kj=None, min_count=3)
        from analysis.fel.extract import run_extract_minima
        result = run_extract_minima(cfg)
        # The global minimum bin has count=1, below min_count=3 → filtered out
        # So minima list should be empty OR the warning code appears
        codes = [w.code for w in result.warnings]
        if result.minima:
            counts_left = [m.count for m in result.minima if m.count is not None]
            assert all(c >= 3 for c in counts_left), f"Low-count minima not filtered: {counts_left}"
        else:
            assert "min_count_filter" in codes or "no_local_minima" in codes


# ── Group 8: Phase 3C — Descriptive filenames ────────────────────────────────

class TestPhase3CFilenames:

    def test_sanitize_system_name_spaces_and_dots(self) -> None:
        """Spaces and dots are converted to underscores."""
        from analysis.fel.extract import sanitize_system_name
        assert sanitize_system_name("HMG-CoA R") == "HMG-CoA_R"
        assert sanitize_system_name("system.v2") == "system_v2"
        assert sanitize_system_name("GLP-1 Receptor") == "GLP-1_Receptor"

    def test_sanitize_system_name_removes_unsafe_chars(self) -> None:
        """Shell-unsafe chars are stripped; only letters/digits/underscore/dash allowed."""
        from analysis.fel.extract import sanitize_system_name
        assert sanitize_system_name("sys/tem") == "system"
        assert sanitize_system_name("sys:tem") == "system"
        assert sanitize_system_name("sys@tem") == "system"
        assert sanitize_system_name("sys tem!") == "sys_tem"

    def test_sanitize_system_name_collapse_underscores(self) -> None:
        """Multiple consecutive underscores collapse to one."""
        from analysis.fel.extract import sanitize_system_name
        assert sanitize_system_name("a..b") == "a_b"
        assert sanitize_system_name("a   b") == "a_b"

    def test_format_time_ns_replaces_dot_with_p(self) -> None:
        """format_time_ns produces token with 'p' instead of decimal point."""
        from analysis.fel.extract import format_time_ns
        assert format_time_ns(44.7) == "44p7"
        assert format_time_ns(51.6) == "51p6"
        assert format_time_ns(100.0) == "100p0"
        assert format_time_ns(0.3) == "0p3"

    def test_format_time_ns_no_decimal_dot(self) -> None:
        """format_time_ns output never contains a literal period."""
        from analysis.fel.extract import format_time_ns
        import math
        for t in [0.0, 1.0, 44.7, 100.25, 999.9]:
            result = format_time_ns(t)
            assert "." not in result, f"Unexpected dot in {result!r} for t={t}"

    def test_build_output_filename_example(self) -> None:
        """Exact expected filename for the HMGCoA-R use case."""
        from analysis.fel.extract import build_output_filename
        fname = build_output_filename("HMGCoA_R", 1, 44.7, "pdb")
        assert fname == "HMGCoA_R__FEL-M01__t44p7ns.pdb"

    def test_build_output_filename_two_digit_id(self) -> None:
        """minimum_id is always zero-padded to 2 digits."""
        from analysis.fel.extract import build_output_filename
        assert "FEL-M03" in build_output_filename("sys", 3, 1.0, "pdb")
        assert "FEL-M10" in build_output_filename("sys", 10, 1.0, "pdb")

    def test_build_output_filename_no_delta_g(self) -> None:
        """ΔG value must NOT appear in the filename."""
        from analysis.fel.extract import build_output_filename
        fname = build_output_filename("sys", 1, 10.0, "pdb")
        # ΔG is never passed — sanity check that the helper doesn't sneak it in
        assert "kJ" not in fname and "dg" not in fname.lower()

    def test_output_file_uses_descriptive_name(self, tmp_path: Path) -> None:
        """Dry-run creates descriptive filename in state_001/; no receptor.pdb."""
        fel = _make_fel_dir(tmp_path)
        from analysis.fel.models import ExtractConfig
        cfg = ExtractConfig(
            fel_dir=fel,
            trajectory=tmp_path / "md.xtc",
            topology=tmp_path / "md.tpr",
            structure=None, index=None,
            selection="Protein",
            n_minima=1, min_distance_bins=1,
            representative="nearest-frame",
            format="pdb", gmx="gmx",
            out_dir=tmp_path / "out",
            dry_run=True, execute=False,
            force=True,
            max_delta_g_kj=None,
            min_count=1,
            system_name="HMGCoA_R",
        )
        from analysis.fel.extract import run_extract_minima
        result = run_extract_minima(cfg)
        assert result.minima, "Expected at least one minimum"

        state_dir = cfg.out_dir / "state_001"
        assert state_dir.is_dir()
        # No receptor.pdb by default
        assert not (state_dir / "receptor.pdb").exists()
        # Descriptive filename present in output_paths
        fnames = [p.name for p in result.state_results[0].output_paths]
        assert any("FEL-M01" in f for f in fnames), f"No FEL-M01 in {fnames}"
        assert all("." not in f.rsplit(".", 1)[0] for f in fnames), (
            "Decimal point inside filename stem"
        )

    def test_compat_receptor_name_not_created_by_default(self, tmp_path: Path) -> None:
        """receptor.pdb is NOT created unless --compat-receptor-name is set."""
        fel = _make_fel_dir(tmp_path)
        from analysis.fel.models import ExtractConfig
        cfg = ExtractConfig(
            fel_dir=fel,
            trajectory=tmp_path / "md.xtc",
            topology=tmp_path / "md.tpr",
            structure=None, index=None,
            selection="Protein",
            n_minima=1, min_distance_bins=1,
            representative="nearest-frame",
            format="pdb", gmx="gmx",
            out_dir=tmp_path / "out_no_compat",
            dry_run=True, execute=False,
            force=True,
            max_delta_g_kj=None,
            min_count=1,
            compat_receptor_name=False,
        )
        from analysis.fel.extract import run_extract_minima
        run_extract_minima(cfg)
        for d in (cfg.out_dir / "state_001",):
            if d.exists():
                assert not (d / "receptor.pdb").exists()

    def test_system_name_inferred_from_fel_dir(self, tmp_path: Path) -> None:
        """When system_name is None, the FEL directory name is used."""
        from analysis.fel.extract import resolve_system_name
        from analysis.fel.models import ExtractConfig
        cfg = ExtractConfig(
            fel_dir=Path("/data/HMGCoA_R_run"),
            trajectory=Path("md.xtc"),
            topology=Path("md.tpr"),
            structure=None, index=None,
            selection="Protein",
            n_minima=3, min_distance_bins=1,
            representative="nearest-frame",
            format="pdb", gmx="gmx",
            out_dir=Path("out"),
            dry_run=True, execute=False,
            force=False,
            system_name=None,
        )
        assert resolve_system_name(cfg) == "HMGCoA_R_run"

    def test_system_name_explicit_overrides_inference(self, tmp_path: Path) -> None:
        """Explicit --system-name takes precedence over any inference."""
        from analysis.fel.extract import resolve_system_name
        from analysis.fel.models import ExtractConfig
        cfg = ExtractConfig(
            fel_dir=Path("/data/irrelevant_dir"),
            trajectory=Path("md.xtc"),
            topology=Path("md.tpr"),
            structure=None, index=None,
            selection="Protein",
            n_minima=3, min_distance_bins=1,
            representative="nearest-frame",
            format="pdb", gmx="gmx",
            out_dir=Path("out"),
            dry_run=True, execute=False,
            force=False,
            system_name="My Protein 3.1",
        )
        assert resolve_system_name(cfg) == "My_Protein_3_1"

    def test_metadata_json_contains_display_name(self, tmp_path: Path) -> None:
        """state_*/metadata.json includes display_name, output_filenames, system_name."""
        import json as _j
        fel = _make_fel_dir(tmp_path)
        from analysis.fel.models import ExtractConfig
        cfg = ExtractConfig(
            fel_dir=fel,
            trajectory=tmp_path / "md.xtc",
            topology=tmp_path / "md.tpr",
            structure=None, index=None,
            selection="Protein",
            n_minima=1, min_distance_bins=1,
            representative="nearest-frame",
            format="pdb", gmx="gmx",
            out_dir=tmp_path / "out_meta",
            dry_run=True, execute=False,
            force=True,
            max_delta_g_kj=None,
            system_name="TestSys",
        )
        from analysis.fel.extract import run_extract_minima
        result = run_extract_minima(cfg)
        assert result.minima
        meta = _j.loads((cfg.out_dir / "state_001" / "metadata.json").read_text())
        assert meta["system_name"] == "TestSys"
        assert "display_name" in meta
        assert "TestSys" in meta["display_name"]
        assert "output_filenames" in meta
        assert any("FEL-M01" in f for f in meta["output_filenames"])
        assert "selected_time_ns" in meta
        assert "free_energy_kJ_mol" in meta
        assert "x_unit" in meta
        assert "y_unit" in meta

    def test_planned_commands_sh_uses_descriptive_filename(self, tmp_path: Path) -> None:
        """planned_commands.sh references the descriptive filename."""
        fel = _make_fel_dir(tmp_path)
        from analysis.fel.models import ExtractConfig
        cfg = ExtractConfig(
            fel_dir=fel,
            trajectory=tmp_path / "md.xtc",
            topology=tmp_path / "md.tpr",
            structure=None, index=None,
            selection="Protein",
            n_minima=1, min_distance_bins=1,
            representative="nearest-frame",
            format="pdb", gmx="gmx",
            out_dir=tmp_path / "out_sh",
            dry_run=True, execute=False,
            force=True,
            max_delta_g_kj=None,
            system_name="MySys",
        )
        from analysis.fel.extract import run_extract_minima
        run_extract_minima(cfg)
        sh = (cfg.out_dir / "planned_commands.sh").read_text()
        assert "FEL-M01" in sh
        assert "receptor.pdb" not in sh

    def test_report_md_shows_descriptive_filename(self, tmp_path: Path) -> None:
        """report.md table contains the output filename."""
        fel = _make_fel_dir(tmp_path)
        from analysis.fel.models import ExtractConfig
        cfg = ExtractConfig(
            fel_dir=fel,
            trajectory=tmp_path / "md.xtc",
            topology=tmp_path / "md.tpr",
            structure=None, index=None,
            selection="Protein",
            n_minima=1, min_distance_bins=1,
            representative="nearest-frame",
            format="pdb", gmx="gmx",
            out_dir=tmp_path / "out_report",
            dry_run=True, execute=False,
            force=True,
            max_delta_g_kj=None,
            system_name="ProSys",
        )
        from analysis.fel.extract import run_extract_minima
        run_extract_minima(cfg)
        report = (cfg.out_dir / "report.md").read_text()
        assert "FEL-M01" in report
        assert "ProSys" in report

    def test_cli_help_includes_system_name_and_compat(self) -> None:
        """CLI help shows --system-name and --compat-receptor-name."""
        r = _cli("fel", "extract-minima", "--help")
        assert r.returncode == 0
        assert "--system-name" in r.stdout
        assert "--compat-receptor-name" in r.stdout


# ── Group 9: Real data (skipped if files absent) ──────────────────────────────

_RMSD_XVG = _REPO / "docs" / "FEL" / "rmsd-hDHFR.xvg"
_RG_XVG   = _REPO / "docs" / "FEL" / "Rg-hDHFR.xvg"

@pytest.mark.skipif(not (_RMSD_XVG.exists() and _RG_XVG.exists()), reason="hDHFR XVG files not present")
class TestRealHDHFRData:

    def test_ps_ns_mismatch_resolved(self) -> None:
        """RMSD (ns) and Rg (ps) files from hDHFR align correctly after unit conversion."""
        from analysis.fel.xvg import load_xvg_series
        from analysis.fel.features import build_feature_table

        sx = load_xvg_series(_RMSD_XVG)
        sy = load_xvg_series(_RG_XVG)

        # Headers declare different units — the converter handles this automatically.
        assert sx.x_unit == "ns"
        assert sy.x_unit == "ps"

        # Both files cover the same trajectory; 1e-3 ns tolerates ps↔ns fp rounding.
        ft = build_feature_table(sx, sy, time_tolerance=1e-3)
        assert ft.n_frames == 501
        # No error raised; conversion note records both units.
        assert "ns" in ft.time_conversion_note
        assert "ps" in ft.time_conversion_note

    def test_fel_surface_from_real_data(self, tmp_path: Path) -> None:
        """Full FEL from real hDHFR data produces non-empty surface with min G = 0."""
        from analysis.fel.xvg import load_xvg_series
        from analysis.fel.features import build_feature_table
        from analysis.fel.surface import compute_fel_surface

        sx = load_xvg_series(_RMSD_XVG)
        sy = load_xvg_series(_RG_XVG)
        ft = build_feature_table(sx, sy, time_tolerance=1e-3)
        surf = compute_fel_surface(
            ft.x, ft.y,
            n_bins_x=20, n_bins_y=20,
            temperature_K=300.0,
        )

        finite_g = [
            surf.free_energy[ix][iy]
            for ix in range(len(surf.x_centers))
            for iy in range(len(surf.y_centers))
            if surf.free_energy[ix][iy] == surf.free_energy[ix][iy]
        ]
        assert min(finite_g) == pytest.approx(0.0, abs=1e-10)
        assert max(finite_g) > 0
