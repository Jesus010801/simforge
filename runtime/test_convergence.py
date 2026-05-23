"""Tests for runtime/convergence_analyzer.py — no GROMACS required."""
from __future__ import annotations

import math
from pathlib import Path

import pytest

from runtime.xvg_parser import XVGData, XVGSeries
from runtime.convergence_analyzer import (
    analyze_rmsd_convergence,
    analyze_energy_stability,
    ConvergenceResult,
    StabilityResult,
    _mean,
    _std,
    _linear_slope,
)


# ─────────────────────────────────────────────────────────────────────────────
# Internal helper tests
# ─────────────────────────────────────────────────────────────────────────────

class TestHelpers:
    def test_mean_empty(self):
        assert _mean([]) == 0.0

    def test_mean_values(self):
        assert _mean([1.0, 2.0, 3.0]) == pytest.approx(2.0)

    def test_std_empty(self):
        assert _std([]) == 0.0

    def test_std_single(self):
        assert _std([5.0]) == 0.0

    def test_std_known(self):
        vals = [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]
        # population std = 2.0; sample std slightly larger
        assert _std(vals) == pytest.approx(2.138, abs=0.01)

    def test_linear_slope_flat(self):
        xs = [0.0, 1.0, 2.0, 3.0]
        ys = [5.0, 5.0, 5.0, 5.0]
        assert _linear_slope(xs, ys) == pytest.approx(0.0)

    def test_linear_slope_positive(self):
        xs = [0.0, 1.0, 2.0]
        ys = [0.0, 1.0, 2.0]
        assert _linear_slope(xs, ys) == pytest.approx(1.0)

    def test_linear_slope_negative(self):
        xs = [0.0, 1.0, 2.0]
        ys = [2.0, 1.0, 0.0]
        assert _linear_slope(xs, ys) == pytest.approx(-1.0)

    def test_linear_slope_degenerate_single_point(self):
        assert _linear_slope([1.0], [1.0]) == 0.0

    def test_linear_slope_degenerate_all_same_x(self):
        assert _linear_slope([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]) == 0.0


# ─────────────────────────────────────────────────────────────────────────────
# Synthetic XVGData factory
# ─────────────────────────────────────────────────────────────────────────────

def _make_xvg(time_ps: list[float], values: list[float],
              title: str = "RMSD", name: str = "Backbone") -> XVGData:
    return XVGData(
        title   = title,
        xlabel  = "Time (ps)",
        ylabel  = "RMSD (nm)",
        time_ps = time_ps,
        series  = [XVGSeries(name=name, values=values)],
        source  = Path("synthetic.xvg"),
    )


def _linspace(start: float, stop: float, n: int) -> list[float]:
    if n == 1:
        return [start]
    step = (stop - start) / (n - 1)
    return [start + i * step for i in range(n)]


# ─────────────────────────────────────────────────────────────────────────────
# RMSD convergence tests
# ─────────────────────────────────────────────────────────────────────────────

class TestAnalyzeRMSDConvergence:

    def test_converged_flat_plateau(self):
        """Flat RMSD in last 20% → converged."""
        n = 100
        time = _linspace(0, 10_000, n)   # 0-10 ns (in ps)
        # Rise for first 60%, then flat at 0.2 nm
        values = [0.1 + 0.1 * (1 - math.exp(-i / 30)) for i in range(n)]
        data = _make_xvg(time, values)
        res = analyze_rmsd_convergence(data, threshold_nm=0.05)
        assert isinstance(res, ConvergenceResult)
        assert res.converged is True
        assert res.std_last20pct < 0.05
        assert res.verdict != ""

    def test_not_converged_rising(self):
        """Monotonically rising RMSD → not converged with strict threshold."""
        n = 100
        time   = _linspace(0, 10_000, n)
        values = _linspace(0.0, 1.0, n)   # linear rise 0→1 nm
        data   = _make_xvg(time, values)
        # Use threshold 0.01 nm — std of last 20% of a linear ramp over 100 points
        # is ~0.06 nm, which is above 0.01 nm → not converged
        res = analyze_rmsd_convergence(data, threshold_nm=0.01)
        assert res.converged is False
        assert "NOT" in res.verdict

    def test_empty_data_not_converged(self):
        """No data → not converged, no crash."""
        data = XVGData(title="", xlabel="", ylabel="", time_ps=[], series=[], source=Path("x.xvg"))
        res  = analyze_rmsd_convergence(data)
        assert res.converged is False
        assert res.verdict != ""

    def test_no_series_not_converged(self):
        """Data with empty series list → not converged."""
        data = XVGData(title="", xlabel="", ylabel="", time_ps=[1.0, 2.0], series=[], source=Path("x.xvg"))
        res  = analyze_rmsd_convergence(data)
        assert res.converged is False

    def test_threshold_respected(self):
        """Strict threshold rejects a slightly noisy plateau."""
        n      = 50
        time   = _linspace(0, 5_000, n)
        values = [0.2 + 0.1 * math.sin(i) for i in range(n)]   # oscillating
        data   = _make_xvg(time, values)

        # With strict threshold → not converged
        res_strict = analyze_rmsd_convergence(data, threshold_nm=0.01)
        # With loose threshold → converged
        res_loose  = analyze_rmsd_convergence(data, threshold_nm=1.0)
        assert res_strict.converged is False
        assert res_loose.converged is True

    def test_plateau_ns_reported_when_converged(self):
        n    = 100
        time = _linspace(0, 10_000, n)
        vals = [0.5 if i < 50 else 0.2 for i in range(n)]  # step drop at halfway
        data = _make_xvg(time, vals)
        res  = analyze_rmsd_convergence(data, threshold_nm=0.05)
        assert res.converged is True
        assert res.plateau_ns is not None
        assert res.plateau_ns >= 0.0

    def test_drift_near_zero_for_flat_series(self):
        n    = 50
        time = _linspace(0, 5_000, n)
        vals = [0.25] * n   # perfectly flat
        data = _make_xvg(time, vals)
        res  = analyze_rmsd_convergence(data)
        assert res.converged is True
        assert abs(res.drift) < 1e-9

    def test_result_fields_present(self):
        """All fields of ConvergenceResult are accessible."""
        data = _make_xvg([0.0, 1000.0], [0.1, 0.15])
        res  = analyze_rmsd_convergence(data)
        _ = res.converged
        _ = res.plateau_ns
        _ = res.mean_last20pct
        _ = res.std_last20pct
        _ = res.drift
        _ = res.verdict


# ─────────────────────────────────────────────────────────────────────────────
# Energy stability tests
# ─────────────────────────────────────────────────────────────────────────────

class TestAnalyzeEnergyStability:

    def test_stable_flat_energy(self):
        """Perfectly flat potential energy → stable."""
        n    = 100
        time = _linspace(0, 10_000, n)
        vals = [-50_000.0] * n
        data = _make_xvg(time, vals, title="Potential Energy", name="E_pot")
        res  = analyze_energy_stability(data)
        assert isinstance(res, StabilityResult)
        assert res.stable is True
        assert abs(res.drift_per_ns) < 1e-6

    def test_unstable_large_drift(self):
        """Large linear drift → not stable."""
        n    = 100
        time = _linspace(0, 10_000, n)          # ps
        vals = _linspace(-50_000, -40_000, n)   # 10% drift over 10 ns
        data = _make_xvg(time, vals, title="Potential Energy")
        res  = analyze_energy_stability(data)
        assert res.stable is False
        assert "NOT" in res.verdict

    def test_empty_data_not_stable(self):
        data = XVGData(title="", xlabel="", ylabel="", time_ps=[], series=[], source=Path("e.xvg"))
        res  = analyze_energy_stability(data)
        assert res.stable is False

    def test_no_series_not_stable(self):
        data = XVGData(title="", xlabel="", ylabel="", time_ps=[1.0], series=[], source=Path("e.xvg"))
        res  = analyze_energy_stability(data)
        assert res.stable is False

    def test_mean_and_std_computed(self):
        vals = [-50_000.0, -50_100.0, -49_900.0]
        time = [0.0, 1000.0, 2000.0]
        data = _make_xvg(time, vals)
        res  = analyze_energy_stability(data)
        assert res.mean == pytest.approx(sum(vals) / 3, rel=1e-6)
        assert res.std  > 0.0

    def test_verdict_not_empty(self):
        vals = [-50_000.0] * 10
        time = _linspace(0, 1_000, 10)
        data = _make_xvg(time, vals)
        res  = analyze_energy_stability(data)
        assert res.verdict != ""

    def test_drift_per_ns_sign_positive(self):
        """Rising energy → positive drift_per_ns."""
        n    = 20
        time = _linspace(0, 2_000, n)
        vals = _linspace(-50_000, -49_000, n)   # +1000 kJ/mol in 2 ns → ~500/ns
        data = _make_xvg(time, vals)
        res  = analyze_energy_stability(data)
        assert res.drift_per_ns > 0

    def test_drift_per_ns_sign_negative(self):
        """Falling energy → negative drift_per_ns."""
        n    = 20
        time = _linspace(0, 2_000, n)
        vals = _linspace(-49_000, -50_000, n)
        data = _make_xvg(time, vals)
        res  = analyze_energy_stability(data)
        assert res.drift_per_ns < 0
