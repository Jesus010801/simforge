"""Phase A — unit tests for runtime/quality_classifier.py."""
from __future__ import annotations

import math
from pathlib import Path

import pytest

from runtime.xvg_parser import XVGData, XVGSeries
from runtime.quality_classifier import (
    RunQuality,
    QualityReport,
    classify_run,
)


# ─────────────────────────────────────────────────────────────────────────────
# Helpers to build synthetic XVGData
# ─────────────────────────────────────────────────────────────────────────────

def _make_xvg(
    values: list[float],
    duration_ns: float = 10.0,
    title: str = "",
    ylabel: str = "",
) -> XVGData:
    """Build a minimal XVGData with evenly-spaced time axis."""
    n = len(values)
    step_ps = (duration_ns * 1000.0) / max(n - 1, 1)
    time_ps = [i * step_ps for i in range(n)]
    return XVGData(
        title   = title,
        xlabel  = "Time (ps)",
        ylabel  = ylabel,
        time_ps = time_ps,
        series  = [XVGSeries(name="col1", values=values)],
        source  = Path("synthetic.xvg"),
    )


def _flat(value: float, n: int = 200, duration_ns: float = 10.0) -> XVGData:
    """Flat signal: converged RMSD or stable energy."""
    import random
    rng = random.Random(42)
    vals = [value + rng.gauss(0, 0.001) for _ in range(n)]
    return _make_xvg(vals, duration_ns=duration_ns)


def _growing(start: float = 0.1, end: float = 0.8,
             n: int = 200, duration_ns: float = 10.0) -> XVGData:
    """Linearly growing signal from start to end nm.

    The last-20% window spans [0.8*range, range], giving enough std
    to hit the NOT_CONVERGED threshold (std >= 0.30 nm) when the
    full ramp is spread across 100% of the time axis.
    The trick: use a wide enough total swing so that the last-20% window
    itself has std >= 0.30 nm via the linear ramp alone (no noise needed).

    For a uniform ramp, std(last 20%) ≈ (0.2 * total_swing) / sqrt(12) ≈ 0.0577 * swing.
    To get std >= 0.30 nm → swing >= 5.2 nm — too large.
    Instead we use a saw-wave in the last 20%: alternate high/low values.
    """
    # Build a ramp for 80%, then oscillate wildly (but below 1 nm) for last 20%
    ramp_n  = int(n * 0.8)
    wave_n  = n - ramp_n
    # Ramp section
    ramp = [start + (0.5 - start) * i / max(ramp_n - 1, 1) for i in range(ramp_n)]
    # Oscillating section in last 20%: alternates between 0.1 and 0.9 nm
    wave = [0.9 if i % 2 == 0 else 0.1 for i in range(wave_n)]
    vals = ramp + wave
    return _make_xvg(vals, duration_ns=duration_ns)


def _short(value: float, duration_ns: float = 0.5) -> XVGData:
    """Short signal with < 1 ns data."""
    return _flat(value, n=50, duration_ns=duration_ns)


def _explosive(n: int = 200, duration_ns: float = 10.0) -> XVGData:
    """RMSD that suddenly jumps above 1 nm."""
    vals = [0.2] * (n // 2) + [1.5] * (n - n // 2)
    return _make_xvg(vals, duration_ns=duration_ns)


def _energy_stable(mean: float = -50000.0, n: int = 200, duration_ns: float = 10.0) -> XVGData:
    """Energy with small fluctuations (< 1% std)."""
    import random
    rng = random.Random(99)
    std = abs(mean) * 0.005   # 0.5% — below threshold
    vals = [mean + rng.gauss(0, std) for _ in range(n)]
    return _make_xvg(vals, duration_ns=duration_ns, ylabel="kJ/mol")


def _energy_drifting(mean: float = -50000.0, n: int = 200, duration_ns: float = 10.0) -> XVGData:
    """Energy with > 5% drift per ns (problematic)."""
    # drift of 10% of |mean| per ns over 10 ns → huge
    drift_total = abs(mean) * 0.10 * duration_ns
    vals = [mean + drift_total * i / (n - 1) for i in range(n)]
    return _make_xvg(vals, duration_ns=duration_ns, ylabel="kJ/mol")


def _energy_noisy(mean: float = -50000.0, n: int = 200, duration_ns: float = 10.0) -> XVGData:
    """Energy with very large fluctuations (std > 10% of mean) → PROBLEMATIC."""
    import random
    rng = random.Random(7)
    std = abs(mean) * 0.15   # 15%
    vals = [mean + rng.gauss(0, std) for _ in range(n)]
    return _make_xvg(vals, duration_ns=duration_ns, ylabel="kJ/mol")


# ─────────────────────────────────────────────────────────────────────────────
# Tests
# ─────────────────────────────────────────────────────────────────────────────

class TestConverged:
    def test_low_std_zero_drift_gives_converged(self):
        rmsd = _flat(0.2)
        report = classify_run(rmsd_data=rmsd)
        assert report.quality == RunQuality.CONVERGED

    def test_both_converged_gives_converged(self):
        rmsd   = _flat(0.2)
        energy = _energy_stable()
        report = classify_run(rmsd_data=rmsd, energy_data=energy)
        assert report.quality == RunQuality.CONVERGED

    def test_converged_confidence_above_threshold(self):
        rmsd = _flat(0.2)
        report = classify_run(rmsd_data=rmsd)
        assert report.confidence >= 0.85

    def test_converged_long_sim_high_confidence(self):
        rmsd = _flat(0.2, n=1000, duration_ns=100.0)
        energy = _energy_stable(n=1000, duration_ns=100.0)
        report = classify_run(rmsd_data=rmsd, energy_data=energy)
        assert report.quality == RunQuality.CONVERGED
        assert report.confidence >= 0.85

    def test_evidence_non_empty_when_converged(self):
        rmsd = _flat(0.2)
        report = classify_run(rmsd_data=rmsd)
        assert len(report.evidence) > 0


class TestNotConverged:
    def test_growing_rmsd_gives_not_converged(self):
        # RMSD grows from 0.1 to 0.8 nm — clear drift, std >= 0.30
        rmsd = _growing(0.1, 0.8)
        report = classify_run(rmsd_data=rmsd)
        assert report.quality == RunQuality.NOT_CONVERGED

    def test_not_converged_confidence_high(self):
        rmsd = _growing(0.1, 0.8)
        report = classify_run(rmsd_data=rmsd)
        assert report.confidence >= 0.7

    def test_evidence_non_empty_when_not_converged(self):
        rmsd = _growing(0.1, 0.8)
        report = classify_run(rmsd_data=rmsd)
        assert len(report.evidence) > 0


class TestProblematic:
    def test_energy_explosion_std_gives_problematic(self):
        energy = _energy_noisy()
        report = classify_run(energy_data=energy)
        assert report.quality == RunQuality.PROBLEMATIC

    def test_energy_large_drift_gives_problematic(self):
        energy = _energy_drifting()
        report = classify_run(energy_data=energy)
        assert report.quality == RunQuality.PROBLEMATIC

    def test_rmsd_above_1nm_gives_problematic(self):
        rmsd = _explosive()
        report = classify_run(rmsd_data=rmsd)
        assert report.quality == RunQuality.PROBLEMATIC

    def test_problematic_has_warnings(self):
        energy = _energy_noisy()
        report = classify_run(energy_data=energy)
        assert len(report.warnings) > 0

    def test_problematic_confidence_high(self):
        energy = _energy_noisy()
        report = classify_run(energy_data=energy)
        assert report.confidence >= 0.8


class TestInsufficientData:
    def test_short_rmsd_gives_insufficient(self):
        rmsd = _short(0.2, duration_ns=0.5)
        report = classify_run(rmsd_data=rmsd)
        assert report.quality == RunQuality.INSUFFICIENT_DATA

    def test_no_data_gives_insufficient(self):
        report = classify_run()
        assert report.quality == RunQuality.INSUFFICIENT_DATA

    def test_none_data_gives_insufficient(self):
        report = classify_run(rmsd_data=None, energy_data=None)
        assert report.quality == RunQuality.INSUFFICIENT_DATA

    def test_insufficient_confidence_is_one(self):
        report = classify_run()
        assert report.confidence == 1.0


class TestPartiallyConverged:
    def test_rmsd_converged_energy_drifting_gives_partially(self):
        # RMSD flat, energy has large std (but not explosion-level drift)
        rmsd   = _flat(0.2)
        # Energy with ~2% std — above 1% threshold but below 10%
        import random
        rng = random.Random(5)
        mean = -50000.0
        vals = [mean + rng.gauss(0, abs(mean) * 0.02) for _ in range(200)]
        energy_xvg = _make_xvg(vals, duration_ns=10.0)
        report = classify_run(rmsd_data=rmsd, energy_data=energy_xvg)
        assert report.quality == RunQuality.PARTIALLY_CONVERGED

    def test_only_rmsd_partial_drift_gives_partially(self):
        # RMSD drifts slowly: std in last 20% is 0.10–0.29 nm with slight drift
        # Use a signal that has moderate fluctuation but stays below 0.30 std
        import random
        rng = random.Random(13)
        # slow rise from 0.1 to 0.25, then noisy plateau
        vals = (
            [0.1 + 0.15 * i / 100 for i in range(100)]
            + [0.25 + rng.gauss(0, 0.04) for _ in range(100)]
        )
        rmsd = _make_xvg(vals, duration_ns=10.0)
        report = classify_run(rmsd_data=rmsd)
        # Should not be CONVERGED (too much drift/std) and not PROBLEMATIC
        assert report.quality in (RunQuality.PARTIALLY_CONVERGED, RunQuality.NOT_CONVERGED)

    def test_partial_evidence_non_empty(self):
        rmsd = _flat(0.2)
        import random
        rng = random.Random(5)
        mean = -50000.0
        vals = [mean + rng.gauss(0, abs(mean) * 0.02) for _ in range(200)]
        energy_xvg = _make_xvg(vals, duration_ns=10.0)
        report = classify_run(rmsd_data=rmsd, energy_data=energy_xvg)
        assert len(report.evidence) > 0


class TestPartialData:
    def test_only_rmsd_converged(self):
        rmsd   = _flat(0.2)
        report = classify_run(rmsd_data=rmsd, energy_data=None)
        assert report.quality == RunQuality.CONVERGED

    def test_only_energy_stable(self):
        energy = _energy_stable()
        report = classify_run(rmsd_data=None, energy_data=energy)
        # With only stable energy and no RMSD, partially converged
        # (can't confirm full convergence without RMSD)
        assert report.quality in (RunQuality.CONVERGED, RunQuality.PARTIALLY_CONVERGED)

    def test_only_rmsd_not_converged(self):
        rmsd = _growing(0.1, 0.8)
        report = classify_run(rmsd_data=rmsd, energy_data=None)
        assert report.quality == RunQuality.NOT_CONVERGED


class TestConfidenceInvariant:
    """Confidence must always be in [0.0, 1.0]."""

    def _check(self, rmsd=None, energy=None):
        r = classify_run(rmsd_data=rmsd, energy_data=energy)
        assert 0.0 <= r.confidence <= 1.0, f"confidence={r.confidence} out of range"

    def test_no_data(self):
        self._check()

    def test_converged(self):
        self._check(rmsd=_flat(0.2))

    def test_not_converged(self):
        self._check(rmsd=_growing(0.1, 0.8))

    def test_problematic_rmsd(self):
        self._check(rmsd=_explosive())

    def test_problematic_energy(self):
        self._check(energy=_energy_noisy())

    def test_both_converged(self):
        self._check(rmsd=_flat(0.2), energy=_energy_stable())

    def test_both_problematic(self):
        self._check(rmsd=_explosive(), energy=_energy_noisy())

    def test_short_data(self):
        self._check(rmsd=_short(0.2))

    def test_100ns(self):
        self._check(rmsd=_flat(0.2, n=1000, duration_ns=100.0))


class TestOutputMethods:
    def test_as_dict_keys(self):
        report = classify_run()
        d = report.as_dict()
        assert "quality" in d
        assert "confidence" in d
        assert "evidence" in d
        assert "warnings" in d
        assert "recommendations" in d
        assert "metrics" in d

    def test_as_markdown_returns_string(self):
        report = classify_run(rmsd_data=_flat(0.2))
        md = report.as_markdown()
        assert isinstance(md, str)
        assert len(md) > 0

    def test_as_dict_quality_is_string(self):
        report = classify_run()
        d = report.as_dict()
        assert isinstance(d["quality"], str)


# ─────────────────────────────────────────────────────────────────────────────
# XVG time-unit handling
# ─────────────────────────────────────────────────────────────────────────────

def _make_xvg_with_unit(
    values: list[float],
    x_values: list[float],
    x_unit: "str | None",
    xlabel: "str | None" = None,
) -> XVGData:
    """Build XVGData with explicit x-axis values and unit declaration.

    xlabel=None → auto-generate from x_unit; xlabel="" → no label at all.
    """
    if xlabel is None:
        xlabel = f"Time ({x_unit})" if x_unit else ""
    return XVGData(
        title   = "RMSD",
        xlabel  = xlabel,
        ylabel  = "RMSD (nm)",
        time_ps = x_values,
        series  = [XVGSeries(name="Backbone", values=values)],
        source  = Path("test.xvg"),
        x_unit  = x_unit,
    )


def _flat_values(n: int = 200, value: float = 0.15) -> list[float]:
    import random
    rng = random.Random(7)
    return [value + rng.gauss(0, 0.005) for _ in range(n)]


class TestXVGUnitHandling:
    """Requirement: x-axis unit read from XVG header; conversions applied correctly."""

    def test_100_ps_is_insufficient_data(self):
        """100 ps = 0.1 ns → below 1 ns threshold → INSUFFICIENT_DATA."""
        x = [i * 1.0 for i in range(101)]   # 0..100 ps
        data = _make_xvg_with_unit(_flat_values(101), x, x_unit="ps")
        report = classify_run(rmsd_data=data)
        assert report.quality == RunQuality.INSUFFICIENT_DATA

    def test_100_ns_passes_duration_threshold(self):
        """100 ns is well above 1 ns threshold → not INSUFFICIENT_DATA."""
        x = [i * 1.0 for i in range(101)]   # 0..100 ns
        data = _make_xvg_with_unit(_flat_values(101), x, x_unit="ns")
        report = classify_run(rmsd_data=data)
        assert report.quality != RunQuality.INSUFFICIENT_DATA

    def test_0_1_us_passes_duration_threshold(self):
        """0.1 µs = 100 ns → well above 1 ns threshold → not INSUFFICIENT_DATA."""
        x = [i * 0.001 for i in range(101)]   # 0..0.1 µs
        data = _make_xvg_with_unit(_flat_values(101), x, x_unit="us")
        report = classify_run(rmsd_data=data)
        assert report.quality != RunQuality.INSUFFICIENT_DATA

    def test_ns_unit_total_ns_metric_correct(self):
        """When x_unit=ns, rmsd_total_ns metric should equal the last x value."""
        x = [float(i) for i in range(101)]   # 0..100 ns
        data = _make_xvg_with_unit(_flat_values(101), x, x_unit="ns")
        report = classify_run(rmsd_data=data)
        assert report.metrics.get("rmsd_total_ns") == pytest.approx(100.0, abs=0.1)

    def test_ps_unit_total_ns_metric_correct(self):
        """When x_unit=ps, 1000 ps → rmsd_total_ns = 1.0."""
        x = [i * 10.0 for i in range(101)]   # 0..1000 ps
        data = _make_xvg_with_unit(_flat_values(101), x, x_unit="ps")
        report = classify_run(rmsd_data=data)
        assert report.metrics.get("rmsd_total_ns") == pytest.approx(1.0, abs=0.05)

    def test_missing_unit_assumes_ps_warns(self):
        """No unit declared (x_unit=None) + undeclared xlabel → warning emitted."""
        x = [i * 10.0 for i in range(101)]   # 0..1000 (assumed ps)
        data = _make_xvg_with_unit(
            _flat_values(101), x,
            x_unit=None,
            xlabel="Time",   # non-empty but no unit in parens → triggers warning
        )
        report = classify_run(rmsd_data=data)
        unit_warnings = [w for w in report.warnings if "assumed ps" in w.lower()]
        assert len(unit_warnings) >= 1

    def test_missing_unit_no_warning_when_xlabel_empty(self):
        """No xlabel at all → no unit warning (can't know intent, silently assume ps)."""
        x = [i * 10.0 for i in range(101)]   # 0..1000 ps
        data = _make_xvg_with_unit(
            _flat_values(101), x,
            x_unit=None,
            xlabel="",
        )
        report = classify_run(rmsd_data=data)
        unit_warnings = [w for w in report.warnings if "assumed ps" in w.lower()]
        assert len(unit_warnings) == 0

    def test_no_warning_when_unit_declared(self):
        """Explicit x_unit set → no unit assumption warning."""
        x = [float(i) for i in range(101)]
        data = _make_xvg_with_unit(_flat_values(101), x, x_unit="ns")
        report = classify_run(rmsd_data=data)
        unit_warnings = [w for w in report.warnings if "assumed ps" in w.lower()]
        assert len(unit_warnings) == 0
