"""Tests for runtime/xvg_parser.py — no GROMACS required."""
from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from runtime.xvg_parser import parse_xvg, XVGData, XVGSeries, _extract_x_unit


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _write_xvg(tmp_path: Path, content: str, name: str = "test.xvg") -> Path:
    p = tmp_path / name
    p.write_text(textwrap.dedent(content))
    return p


# ─────────────────────────────────────────────────────────────────────────────
# Metadata parsing
# ─────────────────────────────────────────────────────────────────────────────

class TestXVGMetadata:
    def test_title_parsed(self, tmp_path):
        p = _write_xvg(tmp_path, """
            @ title "RMSD"
            0.0  0.1
            1.0  0.2
        """)
        data = parse_xvg(p)
        assert data.title == "RMSD"

    def test_xlabel_ylabel_parsed(self, tmp_path):
        p = _write_xvg(tmp_path, """
            @ xaxis label "Time (ps)"
            @ yaxis label "RMSD (nm)"
            0.0  0.1
        """)
        data = parse_xvg(p)
        assert "Time" in data.xlabel
        assert "RMSD" in data.ylabel

    def test_legend_names_assigned_to_series(self, tmp_path):
        p = _write_xvg(tmp_path, """
            @ s0 legend "Backbone"
            @ s1 legend "All-atom"
            0.0  0.10  0.15
            1.0  0.11  0.16
        """)
        data = parse_xvg(p)
        assert len(data.series) == 2
        assert data.series[0].name == "Backbone"
        assert data.series[1].name == "All-atom"

    def test_no_legend_uses_colN_fallback(self, tmp_path):
        p = _write_xvg(tmp_path, """
            0.0  0.1  0.2
            1.0  0.3  0.4
        """)
        data = parse_xvg(p)
        assert data.series[0].name == "col1"
        assert data.series[1].name == "col2"

    def test_missing_metadata_defaults_to_empty_string(self, tmp_path):
        p = _write_xvg(tmp_path, "0.0  0.1\n1.0  0.2\n")
        data = parse_xvg(p)
        assert data.title  == ""
        assert data.xlabel == ""
        assert data.ylabel == ""


# ─────────────────────────────────────────────────────────────────────────────
# Data parsing
# ─────────────────────────────────────────────────────────────────────────────

class TestXVGData:
    def test_time_axis_parsed(self, tmp_path):
        p = _write_xvg(tmp_path, "0.0 0.1\n100.0 0.2\n200.0 0.3\n")
        data = parse_xvg(p)
        assert data.time_ps == pytest.approx([0.0, 100.0, 200.0])

    def test_single_series_values(self, tmp_path):
        p = _write_xvg(tmp_path, "0.0 0.10\n1.0 0.20\n2.0 0.30\n")
        data = parse_xvg(p)
        assert len(data.series) == 1
        assert data.series[0].values == pytest.approx([0.10, 0.20, 0.30])

    def test_multi_series_values(self, tmp_path):
        p = _write_xvg(tmp_path, "0.0 0.1 0.5\n1.0 0.2 0.6\n")
        data = parse_xvg(p)
        assert len(data.series) == 2
        assert data.series[0].values == pytest.approx([0.1, 0.2])
        assert data.series[1].values == pytest.approx([0.5, 0.6])

    def test_times_ns_conversion(self, tmp_path):
        p = _write_xvg(tmp_path, "0.0 0.1\n1000.0 0.2\n2000.0 0.3\n")
        data = parse_xvg(p)
        assert data.times_ns() == pytest.approx([0.0, 1.0, 2.0])

    def test_comment_lines_skipped(self, tmp_path):
        p = _write_xvg(tmp_path, """
            # This is a comment
            # Another comment
            0.0  0.1
            1.0  0.2
        """)
        data = parse_xvg(p)
        assert len(data.time_ps) == 2

    def test_source_path_stored(self, tmp_path):
        p = _write_xvg(tmp_path, "0.0  0.1\n")
        data = parse_xvg(p)
        assert data.source == p


# ─────────────────────────────────────────────────────────────────────────────
# X-axis unit detection
# ─────────────────────────────────────────────────────────────────────────────

class TestXVGUnitDetection:
    """Unit detection from xlabel and unit-aware times_ns() conversion."""

    # ── _extract_x_unit helper ────────────────────────────────────────────────

    def test_extract_ps(self):
        assert _extract_x_unit('Time (ps)') == "ps"

    def test_extract_ns(self):
        assert _extract_x_unit('Time (ns)') == "ns"

    def test_extract_us_ascii(self):
        assert _extract_x_unit('Time (us)') == "us"

    def test_extract_us_unicode_mu(self):
        assert _extract_x_unit('Time (µs)') == "us"

    def test_extract_none_no_parens(self):
        assert _extract_x_unit('Time') is None

    def test_extract_none_empty(self):
        assert _extract_x_unit('') is None

    def test_extract_case_insensitive(self):
        assert _extract_x_unit('TIME (NS)') == "ns"

    # ── parse_xvg sets x_unit from file ──────────────────────────────────────

    def test_parse_sets_x_unit_ps(self, tmp_path):
        p = _write_xvg(tmp_path, '@ xaxis label "Time (ps)"\n0.0 0.1\n1000.0 0.2\n')
        data = parse_xvg(p)
        assert data.x_unit == "ps"

    def test_parse_sets_x_unit_ns(self, tmp_path):
        p = _write_xvg(tmp_path, '@ xaxis label "Time (ns)"\n0.0 0.1\n100.0 0.2\n')
        data = parse_xvg(p)
        assert data.x_unit == "ns"

    def test_parse_sets_x_unit_us(self, tmp_path):
        p = _write_xvg(tmp_path, '@ xaxis label "Time (µs)"\n0.0 0.1\n0.1 0.2\n')
        data = parse_xvg(p)
        assert data.x_unit == "us"

    def test_parse_x_unit_none_when_no_xlabel(self, tmp_path):
        p = _write_xvg(tmp_path, '0.0 0.1\n1000.0 0.2\n')
        data = parse_xvg(p)
        assert data.x_unit is None

    def test_parse_x_unit_none_when_xlabel_has_no_unit(self, tmp_path):
        p = _write_xvg(tmp_path, '@ xaxis label "Time"\n0.0 0.1\n1000.0 0.2\n')
        data = parse_xvg(p)
        assert data.x_unit is None

    # ── times_ns() unit-aware conversion ─────────────────────────────────────

    def test_times_ns_with_ps_label(self, tmp_path):
        p = _write_xvg(tmp_path, '@ xaxis label "Time (ps)"\n0.0 0.1\n1000.0 0.2\n')
        data = parse_xvg(p)
        assert data.times_ns() == pytest.approx([0.0, 1.0])

    def test_times_ns_with_ns_label(self, tmp_path):
        """100 ns in file → 100 ns after conversion (no scaling)."""
        p = _write_xvg(tmp_path, '@ xaxis label "Time (ns)"\n0.0 0.1\n100.0 0.2\n')
        data = parse_xvg(p)
        assert data.times_ns() == pytest.approx([0.0, 100.0])

    def test_times_ns_with_us_label(self, tmp_path):
        """0.1 µs in file → 100 ns after conversion."""
        p = _write_xvg(tmp_path, '@ xaxis label "Time (us)"\n0.0 0.1\n0.1 0.2\n')
        data = parse_xvg(p)
        assert data.times_ns() == pytest.approx([0.0, 100.0])

    def test_times_ns_no_unit_assumes_ps(self, tmp_path):
        """Missing xlabel → assume ps → divide by 1000."""
        p = _write_xvg(tmp_path, '0.0 0.1\n1000.0 0.2\n')
        data = parse_xvg(p)
        assert data.x_unit is None
        assert data.times_ns() == pytest.approx([0.0, 1.0])


# ─────────────────────────────────────────────────────────────────────────────
# Edge cases / robustness
# ─────────────────────────────────────────────────────────────────────────────

class TestXVGEdgeCases:
    def test_empty_file_returns_empty_data(self, tmp_path):
        p = _write_xvg(tmp_path, "")
        data = parse_xvg(p)
        assert data.time_ps == []
        assert data.series  == []

    def test_only_comments_returns_empty_data(self, tmp_path):
        p = _write_xvg(tmp_path, "# nothing here\n# nor here\n")
        data = parse_xvg(p)
        assert data.time_ps == []

    def test_malformed_line_skipped(self, tmp_path):
        p = _write_xvg(tmp_path, "0.0  0.1\nBAD LINE\n1.0  0.2\n")
        data = parse_xvg(p)
        assert len(data.time_ps) == 2
        assert data.time_ps == pytest.approx([0.0, 1.0])

    def test_single_column_line_skipped(self, tmp_path):
        """A line with only one value (no y column) should be skipped."""
        p = _write_xvg(tmp_path, "0.0\n1.0  0.2\n")
        data = parse_xvg(p)
        assert len(data.time_ps) == 1
        assert data.time_ps[0] == pytest.approx(1.0)

    def test_missing_file_returns_empty_data(self, tmp_path):
        p = tmp_path / "nonexistent.xvg"
        data = parse_xvg(p)
        assert data.time_ps == []
        assert data.series  == []

    def test_full_xvg_realistic(self, tmp_path):
        """Parse a realistic GROMACS RMSD XVG snippet."""
        content = """\
# This file was created by gmx rms
# gmx rms is part of G R O M A C S
@ title "RMSD"
@ xaxis label "Time (ps)"
@ yaxis label "RMSD (nm)"
@ TYPE xy
@ s0 legend "Backbone"
     0.000     0.000
   100.000     0.094
   200.000     0.112
   300.000     0.108
   400.000     0.115
"""
        p = _write_xvg(tmp_path, content)
        data = parse_xvg(p)
        assert data.title  == "RMSD"
        assert len(data.time_ps) == 5
        assert data.time_ps[0] == pytest.approx(0.0)
        assert data.time_ps[-1] == pytest.approx(400.0)
        assert len(data.series) == 1
        assert data.series[0].name == "Backbone"
        assert data.series[0].values[0] == pytest.approx(0.0, abs=1e-6)
        assert data.series[0].values[1] == pytest.approx(0.094, abs=1e-6)
