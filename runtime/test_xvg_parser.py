"""Tests for runtime/xvg_parser.py — no GROMACS required."""
from __future__ import annotations

import textwrap
from pathlib import Path

import pytest

from runtime.xvg_parser import parse_xvg, XVGData, XVGSeries


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
