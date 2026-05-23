"""Sprint 2 — Parse GROMACS .xvg files.

XVG format produced by gmx rms, gmx energy, gmx rmsf, etc.
  # comment lines
  @ directive lines (title, xlabel, ylabel, legend names)
  data lines — whitespace-separated floats
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path


# Patterns for @ directives
_TITLE_RE   = re.compile(r'@\s+title\s+"(.+)"',   re.IGNORECASE)
_XLABEL_RE  = re.compile(r'@\s+xaxis\s+label\s+"(.+)"', re.IGNORECASE)
_YLABEL_RE  = re.compile(r'@\s+yaxis\s+label\s+"(.+)"', re.IGNORECASE)
_LEGEND_RE  = re.compile(r'@\s+s(\d+)\s+legend\s+"(.+)"', re.IGNORECASE)


@dataclass
class XVGSeries:
    name:   str             # from @ s0 legend "..." or "col<n>"
    values: list[float]     # data column values


@dataclass
class XVGData:
    title:    str
    xlabel:   str
    ylabel:   str
    time_ps:  list[float]     # first column (time or x-axis)
    series:   list[XVGSeries]
    source:   Path

    def times_ns(self) -> list[float]:
        """Convert time_ps to nanoseconds."""
        return [t / 1000.0 for t in self.time_ps]


def parse_xvg(path: Path) -> XVGData:
    """Parse a GROMACS XVG file into XVGData.

    Gracefully handles:
    - malformed / non-float data lines (skipped)
    - missing directives (empty string defaults)
    - files with no data (returns empty series)
    """
    title:   str = ""
    xlabel:  str = ""
    ylabel:  str = ""
    series_names: dict[int, str] = {}   # s-index → legend name
    time_ps: list[float] = []
    columns: dict[int, list[float]] = {}  # col_idx (1-based) → values

    try:
        with path.open(errors="replace") as f:
            for raw_line in f:
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue

                if line.startswith("@"):
                    # Parse metadata directives
                    m = _TITLE_RE.search(line)
                    if m:
                        title = m.group(1)
                        continue
                    m = _XLABEL_RE.search(line)
                    if m:
                        xlabel = m.group(1)
                        continue
                    m = _YLABEL_RE.search(line)
                    if m:
                        ylabel = m.group(1)
                        continue
                    m = _LEGEND_RE.search(line)
                    if m:
                        series_names[int(m.group(1))] = m.group(2)
                    continue

                # Data line
                parts = line.split()
                if len(parts) < 2:
                    continue
                try:
                    vals = [float(p) for p in parts]
                except ValueError:
                    continue   # skip malformed lines

                time_ps.append(vals[0])
                for i, v in enumerate(vals[1:], 1):
                    columns.setdefault(i, []).append(v)

    except OSError:
        # File not readable — return empty data
        pass

    # Build series list ordered by column index
    series: list[XVGSeries] = []
    for i in sorted(columns):
        name = series_names.get(i - 1, f"col{i}")
        series.append(XVGSeries(name=name, values=columns[i]))

    return XVGData(
        title   = title,
        xlabel  = xlabel,
        ylabel  = ylabel,
        time_ps = time_ps,
        series  = series,
        source  = path,
    )
