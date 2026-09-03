"""FEL-specific XVG reader.

Wraps runtime.xvg_parser but adds:
  - 1-based column selection with clear error messages
  - Strict error when the file is empty or column is out of range
  - Skipped-line counting for QC warnings
  - Returns analysis.fel.models.XVGSeries
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

from analysis.fel.models import XVGSeries, FELWarning

# ── header pattern constants (mirror runtime/xvg_parser.py) ──────────────────
_TITLE_RE  = re.compile(r'@\s+title\s+"(.+)"',        re.IGNORECASE)
_XLABEL_RE = re.compile(r'@\s+xaxis\s+label\s+"(.+)"', re.IGNORECASE)
_YLABEL_RE = re.compile(r'@\s+yaxis\s+label\s+"(.+)"', re.IGNORECASE)
_LEGEND_RE = re.compile(r'@\s+s(\d+)\s+legend\s+"(.+)"', re.IGNORECASE)

_NS_FACTORS: dict[str, float] = {"ps": 1e-3, "ns": 1.0, "us": 1e3}


def _detect_x_unit(xlabel: str) -> Optional[str]:
    lower = xlabel.lower()
    for token, unit in [("(us)", "us"), ("(µs)", "us"), ("(μs)", "us"),
                        ("(ns)", "ns"), ("(ps)", "ps")]:
        if token in lower:
            return unit
    return None


def _detect_data_unit(ylabel: str) -> Optional[str]:
    """Extract the physical unit from a Y-axis label, e.g. 'RMSD (nm)' → 'nm'."""
    m = re.search(r'\(([^)]+)\)', ylabel)
    if m:
        candidate = m.group(1).strip()
        if 0 < len(candidate) < 20:  # skip long parenthetical notes
            return candidate
    return None


def default_unit_for_label(label: str) -> Optional[str]:
    """Return 'nm' for RMSD/Rg/RMSF family labels when no unit was detected in the XVG.

    Uses word-boundary-aware matching so that 'energy' does not match 'rg'.
    """
    lower = label.lower().strip()
    # Exact abbreviated names
    _NM_EXACT = {"rmsd", "rg", "rmsf"}
    # Prefixes / substrings that unambiguously identify nm-scale observables
    _NM_PREFIX = ("rmsd", "rmsf", "radius")
    _NM_CONTAINS = ("gyration",)

    if lower in _NM_EXACT:
        return "nm"
    if any(lower.startswith(p) for p in _NM_PREFIX):
        return "nm"
    if any(kw in lower for kw in _NM_CONTAINS):
        return "nm"
    # Rg as a standalone token (not a substring of longer words like "energy")
    if re.search(r'\brg\b', lower):
        return "nm"
    return None


def load_xvg_series(
    path: Path,
    column: int = 1,
    *,
    override_unit: Optional[str] = None,
) -> XVGSeries:
    """Parse an XVG file and extract one data column.

    Parameters
    ----------
    path:
        Path to the .xvg file.
    column:
        1-based data column index (1 = first column after the time column).
    override_unit:
        Force time unit to "ps", "ns", or "us" instead of auto-detecting.

    Raises
    ------
    FileNotFoundError
        When the file does not exist.
    ValueError
        When the file has no numeric data, or when `column` is out of range.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"XVG file not found: {path}")

    title    = ""
    xlabel   = ""
    ylabel   = ""
    legends: dict[int, str] = {}
    time_raw: list[float]   = []
    col_data: dict[int, list[float]] = {}  # 0-based data col → values
    n_skipped = 0

    with path.open(errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("@"):
                m = _TITLE_RE.search(line)
                if m:
                    title = m.group(1); continue
                m = _XLABEL_RE.search(line)
                if m:
                    xlabel = m.group(1); continue
                m = _YLABEL_RE.search(line)
                if m:
                    ylabel = m.group(1); continue
                m = _LEGEND_RE.search(line)
                if m:
                    legends[int(m.group(1))] = m.group(2)
                continue
            # data row
            parts = line.split()
            if len(parts) < 2:
                n_skipped += 1
                continue
            try:
                vals = [float(p) for p in parts]
            except ValueError:
                n_skipped += 1
                continue
            time_raw.append(vals[0])
            for di, v in enumerate(vals[1:]):
                col_data.setdefault(di, []).append(v)

    if not time_raw:
        raise ValueError(
            f"No numeric data rows found in {path}. "
            "File may be empty, malformed, or contain only header lines."
        )

    n_data_cols = len(col_data)
    if column < 1 or column > n_data_cols:
        available = list(range(1, n_data_cols + 1))
        names = [legends.get(i, f"col{i+1}") for i in range(n_data_cols)]
        raise ValueError(
            f"Column {column} not found in {path}. "
            f"Available columns (1-based): {available} — names: {names}"
        )

    col_idx_0 = column - 1
    col_name  = legends.get(col_idx_0, f"col{column}")
    values    = col_data[col_idx_0]

    x_unit = override_unit or _detect_x_unit(xlabel)
    data_unit = _detect_data_unit(ylabel)

    return XVGSeries(
        path         = path,
        title        = title,
        xlabel       = xlabel,
        ylabel       = ylabel,
        x_unit       = x_unit,
        time_raw     = time_raw,
        data_values  = values,
        column_index = column,
        column_name  = col_name,
        n_skipped    = n_skipped,
        data_unit    = data_unit,
    )


# ── Feature auto-discovery from a run directory ───────────────────────────────

# Maps feature keyword → ordered list of candidate filenames to probe
_FEATURE_CANDIDATES: dict[str, list[str]] = {
    "rmsd":   ["rmsd.xvg"],
    "rg":     ["gyrate.xvg", "Rg.xvg", "rg.xvg"],
    "rmsf":   ["rmsf.xvg"],
    "sasa":   ["sasa.xvg"],
    "hbonds": ["hbnum.xvg"],
    "energy": ["energy.xvg"],
}

_FEATURE_LABELS: dict[str, str] = {
    "rmsd":   "RMSD (nm)",
    "rg":     "Rg (nm)",
    "rmsf":   "RMSF (nm)",
    "sasa":   "SASA (nm²)",
    "hbonds": "H-bonds",
    "energy": "Potential Energy (kJ/mol)",
}


def discover_feature_xvg(run_dir: Path, feature: str) -> Optional[Path]:
    """Return the first matching XVG candidate for a feature keyword in run_dir."""
    for name in _FEATURE_CANDIDATES.get(feature.lower(), []):
        p = run_dir / name
        if p.exists():
            return p
    return None


def feature_label(feature: str) -> str:
    return _FEATURE_LABELS.get(feature.lower(), feature)
