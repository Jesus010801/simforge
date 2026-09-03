"""Data models for FEL analysis.

All models are plain dataclasses to avoid extra dependencies.
Every model that is written to disk implements to_dict() for JSON serialization.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class FELWarning:
    code: str
    message: str
    severity: str = "warn"  # warn | info | error

    def to_dict(self) -> dict:
        return {"code": self.code, "message": self.message, "severity": self.severity}


@dataclass
class XVGSeries:
    """A single data series extracted from a GROMACS XVG file."""
    path: Path
    title: str
    xlabel: str
    ylabel: str
    x_unit: Optional[str]      # "ps" | "ns" | "us" | None (unknown, assumed ps)
    time_raw: list[float]      # raw time values in x_unit
    data_values: list[float]   # extracted data column
    column_index: int          # 1-based user column index
    column_name: str           # series name from @ sN legend or auto-generated
    n_skipped: int = 0         # malformed data lines skipped during parsing
    data_unit: Optional[str] = None  # feature data unit detected from ylabel, e.g. "nm"

    def n_frames(self) -> int:
        return len(self.time_raw)

    def time_ns(self) -> list[float]:
        """Convert raw time to nanoseconds."""
        _NS_FACTORS = {"ps": 1e-3, "ns": 1.0, "us": 1e3}
        factor = _NS_FACTORS.get(self.x_unit or "ps", 1e-3)
        return [t * factor for t in self.time_raw]

    def to_dict(self) -> dict:
        return {
            "path": str(self.path),
            "title": self.title,
            "xlabel": self.xlabel,
            "ylabel": self.ylabel,
            "x_unit": self.x_unit,
            "n_frames": self.n_frames(),
            "column_index": self.column_index,
            "column_name": self.column_name,
            "n_skipped_lines": self.n_skipped,
        }


@dataclass
class FeatureTable:
    """Time-aligned two-feature table ready for histogram computation."""
    time_ns: list[float]
    x: list[float]
    y: list[float]
    x_label: str
    y_label: str
    x_unit_original: Optional[str]
    y_unit_original: Optional[str]
    time_conversion_note: str
    n_frames: int
    warnings: list[FELWarning] = field(default_factory=list)
    x_unit: Optional[str] = None           # feature data unit, e.g. "nm"
    y_unit: Optional[str] = None
    time_unit_normalized: str = "ns"
    time_conversion_x_to_ns: float = 1e-3  # factor applied (1e-3 for ps, 1.0 for ns)
    time_conversion_y_to_ns: float = 1e-3

    def to_csv_lines(self) -> list[str]:
        lines = [f"time_ns,{self.x_label},{self.y_label}"]
        for t, x, y in zip(self.time_ns, self.x, self.y):
            lines.append(f"{t:.6g},{x:.6g},{y:.6g}")
        return lines


@dataclass
class FELConfig:
    """Full resolved configuration for one FEL run."""
    xvg_x: Path
    xvg_y: Path
    x_column: int
    y_column: int
    x_label: str
    y_label: str
    temperature_K: float
    n_bins_x: int
    n_bins_y: int
    out_dir: Path
    time_unit_x: str          # "auto" | "ps" | "ns"
    time_unit_y: str
    time_tolerance: float
    empty_bin_policy: str     # "nan" | "cap"
    energy_cap_kj: float
    dry_run: bool
    x_unit: Optional[str] = None  # populated in workflow after XVG parsing
    y_unit: Optional[str] = None

    def to_dict(self) -> dict:
        return {
            "xvg_x": str(self.xvg_x),
            "xvg_y": str(self.xvg_y),
            "x_column": self.x_column,
            "y_column": self.y_column,
            "x_label": self.x_label,
            "y_label": self.y_label,
            "temperature_K": self.temperature_K,
            "n_bins_x": self.n_bins_x,
            "n_bins_y": self.n_bins_y,
            "out_dir": str(self.out_dir),
            "time_unit_x": self.time_unit_x,
            "time_unit_y": self.time_unit_y,
            "time_tolerance": self.time_tolerance,
            "empty_bin_policy": self.empty_bin_policy,
            "energy_cap_kj": self.energy_cap_kj,
            "dry_run": self.dry_run,
            "x_unit": self.x_unit,
            "y_unit": self.y_unit,
        }


@dataclass
class FELSurface:
    """Computed 2D FEL surface."""
    x_edges: list[float]
    y_edges: list[float]
    x_centers: list[float]
    y_centers: list[float]
    # 2D grids stored row-major [ix][iy]: outer list = x bins, inner = y bins
    counts: list[list[float]]
    probability: list[list[float]]   # NaN sentinel: float('nan')
    free_energy: list[list[float]]   # kJ/mol; NaN for empty bins
    n_frames_used: int
    temperature_K: float
    kB_kj_mol_K: float
    x_label: str
    y_label: str
    formula: str
    warnings: list[FELWarning] = field(default_factory=list)
    x_unit: Optional[str] = None  # feature data unit, e.g. "nm"
    y_unit: Optional[str] = None

    def to_dict(self) -> dict:
        return {
            "n_frames_used": self.n_frames_used,
            "temperature_K": self.temperature_K,
            "kB_kj_mol_K": self.kB_kj_mol_K,
            "formula": self.formula,
            "x_label": self.x_label,
            "y_label": self.y_label,
            "n_bins_x": len(self.x_centers),
            "n_bins_y": len(self.y_centers),
            "x_range": [self.x_edges[0], self.x_edges[-1]],
            "y_range": [self.y_edges[0], self.y_edges[-1]],
            "x_unit": self.x_unit,
            "y_unit": self.y_unit,
        }

    def _flat_csv(self, values: list[list[float]], header: str) -> list[str]:
        """Produce flat x_center,y_center,<header> CSV lines."""
        lines = [f"x_{self.x_label},y_{self.y_label},{header}"]
        for ix, xc in enumerate(self.x_centers):
            for iy, yc in enumerate(self.y_centers):
                v = values[ix][iy]
                v_str = "NaN" if v != v else f"{v:.6g}"  # NaN check: nan != nan
                lines.append(f"{xc:.6g},{yc:.6g},{v_str}")
        return lines

    def counts_csv(self) -> list[str]:
        return self._flat_csv(self.counts, "counts")

    def probability_csv(self) -> list[str]:
        return self._flat_csv(self.probability, "probability")

    def free_energy_csv(self) -> list[str]:
        return self._flat_csv(self.free_energy, "free_energy_kJ_mol")


@dataclass
class FELRunResult:
    config: FELConfig
    series_x: Optional[XVGSeries] = None
    series_y: Optional[XVGSeries] = None
    features: Optional[FeatureTable] = None
    surface: Optional[FELSurface] = None
    warnings: list[FELWarning] = field(default_factory=list)
    output_files: list[Path] = field(default_factory=list)


# ── Phase 3: minima and extraction models ─────────────────────────────────────

@dataclass
class FELMinimum:
    """One detected local free-energy minimum."""
    minimum_id: int           # 1-based, set after deduplication
    x_bin: int                # 0-based bin index in the x direction
    y_bin: int
    x: float                  # bin center value
    y: float
    free_energy_kJ_mol: float
    count: Optional[int] = None          # frames in this bin (from histogram)
    probability: Optional[float] = None  # normalised probability

    def to_dict(self) -> dict:
        d: dict = {
            "minimum_id": self.minimum_id,
            "x_bin": self.x_bin,
            "y_bin": self.y_bin,
            "x": self.x,
            "y": self.y,
            "free_energy_kJ_mol": self.free_energy_kJ_mol,
        }
        if self.count is not None:
            d["count"] = self.count
        if self.probability is not None:
            d["probability"] = self.probability
        return d


@dataclass
class FrameMapping:
    """Nearest trajectory frame mapped to one FEL minimum."""
    minimum_id: int
    target_x: float
    target_y: float
    x_unit: Optional[str]
    y_unit: Optional[str]
    free_energy_kJ_mol: float
    selected_time_ns: float
    selected_frame_index: int
    frame_x: float
    frame_y: float
    distance_in_feature_space: float

    def to_dict(self) -> dict:
        return {
            "minimum_id": self.minimum_id,
            "target_x": self.target_x,
            "target_y": self.target_y,
            "x_unit": self.x_unit,
            "y_unit": self.y_unit,
            "free_energy_kJ_mol": self.free_energy_kJ_mol,
            "selected_time_ns": self.selected_time_ns,
            "selected_frame_index": self.selected_frame_index,
            "frame_x": self.frame_x,
            "frame_y": self.frame_y,
            "distance_in_feature_space": self.distance_in_feature_space,
        }


@dataclass
class StateResult:
    """Outcome of planning or executing extraction for one state."""
    minimum_id: int
    planned_commands: list[str]     # shell-ready command strings
    planned_cmd_lists: list[list[str]]  # subprocess argument lists
    output_paths: list[Path]
    status: str                     # "planned" | "completed" | "failed"
    return_codes: list[int] = field(default_factory=list)
    stdout_lines: list[str] = field(default_factory=list)
    stderr_lines: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "minimum_id": self.minimum_id,
            "status": self.status,
            "planned_commands": self.planned_commands,
            "output_paths": [str(p) for p in self.output_paths],
            "return_codes": self.return_codes,
        }


@dataclass
class ExtractConfig:
    """Configuration for one extract-minima run."""
    fel_dir: Path
    trajectory: Path
    topology: Path
    structure: Optional[Path]
    index: Optional[Path]
    selection: str
    n_minima: int
    min_distance_bins: int
    representative: str       # "nearest-frame"
    format: str               # "pdb" | "gro" | "both"
    gmx: str
    out_dir: Path
    dry_run: bool
    execute: bool
    force: bool
    max_delta_g_kj: Optional[float] = 2.5  # None = no energy filter
    min_count: int = 1                      # bins with fewer frames are excluded
    system_name: Optional[str] = None       # explicit override; inferred when None
    compat_receptor_name: bool = False      # write receptor.{ext} alias alongside descriptive name

    def to_dict(self) -> dict:
        return {
            "fel_dir": str(self.fel_dir),
            "trajectory": str(self.trajectory),
            "topology": str(self.topology),
            "structure": str(self.structure) if self.structure else None,
            "index": str(self.index) if self.index else None,
            "selection": self.selection,
            "n_minima": self.n_minima,
            "min_distance_bins": self.min_distance_bins,
            "representative": self.representative,
            "format": self.format,
            "gmx": self.gmx,
            "out_dir": str(self.out_dir),
            "dry_run": self.dry_run,
            "execute": self.execute,
            "force": self.force,
            "max_delta_g_kj": self.max_delta_g_kj,
            "min_count": self.min_count,
            "system_name": self.system_name,
            "compat_receptor_name": self.compat_receptor_name,
        }


@dataclass
class ExtractResult:
    """Full result of an extract-minima run."""
    config: ExtractConfig
    minima: list[FELMinimum] = field(default_factory=list)
    frame_map: list[FrameMapping] = field(default_factory=list)
    warnings: list[FELWarning] = field(default_factory=list)
    state_results: list[StateResult] = field(default_factory=list)
    output_files: list[Path] = field(default_factory=list)
    rejected_minima: list[FELMinimum] = field(default_factory=list)
