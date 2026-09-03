"""Provenance record for a FEL run."""
from __future__ import annotations

import datetime
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from analysis.fel.models import FELConfig, FELRunResult, FELWarning


@dataclass
class FELProvenance:
    timestamp: str
    simforge_version: str
    argv: list[str]
    config: dict
    series_x: Optional[dict]
    series_y: Optional[dict]
    features: Optional[dict]
    surface: Optional[dict]
    warnings: list[dict]

    def to_dict(self) -> dict:
        return {
            "timestamp":        self.timestamp,
            "simforge_version": self.simforge_version,
            "argv":             self.argv,
            "config":           self.config,
            "series_x":         self.series_x,
            "series_y":         self.series_y,
            "features":         self.features,
            "surface":          self.surface,
            "warnings":         self.warnings,
        }


def build_provenance(result: FELRunResult) -> FELProvenance:
    try:
        from simforge import __version__ as _v
        version = _v
    except Exception:
        version = "unknown"

    all_warnings = list(result.warnings)
    if result.features:
        all_warnings += result.features.warnings
    if result.surface:
        all_warnings += result.surface.warnings

    features_meta: Optional[dict] = None
    if result.features:
        ft = result.features
        features_meta = {
            "n_frames":                ft.n_frames,
            "x_label":                 ft.x_label,
            "y_label":                 ft.y_label,
            "x_unit":                  ft.x_unit,
            "y_unit":                  ft.y_unit,
            "x_unit_original":         ft.x_unit_original,
            "y_unit_original":         ft.y_unit_original,
            "time_unit_normalized":    ft.time_unit_normalized,
            "time_conversion_x_to_ns": ft.time_conversion_x_to_ns,
            "time_conversion_y_to_ns": ft.time_conversion_y_to_ns,
            "time_conversion_note":    ft.time_conversion_note,
        }

    return FELProvenance(
        timestamp        = datetime.datetime.utcnow().isoformat() + "Z",
        simforge_version = version,
        argv             = sys.argv[:],
        config           = result.config.to_dict(),
        series_x         = result.series_x.to_dict() if result.series_x else None,
        series_y         = result.series_y.to_dict() if result.series_y else None,
        features         = features_meta,
        surface          = result.surface.to_dict() if result.surface else None,
        warnings         = [w.to_dict() for w in all_warnings],
    )
