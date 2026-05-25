# runtime/water_gate.py
"""
Gate: water-in-bilayer check after clean_water.

Reads water_report.json written by clean_water/run_clean_water.py.
Blocks if water oxygens remain inside the bilayer hydrophobic core after cleanup.
"""
from __future__ import annotations

import json
from pathlib import Path

from runtime.gate_runner import GateResult

_WARN_THRESHOLD = 5    # warn if 1–5 waters remain (might be at boundary)


def evaluate_water_gate(step_dir: Path) -> GateResult | None:
    """
    Returns None when the report is absent.
    Blocks when n_waters_remaining > _WARN_THRESHOLD (cleanup clearly failed).
    Warns for 1–_WARN_THRESHOLD waters (boundary ambiguity).
    """
    report_path = step_dir / "water_report.json"
    if not report_path.exists():
        return None
    try:
        data = json.loads(report_path.read_text())
    except Exception:
        return None

    errors   = data.get("errors",   [])
    warnings = data.get("warnings", [])

    # Augment with structured threshold logic if report has raw counts
    n_remain = data.get("n_waters_remaining", 0)
    if not errors and not warnings and n_remain > 0:
        msg = f"{n_remain} water oxygen(s) remain inside bilayer core after cleanup"
        if n_remain > _WARN_THRESHOLD:
            errors = [msg]
        else:
            warnings = [msg]

    return GateResult(
        passed     = data.get("passed", False),
        blocked    = len(errors) > 0,
        confidence = data.get("confidence", 1.0),
        errors     = errors,
        warnings   = warnings,
    )
