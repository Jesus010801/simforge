# runtime/overlap_gate.py
"""
Gate: protein-lipid overlap check after embed_in_bilayer.

Reads overlap_report.json written by embed_in_bilayer/run_embed.py.
Blocks if any protein–lipid clashes are detected.
"""
from __future__ import annotations

import json
from pathlib import Path

from runtime.gate_runner import GateResult


def evaluate_overlap_gate(step_dir: Path) -> GateResult | None:
    """
    Returns None when the report is absent (gate not generated = no block).
    Blocks when n_clashes > 0 or protein/lipid atoms not found.
    """
    report_path = step_dir / "overlap_report.json"
    if not report_path.exists():
        return None
    try:
        data = json.loads(report_path.read_text())
    except Exception:
        return None

    errors   = data.get("errors",   [])
    warnings = data.get("warnings", [])

    return GateResult(
        passed     = data.get("passed", False),
        blocked    = len(errors) > 0,
        confidence = data.get("confidence", 1.0),
        errors     = errors,
        warnings   = warnings,
    )
