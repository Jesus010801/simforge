# runtime/tm_aware_gate.py
"""
Gate: TM-aware embedding quality check.
Reads tm_aware_backend_report.json to block or pass downstream stages.
"""
from __future__ import annotations

import json
from pathlib import Path

from runtime.gate_runner import GateResult


def evaluate_tm_aware_gate(step_dir: Path) -> GateResult | None:
    """
    Returns None when the backend report is missing.
    Otherwise, returns a GateResult based on the quality metrics and block status.
    """
    report_path = step_dir / "tm_aware_backend_report.json"
    if not report_path.exists():
        return None
    try:
        data = json.loads(report_path.read_text())
    except Exception:
        return None

    quality_passed = data.get("quality_passed", False)
    warnings = data.get("warnings", [])
    blocking_reason = data.get("blocking_reason")

    errors = []
    if not quality_passed:
        msg = blocking_reason or "TM-aware quality gates failed."
        errors.append(msg)

    return GateResult(
        passed=quality_passed,
        blocked=not quality_passed,
        confidence=1.0,
        errors=errors,
        warnings=warnings,
    )
