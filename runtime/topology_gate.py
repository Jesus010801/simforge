# runtime/topology_gate.py
"""
Gate: topology consistency check after generate_topology.

Reads topology_consistency_report.json written by generate_topology/run_topology.py.
Blocks if topol.top is missing or posre.itp not included.
Warns if strong_posre.itp was expected but not injected.
"""
from __future__ import annotations

import json
from pathlib import Path

from runtime.gate_runner import GateResult


def evaluate_topology_gate(step_dir: Path) -> GateResult | None:
    """
    Returns None when the report is absent.
    Blocks on hard errors (missing topol.top, missing posre.itp).
    Warns when strong_posre.itp was not injected (e.g. source file not found).
    """
    report_path = step_dir / "topology_consistency_report.json"
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
