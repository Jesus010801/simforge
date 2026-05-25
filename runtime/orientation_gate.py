from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class GateResult:
    """Outcome of reading orientation_report.json and applying gate rules."""
    passed: bool
    blocked: bool        # True → downstream steps must be FAILED/BLOCKED
    confidence: float
    errors: list[str]   = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def status_str(self) -> str:
        if self.blocked:
            return "BLOCKED"
        if self.warnings:
            return "PASS (advisory)"
        return "PASS"


def evaluate_orientation_gate(step_dir: Path) -> GateResult | None:
    """
    Read orientation_report.json from step_dir and apply gate rules:
      - errors > 0  → blocked=True  (downstream FAILED)
      - warnings > 0 → blocked=False, caller shows advisory
      - else         → clean pass

    Returns None when the report file is absent or unreadable.
    """
    report_path = step_dir / "orientation_report.json"
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
        confidence = data.get("confidence", 0.0),
        errors     = errors,
        warnings   = warnings,
    )


def read_orientation_report(workspace_path: Path) -> dict | None:
    """
    Scan workspace steps for the first orientation_report.json.
    Returns the raw dict (for CLI display) or None.
    """
    steps_dir = workspace_path / "steps"
    if not steps_dir.exists():
        return None
    for step_dir in sorted(steps_dir.iterdir()):
        if not step_dir.is_dir():
            continue
        p = step_dir / "orientation_report.json"
        if p.exists():
            try:
                return json.loads(p.read_text())
            except Exception:
                return None
    return None
