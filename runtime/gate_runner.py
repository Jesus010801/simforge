# runtime/gate_runner.py
"""
Shared GateResult type and unified gate dispatcher.

All new gate modules return a GateResult.
The executor calls run_gate() instead of hard-coding each gate type.
Legacy gates (orientation_report, box_match_report) are wrapped transparently.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class GateResult:
    """Common gate outcome used by all gate modules."""
    passed:     bool
    blocked:    bool        # True → downstream steps FAILED/BLOCKED
    confidence: float       # 0.0–1.0; higher = more certain
    errors:     list[str] = field(default_factory=list)
    warnings:   list[str] = field(default_factory=list)

    @property
    def status_str(self) -> str:
        if self.blocked:   return "BLOCKED"
        if self.warnings:  return "PASS (advisory)"
        return "PASS"


# Human-readable names shown in executor output
GATE_LABELS: dict[str, str] = {
    "orientation_report":   "Orientation",
    "box_match_report":     "Box-match",
    "overlap_report":       "Overlap",
    "topology_consistency": "Topology consistency",
    "apl_report":           "APL convergence",
    "water_report":         "Water cleanup",
}


def run_gate(gate_type: str, step_dir: Path) -> GateResult | None:
    """
    Dispatch to the appropriate gate evaluator.
    Returns None when the report file is absent or unreadable.
    """
    if gate_type == "orientation_report":
        from runtime.orientation_gate import evaluate_orientation_gate
        return _adapt(evaluate_orientation_gate(step_dir))
    elif gate_type == "box_match_report":
        from runtime.box_match_gate import evaluate_box_match_gate
        return _adapt(evaluate_box_match_gate(step_dir))
    elif gate_type == "overlap_report":
        from runtime.overlap_gate import evaluate_overlap_gate
        return evaluate_overlap_gate(step_dir)
    elif gate_type == "topology_consistency":
        from runtime.topology_gate import evaluate_topology_gate
        return evaluate_topology_gate(step_dir)
    elif gate_type == "apl_report":
        from runtime.apl_gate import evaluate_apl_gate
        return evaluate_apl_gate(step_dir)
    elif gate_type == "water_report":
        from runtime.water_gate import evaluate_water_gate
        return evaluate_water_gate(step_dir)
    return None


def _adapt(result) -> GateResult | None:
    """Wrap legacy gate result (orientation/box_match) to canonical GateResult."""
    if result is None:
        return None
    return GateResult(
        passed=result.passed,
        blocked=result.blocked,
        confidence=result.confidence,
        errors=list(result.errors),
        warnings=list(result.warnings),
    )
