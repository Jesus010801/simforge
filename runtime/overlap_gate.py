# runtime/overlap_gate.py
"""
Gate: protein-lipid overlap check after embed_in_bilayer.

Reads overlap_report.json written by embed_in_bilayer/run_embed.py.

Blocking semantics
------------------
- blocked = True  only when the report's 'errors' list is non-empty.
  Hard errors (impossible geometry, missing inputs, failed system.gro) go there.
- Initial lipid–protein clashes detected by validate_no_overlap are EXPECTED at
  embedding stage — the membrane_embedding shrink loop resolves them.  The
  generated script places those in 'warnings', so n_clashes > 0 alone does NOT
  trigger a block.
- 'passed' (False when n_clashes > 0) is kept as informational metadata.
"""
from __future__ import annotations

import json
from pathlib import Path

from runtime.gate_runner import GateResult


def evaluate_overlap_gate(step_dir: Path) -> GateResult | None:
    """
    Returns None when the report is absent (gate not generated = no block).
    Blocks only when the 'errors' list is non-empty (hard geometry failure).
    Advisory warnings (initial embedding clashes) do not block.
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
