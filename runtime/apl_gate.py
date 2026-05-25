# runtime/apl_gate.py
"""
Gate: APL convergence check after membrane_embedding shrink loop.

Reads shrink_telemetry.json (written by the shrink loop bash script).
The shrink loop already exits 1 on non-convergence, so the gate is
informational — it emits a warning if APL is only barely within tolerance.
"""
from __future__ import annotations

import json
from pathlib import Path

from runtime.gate_runner import GateResult

_NEAR_TOLERANCE_FRACTION = 0.25   # warn if APL is within 25% of the tolerance band


def evaluate_apl_gate(step_dir: Path) -> GateResult | None:
    """
    Returns None when telemetry is absent (shrink loop did not run).
    Never blocks (script exits 1 on non-convergence, executor handles that).
    Warns when APL is very close to the upper tolerance limit.
    """
    telemetry_path = step_dir / "shrink_telemetry.json"
    if not telemetry_path.exists():
        return None
    try:
        data = json.loads(telemetry_path.read_text())
    except Exception:
        return None

    converged    = data.get("converged", False)
    final_apl    = data.get("final_apl_ang2")
    n_iterations = data.get("n_iterations", 0)

    # Read target + tolerance from metadata.json if available
    target_apl  = None
    tolerance   = None
    meta_path   = step_dir / "metadata.json"
    if meta_path.exists():
        try:
            meta = json.loads(meta_path.read_text())
            params = meta.get("params", {})
            target_apl = params.get("apl_target_ang2")
            tolerance  = params.get("apl_tolerance_ang2")
        except Exception:
            pass

    if not converged:
        # Script already failed, but if somehow we reach the gate, block.
        msg = f"Shrink loop did not converge (final APL={final_apl} Å², {n_iterations} iterations)"
        return GateResult(
            passed=False, blocked=True, confidence=1.0,
            errors=[msg],
        )

    warnings: list[str] = []
    if final_apl is not None and target_apl is not None and tolerance is not None:
        cutoff = target_apl + tolerance
        gap    = cutoff - final_apl
        if gap < tolerance * _NEAR_TOLERANCE_FRACTION:
            warnings.append(
                f"APL={final_apl:.1f} Å² is very close to tolerance limit "
                f"({cutoff:.1f} Å²) — consider additional deflation iterations"
            )

    apl_str = f"{final_apl:.1f} Å²" if final_apl is not None else "unknown"
    msg = f"Shrink loop converged after {n_iterations} iterations (APL={apl_str})"

    return GateResult(
        passed=True, blocked=False, confidence=1.0,
        warnings=warnings,
    )
