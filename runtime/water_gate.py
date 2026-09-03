# runtime/water_gate.py
"""
Gate: water-in-bilayer check after clean_water.

Checks clean_water_report.json (primary, written by run_clean_water.py ≥ v2)
or water_report.json (backward compat for workspaces built before the report
was renamed).  Blocks if water oxygens remain inside the bilayer hydrophobic
core after cleanup.
"""
from __future__ import annotations

import json
from pathlib import Path

from runtime.gate_runner import GateResult

_WARN_THRESHOLD = 5    # warn if 1–5 waters remain (might be at boundary)


def evaluate_water_gate(step_dir: Path) -> GateResult | None:
    """
    Returns None when no report is present.
    Blocks when remaining waters > _WARN_THRESHOLD (cleanup clearly failed).
    Warns for 1–_WARN_THRESHOLD waters (boundary ambiguity).

    Report priority:
      1. clean_water_report.json  (primary, contains full audit fields)
      2. water_report.json        (backward compat for older workspaces)
    """
    # Try primary report first
    primary = step_dir / "clean_water_report.json"
    if primary.exists():
        return _evaluate_clean_water_report(primary)

    # Backward compat: old workspaces only wrote water_report.json
    legacy = step_dir / "water_report.json"
    if not legacy.exists():
        return None
    return _evaluate_water_report(legacy)


def _evaluate_clean_water_report(path: Path) -> GateResult | None:
    """Read clean_water_report.json and produce a GateResult.

    Uses 'n_water_oxygens_remaining_in_core' (written by assembly_builder ≥ v3)
    which counts only OW atoms inside the hydrophobic core after deletion.
    Falls back to the legacy 'final_water_count' field with a warning, because
    that field is the total remaining waters in the system (not just core waters)
    and will produce false positives for any large solvated membrane system.
    """
    try:
        data = json.loads(path.read_text())
    except Exception:
        return None

    errors:   list[str] = []
    warnings: list[str] = []

    # Prefer the precise core-water count written by assembly_builder ≥ v3
    if "n_water_oxygens_remaining_in_core" in data:
        n_remain = data["n_water_oxygens_remaining_in_core"]
    else:
        # Legacy: final_water_count = total waters remaining in whole system — not core-only.
        # This will be a very large number for solvated membrane systems and is NOT
        # a reliable gate signal.  Treat as advisory only (never block on it).
        n_remain = 0
        legacy   = data.get("final_water_count", 0)
        if legacy > 0:
            warnings.append(
                f"clean_water_report.json uses legacy 'final_water_count' ({legacy} total "
                "waters in system — not waters in core); gate cannot evaluate cleanup accuracy. "
                "Recompile workspace to regenerate run_clean_water.py with core-water reporting."
            )

    if n_remain > _WARN_THRESHOLD:
        errors.append(
            f"{n_remain} water oxygen(s) remain in bilayer hydrophobic core after cleanup"
        )
    elif n_remain > 0:
        warnings.append(
            f"{n_remain} water oxygen(s) near bilayer core boundary after cleanup"
        )

    # Surface metadata fields for structured reporting
    core_z_min = data.get("core_z_min")
    core_z_max = data.get("core_z_max")
    cleanup_passed = data.get("cleanup_passed", n_remain == 0)

    return GateResult(
        passed     = cleanup_passed and n_remain == 0,
        blocked    = len(errors) > 0,
        confidence = 1.0,
        errors     = errors,
        warnings   = warnings,
    )


def _evaluate_water_report(path: Path) -> GateResult | None:
    """Read legacy water_report.json and produce a GateResult."""
    try:
        data = json.loads(path.read_text())
    except Exception:
        return None

    errors   = data.get("errors",   [])
    warnings = data.get("warnings", [])

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
