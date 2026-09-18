"""Water-cleanup gate based on confirmed region defects and steric clashes."""
from __future__ import annotations
import json
from pathlib import Path
from runtime.gate_runner import GateResult


def evaluate_water_gate(step_dir: Path) -> GateResult | None:
    primary=Path(step_dir)/"clean_water_report.json"
    legacy=Path(step_dir)/"water_report.json"
    path=primary if primary.exists() else legacy
    if not path.exists(): return None
    try: data=json.loads(path.read_text())
    except Exception: return None
    warnings=list(data.get("warnings",[]));errors=[]
    has_atlas=("supported_membrane_defects_remaining" in data or "confirmed_hard_clashes_remaining" in data
               or "schema_version" in data)
    if has_atlas:
        defects=int(data.get("supported_membrane_defects_remaining",0))
        clashes=int(data.get("confirmed_hard_clashes_remaining",0))
        if defects: errors.append(f"{defects} supported membrane-defect waters remain")
        if clashes: errors.append(f"{clashes} confirmed water clashes remain")
        ncore=int(data.get("n_water_oxygens_remaining_in_core",0))
        if ncore: warnings.append(f"{ncore} waters remain in classified membrane-interior regions")
        passed=not errors and bool(data.get("cleanup_passed",True))
    else:
        # Historical reports only measure a Z slab. This cannot establish that
        # retained water is invalid, so it remains an advisory signal.
        has_core_count=("n_water_oxygens_remaining_in_core" in data or "n_waters_remaining" in data)
        ncore=int(data.get("n_water_oxygens_remaining_in_core",data.get("n_waters_remaining",0)))
        if ncore: warnings.append(f"{ncore} water(s) remain according to a legacy slab/count report; review region classification")
        if not has_core_count:
            warnings.append("Legacy water report has no region-atlas counts; membrane-water status is unresolved")
        for item in data.get("errors",[]):
            low=str(item).lower()
            if any(x in low for x in ("core","bilayer","membrane","slab")):
                warnings.append(str(item))
            else: errors.append(str(item))
        passed=not errors
    return GateResult(passed=passed,blocked=bool(errors),confidence=float(data.get("confidence",1.0)),
                      errors=errors,warnings=warnings)
