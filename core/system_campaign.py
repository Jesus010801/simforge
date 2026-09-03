"""
Declarative campaigns (Phase 13) — orchestrate several unified system
specifications that share a protein (and optionally shared components) without
repeating the protein declaration for every system.

    campaign:
      name: hmgcoa_series
      shared:
        protein:
          structure: protein_only.gro
          topology: auto
          restraints: auto
        forcefield: oplsaa
        water_model: spce
        components:                       # shared across every system
          - {id: COA, topology: COA.itp, coordinates: COA.gro, role: cofactor}
      systems:
        - name: A1
          components:
            - {id: A1, topology: A1.itp, coordinates: A1.gro, role: competitive_ligand}
        - name: A3
          components:
            - {id: A3, topology: A3.itp, coordinates: A3.gro, role: competitive_ligand}

Each entry under ``systems`` is expanded into a full ``SystemSpec`` by deep-merging
``shared`` with the entry (the entry wins on scalars; component lists are
concatenated, shared components first). A campaign is *only* an orchestrator —
it calls ``core.system_build.build_system`` per system, never re-implements it.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import yaml

from core.system_spec import SystemSpec, SystemSpecError
from core.system_build import BuildOutcome, build_system


@dataclass
class CampaignResult:
    name: str
    outcomes: dict[str, BuildOutcome] = field(default_factory=dict)
    errors: dict[str, str] = field(default_factory=dict)

    @property
    def total(self) -> int:
        return len(self.outcomes) + len(self.errors)

    @property
    def succeeded(self) -> int:
        return sum(1 for o in self.outcomes.values() if o.success)

    @property
    def success(self) -> bool:
        return not self.errors and all(o.success for o in self.outcomes.values())


def _unwrap_campaign(doc: Any) -> dict:
    if not isinstance(doc, dict):
        raise SystemSpecError("Campaign spec must be a YAML mapping.")
    if "campaign" in doc and isinstance(doc["campaign"], dict):
        return doc["campaign"]
    return doc


def expand_campaign(path: str | Path) -> tuple[str, list[tuple[str, SystemSpec]], Path]:
    """Parse a campaign YAML into (campaign_name, [(system_name, SystemSpec), ...], base_dir).

    Every returned ``SystemSpec`` is already path-resolved and validated.
    """
    p = Path(path).resolve()
    if not p.exists():
        raise SystemSpecError(f"campaign file not found: {p}")
    try:
        doc = yaml.safe_load(p.read_text())
    except yaml.YAMLError as exc:
        raise SystemSpecError(f"invalid YAML in {p}: {exc}")

    body = _unwrap_campaign(doc)
    name = body.get("name") or p.stem
    shared = body.get("shared") or {}
    systems = body.get("systems") or []
    if not isinstance(systems, list) or not systems:
        raise SystemSpecError("campaign.systems must be a non-empty list.")
    if "protein" not in shared:
        raise SystemSpecError("campaign.shared.protein is required.")

    shared_components = list(shared.get("components") or [])

    specs: list[tuple[str, SystemSpec]] = []
    errors: list[str] = []
    seen_names: set[str] = set()
    for i, entry in enumerate(systems):
        entry = entry or {}
        sys_name = entry.get("name") or f"system_{i + 1}"
        if sys_name in seen_names:
            errors.append(f"duplicate system name '{sys_name}' in campaign.")
        seen_names.add(sys_name)

        merged: dict[str, Any] = {
            "name": f"{name}__{sys_name}",
            "protein": copy.deepcopy(entry.get("protein") or shared["protein"]),
            "forcefield": entry.get("forcefield", shared.get("forcefield", "oplsaa")),
            "water_model": entry.get("water_model", shared.get("water_model", "spce")),
            "components": shared_components + list(entry.get("components") or []),
        }
        if "box" in entry or "box" in shared:
            merged["box"] = entry.get("box", shared.get("box"))

        try:
            spec = SystemSpec.model_validate(merged).resolve(base_dir=p.parent, spec_path=p)
            specs.append((sys_name, spec))
        except SystemSpecError as exc:
            errors.append(f"[{sys_name}] " + "; ".join(exc.messages))

    if errors:
        raise SystemSpecError(errors)
    return name, specs, p.parent


def run_campaign(
    path: str | Path,
    out_dir: str | Path,
    *,
    run_grompp: bool = True,
    command: Optional[list[str]] = None,
) -> CampaignResult:
    name, specs, _ = expand_campaign(path)
    out_dir = Path(out_dir).resolve()
    result = CampaignResult(name=name)

    for sys_name, spec in specs:
        sys_out = out_dir / sys_name
        try:
            outcome = build_system(
                spec, sys_out, run_grompp=run_grompp, command=command,
            )
            result.outcomes[sys_name] = outcome
        except Exception as exc:  # noqa: BLE001 - report, don't abort the campaign
            result.errors[sys_name] = str(exc)

    (out_dir / "campaign_summary.yaml").parent.mkdir(parents=True, exist_ok=True)
    (out_dir / "campaign_summary.yaml").write_text(yaml.safe_dump({
        "campaign": name,
        "systems": {
            n: {
                "success": o.success,
                "total_atoms": getattr(o.assembly_result, "total_atom_count", None),
                "validation": o.validation.get("readiness"),
                "out_dir": str(o.out_dir),
            }
            for n, o in result.outcomes.items()
        },
        "errors": result.errors,
    }, sort_keys=False))
    return result
