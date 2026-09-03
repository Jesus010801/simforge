"""
High-level system build: ``SystemSpec`` -> assembled GROMACS system.

This is the orchestration layer behind ``simforge build``. It:

  1. maps a resolved ``SystemSpec`` onto the already-validated assembly engine
     ``ligand.integrate.assemble_system_multi`` (raw-PDB or preparameterized
     protein mode, N parameterized components, atomtype namespace isolation),
  2. runs a deterministic *preflight* — every atomtype a component references
     must be resolvable (its own [ atomtypes ], another component's, or the
     base force field) — and fails BEFORE assembly if not (Phase 6),
  3. writes ``provenance.json`` and ``system_build_report.json``,
  4. runs ``core.system_validate.validate_system_dir`` on the output.

No scientific assembly logic is duplicated here.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

from core.provenance import build_provenance, gromacs_data_prefix
from core.system_spec import SystemSpec
from core.system_validate import validate_system_dir
from ligand.integrate import ParameterizedMolecule, assemble_system_multi
from utils.itp_parser import parse_itp


class SystemBuildError(RuntimeError):
    pass


@dataclass
class BuildOutcome:
    success: bool
    out_dir: Path
    spec: SystemSpec
    assembly_result: Any = None
    provenance: dict[str, Any] = field(default_factory=dict)
    validation: dict[str, Any] = field(default_factory=dict)
    preflight_errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)


# ── Phase 6: deterministic atomtype preflight ────────────────────────────────

def _atomtypes_defined_in(itp_text: str) -> set[str]:
    names: set[str] = set()
    in_block = False
    for line in itp_text.splitlines():
        s = line.split(";")[0].strip()
        if not s:
            continue
        if s.startswith("["):
            in_block = s.strip("[] ").lower() == "atomtypes"
            continue
        if in_block:
            parts = s.split()
            if parts:
                names.add(parts[0])
    return names


def _ff_atomtypes(forcefield: str) -> Optional[set[str]]:
    top = gromacs_data_prefix()
    if top is None:
        return None
    ff_dir = Path(top) / f"{forcefield}.ff"
    if not ff_dir.is_dir():
        return None
    names: set[str] = set()
    atp = ff_dir / "atomtypes.atp"
    if atp.exists():
        for line in atp.read_text().splitlines():
            s = line.split(";")[0].strip()
            if s:
                names.add(s.split()[0])
    nb = ff_dir / "ffnonbonded.itp"
    if nb.exists():
        names |= _atomtypes_defined_in(nb.read_text())
    return names or None


def preflight_component_atomtypes(
    components: list[ParameterizedMolecule], forcefield: str
) -> list[str]:
    """Return a list of actionable error strings (empty == all resolvable).

    A component that references atomtypes with no [ atomtypes ] section at all,
    or references types absent from every component's [ atomtypes ] AND the
    base force field, is reported here — before assembly — with the component,
    the missing types, the source ITP and a remediation hint.
    """
    ff_types = _ff_atomtypes(forcefield)  # None -> GROMACS FF data not inspectable
    per_component_defined: dict[str, set[str]] = {}
    per_component_referenced: dict[str, set[str]] = {}
    per_component_has_block: dict[str, bool] = {}

    for c in components:
        try:
            text = Path(c.topology_path).read_text()
        except Exception as exc:
            return [f"[{c.component_id}] cannot read topology {c.topology_path}: {exc}"]
        defined = _atomtypes_defined_in(text)
        per_component_defined[c.component_id] = defined
        per_component_has_block[c.component_id] = "[ atomtypes ]" in text or "[atomtypes]" in text
        try:
            parsed = parse_itp(c.topology_path)
            per_component_referenced[c.component_id] = {a.atom_type for a in parsed.atoms}
        except Exception as exc:
            return [f"[{c.component_id}] cannot parse topology {c.topology_path}: {exc}"]

    all_component_defined: set[str] = set().union(*per_component_defined.values()) if per_component_defined else set()

    errors: list[str] = []
    for c in components:
        referenced = per_component_referenced.get(c.component_id, set())
        resolvable = set(all_component_defined)
        if ff_types is not None:
            resolvable |= ff_types
        missing = sorted(t for t in referenced if t not in resolvable)
        if not missing:
            continue

        if ff_types is None and not per_component_has_block[c.component_id]:
            # Can't consult FF; but a component with zero [ atomtypes ] that
            # references non-trivial types is the exact known failure mode.
            errors.append(
                f"[{c.component_id}] {Path(c.topology_path).name} references "
                f"atomtype(s) {missing} but contains NO [ atomtypes ] section. "
                f"LigParGen components must keep their original [ atomtypes ] "
                f"block so SimForge can namespace collisions. Re-export the "
                f"component with its [ atomtypes ] definitions included."
            )
        elif ff_types is not None:
            hint = (
                "contains no [ atomtypes ] section — re-export it from LigParGen "
                "with atomtypes included"
                if not per_component_has_block[c.component_id]
                else f"its [ atomtypes ] block defines {sorted(per_component_defined[c.component_id])}"
            )
            errors.append(
                f"[{c.component_id}] {Path(c.topology_path).name} references "
                f"atomtype(s) {missing} that are defined neither by any component "
                f"nor by force field '{forcefield}'. Component {hint}."
            )
    return errors


# ── mapping spec -> assembly kwargs ─────────────────────────────────────────

def _components_from_spec(spec: SystemSpec) -> list[ParameterizedMolecule]:
    return [
        ParameterizedMolecule(
            component_id=c.id,
            topology_path=spec.component_topology_paths[c.id],
            coordinate_path=spec.component_coordinate_paths[c.id],
            role=c.role,
        )
        for c in spec.components
    ]


def build_system(
    spec: SystemSpec,
    out_dir: str | Path,
    *,
    run_grompp: bool = True,
    command: Optional[list[str]] = None,
    skip_preflight: bool = False,
) -> BuildOutcome:
    out_dir = Path(out_dir).resolve()
    components = _components_from_spec(spec)
    outcome = BuildOutcome(success=False, out_dir=out_dir, spec=spec)
    outcome.warnings.extend(spec.warnings)

    # ── Phase 6 preflight ──────────────────────────────────────────────────
    if not skip_preflight:
        pf = preflight_component_atomtypes(components, spec.forcefield)
        if pf:
            outcome.preflight_errors = pf
            outcome.errors = pf
            return outcome

    # ── assembly (validated engine) ───────────────────────────────────────
    kwargs: dict[str, Any] = dict(
        components=components,
        out_dir=out_dir,
        forcefield=spec.forcefield,
        water_model=spec.water_model,
        run_grompp=run_grompp,
    )
    if spec.protein_mode == "raw_pdb":
        kwargs["protein_pdb"] = spec.protein_structure_path
    else:
        kwargs["protein_gro"] = spec.protein_structure_path
        kwargs["protein_topology_itps"] = spec.protein_topology_paths
        kwargs["protein_posre_itps"] = spec.protein_restraint_paths or None

    result = assemble_system_multi(**kwargs)
    outcome.assembly_result = result

    if result.error:
        outcome.errors = [result.error]
        return outcome
    outcome.warnings.extend(result.warnings or [])
    outcome.errors.extend(result.errors or [])

    # ── provenance ────────────────────────────────────────────────────────
    inputs: dict[str, Any] = {"protein_structure": spec.protein_structure_path}
    for i, p in enumerate(spec.protein_topology_paths):
        inputs[f"protein_topology_{i}"] = p
    for i, p in enumerate(spec.protein_restraint_paths):
        inputs[f"protein_restraint_{i}"] = p
    for c in spec.components:
        inputs[f"component_{c.id}_topology"] = spec.component_topology_paths[c.id]
        inputs[f"component_{c.id}_coordinates"] = spec.component_coordinate_paths[c.id]

    outputs = {
        "topol_top": result.topol_top,
        "complex_gro": result.complex_gro,
        "protein_gro": result.protein_gro,
        "assembly_report": result.report_path,
    }

    assembly_meta = {
        "forcefield": spec.forcefield,
        "water_model": spec.water_model,
        "protein_mode": spec.protein_mode,
        "pdb2gmx_used": spec.protein_mode == "raw_pdb",
        "pdb2gmx_ignh_fallback": bool(result.pdb2gmx_fallback_used),
        "discovered_protein_topology": [p.name for p in spec.protein_topology_paths],
        "discovery_notes": spec.discovery_notes,
        "protein_molecules": result.protein_molecules,
        "component_molecule_names": {c["id"]: c["molecule_name"] for c in result.components},
        "atomtype_renames": {
            c["id"]: c.get("atomtype_renames", {}) for c in result.components
            if c.get("atomtype_renames")
        },
        "total_atom_count": result.total_atom_count,
        "warnings": outcome.warnings,
    }

    # ── validation ────────────────────────────────────────────────────────
    vreport = validate_system_dir(
        out_dir, forcefield=spec.forcefield, run_grompp=run_grompp
    )
    outcome.validation = vreport.to_dict()

    provenance = build_provenance(
        spec=spec,
        spec_path=spec._spec_path,
        command=command,
        inputs=inputs,
        outputs=outputs,
        assembly=assembly_meta,
        validation={"readiness": vreport.readiness, "ok": vreport.ok,
                    "checks": [c.to_dict() for c in vreport.checks]},
    )
    outcome.provenance = provenance

    (out_dir / "provenance.json").write_text(json.dumps(provenance, indent=2, default=str))
    (out_dir / "validation_report.json").write_text(
        json.dumps(vreport.to_dict(), indent=2, default=str)
    )
    (out_dir / "system_build_report.json").write_text(json.dumps({
        "spec": spec.normalized_dict(),
        "protein_mode": spec.protein_mode,
        "discovery_notes": spec.discovery_notes,
        "components": result.components,
        "protein_molecules": result.protein_molecules,
        "total_atom_count": result.total_atom_count,
        "warnings": outcome.warnings,
        "errors": outcome.errors,
        "validation_readiness": vreport.readiness,
    }, indent=2, default=str))

    outcome.success = bool(result.success) and not outcome.errors
    return outcome
