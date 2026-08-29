"""Deterministic semantic index (.ndx) generation.

Turns resolved :class:`MolecularComponent` objects into named GROMACS index
groups (``Receptor``, ``Receptor_Backbone``, ``Peptide``, ``Complex``,
``Membrane``, ``Water``, ``Ions`` …) so the analysis layer never references
numeric GROMACS group IDs.

Selections are built with ``gmx select`` from the reference structure:

* PDB with chain IDs  -> ``chain A B`` style selections;
* GRO / no chains     -> contiguous ``resid LO to HI`` ranges derived from the
  parsed residue list; non-contiguous components emit a warning and the group
  is omitted (analyses needing it then become ``review_required``).

The user's own ``.ndx`` is never overwritten — output goes to the analysis
workspace.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

from analysis.campaign.gmx import gmx_available, run_gmx
from analysis.campaign.models import (
    CampaignWarning, ClassificationState, ComponentType, MolecularComponent,
    SemanticIndex, SemanticIndexGroup, Severity,
)
from analysis.campaign.structure.parsers import StructureModel, parse_structure


def _contiguous_ranges(nums: list[int]) -> list[tuple[int, int]]:
    if not nums:
        return []
    s = sorted(set(nums))
    ranges: list[tuple[int, int]] = []
    lo = prev = s[0]
    for n in s[1:]:
        if n == prev + 1:
            prev = n
            continue
        ranges.append((lo, prev))
        lo = prev = n
    ranges.append((lo, prev))
    return ranges


def _component_selection(
    comp: MolecularComponent, structure: StructureModel,
) -> tuple[Optional[str], str, list[str]]:
    """Return (gmx-select-expression, logic-description, warnings)."""
    warns: list[str] = []
    if structure.has_chain_ids and comp.chain_ids and all(
        c not in ("_", " ") for c in comp.chain_ids
    ):
        chains = " ".join(sorted(comp.chain_ids))
        return f'chain {chains}', f"chain(s) {chains}", warns

    # residue-range fallback
    resids: list[int] = []
    for chain in structure.chains:
        if comp.chain_ids and chain.chain_id not in comp.chain_ids:
            continue
        for resnum, rn in chain.residues:
            from analysis.campaign.structure.residue_classes import classify_resname
            if comp.component_type in (
                ComponentType.RECEPTOR, ComponentType.PEPTIDE,
                ComponentType.PROTEIN_PARTNER, ComponentType.COMPLEX,
                ComponentType.NUCLEIC_ACID,
            ) and classify_resname(rn) not in ("amino_acid", "nucleotide"):
                continue
            resids.append(resnum)
    if not resids:
        warns.append(f"could not determine residues for component '{comp.label}'")
        return None, "", warns
    ranges = _contiguous_ranges(resids)
    if len(ranges) > 8:
        warns.append(
            f"component '{comp.label}' residues are split into {len(ranges)} "
            f"non-contiguous blocks; residue-range selection is unreliable"
        )
    expr = " or ".join(f"(resid {lo} to {hi})" for lo, hi in ranges)
    logic = f"resid ranges {ranges}"
    return expr, logic, warns


_STANDARD_BACKBONE = {
    ComponentType.RECEPTOR: "Receptor",
    ComponentType.PEPTIDE: "Peptide",
    ComponentType.PROTEIN_PARTNER: "ProteinPartner",
    ComponentType.COMPLEX: "Complex",
}


def build_semantic_index(
    *,
    structure_path: str | Path,
    components: list[MolecularComponent],
    out_ndx: str | Path,
    gmx: str = "gmx",
    existing_user_index: Optional[str | Path] = None,
) -> SemanticIndex:
    structure_path = Path(structure_path)
    idx = SemanticIndex(source_structure=str(structure_path.resolve()))

    if existing_user_index and Path(existing_user_index).is_file():
        idx.warnings.append(
            f"user index {existing_user_index} left untouched; semantic groups "
            f"written to a separate file"
        )

    try:
        structure = parse_structure(structure_path)
    except ValueError as exc:
        idx.warnings.append(str(exc))
        return idx

    if not gmx_available(gmx):
        idx.warnings.append(
            f"'{gmx}' not available; semantic index not generated. Analyses that "
            f"need named groups will be marked review_required."
        )
        return idx

    # Build the list of (name, selection-expression, logic) to request.
    requests: list[tuple[str, str, str, list[str]]] = []
    for comp in components:
        if comp.component_type in (
            ComponentType.WATER, ComponentType.IONS, ComponentType.MEMBRANE,
        ):
            name = {"water": "Water", "ions": "Ions", "membrane": "Membrane"}[comp.label] \
                if comp.label in ("water", "ions", "membrane") else comp.component_type.capitalize()
            resn = " ".join(comp.resnames) if comp.resnames else ""
            if resn:
                requests.append((name, f"resname {resn}", f"resname {resn}", []))
            continue
        if comp.component_type not in _STANDARD_BACKBONE:
            continue
        if comp.classification_state in (
            ClassificationState.UNKNOWN, ClassificationState.REVIEW_REQUIRED,
        ):
            idx.warnings.append(
                f"component '{comp.label}' is {comp.classification_state}; "
                f"its semantic group was not generated"
            )
            continue
        base = _STANDARD_BACKBONE[comp.component_type]
        expr, logic, w = _component_selection(comp, structure)
        idx.warnings.extend(w)
        if expr is None:
            continue
        prot = "" if comp.component_type == ComponentType.COMPLEX else " and "
        # constrain protein-type groups to protein atoms for safety
        if comp.component_type in (ComponentType.RECEPTOR, ComponentType.PEPTIDE,
                                   ComponentType.PROTEIN_PARTNER):
            full = f'(group "Protein" and ({expr}))'
        else:
            full = f'({expr})'
        requests.append((base, full, logic, []))
        requests.append((f"{base}_Backbone", f'{full} and name CA C N',
                         f"{logic}; protein backbone atoms (N, CA, C — GROMACS 'Backbone' convention)", []))
        requests.append((f"{base}_CA", f'{full} and name CA', f"{logic}; C-alpha atoms", []))

    if not requests:
        idx.warnings.append("no semantic groups could be constructed")
        return idx

    select_arg = "; ".join(f'"{name}" {expr}' for name, expr, _logic, _w in requests)
    out_ndx = Path(out_ndx)
    out_ndx.parent.mkdir(parents=True, exist_ok=True)
    res = run_gmx(
        ["select", "-s", str(structure_path), "-on", str(out_ndx), "-select", select_arg],
        gmx=gmx, timeout=300,
    )
    idx.gmx_commands.append(res.argv)
    if not res.ok:
        idx.warnings.append(
            f"`gmx select` failed (rc={res.returncode}): {res.stderr.strip()[-400:]}"
        )
        return idx

    idx.path = str(out_ndx.resolve())
    _populate_group_sizes(idx, out_ndx, requests, structure)
    return idx


def _populate_group_sizes(
    idx: SemanticIndex, ndx_path: Path,
    requests: list[tuple[str, str, str, list[str]]], structure: StructureModel,
) -> None:
    logic_by_name = {name: logic for name, _e, logic, _w in requests}
    counts: dict[str, int] = {}
    current: Optional[str] = None
    for line in ndx_path.read_text().splitlines():
        line = line.strip()
        if line.startswith("[") and line.endswith("]"):
            current = line.strip("[] ").strip()
            counts[current] = 0
        elif current is not None and line:
            counts[current] += len(line.split())
    for name, n in counts.items():
        idx.groups.append(SemanticIndexGroup(
            name=name, n_atoms=n,
            source="gmx select from reference structure",
            selection_logic=logic_by_name.get(name, ""),
        ))
    for name, n in counts.items():
        if n == 0:
            idx.warnings.append(f"semantic group '{name}' is empty")
