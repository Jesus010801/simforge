"""Does a calculation need a finite-assembly (clustered) view?  (Phase 13.5)

An observable declares *which semantic groups* must be globally coherent
(``assembly_groups``: its selection for Rg / COM / RMSD, the union A ∪ B for a
COM distance).  Whether that needs reconstruction is a property of the
topology, decided here generically:

* the union spans ONE topological molecule → per-molecule ``whole`` already
  gives its global geometry: nothing changes (identical numerics, identity);
* it spans SEVERAL molecules → every molecule must share one periodic image:
  ``cluster_groups`` is requested and judged by the policy;
* the selection, or any molecule it spans, belongs to a periodically extended
  component (membrane, solvent, ions) → a finite assembly is undefined:
  blocked (checked on atoms, independently of the molecule partition).

The decision does not depend on whether a wrap was observed — the metric's
geometry requirement is primary; diagnostics only explain the risk.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from analysis.campaign.models import ComponentType, GeometryRequirement

EXTENDED_GROUPS = {"Membrane": ComponentType.MEMBRANE, "Water": ComponentType.WATER,
                   "Ions": ComponentType.IONS}


@dataclass
class AssemblyDecision:
    groups: tuple[str, ...]
    geometry: str                            # GeometryRequirement.*
    needed: bool = False
    blocked: str = ""
    evidence: dict = field(default_factory=dict)


def decide_assembly(system, semantic_index, groups, gmx: str = "gmx") -> Optional[AssemblyDecision]:
    """None when nothing was requested; otherwise the decision with evidence."""
    groups = tuple(dict.fromkeys(g for g in groups if g))
    if not groups:
        return None
    from analysis.campaign.structure.index_groups import resolve_index_group
    from analysis.campaign.structure.molecules import MoleculePartitionError, molecule_partition
    if semantic_index is None or not semantic_index.path:
        return AssemblyDecision(groups, GeometryRequirement.FINITE_ASSEMBLY,
                                blocked="no semantic index to establish the selection's atoms")
    atoms: set[int] = set()
    for g in groups:
        grp = resolve_index_group(semantic_index.path, g)
        if grp is None:
            return AssemblyDecision(groups, GeometryRequirement.FINITE_ASSEMBLY,
                                    blocked=f"semantic group '{g}' is not available")
        atoms.update(grp.atom_ids)
    # a periodically extended component has no finite global geometry, whatever
    # the topology says (a build may even merge lipids into the protein molecule)
    for name, ctype in EXTENDED_GROUPS.items():
        grp = resolve_index_group(semantic_index.path, name)
        if grp is not None and atoms & set(grp.atom_ids):
            return AssemblyDecision(
                groups, GeometryRequirement.PERIODIC_EXTENDED,
                blocked=(f"the selection includes the periodically extended {ctype} component; "
                         f"a finite global geometry (radius of gyration, centre of mass, RMSD) is "
                         f"undefined for it — needs domain-specific (e.g. membrane-frame) geometry"))
    try:
        part = molecule_partition(system.topology_path, gmx)
    except MoleculePartitionError as exc:
        return AssemblyDecision(groups, GeometryRequirement.FINITE_ASSEMBLY,
                                blocked=f"molecule topology unreadable: {exc}")
    if part is None:
        return AssemblyDecision(groups, GeometryRequirement.SINGLE_MOLECULE_WHOLE,
                                evidence={"molecule_partition": "unknown (no .tpr)"})
    sel = part.select(atoms)
    ev = {"groups": list(groups), "molecules": sel.evidence()}
    if not sel.multi_molecule:
        return AssemblyDecision(groups, GeometryRequirement.SINGLE_MOLECULE_WHOLE, evidence=ev)
    closure = set(part.closure(atoms))
    for name, ctype in EXTENDED_GROUPS.items():
        grp = resolve_index_group(semantic_index.path, name)
        if grp is not None and closure & set(grp.atom_ids):
            return AssemblyDecision(
                groups, GeometryRequirement.PERIODIC_EXTENDED, evidence=ev,
                blocked=(f"the selection spans {sel.n_molecules} topological molecules and "
                         f"they include the periodically extended {ctype} component; a finite "
                         f"global geometry is undefined for it (needs domain-specific geometry)"))
    return AssemblyDecision(groups, GeometryRequirement.FINITE_ASSEMBLY, needed=True, evidence=ev)


def assembly_evidence(view) -> Optional[dict]:
    """Definition evidence of a clustered view (None for any other view)."""
    req = getattr(view, "requirements", None)
    if req is None or not req.cluster_groups:
        return None
    return {"geometry": GeometryRequirement.FINITE_ASSEMBLY, "cluster_groups": list(req.cluster_groups),
            "method": "gmx trjconv -pbc cluster on the complete molecules the groups span"}
