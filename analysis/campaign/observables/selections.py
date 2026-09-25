"""One selection reference for generic observables — no selection language.

    component:<type>     a RESOLVED molecular component (semantic index group)
    annotation:<id>      an ACTIVE persistent annotation (``Ann_<id>`` group)

Resolution never guesses: an unresolved component, an inactive / unresolved
annotation, or a semantic group whose atoms disagree with the annotation
record is reported with its reason.  The resolved evidence carries the exact
atom-set identity, so observable definition tokens change exactly when the
atoms they measure change.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from analysis.campaign.models import ClassificationState, SemanticIndex, SystemRecord

#: component type -> semantic index group (index_builder naming)
COMPONENT_GROUPS = {
    "receptor": "Receptor", "peptide": "Peptide", "protein_partner": "ProteinPartner",
    "complex": "Complex", "ligand": "Ligand", "cofactor": "Cofactor",
    "nucleic_acid": "NucleicAcid", "membrane": "Membrane", "system": "System",
}


@dataclass(frozen=True)
class SelectionRef:
    kind: str                          # "component" | "annotation"
    name: str

    @classmethod
    def parse(cls, text: str) -> "SelectionRef":
        kind, sep, name = str(text).partition(":")
        if not sep or kind not in ("component", "annotation") or not name:
            raise ValueError(f"selection {text!r} must be 'component:<type>' or "
                             f"'annotation:<id>'")
        if kind == "component" and name not in COMPONENT_GROUPS:
            raise ValueError(f"unknown component type {name!r}; one of {sorted(COMPONENT_GROUPS)}")
        return cls(kind, name)

    def __str__(self) -> str:
        return f"{self.kind}:{self.name}"


@dataclass
class ResolvedSelection:
    ref: SelectionRef
    ok: bool
    group: Optional[str] = None
    reason: str = ""
    evidence: dict = field(default_factory=dict)


def resolve_selection(ref: SelectionRef, system: SystemRecord,
                      semantic_index: Optional[SemanticIndex]) -> ResolvedSelection:
    from analysis.campaign.annotations import annotation_evidence
    from analysis.campaign.results import index_group_evidence

    def fail(reason: str) -> ResolvedSelection:
        return ResolvedSelection(ref, False, reason=reason)

    ndx = semantic_index.path if semantic_index is not None and semantic_index.path else None
    if ref.kind == "component":
        group = COMPONENT_GROUPS[ref.name]
        if ref.name != "system":
            comps = [c for c in system.components if c.component_type == ref.name]
            if len(comps) != 1:
                return fail(f"component '{ref.name}': {len(comps)} components of this type")
            if comps[0].classification_state != ClassificationState.RESOLVED:
                return fail(f"component '{ref.name}' is {comps[0].classification_state}, "
                            f"not RESOLVED")
        atoms = index_group_evidence(ndx, group) if ndx else None
        if atoms is None:
            return fail(f"semantic index group '{group}' is not available")
        return ResolvedSelection(ref, True, group, evidence={
            "type": "component", "component": ref.name, "group": group, "atoms": atoms})

    rec = next((a for a in system.annotations if a.annotation_id == ref.name), None)
    if rec is None:
        return fail(f"annotation '{ref.name}' is not declared for this system")
    ev = annotation_evidence(system.annotations, ref.name)
    if ev is None:
        why = "; ".join(rec.reasons) or "no reason recorded"
        return fail(f"annotation '{ref.name}' is {rec.state} ({why}) — only ACTIVE "
                    f"annotations can be measured")
    if rec.category != "residue_set" or not rec.group_name:
        return fail(f"annotation '{ref.name}' ({rec.category}) is not an atom selection")
    atoms = index_group_evidence(ndx, rec.group_name) if ndx else None
    if atoms is None:
        return fail(f"semantic index group '{rec.group_name}' is not available")
    if atoms["atoms_sha256"] != rec.atoms_sha256:
        return fail(f"semantic index group '{rec.group_name}' does not match the resolved "
                    f"annotation atoms (stale index?)")
    return ResolvedSelection(ref, True, rec.group_name, evidence={
        "type": "annotation", "annotation": ev, "group": rec.group_name})
