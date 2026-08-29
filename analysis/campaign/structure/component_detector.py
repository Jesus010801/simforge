"""Evidence-combining molecular-component inference.

Scientific rules enforced here (see the SimForge study-campaign spec):

* **No unconditional size rule.**  "Largest chain = receptor" is at most one
  weak tie-breaker (weight 0.10) and never decides on its own.
* Classification is *explainable*: every :class:`MolecularComponent` carries a
  list of weighted :class:`ComponentEvidence`.
* When receptor vs. peptide cannot be separated confidently the component's
  ``classification_state`` is ``AMBIGUOUS`` (or ``REVIEW_REQUIRED`` once an
  analysis that needs the distinction is requested).  It is never guessed.

Inputs are whatever is cheaply available: the reference structure (chains,
residue sequences, coordinates), the ``.top`` ``[ molecules ]`` section, and the
membrane-embedding signal.  Sequence identity / databases are intentionally
out of scope for Phase 1 (no network dependency).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from analysis.campaign.models import (
    ClassificationState, ComponentEvidence, ComponentType, MolecularComponent,
)
from analysis.campaign.structure.parsers import (
    ChainInfo, StructureModel, TopologyMolecules, parse_structure, parse_top_molecules,
)
from analysis.campaign.structure.residue_classes import classify_resname
from analysis.campaign.structure.spatial import analyse_membrane_embedding

# A chain shorter than this (and clearly not transmembrane) is *peptide-like*.
PEPTIDE_MAX_RESIDUES = 60
# Two protein chains whose residue counts differ by less than this ratio are
# treated as "similar size" -> oligomer hypothesis, size cannot separate them.
OLIGOMER_SIZE_RATIO = 0.80


@dataclass
class ComponentDetectionResult:
    components: list[MolecularComponent]
    structure_model: Optional[StructureModel]
    membrane_present: Optional[bool]
    warnings: list[str]

    def by_type(self, ctype: str) -> Optional[MolecularComponent]:
        for c in self.components:
            if c.component_type == ctype:
                return c
        return None


def _confidence_from_evidence(ev: list[ComponentEvidence]) -> float:
    """Squash summed evidence weight to (0,1)."""
    total = sum(e.weight for e in ev)
    if total <= 0:
        return 0.0
    return round(min(0.99, 1.0 - 0.5 ** (total / 0.4)), 3)


def _composition_vector(chain: ChainInfo) -> dict[str, float]:
    counts: dict[str, int] = {}
    for _, rn in chain.residues:
        counts[rn.upper()] = counts.get(rn.upper(), 0) + 1
    total = sum(counts.values()) or 1
    return {k: v / total for k, v in counts.items()}


def _similar_composition(chains: list[ChainInfo], threshold: float = 0.15) -> bool:
    """True if every pair of chains has a small residue-composition L1 distance.

    Homo-oligomer chains share a sequence -> near-identical composition.  A
    hetero-complex (e.g. GPCR + G-protein) does not, even at similar length.
    """
    vecs = [_composition_vector(c) for c in chains]
    for i in range(len(vecs)):
        for j in range(i + 1, len(vecs)):
            keys = set(vecs[i]) | set(vecs[j])
            l1 = sum(abs(vecs[i].get(k, 0.0) - vecs[j].get(k, 0.0)) for k in keys)
            if l1 > threshold * 2:      # L1 of two distributions is in [0, 2]
                return False
    return True


def _solvent_components(
    structure: StructureModel, topmol: Optional[TopologyMolecules],
) -> list[MolecularComponent]:
    hist: dict[str, dict] = {}
    for chain in structure.chains:
        for resnum, rn in chain.residues:
            cls = classify_resname(rn)
            if cls in ("water", "ion", "lipid", "cofactor"):
                d = hist.setdefault(cls, {"resnames": set(), "residues": 0})
                d["resnames"].add(rn.upper())
                d["residues"] += 1

    mapping = {
        "water": (ComponentType.WATER, "water", "resname SOL HOH WAT TIP3 SPC"),
        "ion":   (ComponentType.IONS, "ions", "resname NA CL K MG CA ZN"),
        "lipid": (ComponentType.MEMBRANE, "membrane", "lipid residues"),
        "cofactor": (ComponentType.COFACTOR, "cofactor", "cofactor residues"),
    }
    out: list[MolecularComponent] = []
    for cls, (ctype, label, sel) in mapping.items():
        if cls not in hist:
            continue
        d = hist[cls]
        out.append(MolecularComponent(
            component_type=ctype, label=label,
            selection=sel, resnames=sorted(d["resnames"]),
            residue_count=d["residues"],
            confidence=0.98,
            classification_state=ClassificationState.RESOLVED,
            evidence=[ComponentEvidence(
                kind="known_residue_class",
                detail=f"{d['residues']} residues match the {cls} residue-name set",
                weight=1.0, value=sorted(d["resnames"]),
            )],
        ))
    return out


def _polymer_chains(structure: StructureModel) -> list[ChainInfo]:
    return [c for c in structure.chains
            if c.residue_count > 0 and c.polymer_fraction() >= 0.5
            and c.dominant_class() in ("amino_acid", "nucleotide")]


def detect_components(
    *,
    structure_path: Optional[str | Path],
    topology_path: Optional[str | Path] = None,
    membrane_present_hint: Optional[bool] = None,
) -> ComponentDetectionResult:
    warnings: list[str] = []
    if not structure_path or not Path(structure_path).is_file():
        return ComponentDetectionResult(
            components=[], structure_model=None,
            membrane_present=membrane_present_hint,
            warnings=["no reference structure available for component inference"],
        )

    try:
        structure = parse_structure(Path(structure_path))
    except ValueError as exc:
        return ComponentDetectionResult(
            components=[], structure_model=None,
            membrane_present=membrane_present_hint, warnings=[str(exc)],
        )
    warnings.extend(structure.warnings)

    topmol: Optional[TopologyMolecules] = None
    if topology_path and Path(topology_path).is_file() and Path(topology_path).suffix.lower() == ".top":
        topmol = parse_top_molecules(Path(topology_path))
        warnings.extend(topmol.warnings)

    embedding = analyse_membrane_embedding(structure_path)
    membrane_present = embedding.membrane_present
    if membrane_present_hint is not None and membrane_present_hint != membrane_present:
        warnings.append(
            f"membrane presence hint ({membrane_present_hint}) disagrees with "
            f"structure evidence ({membrane_present}); using structure evidence"
        )

    components: list[MolecularComponent] = _solvent_components(structure, topmol)

    poly = _polymer_chains(structure)
    if not poly:
        warnings.append("no protein/nucleic polymer chain detected")
        return ComponentDetectionResult(components, structure, membrane_present, warnings)

    tm_chains = set(embedding.transmembrane_chains()) if membrane_present else set()
    max_res = max(c.residue_count for c in poly)

    # ── Case A: single polymer chain ────────────────────────────────────────
    if len(poly) == 1:
        chain = poly[0]
        ev = [ComponentEvidence("sole_polymer_chain",
                                "only one amino-acid/nucleic polymer chain present",
                                0.6)]
        ctype = ComponentType.RECEPTOR
        label = "primary protein"
        if chain.chain_id in tm_chains:
            ev.append(ComponentEvidence(
                "transmembrane_embedding",
                f"chain backbone spans {embedding.chain_span_fraction.get(chain.chain_id):.0%} "
                f"of the bilayer", 0.5))
        elif membrane_present:
            ev.append(ComponentEvidence(
                "membrane_system_soluble_chain",
                "membrane present but this chain does not span it", 0.1))
        if chain.dominant_class() == "nucleic_acid":
            ctype, label = ComponentType.NUCLEIC_ACID, "primary nucleic acid"
        ev.append(ComponentEvidence(
            "study_role_assumed",
            "used as the 'receptor' study role for receptor-based analyses; this is a "
            "study-level role assignment, not a verified biological receptor identity",
            0.0))
        components.append(MolecularComponent(
            component_type=ctype, label=label,
            selection="Receptor", chain_ids=[chain.chain_id],
            resnames=sorted(chain.resname_set), residue_count=chain.residue_count,
            confidence=_confidence_from_evidence(ev),
            classification_state=ClassificationState.RESOLVED,
            study_role="receptor (assumed — sole polymer chain)",
            evidence=ev,
        ))
        # explicit: there is no partner
        components.append(_no_partner_marker("only one polymer chain in the system"))
        # a single protein is still a valid "complex" of one (rmsd-complex == rmsd-receptor here)
        _append_complex(components, [chain])
        return ComponentDetectionResult(components, structure, membrane_present, warnings)

    # ── Case B: multiple polymer chains ────────────────────────────────────
    # classify each chain independently, then reconcile.
    scored: list[tuple[ChainInfo, dict]] = []
    for chain in poly:
        rc = chain.residue_count
        rel = rc / max_res if max_res else 1.0
        similar_peers = [
            o for o in poly if o is not chain
            and min(o.residue_count, rc) / max(o.residue_count, rc) >= OLIGOMER_SIZE_RATIO
        ]
        scored.append((chain, {
            "rel_size": rel,
            "is_tm": chain.chain_id in tm_chains,
            "short": rc <= PEPTIDE_MAX_RESIDUES,
            "similar_peers": [c.chain_id for c in similar_peers],
            "independent_chain": True,
        }))

    tm_present = any(s["is_tm"] for _, s in scored)
    all_similar = all(len(s["similar_peers"]) == len(poly) - 1 for _, s in scored)
    same_composition = _similar_composition(poly)

    # B1: homo-oligomer — every chain similar in SIZE *and* COMPOSITION (so a
    # GPCR + G-protein complex, similar in size but not sequence, is NOT mistaken
    # for a receptor oligomer), and — if membrane — all transmembrane.
    if (all_similar and same_composition and len(poly) >= 2
            and (not membrane_present or all(s["is_tm"] for _, s in scored))):
        ev = [
            ComponentEvidence("similar_chain_sizes",
                              f"all {len(poly)} chains within {OLIGOMER_SIZE_RATIO:.0%} residue count",
                              0.4),
            ComponentEvidence("similar_amino_acid_composition",
                              "chains have near-identical residue-name composition "
                              "(consistent with a homo-oligomer, not a hetero-complex)",
                              0.4),
        ]
        if membrane_present:
            ev.append(ComponentEvidence("all_chains_transmembrane",
                                        "every chain spans the bilayer", 0.4))
        components.append(MolecularComponent(
            component_type=ComponentType.RECEPTOR,
            label=f"receptor ({len(poly)}-mer)",
            selection="Receptor",
            chain_ids=[c.chain_id for c in poly],
            residue_count=sum(c.residue_count for c in poly),
            confidence=_confidence_from_evidence(ev),
            classification_state=ClassificationState.RESOLVED,
            evidence=ev,
            warnings=["treated as a homo-oligomeric receptor; there is no distinct peptide/ligand partner"],
        ))
        components.append(_no_partner_marker("homo-oligomeric receptor; no distinct partner chain"))
        _append_complex(components, poly)
        return ComponentDetectionResult(components, structure, membrane_present, warnings)

    # B2: one clear membrane-embedded chain + short non-TM chain(s)
    tm_poly = [c for c, s in scored if s["is_tm"]]
    non_tm_short = [c for c, s in scored if not s["is_tm"] and s["short"]]
    non_tm_long = [c for c, s in scored if not s["is_tm"] and not s["short"]]

    if len(tm_poly) == 1 and non_tm_short and not non_tm_long:
        receptor_chain = tm_poly[0]
        r_ev = [
            ComponentEvidence("transmembrane_embedding",
                              f"backbone spans "
                              f"{embedding.chain_span_fraction.get(receptor_chain.chain_id):.0%} of the bilayer",
                              0.6),
            ComponentEvidence("large_relative_to_partner",
                              "much larger than the non-transmembrane chain(s)", 0.1),
        ]
        components.append(MolecularComponent(
            component_type=ComponentType.RECEPTOR, label="receptor",
            selection="Receptor", chain_ids=[receptor_chain.chain_id],
            resnames=sorted(receptor_chain.resname_set),
            residue_count=receptor_chain.residue_count,
            confidence=_confidence_from_evidence(r_ev),
            classification_state=ClassificationState.RESOLVED, evidence=r_ev,
        ))
        if len(non_tm_short) == 1:
            pep = non_tm_short[0]
            p_ev = [
                ComponentEvidence("independent_chain", "separate chain from the receptor", 0.25),
                ComponentEvidence("amino_acid_polymer",
                                  f"{pep.residue_count}-residue amino-acid chain", 0.2),
                ComponentEvidence("non_transmembrane", "does not span the bilayer", 0.25),
                ComponentEvidence("small_relative_to_membrane_receptor",
                                  "far shorter than the transmembrane receptor", 0.15),
            ]
            components.append(MolecularComponent(
                component_type=ComponentType.PEPTIDE, label="peptide",
                selection="Peptide", chain_ids=[pep.chain_id],
                resnames=sorted(pep.resname_set), residue_count=pep.residue_count,
                confidence=_confidence_from_evidence(p_ev),
                classification_state=ClassificationState.RESOLVED, evidence=p_ev,
            ))
        else:
            opts = [c.chain_id for c in non_tm_short]
            components.append(MolecularComponent(
                component_type=ComponentType.PEPTIDE, label="peptide (ambiguous)",
                selection="Peptide", chain_ids=opts,
                confidence=0.4,
                classification_state=ClassificationState.AMBIGUOUS,
                evidence=[ComponentEvidence(
                    "multiple_short_non_tm_chains",
                    f"{len(opts)} short non-transmembrane chains are all peptide-compatible: {opts}",
                    0.0, value=opts)],
                warnings=["multiple peptide-compatible chains; assign one in study_manifest.yaml"],
            ))
        _append_complex(components, poly)
        return ComponentDetectionResult(components, structure, membrane_present, warnings)

    # B3: soluble multi-chain, one long + short(s) — weak, keep AMBIGUOUS
    if non_tm_long and non_tm_short and not tm_present:
        long_chain = max(non_tm_long, key=lambda c: c.residue_count)
        r_ev = [
            ComponentEvidence("larger_chain",
                              "largest polymer chain (weak evidence only, no membrane anchor)", 0.1),
            ComponentEvidence("has_short_partner_chain",
                              "co-exists with a much shorter independent chain", 0.15),
        ]
        components.append(MolecularComponent(
            component_type=ComponentType.RECEPTOR, label="receptor (tentative)",
            selection="Receptor", chain_ids=[long_chain.chain_id],
            resnames=sorted(long_chain.resname_set), residue_count=long_chain.residue_count,
            confidence=_confidence_from_evidence(r_ev),
            classification_state=ClassificationState.AMBIGUOUS, evidence=r_ev,
            warnings=["soluble complex: receptor/partner roles inferred from size only — verify"],
        ))
        for pep in non_tm_short:
            p_ev = [
                ComponentEvidence("independent_chain", "separate short chain", 0.2),
                ComponentEvidence("non_transmembrane", "no membrane to anchor to", 0.1),
                ComponentEvidence("small_relative_to_largest_chain",
                                  "much shorter than the largest chain", 0.1),
            ]
            components.append(MolecularComponent(
                component_type=ComponentType.PEPTIDE, label="peptide (tentative)",
                selection="Peptide", chain_ids=[pep.chain_id],
                resnames=sorted(pep.resname_set), residue_count=pep.residue_count,
                confidence=_confidence_from_evidence(p_ev),
                classification_state=ClassificationState.AMBIGUOUS, evidence=p_ev,
                warnings=["peptide role inferred from relative size only — verify"],
            ))
        _append_complex(components, poly)
        return ComponentDetectionResult(components, structure, membrane_present, warnings)

    # B4: genuinely can't tell — every polymer chain a candidate
    opts = [c.chain_id for c in poly]
    components.append(MolecularComponent(
        component_type=ComponentType.RECEPTOR, label="protein chains (unresolved)",
        selection="Receptor", chain_ids=opts,
        residue_count=sum(c.residue_count for c in poly),
        confidence=0.2,
        classification_state=ClassificationState.REVIEW_REQUIRED,
        evidence=[ComponentEvidence(
            "multiple_compatible_protein_chains",
            f"chains {opts} with residue counts "
            f"{[c.residue_count for c in poly]} — no membrane anchor and no dominant signal "
            f"to separate receptor from partner",
            0.0, value={c.chain_id: c.residue_count for c in poly})],
        warnings=["assign receptor / partner chains in study_manifest.yaml before "
                  "running analyses that need the distinction"],
    ))
    _append_complex(components, poly)
    return ComponentDetectionResult(components, structure, membrane_present, warnings)


def _no_partner_marker(reason: str) -> MolecularComponent:
    return MolecularComponent(
        component_type=ComponentType.OTHER, label="no partner",
        selection="", confidence=1.0,
        classification_state=ClassificationState.RESOLVED,
        evidence=[ComponentEvidence("no_partner_chain", reason, 1.0)],
    )


def _append_complex(components: list[MolecularComponent], poly: list[ChainInfo]) -> None:
    components.append(MolecularComponent(
        component_type=ComponentType.COMPLEX, label="complex",
        selection="Complex",
        chain_ids=[c.chain_id for c in poly],
        residue_count=sum(c.residue_count for c in poly),
        confidence=0.95,
        classification_state=ClassificationState.RESOLVED,
        evidence=[ComponentEvidence("all_polymer_chains",
                                    "union of all amino-acid/nucleic chains", 1.0)],
    ))
