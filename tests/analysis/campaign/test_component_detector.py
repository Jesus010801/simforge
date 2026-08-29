"""Molecular-component inference — the scientifically critical path.

These tests are adversarial: the classifier must NEVER silently resolve a
receptor/peptide split on size alone, must always be explainable, and must
degrade to AMBIGUOUS / REVIEW_REQUIRED rather than guess.
"""
from __future__ import annotations

import pytest

from analysis.campaign.models import (
    ClassificationState, ComponentType, MolecularComponent,
)
from analysis.campaign.structure import component_detector as cd
from analysis.campaign.structure.component_detector import detect_components
from analysis.campaign.structure.spatial import MembraneEmbedding


def _by_type(res, ctype):
    return [c for c in res.components if c.component_type == ctype]


# ── Case A: single soluble chain ────────────────────────────────────────────

def test_single_protein_chain_no_membrane(minimal_pdb_single_chain):
    res = detect_components(structure_path=str(minimal_pdb_single_chain))
    receptors = _by_type(res, ComponentType.RECEPTOR)
    assert len(receptors) == 1
    assert receptors[0].classification_state == ClassificationState.RESOLVED
    assert receptors[0].chain_ids == ["A"]
    # explicit "no partner" marker, and NO peptide component
    others = _by_type(res, ComponentType.OTHER)
    assert any(o.label == "no partner" for o in others)
    assert _by_type(res, ComponentType.PEPTIDE) == []
    assert res.membrane_present is False
    assert receptors[0].evidence, "classification must be explainable"


# ── Case B3: soluble long + short => AMBIGUOUS, never silently resolved ──────

def test_receptor_peptide_soluble_is_ambiguous_not_resolved(minimal_pdb_receptor_peptide):
    res = detect_components(structure_path=str(minimal_pdb_receptor_peptide))
    receptor = _by_type(res, ComponentType.RECEPTOR)
    peptide = _by_type(res, ComponentType.PEPTIDE)
    assert len(receptor) == 1 and len(peptide) == 1

    assert receptor[0].classification_state == ClassificationState.AMBIGUOUS
    assert peptide[0].classification_state == ClassificationState.AMBIGUOUS
    # each carries a size-only caveat
    assert any("size only" in w or "relative size only" in w for w in receptor[0].warnings)
    assert peptide[0].warnings
    # evidence is present and explainable
    assert receptor[0].evidence
    kinds = {e.kind for e in receptor[0].evidence}
    assert "has_short_partner_chain" in kinds


def test_no_pure_size_classification_low_confidence(minimal_pdb_unequal):
    """A=200, B=60, both soluble non-TM: must stay AMBIGUOUS w/ combined evidence."""
    res = detect_components(structure_path=str(minimal_pdb_unequal))
    receptor = _by_type(res, ComponentType.RECEPTOR)[0]
    assert receptor.classification_state == ClassificationState.AMBIGUOUS
    assert receptor.confidence < 0.8
    kinds = {e.kind for e in receptor.evidence}
    assert "has_short_partner_chain" in kinds
    # the "larger chain" evidence alone must be weak (<= 0.15)
    larger = [e for e in receptor.evidence if e.kind == "larger_chain"]
    assert larger and larger[0].weight <= 0.15


# ── Case B1: homodimer ─────────────────────────────────────────────────────

def test_homodimer_is_oligomer_not_receptor_peptide(minimal_pdb_homodimer):
    res = detect_components(structure_path=str(minimal_pdb_homodimer))
    receptor = _by_type(res, ComponentType.RECEPTOR)
    assert len(receptor) == 1
    assert receptor[0].classification_state == ClassificationState.RESOLVED
    assert "-mer" in receptor[0].label
    assert set(receptor[0].chain_ids) == {"A", "B"}
    assert any(e.kind == "similar_chain_sizes" for e in receptor[0].evidence)
    assert _by_type(res, ComponentType.PEPTIDE) == []
    assert any(o.label == "no partner" for o in _by_type(res, ComponentType.OTHER))


# ── Case B2: one TM chain + short non-TM chains (monkeypatched embedding) ────

def test_membrane_receptor_plus_two_short_chains(minimal_pdb_three_chains, monkeypatch):
    fake = MembraneEmbedding(
        membrane_present=True, membrane_z_min=1.0, membrane_z_max=5.0,
        chain_span_fraction={"A": 0.92, "B": 0.08, "C": 0.10},
        note="synthetic",
    )
    monkeypatch.setattr(cd, "analyse_membrane_embedding", lambda _p: fake)

    res = detect_components(structure_path=str(minimal_pdb_three_chains))
    receptor = _by_type(res, ComponentType.RECEPTOR)[0]
    peptide = _by_type(res, ComponentType.PEPTIDE)[0]

    assert receptor.chain_ids == ["A"]
    assert receptor.classification_state == ClassificationState.RESOLVED
    assert peptide.classification_state == ClassificationState.AMBIGUOUS
    assert set(peptide.chain_ids) == {"B", "C"}
    assert res.membrane_present is True


# ── Real membrane system via the bundled trajectory's .gro ──────────────────

def test_real_membrane_gro(real_membrane_gro):
    res = detect_components(structure_path=str(real_membrane_gro))
    assert res.membrane_present is True
    receptor = _by_type(res, ComponentType.RECEPTOR)
    assert len(receptor) == 1
    assert receptor[0].classification_state == ClassificationState.RESOLVED
    assert _by_type(res, ComponentType.MEMBRANE)
    assert _by_type(res, ComponentType.WATER)
    assert _by_type(res, ComponentType.IONS)
    for c in res.components:
        if c.component_type in (ComponentType.WATER, ComponentType.IONS, ComponentType.MEMBRANE):
            assert c.classification_state == ClassificationState.RESOLVED


# ── solvent/ion/lipid never enter polymer classification ───────────────────

def test_solvent_excluded_from_polymer(minimal_gro_solvated_protein):
    res = detect_components(structure_path=str(minimal_gro_solvated_protein))
    receptor = _by_type(res, ComponentType.RECEPTOR)
    assert len(receptor) == 1
    # the receptor residue_count is the 50-residue protein, not 50+2000+10
    assert receptor[0].residue_count == 50
    water = _by_type(res, ComponentType.WATER)[0]
    assert water.residue_count == 2000


# ── serialisation round-trip ──────────────────────────────────────────────

def test_component_roundtrip(minimal_pdb_receptor_peptide):
    res = detect_components(structure_path=str(minimal_pdb_receptor_peptide))
    for comp in res.components:
        d = comp.to_dict()
        back = MolecularComponent.from_dict(d)
        assert back.to_dict() == d
        assert [e.kind for e in back.evidence] == [e.kind for e in comp.evidence]


def test_gro_transmembrane_receptor_evidence(real_membrane_gro):
    res = detect_components(structure_path=str(real_membrane_gro))
    receptor = _by_type(res, ComponentType.RECEPTOR)[0]
    kinds = {e.kind for e in receptor.evidence}
    assert "transmembrane_embedding" in kinds
    assert "membrane_system_soluble_chain" not in kinds


def test_missing_structure_returns_empty_not_crash():
    res = detect_components(structure_path=None)
    assert res.components == []
    assert res.warnings
    res2 = detect_components(structure_path="/no/such/file.pdb")
    assert res2.components == []
