# core/test_parser.py
from __future__ import annotations


# ── Parser / SystemState ──────────────────────────────────────────────────────

def test_parse_returns_correct_system_type(state):
    assert state.inferred_system_type == "competitive-inhibition"


def test_parse_has_three_components(state):
    assert len(state.components) == 3


def test_all_components_have_ids(state):
    for comp in state.components:
        assert comp.id


def test_component_ids_non_empty(state):
    assert len(state.component_ids()) == 3


def test_global_reasoning_present(state):
    assert state.global_reasoning is not None


def test_at_least_one_component_validated(state):
    assert state.global_reasoning.n_components_validated > 0


def test_protein_component_is_valid(state):
    protein = next((c for c in state.components if c.role == "protein"), None)
    assert protein is not None
    assert protein.validation is not None
    assert protein.validation.is_valid


def test_ligand_component_has_descriptors(state):
    # role is "competitive_ligand" in hmg_competition.yaml
    ligand = next(
        (c for c in state.components if "ligand" in c.role),
        None,
    )
    assert ligand is not None
    assert ligand.descriptors is not None


def test_substrate_component_has_descriptors(state):
    substrate = next((c for c in state.components if c.role == "substrate"), None)
    assert substrate is not None
    assert substrate.descriptors is not None


def test_collect_all_warnings_returns_list(state):
    warnings = state.collect_all_warnings()
    assert isinstance(warnings, list)


def test_collect_all_risks_returns_list(state):
    risks = state.collect_all_risks()
    assert isinstance(risks, list)


def test_forcefields_set(state):
    assert state.forcefields.protein
    assert state.forcefields.ligands


def test_has_membrane_is_bool(state):
    assert isinstance(state.has_membrane(), bool)
