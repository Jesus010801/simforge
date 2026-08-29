"""Observable registry contract."""
from __future__ import annotations

import pytest

from analysis.campaign.models import ComponentType, ViewKind
from analysis.campaign.observables import registry
from analysis.campaign.observables.base import ObservableSpec

_BUILTINS = [
    "rmsd-complex", "rmsd-peptide-intrinsic", "rmsd-peptide-receptor-frame", "rmsd-receptor",
]


def test_ids_are_the_four_builtins_sorted():
    assert registry.ids() == _BUILTINS


def test_get_returns_spec_with_required_components():
    spec = registry.get("rmsd-receptor")
    assert spec.required_components == (ComponentType.RECEPTOR,)
    assert registry.get("rmsd-peptide-receptor-frame").required_components == (
        ComponentType.RECEPTOR, ComponentType.PEPTIDE,
    )


def test_get_unknown_raises_with_valid_ids():
    with pytest.raises(registry.UnknownObservableError) as ei:
        registry.get("nope")
    msg = str(ei.value)
    assert "rmsd-receptor" in msg


def test_register_duplicate_raises():
    with pytest.raises(registry.DuplicateObservableError):
        registry.register(registry.get("rmsd-receptor"))


def test_register_empty_id_raises():
    class Bad(ObservableSpec):
        id = ""
    with pytest.raises(ValueError):
        registry.register(Bad())


def test_catalog_entry_shape():
    for entry in registry.catalog():
        assert set(entry) >= {
            "id", "display_name", "description", "category",
            "required_components", "supports_replicate_aggregation",
        }
        assert isinstance(entry["required_components"], list)
        assert isinstance(entry["supports_replicate_aggregation"], bool)


def test_register_then_manual_restore_yields_builtins():
    class Fake(ObservableSpec):
        id = "fake-observable"
        required_components = (ComponentType.RECEPTOR,)

    registry.register(Fake())
    assert registry.is_registered("fake-observable")

    # manual restore (what a working _reset_for_tests would achieve)
    registry._REGISTRY.pop("fake-observable")
    assert registry.ids() == _BUILTINS


def test_reset_for_tests_restores_exactly_the_builtins():
    class Fake(ObservableSpec):
        id = "fake-observable"
        required_components = (ComponentType.RECEPTOR,)

    registry.register(Fake())
    registry._reset_for_tests()
    assert registry.ids() == _BUILTINS
    assert not registry.is_registered("fake-observable")


def test_trajectory_requirements_view_kinds():
    frame = registry.get("rmsd-peptide-receptor-frame").trajectory_requirements({})
    assert frame.view_kind() == ViewKind.FITTED
    assert frame.fit_selection == "Receptor_Backbone"

    recep = registry.get("rmsd-receptor").trajectory_requirements({})
    assert recep.view_kind() == ViewKind.WHOLE
    assert recep.fit_selection is None
    assert recep.requires_whole_molecules is True
