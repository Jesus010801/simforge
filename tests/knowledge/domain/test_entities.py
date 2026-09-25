import pytest
from simforge.knowledge.domain import Entity, Protein
from test_identity import base, provenance


@pytest.mark.parametrize('revision', [0, -1, True, '1', 1.0])
def test_internal_revision_is_strict(revision):
    with pytest.raises(ValueError):
        Entity(**(base() | dict(entity_type='protein', entity_revision=revision)))


def test_descriptor_and_concrete_type_are_distinct():
    descriptor = Entity(**base(), entity_type='protein')
    assert not isinstance(descriptor, Protein)
    with pytest.raises(ValueError):
        Protein.model_validate(descriptor)
    with pytest.raises(ValueError):
        Protein(**base(), entity_type='chemical_identity', provenance=[provenance()])


@pytest.mark.parametrize('status', ['ACTIVE', 'DEPRECATED', 'RETRACTED'])
def test_lifecycle_is_preserved(status):
    assert Entity(**base(), entity_type='protein', lifecycle=status).lifecycle == status


def test_extra_fields_are_forbidden():
    with pytest.raises(ValueError):
        Entity(**base(), entity_type='protein', engine='gromacs')
