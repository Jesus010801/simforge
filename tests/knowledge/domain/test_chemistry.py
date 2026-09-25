import pytest
from simforge.knowledge.domain import ChemicalIdentity, ChemicalSpecies
from test_identity import base, ref, provenance


def representation():
    return dict(representation_type='chem:smiles', value='CCO', specificity='PARTIAL', provenance=provenance())


def chemical(n=1):
    return dict(**base(n), representations=[representation()], stereochemistry_status='UNSPECIFIED', isotope_status='UNSPECIFIED')


def species():
    return dict(**base(2), chemical_identity=ref(1, 'chemical_identity'), representations=[representation()], formal_charge=None, protonation_state=None, tautomer_state=None, description_status='PARTIAL')


def test_matching_chemical_strings_do_not_merge_identity():
    a, b = ChemicalIdentity(**chemical()), ChemicalIdentity(**chemical(3))
    assert a.representations == b.representations
    assert a.sf_id != b.sf_id
    assert a != b


def test_species_is_separate_and_unknown_charge_is_not_zero():
    data = species()
    unknown = ChemicalSpecies(**data)
    zero = ChemicalSpecies(**(data | dict(formal_charge=0)))
    assert unknown != zero
    assert not isinstance(unknown, ChemicalIdentity)
    with pytest.raises(ValueError):
        ChemicalSpecies(**(data | dict(chemical_identity=ref(1, 'protein'))))
    with pytest.raises(ValueError):
        ChemicalSpecies(**(data | dict(sf_id=ref(1)['entity_id'])))


def test_specification_status_does_not_compute_chemistry():
    data = species()
    with pytest.raises(ValueError):
        ChemicalSpecies(**(data | dict(description_status='SPECIFIED')))
    data.update(description_status='SPECIFIED', formal_charge=0, protonation_state='source-specified', tautomer_state='not applicable')
    data['representations'][0]['value'] = 'a source assertion, not validated SMILES'
    assert ChemicalSpecies(**data).formal_charge == 0
    with pytest.raises(ValueError):
        ChemicalSpecies(**(data | dict(formal_charge=True)))


@pytest.mark.parametrize('field', ['representation_type', 'value', 'specificity', 'provenance'])
def test_representation_requires_attribution_and_specificity(field):
    data = chemical()
    del data['representations'][0][field]
    with pytest.raises(ValueError):
        ChemicalIdentity(**data)
