from decimal import Decimal
import pytest
from simforge.knowledge.domain import ExperimentalStructure
from test_identity import base, external, ref, provenance, artifact


def structure():
    source = external()
    source.update(entity=ref(10, 'experimental_structure'), mapping_status='CONFIRMED_MAPPING', evidence_refs=[ref(11, 'evidence')])
    return dict(**base(10), origin='experimental', source_identifier=source, method=dict(category='X_RAY_DIFFRACTION', raw_method='X-RAY DIFFRACTION', declared_origin='EXPERIMENTAL', classification_basis='SOURCE_DECLARATION', classification_claim=ref(12, 'claim'), provenance=provenance()), artifact=artifact(), evidence_refs=[ref(11, 'evidence')])


@pytest.mark.parametrize('label,category', [('cryo-EM', 'ELECTRON_MICROSCOPY'), ('cryo-electron microscopy', 'ELECTRON_MICROSCOPY'), ('solution NMR', 'SOLUTION_NMR'), ('neutron diffraction', 'NEUTRON_DIFFRACTION')])
def test_known_method_normalization(label, category):
    data = structure()
    data['method'].update(category=label, raw_method=label)
    result = ExperimentalStructure(**data)
    assert result.method.category == category
    assert result.method.raw_method == label


def test_extensible_experimental_method_requires_explicit_attribution():
    data = structure()
    data['method'].update(category='OTHER_EXPERIMENTAL', raw_method='Novel experimental technique')
    assert ExperimentalStructure(**data).method.category == 'OTHER_EXPERIMENTAL'
    for field in ('classification_claim', 'classification_basis', 'provenance'):
        broken = dict(data, method=dict(data['method']))
        del broken['method'][field]
        with pytest.raises(ValueError):
            ExperimentalStructure(**broken)


@pytest.mark.parametrize('origin', ['UNKNOWN', 'PREDICTED', 'COMPUTATIONAL', 'THEORETICAL'])
def test_nonexperimental_origin_rejected(origin):
    data = structure()
    data['method'].update(category='OTHER_EXPERIMENTAL', raw_method='Unrecognized method', declared_origin=origin)
    with pytest.raises(ValueError):
        ExperimentalStructure(**data)


@pytest.mark.parametrize('label', ['unknown', 'unspecified', 'theoretical model', 'predicted structure', 'computational model', 'AlphaFold prediction', 'homology model'])
def test_explicit_prediction_cannot_use_extension_path(label):
    data = structure()
    data['method'].update(category='OTHER_EXPERIMENTAL', raw_method=label)
    with pytest.raises(ValueError):
        ExperimentalStructure(**data)


def test_source_mapping_and_resolution():
    data = structure()
    result = ExperimentalStructure(**data, resolution_angstrom=Decimal('2.50'))
    assert ExperimentalStructure.from_canonical_bytes(result.to_canonical_bytes()) == result
    for resolution in (Decimal('0'), Decimal('-1'), Decimal('NaN'), 2.5):
        with pytest.raises(ValueError):
            ExperimentalStructure(**data, resolution_angstrom=resolution)
    data['source_identifier']['entity']['entity_revision'] = 2
    with pytest.raises(ValueError):
        ExperimentalStructure(**data)


@pytest.mark.parametrize('label', ['in silico', 'molecular dynamics', 'comparative modeling'])
def test_known_computational_labels_cannot_be_declared_experimental(label):
    data = structure()
    data['method'].update(category='OTHER_EXPERIMENTAL', raw_method=label)
    with pytest.raises(ValueError):
        ExperimentalStructure(**data)


@pytest.mark.parametrize('label', [
    'molecular dynamics simulation', 'molecular-dynamics simulation', 'molecular dynamics model',
    'MD simulation', 'in silico model', 'in-silico model', 'computational model',
    'computational prediction', 'predicted model', 'theoretical model', 'homology model',
    'comparative model', 'comparative modelling', 'AlphaFold prediction', 'AlphaFold model',
    'machine-learning prediction', 'machine-learning model', 'MOLECULAR—DYNAMICS model',
    'experiment followed by computational prediction',
])
def test_computational_phrase_variants_cannot_use_extension(label):
    data = structure()
    data['method'].update(category='OTHER_EXPERIMENTAL', raw_method=label)
    with pytest.raises(ValueError):
        ExperimentalStructure(**data)


@pytest.mark.parametrize('basis', ['USER_ASSERTION', 'CURATED_ASSERTION'])
def test_extension_requires_source_declaration_not_assertion(basis):
    data = structure()
    data['method'].update(category='OTHER_EXPERIMENTAL', raw_method='Serial femtosecond crystallography', classification_basis=basis)
    with pytest.raises(ValueError):
        ExperimentalStructure(**data)


@pytest.mark.parametrize('label', ['unknown technique', 'not determined', 'unclassified method',
    'ambiguous', 'unclear', 'experimental?', 'possibly experimental', 'N/A', 'model', '???'])
def test_ambiguous_extension_labels_rejected(label):
    data = structure()
    data['method'].update(category='OTHER_EXPERIMENTAL', raw_method=label)
    with pytest.raises(ValueError):
        ExperimentalStructure(**data)


def test_legitimate_extension_retains_source_classification_and_provenance():
    data = structure()
    data['method'].update(category='OTHER_EXPERIMENTAL', raw_method='Serial femtosecond crystallography')
    obj = ExperimentalStructure(**data)
    assert obj.method.raw_method == 'Serial femtosecond crystallography'
    assert obj.method.classification_basis == 'SOURCE_DECLARATION'
    assert obj.method.provenance.source_record.record_revision == '2'
    assert obj.method.classification_claim.entity_id.value == ref(12, 'claim')['entity_id']
    assert ExperimentalStructure.from_canonical_bytes(obj.to_canonical_bytes()) == obj


# Independent source-like descriptions, rather than production regex fragments.
_CONTRADICTORY_METHODS = [
    'non-experimental model', 'non experimental model', 'nonexperimental structure',
    'not experimental', 'not an experimental method', 'no experimental method',
    'no experimental structure', 'AlphaFold model', 'AlphaFold2 model', 'AlphaFold3 prediction',
    'computer-generated structural model', 'computer generated structure',
    'computational model', 'computational structure', 'computational prediction',
    'in silico model', 'in-silico structure', 'simulated structure', 'simulation-derived model',
    'predicted structure', 'structure prediction', 'protein structure prediction',
    'AI-predicted structure', 'homology model', 'homology modelling', 'homology modeling',
    'comparative model', 'comparative modelling', 'comparative modeling',
    'theoretical model', 'theoretical structure', 'molecular dynamics simulation',
    'molecular-dynamics model', 'MD simulation', 'MD-derived structure',
    'machine-learning prediction', 'machine learning model', 'model-derived structure',
]


@pytest.mark.parametrize('label', _CONTRADICTORY_METHODS)
@pytest.mark.parametrize('variant', ['original', 'uppercase', 'whitespace', 'punctuation'])
def test_source_declaration_cannot_override_semantic_contradiction(label, variant):
    if variant == 'uppercase':
        label = label.upper()
    elif variant == 'whitespace':
        label = '  ' + label.replace('-', ' ').replace(' ', ' \t  ') + '  '
    elif variant == 'punctuation':
        label = label.replace('-', '—').replace(' ', ' /_ ')
    data = structure()
    data['method'].update(category='OTHER_EXPERIMENTAL', raw_method=label,
                          declared_origin='EXPERIMENTAL', classification_basis='SOURCE_DECLARATION')
    with pytest.raises(ValueError, match='contradicts experimental origin'):
        ExperimentalStructure(**data)


@pytest.mark.parametrize('label,category', [
    ('X-ray diffraction', 'X_RAY_DIFFRACTION'),
    ('electron microscopy', 'ELECTRON_MICROSCOPY'),
    ('cryo-EM', 'ELECTRON_MICROSCOPY'),
    ('solution NMR', 'SOLUTION_NMR'),
    ('solid-state NMR', 'SOLID_STATE_NMR'),
    ('neutron diffraction', 'NEUTRON_DIFFRACTION'),
    ('electron diffraction', 'ELECTRON_DIFFRACTION'),
    ('fiber diffraction', 'FIBER_DIFFRACTION'),
    ('fibre diffraction', 'FIBER_DIFFRACTION'),
])
def test_contradiction_classifier_preserves_known_methods(label, category):
    data = structure()
    data['method'].update(category=label, raw_method=label)
    obj = ExperimentalStructure(**data)
    assert obj.method.category == category
    assert obj.method.raw_method == label
    assert ExperimentalStructure.from_canonical_bytes(obj.to_canonical_bytes()) == obj


@pytest.mark.parametrize('label', ['Serial  femtosecond—crystallography', 'Neutron spin-echo spectroscopy'])
def test_source_declared_extension_preserves_raw_label_and_full_provenance(label):
    data = structure()
    data['method'].update(category='OTHER_EXPERIMENTAL', raw_method=label)
    original_provenance = ExperimentalStructure(**structure()).method.provenance
    obj = ExperimentalStructure(**data)
    assert obj.method.raw_method == label
    assert obj.method.classification_basis == 'SOURCE_DECLARATION'
    assert obj.method.provenance == original_provenance
    assert obj.method.classification_claim.entity_id.value == ref(12, 'claim')['entity_id']
    restored = ExperimentalStructure.from_canonical_bytes(obj.to_canonical_bytes())
    assert restored == obj and restored.method.raw_method == label


@pytest.mark.parametrize('label', ['AlphaFoldv2 coordinates', 'ALPHAFOLDv31 output',
                                   'AlphaFold17.4 coordinates', 'alpha.fold_v8 coordinates'])
def test_numbered_prediction_family_is_not_experimental(label):
    data = structure()
    data['method'].update(category='OTHER_EXPERIMENTAL', raw_method=label)
    with pytest.raises(ValueError, match='contradicts experimental origin: prediction'):
        ExperimentalStructure(**data)
