"""Strict constructors, nested errors, and defensive scalar ownership."""
from dataclasses import replace, FrozenInstanceError
from datetime import datetime, timezone, timedelta, tzinfo
from uuid import UUID
import pytest
from simforge.knowledge.domain import Entity, Protein
from simforge.knowledge.identity import (
    RegistryChange, RegistryWrite, RegistryRecord, ExternalIdentifierKey,
    ExternalIdentifierMatch, InvalidRegistryOperationError, EntityNotFoundError,
    RegistryRevisionNotFoundError, EntityRevisionNotFoundError, RevisionConflictError,
    OperationConflictError,
)
from registry_contract import NOW, eid, descriptor, alias, change, command, invalid


@pytest.mark.parametrize('value', [True, False, 1.0, '1', -1, None])
def test_expected_revision_strict(value):
    invalid(lambda: replace(command(), expected_registry_revision=value), 'INVALID_ARGUMENT', ('expected_registry_revision',))


@pytest.mark.parametrize('value', [True, False, 1.0, '1', -1, 0, None])
@pytest.mark.parametrize('constructor', [lambda n: RegistryRecord(n, descriptor(), (), change()),
    lambda n: ExternalIdentifierMatch(descriptor(), n, alias())])
def test_stored_revision_strict(value, constructor):
    invalid(lambda: constructor(value), 'INVALID_ARGUMENT', ('registry_revision',))


@pytest.mark.parametrize('value', [None, '', ' ', 10, '\ud800'])
@pytest.mark.parametrize('field', ['actor', 'reason'])
def test_change_text(value, field):
    invalid(lambda: replace(change(), **{field: value}), 'INVALID_ARGUMENT', (field,))


@pytest.mark.parametrize('value', [None, '2026-01-02', NOW.replace(tzinfo=None), 1])
def test_time_strict(value):
    invalid(lambda: replace(change(), recorded_at=value), 'INVALID_ARGUMENT', ('recorded_at',))


@pytest.mark.parametrize('value', [None, str(change().operation_id), 1, UUID(int=0),
    UUID('10000000-0000-1000-8000-000000000001'), UUID('10000000-0000-4000-0000-000000000001')])
def test_operation_id_strict(value):
    invalid(lambda: replace(change(), operation_id=value), 'INVALID_ARGUMENT', ('operation_id',))


@pytest.mark.parametrize('value', ['UniProt', ' uniprot', 'uniprot ', 'uni prot', 'uniprot\n', '', 1, None, '\ud800'])
def test_namespace_strict(value):
    invalid(lambda: ExternalIdentifierKey(value, 'P00533'), 'INVALID_ARGUMENT', ('namespace',))


def test_accession_normalization_preserves_case_and_whitespace():
    key = ExternalIdentifierKey('uniprot', ' e\u0301Ab.1 ')
    assert key.accession == ' éAb.1 '
    assert key != ExternalIdentifierKey('uniprot', 'éab.1')
    assert key == ExternalIdentifierKey('uniprot', ' éAb.1 ')


@pytest.mark.parametrize('value', ['', ' ', None, 1, '\udfff'])
def test_accession_invalid(value):
    invalid(lambda: ExternalIdentifierKey('uniprot', value), 'INVALID_ARGUMENT', ('accession',))


def test_unknown_fields_and_binding_errors():
    for cls, args in ((RegistryChange, ()), (RegistryWrite, ()), (RegistryRecord, ()),
                      (ExternalIdentifierKey, ()), (ExternalIdentifierMatch, ())):
        with pytest.raises(TypeError):
            cls(*args)
        with pytest.raises(TypeError):
            cls(unknown=True)


def test_concrete_model_is_not_flattened():
    protein = Protein(sf_id=eid(), entity_revision=1, label='payload', provenance=[alias().provenance])
    invalid(lambda: command(entity=protein), 'UNSUPPORTED_ENTITY_MODEL', ('entity',))
    assert protein.label == 'payload'
    invalid(lambda: command(entity=descriptor().model_dump()), 'INVALID_ARGUMENT', ('entity',))


@pytest.mark.parametrize('field,value', [('entity', None), ('aliases', None), ('aliases', {alias()}),
    ('aliases', iter(())), ('aliases', ({},)), ('change', {})])
def test_wrong_nested_types(field, value):
    error = invalid(lambda: replace(command(), **{field: value}), 'INVALID_ARGUMENT')
    assert error.field_path[0] == field


def test_nested_validation_issues_preserve_all_locations():
    broken = alias()
    # Deliberately corrupt otherwise valid frozen values to exercise revalidation.
    object.__setattr__(broken.provenance.source_record, 'version_status', 'INVALID')
    object.__setattr__(broken.provenance, 'importer_id', '')
    error = invalid(lambda: command(aliases=(broken,)), 'INVALID_ARGUMENT')
    assert (('aliases', 0, 'provenance', 'source_record', 'version_status'), 'literal_error') in error.validation_issues
    assert (('aliases', 0, 'provenance', 'importer_id'), 'value_error') in error.validation_issues
    assert error.field_path == error.validation_issues[0][0]


def test_invalid_nested_entity_revision_retains_validation_code():
    broken = descriptor()
    object.__setattr__(broken, 'entity_revision', True)
    error = invalid(lambda: command(entity=broken), 'INVALID_ARGUMENT', ('entity', 'entity_revision'))
    assert error.validation_issues == ((('entity', 'entity_revision'), 'int_type'),)


def test_scalar_and_timezone_ownership():
    class MutableZone(tzinfo):
        hours = 2
        def utcoffset(self, value):
            return timedelta(hours=self.hours)
        def dst(self, value):
            return timedelta(0)
    class Text(str):
        def __str__(self):
            return 'spoofed'
    class Stamp(datetime):
        def astimezone(self, tz=None):
            raise AssertionError('subclass behavior retained')
    zone = MutableZone()
    stamp = Stamp(2026, 1, 2, 2, tzinfo=zone)
    owned = replace(change(), actor=Text('actor'), reason=Text('reason'), recorded_at=stamp)
    zone.hours = 10
    assert owned.recorded_at == NOW
    assert type(owned.recorded_at) is datetime and owned.recorded_at.tzinfo is timezone.utc
    assert type(owned.actor) is str and owned.actor == 'actor'
    assert owned.reason == 'reason'
    key = ExternalIdentifierKey(Text('uniprot'), Text(' Ab '))
    assert type(key.accession) is str and key.accession == ' Ab '


@pytest.mark.parametrize('error', [EntityNotFoundError(eid()), RegistryRevisionNotFoundError(eid(), 2),
    EntityRevisionNotFoundError(eid(), 2), RevisionConflictError(eid(), 1, 2),
    OperationConflictError(change().operation_id), InvalidRegistryOperationError('INVALID_ARGUMENT', field_path=('entity',))])
def test_error_fields_read_only_and_tracebacks_work(error):
    for field in error._field_names:
        with pytest.raises(AttributeError):
            setattr(error, field, None)
        with pytest.raises(AttributeError):
            delattr(error, field)
    with pytest.raises(type(error)):
        raise error


def test_values_frozen_and_outside_wire_contract():
    values = (change(), command(), RegistryRecord(1, descriptor(), (), change()),
              ExternalIdentifierKey('uniprot', 'P00533'), ExternalIdentifierMatch(descriptor(), 1, alias()))
    for value in values:
        assert not hasattr(value, 'to_canonical_bytes')
        assert not hasattr(value, '__dict__')
        field = next(iter(value.__dataclass_fields__))
        with pytest.raises((FrozenInstanceError, AttributeError)):
            setattr(value, field, None)


def test_invalid_values_retain_available_entity_id():
    for build in (lambda: replace(command(), expected_registry_revision=-1),
                  lambda: replace(command(), aliases=(None,)),
                  lambda: replace(command(), change=None)):
        assert invalid(build, 'INVALID_ARGUMENT').entity_id == eid()
    broken = alias()
    object.__setattr__(broken, 'lifecycle', 'invalid')
    assert invalid(lambda: command(aliases=(broken,)), 'INVALID_ARGUMENT').entity_id == eid()
    entity = descriptor()
    object.__setattr__(entity, 'entity_revision', True)
    assert invalid(lambda: command(entity=entity), 'INVALID_ARGUMENT').entity_id == eid()


def _rich_payloads():
    from simforge.knowledge.domain import (ProteinSequence, ChemicalIdentity, ChemicalSpecies,
        ExperimentalStructure, Claim, Evidence)
    provenance = alias().provenance
    base = dict(sf_id=eid(), entity_revision=1)
    def ref(n, kind):
        return dict(entity_id=eid(n), entity_type=kind, entity_revision=1)
    representation = dict(representation_type='chem:smiles', value='CCO', specificity='PARTIAL', provenance=provenance)
    artifact = dict(artifact_id=eid(20), sha256='a'*64)
    return [
        Protein(**base, label='source payload', provenance=[provenance]),
        ProteinSequence(**base, sequence='ACD', provenance=[provenance]),
        ChemicalIdentity(**base, representations=[representation], stereochemistry_status='UNSPECIFIED', isotope_status='UNSPECIFIED'),
        ChemicalSpecies(**base, chemical_identity=ref(2, 'chemical_identity'), representations=[representation],
            formal_charge=None, protonation_state=None, tautomer_state=None, description_status='PARTIAL'),
        ExperimentalStructure(**base, origin='experimental',
            source_identifier=alias(kind='experimental_structure', mapping_status='CONFIRMED_MAPPING', evidence_refs=[ref(11, 'evidence')]),
            method=dict(category='X_RAY_DIFFRACTION', raw_method='X-RAY DIFFRACTION', declared_origin='EXPERIMENTAL',
                classification_basis='SOURCE_DECLARATION', classification_claim=ref(12, 'claim'), provenance=provenance),
            artifact=artifact, evidence_refs=[ref(11, 'evidence')]),
        Claim(**base, subject=ref(2, 'protein'), predicate='science:property',
            object=dict(kind='literal', value=None), assertion_provenance=[provenance]),
        Evidence(**base, claim_links=[dict(claim=ref(2, 'claim'), relation='SUPPORTS')],
            evidence_class='computational_prediction', provenance=provenance, raw_artifacts=[artifact]),
    ]


@pytest.mark.parametrize('payload', _rich_payloads(), ids=lambda x: type(x).__name__)
def test_every_frozen_scientific_subclass_rejected_without_projection(registry, payload):
    before = payload.to_canonical_bytes()
    invalid(lambda: command(entity=payload), 'UNSUPPORTED_ENTITY_MODEL', ('entity',))
    forged = command()
    object.__setattr__(forged, 'entity', payload)
    error = invalid(lambda: registry.write(forged), 'UNSUPPORTED_ENTITY_MODEL', ('entity',))
    assert error.entity_id == eid() and payload.to_canonical_bytes() == before
    with pytest.raises(EntityNotFoundError):
        registry.get(eid())


def test_closed_failure_codes_are_literal_contract():
    from typing import get_type_hints, get_args
    assert set(get_args(get_type_hints(InvalidRegistryOperationError)['code'])) == {
        'INVALID_ARGUMENT', 'UNSUPPORTED_ENTITY_MODEL', 'ENTITY_TYPE_CHANGE',
        'ENTITY_REVISION_SEQUENCE', 'ENTITY_REVISION_REWRITE', 'ALIAS_TARGET_MISMATCH',
        'ALIAS_TARGET_REVISION_MISSING', 'DUPLICATE_ALIAS_DECLARATION'}


def test_match_validation_retains_available_id():
    error = invalid(lambda: ExternalIdentifierMatch(descriptor(), 1, None), 'INVALID_ARGUMENT', ('alias',))
    assert error.entity_id == eid()


def test_malformed_timezone_is_structured_invalid_argument():
    class InvalidZone(tzinfo):
        def utcoffset(self, value):
            return 'not an offset'
    invalid(lambda: replace(change(), recorded_at=NOW.replace(tzinfo=InvalidZone())),
            'INVALID_ARGUMENT', ('recorded_at',))


@pytest.mark.parametrize('boundary', ['write', 'record', 'match', 'operation'])
@pytest.mark.parametrize('location', ['valid_from', 'valid_until', 'imported_at'])
def test_malformed_nested_timestamp_is_structured(registry, boundary, location):
    class InvalidOffset(tzinfo):
        def utcoffset(self, value):
            return 'invalid offset'
    declaration = alias()
    cmd = command(aliases=(declaration,))
    if boundary == 'operation':
        declaration = cmd.aliases[0]
    owner = declaration.provenance if location == 'imported_at' else declaration
    object.__setattr__(owner, location, NOW.replace(tzinfo=InvalidOffset()))
    if boundary == 'write':
        call = lambda: command(aliases=(declaration,))
    elif boundary == 'record':
        call = lambda: RegistryRecord(1, descriptor(), (declaration,), change())
    elif boundary == 'match':
        call = lambda: ExternalIdentifierMatch(descriptor(), 1, declaration)
    else:
        call = lambda: registry.write(cmd)
    prefix = ('alias',) if boundary == 'match' else ('aliases', 0)
    path = prefix + (('provenance', 'imported_at') if location == 'imported_at' else (location,))
    error = invalid(call, 'INVALID_ARGUMENT', path)
    assert error.entity_id == eid()
    with pytest.raises(EntityNotFoundError):
        registry.get(eid())
    assert registry.write(command()).registry_revision == 1


@pytest.mark.parametrize('model', ['Entity', 'ExternalIdentifier'])
def test_valid_domain_input_internal_typeerror_is_not_reclassified(registry, monkeypatch, model):
    from simforge.knowledge import domain
    cmd = command(aliases=(alias(),))
    def programming_defect(*args, **kwargs):
        raise TypeError('internal validator defect on valid input')
    monkeypatch.setattr(getattr(domain, model), 'model_validate', programming_defect)
    with pytest.raises(TypeError, match='internal validator defect'):
        registry.write(cmd)


@pytest.mark.parametrize('factory', [lambda: EntityNotFoundError(eid()),
    lambda: RegistryRevisionNotFoundError(eid(), 2), lambda: EntityRevisionNotFoundError(eid(), 2),
    lambda: RevisionConflictError(eid(), 1, 2), lambda: OperationConflictError(change().operation_id),
    lambda: InvalidRegistryOperationError('INVALID_ARGUMENT', entity_id=eid(), field_path=('aliases', 0),
        validation_issues=((('aliases', 0, 'lifecycle'), 'literal_error'),))])
def test_every_error_field_is_nonshadowable(factory):
    error = factory()
    for field in type(error)._field_names:
        original = getattr(error, field)
        with pytest.raises(AttributeError):
            setattr(error, field, None)
        with pytest.raises(AttributeError):
            delattr(error, field)
        error.__dict__[field] = object()
        assert getattr(error, field) == original
    error.__dict__['_details'] = ('shadowed',)
    error.__dict__['_field_names'] = ()
    for field in type(error)._field_names:
        with pytest.raises(AttributeError):
            setattr(error, field, None)
    with pytest.raises(type(error)):
        raise error


def test_error_constructor_owns_nested_inputs():
    target = eid()
    path = ['aliases', 0]
    issue_path = ['aliases', 0, 'lifecycle']
    issues = [[issue_path, 'literal_error']]
    error = InvalidRegistryOperationError('INVALID_ARGUMENT', entity_id=target,
        field_path=path, validation_issues=issues)
    object.__setattr__(target, 'value', eid(2).value)
    path.clear(); issue_path.clear(); issues.clear()
    assert error.entity_id == eid()
    assert error.field_path == ('aliases', 0)
    assert error.validation_issues == ((('aliases', 0, 'lifecycle'), 'literal_error'),)
    operation = change().operation_id
    conflict = OperationConflictError(operation)
    object.__setattr__(operation, 'int', 0)
    assert conflict.operation_id == change().operation_id
