"""Identity regressions and deterministic fixture data shared by domain tests."""
from datetime import datetime, timezone
import pytest
from simforge.knowledge.domain import EntityId, ExternalIdentifier

DIGEST = 'a' * 64
NOW = datetime(2026, 1, 2, 3, 4, 5, tzinfo=timezone.utc)


def uid(n=1):
    return f'sf:00000000-0000-4000-8000-{n:012x}'


def ref(n=1, kind='protein', revision=1):
    return dict(entity_id=uid(n), entity_revision=revision, entity_type=kind)


def provenance():
    return dict(source_record=dict(source_id='fixture', record_id='record-1', source_release='2026.1', record_revision='2', version_status='PINNED'), raw_payload_sha256=DIGEST, importer_id='fixture-importer', importer_version='1', imported_at=NOW, extraction_method='fixture')


def external():
    return dict(entity=ref(), namespace='uniprot', accession='P00533', accession_version='7', source_record_revision='2', mapping_status='CANDIDATE_MAPPING', lifecycle='ACTIVE', provenance=provenance())


def base(n=1):
    return dict(sf_id=uid(n), entity_revision=1)


def artifact(n=20, digest=DIGEST):
    return dict(artifact_id=uid(n), sha256=digest)


def policy():
    return dict(policy_id=uid(30), policy_version='1', definition_sha256=DIGEST)


def operation():
    return dict(operation_key='protonation', target_stage='preparation', action='science:protonation', method_id='science:method', method_version='1', method_sha256=DIGEST, required_inputs=[dict(name='molecule', kind='entity', value=ref())], expected_outputs=[dict(name='species', scientific_type='science:species')], success_criteria=['charge_consistency'], on_failure='ABORT')


def decision(status='RESOLVED'):
    d = dict(decision_id=uid(40), decision_revision=1, question='Select input', status=status, policy=policy(), rationale='Pinned fixture policy', entity_inputs=[ref()])
    if status in ('RESOLVED', 'RESOLVED_WITH_WARNING'):
        d['result'] = dict(selected_entities=[ref()])
    if status == 'RESOLVED_WITH_WARNING':
        d['warnings'] = ['Limited source coverage']
    if status == 'DEFERRED':
        d['deferred_operation'] = operation()
    return d


def context():
    return dict(context_id=uid(50), context_revision=1, request_sha256=DIGEST, entity_refs=[ref()], decisions=[decision()])


def manifest():
    return dict(bundle_id=uid(60), created_at=NOW, knowledge_snapshot=dict(snapshot_id=uid(61), sha256=DIGEST), policy_set=dict(policy_set_id=uid(62), sha256=DIGEST, policies=[policy()]), normalized_input_sha256=DIGEST, decision_context=context())


@pytest.mark.parametrize('value', ['P00533', 'uniprot:P00533', 'sf:protein:P00533', 'sf:protein:'+uid()[3:], uid().replace('-4000-', '-7000-'), uid().replace('-8000-', '-0000-'), uid().upper(), '01ARZ3NDEKTSV4RRFFQ69G5FAV', 1, None])
def test_invalid_internal_identity(value):
    with pytest.raises(ValueError):
        EntityId.model_validate(value)


def test_kind_free_identity_and_round_trip():
    value = EntityId.model_validate(uid())
    assert value.value == uid()
    assert EntityId.from_canonical_bytes(value.to_canonical_bytes()) == value
    with pytest.raises(ValueError):
        EntityId()


def test_external_identity_and_source_revision_are_separate():
    data = external()
    a = ExternalIdentifier(**data)
    data['entity'] = ref(2, revision=9)
    b = ExternalIdentifier(**data)
    assert a.accession == b.accession
    assert a.entity.entity_id != b.entity.entity_id
    assert b.entity.entity_revision == 9
    assert b.source_record_revision == '2'
    assert b.accession_version == '7'
    with pytest.raises(ValueError):
        EntityId.model_validate(a)


def test_confirmation_is_explicit():
    data = external()
    data['mapping_status'] = 'CONFIRMED_MAPPING'
    with pytest.raises(ValueError):
        ExternalIdentifier(**data)
    data['evidence_refs'] = [ref(5, 'evidence')]
    assert ExternalIdentifier(**data).mapping_status == 'CONFIRMED_MAPPING'
    data['mapping_status'] = 'CANDIDATE_MAPPING'
    assert ExternalIdentifier(**data).mapping_status == 'CANDIDATE_MAPPING'
    data['mapping_status'] = 'CONFIRMED_MAPPING'
    data['evidence_refs'] = []
    data['confirmation_rule'] = dict(rule_id='rule', version='1', definition_sha256=DIGEST)
    assert ExternalIdentifier(**data).confirmation_rule.version == '1'


@pytest.mark.parametrize('change', [dict(namespace='UniProt'), dict(valid_from=NOW, valid_until=NOW), dict(valid_from=NOW.replace(tzinfo=None)), dict(lifecycle='DELETED')])
def test_mapping_rejects_invalid_metadata(change):
    with pytest.raises(ValueError):
        ExternalIdentifier(**(external() | change))


def test_source_revision_consistency_and_unversioned_content_pin():
    data = external()
    data['accession'] = 'record-1'
    data['source_record_revision'] = '3'
    with pytest.raises(ValueError):
        ExternalIdentifier(**data)
    data = external()
    source = data['provenance']['source_record']
    source.update(version_status='UNVERSIONED', source_release=None, record_revision=None)
    assert ExternalIdentifier(**data).provenance.raw_payload_sha256.value == DIGEST
    source['source_release'] = '1'
    with pytest.raises(ValueError):
        ExternalIdentifier(**data)
