import pytest
from simforge.knowledge.domain import KnowledgeBundleManifest
from test_identity import manifest, artifact, decision, policy


def with_artifact(path='artifacts/structures/a.cif'):
    data = manifest()
    data['artifacts'] = [dict(artifact=artifact(), role='structure', path=path)]
    data['decision_context']['artifact_refs'] = [artifact()]
    return data


@pytest.mark.parametrize('path', ['/a', 'C:/a', 'C:a', 'a\\b', 'a//b', './a', 'a/../b', 'a/./b', '../a', 'a/', 'a/\x00b', ''])
def test_path_rejections(path):
    with pytest.raises(ValueError):
        KnowledgeBundleManifest(**with_artifact(path))


def test_manifest_consistency_and_roundtrip():
    data = with_artifact()
    value = KnowledgeBundleManifest(**data)
    assert KnowledgeBundleManifest.from_canonical_bytes(value.to_canonical_bytes()) == value
    data['artifacts'].append(dict(artifact=artifact(), role='copy', path='artifacts/copy.cif'))
    assert len(KnowledgeBundleManifest(**data).artifacts) == 2


@pytest.mark.parametrize('mismatch', ['request', 'policy', 'policy_hash', 'artifact', 'artifact_hash', 'duplicate_path', 'conflicting_artifact', 'conflicting_policy'])
def test_manifest_invariant_mismatch(mismatch):
    data = with_artifact()
    if mismatch == 'request': data['normalized_input_sha256'] = 'b'*64
    elif mismatch == 'policy': data['policy_set']['policies'] = []
    elif mismatch == 'policy_hash': data['policy_set']['policies'][0]['definition_sha256'] = 'b'*64
    elif mismatch == 'artifact': data['artifacts'] = []
    elif mismatch == 'artifact_hash': data['artifacts'][0]['artifact']['sha256'] = 'b'*64
    elif mismatch == 'duplicate_path': data['artifacts'].append(dict(artifact=artifact(21), role='other', path=data['artifacts'][0]['path']))
    elif mismatch == 'conflicting_artifact': data['artifacts'].append(dict(artifact=artifact(digest='b'*64), role='other', path='artifacts/other'))
    elif mismatch == 'conflicting_policy': data['policy_set']['policies'].append(policy() | dict(definition_sha256='b'*64))
    with pytest.raises(ValueError):
        KnowledgeBundleManifest(**data)


def test_optional_abort_retained_mandatory_abort_rejected():
    data = manifest()
    data['decision_context']['decisions'] = [decision('ABORT')]
    with pytest.raises(ValueError):
        KnowledgeBundleManifest(**data)
    data['decision_context']['decisions'][0]['mandatory'] = False
    assert KnowledgeBundleManifest(**data).decision_context.decisions[0].status == 'ABORT'


@pytest.mark.parametrize('where', ['snapshot', 'policy_set', 'input'])
def test_all_manifest_digests_are_strict(where):
    data = manifest()
    if where == 'snapshot': data['knowledge_snapshot']['sha256'] = 'bad'
    elif where == 'policy_set': data['policy_set']['sha256'] = 'bad'
    else: data['normalized_input_sha256'] = 'bad'
    with pytest.raises(ValueError):
        KnowledgeBundleManifest(**data)


@pytest.mark.parametrize('value', [False, 1, 'true'])
def test_immutable_is_literal_boolean(value):
    with pytest.raises(ValueError):
        KnowledgeBundleManifest(**(manifest() | dict(immutable=value)))


@pytest.mark.parametrize('role', ['snapshot', 'policy_set', 'context', 'policy', 'artifact'])
def test_generic_ids_cannot_acquire_conflicting_roles_inside_manifest(role):
    data = with_artifact()
    same = data['bundle_id']
    if role == 'snapshot': data['knowledge_snapshot']['snapshot_id'] = same
    elif role == 'policy_set': data['policy_set']['policy_set_id'] = same
    elif role == 'context': data['decision_context']['context_id'] = same
    elif role == 'policy':
        data['policy_set']['policies'][0]['policy_id'] = same
        data['decision_context']['decisions'][0]['policy']['policy_id'] = same
    else:
        data['artifacts'][0]['artifact']['artifact_id'] = same
        data['decision_context']['artifact_refs'][0]['artifact_id'] = same
    with pytest.raises(ValueError, match='conflicting'):
        KnowledgeBundleManifest(**data)
