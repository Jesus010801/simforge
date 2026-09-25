import pytest
from simforge.knowledge.domain import Claim, Evidence
from test_identity import base, ref, provenance, artifact


def claim():
    return dict(**base(70), subject=ref(), predicate='science:property', object=dict(kind='literal', value=None), assertion_provenance=[provenance()])


def evidence():
    return dict(**base(71), claim_links=[dict(claim=ref(70, 'claim'), relation='SUPPORTS')], evidence_class='curated', provenance=provenance(), raw_artifacts=[artifact()])


def test_claim_evidence_and_explicit_null():
    a, b = Claim(**claim()), Evidence(**evidence())
    assert not isinstance(a, Evidence)
    assert not isinstance(b, Claim)
    assert a.object.value.kind == 'null'
    assert b.claim_links[0].claim.entity_id == a.sf_id
    assert not hasattr(b, 'confidence')
    data = claim()
    del data['object']['value']
    with pytest.raises(ValueError):
        Claim(**data)


@pytest.mark.parametrize('kind', ['experimental', 'curated', 'computational_prediction', 'computational_simulation', 'inferred', 'user_supplied'])
def test_evidence_classes_are_preserved(kind):
    obj = Evidence(**(evidence() | dict(evidence_class=kind)))
    assert Evidence.from_canonical_bytes(obj.to_canonical_bytes()).evidence_class == kind


@pytest.mark.parametrize('relation', ['SUPPORTS', 'REFUTES', 'CONTEXT'])
def test_evidence_relation(relation):
    data = evidence()
    data['claim_links'][0]['relation'] = relation
    assert Evidence(**data).claim_links[0].relation == relation


def test_duplicate_conflicting_links_and_raw_hash_rejected():
    data = evidence()
    data['claim_links'].append(dict(claim=ref(70, 'claim'), relation='REFUTES'))
    with pytest.raises(ValueError):
        Evidence(**data)
    data = evidence()
    data['raw_artifacts'][0]['sha256'] = 'b'*64
    with pytest.raises(ValueError):
        Evidence(**data)


def test_contradictions_survive_without_resolution():
    data = claim()
    data['contradicts'] = [ref(72, 'claim')]
    result = Claim(**data)
    assert Claim.from_canonical_bytes(result.to_canonical_bytes()) == result
    data['contradicts'] = [ref(70, 'claim')]
    with pytest.raises(ValueError):
        Claim(**data)
