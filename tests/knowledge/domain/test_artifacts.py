import pytest
from simforge.knowledge.domain import ArtifactRef, PolicyRef, DecisionContext
from test_identity import artifact, policy, context


@pytest.mark.parametrize('digest', ['', 'a'*63, 'a'*65, 'g'*64, ' '+ 'a'*64, 'sha256:'+'a'*64, b'a'*64, 1, None])
def test_strict_sha256(digest):
    with pytest.raises(ValueError):
        ArtifactRef(**artifact(digest=digest))
    with pytest.raises(ValueError):
        PolicyRef(**(policy() | dict(definition_sha256=digest)))
    with pytest.raises(ValueError):
        DecisionContext(**(context() | dict(request_sha256=digest)))


def test_digest_case_and_zero_semantics():
    assert ArtifactRef(**artifact(digest='Ab'*32)).sha256.value == 'ab'*32
    assert ArtifactRef(**artifact(digest='0'*64)).sha256.value == '0'*64


@pytest.mark.parametrize('size', [-1, True, 1.0, '1'])
def test_invalid_size(size):
    with pytest.raises(ValueError):
        ArtifactRef(**artifact(), size_bytes=size)


def test_absent_digest_and_unknown_size():
    data = artifact()
    del data['sha256']
    with pytest.raises(ValueError):
        ArtifactRef(**data)
    assert ArtifactRef(**artifact()).size_bytes is None
