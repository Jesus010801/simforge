import pytest
from simforge.knowledge.domain import Protein, ProteinSequence
from test_identity import base, provenance


def seq(value, n=1):
    return ProteinSequence(**base(n), sequence=value, provenance=[provenance()])


def test_normalization_and_ambiguity_preserved():
    a = seq(' \tac\nd\re\vf\fuobzjx ')
    assert a.sequence == 'ACDEFUOBZJX'
    assert a.length == 11
    assert a.contains_ambiguous_symbols
    assert a.sequence_sha256 == seq('ACDEFUOBZJX', 2).sequence_sha256
    assert a.sf_id != seq('ACDEFUOBZJX', 2).sf_id


@pytest.mark.parametrize('value', ['', ' \n', 'A-C', 'AC*', '*AC', 'A*C', 'AéC', 'A\u00a0C', 'A1C', 'AßC'])
def test_rejected_sequences(value):
    with pytest.raises(ValueError):
        seq(value)


def test_protein_is_not_sequence():
    protein = Protein(**base(), provenance=[provenance()])
    assert not hasattr(protein, 'sequence')
    assert not isinstance(seq('ACD'), Protein)


@pytest.mark.parametrize('sequence,expected', [
    ('ACD', 'a1f16901979311ac97bd78aac7809fbd630fa14723749b2dd993bb4fc413f985'),
    ('ACDEFGHIKLMNPQRSTVWYUOBZJX', '21775e157676b6cd3fb9f39624a832942e76415e820ba3ebe1b9192c1bbd9e8a'),
    ('X', '25dc25b0aa6137291b0936b81fd1c21857a9bcccc44fc66aa4aa088ca7b1c24a'),
])
def test_sequence_golden_vectors(sequence, expected):
    assert seq(sequence).sequence_sha256.value == expected
