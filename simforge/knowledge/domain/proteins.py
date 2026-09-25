"""Protein concepts and provider-independent sequence fingerprints."""
import hashlib
from typing import Annotated, Literal
from pydantic import Field, field_validator
from .entities import Entity
from .identity import _ImportProvenance
from .serialization import _Items, _Text, _SHA256


class Protein(Entity):
    entity_type: Literal['protein'] = 'protein'
    label: _Text | None = None
    provenance: Annotated[_Items[_ImportProvenance], Field(min_length=1)]
    _unordered_fields = ('provenance',)


class ProteinSequence(Entity):
    entity_type: Literal['protein_sequence'] = 'protein_sequence'
    sequence: str
    canonicalization_version: Literal['sski-protein-sequence-v1'] = 'sski-protein-sequence-v1'
    provenance: Annotated[_Items[_ImportProvenance], Field(min_length=1)]
    _unordered_fields = ('provenance',)

    @field_validator('sequence', mode='before')
    @classmethod
    def _sequence(cls, value):
        if not isinstance(value, str) or not value.isascii():
            raise ValueError('ASCII protein sequence required')
        value = value.translate(str.maketrans('', '', ' \t\n\r\v\f')).upper()
        if not value or set(value) - set('ACDEFGHIKLMNPQRSTVWYUOBZJX'):
            raise ValueError('invalid protein sequence alphabet')
        return value

    @property
    def length(self):
        return len(self.sequence)

    @property
    def contains_ambiguous_symbols(self):
        return bool(set(self.sequence) & set('BZJX'))

    @property
    def sequence_sha256(self):
        return _SHA256(value=hashlib.sha256(b'sski-protein-sequence-v1\x00' + self.sequence.encode('ascii')).hexdigest())
