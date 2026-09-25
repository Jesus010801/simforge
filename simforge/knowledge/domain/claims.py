"""Scientific propositions, separate from their evidential support."""
from typing import Annotated, Literal
from pydantic import Field, model_validator
from .entities import Entity
from .identity import _Ref, _ClaimRef, _ImportProvenance, _consistent_types
from .serialization import _DomainModel, _Items, _Symbol, _Structured, _Object


class _EntityObject(_DomainModel):
    kind: Literal['entity']
    value: _Ref


class _LiteralObject(_DomainModel):
    kind: Literal['literal']
    value: _Structured


class Claim(Entity):
    entity_type: Literal['claim'] = 'claim'
    subject: _Ref
    predicate: _Symbol
    object: _EntityObject | _LiteralObject
    qualifiers: _Object = Field(default_factory=dict)
    scope: _Object = Field(default_factory=dict)
    assertion_provenance: Annotated[_Items[_ImportProvenance], Field(min_length=1)]
    contradicts: _Items[_ClaimRef] = ()
    _unordered_fields = ('assertion_provenance', 'contradicts')

    @model_validator(mode='after')
    def _relationships(self):
        own = _Ref(entity_id=self.sf_id, entity_revision=self.entity_revision, entity_type='claim')
        if any(r.entity_id == self.sf_id and r.entity_revision == self.entity_revision for r in self.contradicts):
            raise ValueError('claim cannot contradict itself at the same revision')
        refs = (own, self.subject, *self.contradicts)
        if isinstance(self.object, _EntityObject):
            refs += (self.object.value,)
        _consistent_types(refs)
        return self
