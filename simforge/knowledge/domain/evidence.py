"""Source-attributed support, refutation and context for claim revisions."""
from typing import Annotated, Literal
from pydantic import Field, model_validator
from .entities import Entity
from .identity import _ClaimRef, _ImportProvenance, _consistent_types, _Ref
from .artifacts import ArtifactRef, _artifact_consistency
from .serialization import _DomainModel, _Items, _Text


class _EvidenceLink(_DomainModel):
    claim: _ClaimRef
    relation: Literal['SUPPORTS', 'REFUTES', 'CONTEXT']


class Evidence(Entity):
    entity_type: Literal['evidence'] = 'evidence'
    claim_links: Annotated[_Items[_EvidenceLink], Field(min_length=1)]
    evidence_class: Literal['experimental', 'curated', 'computational_prediction', 'computational_simulation', 'inferred', 'user_supplied']
    provenance: _ImportProvenance
    raw_artifacts: Annotated[_Items[ArtifactRef], Field(min_length=1)]
    method: _Text | None = None
    references: _Items[_Text] = ()
    limitations: _Items[_Text] = ()
    _unordered_fields = ('claim_links', 'raw_artifacts', 'references')

    @model_validator(mode='after')
    def _links(self):
        keys = [(x.claim.entity_id.value, x.claim.entity_revision) for x in self.claim_links]
        if len(set(keys)) != len(keys):
            raise ValueError('duplicate or conflicting claim links')
        if not any(a.sha256 == self.provenance.raw_payload_sha256 for a in self.raw_artifacts):
            raise ValueError('raw provenance digest is absent from artifacts')
        _artifact_consistency(self.raw_artifacts)
        _consistent_types((_Ref(entity_id=self.sf_id, entity_revision=self.entity_revision, entity_type='evidence'), *(x.claim for x in self.claim_links)))
        return self
