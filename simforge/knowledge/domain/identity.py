"""Opaque internal identity, external mappings and source attribution."""
from __future__ import annotations
import re
from typing import Literal
from uuid import UUID, RFC_4122
from pydantic import model_validator
from .serialization import _DomainModel, _Text, _Token, _PositiveInt, _UTC, _SHA256, _Items

_EntityType = Literal['protein', 'protein_sequence', 'chemical_identity', 'chemical_species', 'experimental_structure', 'claim', 'evidence']
_Lifecycle = Literal['ACTIVE', 'DEPRECATED', 'RETRACTED']


class EntityId(_DomainModel):
    """An opaque supplied UUIDv4; no kind, provider or generation behavior."""
    _scalar = True
    value: str

    @model_validator(mode='before')
    @classmethod
    def _identity(cls, value):
        if isinstance(value, str):
            value = {'value': value}
        if isinstance(value, dict):
            raw = value.get('value')
            if not isinstance(raw, str) or not re.fullmatch(r'sf:[0-9a-f]{8}(?:-[0-9a-f]{4}){3}-[0-9a-f]{12}', raw):
                raise ValueError('kind-free canonical sf:UUIDv4 required')
            uid = UUID(raw[3:])
            if uid.version != 4 or uid.variant != RFC_4122:
                raise ValueError('UUIDv4 required')
        return value


class _Ref(_DomainModel):

    def _identity_roles(self) -> tuple[tuple[str, str], ...]:
        return ((self.entity_id.value, self.entity_type),)
    entity_id: EntityId
    entity_revision: _PositiveInt
    entity_type: _EntityType


class _ClaimRef(_Ref):
    entity_type: Literal['claim']


class _EvidenceRef(_Ref):
    entity_type: Literal['evidence']


class _ChemicalRef(_Ref):
    entity_type: Literal['chemical_identity']


def _consistent_types(refs):
    types = {}
    for ref in refs:
        key = ref.entity_id.value
        if key in types and types[key] != ref.entity_type:
            raise ValueError('same internal ID has conflicting entity types')
        types[key] = ref.entity_type


class _SourceRecord(_DomainModel):
    source_id: _Token
    record_id: _Text
    source_release: _Text | None = None
    record_revision: _Text | None = None
    version_status: Literal['PINNED', 'UNVERSIONED']

    @model_validator(mode='after')
    def _pin(self):
        pinned = self.source_release is not None or self.record_revision is not None
        if pinned != (self.version_status == 'PINNED'):
            raise ValueError('source version status disagrees with source versions')
        return self


class _ImportProvenance(_DomainModel):
    source_record: _SourceRecord
    raw_payload_sha256: _SHA256
    importer_id: _Token
    importer_version: _Token
    imported_at: _UTC
    extraction_method: _Token


class _RuleRef(_DomainModel):
    rule_id: _Token
    version: _Token
    definition_sha256: _SHA256


class ExternalIdentifier(_DomainModel):
    """A declared mapping, never an equality operation or primary key."""
    entity: _Ref
    namespace: _Token
    accession: _Text
    accession_version: _Text | None = None
    source_record_revision: _Text | None = None
    valid_from: _UTC | None = None
    valid_until: _UTC | None = None
    mapping_status: Literal['CANDIDATE_MAPPING', 'CONFIRMED_MAPPING']
    lifecycle: _Lifecycle
    evidence_refs: _Items[_EvidenceRef] = ()
    confirmation_rule: _RuleRef | None = None
    provenance: _ImportProvenance
    _unordered_fields = ('evidence_refs',)

    @model_validator(mode='after')
    def _mapping(self):
        if not re.fullmatch(r'[a-z][a-z0-9_.-]*', self.namespace):
            raise ValueError('lowercase namespace token required')
        if self.valid_from is not None and self.valid_until is not None and self.valid_until <= self.valid_from:
            raise ValueError('invalid validity interval')
        if self.mapping_status == 'CONFIRMED_MAPPING' and not (self.evidence_refs or self.confirmation_rule):
            raise ValueError('confirmed mapping needs evidence or a pinned rule')
        source = self.provenance.source_record
        if (source.record_id == self.accession and source.record_revision is not None
                and self.source_record_revision is not None and source.record_revision != self.source_record_revision):
            raise ValueError('conflicting source record revisions')
        _consistent_types((self.entity, *self.evidence_refs))
        return self
