"""Manifest values only: no availability, filesystem or compilation guarantees."""
import re
from typing import Literal
from pydantic import field_validator, model_validator
from .identity import EntityId
from .artifacts import ArtifactRef, _artifact_consistency, _artifact_key
from .decisions import PolicyRef, DecisionContext, _require_members
from .serialization import _DomainModel, _SHA256, _Token, _Text, _UTC, _Items


def _policy_key(ref):
    return ref.policy_id.value, ref.policy_version, ref.definition_sha256.value


class _SnapshotRef(_DomainModel):

    def _identity_roles(self) -> tuple[tuple[str, str], ...]:
        return ((self.snapshot_id.value, 'snapshot'),)
    snapshot_id: EntityId
    sha256: _SHA256


class _PolicySetRef(_DomainModel):

    def _identity_roles(self) -> tuple[tuple[str, str], ...]:
        return ((self.policy_set_id.value, 'policy_set'),)
    policy_set_id: EntityId
    sha256: _SHA256
    policies: _Items[PolicyRef] = ()
    _unordered_fields = ('policies',)

    @model_validator(mode='after')
    def _unique(self):
        seen = {}
        for p in self.policies:
            key = p.policy_id.value, p.policy_version
            if key in seen and seen[key] != p.definition_sha256:
                raise ValueError('conflicting policy digests')
            seen[key] = p.definition_sha256
        return self


class _ManifestArtifact(_DomainModel):
    artifact: ArtifactRef
    role: _Token
    path: _Text

    @field_validator('path')
    @classmethod
    def _path(cls, value):
        if (value.startswith('/') or '\\' in value or re.match(r'^[A-Za-z]:', value)
                or any(s in ('', '.', '..') for s in value.split('/'))
                or any(ord(c) < 32 or ord(c) == 127 for c in value)):
            raise ValueError('canonical relative POSIX path required')
        return value


class KnowledgeBundleManifest(_DomainModel):

    def _identity_roles(self) -> tuple[tuple[str, str], ...]:
        return ((self.bundle_id.value, 'bundle'),)
    bundle_id: EntityId
    schema_version: Literal['0.1.0'] = '0.1.0'
    created_at: _UTC
    immutable: Literal[True] = True
    knowledge_snapshot: _SnapshotRef
    policy_set: _PolicySetRef
    normalized_input_sha256: _SHA256
    decision_context: DecisionContext
    artifacts: _Items[_ManifestArtifact] = ()
    _unordered_fields = ('artifacts',)

    @field_validator('immutable', mode='before')
    @classmethod
    def _true(cls, value):
        if value is not True:
            raise ValueError('immutable must be true')
        return value

    @model_validator(mode='after')
    def _manifest(self):
        context = self.decision_context
        if self.normalized_input_sha256 != context.request_sha256:
            raise ValueError('request digests disagree')
        if context.has_mandatory_abort:
            raise ValueError('mandatory ABORT decision prevents manifest construction')
        _require_members((d.policy for d in context.decisions), self.policy_set.policies, _policy_key)
        artifacts = tuple(x.artifact for x in self.artifacts)
        _artifact_consistency(artifacts)
        _require_members(context.artifact_refs, artifacts, _artifact_key)
        paths = [a.path for a in self.artifacts]
        if len(set(paths)) != len(paths):
            raise ValueError('duplicate manifest path')
        return self
