"""Content-pinned references; no filesystem access or registry."""
from typing import Annotated
from pydantic import Field
from .identity import EntityId
from .serialization import _DomainModel, _SHA256, _Token


class ArtifactRef(_DomainModel):

    def _identity_roles(self) -> tuple[tuple[str, str], ...]:
        return ((self.artifact_id.value, 'artifact'),)
    artifact_id: EntityId
    sha256: _SHA256
    media_type: _Token | None = None
    size_bytes: Annotated[int, Field(ge=0)] | None = None


def _artifact_consistency(refs):
    seen = {}
    for ref in refs:
        key = ref.artifact_id.value
        if key in seen and seen[key] != ref.sha256:
            raise ValueError('conflicting digests for artifact ID')
        seen[key] = ref.sha256


def _artifact_key(ref):
    return ref.artifact_id.value, ref.sha256.value
