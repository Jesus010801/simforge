"""Backend-independent, per-call linearizable registry contract."""
from typing import Protocol as _Protocol
from simforge.knowledge.domain import Entity as _Entity, EntityId as _EntityId
from .records import RegistryWrite as _Write, RegistryRecord as _Record
from .records import ExternalIdentifierKey as _Key, ExternalIdentifierMatch as _Match

__all__ = ['EntityRegistry']


class EntityRegistry(_Protocol):
    """Atomic writes and immutable reads over one logical registry.

    Revisions are per EntityId. Successful operation IDs span the entire logical
    registry and exact retries return the original result before checking heads.
    Separate calls never imply a cross-entity transaction or a pinned snapshot.
    """

    def write(self, command: _Write) -> _Record:
        """Validate, replay/check operation ID, CAS, then atomically commit."""
        ...

    def get(self, entity_id: _EntityId, *, registry_revision: int | None = None) -> _Record:
        """Return current head or exactly the requested positive revision."""
        ...

    def get_entity_revision(self, entity_id: _EntityId, entity_revision: int) -> _Entity:
        """Return the exact registered descriptor revision, never fall forward."""
        ...

    def history(self, entity_id: _EntityId) -> tuple[_Record, ...]:
        """Return retained snapshots in ascending registry revision order."""
        ...

    def lookup_external(self, key: _Key) -> tuple[_Match, ...]:
        """Current declarations ordered by ID, then canonical alias bytes.

        Includes every lifecycle/status/validity; performs no resolution.
        """
        ...
