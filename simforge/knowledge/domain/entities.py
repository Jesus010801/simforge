"""Scientific entity descriptors, independent of simulation models."""
from .identity import EntityId, _EntityType, _Lifecycle
from .serialization import _DomainModel, _PositiveInt


class Entity(_DomainModel):

    def _identity_roles(self) -> tuple[tuple[str, str], ...]:
        return ((self.sf_id.value, self.entity_type),)
    sf_id: EntityId
    entity_type: _EntityType
    entity_revision: _PositiveInt
    lifecycle: _Lifecycle = 'ACTIVE'
