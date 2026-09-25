"""Volatile semantic reference registry. Each instance owns independent state."""
from threading import RLock as _RLock
from dataclasses import dataclass as _dataclass
from uuid import UUID as _UUID
from simforge.knowledge.domain import Entity as _Entity, EntityId as _EntityId
from simforge.knowledge.identity.registry import EntityRegistry as _EntityRegistry
from simforge.knowledge.identity.records import (
    RegistryWrite as _Write, RegistryRecord as _Record,
    ExternalIdentifierKey as _Key, ExternalIdentifierMatch as _Match,
    _domain, _revision, _invalid,
)
from simforge.knowledge.identity.errors import (
    EntityNotFoundError as _EntityMissing, RegistryRevisionNotFoundError as _RegistryMissing,
    EntityRevisionNotFoundError as _RevisionMissing, RevisionConflictError as _Conflict,
    OperationConflictError as _OperationConflict, InvalidRegistryOperationError as _Invalid,
)

__all__ = ['InMemoryEntityRegistry']


def _copy(record: _Record) -> _Record:
    return _Record(record.registry_revision, record.entity, record.aliases, record.change)


def _command_key(command: _Write) -> tuple:
    # This tuple is a process-local comparison, NOT a wire format or content ID.
    return (command.entity.to_canonical_bytes(),
            tuple(a.to_canonical_bytes() for a in command.aliases),
            command.expected_registry_revision, command.change)


@_dataclass(frozen=True, slots=True)
class _State:
    """Owned maps, never mutated after publication; no public storage interface."""
    history: dict[_EntityId, tuple[_Record, ...]]
    heads: dict[_EntityId, _Record]
    entities: dict[_EntityId, dict[int, _Entity]]
    index: dict[_Key, tuple[_Match, ...]]
    operations: dict[_UUID, tuple[tuple, _Record]]


class InMemoryEntityRegistry(_EntityRegistry):
    """Thread-safe, descriptor-only registry with unbounded immutable history.

    The lock encloses each entire operation; correctness never depends on dict
    atomicity. A complete replacement state is prepared before one publication;
    no published map is mutated and allocation failures leave it untouched.
    No clocks, persistence, scientific payloads, or global registry state.
    """

    def __init__(self) -> None:
        self._lock = _RLock()
        self._state = _State({}, {}, {}, {}, {})

    def write(self, command: _Write) -> _Record:
        with self._lock:
            if type(command) is not _Write:
                _invalid(('command',))
            command = _Write(command.entity, command.aliases, command.expected_registry_revision, command.change)
            state = self._state
            key = _command_key(command)
            operation = command.change.operation_id
            if operation in state.operations:
                original, result = state.operations[operation]
                if original != key:
                    raise _OperationConflict(operation)
                return _copy(result)
            entity = command.entity
            entity_id = entity.sf_id
            head = state.heads.get(entity_id)
            actual = head.registry_revision if head else 0
            if command.expected_registry_revision != actual:
                raise _Conflict(entity_id, command.expected_registry_revision, actual)
            self._validate_transition(entity, head)
            prepared, result = self._prepare(state, command, key, actual)
            # The sole publication point. All potentially failing preparation,
            # including the journal and detached return value, precedes this.
            self._state = prepared
            return result

    @staticmethod
    def _prepare(state: _State, command: _Write, key: tuple, actual: int) -> tuple[_State, _Record]:
        entity = command.entity
        entity_id = entity.sf_id
        descriptors = dict(state.entities.get(entity_id, {}))
        descriptors[entity.entity_revision] = entity
        seen: set[bytes] = set()
        for i, alias in enumerate(command.aliases):
            path = ('aliases', i)
            if alias.entity.entity_id != entity_id or alias.entity.entity_type != entity.entity_type:
                raise _Invalid('ALIAS_TARGET_MISMATCH', entity_id=entity_id, field_path=path + ('entity',))
            if alias.entity.entity_revision not in descriptors:
                raise _Invalid('ALIAS_TARGET_REVISION_MISSING', entity_id=entity_id, field_path=path + ('entity', 'entity_revision'))
            encoded = alias.to_canonical_bytes()
            if encoded in seen:
                raise _Invalid('DUPLICATE_ALIAS_DECLARATION', entity_id=entity_id, field_path=path)
            seen.add(encoded)
        record = _Record(actual + 1, entity, command.aliases, command.change)
        history = state.history.get(entity_id, ()) + (record,)
        # Prepare a new projection without exposing partially replaced sets.
        index = {k: tuple(m for m in matches if m.entity.sf_id != entity_id)
                 for k, matches in state.index.items()}
        for alias in record.aliases:
            alias_key = _Key(alias.namespace, alias.accession)
            match = _Match(record.entity, record.registry_revision, alias)
            index[alias_key] = index.get(alias_key, ()) + (match,)
        index = {k: tuple(sorted(matches, key=lambda m: (m.entity.sf_id.value, m.alias.to_canonical_bytes())))
                 for k, matches in index.items() if matches}
        result = _copy(record)
        histories = dict(state.history)
        histories[entity_id] = history
        heads = dict(state.heads)
        heads[entity_id] = record
        entities = dict(state.entities)
        entities[entity_id] = descriptors
        operations = dict(state.operations)
        operations[command.change.operation_id] = (key, record)
        prepared = _State(histories, heads, entities, index, operations)
        return prepared, result

    @staticmethod
    def _validate_transition(entity: _Entity, head: _Record | None) -> None:
        def fail(code: str, field: str) -> None:
            raise _Invalid(code, entity_id=entity.sf_id, field_path=('entity', field))
        if head is None:
            if entity.entity_revision != 1:
                fail('ENTITY_REVISION_SEQUENCE', 'entity_revision')
            return
        if entity.entity_type != head.entity.entity_type:
            fail('ENTITY_TYPE_CHANGE', 'entity_type')
        current = head.entity.entity_revision
        if entity.entity_revision not in (current, current + 1):
            fail('ENTITY_REVISION_SEQUENCE', 'entity_revision')
        if entity.entity_revision == current and entity.to_canonical_bytes() != head.entity.to_canonical_bytes():
            fail('ENTITY_REVISION_REWRITE', 'entity_revision')

    def _require(self, entity_id: _EntityId) -> _Record:
        if entity_id not in self._state.heads:
            raise _EntityMissing(entity_id)
        return self._state.heads[entity_id]

    def get(self, entity_id: _EntityId, *, registry_revision: int | None = None) -> _Record:
        with self._lock:
            entity_id = _domain(entity_id, _EntityId, ('entity_id',))
            if registry_revision is not None:
                _revision(registry_revision, ('registry_revision',), entity_id=entity_id)
            head = self._require(entity_id)
            if registry_revision is None:
                return _copy(head)
            if registry_revision > head.registry_revision:
                raise _RegistryMissing(entity_id, registry_revision)
            return _copy(self._state.history[entity_id][registry_revision - 1])

    def get_entity_revision(self, entity_id: _EntityId, entity_revision: int) -> _Entity:
        with self._lock:
            entity_id = _domain(entity_id, _EntityId, ('entity_id',))
            _revision(entity_revision, ('entity_revision',), entity_id=entity_id)
            self._require(entity_id)
            if entity_revision not in self._state.entities[entity_id]:
                raise _RevisionMissing(entity_id, entity_revision)
            return _Entity.model_validate(self._state.entities[entity_id][entity_revision])

    def history(self, entity_id: _EntityId) -> tuple[_Record, ...]:
        with self._lock:
            entity_id = _domain(entity_id, _EntityId, ('entity_id',))
            self._require(entity_id)
            return tuple(_copy(record) for record in self._state.history[entity_id])

    def lookup_external(self, key: _Key) -> tuple[_Match, ...]:
        with self._lock:
            if type(key) is not _Key:
                _invalid(('key',))
            try:
                key = _Key(key.namespace, key.accession)
            except _Invalid as error:
                raise _Invalid(error.code, field_path=('key',) + error.field_path,
                               validation_issues=error.validation_issues) from error
            return tuple(_Match(m.entity, m.registry_revision, m.alias) for m in self._state.index.get(key, ()))
