"""Owned immutable registry values, deliberately outside sski-json-v1."""
from dataclasses import dataclass as _dataclass
from datetime import datetime as _datetime, timezone as _timezone
import re as _re
import unicodedata as _unicode
from uuid import UUID as _UUID
from typing import TypeVar as _TypeVar
from simforge.knowledge.domain import Entity as _Entity, EntityId as _EntityId, ExternalIdentifier as _ExternalIdentifier
from .errors import InvalidRegistryOperationError as _Invalid

__all__ = ['RegistryChange', 'RegistryWrite', 'RegistryRecord',
           'ExternalIdentifierKey', 'ExternalIdentifierMatch']
_T = _TypeVar('_T', _Entity, _EntityId, _ExternalIdentifier)


def _invalid(path: tuple[str | int, ...], entity_id: _EntityId | None = None) -> None:
    raise _Invalid('INVALID_ARGUMENT', entity_id=entity_id, field_path=path)


def _revision(value: int, path: tuple[str | int, ...], minimum: int = 1,
              entity_id: _EntityId | None = None) -> int:
    if type(value) is not int or value < minimum:
        _invalid(path, entity_id)
    return value


def _known_id(value: object) -> _EntityId | None:
    candidate = getattr(value, 'sf_id', None) if isinstance(value, _Entity) else None
    if type(candidate) is not _EntityId:
        return None
    try:
        return _EntityId.model_validate(candidate)
    except ValueError as error:
        if not callable(getattr(error, 'errors', None)):
            raise
        return None


def _domain(value: _T, cls: type[_T], path: tuple[str | int, ...]) -> _T:
    if cls is _Entity and isinstance(value, _Entity) and type(value) is not _Entity:
        raise _Invalid('UNSUPPORTED_ENTITY_MODEL', entity_id=_known_id(value), field_path=path)
    if type(value) is not cls:
        _invalid(path)
    try:
        # Public validation revalidates instances and recursively owns values.
        source = _alias_input(value, path) if cls is _ExternalIdentifier else value
        return cls.model_validate(source)
    except ValueError as error:
        # Pydantic validation exposes structured errors via its public method.
        # No pydantic/private-domain dependency, no message parsing, and ordinary
        # implementation ValueErrors are not disguised as input failures.
        errors = getattr(error, 'errors', None)
        if not callable(errors):
            raise
        issues = tuple((path + tuple(issue['loc']), issue['type']) for issue in errors())
        raise _Invalid('INVALID_ARGUMENT', entity_id=_known_id(value),
                       field_path=issues[0][0] if issues else path,
                       validation_issues=issues) from error


def _text(value: str, path: tuple[str | int, ...]) -> str:
    if not isinstance(value, str):
        _invalid(path)
    value = str.__str__(value)
    if not value.strip() or any(0xD800 <= ord(c) <= 0xDFFF for c in value):
        _invalid(path)
    return _unicode.normalize('NFC', value)


def _instant(value: _datetime, path: tuple[str | int, ...]) -> _datetime:
    if not isinstance(value, _datetime):
        _invalid(path)
    fields = {name: getattr(_datetime, name).__get__(value) for name in
              ('year', 'month', 'day', 'hour', 'minute', 'second', 'microsecond', 'tzinfo', 'fold')}
    owned = _datetime(**fields)
    try:
        offset = owned.utcoffset()
        if offset is None:
            _invalid(path)
        return owned.replace(tzinfo=_timezone(offset)).astimezone(_timezone.utc)
    except (TypeError, ValueError, OverflowError):
        _invalid(path)


def _alias_input(value: _ExternalIdentifier, path: tuple[str | int, ...]) -> dict:
    """Snapshot the three approved timestamp inputs before domain revalidation.

    Only datetime primitive operations in _instant classify TypeError as bad
    input. Domain validation/serialization TypeErrors remain programming errors.
    The public Python payload is used for validation, never command comparison
    or a registry wire schema. Every other field retains Phase 1 validation.
    """
    payload = _ExternalIdentifier.model_dump(value, mode='python')
    for name in ('valid_from', 'valid_until'):
        if isinstance(payload.get(name), _datetime):
            payload[name] = _instant(payload[name], path + (name,))
    provenance = payload.get('provenance')
    if type(provenance) is dict and isinstance(provenance.get('imported_at'), _datetime):
        provenance['imported_at'] = _instant(provenance['imported_at'], path + ('provenance', 'imported_at'))
    return payload


@_dataclass(frozen=True, slots=True)
class RegistryChange:
    """Caller-supplied administrative attribution, never a scientific claim."""
    operation_id: _UUID
    actor: str
    reason: str
    recorded_at: _datetime

    def __post_init__(self) -> None:
        if not isinstance(self.operation_id, _UUID):
            _invalid(('operation_id',))
        owned = _UUID(int=_UUID.int.__get__(self.operation_id))
        if owned.version != 4:
            _invalid(('operation_id',))
        object.__setattr__(self, 'operation_id', owned)
        object.__setattr__(self, 'actor', _text(self.actor, ('actor',)))
        object.__setattr__(self, 'reason', _text(self.reason, ('reason',)))
        object.__setattr__(self, 'recorded_at', _instant(self.recorded_at, ('recorded_at',)))


def _change(value: RegistryChange, path: tuple[str | int, ...]) -> RegistryChange:
    if type(value) is not RegistryChange:
        _invalid(path)
    try:
        return RegistryChange(value.operation_id, value.actor, value.reason, value.recorded_at)
    except _Invalid as error:
        raise _Invalid(error.code, field_path=path + error.field_path,
                       validation_issues=tuple((path + p, c) for p, c in error.validation_issues)) from error


def _snapshot(value: 'RegistryWrite | RegistryRecord') -> None:
    object.__setattr__(value, 'entity', _domain(value.entity, _Entity, ('entity',)))
    entity_id = value.entity.sf_id
    if type(value.aliases) not in (tuple, list):
        _invalid(('aliases',), entity_id)
    try:
        aliases = tuple(_domain(alias, _ExternalIdentifier, ('aliases', i)) for i, alias in enumerate(value.aliases))
        owned_change = _change(value.change, ('change',))
    except _Invalid as error:
        raise _Invalid(error.code, entity_id=entity_id, field_path=error.field_path,
                       validation_issues=error.validation_issues) from error
    object.__setattr__(value, 'aliases', tuple(sorted(aliases, key=lambda a: a.to_canonical_bytes())))
    object.__setattr__(value, 'change', owned_change)


@_dataclass(frozen=True, slots=True)
class RegistryWrite:
    """Complete proposed current state; target/transition checks occur in write."""
    entity: _Entity
    aliases: tuple[_ExternalIdentifier, ...]
    expected_registry_revision: int
    change: RegistryChange

    def __post_init__(self) -> None:
        _revision(self.expected_registry_revision, ('expected_registry_revision',), 0, _known_id(self.entity))
        _snapshot(self)


@_dataclass(frozen=True, slots=True)
class RegistryRecord:
    """One immutable snapshot in a per-EntityId registry revision sequence."""
    registry_revision: int
    entity: _Entity
    aliases: tuple[_ExternalIdentifier, ...]
    change: RegistryChange

    def __post_init__(self) -> None:
        _revision(self.registry_revision, ('registry_revision',), entity_id=_known_id(self.entity))
        _snapshot(self)


@_dataclass(frozen=True, slots=True)
class ExternalIdentifierKey:
    namespace: str
    accession: str

    def __post_init__(self) -> None:
        namespace = _text(self.namespace, ('namespace',))
        if not _re.fullmatch(r'[a-z][a-z0-9_.-]*', namespace):
            _invalid(('namespace',))
        object.__setattr__(self, 'namespace', namespace)
        object.__setattr__(self, 'accession', _text(self.accession, ('accession',)))


@_dataclass(frozen=True, slots=True)
class ExternalIdentifierMatch:
    """Current descriptor and exact supplied alias target are distinct."""
    entity: _Entity
    registry_revision: int
    alias: _ExternalIdentifier

    def __post_init__(self) -> None:
        _revision(self.registry_revision, ('registry_revision',), entity_id=_known_id(self.entity))
        object.__setattr__(self, 'entity', _domain(self.entity, _Entity, ('entity',)))
        try:
            owned_alias = _domain(self.alias, _ExternalIdentifier, ('alias',))
        except _Invalid as error:
            raise _Invalid(error.code, entity_id=self.entity.sf_id, field_path=error.field_path,
                           validation_issues=error.validation_issues) from error
        object.__setattr__(self, 'alias', owned_alias)
