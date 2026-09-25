"""Structured registry failures. Message text is never a handling contract."""
from typing import Literal as _Literal
from uuid import UUID as _UUID
from simforge.knowledge.domain import EntityId as _EntityId

__all__ = ['RegistryError', 'EntityNotFoundError', 'RegistryRevisionNotFoundError',
           'EntityRevisionNotFoundError', 'RevisionConflictError',
           'OperationConflictError', 'InvalidRegistryOperationError']

_Code = _Literal['INVALID_ARGUMENT', 'UNSUPPORTED_ENTITY_MODEL', 'ENTITY_TYPE_CHANGE',
                 'ENTITY_REVISION_SEQUENCE', 'ENTITY_REVISION_REWRITE',
                 'ALIAS_TARGET_MISMATCH', 'ALIAS_TARGET_REVISION_MISSING',
                 'DUPLICATE_ALIAS_DECLARATION']
_Path = tuple[str | int, ...]


def _own(value: object) -> object:
    """Own structured error inputs without retaining caller scalar behavior."""
    if value is None:
        return None
    if isinstance(value, _EntityId):
        return _EntityId.model_validate(value)
    if isinstance(value, _UUID):
        return _UUID(int=_UUID.int.__get__(value))
    if isinstance(value, str):
        return str.__str__(value)
    if type(value) is int:
        return value
    if isinstance(value, int) and not isinstance(value, bool):
        return int.__int__(value)
    if type(value) in (tuple, list):
        return tuple(_own(item) for item in value)
    raise TypeError('unsupported structured error field value')


def _field(index: int) -> property:
    # A property is a data descriptor, so Exception.__dict__ cannot shadow it.
    return property(lambda self: self._details[index])


class RegistryError(Exception):
    """Read-only structured data; ordinary exception traceback remains writable."""
    __slots__ = ('_details',)
    _field_names: tuple[str, ...] = ()

    def _initialize(self, *values: object) -> None:
        object.__setattr__(self, '_details', tuple(_own(value) for value in values))
        Exception.__init__(self, type(self).__name__)

    def __setattr__(self, name: str, value: object) -> None:
        if name in type(self)._field_names or name in ('_details', '_field_names'):
            raise AttributeError('registry error fields are read-only')
        super().__setattr__(name, value)

    def __delattr__(self, name: str) -> None:
        if name in type(self)._field_names or name in ('_details', '_field_names'):
            raise AttributeError('registry error fields are read-only')
        super().__delattr__(name)


class EntityNotFoundError(RegistryError):
    entity_id: _EntityId
    _field_names = ('entity_id',)
    entity_id = _field(0)

    def __init__(self, entity_id: _EntityId) -> None:
        self._initialize(entity_id)


class RegistryRevisionNotFoundError(RegistryError):
    entity_id: _EntityId
    registry_revision: int
    _field_names = ('entity_id', 'registry_revision')
    entity_id = _field(0)
    registry_revision = _field(1)

    def __init__(self, entity_id: _EntityId, registry_revision: int) -> None:
        self._initialize(entity_id, registry_revision)


class EntityRevisionNotFoundError(RegistryError):
    entity_id: _EntityId
    entity_revision: int
    _field_names = ('entity_id', 'entity_revision')
    entity_id = _field(0)
    entity_revision = _field(1)

    def __init__(self, entity_id: _EntityId, entity_revision: int) -> None:
        self._initialize(entity_id, entity_revision)


class RevisionConflictError(RegistryError):
    entity_id: _EntityId
    expected: int
    actual: int
    _field_names = ('entity_id', 'expected', 'actual')
    entity_id = _field(0)
    expected = _field(1)
    actual = _field(2)

    def __init__(self, entity_id: _EntityId, expected: int, actual: int) -> None:
        self._initialize(entity_id, expected, actual)


class OperationConflictError(RegistryError):
    operation_id: _UUID
    _field_names = ('operation_id',)
    operation_id = _field(0)

    def __init__(self, operation_id: _UUID) -> None:
        self._initialize(operation_id)


class InvalidRegistryOperationError(RegistryError):
    code: _Code
    entity_id: _EntityId | None
    field_path: _Path
    validation_issues: tuple[tuple[_Path, str], ...]
    _field_names = ('code', 'entity_id', 'field_path', 'validation_issues')
    code = _field(0)
    entity_id = _field(1)
    field_path = _field(2)
    validation_issues = _field(3)

    def __init__(self, code: _Code, *, entity_id: _EntityId | None = None,
                 field_path: _Path = (),
                 validation_issues: tuple[tuple[_Path, str], ...] = ()) -> None:
        self._initialize(code, entity_id, tuple(field_path),
                         tuple((tuple(path), category) for path, category in validation_issues))
