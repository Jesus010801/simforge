"""Revisioned descriptor registry contracts; no backend initialization."""
from .records import RegistryChange, RegistryWrite, RegistryRecord, ExternalIdentifierKey, ExternalIdentifierMatch
from .registry import EntityRegistry
from .errors import (RegistryError, EntityNotFoundError, RegistryRevisionNotFoundError,
                     EntityRevisionNotFoundError, RevisionConflictError,
                     OperationConflictError, InvalidRegistryOperationError)

__all__ = ['RegistryChange', 'RegistryWrite', 'RegistryRecord', 'ExternalIdentifierKey',
           'ExternalIdentifierMatch', 'EntityRegistry', 'RegistryError',
           'EntityNotFoundError', 'RegistryRevisionNotFoundError', 'EntityRevisionNotFoundError',
           'RevisionConflictError', 'OperationConflictError', 'InvalidRegistryOperationError']
