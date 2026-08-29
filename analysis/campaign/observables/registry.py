"""Observable registry.

Adding ``observables/sasa.py`` later only requires calling :func:`register`
(or importing a module that does).  The orchestration layer resolves everything
through :func:`get` / :func:`all_specs` and never contains an
``if analysis == "rmsd": ...`` ladder.
"""
from __future__ import annotations

from typing import Iterable

from analysis.campaign.observables.base import ObservableSpec

_REGISTRY: dict[str, ObservableSpec] = {}


class DuplicateObservableError(ValueError):
    pass


class UnknownObservableError(KeyError):
    pass


def register(observable: ObservableSpec) -> None:
    if not observable.id:
        raise ValueError("observable id must be non-empty")
    if observable.id in _REGISTRY:
        raise DuplicateObservableError(f"observable id already registered: {observable.id!r}")
    _REGISTRY[observable.id] = observable


def get(analysis_id: str) -> ObservableSpec:
    try:
        return _REGISTRY[analysis_id]
    except KeyError:
        raise UnknownObservableError(
            f"unknown analysis id: {analysis_id!r}. Registered: {sorted(_REGISTRY)}"
        ) from None


def is_registered(analysis_id: str) -> bool:
    return analysis_id in _REGISTRY


def ids() -> list[str]:
    return sorted(_REGISTRY)


def all_specs() -> list[ObservableSpec]:
    return [_REGISTRY[k] for k in sorted(_REGISTRY)]


def catalog() -> list[dict]:
    return [s.to_dict() for s in all_specs()]


def _reset_for_tests() -> None:
    _REGISTRY.clear()
    _load_builtins()


def _load_builtins() -> None:
    # register_builtins() is idempotent and re-registers even after a clear(),
    # unlike a one-shot import side-effect.
    from analysis.campaign.observables.rmsd import register_builtins
    register_builtins()


def ensure_loaded() -> None:
    if not _REGISTRY:
        _load_builtins()


# Load built-ins at import time.
_load_builtins()
