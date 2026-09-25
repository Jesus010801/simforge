"""The single intermediate registry (mirrors ``observables.registry``)."""
from __future__ import annotations

from contextlib import contextmanager

from analysis.campaign.intermediates.base import IntermediateSpec

_REGISTRY: dict[str, IntermediateSpec] = {}
_LOADED = False


def register_intermediate(spec: IntermediateSpec) -> IntermediateSpec:
    if not spec.id:
        raise ValueError("an intermediate needs an id")
    if spec.id in _REGISTRY and _REGISTRY[spec.id] is not spec:
        raise ValueError(f"intermediate {spec.id!r} is already registered")
    _REGISTRY[spec.id] = spec
    return spec


def ensure_loaded() -> None:
    global _LOADED
    if not _LOADED:
        _LOADED = True
        import analysis.campaign.intermediates.selection_centres  # noqa: F401


def get(intermediate_id: str) -> IntermediateSpec:
    ensure_loaded()
    if intermediate_id not in _REGISTRY:
        raise KeyError(f"unknown intermediate {intermediate_id!r}; registered: {ids()}")
    return _REGISTRY[intermediate_id]


def is_registered(intermediate_id: str) -> bool:
    ensure_loaded()
    return intermediate_id in _REGISTRY


def ids() -> list[str]:
    ensure_loaded()
    return sorted(_REGISTRY)


@contextmanager
def temporarily_registered(*specs: IntermediateSpec):
    """Tests / experiments: register specs for the duration of a block."""
    ensure_loaded()
    added = []
    for s in specs:
        if s.id not in _REGISTRY:
            _REGISTRY[s.id] = s
            added.append(s.id)
    try:
        yield
    finally:
        for sid in added:
            _REGISTRY.pop(sid, None)
