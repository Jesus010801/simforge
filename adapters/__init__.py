"""
adapters/ — external-tool adapter layer for SimForge.

Each adapter wraps a single external tool (Perl script, Fortran binary,
web service) behind a stable contract defined in adapters/base.py.

Available adapters:
    InflateGROAdapter   — inflategro-Jorge.pl  (lipid scaling / APL measurement)
    WaterDeletorAdapter — water_deletor logic  (pure Python reimplementation)
    MoveMembAdapter     — MoveMemb.f logic     (pure Python reimplementation)
"""

from adapters.base import (
    AdapterError,
    AdapterResult,
    AvailabilityResult,
    ExternalToolAdapter,
    PreconditionError,
    PreconditionViolation,
    ToolNotAvailableError,
)
from adapters.inflategro_adapter import InflateGROAdapter
from adapters.water_deletor_adapter import WaterDeletorAdapter
from adapters.movememb_adapter import MoveMembAdapter

__all__ = [
    "AdapterError",
    "AdapterResult",
    "AvailabilityResult",
    "ExternalToolAdapter",
    "InflateGROAdapter",
    "MoveMembAdapter",
    "PreconditionError",
    "PreconditionViolation",
    "ToolNotAvailableError",
    "WaterDeletorAdapter",
]
