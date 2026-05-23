"""
adapters/base.py — contracts for external-tool adapters.

External tools (Perl scripts, Fortran binaries, web services) are the
most fragile boundary in SimForge.  This module defines the contracts
that every adapter must honour so the rest of the system can reason
about external tools uniformly.

Design rules:
  - AdapterResult is always returned (never swallowed).
  - ToolNotAvailableError / PreconditionError are raised immediately —
    they represent programmer/environment errors, not runtime failures.
  - Subprocess failures become AdapterResult(success=False).
  - Adapters never attempt remediation — they report faithfully.
  - All I/O paths are absolute.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Optional

from pydantic import BaseModel, Field


# ── Error hierarchy ───────────────────────────────────────────────────────────

class AdapterError(Exception):
    """Base for all adapter errors."""


class ToolNotAvailableError(AdapterError):
    """External tool binary/script cannot be found or executed.

    Raised during check_availability().  Callers should catch this and
    surface it as a configuration problem, not a simulation failure.
    """


class PreconditionError(AdapterError):
    """One or more preconditions for running the tool were violated.

    Raised during validate_preconditions() when violations are fatal.
    Contains the full list of violations for structured reporting.
    """

    def __init__(self, violations: list[PreconditionViolation]) -> None:
        self.violations = violations
        super().__init__(
            "Precondition violations:\n"
            + "\n".join(f"  [{v.field}] {v.message}" for v in violations)
        )


# ── Supporting data types ─────────────────────────────────────────────────────

@dataclass(frozen=True)
class AvailabilityResult:
    """Result of checking whether an external tool can be invoked."""
    available: bool
    tool_name: str
    reason: Optional[str] = None           # why unavailable; None when available
    binary_path: Optional[str] = None      # resolved binary path if found


@dataclass(frozen=True)
class PreconditionViolation:
    """A single violated precondition."""
    field: str       # input parameter name or file path label
    message: str     # human-readable reason
    is_fatal: bool = True   # if True, run() must not proceed


# ── AdapterResult — structured telemetry ─────────────────────────────────────

class AdapterResult(BaseModel):
    """
    Structured result returned by every adapter.run() call.

    Always returned — never None.  Callers inspect .success to decide
    what to do next.  Validators consume .metadata for convergence checks.
    Telemetry systems log the full record to disk.
    """
    # Identity
    tool_name:    str
    adapter_type: str                       # concrete adapter class name

    # Outcome
    success:      bool
    exit_code:    Optional[int]   = None
    error_message: Optional[str]  = None

    # I/O transcript
    stdout:       str             = ""
    stderr:       str             = ""

    # Timing
    started_at:   datetime        = Field(default_factory=datetime.now)
    finished_at:  Optional[datetime] = None
    elapsed_s:    Optional[float] = None

    # Resolved outputs: logical_name → absolute_path_str
    # e.g. {"gro_out": "/path/system.gro", "area_dat": "/path/area_2.dat"}
    outputs:      dict[str, str]  = Field(default_factory=dict)

    # Adapter-specific structured metadata.
    # Validators and telemetry systems consume these fields.
    # Each adapter documents the keys it guarantees.
    metadata:     dict[str, Any]  = Field(default_factory=dict)

    def output_path(self, name: str) -> Optional[Path]:
        """Return a named output as a Path, or None if not present."""
        raw = self.outputs.get(name)
        return Path(raw) if raw else None


# ── ExternalToolAdapter — abstract base ──────────────────────────────────────

class ExternalToolAdapter(ABC):
    """
    Abstract base for adapters that wrap external tools.

    Subclasses must implement:
        tool_name         — class-level string identifier
        check_availability()   — probe whether the tool exists
        validate_preconditions(**kwargs) — validate inputs before run
        run(**kwargs) → AdapterResult    — execute and return structured result

    Typical call sequence:
        adapter.assert_available()          # raises ToolNotAvailableError
        violations = adapter.validate_preconditions(...)
        if violations: handle_or_raise()
        result = adapter.run(...)
        if not result.success: ...
    """

    tool_name: str = "external_tool"    # override in subclass

    # ── Availability ──────────────────────────────────────────────────────────

    @abstractmethod
    def check_availability(self) -> AvailabilityResult:
        """Probe whether the tool can be invoked on this machine.

        Must not modify any files.  Should be fast (< 1s).
        """
        ...

    def assert_available(self) -> None:
        """Raise ToolNotAvailableError if the tool is not available."""
        result = self.check_availability()
        if not result.available:
            raise ToolNotAvailableError(
                f"{self.tool_name} is not available: {result.reason}"
            )

    # ── Preconditions ─────────────────────────────────────────────────────────

    @abstractmethod
    def validate_preconditions(self, **kwargs) -> list[PreconditionViolation]:
        """Validate inputs without running the tool.

        Returns a list of violations.  Empty list = all preconditions met.
        Does not raise — callers decide how to handle non-fatal violations.
        """
        ...

    def assert_preconditions(self, **kwargs) -> None:
        """Raise PreconditionError if any fatal precondition is violated."""
        violations = self.validate_preconditions(**kwargs)
        fatal = [v for v in violations if v.is_fatal]
        if fatal:
            raise PreconditionError(fatal)

    # ── Execution ─────────────────────────────────────────────────────────────

    @abstractmethod
    def run(self, **kwargs) -> AdapterResult:
        """Run the external tool and return a structured result.

        Implementations must:
          1. Record started_at before subprocess launch.
          2. Capture stdout + stderr.
          3. Populate outputs dict with resolved absolute paths.
          4. Populate metadata with adapter-specific convergence data.
          5. Set success=False on non-zero exit or missing outputs.
          6. Never raise for execution failures — return them in AdapterResult.

        May raise ToolNotAvailableError or PreconditionError if callers
        skip the assert_available / assert_preconditions steps.
        """
        ...

    # ── Helpers available to subclasses ──────────────────────────────────────

    @staticmethod
    def _which(binary: str) -> Optional[str]:
        """Return the resolved path of a binary, or None if not found."""
        import shutil
        return shutil.which(binary)

    @staticmethod
    def _make_result(
        *,
        tool_name: str,
        adapter_type: str,
        success: bool,
        started_at: datetime,
        exit_code: Optional[int] = None,
        stdout: str = "",
        stderr: str = "",
        error_message: Optional[str] = None,
        outputs: Optional[dict[str, str]] = None,
        metadata: Optional[dict[str, Any]] = None,
    ) -> AdapterResult:
        """Convenience constructor that fills in finished_at and elapsed_s."""
        now = datetime.now()
        return AdapterResult(
            tool_name=tool_name,
            adapter_type=adapter_type,
            success=success,
            exit_code=exit_code,
            error_message=error_message,
            stdout=stdout,
            stderr=stderr,
            started_at=started_at,
            finished_at=now,
            elapsed_s=(now - started_at).total_seconds(),
            outputs=outputs or {},
            metadata=metadata or {},
        )
