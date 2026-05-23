"""Sprint 1b — ExecutionBackend abstraction.

Decouples the "how to run a command" concern from RemediationExecutor,
allowing injection of custom backends in tests and future SLURM/cloud runners.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path

from runtime.events import BoundPublisher
from runtime.stream import ProcessResult


class ExecutionBackend(ABC):
    """Abstract backend: submit a shell command and return its result."""

    @abstractmethod
    async def submit(
        self,
        cmd:       list[str],
        cwd:       Path,
        publisher: BoundPublisher,
    ) -> ProcessResult:
        ...


class LocalSubprocessBackend(ExecutionBackend):
    """Wraps AsyncProcessRunner — current local subprocess implementation."""

    async def submit(
        self,
        cmd:       list[str],
        cwd:       Path,
        publisher: BoundPublisher,
    ) -> ProcessResult:
        from runtime.stream import AsyncProcessRunner
        runner = AsyncProcessRunner(publisher)
        return await runner.run(cmd=cmd, cwd=cwd)
