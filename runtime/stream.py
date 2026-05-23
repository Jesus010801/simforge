from __future__ import annotations

import asyncio
import re
import time
from dataclasses import dataclass, field
from typing import AsyncIterator

from runtime.events import BoundPublisher, EventSeverity, EventType

# GROMACS interactive-prompt patterns that need stdin response
_INTERACTIVE_PATTERNS: list[re.Pattern] = [
    re.compile(r"Select a group:", re.IGNORECASE),
    re.compile(r"Select group for"),
    re.compile(r"Enter a selection"),
    re.compile(r"choice\s*:"),
]

# Parses "Performance:  X ns/day  Y hours/ns"
_PERFORMANCE_RE = re.compile(
    r"Performance:\s+([\d.]+)\s+ns/day\s+([\d.]+)\s+hours/ns"
)

# Hang: no stdout for this many seconds → WARNING event
_HANG_TIMEOUT_S = 120

# 10 MB per stream — well above any real GROMACS output line.
# asyncio.StreamReader default is 64 KB, which mdrun log blocks can exceed.
_STREAM_LIMIT = 10 * 1024 * 1024


@dataclass
class ProcessResult:
    returncode:   int
    ns_per_day:   float | None       = None
    hours_per_ns: float | None       = None
    wall_time_s:  float              = 0.0
    stdout_lines: list[str]          = field(default_factory=list)
    stderr_lines: list[str]          = field(default_factory=list)


class AsyncProcessRunner:
    """
    Runs an external command asynchronously with line-by-line streaming.

    Design principles:
      - Process supervision is decoupled from observability: even if all
        streaming fails, proc.wait() is always called and returncode captured.
      - Buffer limit is 10 MB — prevents LimitOverrunError for all real GROMACS
        output. A safety-net handler discards oversized lines rather than crashing.
      - The heartbeat runs as a separate cancellable task, not inside gather,
        so it cannot prevent the run from completing when drains finish.
      - Streaming exceptions are reported as WARNING events and never propagated
        to the caller.
    """

    def __init__(
        self,
        publisher:       BoundPublisher,
        hang_timeout_s:  int = _HANG_TIMEOUT_S,
        heartbeat_s:     int = 30,
    ) -> None:
        self._pub           = publisher
        self._hang_timeout  = hang_timeout_s
        self._heartbeat_s   = heartbeat_s

    async def run(
        self,
        cmd:         list[str],
        cwd:         str | None = None,
        env:         dict[str, str] | None = None,
        stdin_text:  str | None = None,
    ) -> ProcessResult:
        t0 = time.monotonic()

        proc = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            stdin=asyncio.subprocess.PIPE if stdin_text else None,
            cwd=cwd,
            env=env,
            limit=_STREAM_LIMIT,
        )

        if stdin_text and proc.stdin:
            proc.stdin.write(stdin_text.encode())
            await proc.stdin.drain()
            proc.stdin.close()

        stdout_lines: list[str] = []
        stderr_lines: list[str] = []
        ns_per_day:   float | None = None
        hours_per_ns: float | None = None

        async def _drain_stdout() -> None:
            nonlocal ns_per_day, hours_per_ns
            assert proc.stdout
            async for line in self._iter_lines(proc.stdout):
                text = line.rstrip("\n")
                stdout_lines.append(text)
                self._pub.emit(EventType.STDOUT, message=text)
                self._check_interactive(text)
                perf = _PERFORMANCE_RE.search(text)
                if perf:
                    ns_per_day   = float(perf.group(1))
                    hours_per_ns = float(perf.group(2))
                    self._pub.emit(
                        EventType.PERFORMANCE,
                        message=f"{ns_per_day:.2f} ns/day",
                        ns_per_day=ns_per_day,
                        hours_per_ns=hours_per_ns,
                    )

        async def _drain_stderr() -> None:
            assert proc.stderr
            async for line in self._iter_lines(proc.stderr):
                text = line.rstrip("\n")
                stderr_lines.append(text)
                sev = EventSeverity.WARNING if "warning" in text.lower() else EventSeverity.INFO
                self._pub.emit(EventType.STDERR, message=text, severity=sev)

        # Heartbeat runs as an independent task cancelled when draining is done.
        # Keeping it outside the drain gather prevents it from blocking completion.
        async def _heartbeat() -> None:
            try:
                while True:
                    await asyncio.sleep(self._heartbeat_s)
                    self._pub.emit(EventType.HEARTBEAT, message="process alive")
            except asyncio.CancelledError:
                pass

        heartbeat_task = asyncio.create_task(_heartbeat())
        try:
            drain_results = await asyncio.gather(
                _drain_stdout(),
                _drain_stderr(),
                return_exceptions=True,
            )
        finally:
            heartbeat_task.cancel()
            try:
                await heartbeat_task
            except asyncio.CancelledError:
                pass

        # Report streaming failures as warnings — they don't affect correctness.
        for exc in drain_results:
            if isinstance(exc, Exception):
                self._pub.emit(
                    EventType.WARNING,
                    message=f"Stream error (output may be incomplete): {exc}",
                    severity=EventSeverity.WARNING,
                )

        # Always wait for the process — returncode must be captured regardless
        # of what happened in the observability layer.
        await proc.wait()

        return ProcessResult(
            returncode   = proc.returncode or 0,
            ns_per_day   = ns_per_day,
            hours_per_ns = hours_per_ns,
            wall_time_s  = time.monotonic() - t0,
            stdout_lines = stdout_lines,
            stderr_lines = stderr_lines,
        )

    # ── Helpers ───────────────────────────────────────────────────────────────

    @staticmethod
    async def _iter_lines(stream: asyncio.StreamReader) -> AsyncIterator[str]:
        """
        Yield lines from *stream* with resilient error handling.

        asyncio.LimitOverrunError is raised when readline() cannot find a
        newline within the StreamReader buffer (_STREAM_LIMIT = 10 MB).
        When that happens the oversized segment is discarded so the stream
        can continue — that one line is lost from observability but the
        process is unaffected.
        """
        while True:
            try:
                line = await stream.readline()
            except asyncio.LimitOverrunError as exc:
                # Line exceeds the 10 MB buffer — extremely unusual.
                # Read and discard the stuck segment to unblock readline.
                try:
                    await stream.read(exc.consumed)
                except Exception:
                    return
                continue
            except Exception:
                return
            if not line:
                return
            yield line.decode("utf-8", errors="replace")

    def _check_interactive(self, text: str) -> None:
        for pat in _INTERACTIVE_PATTERNS:
            if pat.search(text):
                self._pub.emit(
                    EventType.INTERACTIVE_PROMPT,
                    message=text,
                    severity=EventSeverity.WARNING,
                    pattern=pat.pattern,
                )
                break
