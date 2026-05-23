# executors/stream.py
"""
Streaming subprocess helper.

Replaces subprocess.run(capture_output=True) for long-running GROMACS commands.

Key behaviors:
  - Lines emitted to stdout/stderr appear in the terminal in real time
  - All output is also captured in StreamResult for post-run analysis
  - Heartbeat is emitted if no output has appeared for heartbeat_interval_s
  - Hard timeout kills the process (configurable, default 24 h)
  - Interactive prompt detection: warns if the process stops emitting output
    unexpectedly (likely waiting for TTY input)
"""

from __future__ import annotations

import subprocess
import threading
import time
from pathlib import Path
from typing import Callable, Optional


# ═══════════════════════════════════════════════════════════════════════════════
# Result model
# ═══════════════════════════════════════════════════════════════════════════════

class StreamResult:
    __slots__ = (
        "returncode",
        "stdout",
        "stderr",
        "timed_out",
        "hang_killed",
        "elapsed_s",
        "likely_interactive",
    )

    def __init__(self) -> None:
        self.returncode:         Optional[int] = None
        self.stdout:             str           = ""
        self.stderr:             str           = ""
        self.timed_out:          bool          = False
        self.hang_killed:        bool          = False  # killed by hang_timeout_s
        self.elapsed_s:          float         = 0.0
        # True if process appeared to hang waiting for TTY input
        self.likely_interactive: bool          = False


# ═══════════════════════════════════════════════════════════════════════════════
# Core streaming runner
# ═══════════════════════════════════════════════════════════════════════════════

def run_streaming(
    cmd:                  list[str],
    cwd:                  str | Path,
    on_stdout:            Optional[Callable[[str], None]] = None,
    on_stderr:            Optional[Callable[[str], None]] = None,
    on_heartbeat:         Optional[Callable[[float], None]] = None,
    timeout_s:            float = 86400.0,
    heartbeat_interval_s: float = 30.0,
    hang_timeout_s:       float = 600.0,
) -> StreamResult:
    """
    Run *cmd* in *cwd*, streaming output line-by-line.

    Args:
        cmd:                  Command + arguments list.
        cwd:                  Working directory.
        on_stdout:            Called for each stdout line (stripped of newline).
        on_stderr:            Called for each stderr line (stripped of newline).
        on_heartbeat:         Called with seconds-since-last-output when silent.
        timeout_s:            Hard kill timeout in seconds (default 24 h).
        heartbeat_interval_s: Seconds of silence before heartbeat fires.
        hang_timeout_s:       Kill if no output for this many seconds (0 = disabled).

    Returns:
        StreamResult with returncode, captured stdout/stderr, timing info.
    """
    result = StreamResult()
    stdout_lines: list[str] = []
    stderr_lines: list[str] = []

    last_output_time: list[float] = [time.monotonic()]

    try:
        proc = subprocess.Popen(
            cmd,
            cwd=str(cwd),
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )
    except FileNotFoundError as exc:
        result.returncode = 127
        result.stderr = str(exc)
        return result

    start_time = time.monotonic()

    def _drain(stream, sink: list[str], cb: Optional[Callable[[str], None]]) -> None:
        try:
            for raw_line in stream:
                line = raw_line.rstrip("\n")
                sink.append(line)
                last_output_time[0] = time.monotonic()
                if cb:
                    try:
                        cb(line)
                    except Exception:
                        pass
        except Exception:
            pass

    t_out = threading.Thread(
        target=_drain, args=(proc.stdout, stdout_lines, on_stdout), daemon=True
    )
    t_err = threading.Thread(
        target=_drain, args=(proc.stderr, stderr_lines, on_stderr), daemon=True
    )
    t_out.start()
    t_err.start()

    poll_interval = min(heartbeat_interval_s / 4.0, 2.0)

    while True:
        try:
            proc.wait(timeout=poll_interval)
            break
        except subprocess.TimeoutExpired:
            elapsed_total = time.monotonic() - start_time
            if elapsed_total >= timeout_s:
                proc.kill()
                result.timed_out = True
                break

            elapsed_silent = time.monotonic() - last_output_time[0]
            if elapsed_silent >= heartbeat_interval_s and on_heartbeat:
                try:
                    on_heartbeat(elapsed_silent)
                except Exception:
                    pass

            # Semantic hang detection: kill if silent beyond hang_timeout_s
            if hang_timeout_s > 0 and elapsed_silent >= hang_timeout_s:
                proc.kill()
                result.hang_killed        = True
                result.likely_interactive = True
                break

            # Heuristic: if silent for > 2× heartbeat, process may need TTY
            if elapsed_silent >= heartbeat_interval_s * 2:
                result.likely_interactive = True

    t_out.join(timeout=5.0)
    t_err.join(timeout=5.0)

    result.returncode = proc.returncode if not result.timed_out else -1
    result.stdout     = "\n".join(stdout_lines)
    result.stderr     = "\n".join(stderr_lines)
    result.elapsed_s  = time.monotonic() - start_time

    return result
