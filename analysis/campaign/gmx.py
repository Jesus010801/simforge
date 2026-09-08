"""Thin, auditable GROMACS command runner for the study-campaign layer.

There is no shared "run an arbitrary gmx command" abstraction in SimForge
(``executors/`` only runs generated step scripts).  ``analysis/fel/extract.py``
established the pattern used here: ``subprocess.run`` with ``input=`` for group
selections, full capture of argv / stdin / stdout / stderr / return code.

This module adds nothing scientific — it only executes and records.
"""
from __future__ import annotations

import shutil
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

_GMX_VERSION_CACHE: dict[str, Optional[str]] = {}


@dataclass
class GmxResult:
    argv: list[str]
    returncode: int
    stdout: str
    stderr: str
    stdin_text: Optional[str] = None
    duration_s: float = 0.0
    gmx_version: Optional[str] = None
    available: bool = True             # False => binary not found
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.available and self.returncode == 0

    def tail(self, stream: str, n: int = 2000) -> str:
        text = self.stdout if stream == "stdout" else self.stderr
        return text[-n:]

    def to_dict(self) -> dict:
        return {
            "argv": self.argv,
            "returncode": self.returncode,
            "stdin_text": self.stdin_text,
            "duration_s": round(self.duration_s, 3),
            "gmx_version": self.gmx_version,
            "available": self.available,
            "stdout_tail": self.tail("stdout"),
            "stderr_tail": self.tail("stderr"),
            "warnings": self.warnings,
        }


def gmx_available(gmx: str = "gmx") -> bool:
    return shutil.which(gmx) is not None or Path(gmx).is_file()


def gmx_version(gmx: str = "gmx") -> Optional[str]:
    if gmx in _GMX_VERSION_CACHE:
        return _GMX_VERSION_CACHE[gmx]
    version: Optional[str] = None
    if gmx_available(gmx):
        try:
            out = subprocess.run(
                [gmx, "--version"], capture_output=True, text=True, timeout=15,
            )
            text = (out.stdout or "") + (out.stderr or "")
            for line in text.splitlines():
                if "GROMACS version" in line:
                    version = line.split(":")[-1].strip() or None
                    break
        except (OSError, subprocess.SubprocessError):
            version = None
    _GMX_VERSION_CACHE[gmx] = version
    return version


def run_gmx(
    args: list[str],
    *,
    stdin: Optional[str] = None,
    gmx: str = "gmx",
    cwd: Optional[Path] = None,
    timeout: float = 3600.0,
    extra_env: Optional[dict[str, str]] = None,
) -> GmxResult:
    """Run ``gmx <args>`` capturing everything.  Never raises for tool failure.

    ``stdin`` is fed verbatim (append your own trailing newline for group
    selections, matching GROMACS's interactive prompt behaviour).
    """
    argv = [gmx, *args]
    if not gmx_available(gmx):
        return GmxResult(
            argv=argv, returncode=127, stdout="", stderr=f"{gmx}: not found",
            stdin_text=stdin, available=False,
            warnings=[f"GROMACS binary '{gmx}' not found on PATH"],
        )

    env = None
    if extra_env:
        import os
        env = {**os.environ, **extra_env}

    started = time.time()
    try:
        proc = subprocess.run(
            argv,
            input=stdin,
            text=True,
            capture_output=True,
            cwd=str(cwd) if cwd else None,
            timeout=timeout,
            env=env,
        )
        return GmxResult(
            argv=argv, returncode=proc.returncode,
            stdout=proc.stdout or "", stderr=proc.stderr or "",
            stdin_text=stdin, duration_s=time.time() - started,
            gmx_version=gmx_version(gmx),
        )
    except subprocess.TimeoutExpired as exc:
        # TimeoutExpired can carry bytes even with text=True. Metadata inspection
        # must report a failed check rather than crash while joining its output.
        def decoded(value):
            return value.decode('utf-8', errors='replace') if isinstance(value, bytes) else value or ''
        return GmxResult(
            argv=argv, returncode=124, stdout=decoded(exc.stdout), stderr=decoded(exc.stderr),
            stdin_text=stdin, duration_s=time.time() - started,
            gmx_version=gmx_version(gmx),
            warnings=[f"gmx command timed out after {timeout:.0f}s"],
        )
    except OSError as exc:
        return GmxResult(
            argv=argv, returncode=126, stdout="", stderr=str(exc),
            stdin_text=stdin, duration_s=time.time() - started,
            warnings=[f"failed to launch gmx: {exc}"],
        )
