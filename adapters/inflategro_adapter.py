"""
adapters/inflategro_adapter.py — adapter for inflategro-Jorge.pl

InflateGRO (Kandt et al. 2007) scales lipid XY coordinates radially
to create space around an embedded membrane protein, or deflates them
iteratively until the system reaches target APL convergence.

This adapter:
  - wraps a single inflategro Perl script invocation
  - extracts APL from area_2.dat in Python (no Fortran AperR needed)
  - returns structured telemetry including current APL in both nm² and Å²
  - does NOT implement the shrink loop — that lives in EmbeddingBuilder

Guaranteed metadata keys on success:
    apl_nm2        float   area per lipid in nm² (raw from area_2.dat)
    apl_ang2       float   area per lipid in Å²  (apl_nm2 × 100)
    scale_applied  float   the scale factor that was used
    lipid_name     str     the lipid residue name that was used
    n_lines_gro    int     number of lines in the output .gro (sanity check)
"""

from __future__ import annotations

import subprocess
from datetime import datetime
from pathlib import Path
from typing import Optional

from adapters.base import (
    AdapterResult,
    AvailabilityResult,
    ExternalToolAdapter,
    PreconditionViolation,
    ToolNotAvailableError,
)


class InflateGROAdapter(ExternalToolAdapter):
    """
    Adapter for inflategro-Jorge.pl.

    Args:
        script_path: Absolute path to the inflategro Perl script.
        timeout_s:   Max seconds to allow the Perl process to run.
    """

    tool_name = "inflategro"

    def __init__(
        self,
        script_path: Path | str,
        timeout_s: int = 300,
    ) -> None:
        self._script = Path(script_path).resolve()
        self._timeout_s = timeout_s

    # ── Availability ──────────────────────────────────────────────────────────

    def check_availability(self) -> AvailabilityResult:
        perl = self._which("perl")
        if perl is None:
            return AvailabilityResult(
                available=False,
                tool_name=self.tool_name,
                reason="perl interpreter not found on PATH",
            )
        if not self._script.exists():
            return AvailabilityResult(
                available=False,
                tool_name=self.tool_name,
                reason=f"inflategro script not found: {self._script}",
            )
        return AvailabilityResult(
            available=True,
            tool_name=self.tool_name,
            binary_path=str(self._script),
        )

    # ── Preconditions ─────────────────────────────────────────────────────────

    def validate_preconditions(self, **kwargs) -> list[PreconditionViolation]:  # type: ignore[override]
        return self._check(
            gro_in=kwargs.get("gro_in"),
            scale=kwargs.get("scale"),
            lipid_name=kwargs.get("lipid_name"),
            cutoff=kwargs.get("cutoff"),
            gro_out=kwargs.get("gro_out"),
            gridsize=kwargs.get("gridsize", 5),
        )

    def _check(
        self,
        gro_in,
        scale,
        lipid_name,
        cutoff,
        gro_out,
        gridsize,
    ) -> list[PreconditionViolation]:
        violations: list[PreconditionViolation] = []

        if gro_in is None:
            violations.append(PreconditionViolation("gro_in", "required, got None"))
        else:
            p = Path(gro_in)
            if not p.exists():
                violations.append(PreconditionViolation("gro_in", f"file not found: {p}"))
            elif p.suffix != ".gro":
                violations.append(
                    PreconditionViolation("gro_in", f"expected .gro file, got {p.suffix}", is_fatal=False)
                )

        if scale is None:
            violations.append(PreconditionViolation("scale", "required, got None"))
        elif not (0.1 <= float(scale) <= 20.0):
            violations.append(
                PreconditionViolation("scale", f"unreasonable scale factor {scale!r}; expected 0.1–20.0")
            )

        if not lipid_name:
            violations.append(PreconditionViolation("lipid_name", "required, must not be empty"))

        if cutoff is None:
            violations.append(PreconditionViolation("cutoff", "required, got None"))
        elif float(cutoff) < 0.0:
            violations.append(PreconditionViolation("cutoff", f"cutoff must be ≥ 0.0, got {cutoff!r}"))

        if gro_out is None:
            violations.append(PreconditionViolation("gro_out", "required, got None"))
        else:
            parent = Path(gro_out).parent
            if not parent.exists():
                violations.append(PreconditionViolation("gro_out", f"parent directory does not exist: {parent}"))

        if gridsize is not None and int(gridsize) <= 0:
            violations.append(PreconditionViolation("gridsize", f"must be > 0, got {gridsize!r}"))

        return violations

    # ── Execution ─────────────────────────────────────────────────────────────

    def run(  # type: ignore[override]
        self,
        *,
        gro_in: Path | str,
        scale: float,
        lipid_name: str,
        cutoff: float,
        gro_out: Path | str,
        gridsize: int = 5,
        area_dat: Optional[Path | str] = None,
        include_protein: bool = False,
    ) -> AdapterResult:
        """
        Run inflategro-Jorge.pl once and return structured telemetry.

        Args:
            gro_in:          Input .gro (protein + bilayer system).
            scale:           Radial scaling factor for lipid XY coords.
                             > 1.0 = inflate; < 1.0 = deflate (shrink step).
            lipid_name:      Residue name of the lipid (e.g. "DPP" for DPPC OPLS-AA).
            cutoff:          Exclusion cutoff around protein in nm; 0 to disable.
            gro_out:         Output .gro path.
            gridsize:        Grid resolution for area calculation (default 5).
            area_dat:        Where to write the APL dat file.
                             Defaults to gro_out.parent/area_2.dat.
            include_protein: Pass "protein" flag to inflategro (only XY protein coords).
        """
        started_at = datetime.now()

        # Normalize paths
        gro_in  = Path(gro_in).resolve()
        gro_out = Path(gro_out).resolve()
        if area_dat is None:
            area_dat = gro_out.parent / "area_2.dat"
        area_dat = Path(area_dat).resolve()

        # Guard: availability + preconditions
        avail = self.check_availability()
        if not avail.available:
            raise ToolNotAvailableError(f"{self.tool_name}: {avail.reason}")

        violations = self._check(
            gro_in=gro_in,
            scale=scale,
            lipid_name=lipid_name,
            cutoff=cutoff,
            gro_out=gro_out,
            gridsize=gridsize,
        )
        fatal = [v for v in violations if v.is_fatal]
        if fatal:
            from adapters.base import PreconditionError
            raise PreconditionError(fatal)

        # Build command
        # inflategro-Jorge.pl takes positional args:
        # bilayer.gro scale_factor lipid_name cutoff_nm gro_out gridsize area_dat [protein]
        cmd = [
            "perl",
            str(self._script),
            str(gro_in),
            str(scale),
            lipid_name,
            str(cutoff),
            str(gro_out),
            str(gridsize),
            str(area_dat),
        ]
        if include_protein:
            cmd.append("protein")

        # Execute
        cwd = gro_out.parent
        try:
            proc = subprocess.run(
                cmd,
                cwd=str(cwd),
                capture_output=True,
                text=True,
                timeout=self._timeout_s,
            )
        except subprocess.TimeoutExpired as exc:
            return self._make_result(
                tool_name=self.tool_name,
                adapter_type=type(self).__name__,
                success=False,
                started_at=started_at,
                error_message=f"Timeout after {self._timeout_s}s",
                stdout=exc.stdout or "",
                stderr=exc.stderr or "",
            )
        except Exception as exc:
            return self._make_result(
                tool_name=self.tool_name,
                adapter_type=type(self).__name__,
                success=False,
                started_at=started_at,
                error_message=f"Subprocess error: {exc}",
            )

        stdout = proc.stdout or ""
        stderr = proc.stderr or ""

        # Non-zero exit
        if proc.returncode != 0:
            return self._make_result(
                tool_name=self.tool_name,
                adapter_type=type(self).__name__,
                success=False,
                started_at=started_at,
                exit_code=proc.returncode,
                stdout=stdout,
                stderr=stderr,
                error_message=f"inflategro exited with code {proc.returncode}",
            )

        # Output validation
        missing = []
        if not gro_out.exists():
            missing.append(str(gro_out))
        if not area_dat.exists():
            missing.append(str(area_dat))

        if missing:
            return self._make_result(
                tool_name=self.tool_name,
                adapter_type=type(self).__name__,
                success=False,
                started_at=started_at,
                exit_code=proc.returncode,
                stdout=stdout,
                stderr=stderr,
                error_message=f"Expected output(s) not found: {missing}",
            )

        # Parse APL from area_2.dat
        apl_nm2, apl_parse_error = self._parse_apl(area_dat)
        n_lines_gro = _count_lines(gro_out)

        metadata: dict = {
            "scale_applied": scale,
            "lipid_name":    lipid_name,
            "n_lines_gro":   n_lines_gro,
        }
        if apl_nm2 is not None:
            metadata["apl_nm2"]  = apl_nm2
            metadata["apl_ang2"] = round(apl_nm2 * 100, 2)
        if apl_parse_error:
            metadata["apl_parse_error"] = apl_parse_error

        return self._make_result(
            tool_name=self.tool_name,
            adapter_type=type(self).__name__,
            success=True,
            started_at=started_at,
            exit_code=proc.returncode,
            stdout=stdout,
            stderr=stderr,
            outputs={
                "gro_out":  str(gro_out),
                "area_dat": str(area_dat),
            },
            metadata=metadata,
        )

    # ── Internal helpers ──────────────────────────────────────────────────────

    @staticmethod
    def _parse_apl(area_dat: Path) -> tuple[Optional[float], Optional[str]]:
        """Read APL in nm² from area_2.dat.

        The Fortran AperR program reads this file with format f5.3,
        meaning a 5-char float with 3 decimal places (e.g. "0.620").
        We parse it directly so AperR compilation is not required.

        Returns (apl_nm2, error_message).  error_message is None on success.
        """
        try:
            raw = area_dat.read_text().strip()
            if not raw:
                return None, "area_2.dat is empty"
            return float(raw), None
        except ValueError:
            return None, f"Could not parse APL from area_2.dat: {raw!r}"
        except Exception as exc:
            return None, f"Could not read area_2.dat: {exc}"


def _count_lines(path: Path) -> int:
    try:
        return sum(1 for _ in path.open())
    except Exception:
        return -1
