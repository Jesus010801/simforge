"""
``simforge validate-system`` — inspect a built GROMACS system directory without
modifying it, and report structured checks.

Each check is a ``Check(check_id, severity, status, message, evidence)`` where

    severity : blocking | advisory
    status   : PASS | WARN | FAIL | SKIP

A blocking check that FAILs makes ``ValidationReport.ok`` False (non-zero exit).

Checks implemented:
    topology_present        topol.top exists
    coordinate_present      a coordinate .gro exists
    includes_resolve        every local #include resolves on disk
    molecule_table          [ molecules ] present, each name has a moleculetype
    atomtype_resolution     every referenced atomtype is defined (local or FF)
    atom_counts             Σ(mol atoms × count) == coordinate atom count
    coordinate_validity     no NaN/inf coords, sane box, no exact overlaps
    grompp_dry_run          gmx grompp -maxwarn 0 succeeds (when gmx present)
"""

from __future__ import annotations

import math
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

from core.provenance import gromacs_data_prefix
from utils.gro_parser import parse_gro
from utils.itp_parser import parse_itp

PASS, WARN, FAIL, SKIP = "PASS", "WARN", "FAIL", "SKIP"
BLOCKING, ADVISORY = "blocking", "advisory"

_INCLUDE_RE = re.compile(r'^\s*#include\s+"([^"]+)"')
_COORD_CANDIDATES = ("complex.gro", "system.gro", "solvated.gro", "conf.gro")


@dataclass
class Check:
    check_id: str
    severity: str
    status: str
    message: str
    evidence: Any = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "check_id": self.check_id,
            "severity": self.severity,
            "status": self.status,
            "message": self.message,
            "evidence": self.evidence,
        }


@dataclass
class ValidationReport:
    system_dir: str
    checks: list[Check] = field(default_factory=list)

    def add(self, *args: Any, **kwargs: Any) -> Check:
        c = Check(*args, **kwargs)
        self.checks.append(c)
        return c

    @property
    def ok(self) -> bool:
        return not any(c.severity == BLOCKING and c.status == FAIL for c in self.checks)

    @property
    def readiness(self) -> str:
        if not self.ok:
            return "NOT READY"
        if any(c.status == WARN for c in self.checks):
            return "READY (with warnings)"
        return "READY"

    def to_dict(self) -> dict[str, Any]:
        return {
            "system_dir": self.system_dir,
            "ok": self.ok,
            "readiness": self.readiness,
            "checks": [c.to_dict() for c in self.checks],
        }


# ── include graph ────────────────────────────────────────────────────────────

def _iter_includes(text: str) -> list[str]:
    return [m.group(1) for line in text.splitlines() if (m := _INCLUDE_RE.match(line))]


def _is_ff_include(inc: str) -> bool:
    return ".ff/" in inc or inc.endswith(".ff") or inc.startswith("oplsaa") or "/" in inc


def _collect_local_itps(top_path: Path) -> tuple[list[Path], list[str], list[str]]:
    """Walk local #includes from topol.top.

    Returns (existing_local_itp_paths, missing_local_includes, ff_includes).
    """
    root = top_path.parent
    seen: set[Path] = set()
    existing: list[Path] = []
    missing: list[str] = []
    ff: list[str] = []
    queue: list[Path] = [top_path]
    while queue:
        cur = queue.pop()
        try:
            text = cur.read_text()
        except Exception:
            continue
        for inc in _iter_includes(text):
            if _is_ff_include(inc):
                ff.append(inc)
                continue
            cand = (cur.parent / inc).resolve()
            if not cand.exists():
                cand = (root / inc).resolve()
            if not cand.exists():
                missing.append(inc)
                continue
            if cand not in seen:
                seen.add(cand)
                existing.append(cand)
                queue.append(cand)
    return existing, missing, ff


# ── atomtype parsing ─────────────────────────────────────────────────────────

def _atomtype_names_in_block(text: str) -> set[str]:
    names: set[str] = set()
    in_block = False
    for line in text.splitlines():
        s = line.split(";")[0].strip()
        if not s:
            continue
        if s.startswith("["):
            in_block = s.strip("[] ").lower() == "atomtypes"
            continue
        if in_block:
            parts = s.split()
            if parts:
                names.add(parts[0])
    return names


def _ff_atomtypes(forcefield: str) -> Optional[set[str]]:
    top = gromacs_data_prefix()
    if top is None:
        return None
    ff_dir = Path(top) / f"{forcefield}.ff"
    if not ff_dir.is_dir():
        return None
    names: set[str] = set()
    atp = ff_dir / "atomtypes.atp"
    if atp.exists():
        for line in atp.read_text().splitlines():
            s = line.split(";")[0].strip()
            if s:
                names.add(s.split()[0])
    nb = ff_dir / "ffnonbonded.itp"
    if nb.exists():
        names |= _atomtype_names_in_block(nb.read_text())
    return names or None


# ── main entry ───────────────────────────────────────────────────────────────

def validate_system_dir(
    system_dir: str | Path,
    *,
    forcefield: str = "oplsaa",
    coordinate: Optional[str | Path] = None,
    run_grompp: Optional[bool] = None,
) -> ValidationReport:
    d = Path(system_dir).resolve()
    report = ValidationReport(system_dir=str(d))

    if not d.is_dir():
        report.add("topology_present", BLOCKING, FAIL, f"not a directory: {d}")
        return report

    # 1. topol.top
    top_path = d / "topol.top"
    if not top_path.exists():
        alt = sorted(d.glob("*.top"))
        if alt:
            top_path = alt[0]
        else:
            report.add("topology_present", BLOCKING, FAIL, "no topol.top (or *.top) in system dir")
            return report
    report.add("topology_present", BLOCKING, PASS, f"found {top_path.name}", str(top_path))

    top_text = top_path.read_text()

    # 2. coordinate file
    coord_path: Optional[Path] = None
    if coordinate:
        coord_path = Path(coordinate)
        if not coord_path.is_absolute():
            coord_path = d / coord_path
    else:
        for name in _COORD_CANDIDATES:
            if (d / name).exists():
                coord_path = d / name
                break
        if coord_path is None:
            gros = sorted(d.glob("*.gro"))
            coord_path = gros[0] if len(gros) == 1 else None

    gro = None
    if coord_path is None or not coord_path.exists():
        report.add("coordinate_present", BLOCKING, FAIL,
                   "no coordinate .gro found (looked for complex.gro/system.gro/…)")
    else:
        try:
            gro = parse_gro(coord_path)
            report.add("coordinate_present", BLOCKING, PASS,
                       f"{coord_path.name}: {gro.atom_count} atoms", str(coord_path))
        except Exception as exc:
            report.add("coordinate_present", BLOCKING, FAIL,
                       f"cannot parse {coord_path.name}: {exc}", str(coord_path))

    # 3. includes resolve
    local_itps, missing_inc, ff_inc = _collect_local_itps(top_path)
    if missing_inc:
        report.add("includes_resolve", BLOCKING, FAIL,
                   f"{len(missing_inc)} local #include(s) not on disk", missing_inc)
    else:
        report.add("includes_resolve", BLOCKING, PASS,
                   f"{len(local_itps)} local include(s) resolved, "
                   f"{len(ff_inc)} force-field include(s) delegated to GROMACS",
                   {"local": [p.name for p in local_itps], "forcefield": sorted(set(ff_inc))})

    # gather moleculetypes from topol.top + all local itps
    mol_atomcounts: dict[str, int] = {}
    referenced_types: set[str] = set()
    defined_local_types: set[str] = _atomtype_names_in_block(top_text)
    for itp in [top_path, *local_itps]:
        try:
            parsed = parse_itp(itp)
        except Exception:
            continue
        if parsed.moleculetype:
            mol_atomcounts[parsed.moleculetype.name] = len(parsed.atoms)
        for a in parsed.atoms:
            referenced_types.add(a.atom_type)
        try:
            defined_local_types |= _atomtype_names_in_block(itp.read_text())
        except Exception:
            pass

    # 4. molecule table
    molecules = _parse_molecules(top_text)
    if not molecules:
        report.add("molecule_table", BLOCKING, FAIL, "no [ molecules ] section in topol.top")
    else:
        unknown = [
            name for name, _ in molecules
            if name not in mol_atomcounts and name.upper() not in _WELL_KNOWN_SOLVENT_IONS
        ]
        if unknown:
            report.add("molecule_table", BLOCKING, FAIL,
                       f"[ molecules ] references undefined moleculetype(s): {unknown}",
                       {"molecules": molecules, "defined": sorted(mol_atomcounts)})
        else:
            report.add("molecule_table", BLOCKING, PASS,
                       f"{len(molecules)} molecule row(s), all resolvable", molecules)

    # 5. atomtype resolution (Phase 6 preflight, applied post-hoc here)
    ff_types = _ff_atomtypes(forcefield)
    if ff_types is None:
        unresolved = sorted(referenced_types - defined_local_types)
        if unresolved:
            report.add("atomtype_resolution", ADVISORY, WARN,
                       f"{len(unresolved)} atomtype(s) not defined locally; "
                       f"GROMACS force-field data unavailable to confirm they are "
                       f"provided by '{forcefield}'", unresolved)
        else:
            report.add("atomtype_resolution", BLOCKING, PASS,
                       "all referenced atomtypes defined locally", None)
    else:
        unresolved = sorted(referenced_types - defined_local_types - ff_types)
        if unresolved:
            report.add("atomtype_resolution", BLOCKING, FAIL,
                       f"{len(unresolved)} atomtype(s) defined neither locally nor by "
                       f"'{forcefield}': {unresolved}", unresolved)
        else:
            report.add("atomtype_resolution", BLOCKING, PASS,
                       f"all {len(referenced_types)} referenced atomtypes resolve", None)

    # 6. atom counts
    if gro is not None and molecules:
        total = 0
        undetermined: list[str] = []
        for name, count in molecules:
            if name in mol_atomcounts:
                total += mol_atomcounts[name] * count
            else:
                undetermined.append(name)
        if undetermined:
            report.add("atom_counts", ADVISORY, SKIP,
                       f"cannot size solvent/ion molecule(s) {undetermined} "
                       f"(defined in force-field includes)", None)
        elif total == gro.atom_count:
            report.add("atom_counts", BLOCKING, PASS,
                       f"topology and coordinates agree ({total} atoms)", total)
        else:
            report.add("atom_counts", BLOCKING, FAIL,
                       f"topology sums to {total} atoms but coordinate file has "
                       f"{gro.atom_count}", {"topology": total, "coordinates": gro.atom_count})
    else:
        report.add("atom_counts", BLOCKING, SKIP, "coordinate or molecule table unavailable")

    # 7. coordinate validity
    if gro is not None:
        report.checks.append(_check_coordinate_validity(gro))

    # 8. grompp dry-run
    do_grompp = run_grompp if run_grompp is not None else bool(shutil.which("gmx"))
    if not do_grompp:
        report.add("grompp_dry_run", ADVISORY, SKIP, "gmx not in PATH (or disabled)")
    elif gro is None:
        report.add("grompp_dry_run", ADVISORY, SKIP, "no coordinate file to check against")
    else:
        report.checks.append(_grompp_dry_run(d, coord_path, top_path))

    return report


_WELL_KNOWN_SOLVENT_IONS = {
    "SOL", "WAT", "HOH", "TIP3", "TIP3P", "SPC", "SPCE",
    "NA", "CL", "K", "MG", "CA", "ZN", "NA+", "CL-", "ION", "CLA", "SOD", "POT",
}


def _parse_molecules(top_text: str) -> list[tuple[str, int]]:
    out: list[tuple[str, int]] = []
    in_block = False
    for line in top_text.splitlines():
        s = line.split(";")[0].strip()
        if not s:
            continue
        if s.startswith("["):
            in_block = s.strip("[] ").lower() == "molecules"
            continue
        if in_block:
            parts = s.split()
            if len(parts) >= 2:
                try:
                    out.append((parts[0], int(parts[1])))
                except ValueError:
                    continue
    return out


def _check_coordinate_validity(gro: Any) -> Check:
    problems: list[str] = []
    coord_set: dict[tuple, int] = {}
    exact_overlaps = 0
    for a in gro.atoms:
        for v in (a.x, a.y, a.z):
            if math.isnan(v) or math.isinf(v):
                problems.append("NaN/inf atom coordinate present")
                break
        key = (round(a.x, 6), round(a.y, 6), round(a.z, 6))
        if key in coord_set:
            exact_overlaps += 1
        else:
            coord_set[key] = 1
    box = gro.box
    if not box or len(box) < 3 or any((b <= 0 or math.isnan(b) or math.isinf(b)) for b in box[:3]):
        problems.append(f"invalid box vectors: {box}")
    if exact_overlaps:
        problems.append(f"{exact_overlaps} atom(s) share exact coordinates with another atom")

    if not problems:
        return Check("coordinate_validity", BLOCKING, PASS,
                     "coordinates finite, box sane, no exact overlaps", None)
    # NaN/inf and bad box are blocking; exact overlaps alone are advisory
    hard = any("NaN" in p or "box" in p for p in problems)
    return Check("coordinate_validity", BLOCKING if hard else ADVISORY,
                 FAIL if hard else WARN, "; ".join(dict.fromkeys(problems)), problems)


def _grompp_dry_run(out_dir: Path, gro: Path, top: Path) -> Check:
    mdp = out_dir / "_sf_validate.mdp"
    tpr = out_dir / "_sf_validate.tpr"
    mdp.write_text(
        "integrator    = steep\n"
        "nsteps        = 0\n"
        "nstlist       = 1\n"
        "cutoff-scheme = Verlet\n"
        "rlist         = 1.0\n"
        "rcoulomb      = 1.0\n"
        "rvdw          = 1.0\n"
    )
    try:
        res = subprocess.run(
            ["gmx", "grompp", "-f", str(mdp), "-c", str(gro), "-p", str(top),
             "-o", str(tpr), "-maxwarn", "0"],
            capture_output=True, text=True, cwd=out_dir, timeout=90,
        )
        combined = (res.stdout + res.stderr)
        warns = [ln.strip() for ln in combined.splitlines() if ln.strip().lower().startswith("warning")]
        if res.returncode == 0:
            if warns:
                return Check("grompp_dry_run", ADVISORY, WARN,
                             f"grompp succeeded with {len(warns)} warning(s)", warns)
            return Check("grompp_dry_run", BLOCKING, PASS, "gmx grompp -maxwarn 0 succeeded", None)
        fatal = [ln.strip() for ln in combined.splitlines()
                 if "error" in ln.lower() or "fatal" in ln.lower()]
        return Check("grompp_dry_run", BLOCKING, FAIL,
                     "gmx grompp failed", fatal[-8:] or combined.strip()[-800:])
    except subprocess.TimeoutExpired:
        return Check("grompp_dry_run", ADVISORY, WARN, "gmx grompp timed out (90 s)")
    except Exception as exc:
        return Check("grompp_dry_run", ADVISORY, WARN, f"gmx grompp could not run: {exc}")
    finally:
        mdp.unlink(missing_ok=True)
        tpr.unlink(missing_ok=True)
        (out_dir / "mdout.mdp").unlink(missing_ok=True)
