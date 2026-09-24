"""Helpers for the labelled-array result contract (see ``models.ResultArray``).

* :func:`definition_token` — deterministic identity of an observable
  *definition* from explicit evidence (canonical JSON → sha256).  No
  timestamps, no UUIDs; evidence must not contain machine-local output paths.
* :func:`atom_set_hash` / :func:`index_group_evidence` — content identity of a
  named index group, so a token changes when a group's atoms change even if
  its name does not.
* :func:`read_column` — minimal readers for the tabular formats in use
  today (XVG, CSV).  NPY/NPZ references are declarable and validated but not
  read here yet.
* :func:`xvg_units` — units as *declared* in an XVG header.
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Any, Optional

from analysis.campaign.models import StorageFormat, StorageRef

DEFINITION_TOKEN_SCHEMA = "simforge/definition-token/v1"

_UNIT_RE = re.compile(r"\(([^()]+)\)\s*$")


def definition_token(evidence: dict[str, Any]) -> str:
    """sha256 over canonical JSON of *evidence* (sorted keys, no whitespace)."""
    payload = json.dumps({"schema": DEFINITION_TOKEN_SCHEMA, "evidence": evidence},
                         sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(payload.encode()).hexdigest()[:32]


def atom_set_hash(atom_ids) -> str:
    """Order-independent content hash of a set of 1-based atom ids."""
    ids = sorted(set(int(a) for a in atom_ids))
    return hashlib.sha256(",".join(map(str, ids)).encode()).hexdigest()[:24]


def atom_multiset_hash(atom_ids) -> str:
    """Like :func:`atom_set_hash` but keeps repeated ids (sorted, with multiplicity).

    Order within a GROMACS group does not affect centering or fitting (both
    reduce over the group), but a repeated atom is weighted twice — so
    multiplicity, unlike order, can matter.
    """
    ids = sorted(int(a) for a in atom_ids)
    return hashlib.sha256(",".join(map(str, ids)).encode()).hexdigest()[:24]


def index_group_evidence(index_path: str | Path, group_name: str) -> Optional[dict]:
    """``{"name", "n_atoms", "atoms_sha256"}`` for a group, or None if absent/ambiguous.

    Groups that list an atom more than once additionally carry
    ``n_entries`` and ``multiset_sha256``; duplicate-free groups keep exactly
    the original three keys (so existing tokens are unchanged).
    """
    from analysis.campaign.structure.index_groups import resolve_index_group
    try:
        grp = resolve_index_group(index_path, group_name)
    except (OSError, ValueError):
        return None
    if grp is None:
        return None
    ev = {"name": group_name, "n_atoms": len(set(grp.atom_ids)),
          "atoms_sha256": atom_set_hash(grp.atom_ids)}
    if len(grp.atom_ids) != len(set(grp.atom_ids)):
        ev["n_entries"] = len(grp.atom_ids)
        ev["multiset_sha256"] = atom_multiset_hash(grp.atom_ids)
    return ev


def xvg_units(path: str | Path) -> tuple[Optional[str], Optional[str]]:
    """``(x_unit, y_unit)`` as declared by the XVG axis labels; None if undeclared."""
    from runtime.xvg_parser import parse_xvg
    data = parse_xvg(Path(path))
    m = _UNIT_RE.search(data.ylabel or "")
    return data.x_unit, (m.group(1).strip() if m else None)


def read_column(ref: StorageRef) -> list[float]:
    """Read one numeric column referenced by a tabular :class:`StorageRef`.

    Values are returned raw, in file order, in the file's own units — no
    conversion, sorting or de-duplication.
    """
    problems = ref.validate()
    if problems:
        raise ValueError("; ".join(problems))
    if ref.format == StorageFormat.XVG:
        return _read_xvg_column(Path(ref.path), ref.column)
    if ref.format == StorageFormat.CSV:
        return _read_csv_column(Path(ref.path), ref.column)
    raise NotImplementedError(f"reading {ref.format!r} storage is not implemented yet")


def _read_xvg_column(path: Path, column: int) -> list[float]:
    from runtime.xvg_parser import parse_xvg
    data = parse_xvg(path)
    if column == 0:
        return list(data.time_ps)            # raw x values (unit per the header)
    if column - 1 >= len(data.series):
        raise ValueError(f"{path.name}: column {column} out of range "
                         f"(file has {len(data.series) + 1} columns)")
    return list(data.series[column - 1].values)


def _read_csv_column(path: Path, column: int) -> list[float]:
    out: list[float] = []
    with path.open(newline="") as fh:
        for row in csv.reader(line for line in fh if not line.lstrip().startswith("#")):
            if not row:
                continue
            try:
                out.append(float(row[column]))
            except ValueError:
                if out:
                    raise ValueError(f"{path.name}: non-numeric value in column {column}") from None
                continue                          # header row
            except IndexError:
                raise ValueError(f"{path.name}: column {column} out of range") from None
    return out
