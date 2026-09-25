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
    if ref.format == StorageFormat.NPY:
        import numpy as np
        arr = np.load(Path(ref.path), allow_pickle=False)
        if arr.ndim == 1 and ref.column in (None, 0):
            return arr.tolist()
        if arr.ndim == 2 and ref.column is not None and 0 <= ref.column < arr.shape[1]:
            return arr[:, ref.column].tolist()
        raise ValueError(f"{Path(ref.path).name}: cannot address column {ref.column} of a "
                         f"{arr.ndim}-D array")
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


_GRACE_MARKUP = {"\\S2\\N": "^2", "\\S3\\N": "^3", "\\S-1\\N": "^-1"}


def normalize_unit(unit: Optional[str]) -> Optional[str]:
    """Normalise xmgrace markup in a unit label (``nm\\S2\\N`` → ``nm^2``)."""
    if not unit:
        return None
    for k, v in _GRACE_MARKUP.items():
        unit = unit.replace(k, v)
    return unit.strip() or None


def import_external_series(path: str | Path, *, quantity: str, unit: str, time_unit: str,
                           value_column: int = 1, fmt: Optional[str] = None,
                           description: str = "") -> "ResultArray":
    """Explicitly declare an existing XVG/CSV series as a ResultArray.

    Meaning comes only from the caller's declaration — never from the filename.
    A unit/time unit that contradicts the file header is refused.  The result is
    marked externally supplied, is fingerprinted, and is never reused
    automatically (compatibility is a separate, later decision).
    """
    from analysis.campaign.fingerprint import fingerprint_file
    from analysis.campaign.models import Axis, AxisKind, MissingSemantics, ResultArray
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(p)
    fmt = fmt or p.suffix.lower().lstrip(".")
    if fmt not in (StorageFormat.XVG, StorageFormat.CSV):
        raise ValueError(f"unsupported format {fmt!r} (xvg or csv)")
    header = {"x": None, "y": None}
    if fmt == StorageFormat.XVG:
        x_unit, y_unit = xvg_units(p)
        header = {"x": x_unit, "y": normalize_unit(y_unit)}
        if header["x"] and header["x"] != time_unit:
            raise ValueError(f"declared time unit {time_unit!r} contradicts the header ({header['x']!r})")
        if header["y"] and header["y"] != unit:
            raise ValueError(f"declared unit {unit!r} contradicts the header ({header['y']!r})")
    values = StorageRef(path=str(p.resolve()), format=fmt, column=value_column,
                        fingerprint=fingerprint_file(p))
    read_column(values)                               # must be readable as declared
    return ResultArray(
        name=quantity, quantity=quantity, unit=unit,
        axes=[Axis(name="time", kind=AxisKind.TIME, unit=time_unit,
                   values_ref=StorageRef(path=str(p.resolve()), format=fmt, column=0),
                   attrs={"unit_source": "xvg header" if header["x"] else "declared"})],
        storage=values, missing=MissingSemantics.UNSPECIFIED,
        attrs={"externally_supplied": True, "declared_by": "caller",
               "unit_source": "xvg header" if header["y"] else "declared",
               "description": description})
