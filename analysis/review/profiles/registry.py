"""The one place review profiles are found and composed.

Built-in profiles are YAML resources in ``builtin/``; user profiles are YAML
files.  ``extends`` names profile ids (built-ins, or ids of other profiles
passed in ``extra``) and is resolved as a DAG: parents first in the listed
order, each profile once, cycles refused.  Merging is explicit:

* observables / optional requests are concatenated in that order; a request
  already present is dropped and listed in ``duplicates_dropped``;
* requirements are unioned; a component both required and absent is an error;
* ``display``: the profile's own value, else the single value its parents
  agree on (parents that disagree are an error, never a silent choice).
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import yaml

from analysis.review.profiles.model import (
    ProfileError, ReviewProfile, load_profile_file, parse_profile,
)

BUILTIN_DIR = Path(__file__).resolve().parent / "builtin"


def builtin_profiles(directory: Path = BUILTIN_DIR) -> dict[str, ReviewProfile]:
    out: dict[str, ReviewProfile] = {}
    for f in sorted(directory.glob("*.yaml")):
        p = parse_profile(yaml.safe_load(f.read_text()), source="builtin")
        if p.id in out:
            raise ProfileError(f"duplicate built-in profile id {p.id!r} ({f.name})")
        out[p.id] = p
    return out


def _is_path(ref: str) -> bool:
    return "/" in ref or ref.endswith((".yaml", ".yml")) or Path(ref).is_file()


def get_profile(ref: str, *, extra: Optional[dict[str, ReviewProfile]] = None,
                builtins: Optional[dict[str, ReviewProfile]] = None) -> ReviewProfile:
    """A built-in id or a YAML path → the composed (resolved) profile."""
    lookup = dict(builtins if builtins is not None else builtin_profiles())
    for pid, p in (extra or {}).items():
        if pid in lookup:
            raise ProfileError(f"profile id {pid!r} is defined twice")
        lookup[pid] = p
    if _is_path(ref):
        root = load_profile_file(ref)
    elif ref in lookup:
        root = lookup[ref]
    else:
        raise ProfileError(f"unknown profile {ref!r}; built-in: {sorted(lookup)} "
                           f"(or pass a YAML file)")
    return compose(root, lookup)


def compose(root: ReviewProfile, lookup: dict[str, ReviewProfile]) -> ReviewProfile:
    order: list[ReviewProfile] = []
    seen: set[str] = set()

    def visit(p: ReviewProfile, stack: tuple[str, ...]) -> None:
        if p.id in stack:
            raise ProfileError("profile composition cycle: " + " -> ".join(stack + (p.id,)))
        if p.id in seen:
            return
        for parent_id in p.extends:
            parent = lookup.get(parent_id)
            if parent is None or parent is p:
                raise ProfileError(f"{p.id} extends unknown profile {parent_id!r}")
            visit(parent, stack + (p.id,))
        seen.add(p.id)
        order.append(p)

    visit(root, ())
    obs, opt, dropped = [], [], []
    comps, absent = [], []
    for p in order:
        for o in p.observables:
            if o in obs:
                dropped.append(f"{p.id}: {o.instance_id}")
            else:
                obs.append(o)
        for o in p.optional:
            if o in opt:
                dropped.append(f"{p.id}: optional {o.request.instance_id}")
            else:
                opt.append(o)
        comps += [c for c in p.components if c not in comps]
        absent += [c for c in p.absent_components if c not in absent]
    clash = set(comps) & set(absent)
    if clash:
        raise ProfileError(f"{root.id}: component(s) {sorted(clash)} both required and absent")
    display = root.display
    if display is None:
        parents = {lookup[x].id: compose(lookup[x], lookup).display for x in root.extends}
        values = {v for v in parents.values() if v is not None}
        if len(values) > 1:
            raise ProfileError(f"{root.id}: parents disagree on display {parents}; "
                               f"set 'display' explicitly")
        display = values.pop() if values else None
    return ReviewProfile(id=root.id, version=root.version, description=root.description,
                         extends=list(root.extends), display=display, components=comps,
                         absent_components=absent, observables=obs, optional=opt,
                         source=root.source, composition=[p.id for p in order],
                         duplicates_dropped=dropped)


def list_profiles() -> list[ReviewProfile]:
    b = builtin_profiles()
    return [compose(p, b) for p in b.values()]
