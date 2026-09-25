"""Declarative review profiles (``simforge/review-profile/v1``).

A profile is *data*: which existing observable instances would be useful to
review, an optional display preference, component requirements, and
optional requests conditioned on ACTIVE annotations.  It never contains
commands, callbacks or expressions, and it never carries policy intent:
its display preference is judged by the Phase 5 policy exactly like an
unflagged request (``IntentSource.AUTO``) — choosing a profile is not an
explicit authorisation of interpretation-changing transformations.

Schema::

    schema: simforge/review-profile/v1
    id: protein-ligand                 # [a-z0-9][a-z0-9-]*
    version: 1                         # integer, bumped when content changes
    description: free text             # not part of the identity
    extends: [general]                 # optional; profile ids, resolved as a DAG
    display: raw                       # optional; raw | whole | center:<Group> | fit:<Group>
    requirements:
      components: [receptor, ligand]   # must be RESOLVED components
      absent_components: [membrane]    # must NOT be resolved
    observables:                       # always requested
      - observable: rg
        parameters: {selection: "component:receptor"}
    optional_observables:              # requested only when their condition holds
      - when: {annotation: selected_subunit}          # a named ACTIVE annotation
        observable: rg
        parameters: {selection: "annotation:selected_subunit"}
      - when: {annotation_kind: binding_site}         # exactly one ACTIVE of that kind
        observable: min-distance
        parameters: {selection_a: "component:ligand", selection_b: "annotation:{match}"}

``annotation:{match}`` is the only substitution, allowed only with
``annotation_kind`` — it names the single matching annotation id.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import yaml

from analysis.review.request import ObservableRequest

PROFILE_SCHEMA = "simforge/review-profile/v1"
MATCH_TOKEN = "annotation:{match}"
_ID = re.compile(r"^[a-z0-9][a-z0-9-]{0,63}$")
_TOP_KEYS = {"schema", "id", "version", "description", "extends", "display", "requirements",
             "observables", "optional_observables"}
_REQ_KEYS = {"components", "absent_components"}
_OBS_KEYS = {"observable", "parameters"}
_OPT_KEYS = {"when", "observable", "parameters"}
_WHEN_KEYS = {"annotation", "annotation_kind"}


class ProfileError(ValueError):
    """A profile that is not valid declarative data of this schema."""


@dataclass(frozen=True)
class OptionalRequest:
    request: ObservableRequest
    annotation: Optional[str] = None          # a named annotation id …
    annotation_kind: Optional[str] = None     # … or exactly one ACTIVE of this kind

    def to_dict(self) -> dict:
        when = ({"annotation": self.annotation} if self.annotation
                else {"annotation_kind": self.annotation_kind})
        return {"when": when, "observable": self.request.observable,
                "parameters": self.request.params}


@dataclass
class ReviewProfile:
    id: str
    version: int
    description: str = ""
    extends: list[str] = field(default_factory=list)
    display: Optional[str] = None
    components: list[str] = field(default_factory=list)
    absent_components: list[str] = field(default_factory=list)
    observables: list[ObservableRequest] = field(default_factory=list)
    optional: list[OptionalRequest] = field(default_factory=list)
    source: str = "file"                      # "builtin" | "file" — provenance, not identity
    composition: list[str] = field(default_factory=list)   # resolved ids, parents first
    duplicates_dropped: list[str] = field(default_factory=list)

    # ── identity ───────────────────────────────────────────────────────────
    def definition_content(self) -> dict:
        """Everything that defines what the profile requests (no description,
        no source path, no file formatting)."""
        return {"schema": PROFILE_SCHEMA, "id": self.id, "version": self.version,
                "composition": self.composition or [self.id],
                "display": self.display,
                "requirements": {"components": sorted(self.components),
                                 "absent_components": sorted(self.absent_components)},
                "observables": [[o.observable, sorted(o.params.items())] for o in self.observables],
                "optional": [[o.annotation, o.annotation_kind, o.request.observable,
                              sorted(o.request.params.items())] for o in self.optional]}

    @property
    def definition_identity(self) -> str:
        from analysis.campaign.results import definition_token
        return definition_token(self.definition_content())

    def to_dict(self) -> dict:
        return {"schema": PROFILE_SCHEMA, "id": self.id, "version": self.version,
                "description": self.description, "extends": self.extends,
                "display": self.display,
                "requirements": {"components": self.components,
                                 "absent_components": self.absent_components},
                "observables": [{"observable": o.observable, "parameters": o.params}
                                for o in self.observables],
                "optional_observables": [o.to_dict() for o in self.optional],
                "source": self.source, "composition": self.composition,
                "duplicates_dropped": self.duplicates_dropped,
                "definition_identity": self.definition_identity}


# ═══════════════════════════════════════════════════════════════════════════════
# Parsing / validation (strict)
# ═══════════════════════════════════════════════════════════════════════════════

def _unknown(where: str, d: dict, allowed: set) -> None:
    extra = sorted(set(d) - allowed)
    if extra:
        raise ProfileError(f"{where}: unknown key(s) {extra}; allowed: {sorted(allowed)}")


def _str_list(where: str, v) -> list[str]:
    if v is None:
        return []
    if not isinstance(v, list) or not all(isinstance(x, str) for x in v):
        raise ProfileError(f"{where} must be a list of strings")
    return list(v)


def _request(where: str, d: Any, *, allow_match: bool) -> ObservableRequest:
    from analysis.campaign.observables import registry
    from analysis.campaign.observables.selections import SelectionRef
    from analysis.review.request import RequestError, request_problems
    if not isinstance(d, dict):
        raise ProfileError(f"{where} must be a mapping")
    _unknown(where, d, _OPT_KEYS if allow_match else _OBS_KEYS)
    oid = d.get("observable")
    if not isinstance(oid, str):
        raise ProfileError(f"{where}: 'observable' is required")
    params = d.get("parameters") or {}
    if not isinstance(params, dict):
        raise ProfileError(f"{where}: 'parameters' must be a mapping")
    try:
        req = ObservableRequest.make(oid, params)
    except RequestError as exc:
        raise ProfileError(f"{where}: {exc}") from exc
    problems = request_problems(req)
    if problems:
        raise ProfileError(f"{where}: " + "; ".join(problems))
    for name in registry.get(oid).parameters_schema():
        if name not in req.params:
            raise ProfileError(f"{where}: {oid} needs parameter {name!r}")
        value = str(req.params[name])
        if value == MATCH_TOKEN:
            if not allow_match:
                raise ProfileError(f"{where}: {MATCH_TOKEN!r} is only valid with annotation_kind")
            continue
        try:
            SelectionRef.parse(value)
        except ValueError as exc:
            raise ProfileError(f"{where}: {exc}") from exc
    for name, value in req.params.items():
        if str(value) == MATCH_TOKEN and name not in registry.get(oid).parameters_schema():
            raise ProfileError(f"{where}: {MATCH_TOKEN!r} used in a non-selection parameter")
    return req


def parse_profile(data: Any, *, source: str = "file") -> ReviewProfile:
    """Validate one profile document (not yet composed)."""
    from analysis.campaign.models import ComponentType
    from analysis.review.request import DisplayRequest, RequestError
    if not isinstance(data, dict):
        raise ProfileError("a profile must be a mapping")
    _unknown("profile", data, _TOP_KEYS)
    if data.get("schema") != PROFILE_SCHEMA:
        raise ProfileError(f"profile schema must be {PROFILE_SCHEMA!r}, got {data.get('schema')!r}")
    pid = data.get("id")
    if not isinstance(pid, str) or not _ID.match(pid):
        raise ProfileError(f"profile id {pid!r} must match {_ID.pattern}")
    version = data.get("version")
    if not isinstance(version, int) or isinstance(version, bool) or version < 1:
        raise ProfileError(f"{pid}: version must be a positive integer")
    display = data.get("display")
    if display is not None:
        try:
            DisplayRequest.parse(str(display))
        except RequestError as exc:
            raise ProfileError(f"{pid}: display: {exc}") from exc
    req = data.get("requirements") or {}
    if not isinstance(req, dict):
        raise ProfileError(f"{pid}: requirements must be a mapping")
    _unknown(f"{pid}.requirements", req, _REQ_KEYS)
    known = {v for k, v in vars(ComponentType).items() if k.isupper() and isinstance(v, str)}
    comps = _str_list(f"{pid}.requirements.components", req.get("components"))
    absent = _str_list(f"{pid}.requirements.absent_components", req.get("absent_components"))
    for c in comps + absent:
        if c not in known:
            raise ProfileError(f"{pid}: unknown component type {c!r}; one of {sorted(known)}")
    obs = [_request(f"{pid}.observables[{i}]", d, allow_match=False)
           for i, d in enumerate(data.get("observables") or [])]
    optional = []
    for i, d in enumerate(data.get("optional_observables") or []):
        where = f"{pid}.optional_observables[{i}]"
        r = _request(where, d, allow_match=True)
        when = d.get("when")
        if not isinstance(when, dict) or len(when) != 1:
            raise ProfileError(f"{where}: 'when' needs exactly one of {sorted(_WHEN_KEYS)}")
        _unknown(f"{where}.when", when, _WHEN_KEYS)
        (key, value), = when.items()
        if not isinstance(value, str) or not value:
            raise ProfileError(f"{where}: when.{key} must be a non-empty string")
        uses_match = any(str(v) == MATCH_TOKEN for v in r.params.values())
        if key == "annotation_kind" and not uses_match:
            raise ProfileError(f"{where}: annotation_kind requests must reference "
                               f"{MATCH_TOKEN!r} (the single matching annotation)")
        if key == "annotation" and uses_match:
            raise ProfileError(f"{where}: {MATCH_TOKEN!r} is only valid with annotation_kind")
        optional.append(OptionalRequest(r, annotation=value if key == "annotation" else None,
                                        annotation_kind=value if key == "annotation_kind" else None))
    extends = _str_list(f"{pid}.extends", data.get("extends"))
    if pid in extends:
        raise ProfileError(f"{pid} extends itself")
    return ReviewProfile(id=pid, version=version, description=str(data.get("description") or ""),
                         extends=extends, display=str(display) if display is not None else None,
                         components=comps, absent_components=absent, observables=obs,
                         optional=optional, source=source)


def load_profile_file(path: str | Path) -> ReviewProfile:
    p = Path(path)
    try:
        data = yaml.safe_load(p.read_text())
    except (OSError, yaml.YAMLError) as exc:
        raise ProfileError(f"cannot read profile {p}: {exc}") from exc
    return parse_profile(data, source="file")
