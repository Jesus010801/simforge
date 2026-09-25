"""What a review session asks for: observable instances and a display request.

Observable syntax (CLI ``--show``, repeatable, comma-separated at top level)::

    rmsd-receptor
    rg(selection=component:receptor)
    com-distance(selection_a=component:ligand,selection_b=annotation:selected_subunit)

An *instance* is one observable id with one parameter set, so the same
observable can be requested several times (Rg of the receptor and of a
subunit).  Its id is the canonical call string — sorted ``key=value`` pairs —
so it is stable and human readable.  Values are opaque strings (selection
references never contain ``,``, ``(`` or ``)``); anything richer goes through a
YAML/JSON request file::

    observables:
      - observable: rg
        parameters: {selection: "annotation:selected_subunit"}

Display request (``--display``): ``raw`` (default: the source coordinates,
nothing transformed), ``whole``, ``center:<Group>`` or ``fit:<Group>``.  It is
planned by the Phase 5 policy with purpose DISPLAY; interpretation-changing
operations need ``--display-intent`` (recorded as ``user_flag`` intent for the
display view only — never for the observables).
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import yaml

from analysis.campaign.models import IntentSource, TrajectoryRequirements


class RequestError(ValueError):
    """An observable / display request that cannot be parsed or is unknown."""


# ═══════════════════════════════════════════════════════════════════════════════
# Observable instances
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class ObservableRequest:
    observable: str
    parameters: tuple[tuple[str, Any], ...] = ()

    @property
    def params(self) -> dict[str, Any]:
        return dict(self.parameters)

    @property
    def instance_id(self) -> str:
        if not self.parameters:
            return self.observable
        return f"{self.observable}(" + ",".join(f"{k}={v}" for k, v in self.parameters) + ")"

    def to_dict(self) -> dict:
        return {"instance_id": self.instance_id, "observable": self.observable,
                "parameters": self.params}

    @classmethod
    def make(cls, observable: str, parameters: Optional[dict] = None) -> "ObservableRequest":
        params = parameters or {}
        for k, v in params.items():
            if isinstance(v, (dict, list)):
                raise RequestError(f"{observable}: parameter {k!r} must be a scalar")
        return cls(observable.strip(), tuple(sorted((str(k), v) for k, v in params.items())))

    @classmethod
    def from_dict(cls, d: dict) -> "ObservableRequest":
        return cls.make(d["observable"], d.get("parameters"))


def _split_top_level(text: str) -> list[str]:
    out, depth, cur = [], 0, []
    for ch in text:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth < 0:
                raise RequestError(f"unbalanced ')' in {text!r}")
        if ch == "," and depth == 0:
            out.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    if depth:
        raise RequestError(f"unbalanced '(' in {text!r}")
    out.append("".join(cur))
    return [p.strip() for p in out if p.strip()]


def parse_observable(item: str) -> ObservableRequest:
    item = item.strip()
    if "(" not in item:
        if ")" in item or "=" in item:
            raise RequestError(f"cannot parse observable {item!r}")
        return ObservableRequest.make(item)
    if not item.endswith(")"):
        raise RequestError(f"observable {item!r} must end with ')'")
    oid, _, body = item[:-1].partition("(")
    params: dict[str, str] = {}
    for pair in _split_top_level(body):
        key, sep, value = pair.partition("=")
        if not sep or not key.strip() or not value.strip():
            raise RequestError(f"parameter {pair!r} of {oid!r} must be key=value")
        if key.strip() in params:
            raise RequestError(f"parameter {key.strip()!r} given twice for {oid!r}")
        params[key.strip()] = value.strip()
    return ObservableRequest.make(oid, params)


def parse_show(values: list[str]) -> list[ObservableRequest]:
    """``--show`` values → instances, de-duplicated, first occurrence order."""
    out: list[ObservableRequest] = []
    for v in values or []:
        for item in _split_top_level(v):
            r = parse_observable(item)
            if r not in out:
                out.append(r)
    return out


def load_request_file(path: str | Path) -> list[ObservableRequest]:
    p = Path(path)
    text = p.read_text()
    data = json.loads(text) if p.suffix.lower() == ".json" else yaml.safe_load(text)
    items = (data or {}).get("observables", []) if isinstance(data, dict) else data
    if not isinstance(items, list):
        raise RequestError(f"{p}: expected an 'observables' list")
    return [ObservableRequest.from_dict(i) if isinstance(i, dict) else parse_observable(str(i))
            for i in items]


def request_problems(r: ObservableRequest) -> list[str]:
    """Unknown id, or parameters given to an observable that takes none.
    (Selection parameters are validated by the observable itself; metric
    parameters such as ``probe`` are not in ``parameters_schema``.)"""
    from analysis.campaign.observables import registry
    registry.ensure_loaded()
    if not registry.is_registered(r.observable):
        return [f"'{r.observable}' is not a registered observable (valid: {registry.ids()})"]
    if r.params and not registry.get(r.observable).parameters_schema():
        return [f"'{r.observable}' takes no parameters"]
    return []


# ═══════════════════════════════════════════════════════════════════════════════
# Display request
# ═══════════════════════════════════════════════════════════════════════════════

DISPLAY_MODES = ("raw", "whole", "center", "fit")


@dataclass(frozen=True)
class DisplayRequest:
    mode: str = "raw"
    target: Optional[str] = None           # semantic group for center / fit
    intent: bool = False                   # --display-intent (user_flag)

    @classmethod
    def parse(cls, text: Optional[str], *, intent: bool = False) -> "DisplayRequest":
        text = (text or "raw").strip()
        mode, _, target = text.partition(":")
        if mode not in DISPLAY_MODES:
            raise RequestError(f"--display must be one of raw, whole, center:<Group>, "
                               f"fit:<Group>; got {text!r}")
        if mode in ("center", "fit") and not target:
            raise RequestError(f"--display {mode} needs a semantic group, e.g. {mode}:Receptor")
        if mode in ("raw", "whole") and target:
            raise RequestError(f"--display {mode} takes no group")
        return cls(mode, target or None, intent)

    def requirements(self) -> TrajectoryRequirements:
        """Display modes as typed requests; center/fit keep molecules whole."""
        if self.mode == "raw":
            return TrajectoryRequirements(rationale="review display: source coordinates")
        if self.mode == "whole":
            return TrajectoryRequirements(requires_whole_molecules=True,
                                          rationale="review display: whole molecules")
        if self.mode == "center":
            return TrajectoryRequirements(requires_whole_molecules=True,
                                          centering_target=self.target,
                                          rationale=f"review display: centred on {self.target}")
        return TrajectoryRequirements(requires_whole_molecules=True, fit_selection=self.target,
                                      rationale=f"review display: fitted to {self.target}")

    def policy_intent(self):
        from analysis.campaign.trajectory.policy import PolicyIntent
        if not self.intent:
            return PolicyIntent()
        ops = [op for op, on in (("make_whole", self.mode != "raw"),
                                 ("center", self.mode == "center"),
                                 ("fit", self.mode == "fit")) if on]
        return PolicyIntent(source=IntentSource.USER_FLAG, operations=tuple(ops))

    def to_dict(self) -> dict:
        return {"mode": self.mode, "target": self.target,
                "intent_source": IntentSource.USER_FLAG if self.intent else IntentSource.AUTO}


@dataclass
class ReviewRequest:
    observables: list[ObservableRequest] = field(default_factory=list)
    #: ``None`` = no explicit display choice (a profile's preference, else raw)
    display: Optional[DisplayRequest] = field(default_factory=DisplayRequest)

    def to_dict(self) -> dict:
        return {"observables": [o.to_dict() for o in self.observables],
                "display": self.display.to_dict() if self.display else None}
