"""Resolve a profile against one discovered system → standard ReviewRequest.

Uses only existing resolved scientific state:

* requirements: ``MolecularComponent.classification_state`` (RESOLVED counts;
  ambiguous / review-required components make the profile REVIEW_REQUIRED;
  nothing is inferred from names or paths);
* optional requests: ``AnnotationRecord.state`` — only ACTIVE annotations
  count; PROPOSED ones are never activated; a kind-only condition with more
  than one ACTIVE candidate is REVIEW_REQUIRED (no ranking, no choice).

The output is ordinary :class:`ObservableRequest` instances; every one still
goes through ``ObservableSpec.applicability``, the policy and Phase 8 reuse
in generic preparation.  Nothing here executes anything.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from analysis.campaign.models import AnnotationState, ClassificationState
from analysis.review.profiles.model import MATCH_TOKEN, ProfileError, ReviewProfile
from analysis.review.request import DisplayRequest, ObservableRequest, ReviewRequest


class ProfileStatus:
    APPLICABLE = "applicable"
    PARTIALLY_APPLICABLE = "partially_applicable"
    NOT_APPLICABLE = "not_applicable"
    REVIEW_REQUIRED = ClassificationState.REVIEW_REQUIRED


class OptionalOutcome:
    INCLUDED = "included"
    UNAVAILABLE = "unavailable"               # condition not met — simply not requested
    REVIEW_REQUIRED = ClassificationState.REVIEW_REQUIRED


@dataclass
class ProfileResolution:
    profile: ReviewProfile
    status: str
    requirements: list[dict] = field(default_factory=list)
    requests: list[ObservableRequest] = field(default_factory=list)   # profile-derived, in order
    optional: list[dict] = field(default_factory=list)
    reasons: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {"profile_id": self.profile.id, "status": self.status,
                "requirements": self.requirements,
                "requests": [r.instance_id for r in self.requests],
                "optional": self.optional, "reasons": self.reasons}


def _requirement(rec, ctype: str, *, absent: bool) -> dict:
    comps = [c for c in rec.components if c.component_type == ctype]
    resolved = any(c.classification_state == ClassificationState.RESOLVED for c in comps)
    states = sorted({c.classification_state for c in comps}) or ["absent"]
    if absent:
        ok = not resolved
        return {"component": ctype, "requirement": "absent", "states": states, "ok": ok,
                "status": ProfileStatus.APPLICABLE if ok else ProfileStatus.NOT_APPLICABLE,
                "reason": "" if ok else f"a resolved {ctype} component is present"}
    if resolved:
        return {"component": ctype, "requirement": "resolved", "states": states, "ok": True,
                "status": ProfileStatus.APPLICABLE, "reason": ""}
    unsure = [s for s in states if s in (ClassificationState.AMBIGUOUS,
                                         ClassificationState.REVIEW_REQUIRED)]
    return {"component": ctype, "requirement": "resolved", "states": states, "ok": False,
            "status": ProfileStatus.REVIEW_REQUIRED if unsure else ProfileStatus.NOT_APPLICABLE,
            "reason": (f"{ctype} component is {unsure[0]}" if unsure
                       else f"no resolved {ctype} component")}


def _optional(rec, opt) -> tuple[dict, Optional[ObservableRequest]]:
    anns = rec.annotations
    template = opt.request
    out = {"request": template.instance_id,
           "condition": ({"annotation": opt.annotation} if opt.annotation
                         else {"annotation_kind": opt.annotation_kind})}
    if opt.annotation:
        rec_a = next((a for a in anns if a.annotation_id == opt.annotation), None)
        if rec_a is None:
            return {**out, "outcome": OptionalOutcome.UNAVAILABLE,
                    "reason": f"annotation '{opt.annotation}' is not declared for this system"}, None
        if rec_a.state != AnnotationState.ACTIVE:
            return {**out, "outcome": OptionalOutcome.UNAVAILABLE,
                    "reason": f"annotation '{opt.annotation}' is {rec_a.state}, not ACTIVE"}, None
        return {**out, "outcome": OptionalOutcome.INCLUDED, "annotation": opt.annotation,
                "instance_id": template.instance_id}, template
    kind = opt.annotation_kind
    active = sorted(a.annotation_id for a in anns
                    if a.kind == kind and a.state == AnnotationState.ACTIVE)
    inactive = sorted(f"{a.annotation_id} ({a.state})" for a in anns
                      if a.kind == kind and a.state != AnnotationState.ACTIVE)
    if not active:
        why = f"no ACTIVE {kind} annotation"
        if inactive:
            why += f" (present but not active: {', '.join(inactive)})"
        return {**out, "outcome": OptionalOutcome.UNAVAILABLE, "reason": why}, None
    if len(active) > 1:
        return {**out, "outcome": OptionalOutcome.REVIEW_REQUIRED, "candidates": active,
                "reason": f"{len(active)} ACTIVE {kind} annotations {active}; the profile must "
                          f"name one (no automatic choice)"}, None
    match = active[0]
    params = {k: (f"annotation:{match}" if str(v) == MATCH_TOKEN else v)
              for k, v in template.params.items()}
    req = ObservableRequest.make(template.observable, params)
    return {**out, "outcome": OptionalOutcome.INCLUDED, "annotation": match,
            "instance_id": req.instance_id}, req


def resolve_profile(profile: ReviewProfile, rec) -> ProfileResolution:
    reqs = ([_requirement(rec, c, absent=False) for c in profile.components]
            + [_requirement(rec, c, absent=True) for c in profile.absent_components])
    failed = [r for r in reqs if not r["ok"]]
    if any(r["status"] == ProfileStatus.NOT_APPLICABLE for r in failed):
        status = ProfileStatus.NOT_APPLICABLE
    elif failed:
        status = ProfileStatus.REVIEW_REQUIRED
    else:
        status = ProfileStatus.APPLICABLE
    requests = list(profile.observables)
    optional = []
    for opt in profile.optional:
        outcome, req = _optional(rec, opt)
        optional.append(outcome)
        if req is not None and req not in requests:
            requests.append(req)
    return ProfileResolution(profile, status, reqs, requests, optional,
                             [r["reason"] for r in failed])


def merge_request(resolution: ProfileResolution, user: ReviewRequest, hide=()) -> tuple:
    """Profile requests + explicit user requests (user wins on display).

    ``hide`` entries are instance ids or observable ids; an entry that
    matches nothing is an error rather than silently ignored."""
    hide = list(hide or [])
    matched = {h: False for h in hide}

    def hidden(r):
        hit = [h for h in hide if h in (r.instance_id, r.observable)]
        for h in hit:
            matched[h] = True
        return bool(hit)
    from_profile = [r for r in resolution.requests if not hidden(r)]
    added = [r for r in user.observables if r not in from_profile]
    unused = [h for h, ok in matched.items() if not ok]
    if unused:
        raise ProfileError(f"--hide {unused} matches no request of profile "
                           f"{resolution.profile.id!r}")
    if user.display is not None:
        display, source = user.display, "user"
    elif resolution.profile.display is not None:
        display, source = DisplayRequest.parse(resolution.profile.display), "profile"
    else:
        display, source = DisplayRequest(), "default"
    overrides = {"added": [r.instance_id for r in added],
                 "hidden": [r.instance_id for r in resolution.requests if r not in from_profile],
                 "display_source": source}
    return ReviewRequest(from_profile + added, display), overrides


def final_status(resolution: ProfileResolution, entries) -> str:
    """Requirement status, refined by what preparation actually delivered."""
    if resolution.status != ProfileStatus.APPLICABLE:
        return resolution.status
    profile_ids = {r.instance_id for r in resolution.requests}
    missing = [e for e in entries if e.instance_id in profile_ids
               and e.state not in ("available", "planned")]
    review = any(o["outcome"] == OptionalOutcome.REVIEW_REQUIRED for o in resolution.optional)
    return ProfileStatus.PARTIALLY_APPLICABLE if (missing or review) else ProfileStatus.APPLICABLE
