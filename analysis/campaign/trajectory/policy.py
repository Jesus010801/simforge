"""Scientific preprocessing policy — is a *requested* coordinate operation
admissible for this purpose, in this system, given this evidence and intent?

    requirements + purpose + system semantics + diagnostics + intent
                               ↓  plan_preprocessing()
                      PreprocessingDecision[]  (+ resolved requirements)

Principles (enforced here, and only here):

* The policy never invents operations: it judges exactly what was requested.
  A blocked request is never silently dropped — the plan becomes
  non-executable and the observable must not run on other coordinates.
* Context comes from ``MolecularComponent.classification_state`` — never from
  the mere existence of an index group, and never from ``Complex`` membership.
* Diagnostics are evidence, not commands: they can only make a decision more
  conservative, never authorise an operation by themselves.
* ``VALID_WITH_INTENT`` is applied only with explicit (non-auto) intent.
* No functional annotations (binding site, pocket, subunit …) are assumed.
* Pure evaluation: no GROMACS, no files written, no ``build_view()``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

from analysis.campaign.models import (
    ClassificationState, ComponentType, DecisionClass, DiagnosticReport, IntentSource,
    MolecularComponent, ObservablePurpose as P, PreprocessingDecision, SystemRecord,
    TrajectoryRequirements,
)

RULE_VERSION = "1"
POLICY_ENGINE = "simforge/preprocessing-policy/v1"

_MEMBRANE_FRAME_PURPOSES = (P.MEMBRANE_FRAME_PROPERTY, P.DENSITY_PROFILE, P.SOLVENT_OCCUPANCY)
_TIMELINE_INVALID = ("non_monotonic_time", "time_index_invalid", "empty_trajectory")
_PBC_PAIR_SUFFIXES = ("_receptor_periodic_separation", "_receptor_image_discrepancy")

#: semantic group name -> component type it stands for
_GROUP_COMPONENT = {
    "Receptor": ComponentType.RECEPTOR, "Peptide": ComponentType.PEPTIDE,
    "ProteinPartner": ComponentType.PROTEIN_PARTNER, "Complex": ComponentType.COMPLEX,
    "Ligand": ComponentType.LIGAND, "Cofactor": ComponentType.COFACTOR,
    "NucleicAcid": ComponentType.NUCLEIC_ACID, "Membrane": ComponentType.MEMBRANE,
    "Water": ComponentType.WATER, "Ions": ComponentType.IONS, "System": "system",
}


# ═══════════════════════════════════════════════════════════════════════════════
# Context and intent
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class PolicyIntent:
    """How intent entered.  ``operations=None`` with an explicit source means
    'every requested operation'; AUTO authorises nothing."""
    source: str = IntentSource.AUTO
    operations: Optional[tuple[str, ...]] = None

    def authorizes(self, operation: str) -> bool:
        return self.source in IntentSource.EXPLICIT and (
            self.operations is None or operation in self.operations)

    def to_dict(self) -> dict:
        return {"source": self.source,
                "operations": list(self.operations) if self.operations is not None else None}


@dataclass
class PolicyContext:
    components: list[MolecularComponent] = field(default_factory=list)
    topology_path: Optional[str] = None
    structure_path: Optional[str] = None
    index_path: Optional[str] = None
    diagnostics: Optional[DiagnosticReport] = None

    @classmethod
    def from_system(cls, rec: SystemRecord, *, diagnostics: Optional[DiagnosticReport] = None
                    ) -> "PolicyContext":
        idx = rec.semantic_index
        return cls(components=list(rec.components), topology_path=rec.topology_path,
                   structure_path=rec.structure_path,
                   index_path=idx.path if idx is not None and idx.path else None,
                   diagnostics=diagnostics)

    # ── semantics (component STATE, not group existence) ───────────────────
    def components_of(self, ctype: str) -> list[MolecularComponent]:
        return [c for c in self.components if c.component_type == ctype]

    def state(self, ctype: str) -> str:
        if ctype == "system":
            return ClassificationState.RESOLVED          # a structural primitive
        comps = self.components_of(ctype)
        if not comps:
            return "absent"
        if len(comps) > 1:
            return ClassificationState.AMBIGUOUS
        return comps[0].classification_state

    def resolved(self, ctype: str) -> bool:
        return self.state(ctype) == ClassificationState.RESOLVED

    def group_names(self) -> Optional[set[str]]:
        if not self.index_path or not Path(self.index_path).is_file():
            return None
        from analysis.campaign.structure.index_groups import parse_index_groups
        try:
            return {g.name for g in parse_index_groups(self.index_path)}
        except ValueError:
            return None

    def group_evidence(self, group: str) -> Optional[dict]:
        if not self.index_path:
            return None
        from analysis.campaign.results import index_group_evidence
        return index_group_evidence(self.index_path, group)

    def component_summary(self, ctype: str) -> dict:
        comps = self.components_of(ctype)
        out = {"state": self.state(ctype)}
        if len(comps) == 1 and comps[0].atom_ids:
            from analysis.campaign.results import atom_set_hash
            out["atoms_sha256"] = atom_set_hash(comps[0].atom_ids)
        return out

    def connectivity_source(self) -> Optional[str]:
        t = self.topology_path
        return "tpr" if t and t.lower().endswith(".tpr") else None

    def fit_reference(self) -> dict:
        """The -s file build_view uses for center/fit (structure first — recorded, not changed)."""
        ref = self.structure_path or self.topology_path
        kind = Path(ref).suffix.lower().lstrip(".") if ref else None
        alternatives = sorted({Path(p).suffix.lower().lstrip(".") for p in
                               (self.structure_path, self.topology_path) if p} - {kind})
        return {"reference_file_type": kind, "alternative_reference_types": alternatives,
                "ambiguous": kind != "tpr" and "tpr" in alternatives}

    # ── diagnostics as evidence ────────────────────────────────────────────
    def diagnostics_available(self) -> bool:
        if self.diagnostics is None:
            return False
        return any(r.detector_id == "timeline_integrity" and r.status == "ran"
                   for r in self.diagnostics.detectors)

    def evidence(self, pred: Callable) -> list[dict]:
        if self.diagnostics is None:
            return []
        out = []
        for d in self.diagnostics.diagnostics:
            if pred(d):
                out.append({
                    "code": d.code, "severity": d.severity, "scope": d.scope,
                    "interpretation": d.interpretation, "confidence": d.confidence,
                    "detector": d.detector,
                    "detector_identity": d.provenance.get("detector_identity"),
                    "centre_convention": d.provenance.get("centre_convention"),
                })
        return out


def _timeline_invalid(ctx):
    return ctx.evidence(lambda d: d.code in _TIMELINE_INVALID)


def _timeline_review(ctx):
    return ctx.evidence(lambda d: d.code == "duplicate_timestamps")


def _pbc_pair_evidence(ctx):
    return ctx.evidence(lambda d: any(d.code.endswith(s) for s in _PBC_PAIR_SUFFIXES))


def _wrap_evidence(ctx):
    return ctx.evidence(lambda d: d.code in ("component_periodic_wrap",
                                             "system_periodic_translation"))


# ═══════════════════════════════════════════════════════════════════════════════
# Rules — one function per operation; each returns a Verdict
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class Verdict:
    classification: str
    rule_id: str
    reason: str
    context: dict = field(default_factory=dict)
    evidence: list[dict] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)


V, VWI, REF, UNS = (DecisionClass.VALID, DecisionClass.VALID_WITH_INTENT,
                    DecisionClass.REFUSED, DecisionClass.UNSUPPORTED)


def rule_make_whole(req, purpose, ctx) -> Verdict:
    src = ctx.connectivity_source()
    context = {"connectivity_source": src, "topology_file_type":
               Path(ctx.topology_path).suffix.lower().lstrip(".") if ctx.topology_path else None}
    if src == "tpr":
        return Verdict(V, "make_whole.tpr_connectivity",
                       "molecules are reassembled from .tpr bond connectivity; no frame change",
                       context)
    return Verdict(UNS, "make_whole.no_connectivity",
                   "no .tpr — GROMACS cannot reconstruct molecular connectivity reliably", context)


def rule_nojump(req, purpose, ctx) -> Verdict:
    bad = _timeline_invalid(ctx)
    if bad:
        return Verdict(REF, "nojump.timeline_invalid",
                       "nojump depends on sequential frame history; the timeline is invalid or "
                       "non-monotonic — resolve segments first", evidence=bad)
    wraps = _wrap_evidence(ctx)
    if purpose in (P.MEMBRANE_FRAME_PROPERTY, P.DENSITY_PROFILE):
        return Verdict(REF, "nojump.bounded_frame_required",
                       "globally unwrapped coordinates leave the periodic reference box that "
                       "membrane-frame / spatial-binning quantities are defined in",
                       evidence=wraps)
    if purpose == P.SOLVENT_OCCUPANCY:
        return Verdict(REF, "nojump.bounded_frame_required",
                       "unwrapped solvent leaves the region occupancy is counted in",
                       evidence=wraps)
    if purpose == P.DIFFUSION:
        if req.centering_target:
            return Verdict(REF, "nojump.rewrapped_by_center",
                           "the requested centring step re-wraps molecules (-pbc mol) after "
                           "nojump, destroying the unwrapped history diffusion needs")
        review = _timeline_review(ctx)
        if review:
            return Verdict(VWI, "nojump.timeline_review",
                           "duplicate timestamps present — confirm segment handling before "
                           "unwrapping", evidence=review)
        if not ctx.diagnostics_available():
            return Verdict(VWI, "nojump.timeline_unverified",
                           "timeline integrity has not been diagnosed for this trajectory")
        return Verdict(V, "nojump.diffusion",
                       "diffusion requires unwrapped, sequential positions; timeline verified",
                       evidence=wraps)
    return Verdict(VWI, "nojump.representation_choice",
                   f"unwrapping is not required for '{purpose}' and changes the coordinate "
                   f"representation — explicit intent needed", evidence=wraps)


def _target_check(group, ctx) -> Optional[Verdict]:
    names = ctx.group_names()
    if names is None or group not in names:
        return Verdict(UNS, "group.missing", f"semantic group '{group}' is not available",
                       {"group": group})
    return None


def rule_center(req, purpose, ctx) -> Verdict:
    target = req.centering_target
    missing = _target_check(target, ctx)
    if missing:
        return missing
    ctype = _GROUP_COMPONENT.get(target)
    context = {"target_group": target, "target_component": ctype,
               "target_state": ctx.component_summary(ctype) if ctype else {"state": "unmapped"},
               "membrane": ctx.state(ComponentType.MEMBRANE)}
    if purpose == P.DIFFUSION:
        return Verdict(REF, "center.removes_drift",
                       "centring removes centre-of-mass drift and re-wraps molecules", context)
    if ctx.resolved(ComponentType.MEMBRANE):
        if purpose in _MEMBRANE_FRAME_PURPOSES:
            if ctype == ComponentType.MEMBRANE:
                return Verdict(VWI, "center.membrane_on_membrane",
                               "centring the bilayer keeps the membrane frame, but the "
                               "reference choice is explicit", context)
            return Verdict(REF, "center.membrane_frame_foreign_target",
                           "centring on a non-membrane group re-wraps the bilayer and shifts "
                           "the membrane frame used by this quantity", context)
        if purpose == P.DISPLAY:
            return Verdict(VWI, "center.membrane_display",
                           "re-wrapping around the target can split the bilayer on display", context)
    if ctype is None or not ctx.resolved(ctype):
        return Verdict(VWI, "center.target_unresolved",
                       f"centring target '{target}' is not a RESOLVED component "
                       f"(state: {context['target_state']['state']})", context)
    if purpose in (P.DISPLAY, P.INTRAMOLECULAR_SHAPE):
        return Verdict(V, "center.whole_translation",
                       "translation + whole-molecule re-wrap; internal geometry unchanged", context)
    if purpose == P.INTER_COMPONENT_GEOMETRY:
        return Verdict(VWI, "center.inter_component_representation",
                       "-pbc mol re-wrapping chooses periodic images of the components — "
                       "a representation choice for inter-component geometry", context,
                       evidence=_pbc_pair_evidence(ctx))
    return Verdict(VWI, "center.purpose_unspecified",
                   f"centring for '{purpose}' needs an explicit reference choice", context)


def rule_fit(req, purpose, ctx) -> Verdict:
    group = req.fit_selection
    if req.fit_mode != "rot+trans":
        return Verdict(UNS, "fit.mode_unsupported",
                       f"fit mode '{req.fit_mode}' has no executor in build_view (only rot+trans)",
                       {"fit_mode": req.fit_mode})
    missing = _target_check(group, ctx)
    if missing:
        return missing
    base = group.split("_")[0]
    ctype = _GROUP_COMPONENT.get(base)
    context = {"fit_group": group, "fit_component": ctype,
               "component": ctx.component_summary(ctype) if ctype else {"state": "unmapped"},
               "membrane": ctx.state(ComponentType.MEMBRANE), "reference": ctx.fit_reference()}
    if purpose == P.DIFFUSION:
        return Verdict(REF, "fit.removes_motion",
                       "least-squares fitting removes the translation/rotation diffusion measures",
                       context)
    if ctx.resolved(ComponentType.MEMBRANE) and purpose in (P.DISPLAY, *_MEMBRANE_FRAME_PURPOSES):
        return Verdict(REF, "fit.membrane_frame",
                       "a rot+trans fit to a protein group rotates the membrane normal away from "
                       "z and changes leaflet / z / thickness / density meaning", context)
    if purpose in _MEMBRANE_FRAME_PURPOSES:
        return Verdict(VWI, "fit.rotating_spatial_frame",
                       "fitting rotates the spatial frame this quantity is binned in", context)
    if ctype is None or not ctx.resolved(ctype):
        return Verdict(VWI, "fit.group_unresolved",
                       f"fit group '{group}' is not a RESOLVED component "
                       f"(state: {context['component']['state']})", context)
    if context["reference"]["ambiguous"]:
        return Verdict(VWI, "fit.reference_ambiguous",
                       "the fit reference is the structure file while a .tpr also exists — "
                       "two plausible reference coordinate sets; choose explicitly", context)
    return Verdict(V, "fit.frame_alignment",
                   f"aligning to '{group}' defines the frame for '{purpose}'", context)


def rule_minimum_image(req, purpose, ctx) -> Verdict:
    return Verdict(V, "minimum_image.no_coordinate_change",
                   "the observable applies minimum-image geometry itself; coordinates are not "
                   "transformed", evidence=_pbc_pair_evidence(ctx))


def rule_reconstruction(req, purpose, ctx) -> Verdict:
    return Verdict(UNS, "reconstruction.no_executor",
                   f"component reconstruction '{req.reconstruction}' has no executor in "
                   f"build_view; scientifically it would be VALID_WITH_INTENT (a representation "
                   f"choice among cluster / nojump+center / minimum-image)",
                   {"requested": req.reconstruction,
                    "ligand": ctx.component_summary(ComponentType.LIGAND)},
                   evidence=_pbc_pair_evidence(ctx))


#: (operation, is-requested, parameters, rule) — evaluated in build_view order
_OPERATIONS: list[tuple[str, Callable, Callable, Callable]] = [
    ("make_whole", lambda r: r.requires_whole_molecules, lambda r: {}, rule_make_whole),
    ("nojump", lambda r: r.requires_nojump, lambda r: {}, rule_nojump),
    ("center", lambda r: bool(r.centering_target),
     lambda r: {"target": r.centering_target}, rule_center),
    ("fit", lambda r: bool(r.fit_selection),
     lambda r: {"group": r.fit_selection, "mode": r.fit_mode}, rule_fit),
    ("reconstruction", lambda r: bool(r.reconstruction),
     lambda r: {"method": r.reconstruction}, rule_reconstruction),
    ("minimum_image", lambda r: r.minimum_image_distances, lambda r: {}, rule_minimum_image),
]
_TRANSFORMATIONS = ("make_whole", "nojump", "center", "fit", "reconstruction")

#: detectors whose outcome (including "ran, found nothing") bears on each operation
_COVERAGE: dict[str, tuple[str, ...]] = {
    "nojump": ("timeline_integrity", "duplicate_timestamps", "non_monotonic_time",
               "time_gaps", "component_periodic_jump"),
    "center": ("box_behaviour", "component_periodic_jump", "partner_receptor_separation"),
    "fit": ("box_behaviour",),
    "minimum_image": ("partner_receptor_separation",),
    "reconstruction": ("component_periodic_jump", "partner_receptor_separation"),
}


def diagnostic_coverage(operation: str, ctx: PolicyContext) -> list[dict]:
    """Pure: the status of each relevant detector in the supplied report."""
    wanted = _COVERAGE.get(operation, ())
    if ctx.diagnostics is None:
        return [{"detector": d, "status": "not_supplied"} for d in wanted]
    runs = {r.detector_id: r for r in ctx.diagnostics.detectors}
    out = []
    for d in wanted:
        r = runs.get(d)
        if r is None:
            out.append({"detector": d, "status": "absent_from_report"})
            continue
        out.append({"detector": d, "version": r.version, "status": r.status,
                    "identity": r.identity, "n_diagnostics": r.n_diagnostics,
                    "reason": r.reason,
                    "sampling": {k: r.sampling.get(k) for k in
                                 ("strategy", "frames_examined", "total_frames")}
                    if r.sampling else {},
                    "finding": ("none" if r.status == "ran" and r.n_diagnostics == 0 else
                                "reported" if r.status == "ran" else "no_evidence")})
    return out


# ═══════════════════════════════════════════════════════════════════════════════
# Plan
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class PolicyPlan:
    purpose: str
    requested: TrajectoryRequirements
    decisions: list[PreprocessingDecision]
    intent: PolicyIntent
    resolved: Optional[TrajectoryRequirements] = None     # None => not executable
    blocked: list[str] = field(default_factory=list)
    policy_name: Optional[str] = None
    execution: str = "build_view"                           # or "external_protocol"
    diagnostics: dict = field(default_factory=dict)         # identity/status of the report used

    @property
    def executable(self) -> bool:
        return self.resolved is not None

    def to_dict(self) -> dict:
        return {"engine": POLICY_ENGINE, "purpose": self.purpose,
                "policy_name": self.policy_name, "execution": self.execution,
                "intent": self.intent.to_dict(), "requested": self.requested.to_dict(),
                "executable": self.executable, "blocked": self.blocked,
                "resolved": self.resolved.to_dict() if self.resolved else None,
                "diagnostics": self.diagnostics,
                "decisions": [d.to_dict() for d in self.decisions]}


def _report_summary(report: Optional[DiagnosticReport]) -> dict:
    from analysis.campaign.diagnostics.registry import diagnostics_status, report_identity
    if report is None:
        return {"supplied": False, "status": "unavailable"}
    fp = report.trajectory_fingerprint
    return {"supplied": True, "status": diagnostics_status(report),
            "report_identity": report_identity(report),
            "trajectory_digest": fp.digest if fp else None}


def plan_preprocessing(requirements: TrajectoryRequirements, *, purpose: Optional[str],
                       context: PolicyContext, intent: Optional[PolicyIntent] = None,
                       policy_name: Optional[str] = None) -> PolicyPlan:
    """Judge every requested operation.  Pure: runs nothing, writes nothing."""
    from analysis.campaign.results import definition_token
    intent = intent or PolicyIntent()
    purpose = purpose or "unspecified"
    decisions: list[PreprocessingDecision] = []
    for op, requested, params, rule in _OPERATIONS:
        if not requested(requirements):
            continue
        v = rule(requirements, purpose, context)
        applied = v.classification == V or (v.classification == VWI and intent.authorizes(op))
        d = PreprocessingDecision(
            operation=op, purpose=purpose, classification=v.classification, reason=v.reason,
            rule_id=v.rule_id, rule_version=RULE_VERSION, intent_source=intent.source,
            parameters=params(requirements), requested=True, applied=applied,
            system_context=v.context, diagnostic_evidence=v.evidence, notes=v.notes,
            diagnostic_coverage=diagnostic_coverage(op, context))
        if op == "nojump" and requirements.centering_target and v.classification != REF:
            d.notes.append("the following centring step re-wraps molecules (-pbc mol): the "
                           "resulting view is NOT globally unwrapped")
        d.identity = definition_token({
            "engine": POLICY_ENGINE, "rule": d.rule_id, "version": d.rule_version,
            "operation": op, "parameters": d.parameters, "purpose": purpose,
            "intent": intent.source, "context": d.system_context,
            "evidence": [(e["code"], e["detector"], e["detector_identity"])
                         for e in d.diagnostic_evidence],
            "coverage": [(c["detector"], c["status"], c.get("identity"))
                         for c in d.diagnostic_coverage]})
        decisions.append(d)
    blocked = [d.operation for d in decisions if not d.applied]
    plan = PolicyPlan(purpose=purpose, requested=requirements, decisions=decisions,
                      intent=intent, blocked=blocked, policy_name=policy_name,
                      diagnostics=_report_summary(context.diagnostics))
    if not blocked:
        plan.resolved = requirements
    return plan


# ═══════════════════════════════════════════════════════════════════════════════
# Named policies — presets of *requests*, always judged by the same rules
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class NamedPolicy:
    name: str
    description: str
    requirements: Callable[[], TrajectoryRequirements]
    execution: str = "build_view"
    recorded_chain: tuple[str, ...] = ()


def _soluble_complex_chain() -> tuple[str, ...]:
    from analysis.campaign.xanthone_short import PREPROCESSING_CHAIN
    return tuple(PREPROCESSING_CHAIN)


NAMED_POLICIES: dict[str, NamedPolicy] = {
    "raw": NamedPolicy("raw", "stored coordinates, no transformation",
                       lambda: TrajectoryRequirements(rationale="named policy: raw")),
    "whole-only": NamedPolicy("whole-only", "molecules made whole",
                              lambda: TrajectoryRequirements(requires_whole_molecules=True,
                                                             rationale="named policy: whole-only")),
    "soluble-fit": NamedPolicy("soluble-fit", "whole + rot+trans fit to the receptor backbone",
                               lambda: TrajectoryRequirements(
                                   requires_whole_molecules=True, fit_selection="Receptor_Backbone",
                                   rationale="named policy: soluble-fit")),
    "membrane-frame": NamedPolicy("membrane-frame", "whole + bilayer centred, no rotation",
                                  lambda: TrajectoryRequirements(
                                      requires_whole_molecules=True, centering_target="Membrane",
                                      rationale="named policy: membrane-frame")),
    # The xanthone_short / mechanistic protocol: executed OUTSIDE SimForge (mdfit.xtc is
    # consumed, never produced here).  Its operations are judged as center(Complex) +
    # fit(Complex_Backbone); `-pbc res -ur compact` has no build_view executor.
    "soluble-complex-fit": NamedPolicy(
        "soluble-complex-fit",
        "external protocol: trjconv -pbc res -ur compact -center, then -fit rot+trans",
        lambda: TrajectoryRequirements(
            requires_whole_molecules=True, centering_target="Complex",
            fit_selection="Complex_Backbone", rationale="named policy: soluble-complex-fit"),
        execution="external_protocol", recorded_chain=()),
}


def named_policy(name: str) -> NamedPolicy:
    pol = NAMED_POLICIES[name]
    if name == "soluble-complex-fit":
        return NamedPolicy(pol.name, pol.description, pol.requirements, pol.execution,
                           _soluble_complex_chain())
    return pol


def plan_named_policy(name: str, *, purpose: str, context: PolicyContext,
                      intent: Optional[PolicyIntent] = None) -> PolicyPlan:
    """A preset supplies the *request*; the rules still judge it.  Choosing a
    preset is profile intent unless another explicit source is given."""
    pol = named_policy(name)
    plan = plan_preprocessing(pol.requirements(), purpose=purpose, context=context,
                              intent=intent or PolicyIntent(IntentSource.PROFILE),
                              policy_name=name)
    plan.execution = pol.execution
    if pol.execution != "build_view":
        plan.resolved = None                  # never executed by build_view
        plan.blocked = plan.blocked or ["external_protocol"]
    return plan


def main(argv=None) -> int:
    """Developer dry run: ``python -m analysis.campaign.trajectory.policy -s STRUCT
    [-p TOPOL.tpr] [-n INDEX] [--diagnostics trajectory_diagnostics.json]
    --purpose P (--policy NAME | --whole --nojump --center G --fit G) [--intent SRC]``"""
    import argparse
    import json
    from analysis.campaign.structure.component_detector import detect_components
    ap = argparse.ArgumentParser(description="Evaluate preprocessing policy (no execution).")
    ap.add_argument("-s", "--structure", required=True)
    ap.add_argument("-p", "--topology")
    ap.add_argument("-n", "--index")
    ap.add_argument("--diagnostics")
    ap.add_argument("--purpose", required=True, choices=list(P.ALL))
    ap.add_argument("--policy", choices=sorted(NAMED_POLICIES))
    ap.add_argument("--whole", action="store_true")
    ap.add_argument("--nojump", action="store_true")
    ap.add_argument("--center")
    ap.add_argument("--fit")
    ap.add_argument("--intent", default=IntentSource.AUTO)
    ns = ap.parse_args(argv)
    diags = None
    if ns.diagnostics:
        diags = DiagnosticReport.from_dict(json.loads(Path(ns.diagnostics).read_text()))
    ctx = PolicyContext(components=detect_components(structure_path=ns.structure).components,
                        topology_path=ns.topology, structure_path=ns.structure,
                        index_path=ns.index, diagnostics=diags)
    intent = PolicyIntent(ns.intent)
    if ns.policy:
        plan = plan_named_policy(ns.policy, purpose=ns.purpose, context=ctx, intent=intent)
    else:
        req = TrajectoryRequirements(requires_whole_molecules=ns.whole, requires_nojump=ns.nojump,
                                     centering_target=ns.center, fit_selection=ns.fit)
        plan = plan_preprocessing(req, purpose=ns.purpose, context=ctx, intent=intent)
    print(json.dumps(plan.to_dict(), indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
