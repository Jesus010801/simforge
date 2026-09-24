"""Phase 5 — scientific preprocessing decision/policy layer."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from analysis.campaign.models import (
    ClassificationState as CS, ComponentType as CT, DecisionClass as DC, DetectorRun,
    Diagnostic, DiagnosticReport, IntentSource, MolecularComponent,
    ObservablePurpose as P, TrajectoryRequirements as TR,
)
from analysis.campaign.trajectory import policy as pol
from analysis.campaign.trajectory.policy import (
    PolicyContext, PolicyIntent, plan_named_policy, plan_preprocessing,
)

_GROUPS = ["System", "Receptor", "Receptor_Backbone", "Peptide_Backbone", "Complex",
           "Complex_Backbone", "Ligand", "Membrane"]


def comp(ctype, state=CS.RESOLVED, atom_ids=None):
    return MolecularComponent(component_type=ctype, label=ctype, classification_state=state,
                              atom_ids=atom_ids)


def ctx(tmp_path, components, *, structure=None, topology="md.tpr", diagnostics=None,
        groups=_GROUPS) -> PolicyContext:
    ndx = tmp_path / "semantic.ndx"
    ndx.write_text("".join(f"[ {g} ]\n{i + 1}\n" for i, g in enumerate(groups)))
    return PolicyContext(components=components, topology_path=str(tmp_path / topology),
                         structure_path=str(tmp_path / structure) if structure else None,
                         index_path=str(ndx), diagnostics=diagnostics)


def report(*diags: Diagnostic, timeline_ran=True) -> DiagnosticReport:
    runs = [DetectorRun("timeline_integrity", "1", 0, "ran")] if timeline_ran else []
    return DiagnosticReport(trajectory_path="x.xtc", detectors=runs, diagnostics=list(diags))


def diag(code, confidence="strong", severity="warn", convention="center_of_mass"):
    return Diagnostic(code=code, severity=severity, scope="x", message=code,
                      confidence=confidence, detector={"id": "det", "version": "1"},
                      provenance={"detector_identity": f"id-{code}-{confidence}",
                                  "centre_convention": convention})


SOLUBLE = [comp(CT.RECEPTOR), comp(CT.COMPLEX)]
MEMBRANE = [comp(CT.RECEPTOR), comp(CT.COMPLEX), comp(CT.MEMBRANE)]
LIGAND = [comp(CT.RECEPTOR), comp(CT.COMPLEX), comp(CT.LIGAND, atom_ids=[5, 6, 7])]


def decision(plan, op):
    (d,) = [d for d in plan.decisions if d.operation == op]
    return d


# ═══════════════════════════════════════════════════════════════════════════════
# Context dependence — same operation, different purpose
# ═══════════════════════════════════════════════════════════════════════════════

def test_same_nojump_same_diagnostic_different_purpose(tmp_path):
    rep = report(diag("component_periodic_wrap"))
    c = ctx(tmp_path, MEMBRANE, diagnostics=rep)
    req = TR(requires_whole_molecules=True, requires_nojump=True)
    diff = decision(plan_preprocessing(req, purpose=P.DIFFUSION, context=c), "nojump")
    memb = decision(plan_preprocessing(req, purpose=P.MEMBRANE_FRAME_PROPERTY, context=c), "nojump")
    assert diff.classification == DC.VALID and diff.applied
    assert memb.classification == DC.REFUSED and not memb.applied
    assert diff.diagnostic_evidence[0]["code"] == "component_periodic_wrap"


def test_display_and_analysis_get_different_decisions(tmp_path):
    c = ctx(tmp_path, MEMBRANE)
    req = TR(requires_whole_molecules=True, fit_selection="Receptor_Backbone")
    display = plan_preprocessing(req, purpose=P.DISPLAY, context=c)
    rmsd = plan_preprocessing(req, purpose=P.INTRAMOLECULAR_SHAPE, context=c)
    assert decision(display, "fit").classification == DC.REFUSED and not display.executable
    assert decision(rmsd, "fit").classification == DC.VALID and rmsd.executable


# ═══════════════════════════════════════════════════════════════════════════════
# Operation semantics
# ═══════════════════════════════════════════════════════════════════════════════

def test_make_whole(tmp_path):
    ok = plan_preprocessing(TR(requires_whole_molecules=True), purpose=P.INTRAMOLECULAR_SHAPE,
                            context=ctx(tmp_path, SOLUBLE))
    assert decision(ok, "make_whole").classification == DC.VALID
    assert decision(ok, "make_whole").system_context["connectivity_source"] == "tpr"
    no = plan_preprocessing(TR(requires_whole_molecules=True), purpose=P.INTRAMOLECULAR_SHAPE,
                            context=ctx(tmp_path, SOLUBLE, topology="md.gro"))
    assert decision(no, "make_whole").classification == DC.UNSUPPORTED and not no.executable


@pytest.mark.parametrize("purpose,expected", [
    (P.DIFFUSION, DC.VALID), (P.MEMBRANE_FRAME_PROPERTY, DC.REFUSED),
    (P.DENSITY_PROFILE, DC.REFUSED), (P.SOLVENT_OCCUPANCY, DC.REFUSED),
    (P.DISPLAY, DC.VALID_WITH_INTENT), (P.INTRAMOLECULAR_SHAPE, DC.VALID_WITH_INTENT),
])
def test_nojump_by_purpose(tmp_path, purpose, expected):
    plan = plan_preprocessing(TR(requires_whole_molecules=True, requires_nojump=True),
                              purpose=purpose, context=ctx(tmp_path, SOLUBLE, diagnostics=report()))
    assert decision(plan, "nojump").classification == expected


def test_nojump_then_center(tmp_path):
    c = ctx(tmp_path, SOLUBLE, diagnostics=report())
    req = TR(requires_whole_molecules=True, requires_nojump=True, centering_target="Receptor")
    d = decision(plan_preprocessing(req, purpose=P.DIFFUSION, context=c), "nojump")
    assert d.classification == DC.REFUSED and d.rule_id == "nojump.rewrapped_by_center"
    disp = decision(plan_preprocessing(req, purpose=P.DISPLAY, context=c,
                                       intent=PolicyIntent(IntentSource.USER_FLAG)), "nojump")
    assert any("NOT globally unwrapped" in n for n in disp.notes)


@pytest.mark.parametrize("components,purpose,target,expected", [
    (SOLUBLE, P.DISPLAY, "Receptor", DC.VALID),
    (SOLUBLE, P.INTRAMOLECULAR_SHAPE, "Receptor", DC.VALID),
    (SOLUBLE, P.DIFFUSION, "Receptor", DC.REFUSED),
    (SOLUBLE, P.INTER_COMPONENT_GEOMETRY, "Complex", DC.VALID_WITH_INTENT),
    (MEMBRANE, P.DISPLAY, "Receptor", DC.VALID_WITH_INTENT),
    (MEMBRANE, P.MEMBRANE_FRAME_PROPERTY, "Receptor", DC.REFUSED),
    (MEMBRANE, P.MEMBRANE_FRAME_PROPERTY, "Membrane", DC.VALID_WITH_INTENT),
    (SOLUBLE, P.DISPLAY, "Water", DC.UNSUPPORTED),            # group not in the index
])
def test_center_rules(tmp_path, components, purpose, target, expected):
    plan = plan_preprocessing(TR(requires_whole_molecules=True, centering_target=target),
                              purpose=purpose, context=ctx(tmp_path, components))
    assert decision(plan, "center").classification == expected


@pytest.mark.parametrize("components,purpose,structure,mode,expected", [
    (SOLUBLE, P.INTRAMOLECULAR_SHAPE, None, "rot+trans", DC.VALID),
    (SOLUBLE, P.INTRAMOLECULAR_SHAPE, "md.gro", "rot+trans", DC.VALID_WITH_INTENT),
    (SOLUBLE, P.DIFFUSION, None, "rot+trans", DC.REFUSED),
    (MEMBRANE, P.MEMBRANE_FRAME_PROPERTY, None, "rot+trans", DC.REFUSED),
    (MEMBRANE, P.DISPLAY, None, "rot+trans", DC.REFUSED),
    (MEMBRANE, P.DENSITY_PROFILE, None, "rot+trans", DC.REFUSED),
    (MEMBRANE, P.SOLVENT_OCCUPANCY, None, "rot+trans", DC.REFUSED),
    (MEMBRANE, P.DISPLAY, None, "transxy", DC.UNSUPPORTED),   # membrane-preserving: no executor
])
def test_fit_rules(tmp_path, components, purpose, structure, mode, expected):
    req = TR(requires_whole_molecules=True, fit_selection="Receptor_Backbone", fit_mode=mode)
    plan = plan_preprocessing(req, purpose=purpose,
                              context=ctx(tmp_path, components, structure=structure))
    d = decision(plan, "fit")
    assert d.classification == expected
    if structure:
        assert d.system_context["reference"]["reference_file_type"] == "gro"


def test_membrane_protein_fit_not_even_with_intent(tmp_path):
    plan = plan_preprocessing(TR(requires_whole_molecules=True, fit_selection="Receptor_Backbone"),
                              purpose=P.MEMBRANE_FRAME_PROPERTY, context=ctx(tmp_path, MEMBRANE),
                              intent=PolicyIntent(IntentSource.USER_FLAG))
    assert decision(plan, "fit").classification == DC.REFUSED and not plan.executable


def test_minimum_image_is_not_a_transformation(tmp_path):
    plan = plan_preprocessing(TR(requires_whole_molecules=True, minimum_image_distances=True),
                              purpose=P.INTER_COMPONENT_GEOMETRY, context=ctx(tmp_path, LIGAND))
    d = decision(plan, "minimum_image")
    assert d.classification == DC.VALID and plan.executable
    assert plan.resolved.view_kind() == "whole"                  # no extra operation added


# ═══════════════════════════════════════════════════════════════════════════════
# Intent
# ═══════════════════════════════════════════════════════════════════════════════

def test_valid_with_intent_needs_explicit_intent(tmp_path):
    c = ctx(tmp_path, SOLUBLE, diagnostics=report())
    req = TR(requires_whole_molecules=True, requires_nojump=True)
    auto = plan_preprocessing(req, purpose=P.DISPLAY, context=c)
    assert decision(auto, "nojump").classification == DC.VALID_WITH_INTENT
    assert not decision(auto, "nojump").applied and not auto.executable
    assert auto.blocked == ["nojump"]
    for src in (IntentSource.USER_FLAG, IntentSource.MANIFEST_RESOLUTION,
                IntentSource.EXPLICIT_API, IntentSource.PROFILE):
        ok = plan_preprocessing(req, purpose=P.DISPLAY, context=c, intent=PolicyIntent(src))
        assert ok.executable and decision(ok, "nojump").applied
    narrow = plan_preprocessing(req, purpose=P.DISPLAY, context=c,
                                intent=PolicyIntent(IntentSource.USER_FLAG, ("fit",)))
    assert not narrow.executable                                  # intent was for another op


def test_policy_never_invents_operations(tmp_path):
    rep = report(diag("component_periodic_wrap"), diag("ligand_receptor_periodic_separation"))
    plan = plan_preprocessing(TR(requires_whole_molecules=True), purpose=P.DIFFUSION,
                              context=ctx(tmp_path, LIGAND, diagnostics=rep))
    assert [d.operation for d in plan.decisions] == ["make_whole"]
    assert plan.resolved.cache_token() == TR(requires_whole_molecules=True).cache_token()


# ═══════════════════════════════════════════════════════════════════════════════
# Component state, not group existence
# ═══════════════════════════════════════════════════════════════════════════════

def test_ligand_group_without_resolved_ligand(tmp_path):
    comps = [comp(CT.RECEPTOR), comp(CT.COMPLEX), comp(CT.LIGAND, CS.AMBIGUOUS)]
    plan = plan_preprocessing(TR(requires_whole_molecules=True, centering_target="Ligand"),
                              purpose=P.DISPLAY, context=ctx(tmp_path, comps))
    d = decision(plan, "center")
    assert d.classification == DC.VALID_WITH_INTENT and d.rule_id == "center.target_unresolved"
    assert d.system_context["target_state"]["state"] == CS.AMBIGUOUS


def test_ambiguous_receptor_group_not_trusted(tmp_path):
    comps = [comp(CT.RECEPTOR, CS.AMBIGUOUS), comp(CT.COMPLEX)]
    plan = plan_preprocessing(TR(requires_whole_molecules=True, centering_target="Receptor"),
                              purpose=P.DISPLAY, context=ctx(tmp_path, comps))
    assert decision(plan, "center").classification == DC.VALID_WITH_INTENT


# ═══════════════════════════════════════════════════════════════════════════════
# Ligand / PBC evidence
# ═══════════════════════════════════════════════════════════════════════════════

def test_strong_ligand_pbc_never_picks_a_recipe(tmp_path):
    rep = report(diag("ligand_receptor_periodic_separation"),
                 diag("ligand_receptor_image_discrepancy", confidence="moderate"))
    c = ctx(tmp_path, LIGAND, diagnostics=rep)
    center = decision(plan_preprocessing(TR(requires_whole_molecules=True, centering_target="Complex"),
                                         purpose=P.INTER_COMPONENT_GEOMETRY, context=c), "center")
    assert center.classification == DC.VALID_WITH_INTENT and not center.applied
    assert {e["code"] for e in center.diagnostic_evidence} == {
        "ligand_receptor_periodic_separation", "ligand_receptor_image_discrepancy"}
    rec = decision(plan_preprocessing(TR(requires_whole_molecules=True, reconstruction="cluster"),
                                      purpose=P.INTER_COMPONENT_GEOMETRY, context=c,
                                      intent=PolicyIntent(IntentSource.USER_FLAG)),
                   "reconstruction")
    assert rec.classification == DC.UNSUPPORTED and not rec.applied   # no executor, no faking
    assert "VALID_WITH_INTENT" in rec.reason


def test_physical_separation_is_not_pbc_evidence(tmp_path):
    rep = report(diag("ligand_receptor_separation_change", confidence="weak", severity="info"))
    plan = plan_preprocessing(TR(requires_whole_molecules=True, minimum_image_distances=True),
                              purpose=P.INTER_COMPONENT_GEOMETRY,
                              context=ctx(tmp_path, LIGAND, diagnostics=rep))
    assert decision(plan, "minimum_image").diagnostic_evidence == []


def test_weak_evidence_never_stronger_than_strong(tmp_path):
    req = TR(requires_whole_molecules=True, centering_target="Complex")
    strong = plan_preprocessing(req, purpose=P.INTER_COMPONENT_GEOMETRY, context=ctx(
        tmp_path, LIGAND, diagnostics=report(diag("ligand_receptor_periodic_separation"))))
    weak = plan_preprocessing(req, purpose=P.INTER_COMPONENT_GEOMETRY, context=ctx(
        tmp_path, LIGAND, diagnostics=report(diag("ligand_receptor_periodic_separation",
                                                  confidence="weak", convention="center_of_geometry"))))
    order = {DC.REFUSED: 0, DC.UNSUPPORTED: 0, DC.VALID_WITH_INTENT: 1, DC.VALID: 2}
    assert order[decision(weak, "center").classification] <= order[decision(strong, "center").classification]
    assert decision(weak, "center").diagnostic_evidence[0]["centre_convention"] == "center_of_geometry"
    # evidence can never upgrade: a strong wrap does not make DISPLAY nojump VALID
    disp = plan_preprocessing(TR(requires_whole_molecules=True, requires_nojump=True),
                              purpose=P.DISPLAY,
                              context=ctx(tmp_path, SOLUBLE, diagnostics=report(diag("component_periodic_wrap"))))
    assert decision(disp, "nojump").classification == DC.VALID_WITH_INTENT


# ═══════════════════════════════════════════════════════════════════════════════
# Timeline integrity
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize("rep,expected,rule", [
    (report(diag("non_monotonic_time", severity="review_required")), DC.REFUSED, "nojump.timeline_invalid"),
    (report(diag("time_index_invalid", severity="error")), DC.REFUSED, "nojump.timeline_invalid"),
    (report(diag("duplicate_timestamps")), DC.VALID_WITH_INTENT, "nojump.timeline_review"),
    (None, DC.VALID_WITH_INTENT, "nojump.timeline_unverified"),
])
def test_timeline_gates_nojump(tmp_path, rep, expected, rule):
    plan = plan_preprocessing(TR(requires_whole_molecules=True, requires_nojump=True),
                              purpose=P.DIFFUSION, context=ctx(tmp_path, SOLUBLE, diagnostics=rep))
    d = decision(plan, "nojump")
    assert (d.classification, d.rule_id) == (expected, rule)


# ═══════════════════════════════════════════════════════════════════════════════
# Identity / provenance
# ═══════════════════════════════════════════════════════════════════════════════

def _ident(tmp_path, *, purpose=P.INTER_COMPONENT_GEOMETRY, intent=IntentSource.AUTO,
           lig_atoms=(5, 6, 7), diags=(), groups=_GROUPS, sub="a"):
    d = tmp_path / sub
    d.mkdir(exist_ok=True)
    comps = [comp(CT.RECEPTOR), comp(CT.COMPLEX), comp(CT.LIGAND, atom_ids=list(lig_atoms))]
    plan = plan_preprocessing(TR(requires_whole_molecules=True, reconstruction="cluster"),
                              purpose=purpose, intent=PolicyIntent(intent),
                              context=ctx(d, comps, diagnostics=report(*diags), groups=groups))
    return decision(plan, "reconstruction").identity


def test_decision_identity(tmp_path, monkeypatch):
    base = _ident(tmp_path)
    assert _ident(tmp_path, sub="elsewhere") == base                  # paths never enter
    assert _ident(tmp_path, groups=_GROUPS + ["Unused"], sub="u") == base
    assert _ident(tmp_path, purpose=P.DISPLAY) != base
    assert _ident(tmp_path, intent=IntentSource.USER_FLAG) != base
    assert _ident(tmp_path, lig_atoms=(5, 6, 8)) != base
    assert _ident(tmp_path, diags=[diag("ligand_receptor_periodic_separation")]) != base
    monkeypatch.setattr(pol, "RULE_VERSION", pol.RULE_VERSION + "-test")
    assert _ident(tmp_path) != base


# ═══════════════════════════════════════════════════════════════════════════════
# No automatic remediation
# ═══════════════════════════════════════════════════════════════════════════════

def test_policy_evaluation_runs_nothing(tmp_path, monkeypatch):
    import analysis.campaign.gmx as gmxmod
    import analysis.campaign.trajectory.preprocessor as pre

    def boom(*a, **k):
        raise AssertionError("policy evaluation must not execute anything")
    monkeypatch.setattr(gmxmod, "run_gmx", boom)
    monkeypatch.setattr(pre, "build_view", boom)
    rep = report(diag("component_periodic_wrap"), diag("ligand_receptor_periodic_separation"))
    c = ctx(tmp_path, LIGAND, diagnostics=rep)
    before = sorted(p.name for p in tmp_path.rglob("*"))
    for purpose in P.ALL:
        plan_preprocessing(TR(requires_whole_molecules=True, requires_nojump=True,
                              centering_target="Complex", fit_selection="Complex_Backbone"),
                           purpose=purpose, context=c)
    assert sorted(p.name for p in tmp_path.rglob("*")) == before


# ═══════════════════════════════════════════════════════════════════════════════
# Named policies and the existing protocol anchor
# ═══════════════════════════════════════════════════════════════════════════════

def test_named_policy_does_not_bypass_rules(tmp_path):
    plan = plan_named_policy("soluble-fit", purpose=P.DISPLAY, context=ctx(tmp_path, MEMBRANE))
    assert plan.intent.source == IntentSource.PROFILE
    assert decision(plan, "fit").classification == DC.REFUSED and not plan.executable
    mf = plan_named_policy("membrane-frame", purpose=P.MEMBRANE_FRAME_PROPERTY,
                           context=ctx(tmp_path, MEMBRANE))
    assert decision(mf, "center").classification == DC.VALID_WITH_INTENT
    assert mf.executable                                          # profile intent authorises it


def test_soluble_complex_fit_records_the_existing_protocol(tmp_path):
    from analysis.campaign.xanthone_short import PREPROCESSING_CHAIN
    assert PREPROCESSING_CHAIN == (                               # byte-for-byte anchor
        "md.xtc",
        "gmx trjconv -pbc res -ur compact -center -> mdcenter.xtc",
        "gmx trjconv -fit rot+trans -> mdfit.xtc",
    )
    named = pol.named_policy("soluble-complex-fit")
    assert named.recorded_chain == PREPROCESSING_CHAIN and named.execution == "external_protocol"
    plan = plan_named_policy("soluble-complex-fit", purpose=P.INTER_COMPONENT_GEOMETRY,
                             context=ctx(tmp_path, LIGAND))
    assert plan.execution == "external_protocol" and not plan.executable
    assert decision(plan, "fit").classification == DC.VALID
    assert decision(plan, "center").classification == DC.VALID_WITH_INTENT
    assert decision(plan, "center").applied                       # profile = explicit intent


def test_mechanistic_chain_unchanged():
    src = Path("analysis/campaign/mechanistic_exec.py").read_text()
    assert '"trjconv -pbc res -ur compact -center -> mdcenter.xtc",' in src
    assert '"trjconv -fit rot+trans -> mdfit.xtc",' in src


# ═══════════════════════════════════════════════════════════════════════════════
# build_view integration (fake gmx from the Phase 3 suite)
# ═══════════════════════════════════════════════════════════════════════════════

from tests.analysis.campaign.test_view_hardening import build, env  # noqa: E402,F401


def test_policy_planned_view_records_decisions(env):
    from analysis.campaign.models import SemanticIndex, SemanticIndexGroup
    idx = SemanticIndex(path=str(env["ndx"]), groups=[SemanticIndexGroup(n, 0) for n in
                                                      ("System", "Protein", "Unused")])
    c = PolicyContext(components=SOLUBLE, topology_path=str(env["tpr"]), structure_path=None,
                      index_path=str(env["ndx"]))
    req = TR(requires_whole_molecules=True)
    plan = plan_preprocessing(req, purpose=P.INTRAMOLECULAR_SHAPE, context=c)
    v = build(env, req, decisions=plan.decisions, semantic_index=idx)
    assert v.policy_planned and v.decisions[0]["classification"] == DC.VALID
    m = json.loads(Path(v.manifest_path).read_text())
    assert m["policy"]["planned"] is True and m["policy"]["decisions"][0]["rule_id"] == \
        "make_whole.tpr_connectivity"
    legacy = build(env, req, work="legacy")
    assert legacy.policy_planned is False and legacy.decisions == []


def test_build_view_refuses_unadmitted_or_unsupported(env):
    from analysis.campaign.models import PreprocessingDecision
    blocked = PreprocessingDecision(operation="fit", purpose="display", classification=DC.REFUSED,
                                    reason="r", rule_id="x", rule_version="1", applied=False)
    v = build(env, decisions=[blocked])
    assert not v.safe and v.path is None and any(w.code == "policy_blocked" for w in v.warnings)
    u = build(env, TR(requires_whole_molecules=True, fit_selection="Protein", fit_mode="transxy"))
    assert not u.safe and any(w.code == "unsupported_operation" for w in u.warnings)
    r = build(env, TR(reconstruction="cluster"))
    assert not r.safe and any(w.code == "unsupported_operation" for w in r.warnings)


# ═══════════════════════════════════════════════════════════════════════════════
# Real bundled membrane system (read-only)
# ═══════════════════════════════════════════════════════════════════════════════

from tests.analysis.campaign.conftest import (  # noqa: E402
    REAL_GRO, REAL_TPR, requires_gmx, requires_real_traj,
)


@requires_gmx
@requires_real_traj
def test_real_membrane_decisions(tmp_path):
    from analysis.campaign.structure.component_detector import detect_components
    from analysis.campaign.structure.index_builder import build_semantic_index
    det = detect_components(structure_path=REAL_GRO)
    idx = build_semantic_index(structure_path=REAL_GRO, components=det.components,
                               out_ndx=tmp_path / "s.ndx")
    c = PolicyContext(components=det.components, topology_path=str(REAL_TPR),
                      structure_path=str(REAL_GRO), index_path=idx.path)
    fit = TR(requires_whole_molecules=True, fit_selection="Receptor_Backbone")
    assert decision(plan_preprocessing(fit, purpose=P.MEMBRANE_FRAME_PROPERTY, context=c),
                    "fit").classification == DC.REFUSED
    assert decision(plan_preprocessing(fit, purpose=P.DISPLAY, context=c),
                    "fit").classification == DC.REFUSED
    assert decision(plan_preprocessing(fit, purpose=P.INTRAMOLECULAR_SHAPE, context=c),
                    "fit").classification == DC.VALID_WITH_INTENT          # .gro vs .tpr reference
    nj = TR(requires_whole_molecules=True, requires_nojump=True)
    assert decision(plan_preprocessing(nj, purpose=P.DIFFUSION, context=c),
                    "nojump").classification == DC.VALID_WITH_INTENT       # timeline not diagnosed
    whole = plan_preprocessing(TR(requires_whole_molecules=True), purpose=P.INTRAMOLECULAR_SHAPE,
                               context=c)
    assert whole.executable
