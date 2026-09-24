"""Phase 5.1 — run_analyze feeds the read-only Phase 4 diagnostics to the Phase 5
policy (once per system), without new detectors, rules or remediation."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from analysis.campaign.gmx import run_gmx
from analysis.campaign.models import (
    AnalysisResult, AnalysisStatus, ComponentType, DecisionClass as DC,
    ObservablePurpose as P, TrajectoryRequirements as TR,
)
from analysis.campaign.observables import registry
from analysis.campaign.observables.base import ObservableSpec
from analysis.campaign.orchestration import study_analyzer as sa
from analysis.campaign.orchestration.study_analyzer import clear_view_cache, run_analyze
from tests.analysis.campaign.conftest import (
    REAL_GRO, REAL_TPR, REAL_XTC, requires_gmx, requires_real_traj,
)

pytestmark = requires_gmx


def probe(pid, purpose, req, needs=(ComponentType.RECEPTOR,)):
    class _Probe(ObservableSpec):
        id = pid
        required_components = tuple(needs)

        def trajectory_requirements(self, params):
            return req

        def execute(self, ctx):
            v = ctx.trajectory_view
            ok = v.safe and v.path is not None
            return AnalysisResult(self.id, ctx.system.system_id,
                                  AnalysisStatus.SUCCESS if ok else AnalysisStatus.FAILED,
                                  message="; ".join(w.message for w in v.warnings))
    _Probe.purpose = purpose
    registry.register(_Probe())
    return pid


def view_of(result):
    return json.loads(Path(result.provenance_path).read_text())["trajectory_view"]


def decision(view, op):
    (d,) = [d for d in view["decisions"] if d["operation"] == op]
    return d


def coverage(d, det):
    (c,) = [c for c in d["diagnostic_coverage"] if c["detector"] == det]
    return c


def sha(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ═══════════════════════════════════════════════════════════════════════════════
# Study fixtures
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.fixture
def real_study(tmp_path):
    d = tmp_path / "campaign" / "sys" / "rep1"
    d.mkdir(parents=True)
    (d / "production_run.xtc").symlink_to(REAL_XTC)
    (d / "system_topol.tpr").symlink_to(REAL_TPR)
    (d / "reference.gro").symlink_to(REAL_GRO)
    return tmp_path / "campaign"


AA = ["ALA", "GLY", "SER", "LEU", "VAL", "THR", "LYS", "GLU", "ASP", "ILE"]


def synthetic_study(root: Path, times, lig_x) -> Path:
    """Protein (30 CA) + LIG (10 heavy atoms) + 1 water, 10 nm box; no .tpr."""
    d = root / "campaign" / "sys" / "rep1"
    d.mkdir(parents=True)
    frames = []
    for f, (t, lx) in enumerate(zip(times, lig_x)):
        atoms = [(i + 1, AA[i % 10], "CA", 8.9 + 0.2 * (i % 5), 4.0 + 0.1 * i, 5.0) for i in range(30)]
        atoms += [(31, "LIG", f"C{k + 1}", lx, 5.0 + 0.05 * k, 5.2) for k in range(10)]
        atoms += [(32, "SOL", "OW", 3.0, 3.0, 3.0)]
        lines = [f"synthetic t= {t:.5f} step= {f}", str(len(atoms))]
        lines += [f"{rs:>5}{rn:<5}{an:>5}{i:>5}{x:8.3f}{y:8.3f}{z:8.3f}"
                  for i, (rs, rn, an, x, y, z) in enumerate(atoms, 1)]
        lines.append("  10.00000  10.00000  10.00000")
        frames.append("\n".join(lines))
    gro_all = root / "frames.gro"
    gro_all.write_text("\n".join(frames) + "\n")
    (d / "reference.gro").write_text(frames[0] + "\n")
    (d / "topol.top").write_text("[ molecules ]\nProtein 1\nLIG 1\nSOL 1\n")
    res = run_gmx(["trjconv", "-f", str(gro_all), "-s", str(gro_all),
                   "-o", str(d / "production_run.xtc")], stdin="0\n")
    assert res.ok, res.stderr[-300:]
    return root / "campaign"


STATIC_LIG = [9.7] * 5
TIMES_OK = [0, 20, 40, 60, 80]


# ═══════════════════════════════════════════════════════════════════════════════
# The integration gap — before (diagnostics off) and after (auto)
# ═══════════════════════════════════════════════════════════════════════════════

@requires_real_traj
def test_diffusion_nojump_gets_diagnosed_timeline(real_study):
    pid = probe("probe-diffusion", P.DIFFUSION,
                TR(requires_whole_molecules=True, requires_nojump=True))
    # previous behaviour (still what "off" gives): nothing diagnosed -> needs intent
    off = run_analyze(real_study, [pid], inspect_trajectories=True, diagnostics="off").results[0]
    d_off = decision(view_of(off), "nojump")
    assert (d_off["classification"], d_off["rule_id"]) == (DC.VALID_WITH_INTENT,
                                                           "nojump.timeline_unverified")
    assert off.status == AnalysisStatus.FAILED
    assert view_of(off)["policy_ref"]["diagnostics"]["execution"] == "off"
    clear_view_cache()

    r = run_analyze(real_study, [pid], inspect_trajectories=True).results[0]
    v = view_of(r)
    d = decision(v, "nojump")
    assert (d["classification"], d["rule_id"]) == (DC.VALID, "nojump.diffusion")
    assert r.status == AnalysisStatus.SUCCESS
    cov = coverage(d, "timeline_integrity")
    assert cov["status"] == "ran" and cov["finding"] == "none"
    assert cov["sampling"]["frames_examined"] == 51
    ref = v["policy_ref"]["diagnostics"]
    assert ref["execution"] == "auto" and ref["status"] == "complete"
    report = json.loads(Path(ref["report_path"]).read_text())
    from analysis.campaign.diagnostics import report_identity
    from analysis.campaign.models import DiagnosticReport
    assert report_identity(DiagnosticReport.from_dict(report)) == ref["report_identity"]
    manifest = json.loads(Path(v["manifest_path"]).read_text())
    assert manifest["policy"]["context"]["diagnostics"]["report_identity"] == ref["report_identity"]


def test_non_monotonic_timeline_refuses_nojump(tmp_path):
    root = synthetic_study(tmp_path, [0, 20, 40, 10, 30], STATIC_LIG)
    pid = probe("probe-nojump", P.DIFFUSION, TR(requires_nojump=True))
    r = run_analyze(root, [pid], inspect_trajectories=False).results[0]
    d = decision(view_of(r), "nojump")
    assert (d["classification"], d["rule_id"]) == (DC.REFUSED, "nojump.timeline_invalid")
    assert d["diagnostic_evidence"][0]["code"] == "non_monotonic_time"
    assert r.status == AnalysisStatus.FAILED
    assert not list(root.rglob("nojump.xtc"))                      # nothing built


def test_good_synthetic_timeline_allows_nojump(tmp_path):
    root = synthetic_study(tmp_path, TIMES_OK, STATIC_LIG)
    pid = probe("probe-nojump", P.DIFFUSION, TR(requires_nojump=True))
    r = run_analyze(root, [pid], inspect_trajectories=False).results[0]
    d = decision(view_of(r), "nojump")
    # admitted by policy (execution itself is covered on the real .tpr system; this
    # synthetic study only has a .top, which trjconv cannot use as -s)
    assert d["classification"] == DC.VALID and d["applied"] is True


# ═══════════════════════════════════════════════════════════════════════════════
# "Ran, found nothing" vs "not applicable" vs "failed"
# ═══════════════════════════════════════════════════════════════════════════════

def test_ran_no_finding_differs_from_not_applicable_and_failed(tmp_path, monkeypatch, request):
    pid = probe("probe-mi", P.INTER_COMPONENT_GEOMETRY, TR(minimum_image_distances=True))
    root = synthetic_study(tmp_path / "lig", TIMES_OK, STATIC_LIG)
    ran = decision(view_of(run_analyze(root, [pid], inspect_trajectories=False).results[0]),
                   "minimum_image")
    c = coverage(ran, "partner_receptor_separation")
    assert (c["status"], c["finding"]) == ("ran", "none")

    # not applicable: bundled membrane system has no ligand/peptide/cofactor
    if REAL_XTC.is_file():
        clear_view_cache()
        real = request.getfixturevalue("real_study")
        na = decision(view_of(run_analyze(real, [pid], inspect_trajectories=True).results[0]),
                      "minimum_image")
        c = coverage(na, "partner_receptor_separation")
        assert (c["status"], c["finding"]) == ("not_applicable", "no_evidence")

    # failed: the motion probe breaks -> detector FAILED, report PARTIAL
    from analysis.campaign.diagnostics import context as dctx
    def broken(self, roles):
        raise dctx.ProbeError("simulated probe failure")
    monkeypatch.setattr(dctx.DiagnosticContext, "_probe_motion", broken)
    clear_view_cache()
    root2 = synthetic_study(tmp_path / "fail", TIMES_OK, STATIC_LIG)
    v = view_of(run_analyze(root2, [pid], inspect_trajectories=False).results[0])
    c = coverage(decision(v, "minimum_image"), "partner_receptor_separation")
    assert (c["status"], c["finding"]) == ("failed", "no_evidence")
    assert v["policy_ref"]["diagnostics"]["status"] == "partial"


def test_unavailable_diagnostics_stay_conservative(tmp_path):
    root = synthetic_study(tmp_path, TIMES_OK, STATIC_LIG)
    pid = probe("probe-nojump", P.DIFFUSION, TR(requires_nojump=True))
    r = run_analyze(root, [pid], inspect_trajectories=False, diagnostics="off").results[0]
    d = decision(view_of(r), "nojump")
    assert d["classification"] == DC.VALID_WITH_INTENT
    assert {c["status"] for c in d["diagnostic_coverage"]} == {"not_supplied"}


# ═══════════════════════════════════════════════════════════════════════════════
# PBC evidence reaches policy — and still chooses no recipe
# ═══════════════════════════════════════════════════════════════════════════════

def test_ligand_pbc_evidence_reaches_policy_without_recipe(tmp_path):
    root = synthetic_study(tmp_path, TIMES_OK, [9.70, 9.80, 0.10, 0.15, 0.20])
    center = probe("probe-center", P.INTER_COMPONENT_GEOMETRY, TR(centering_target="Complex"))
    mi = probe("probe-mi", P.INTER_COMPONENT_GEOMETRY, TR(minimum_image_distances=True))
    res = {r.analysis_id: r for r in run_analyze(root, [center, mi],
                                                 inspect_trajectories=False).results}
    dc = decision(view_of(res[center]), "center")
    codes = {e["code"] for e in dc["diagnostic_evidence"]}
    assert "ligand_receptor_periodic_separation" in codes
    assert dc["classification"] == DC.VALID_WITH_INTENT and not dc["applied"]
    assert [d["operation"] for d in view_of(res[center])["decisions"]] == ["center"]   # no recipe added
    assert res[center].status == AnalysisStatus.FAILED
    dm = decision(view_of(res[mi]), "minimum_image")
    assert dm["classification"] == DC.VALID and res[mi].status == AnalysisStatus.SUCCESS
    assert "ligand_receptor_periodic_separation" in {e["code"] for e in dm["diagnostic_evidence"]}


def test_physical_separation_not_pbc_evidence(tmp_path):
    root = synthetic_study(tmp_path, TIMES_OK, [9.70, 9.70, 5.80, 5.70, 5.70])
    mi = probe("probe-mi", P.INTER_COMPONENT_GEOMETRY, TR(minimum_image_distances=True))
    v = view_of(run_analyze(root, [mi], inspect_trajectories=False).results[0])
    report = json.loads(Path(v["policy_ref"]["diagnostics"]["report_path"]).read_text())
    assert "ligand_receptor_separation_change" in {d["code"] for d in report["diagnostics"]}
    dm = decision(v, "minimum_image")
    assert dm["diagnostic_evidence"] == []                     # not treated as PBC evidence
    c = coverage(dm, "partner_receptor_separation")
    assert (c["status"], c["finding"]) == ("ran", "reported")


@requires_real_traj
def test_membrane_fit_refusal_is_semantic_not_diagnostic(real_study):
    pid = probe("probe-membrane-fit", P.MEMBRANE_FRAME_PROPERTY,
                TR(requires_whole_molecules=True, fit_selection="Receptor_Backbone"))
    v = view_of(run_analyze(real_study, [pid], inspect_trajectories=True).results[0])
    d = decision(v, "fit")
    assert (d["classification"], d["rule_id"]) == (DC.REFUSED, "fit.membrane_frame")
    assert d["diagnostic_evidence"] == []                      # box_shape_drift did not decide it
    c = coverage(d, "box_behaviour")
    assert c["status"] == "ran" and c["n_diagnostics"] == 1    # the INFO drift is visible


# ═══════════════════════════════════════════════════════════════════════════════
# Once per system; cache reuse; supplied reports
# ═══════════════════════════════════════════════════════════════════════════════

def test_diagnostics_once_per_system_and_cached_across_runs(tmp_path, monkeypatch):
    import analysis.campaign.diagnostics as dmod
    from analysis.campaign.diagnostics import context as dctx
    root = synthetic_study(tmp_path, TIMES_OK, STATIC_LIG)
    a = probe("probe-a", P.INTER_COMPONENT_GEOMETRY, TR(minimum_image_distances=True))
    b = probe("probe-b", P.DIFFUSION, TR(requires_nojump=True))
    calls, gmx_calls = [], []
    real_diag, real_run = dmod.diagnose_system, dctx.run_gmx
    monkeypatch.setattr(dmod, "diagnose_system", lambda rec, **kw: calls.append(1) or real_diag(rec, **kw))
    monkeypatch.setattr(dctx, "run_gmx", lambda args, **kw: gmx_calls.append(args[0]) or real_run(args, **kw))
    res = run_analyze(root, [a, b], inspect_trajectories=False)
    assert len(calls) == 1                                       # one report, two observables
    ids = {view_of(r)["policy_ref"]["diagnostics"]["report_identity"] for r in res.results}
    assert len(ids) == 1
    assert gmx_calls.count("trajectory") == 1
    clear_view_cache()
    res2 = run_analyze(root, [a, b], inspect_trajectories=False)
    assert gmx_calls.count("trajectory") == 1                    # probe evidence reused from disk
    ref2 = view_of(res2.results[0])["policy_ref"]["diagnostics"]
    assert ref2["probes_from_cache"] is True and ref2["report_identity"] in ids


def test_supplied_report_validated(tmp_path):
    from analysis.campaign.diagnostics import diagnose_system
    from analysis.campaign.manifest import build_manifest
    root = synthetic_study(tmp_path / "a", TIMES_OK, STATIC_LIG)
    other = synthetic_study(tmp_path / "b", [0, 10, 20, 30, 40], STATIC_LIG)
    pid = probe("probe-mi", P.INTER_COMPONENT_GEOMETRY, TR(minimum_image_distances=True))

    # a report for the same trajectory (built the way run_analyze would) is accepted
    first = run_analyze(root, [pid], inspect_trajectories=False)
    ref = view_of(first.results[0])["policy_ref"]["diagnostics"]
    from analysis.campaign.models import DiagnosticReport
    good = DiagnosticReport.from_dict(json.loads(Path(ref["report_path"]).read_text()))
    clear_view_cache()
    ok = run_analyze(root, [pid], inspect_trajectories=False, diagnostics=good)
    assert view_of(ok.results[0])["policy_ref"]["diagnostics"]["execution"] == "supplied"

    # a report for another trajectory is rejected (and a fresh one generated)
    clear_view_cache()
    bad = DiagnosticReport.from_dict(json.loads(Path(view_of(run_analyze(
        other, [pid], inspect_trajectories=False).results[0])["policy_ref"]["diagnostics"]
        ["report_path"]).read_text()))
    clear_view_cache()
    rej = run_analyze(root, [pid], inspect_trajectories=False, diagnostics=bad)
    assert any(w.code == "diagnostics_report_rejected" and "fingerprint" in w.message
               for w in rej.warnings)
    assert view_of(rej.results[0])["policy_ref"]["diagnostics"]["execution"] == "auto"


# ═══════════════════════════════════════════════════════════════════════════════
# Read-only, no remediation, policy purity
# ═══════════════════════════════════════════════════════════════════════════════

def test_integration_is_read_only_and_builds_nothing_when_blocked(tmp_path, monkeypatch):
    root = synthetic_study(tmp_path, TIMES_OK, [9.70, 9.80, 0.10, 0.15, 0.20])
    src = list(root.rglob("production_run.xtc")) + list(root.rglob("reference.gro"))
    before = {p: sha(p) for p in src}
    def no_build(**kw):
        raise AssertionError("blocked plans must not reach build_view")
    monkeypatch.setattr(sa, "build_view", no_build)
    pid = probe("probe-center", P.INTER_COMPONENT_GEOMETRY, TR(centering_target="Complex"))
    r = run_analyze(root, [pid], inspect_trajectories=False).results[0]
    assert r.status == AnalysisStatus.FAILED
    assert {p: sha(p) for p in src} == before
    assert not list(root.rglob("*.xtc.building*")) and not list(root.rglob("centered.xtc"))


def test_plan_preprocessing_stays_pure(tmp_path, monkeypatch):
    import analysis.campaign.gmx as gmxmod
    from analysis.campaign.models import DetectorRun, DiagnosticReport
    from analysis.campaign.trajectory.policy import PolicyContext, plan_preprocessing
    monkeypatch.setattr(gmxmod, "run_gmx", lambda *a, **k: (_ for _ in ()).throw(AssertionError))
    rep = DiagnosticReport(trajectory_path="x", detectors=[
        DetectorRun("timeline_integrity", "1", 0, "ran", identity="t"),
        DetectorRun("partner_receptor_separation", "1", 2, "failed", reason="probe")])
    plan = plan_preprocessing(TR(requires_nojump=True, minimum_image_distances=True),
                              purpose=P.DIFFUSION, context=PolicyContext(diagnostics=rep))
    assert plan.diagnostics["supplied"] and plan.diagnostics["status"] == "partial"
    mi = [d for d in plan.decisions if d.operation == "minimum_image"][0]
    assert mi.diagnostic_coverage[0]["status"] == "failed"
