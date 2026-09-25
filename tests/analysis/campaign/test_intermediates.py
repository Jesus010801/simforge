"""Phase 13 — shared scientific intermediates and the dependency DAG."""
from __future__ import annotations

import json
import threading
from pathlib import Path

import numpy as np
import pytest

from analysis.campaign.intermediates import (
    DependencyCycleError, IntermediateConflict, IntermediateRequest as R, IntermediateScope,
    IntermediateSpec, IntermediateStatus as S, IntermediateStore, plan_dependencies, resolve_plan,
)
from analysis.campaign.intermediates import registry as ireg
from analysis.campaign.intermediates.executor import consumer_dependencies, identities
from analysis.campaign.models import (
    AnalysisResult, AnalysisStatus, Axis, AxisKind, FrameAlignmentMode, ResultArray,
    SourceFingerprint, StorageRef, SystemRecord, TrajectoryRequirements,
)
from tests.analysis.campaign.conftest import requires_gmx, requires_real_traj

# ═══════════════════════════════════════════════════════════════════════════════
# Test-only specs (no scientific meaning; architecture only)
# ═══════════════════════════════════════════════════════════════════════════════

CALLS: dict[str, int] = {}


class Dummy(IntermediateSpec):
    scope = IntermediateScope.STRUCTURE_STATIC

    def __init__(self, id, deps=(), version=1, fail=False, noise=False, structured=False):
        self.id, self._deps, self.version = id, list(deps), version
        self.fail, self.noise, self.structured = fail, noise, structured

    def parameters_schema(self):
        return {"x": "any"}

    def dependencies(self, params):
        return [d(params) if callable(d) else d for d in self._deps]

    def definition_evidence(self, ctx):
        return {"x": ctx.parameters.get("x")}

    def compute(self, ctx):
        CALLS[self.id] = CALLS.get(self.id, 0) + 1
        if self.fail:
            raise RuntimeError("deliberate failure")
        seed = [len(ctx.dependencies), ctx.parameters.get("x", 0) or 0]
        if self.noise:
            seed.append(np.random.random())
        f = ctx.work_dir / "values.npy"
        np.save(f, np.asarray(seed, dtype=float))
        arrays = [ResultArray("v", "dummy", None, axes=[Axis("i", AxisKind.INDEX)],
                              storage=StorageRef(str(f), "npy"))]
        artifacts = []
        if self.structured:
            m = ctx.work_dir / "membership.json"
            m.write_text(json.dumps({"frame0": {"atom1": "upper"}}, sort_keys=True))
            artifacts.append(StorageRef(str(m), "npz" if False else "npy", member=None,
                                        attrs={"schema": "frame->atom->region", "media": "json"}))
        return arrays, artifacts, {"note": "test"}


def rec() -> SystemRecord:
    r = SystemRecord("s", "c", "r")
    r.source_fingerprints = {"structure": SourceFingerprint("/s.gro", 1, 1, "strong", "gro-1")}
    return r


def run(plan_consumers, store_root, system=None):
    plan = plan_dependencies(plan_consumers)
    return plan, resolve_plan(plan, system=system or rec(), semantic_index=None,
                              topology_path="", structure_path=None,
                              store=IntermediateStore(store_root))


@pytest.fixture(autouse=True)
def _reset_calls():
    CALLS.clear()


# ═══════════════════════════════════════════════════════════════════════════════
# Planner
# ═══════════════════════════════════════════════════════════════════════════════

def test_fan_out_nested_and_deterministic_order():
    a, b, x = Dummy("t-a"), Dummy("t-b", deps=[R.make("t-a", {"x": 1})]), Dummy("t-x")
    with ireg.temporarily_registered(a, b, x):
        consumers = [("obs-1", [R.make("t-b", {"x": 1}), R.make("t-x")]),
                     ("obs-2", [R.make("t-b", {"x": 1})])]
        plan = plan_dependencies(consumers)
        assert plan.order == ["t-a(x=1)", "t-b(x=1)", "t-x"]         # deps first, stable
        assert plan.nodes["t-b(x=1)"].consumers == ["obs-1", "obs-2"]  # one node, two consumers
        assert plan.closure("t-b(x=1)") == ["t-a(x=1)", "t-b(x=1)"]
        assert plan_dependencies(consumers).order == plan.order
        text = plan.explain()
        assert "shared with obs-1" in text and text.count("t-a(x=1)") == 2


def test_cycle_is_detected_with_its_path():
    a = Dummy("t-a", deps=[R.make("t-b")])
    b = Dummy("t-b", deps=[R.make("t-c")])
    c = Dummy("t-c", deps=[R.make("t-a")])
    with ireg.temporarily_registered(a, b, c):
        with pytest.raises(DependencyCycleError) as e:
            plan_dependencies([("obs", [R.make("t-a")])])
    assert e.value.path == ["t-a", "t-b", "t-c", "t-a"]


def test_contract_errors_block_only_their_consumers(tmp_path):
    a = Dummy("t-a", version=2)
    with ireg.temporarily_registered(a):
        plan, res = run([("pinned", [R.make("t-a", version=1)]),
                         ("unknown", [R.make("t-nope")]),
                         ("badparam", [R.make("t-a", {"y": 1})]),
                         ("fine", [R.make("t-a", version=2)])], tmp_path)
    assert "contract mismatch" in res["t-a"].reason or res["t-a"].status == S.COMPUTED
    # one node t-a: the first consumer pinned v1 → node error, shared by the key
    assert plan.nodes["t-nope"].error.startswith("unknown intermediate")
    assert "unknown parameter" in plan.nodes["t-a(y=1)"].error
    for consumer in ("unknown", "badparam"):
        _, status, why = consumer_dependencies(plan, res, consumer)
        assert status == S.BLOCKED and why


# ═══════════════════════════════════════════════════════════════════════════════
# Executor + store
# ═══════════════════════════════════════════════════════════════════════════════

def test_compute_once_reuse_across_invocations_and_identity_propagation(tmp_path):
    a = Dummy("t-a")
    b = Dummy("t-b", deps=[lambda p: R.make("t-a", {"x": p.get("x")})])
    with ireg.temporarily_registered(a, b):
        consumers = [("c1", [R.make("t-b", {"x": 1})]), ("c2", [R.make("t-b", {"x": 1})])]
        _, r1 = run(consumers, tmp_path)
        assert CALLS == {"t-a": 1, "t-b": 1}                          # once each
        assert r1["t-b(x=1)"].status == S.COMPUTED
        _, r2 = run(consumers, tmp_path)                              # a new invocation
        assert CALLS == {"t-a": 1, "t-b": 1}
        assert r2["t-b(x=1)"].status == S.CACHED
        assert r2["t-b(x=1)"].artifact_identity == r1["t-b(x=1)"].artifact_identity
        assert r2["t-b(x=1)"].dependencies[0]["artifact_identity"] == \
            r1["t-a(x=1)"].artifact_identity
        _, r3 = run([("c1", [R.make("t-b", {"x": 2})])], tmp_path)    # A changes → B changes
        assert r3["t-b(x=2)"].definition_identity != r1["t-b(x=1)"].definition_identity
    with ireg.temporarily_registered(Dummy("t-a", version=2), b):     # version bump
        _, r4 = run([("c1", [R.make("t-b", {"x": 1})])], tmp_path)
        assert r4["t-a(x=1)"].status == S.COMPUTED
        assert r4["t-b(x=1)"].definition_identity != r1["t-b(x=1)"].definition_identity


def test_identity_excludes_consumer_path_and_time_but_separates_inputs(tmp_path):
    from analysis.campaign.intermediates.base import IntermediateContext
    from analysis.campaign.models import TrajectoryView
    spec = Dummy("t-v")
    spec.scope = IntermediateScope.VIEW

    def ctx(view_ref, digest):
        v = TrajectoryView(kind="raw", requirements=TrajectoryRequirements(), cache_key=view_ref,
                           source_fingerprints=[SourceFingerprint("/t.xtc", 1, 1, "strong", digest)])
        return IntermediateContext(system=rec(), semantic_index=None, trajectory_view=v,
                                   topology_path="", structure_path=None, parameters={},
                                   work_dir=tmp_path)
    one = identities(spec, {"x": 1}, [], ctx("view-1", "d1"))
    other_view = identities(spec, {"x": 1}, [], ctx("view-2", "d1"))
    other_source = identities(spec, {"x": 1}, [], ctx("view-1", "d2"))
    assert one["definition_identity"] == other_view["definition_identity"]
    assert len({one["input_identity"], other_view["input_identity"],
                other_source["input_identity"]}) == 3
    assert len({one["artifact_identity"], other_view["artifact_identity"],
                other_source["artifact_identity"]}) == 3
    text = json.dumps(one["definition"])
    assert "/t.xtc" not in text and str(tmp_path) not in text          # no paths


def test_corruption_prevents_reuse_and_is_quarantined(tmp_path):
    a = Dummy("t-a")
    with ireg.temporarily_registered(a):
        _, r1 = run([("c", [R.make("t-a")])], tmp_path)
        f = Path(r1["t-a"].arrays[0].storage.path)
        f.write_bytes(f.read_bytes()[:-8] + b"\xff" * 8)              # tamper
        store = IntermediateStore(tmp_path)
        assert store.load(r1["t-a"].artifact_identity)[0] is None
        _, r2 = run([("c", [R.make("t-a")])], tmp_path)
        assert r2["t-a"].status == S.COMPUTED and CALLS["t-a"] == 2
        assert "changed after publication" in r2["t-a"].provenance["replaced_invalid_entry"]
        assert list((tmp_path / ".quarantine").iterdir())
        assert store.load(r1["t-a"].artifact_identity)[0] is not None  # clean again


def test_failed_dependency_blocks_only_dependents(tmp_path):
    bad = Dummy("t-bad", fail=True)
    user = Dummy("t-user", deps=[R.make("t-bad")])
    ok = Dummy("t-ok")
    with ireg.temporarily_registered(bad, user, ok):
        plan, res = run([("A", [R.make("t-user")]), ("B", [R.make("t-ok")])], tmp_path)
    assert res["t-bad"].status == S.FAILED and "deliberate failure" in res["t-bad"].reason
    assert res["t-user"].status == S.BLOCKED and "t-bad failed" in res["t-user"].reason
    assert res["t-ok"].status == S.COMPUTED
    assert consumer_dependencies(plan, res, "A")[1] == S.BLOCKED
    assert consumer_dependencies(plan, res, "B")[1] is None
    assert not [p for p in tmp_path.iterdir() if p.name.startswith(".building")]  # no residue


def test_republish_must_be_reproducible(tmp_path):
    noisy = Dummy("t-noise", noise=True)
    with ireg.temporarily_registered(noisy):
        _, r1 = run([("c", [R.make("t-noise")])], tmp_path)
        store = IntermediateStore(tmp_path)
        build = store.new_build_dir(r1["t-noise"].artifact_identity)
        np.save(build / "values.npy", np.asarray([9.0]))
        fake = r1["t-noise"]
        fake.arrays = [ResultArray("v", "dummy", axes=[Axis("i", AxisKind.INDEX)],
                                   storage=StorageRef(str(build / "values.npy"), "npy"))]
        with pytest.raises(IntermediateConflict, match="non-reproducible"):
            store.publish(build, fake)
        assert store.load(r1["t-noise"].artifact_identity)[0] is not None   # original intact


def test_concurrent_builders_publish_one_entry(tmp_path):
    slow = Dummy("t-slow")
    with ireg.temporarily_registered(slow):
        out = []
        threads = [threading.Thread(target=lambda: out.append(
            run([("c", [R.make("t-slow")])], tmp_path)[1]["t-slow"])) for _ in range(4)]
        [t.start() for t in threads]
        [t.join() for t in threads]
    assert CALLS["t-slow"] == 1                                       # the lock serialised them
    assert sorted(r.status for r in out) == [S.CACHED] * 3 + [S.COMPUTED]
    entries = [p for p in tmp_path.glob("*/*") if (p / "manifest.json").is_file()]
    assert len(entries) == 1 and not list(tmp_path.glob(".building-*"))


def test_structured_artifacts_round_trip(tmp_path):
    s = Dummy("t-map", structured=True)
    with ireg.temporarily_registered(s):
        _, r1 = run([("c", [R.make("t-map")])], tmp_path)
        _, r2 = run([("c", [R.make("t-map")])], tmp_path)
    (art,) = r2["t-map"].artifacts
    assert art.attrs["schema"] == "frame->atom->region" and Path(art.path).is_file()
    assert json.loads(Path(art.path).read_text()) == {"frame0": {"atom1": "upper"}}


def test_future_domain_graphs_need_no_planner_changes():
    """Architecture only: shapes future membrane / pore dependencies take."""
    frame = Dummy("test-membrane-reference-frame")
    leaflets = Dummy("test-leaflet-assignment", deps=[R.make("test-membrane-reference-frame")])
    axis = Dummy("test-pore-axis")
    with ireg.temporarily_registered(frame, leaflets, axis):
        plan = plan_dependencies([
            ("future-apl", [R.make("test-leaflet-assignment")]),
            ("future-thickness", [R.make("test-leaflet-assignment"),
                                  R.make("test-membrane-reference-frame")]),
            ("future-pore-radius", [R.make("test-pore-axis")])])
    assert plan.order == ["test-membrane-reference-frame", "test-leaflet-assignment",
                          "test-pore-axis"]
    assert plan.nodes["test-leaflet-assignment"].consumers == ["future-apl", "future-thickness"]


def test_no_placeholder_domain_intermediates_and_existing_observables_unchanged():
    from analysis.campaign.observables import registry as oreg
    assert ireg.ids() == ["selection-centre-series"]
    oreg.ensure_loaded()
    for oid in oreg.ids():
        assert oreg.get(oid).dependencies({}) == [], oid


def test_selection_centre_parameter_validation():
    spec = ireg.get("selection-centre-series")
    assert spec.validate({"selection": "component:ligand", "convention": "com"}) == []
    assert any("no default" in e for e in spec.validate({"selection": "component:ligand"}))
    assert spec.validate({"selection": "ligand", "convention": "com"})
    assert spec.validate({"selection": "component:ligand", "convention": "com", "pbc": "unwrap"})
    (schema,) = spec.output_schema({})
    assert schema.axis_signature() == ("time", "coordinate")


# ═══════════════════════════════════════════════════════════════════════════════
# Real GROMACS: consumers of selection-centre-series on the bundled membrane run
# ═══════════════════════════════════════════════════════════════════════════════

from analysis.campaign.observables.base import ObservableSpec  # noqa: E402


class CentreConsumer(ObservableSpec):
    """Test consumer: the z coordinate of a selection centre (reads the dependency)."""
    display_name = "test consumer"
    description = "test only"
    category = "test"
    required_components = ()
    purpose = "inter_component_geometry"

    def __init__(self, oid, selection="component:receptor", convention="com"):
        self.id, self.selection, self.convention = oid, selection, convention

    def _req(self):
        return R.make("selection-centre-series", {"selection": self.selection,
                                                  "convention": self.convention}, version=1)

    def dependencies(self, params=None):
        return [self._req()]

    def trajectory_requirements(self, ctx_params):
        return TrajectoryRequirements(rationale="reads only its dependency")

    def output_schema(self, params=None):
        return [ResultArray("z", "centre_z", "nm", axes=[Axis("time", AxisKind.TIME, "ps")])]

    def definition_evidence(self, ctx):
        return {"observable": self.id, "dependencies": self.dependency_evidence(ctx),
                "trajectory_view": {"cache_key": ctx.trajectory_view.cache_key}}

    def execute(self, ctx):
        from analysis.campaign.results import definition_token, read_column
        dep = ctx.intermediates[self._req().key]
        arr = dep.array("selection_center")
        z = np.load(arr.storage.path)[:, 2]
        t = read_column(arr.axes[0].values_ref)
        ctx.output_dir.mkdir(parents=True, exist_ok=True)
        out = ctx.output_dir / "z.xvg"
        out.write_text("".join(f"{a:.3f} {b:.5f}\n" for a, b in zip(t, z)))
        ev = self.definition_evidence(ctx)
        from analysis.campaign.fingerprint import fingerprint_file
        res = AnalysisResult(self.id, ctx.system.system_id, AnalysisStatus.SUCCESS,
                             definition_evidence=ev, definition_token=definition_token(ev))
        res.arrays = [ResultArray("z", "centre_z", "nm", axes=[Axis(
            "time", AxisKind.TIME, "ps", values_ref=StorageRef(str(out), "xvg", column=0),
            alignment=arr.axes[0].alignment)], storage=StorageRef(
            str(out), "xvg", column=1, fingerprint=fingerprint_file(out)),
            view_ref=ctx.trajectory_view.cache_key, definition_token=res.definition_token)]
        res.output_files = [str(out)]
        return res


SITE = """schema_version: simforge/annotations/v1
structural_annotation:
  residue_sets:
  - id: site
    kind: {kind}
    selection: {{numbering: topology, residues: "{res}", polymer_only: true}}
{prov}  - id: other
    kind: custom
    selection: {{numbering: topology, residues: "{other}", polymer_only: true}}
"""


def write_ann(run, res="177-190", other="140-150", kind="selected_subunit", derived=False):
    prov = "    provenance: {origin: derived}\n" if derived else ""
    (run / "simforge_annotations.yaml").write_text(
        SITE.format(res=res, other=other, kind=kind, prov=prov))


@pytest.fixture
def glp(tmp_path, monkeypatch):
    from analysis.campaign.observables import registry as oreg
    from analysis.campaign.orchestration.study_analyzer import clear_view_cache
    from tests.analysis.review.test_review import _run_dir
    oreg.ensure_loaded()
    for oid, sel in (("test-centre-a", "component:receptor"), ("test-centre-b", "component:receptor"),
                     ("test-centre-site", "annotation:site"),
                     ("test-centre-cog", "component:receptor")):
        monkeypatch.setitem(oreg._REGISTRY, oid, CentreConsumer(
            oid, sel, "cog" if oid.endswith("cog") else "com"))
    clear_view_cache()
    run = _run_dir(tmp_path)
    write_ann(run)
    return run


def _analyze(run, ids, **kw):
    from analysis.campaign.orchestration.study_analyzer import clear_view_cache, run_analyze
    clear_view_cache()
    res = run_analyze(run, ids, inspect_trajectories=False, **kw)
    return res, {r.analysis_id: r for r in res.results}


def _count_computes(monkeypatch):
    spec = ireg.get("selection-centre-series")
    calls = []
    real = type(spec).compute
    monkeypatch.setattr(type(spec), "compute",
                        lambda self, ctx: calls.append(ctx.parameters) or real(self, ctx))
    return calls


@requires_gmx
@requires_real_traj
def test_shared_intermediate_computed_once_and_reused_across_runs(glp, monkeypatch):
    calls = _count_computes(monkeypatch)
    res, by = _analyze(glp, ["test-centre-a", "test-centre-b", "rg"],
                       parameters={"rg": {"selection": "component:receptor"}})
    assert len(calls) == 1                                           # one trajectory computation
    a, b = by["test-centre-a"], by["test-centre-b"]
    assert a.status == b.status == AnalysisStatus.SUCCESS, (a.message, b.message)
    assert a.dependencies[0]["artifact_identity"] == b.dependencies[0]["artifact_identity"]
    assert a.dependencies[0]["status"] == "computed" and by["rg"].dependencies == []
    plan = next(iter(res.dependency_plans.values()))
    assert "shared with test-centre-a" in plan["explain"]
    entry = Path(glp / "simforge_analysis" / "intermediates")
    manifests = list(entry.glob("*/*/manifest.json"))
    assert len(manifests) == 1                                       # one persisted artifact
    m = json.loads(manifests[0].read_text())
    (arr,) = m["arrays"]
    assert [x["kind"] for x in arr["axes"]] == ["time", "coordinate"] and arr["unit"] == "nm"
    assert arr["attrs"]["convention"] == "com" and arr["attrs"]["shape"] == [51, 3]
    assert arr["axes"][0]["alignment"]["mode"] == FrameAlignmentMode.ONE_ROW_PER_FRAME
    assert m["definition_evidence"]["mass_source"] == "tpr masses"
    assert m["definition_evidence"]["gmx_flags"] == ["-rmpbc", "-nopbc"]
    assert not list(glp.rglob("whole.xtc"))                          # no derived trajectory
    calls.clear()
    res2, by2 = _analyze(glp, ["test-centre-a", "test-centre-b"])    # a new invocation
    assert calls == []                                               # intermediate reused
    assert {r.status for r in by2.values()} == {AnalysisStatus.CACHED}   # Phase 8, via deps
    assert by2["test-centre-a"].dependencies[0]["status"] == "cached"
    _, dry = _analyze(glp, ["test-centre-a"], dry_run=True)
    assert dry["test-centre-a"].dependencies[0]["status"] == "cached" and calls == []


@requires_gmx
@requires_real_traj
def test_parameter_and_annotation_identity(glp, monkeypatch):
    calls = _count_computes(monkeypatch)
    _, by = _analyze(glp, ["test-centre-a", "test-centre-cog", "test-centre-site"])
    ids = {k: r.dependencies[0]["artifact_identity"] for k, r in by.items()}
    assert len(set(ids.values())) == 3 and len(calls) == 3           # COM ≠ COG ≠ site
    write_ann(glp, other="140-151")                                  # unrelated annotation
    calls.clear()
    _, by2 = _analyze(glp, ["test-centre-a", "test-centre-site"])
    assert calls == [] and by2["test-centre-site"].dependencies[0]["artifact_identity"] == ids[
        "test-centre-site"]
    write_ann(glp, res="177-191", other="140-151")                   # the used atom set
    _, by3 = _analyze(glp, ["test-centre-a", "test-centre-site"])
    assert len(calls) == 1 and calls[0]["selection"] == "annotation:site"
    assert by3["test-centre-site"].dependencies[0]["artifact_identity"] != ids["test-centre-site"]
    assert by3["test-centre-a"].dependencies[0]["status"] == "cached"


@requires_gmx
@requires_real_traj
def test_inactive_annotation_dependency_blocks_only_its_consumer(glp):
    write_ann(glp, derived=True)                                     # site is PROPOSED
    _, by = _analyze(glp, ["test-centre-site", "test-centre-a"])
    site = by["test-centre-site"]
    assert site.status == AnalysisStatus.REVIEW_REQUIRED
    assert "annotation 'site' is proposed" in site.message and "blocked" in site.message
    assert by["test-centre-a"].status == AnalysisStatus.SUCCESS


@requires_gmx
@requires_real_traj
def test_review_session_records_dependency_references(glp):
    from analysis.review import DisplayRequest, ReviewRequest, parse_show, prepare_review
    prep = prepare_review(glp, ReviewRequest(parse_show(["test-centre-a", "test-centre-b"]),
                                             DisplayRequest.parse("raw")))
    deps = [o.dependencies for o in prep.dataset.observables]
    assert deps[0][0]["artifact_identity"] == deps[1][0]["artifact_identity"]
    assert {d[0]["status"] for d in deps} <= {"computed", "cached"}
    stores = list((glp / "simforge_analysis" / "review" / "intermediates").glob("*/*/manifest.json"))
    assert len(stores) == 1                                          # shared by both instances
