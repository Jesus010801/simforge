"""Phase 12 — declarative review profiles resolved through generic preparation."""
from __future__ import annotations

import ast
import json
from pathlib import Path

import pytest
import yaml

from analysis.campaign.models import AnnotationRecord, MolecularComponent, SystemRecord
from analysis.review import DisplayRequest, ReviewRequest, open_review_session, parse_show
from analysis.review.profiles import (
    ProfileError, builtin_profiles, get_profile, list_profiles, merge_request, parse_profile,
    resolve_profile,
)
from tests.analysis.campaign.conftest import REPO_ROOT, requires_gmx, requires_real_traj

PROFILES_DIR = REPO_ROOT / "analysis" / "review" / "profiles"


def doc(**over):
    d = {"schema": "simforge/review-profile/v1", "id": "custom", "version": 1,
         "description": "a test profile",
         "observables": [{"observable": "rg", "parameters": {"selection": "component:receptor"}}]}
    d.update(over)
    return d


def write(path: Path, d: dict) -> Path:
    path.write_text(yaml.safe_dump(d, sort_keys=False))
    return path


# ═══════════════════════════════════════════════════════════════════════════════
# Schema
# ═══════════════════════════════════════════════════════════════════════════════

def test_builtin_profiles_are_valid_data():
    ps = {p.id: p for p in list_profiles()}
    assert set(ps) == {"general", "soluble-protein", "protein-ligand", "membrane-protein"}
    assert "membrane-channel" not in ps                     # deferred: would be hollow
    for p in ps.values():
        assert p.display == "raw"                           # no transformation by default
    pl = ps["protein-ligand"]
    assert "not binding-site metrics" in " ".join(pl.description.split())
    assert [o.request.params.get("selection_b") for o in pl.optional
            if o.annotation_kind] == ["annotation:{match}"] * 2


@pytest.mark.parametrize("change,needle", [
    ({"colour": "red"}, "unknown key"),
    ({"schema": "simforge/review-profile/v2"}, "schema"),
    ({"id": "Bad Id"}, "must match"),
    ({"version": 0}, "version"),
    ({"display": "spin:Receptor"}, "display"),
    ({"requirements": {"components": ["ligandz"]}}, "unknown component"),
    ({"requirements": {"components": ["ligand"], "absent": ["x"]}}, "unknown key"),
    ({"observables": [{"observable": "pore-radius"}]}, "not a registered observable"),
    ({"observables": [{"observable": "rg", "parameters": {"selection": "ligand"}}]}, "must be"),
    ({"observables": [{"observable": "rg"}]}, "needs parameter"),
    ({"observables": [{"observable": "rmsd-receptor", "parameters": {"x": 1}}]}, "takes no"),
    ({"observables": [{"observable": "rg", "parameters": {"selection": "component:receptor"},
                       "command": "gmx trjconv"}]}, "unknown key"),
    ({"optional_observables": [{"when": {"annotation_kind": "binding_site"}, "observable": "rg",
                                "parameters": {"selection": "annotation:site"}}]}, "{match}"),
    ({"optional_observables": [{"when": {"annotation": "x", "annotation_kind": "y"},
                                "observable": "rg",
                                "parameters": {"selection": "annotation:x"}}]}, "exactly one"),
    ({"observables": [{"observable": "rg", "parameters": {"selection": "annotation:{match}"}}]},
     "only valid with annotation_kind"),
])
def test_strict_schema(change, needle):
    with pytest.raises(ProfileError, match=needle):
        parse_profile(doc(**change))


def test_custom_file_and_composition(tmp_path):
    base = parse_profile(doc(id="base", display="raw",
                             requirements={"components": ["receptor"]}), source="builtin")
    child = write(tmp_path / "child.yaml", doc(
        id="child", extends=["base"],
        observables=[{"observable": "rg", "parameters": {"selection": "component:receptor"}},
                     {"observable": "rmsd-receptor"}]))
    p = get_profile(str(child), builtins={"base": base})
    assert p.composition == ["base", "child"] and p.display == "raw"
    assert [o.instance_id for o in p.observables] == ["rg(selection=component:receptor)",
                                                      "rmsd-receptor"]
    assert p.duplicates_dropped == ["child: rg(selection=component:receptor)"]   # explicit
    assert p.components == ["receptor"] and p.source == "file"
    a = parse_profile(doc(id="a", extends=["b"]))
    b = parse_profile(doc(id="b", extends=["a"]))
    with pytest.raises(ProfileError, match="cycle"):
        get_profile("a", builtins={"a": a, "b": b})
    with pytest.raises(ProfileError, match="unknown profile"):
        get_profile("a", builtins={"a": a})
    x = parse_profile(doc(id="x", display="raw"))
    y = parse_profile(doc(id="y", display="whole"))
    z = parse_profile(doc(id="z", extends=["x", "y"]))
    with pytest.raises(ProfileError, match="disagree"):
        get_profile("z", builtins={"x": x, "y": y, "z": z})
    with pytest.raises(ProfileError, match="defined twice"):
        get_profile("x", builtins={"x": x}, extra={"x": y})
    both = parse_profile(doc(id="w", requirements={"components": ["membrane"],
                                                   "absent_components": ["membrane"]}))
    with pytest.raises(ProfileError, match="both required and absent"):
        get_profile("w", builtins={"w": both})


def test_duplicate_builtin_ids_are_refused(tmp_path):
    write(tmp_path / "a.yaml", doc(id="same"))
    write(tmp_path / "b.yaml", doc(id="same"))
    with pytest.raises(ProfileError, match="duplicate built-in"):
        builtin_profiles(tmp_path)


def test_identity_is_content_not_path_or_wording(tmp_path):
    one = write(tmp_path / "one.yaml", doc(description="first wording"))
    (tmp_path / "elsewhere").mkdir()
    two = tmp_path / "elsewhere" / "copy.yaml"
    two.write_text("# a comment\n" + yaml.safe_dump(doc(description="other words"),
                                                  sort_keys=True, default_flow_style=True))
    p1, p2 = get_profile(str(one)), get_profile(str(two))
    assert p1.definition_identity == p2.definition_identity
    changed = write(tmp_path / "v.yaml", doc(version=2))
    assert get_profile(str(changed)).definition_identity != p1.definition_identity
    param = write(tmp_path / "p.yaml", doc(observables=[
        {"observable": "rg", "parameters": {"selection": "component:ligand"}}]))
    assert get_profile(str(param)).definition_identity != p1.definition_identity


# ═══════════════════════════════════════════════════════════════════════════════
# Resolution against resolved scientific state (synthetic records)
# ═══════════════════════════════════════════════════════════════════════════════

def rec(*components, annotations=()):
    r = SystemRecord("s", "c", "r")
    r.components = [MolecularComponent(t, t, classification_state=st) for t, st in components]
    r.annotations = [AnnotationRecord(aid, "residue_set", kind, "user_yaml", state)
                     for aid, kind, state in annotations]
    return r


RES = "resolved"


def test_applicability_from_component_states():
    B = {p.id: p for p in list_profiles()}
    soluble = rec(("receptor", RES))
    assert resolve_profile(B["soluble-protein"], soluble).status == "applicable"
    assert resolve_profile(B["protein-ligand"], rec(("receptor", RES), ("ligand", RES))).status \
        == "applicable"
    none = resolve_profile(B["protein-ligand"], soluble)
    assert none.status == "not_applicable" and "no resolved ligand" in none.reasons[0]
    amb = resolve_profile(B["protein-ligand"], rec(("receptor", RES), ("ligand", "ambiguous")))
    assert amb.status == "review_required"
    mem = rec(("receptor", RES), ("membrane", RES))
    assert resolve_profile(B["membrane-protein"], mem).status == "applicable"
    assert resolve_profile(B["soluble-protein"], mem).status == "not_applicable"
    assert resolve_profile(B["general"], rec()).status == "applicable"   # no requirements;
    # its requests still face ObservableSpec.applicability during preparation


def test_optional_annotation_requests():
    pl = {p.id: p for p in list_profiles()}["protein-ligand"]
    base = [("receptor", RES), ("ligand", RES)]

    def resolved(*anns):
        r = resolve_profile(pl, rec(*base, annotations=anns))
        return r, {o["request"]: o for o in r.optional}, [q.instance_id for q in r.requests]
    r, opt, ids = resolved(("site1", "binding_site", "active"))
    assert "min-distance(selection_a=component:ligand,selection_b=annotation:site1)" in ids
    assert opt["min-distance(selection_a=component:ligand,selection_b=annotation:{match})"][
        "annotation"] == "site1"
    assert not any("selected_subunit" in i for i in ids)
    r, opt, ids = resolved()
    assert all(o["outcome"] == "unavailable" for o in r.optional) and len(ids) == 5
    assert "no ACTIVE binding_site" in r.optional[2]["reason"]
    r, opt, ids = resolved(("pocket", "binding_site", "proposed"))
    assert not any("pocket" in i for i in ids)                         # never activated
    assert "pocket (proposed)" in r.optional[2]["reason"]
    r, opt, ids = resolved(("siteA", "binding_site", "active"), ("siteB", "binding_site", "active"))
    kind = [o for o in r.optional if o["condition"].get("annotation_kind")]
    assert {o["outcome"] for o in kind} == {"review_required"}
    assert kind[0]["candidates"] == ["siteA", "siteB"] and len(ids) == 5   # no choice made
    r, opt, ids = resolved(("selected_subunit", "selected_subunit", "active"))
    assert "com-distance(selection_a=component:ligand,selection_b=annotation:selected_subunit)" \
        in ids and "rg(selection=annotation:selected_subunit)" in ids
    r, opt, ids = resolved(("selected_subunit", "selected_subunit", "unresolved"))
    assert "is unresolved, not ACTIVE" in r.optional[0]["reason"]


def test_merge_user_precedence_and_hide():
    p = {p.id: p for p in list_profiles()}["membrane-protein"]
    res = resolve_profile(p, rec(("receptor", RES), ("membrane", RES)))
    user = ReviewRequest(parse_show(["rg(selection=annotation:tm_1)"]),
                         DisplayRequest.parse("center:Receptor"))
    merged, ov = merge_request(res, user, hide=["rmsd-receptor"])
    ids = [o.instance_id for o in merged.observables]
    assert ids[-1] == "rg(selection=annotation:tm_1)" and "rmsd-receptor" not in ids
    assert merged.display.mode == "center" and not merged.display.intent  # still policy-judged
    assert ov == {"added": ["rg(selection=annotation:tm_1)"], "hidden": ["rmsd-receptor"],
                  "display_source": "user"}
    merged, ov = merge_request(res, ReviewRequest([], None))
    assert merged.display.mode == "raw" and ov["display_source"] == "profile"
    with pytest.raises(ProfileError, match="matches no request"):
        merge_request(res, ReviewRequest([], None), hide=["sasa-typo"])


# ═══════════════════════════════════════════════════════════════════════════════
# Structural guards: data → ReviewRequest → generic preparation
# ═══════════════════════════════════════════════════════════════════════════════

def test_no_profile_specific_code_paths():
    ids = {p.id for p in list_profiles()}
    for py in (REPO_ROOT / "analysis").rglob("*.py"):
        text = py.read_text()
        for pid in ids - {"general"}:                 # "general" is an ordinary word
            assert f'"{pid}"' not in text and f"'{pid}'" not in text, (py, pid)
    for py in PROFILES_DIR.glob("*.py"):
        tree = ast.parse(py.read_text())
        names = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name)} | \
                {n.attr for n in ast.walk(tree) if isinstance(n, ast.Attribute)}
        mods = {n.module for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)}
        assert not names & {"run_gmx", "run_analyze", "build_view", "plan_view", "eval", "exec",
                            "subprocess", "prepare_system", "diagnose_system"}, py
        assert not {m for m in mods if m and ("gmx" in m or "orchestration" in m)}, py
    obs_dir = REPO_ROOT / "analysis" / "campaign" / "observables"
    assert not any("profile" in p.read_text().lower().replace("density_profile", "")
                   for p in obs_dir.glob("*.py"))


# ═══════════════════════════════════════════════════════════════════════════════
# Real preparation (bundled GLP-1R membrane system)
# ═══════════════════════════════════════════════════════════════════════════════

def _prep(run, profile, show=(), display=None, **kw):
    from analysis.review import prepare_review
    return prepare_review(run, ReviewRequest(parse_show(list(show)), display), profile=profile, **kw)


@pytest.fixture(scope="module")
def glp(tmp_path_factory):
    import shutil
    from tests.analysis.campaign.conftest import REAL_XTC
    if not (REAL_XTC.is_file() and shutil.which("gmx")):
        pytest.skip("needs gmx and the bundled trajectory")
    from tests.analysis.review.test_review import _run_dir
    return _run_dir(tmp_path_factory.mktemp("prof"))


@requires_gmx
@requires_real_traj
def test_membrane_profile_prepares_generically(glp, monkeypatch):
    import analysis.campaign.orchestration.study_analyzer as sa
    import analysis.review.prepare as rp
    seen = []
    real = sa.run_analyze
    monkeypatch.setattr(sa, "run_analyze",
                        lambda root, ids, **kw: seen.append((ids, kw["parameters"])) or real(root, ids, **kw))
    prep = _prep(glp, get_profile("membrane-protein"))
    ds = prep.dataset
    # the profile became ordinary instances, prepared by the generic path
    assert [ids for ids, _ in seen] == [["rg"], ["sasa"], ["rmsd-receptor"]]
    assert seen[0][1] == {"rg": {"selection": "component:receptor"}}
    p = ds.provenance["profile"]
    assert (p["id"], p["version"], p["source"], p["status"]) == \
        ("membrane-protein", 1, "builtin", "applicable")
    assert p["overrides"]["display_source"] == "profile"
    assert ds.display_view["kind"] == "raw" and ds.display_view["decisions"] == []
    assert {o.state for o in ds.observables} == {"available"}
    assert ds.identity_evidence["profile"] == p["definition_identity"]
    assert "/" not in json.dumps(p["definition_identity"])
    assert not any("path" in k for k in p)                            # no source path
    # the dashboard shows it as metadata only
    from analysis.review.runtime.server import ReviewServer
    from analysis.review.runtime.viewer import FakeViewerAdapter
    from tests.analysis.review.test_runtime import _req
    with ReviewServer(prep.dataset_path, viewer=FakeViewerAdapter()) as srv:
        _, s = _req(srv, "/api/session")
        assert s["profile"]["id"] == "membrane-protein"
        _req(srv, "/api/sync/frame", {"frame": 10})
        assert srv.viewer.frame == 10                                   # sync unchanged


@requires_gmx
@requires_real_traj
def test_profile_session_and_result_reuse(glp, monkeypatch, tmp_path):
    first = _prep(glp, get_profile("membrane-protein"))
    import analysis.campaign.observables.generic as gen
    import analysis.campaign.observables.rmsd as rmsd

    def boom(*a, **k):
        raise AssertionError("no observable GROMACS execution expected")
    monkeypatch.setattr(gen, "run_gmx", boom)
    monkeypatch.setattr(rmsd, "run_gmx", boom)
    again = _prep(glp, get_profile("membrane-protein"))
    assert again.reused_session and again.dataset.session_id == first.dataset.session_id
    forced = _prep(glp, get_profile("membrane-protein"), force=True)
    assert {o.availability for o in forced.dataset.observables} == {"cached"}   # Phase 8
    # the same requests from a custom file with other wording: new profile identity,
    # but every result is still reused (the profile is not part of any observable definition)
    custom = write(tmp_path / "mine.yaml", doc(
        id="my-membrane", description="my words", requirements={"components": ["receptor"]},
        observables=[{"observable": "rg", "parameters": {"selection": "component:receptor"}},
                     {"observable": "sasa", "parameters": {"selection": "component:receptor",
                                                           "surface": "component:receptor"}}]))
    mine = _prep(glp, get_profile(str(custom)))
    assert {o.availability for o in mine.dataset.observables} == {"cached"}
    monkeypatch.undo()
    # a changed scientific parameter changes exactly that result and the session
    probe = write(tmp_path / "probe.yaml", {**yaml.safe_load(custom.read_text()), "observables": [
        {"observable": "rg", "parameters": {"selection": "component:receptor"}},
        {"observable": "sasa", "parameters": {"selection": "component:receptor",
                                              "surface": "component:receptor", "probe": 0.16}}]})
    changed = _prep(glp, get_profile(str(probe)))
    assert changed.dataset.session_id != mine.dataset.session_id
    by = {o.observable: o.availability for o in changed.dataset.observables}
    assert by == {"rg": "cached", "sasa": "computed"}
    assert open_review_session(mine.dataset_path)[1].valid


@requires_gmx
@requires_real_traj
def test_profile_cannot_override_policy_and_user_overrides_are_judged(glp, tmp_path):
    unsafe = write(tmp_path / "unsafe.yaml", doc(id="unsafe-display", display="fit:Receptor"))
    prep = _prep(glp, get_profile(str(unsafe)))
    dv = prep.dataset.display_view
    assert dv["status"] == "refused" and "fit.membrane_frame" in dv["reason"]
    assert prep.dataset.provenance["profile"]["overrides"]["display_source"] == "profile"
    centred = write(tmp_path / "centred.yaml", doc(id="centred-display", display="center:Receptor"))
    dv = _prep(glp, get_profile(str(centred))).dataset.display_view
    assert dv["status"] == "review_required"          # profile ≠ explicit intent
    assert {d["intent_source"] for d in dv["decisions"]} == {"auto"}
    user = _prep(glp, get_profile("membrane-protein"), display=DisplayRequest.parse("center:Receptor"))
    dv = user.dataset.display_view
    assert user.dataset.provenance["profile"]["overrides"]["display_source"] == "user"
    assert dv["status"] == "review_required"          # the user's choice is still judged
    both = _prep(glp, get_profile("membrane-protein"),
                 display=DisplayRequest.parse("center:Receptor", intent=True))
    assert both.dataset.display_view["status"] == "available"   # explicit user intent
    assert {d["intent_source"] for d in both.dataset.display_view["decisions"]} == {"user_flag"}


@requires_gmx
@requires_real_traj
def test_not_applicable_profile_is_refused_and_cli_explains(glp):
    import subprocess
    import sys
    from analysis.review.prepare import ReviewError
    with pytest.raises(ReviewError, match="no resolved ligand"):
        _prep(glp, get_profile("protein-ligand"))
    r = subprocess.run([sys.executable, "-m", "cli", "trajectory", "review", str(glp),
                        "--profile", "membrane-protein", "--hide", "sasa", "--dry-run"],
                       cwd=REPO_ROOT, capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "membrane-protein v1 (builtin)" in r.stdout and "hidden" in r.stdout
    assert "sasa" not in r.stdout.split("Observables")[1]
    r = subprocess.run([sys.executable, "-m", "cli", "trajectory", "review-profiles", str(glp),
                        "--json"], cwd=REPO_ROOT, capture_output=True, text=True, timeout=300)
    rows = {x["id"]: x["system"]["requirements_status"] for x in json.loads(r.stdout)}
    assert rows == {"general": "applicable", "membrane-protein": "applicable",
                    "protein-ligand": "not_applicable", "soluble-protein": "not_applicable"}
