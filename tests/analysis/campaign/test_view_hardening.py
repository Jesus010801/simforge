"""Phase 3 — reproducible, identifiable, safely reusable derived trajectory views.

Most tests drive ``build_view`` with a tiny *fake* ``gmx`` executable so that
identity / cache / atomicity / cleanup logic is tested deterministically; the
real-GROMACS test at the end validates the same guarantees on the bundled
trajectory.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from analysis.campaign.models import (
    SemanticIndex, SemanticIndexGroup, TrajectoryRequirements, ValidationState,
    ViewBuildStatus, ViewCacheStatus,
)
from analysis.campaign.trajectory import requirements as reqs
from analysis.campaign.trajectory.preprocessor import (
    MANIFEST_NAME, VIEW_MANIFEST_SCHEMA, _identity_evidence, _timeline_args,
    build_view, load_view_time_index, view_identity,
)
from tests.analysis.campaign.conftest import (
    REAL_GRO, REAL_TPR, REAL_XTC, requires_gmx, requires_real_traj,
)

REPO = Path(__file__).resolve().parents[3]

# ═══════════════════════════════════════════════════════════════════════════════
# Fake gmx
# ═══════════════════════════════════════════════════════════════════════════════
#
# Fake trajectories are text files with one "TIME <t>" line per frame.
# `trjconv -o X.xtc` copies the input and appends an "OP ..." line (so each
# step's output differs); `trjconv -o X.gro` (the Phase 1 time-index probe)
# writes one-atom GRO frames with the TIME values.  Env knobs:
#   FAKE_GMX_FAIL_ON=<flag>   write a partial output then exit 1 when <flag> is in argv
#   FAKE_GMX_SHIFT_ON=<flag>  add 1000 ps to every time when <flag> is in argv
#   FAKE_GMX_LOG=<file>       append one JSON line per invocation

_FAKE = r'''#!{python}
import json, os, sys
args = sys.argv[1:]
if args and args[0] == "--version":
    print("GROMACS version:    {version}")
    sys.exit(0)
stdin = sys.stdin.read() if not sys.stdin.isatty() else ""
log = os.environ.get("FAKE_GMX_LOG")
if log:
    with open(log, "a") as fh:
        fh.write(json.dumps({{"args": args, "stdin": stdin}}) + "\n")
if not args or args[0] != "trjconv":
    sys.exit(0)
f = args[args.index("-f") + 1]
o = args[args.index("-o") + 1]
fail_on = os.environ.get("FAKE_GMX_FAIL_ON")
if fail_on and fail_on in args and not o.endswith(".gro"):
    open(o, "w").write("PARTIAL")
    sys.exit(1)
lines = open(f).read().splitlines()
times = [float(l.split()[1]) for l in lines if l.startswith("TIME")]
if o.endswith(".gro"):
    with open(o, "w") as out:
        for i, t in enumerate(times):
            out.write(f"probe t= {{t:.5f}} step= {{i}}\n1\n    1PRB     PX    1   0.000   0.000   0.000\n   3.00000   3.00000   3.00000\n")
    if times:
        sys.stderr.write(f"Last written: frame {{len(times) - 1}} time {{times[-1]:.3f}}\n")
    sys.exit(0)
shift = os.environ.get("FAKE_GMX_SHIFT_ON")
delta = 1000.0 if shift and shift in args else 0.0
with open(o, "w") as out:
    for t in times:
        out.write(f"TIME {{t + delta}}\n")
    for l in lines:
        if l.startswith("OP"):
            out.write(l + "\n")
    out.write("OP " + " ".join(a for a in args if a not in (f, o)) + "\n")
'''


def make_fake_gmx(dirpath: Path, version: str = "fake-1.0") -> str:
    p = dirpath / f"gmx_{version}"
    p.write_text(_FAKE.format(python=sys.executable, version=version))
    p.chmod(0o755)
    return str(p)


_NDX = """[ System ]
1 2 3 4 5 6 7 8
[ Protein ]
1 2 3 4
[ Unused ]
7 8
"""


@pytest.fixture
def env(tmp_path, monkeypatch):
    src = tmp_path / "inputs"
    src.mkdir()
    traj = src / "md.xtc"
    traj.write_text("TIME 50000\nTIME 50020\nTIME 50040\n")
    tpr = src / "md.tpr"
    tpr.write_bytes(b"TPR-v1")
    gro = src / "md.gro"
    gro.write_text("ref\n1\n    1A     B    1   0 0 0\n 1 1 1\n")
    ndx = src / "semantic.ndx"
    ndx.write_text(_NDX)
    log = tmp_path / "gmx.log"
    monkeypatch.setenv("FAKE_GMX_LOG", str(log))
    monkeypatch.delenv("FAKE_GMX_FAIL_ON", raising=False)
    monkeypatch.delenv("FAKE_GMX_SHIFT_ON", raising=False)
    return {"traj": traj, "tpr": tpr, "gro": gro, "ndx": ndx, "log": log,
            "gmx": make_fake_gmx(tmp_path), "tmp": tmp_path}


def _index(ndx: Path) -> SemanticIndex:
    return SemanticIndex(path=str(ndx), groups=[
        SemanticIndexGroup(name=n, n_atoms=0) for n in ("System", "Protein", "Unused")])


def build(env, req=None, *, work="work", **kw):
    params = dict(
        requirements=req or reqs.fit_to("Protein", "t"),
        trajectory_paths=[str(env["traj"])], topology_path=str(env["tpr"]),
        structure_path=str(env["gro"]), semantic_index=_index(env["ndx"]),
        source_fingerprints=[], topology_fingerprint=None,
        work_dir=env["tmp"] / work, gmx=env["gmx"])
    params.update(kw)
    return build_view(**params)


def _calls(env) -> list[dict]:
    if not env["log"].is_file():
        return []
    return [json.loads(l) for l in env["log"].read_text().splitlines()]


def _sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ═══════════════════════════════════════════════════════════════════════════════
# Identity
# ═══════════════════════════════════════════════════════════════════════════════

def test_identity_stable_and_output_dir_independent(env):
    a = build(env, work="w1")
    b = build(env, work="w2")
    assert a.build_status == b.build_status == ViewBuildStatus.BUILT
    assert a.cache_key == b.cache_key and len(a.cache_key) == 32


def test_changed_used_group_atoms_changes_identity(env):
    a = build(env)
    env["ndx"].write_text(_NDX.replace("[ Protein ]\n1 2 3 4", "[ Protein ]\n1 2 3 4 5"))
    b = build(env)
    assert a.cache_key != b.cache_key
    assert b.cache_status == ViewCacheStatus.MISS and not b.reused


def test_changed_unused_group_keeps_identity(env):
    a = build(env)
    env["ndx"].write_text(_NDX.replace("[ Unused ]\n7 8", "[ Unused ]\n6 7 8"))
    b = build(env)
    assert a.cache_key == b.cache_key and b.reused


def test_changed_source_topology_reference_change_identity(env):
    base = build(env).cache_key
    env["traj"].write_text("TIME 50000\nTIME 50020\nTIME 50060\n")
    k_traj = build(env).cache_key
    env["tpr"].write_bytes(b"TPR-v2")
    k_tpr = build(env).cache_key
    env["gro"].write_text("ref2\n1\n    1A     B    1   0 0 0\n 1 1 1\n")
    k_gro = build(env).cache_key
    assert len({base, k_traj, k_tpr, k_gro}) == 4


def test_changed_operations_change_identity(env):
    keys = {build(env, r).cache_key for r in (
        reqs.whole_only("w"), reqs.nojump("n"), reqs.fit_to("Protein", "f"),
        reqs.centered_on("Protein", "c"))}
    assert len(keys) == 4


def test_gromacs_version_changes_identity(env):
    a = build(env)
    b = build(env, gmx=make_fake_gmx(env["tmp"], "fake-2.0"))
    assert a.cache_key != b.cache_key


def test_timeline_args_are_part_of_identity(env):
    from analysis.campaign.fingerprint import fingerprint_file

    def step(dt):
        args = ["trjconv", "-s", str(env["tpr"]), "-f", str(env["traj"]),
                "-o", "/any/where/out.xtc", "-dt", dt]
        return {"op": "stride", "full_args": args, "stdin": "System", "uses_index": False}

    fps = {str(env["tpr"]): fingerprint_file(env["tpr"])}
    src = [fingerprint_file(env["traj"])]
    r = TrajectoryRequirements()
    e40, _ = _identity_evidence(r, "whole", [step("40")], src, fps, None, "v")
    e80, _ = _identity_evidence(r, "whole", [step("80")], src, fps, None, "v")
    assert view_identity(e40) != view_identity(e80)
    assert "/any/where" not in json.dumps(e40)                   # no output paths
    assert str(env["tmp"]) not in json.dumps(e40)                # no machine-local paths
    assert _timeline_args([step("40")]) == [["stride", "-dt", "40"]]


# ═══════════════════════════════════════════════════════════════════════════════
# Command order (locked) and manifest
# ═══════════════════════════════════════════════════════════════════════════════

def test_command_order_and_arguments_locked(env):
    req = TrajectoryRequirements(requires_whole_molecules=True, requires_nojump=True,
                                 centering_target="Protein", fit_selection="Protein")
    v = build(env, req, dry_run=True)
    entry = env["tmp"] / "work" / v.cache_key
    t, g, n = str(env["tpr"]), str(env["gro"]), str(env["ndx"])
    assert [o.command for o in v.operations] == [
        ["gmx", "trjconv", "-s", t, "-f", str(env["traj"]), "-o", str(entry / "whole.xtc"),
         "-pbc", "whole"],
        ["gmx", "trjconv", "-s", t, "-f", str(entry / "whole.xtc"),
         "-o", str(entry / "nojump.xtc"), "-pbc", "nojump"],
        # -pbc mol needs a .tpr (GROMACS: "Option -pbc mol requires a .tpr file")
        ["gmx", "trjconv", "-s", t, "-f", str(entry / "nojump.xtc"),
         "-o", str(entry / "centered.xtc"), "-pbc", "mol", "-center", "-n", n],
        ["gmx", "trjconv", "-s", g, "-f", str(entry / "centered.xtc"),
         "-o", str(entry / "fitted.xtc"), "-fit", "rot+trans", "-n", n],
    ]
    assert [o.stdin for o in v.operations] == [
        "System", "System", "Protein\nSystem", "Protein\nSystem"]
    assert v.build_status == ViewBuildStatus.PLANNED and not entry.exists()


def test_center_uses_tpr_for_pbc_mol_and_structure_only_without_one(env):
    """GROMACS refuses ``-pbc mol`` without a .tpr; fitting keeps its reference."""
    req = TrajectoryRequirements(centering_target="Protein")
    v = build(env, req, dry_run=True)
    (op,) = v.operations
    assert op.command[op.command.index("-s") + 1] == str(env["tpr"])
    top = env["tmp"] / "inputs" / "topol.top"
    top.write_text("; no tpr\n")
    v = build(env, req, dry_run=True, topology_path=str(top), work="w2")
    (op,) = v.operations
    assert op.command[op.command.index("-s") + 1] == str(env["gro"])   # nothing better exists


def test_executed_commands_match_plan(env):
    req = TrajectoryRequirements(requires_whole_molecules=True, requires_nojump=True,
                                 centering_target="Protein", fit_selection="Protein")
    build(env, req)
    tr = [c for c in _calls(env) if c["args"][0] == "trjconv"
          and not c["args"][c["args"].index("-o") + 1].endswith(".gro")]
    n = str(env["ndx"])
    assert [_tail(c) for c in tr] == [
        ["-pbc", "whole"], ["-pbc", "nojump"],
        ["-pbc", "mol", "-center", "-n", n], ["-fit", "rot+trans", "-n", n]]
    assert [c["args"][1:3] for c in tr] == [["-s", str(env["tpr"])], ["-s", str(env["tpr"])],
                                            ["-s", str(env["tpr"])], ["-s", str(env["gro"])]]
    assert tr[0]["args"][4] == str(env["traj"])
    assert [c["stdin"] for c in tr] == ["System\n", "System\n", "Protein\nSystem\n",
                                        "Protein\nSystem\n"]


def _tail(call):
    a = call["args"]
    return a[a.index("-o") + 2:]


def test_manifest_contents(env):
    v = build(env)
    m = json.loads(Path(v.manifest_path).read_text())
    assert m["schema_version"] == VIEW_MANIFEST_SCHEMA and m["status"] == "complete"
    assert m["view_identity"] == v.cache_key
    assert m["gromacs_version"] == "fake-1.0"
    assert [o["operation"] for o in m["operations"]] == ["make_whole", "fit"]
    assert [o["stdin"] for o in m["operations"]] == ["System", "Protein\nSystem"]
    assert all(o["returncode"] == 0 for o in m["operations"])
    assert m["sources"]["trajectories"][0]["digest"]
    assert set(m["sources"]["structures"]) == {"md.tpr", "md.gro"}
    assert m["sources"]["index"]["digest"]
    prot = [g for g in m["groups"] if g["name"] == "Protein"][0]
    from analysis.campaign.results import atom_set_hash
    assert prot["atoms_sha256"] == atom_set_hash([1, 2, 3, 4])
    assert m["output"]["relative_path"] == "fitted.xtc"
    assert m["output"]["fingerprint"]["digest"] == v.output_fingerprint.digest
    assert m["timeline_args"] == []
    ti = m["time_index"]
    assert ti["n_frames"] == 3 and ti["start_time_ps"] == 50000.0
    assert "times_ps" not in json.dumps(m)                        # no arrays copied
    # round-trip
    assert json.loads(json.dumps(m)) == m


# ═══════════════════════════════════════════════════════════════════════════════
# Cache validation
# ═══════════════════════════════════════════════════════════════════════════════

def _trjconv_count(env) -> int:
    return sum(1 for c in _calls(env) if c["args"][0] == "trjconv")


def test_valid_entry_is_reused_without_running_gmx(env):
    a = build(env)
    n = _trjconv_count(env)
    b = build(env)
    assert b.reused and b.build_status == ViewBuildStatus.REUSED
    assert b.cache_status == ViewCacheStatus.HIT and _trjconv_count(env) == n
    assert b.path == a.path and b.output_fingerprint.digest == a.output_fingerprint.digest
    assert [o.operation for o in b.operations] == ["make_whole", "fit"]
    assert b.time_index_ref["n_frames"] == 3


@pytest.mark.parametrize("damage,reason", [
    (lambda e: (e / MANIFEST_NAME).unlink(), "manifest missing"),
    (lambda e: (e / MANIFEST_NAME).write_text("{not json"), "corrupt"),
    (lambda e: (e / "fitted.xtc").write_text("TIME 1\nTAMPERED\n"), "changed after"),
    (lambda e: (e / "fitted.xtc").unlink(), "missing"),
    (lambda e: _edit_manifest(e, schema_version="simforge/trajectory-view-manifest/v0"),
     "schema"),
    (lambda e: _edit_manifest(e, view_identity="0" * 32), "identity"),
    (lambda e: _edit_manifest(e, status="building"), "status"),
])
def test_invalid_entry_is_rejected_and_regenerated(env, damage, reason):
    a = build(env)
    entry = Path(a.path).parent
    damage(entry)
    b = build(env)
    assert not b.reused and b.build_status == ViewBuildStatus.BUILT
    assert b.cache_status == ViewCacheStatus.INVALID
    assert any(reason in r for r in b.cache_rejections), b.cache_rejections
    assert any(w.code == "view_cache_rejected" for w in b.warnings)
    c = build(env)                                   # regenerated entry is valid again
    assert c.reused


def _edit_manifest(entry: Path, **changes):
    m = json.loads((entry / MANIFEST_NAME).read_text())
    m.update(changes)
    (entry / MANIFEST_NAME).write_text(json.dumps(m))


def test_changed_source_does_not_reuse_old_entry(env):
    a = build(env)
    env["traj"].write_text("TIME 1\nTIME 2\n")
    b = build(env)
    assert b.cache_key != a.cache_key and not b.reused
    assert Path(a.path).is_file()                    # old entry untouched, just not used


def test_force_rebuilds(env):
    build(env)
    b = build(env, force=True)
    assert b.build_status == ViewBuildStatus.BUILT and b.cache_status == ViewCacheStatus.BYPASSED


# ═══════════════════════════════════════════════════════════════════════════════
# Atomic completion and intermediates
# ═══════════════════════════════════════════════════════════════════════════════

def _entries(work: Path) -> list[str]:
    return sorted(p.name for p in work.iterdir()) if work.exists() else []


def test_failure_leaves_no_entry_and_no_temp(env, monkeypatch):
    monkeypatch.setenv("FAKE_GMX_FAIL_ON", "-fit")
    v = build(env)
    assert not v.safe and v.path is None and v.build_status == ViewBuildStatus.FAILED
    assert any(w.code == "preprocessing_failed" for w in v.warnings)
    work = env["tmp"] / "work"
    assert v.cache_key not in _entries(work)
    assert not [n for n in _entries(work) if n.startswith(".building-")]
    monkeypatch.delenv("FAKE_GMX_FAIL_ON")
    ok = build(env)
    assert ok.build_status == ViewBuildStatus.BUILT and ok.cache_status == ViewCacheStatus.MISS


def test_entry_holds_only_final_product_and_manifest(env):
    v = build(env, TrajectoryRequirements(requires_whole_molecules=True, requires_nojump=True,
                                          fit_selection="Protein"))
    entry = Path(v.path).parent
    assert sorted(p.name for p in entry.iterdir()) == ["fitted.xtc", MANIFEST_NAME]
    m = json.loads((entry / MANIFEST_NAME).read_text())
    assert set(m["intermediates_removed"]) >= {"whole.xtc", "nojump.xtc"}


def test_abandoned_build_dir_is_never_a_hit_and_is_swept(env):
    v = build(env, dry_run=True)
    work = env["tmp"] / "work"
    work.mkdir(parents=True, exist_ok=True)
    stale = work / f".building-{v.cache_key}-999999999-abc"
    stale.mkdir()
    (stale / "fitted.xtc").write_text("TIME 50000\n")
    (stale / MANIFEST_NAME).write_text("{}")
    b = build(env)
    assert b.build_status == ViewBuildStatus.BUILT and b.cache_status == ViewCacheStatus.MISS
    assert not stale.exists()


def test_other_keys_build_dirs_are_not_swept(env):
    work = env["tmp"] / "work"
    work.mkdir(parents=True)
    other = work / ".building-deadbeef-999999999-abc"
    other.mkdir()
    build(env)
    assert other.exists()


def test_legacy_entry_without_manifest_is_not_trusted(env):
    v = build(env, dry_run=True)
    work = env["tmp"] / "work"
    entry = work / v.cache_key
    entry.mkdir(parents=True)
    (entry / "fitted.xtc").write_text("TIME 50000\nTIME 50020\nTIME 50040\n")
    legacy_old_key = work / "0123456789abcdef"            # pre-Phase-3 16-hex key
    legacy_old_key.mkdir()
    (legacy_old_key / "fitted.xtc").write_text("old")
    b = build(env)
    assert not b.reused and b.cache_status == ViewCacheStatus.INVALID
    assert any("manifest missing" in r for r in b.cache_rejections)
    assert (legacy_old_key / "fitted.xtc").read_text() == "old"   # never deleted


# ═══════════════════════════════════════════════════════════════════════════════
# Source immutability
# ═══════════════════════════════════════════════════════════════════════════════

def test_sources_never_modified(env):
    before = {k: _sha(env[k]) for k in ("traj", "tpr", "gro", "ndx")}
    build(env, TrajectoryRequirements(requires_whole_molecules=True, requires_nojump=True,
                                      centering_target="Protein", fit_selection="Protein"))
    build(env)
    assert {k: _sha(env[k]) for k in ("traj", "tpr", "gro", "ndx")} == before
    inputs = {str(env[k]) for k in ("traj", "tpr", "gro", "ndx")}
    for c in _calls(env):
        a = c["args"]
        if "-o" in a:
            assert a[a.index("-o") + 1] not in inputs


# ═══════════════════════════════════════════════════════════════════════════════
# Per-view FrameTimeIndex
# ═══════════════════════════════════════════════════════════════════════════════

def test_view_has_its_own_time_index(env):
    v = build(env)
    ref = v.time_index_ref
    assert ref["n_frames"] == 3 and ref["state"] == ValidationState.VALID
    assert ref["trajectory_digest"] == v.output_fingerprint.digest
    assert ref["same_timeline_as_source"] is True                 # demonstrated, not assumed
    idx = load_view_time_index(v)
    assert idx.times_ps == [50000.0, 50020.0, 50040.0]


def test_changed_timeline_is_detected_not_inherited(env, monkeypatch):
    monkeypatch.setenv("FAKE_GMX_SHIFT_ON", "-fit")
    v = build(env)
    assert v.time_index_ref["same_timeline_as_source"] is False
    assert load_view_time_index(v).times_ps == [51000.0, 51020.0, 51040.0]


def test_raw_view_has_identity_but_no_derived_manifest(env):
    v = build(env, reqs.raw())
    assert v.build_status == ViewBuildStatus.RAW and v.path == str(env["traj"])
    assert v.cache_key and v.manifest_path is None


# ═══════════════════════════════════════════════════════════════════════════════
# In-process view cache and import hygiene
# ═══════════════════════════════════════════════════════════════════════════════

def test_in_process_view_cache_keyed_by_index_content(tmp_path, monkeypatch):
    from analysis.campaign.models import CampaignRunResult, SystemRecord
    from analysis.campaign.orchestration import study_analyzer as sa

    calls = []
    from types import SimpleNamespace
    monkeypatch.setattr(sa, "build_view",
                        lambda **kw: calls.append(kw) or SimpleNamespace(warnings=[]))
    ndx = tmp_path / "i.ndx"
    ndx.write_text("[ Protein ]\n1 2 3\n")
    rec = SystemRecord(system_id="s", condition_id="c", replicate_id="r",
                       topology_path=str(tmp_path / "md.tpr"),
                       semantic_index=SemanticIndex(path=str(ndx)))
    res = CampaignRunResult(study_root="x", output_dir="y")
    req = reqs.fit_to("Protein", "f")
    from analysis.campaign.models import IntentSource
    from analysis.campaign.trajectory.policy import PolicyIntent
    intent = PolicyIntent(IntentSource.EXPLICIT_API)     # Phase 5: fit needs explicit intent here
    a = sa._resolve_view(rec, req, tmp_path, "gmx", False, res, intent=intent)
    assert sa._resolve_view(rec, req, tmp_path, "gmx", False, res, intent=intent) is a
    ndx.write_text("[ Protein ]\n1 2 3 4\n")
    b = sa._resolve_view(rec, req, tmp_path, "gmx", False, res, intent=intent)
    assert b is not a and len(calls) == 2


def test_rmsd_module_imports_directly():
    out = subprocess.run([sys.executable, "-c", "import analysis.campaign.observables.rmsd"],
                         cwd=REPO, capture_output=True, text=True)
    assert out.returncode == 0, out.stderr


# ═══════════════════════════════════════════════════════════════════════════════
# Real GROMACS
# ═══════════════════════════════════════════════════════════════════════════════

@requires_gmx
@requires_real_traj
def test_real_derived_view(tmp_path):
    natoms = int(REAL_GRO.read_text().splitlines()[1])
    ndx = tmp_path / "semantic.ndx"

    def write_ndx(n_prot):
        ndx.write_text("[ System ]\n" + " ".join(map(str, range(1, natoms + 1))) + "\n"
                       "[ Protein ]\n" + " ".join(map(str, range(1, n_prot + 1))) + "\n")

    write_ndx(200)
    before = {p: _sha(p) for p in (REAL_XTC, REAL_TPR, REAL_GRO)}
    idx = SemanticIndex(path=str(ndx), groups=[SemanticIndexGroup("System", natoms),
                                               SemanticIndexGroup("Protein", 200)])
    kw = dict(requirements=reqs.fit_to("Protein", "f"), trajectory_paths=[str(REAL_XTC)],
              topology_path=str(REAL_TPR), structure_path=str(REAL_GRO),
              semantic_index=idx, source_fingerprints=[], topology_fingerprint=None,
              work_dir=tmp_path / "views")
    a = build_view(**kw)
    assert a.build_status == ViewBuildStatus.BUILT, [w.message for w in a.warnings]
    assert Path(a.path).is_file() and Path(a.manifest_path).is_file()
    m = json.loads(Path(a.manifest_path).read_text())
    assert m["gromacs_version"] and [o["operation"] for o in m["operations"]] == ["make_whole", "fit"]
    assert a.output_fingerprint.digest == m["output"]["fingerprint"]["digest"]
    ti = load_view_time_index(a)
    assert ti.n_frames == 51 and ti.start_time_ps == 0.0 and ti.end_time_ps == 1000.0
    assert a.time_index_ref["same_timeline_as_source"] is True

    b = build_view(**kw)
    assert b.reused and b.cache_status == ViewCacheStatus.HIT and b.cache_key == a.cache_key

    write_ndx(201)                                   # used group atoms changed
    c = build_view(**kw)
    assert c.cache_key != a.cache_key and c.cache_status == ViewCacheStatus.MISS

    assert {p: _sha(p) for p in (REAL_XTC, REAL_TPR, REAL_GRO)} == before
