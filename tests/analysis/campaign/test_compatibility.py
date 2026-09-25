"""Phase 8 — scientific compatibility and reuse of existing results."""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import pytest

from analysis.campaign.compatibility import (
    Candidate, CheckOutcome, CompatibilityStatus as CS, ProvenanceTier, definition_identity,
    evaluate, output_resolution, storage_integrity,
)
from analysis.campaign.gmx import run_gmx
from analysis.campaign.models import (
    AnalysisResult, AnalysisStatus, Axis, AxisKind, FrameAlignment, ResultArray, StorageRef,
)
from analysis.campaign.orchestration.study_analyzer import clear_view_cache, run_analyze
from tests.analysis.campaign.conftest import (
    REAL_GRO, REAL_TPR, REAL_XTC, requires_gmx, requires_real_traj,
)

GENERIC = ["rg", "sasa", "com-distance", "min-distance", "hbond-count"]


def sha(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ═══════════════════════════════════════════════════════════════════════════════
# Unit level: evaluation dimensions
# ═══════════════════════════════════════════════════════════════════════════════

def _native(tmp_path, times=(0.0, 20.0, 40.0), values=(1.0, 1.1, 1.2), *, evidence=None,
            view="v1", name="d.xvg"):
    from analysis.campaign.fingerprint import fingerprint_file
    f = tmp_path / name
    f.write_text('@    xaxis  label "Time (ps)"\n@    yaxis  label "Distance (nm)"\n'
                 + "".join(f"{t:10.3f} {v:8.3f}\n" for t, v in zip(times, values)))
    ev = evidence or {"observable": "com-distance", "selections": {"a": 1},
                      "backend": {"tool": "gmx distance", "version": "2025.2"},
                      "trajectory_view": {"kind": "raw", "cache_key": view}}
    arr = ResultArray(
        name="com_distance", quantity="com_distance", unit="nm",
        axes=[Axis("time", AxisKind.TIME, unit="ps",
                   values_ref=StorageRef(str(f), "xvg", column=0),
                   alignment=FrameAlignment("one_row_per_frame", trajectory_digest="dig",
                                            view_ref=view))],
        storage=StorageRef(str(f), "xvg", column=1, fingerprint=fingerprint_file(f)))
    res = AnalysisResult("com-distance", "s", "success", arrays=[arr], definition_evidence=ev,
                         definition_token="tok")
    return Candidate(source=str(f), tier=ProvenanceTier.NATIVE, result=res), ev


def _schema():
    return [ResultArray(name="com_distance", quantity="com_distance", unit="nm",
                        axes=[Axis("time", AxisKind.TIME, unit="ps")])]


def _eval(cand, ev, *, view="v1", timeline=(0.0, 20.0, 40.0)):
    return evaluate(cand, requested_evidence=ev, requested_view_ref=view,
                    requested_schema=_schema(), requested_timeline=list(timeline))


def test_exact_match_is_compatible_with_structured_report(tmp_path):
    cand, ev = _native(tmp_path)
    r = _eval(cand, ev)
    assert r.status == CS.COMPATIBLE
    dims = {c.dimension: c.outcome for c in r.checks}
    assert set(dims) == {"provenance", "definition", "backend", "view", "source", "timeline",
                         "output_schema", "output_integrity", "output_resolution"}
    assert set(dims.values()) == {CheckOutcome.MATCH}
    reso = next(c for c in r.checks if c.dimension == "output_resolution").evidence
    assert reso["measured"]["com_distance"] == {"format": "fixed", "decimals": 3,
                                                "resolution": 0.001}
    assert "DEFINITION" in r.summary().upper() and "COMPATIBLE" in r.summary()


def test_definition_identity_excludes_input_identity():
    a = {"observable": "x", "p": 1, "trajectory_view": {"cache_key": "v1"}}
    b = {"observable": "x", "p": 1, "trajectory_view": {"cache_key": "v2"}}
    assert definition_identity(a) == definition_identity(b)
    assert definition_identity(a) != definition_identity({**a, "p": 2})


@pytest.mark.parametrize("change,dimension", [
    (lambda ev: {**ev, "selections": {"a": 2}}, "definition"),
    (lambda ev: {**ev, "backend": {"tool": "gmx distance", "version": "2026.0"}}, "backend"),
])
def test_definition_and_backend_changes(tmp_path, change, dimension):
    cand, ev = _native(tmp_path)
    r = _eval(cand, change(ev))
    assert r.status == CS.INCOMPATIBLE
    assert next(c for c in r.checks if c.dimension == dimension).outcome == CheckOutcome.MISMATCH


def test_view_source_and_timeline_changes(tmp_path):
    cand, ev = _native(tmp_path)
    assert _eval(cand, ev, view="v2").status == CS.INCOMPATIBLE             # different view
    r = _eval(cand, ev, timeline=(0.0, 40.0, 80.0))                          # coarser sampling
    assert r.status == CS.INCOMPATIBLE and "sample 1" in r.reason
    assert _eval(cand, ev, timeline=(0.0, 20.0)).status == CS.INCOMPATIBLE    # no slicing
    dup, _ = _native(tmp_path, times=(0.0, 20.0, 20.0), name="dup.xvg")
    assert _eval(dup, ev, timeline=(0.0, 20.0, 40.0)).status == CS.INCOMPATIBLE
    assert _eval(dup, ev, timeline=(0.0, 20.0, 20.0)).status == CS.COMPATIBLE  # order kept
    unknown = evaluate(cand, requested_evidence=ev, requested_view_ref="v1",
                       requested_schema=_schema(), requested_timeline=None)
    assert unknown.status == CS.INSUFFICIENT_EVIDENCE                        # never guessed


def test_view_invariance_declaration_is_not_enough_without_source_proof(tmp_path):
    """v1 is conservative: a declared view invariance is recorded in the reason but
    a different view is still refused (no cross-view source equality proof yet)."""
    cand, ev = _native(tmp_path)
    r = evaluate(cand, requested_evidence=ev, requested_view_ref="v2", requested_schema=_schema(),
                 requested_timeline=[0.0, 20.0, 40.0], view_invariant=True)
    assert r.status == CS.INCOMPATIBLE and "view invariance" in r.reason


def test_modified_output_file_blocks_reuse(tmp_path):
    cand, ev = _native(tmp_path)
    Path(cand.source).write_text(Path(cand.source).read_text().replace("1.100", "9.100"))
    r = _eval(cand, ev)
    assert r.status == CS.INCOMPATIBLE and "changed after the analysis" in r.reason
    Path(cand.source).unlink()
    assert "missing" in _eval(cand, ev).reason


def test_external_and_name_only_candidates_never_auto(tmp_path):
    from analysis.campaign.results import import_external_series
    f = tmp_path / "rmsd.xvg"
    f.write_text('@    xaxis  label "Time (ps)"\n@    yaxis  label "RMSD (nm)"\n0 0.1\n20 0.2\n')
    name_only = Candidate(source=str(f), tier=ProvenanceTier.NAME_ONLY)
    r = evaluate(name_only, requested_evidence={"observable": "rmsd-receptor"},
                 requested_view_ref="v", requested_schema=[], requested_timeline=[0.0, 20.0])
    assert r.status == CS.INSUFFICIENT_EVIDENCE and "file name" in r.reason
    arr = import_external_series(f, quantity="rmsd", unit="nm", time_unit="ps")
    imp = Candidate(source=str(f), tier=ProvenanceTier.DECLARED,
                    result=AnalysisResult("rmsd-receptor", "s", "success", arrays=[arr]))
    r = evaluate(imp, requested_evidence={"observable": "rmsd-receptor"}, requested_view_ref="v",
                 requested_schema=[], requested_timeline=[0.0, 20.0])
    assert r.status == CS.EXPLICIT_IMPORT_ONLY


def test_output_resolution_formats(tmp_path):
    f = tmp_path / "m.xvg"
    f.write_text("0.000000e+00  1.723542e-01\n1.000000e+03  1.596403e-01\n")
    assert output_resolution(f, 1) == {"format": "exponent", "significant_digits": 7}
    g = tmp_path / "h.xvg"
    g.write_text("      0.000         11\n")
    assert output_resolution(g, 1)["decimals"] == 0


# ═══════════════════════════════════════════════════════════════════════════════
# Real GROMACS end-to-end (bundled GLP-1R membrane system)
# ═══════════════════════════════════════════════════════════════════════════════

SITE_YAML = """schema_version: simforge/annotations/v1
structural_annotation:
  residue_sets:
  - id: site
    kind: binding_site
    selection: {{numbering: topology, residues: "{res}", polymer_only: true}}
    provenance: {{origin: user_yaml}}
"""

PARAMS = {
    "rg": {"selection": "component:receptor"},
    "sasa": {"selection": "annotation:site", "surface": "component:receptor"},
    "com-distance": {"selection_a": "annotation:tm_1", "selection_b": "annotation:site"},
    "min-distance": {"selection_a": "component:receptor", "selection_b": "component:membrane"},
    "hbond-count": {"selection_a": "component:receptor", "selection_b": "component:membrane"},
}


def _study(root: Path, res="177-190", xtc=None) -> Path:
    """Phase 7 real layout (build-spec annotations tm_1..) + a user ``site``."""
    run = root / "run"
    (run / "metadata").mkdir(parents=True)
    spec = json.loads((REAL_GRO.parents[2] / "metadata" / "run_info.json").read_text())["yaml_source"]
    (run / "metadata" / "run_info.json").write_text(json.dumps({"yaml_source": spec}))
    step = run / "steps" / "11_production_md"
    step.mkdir(parents=True)
    (step / "md.xtc").symlink_to(xtc or REAL_XTC)
    for f in (REAL_TPR, REAL_GRO):
        (step / f.name).symlink_to(f)
    (run / "simforge_annotations.yaml").write_text(SITE_YAML.format(res=res))
    return run


def _guard(monkeypatch):
    import analysis.campaign.observables.generic as gen
    import analysis.campaign.observables.rmsd as rmsd
    def boom(*a, **k):
        raise AssertionError("a compatible reuse must not run GROMACS observables")
    monkeypatch.setattr(gen, "run_gmx", boom)
    monkeypatch.setattr(rmsd, "run_gmx", boom)


def _run(study, analyses=GENERIC + ["rmsd-receptor"], **kw):
    clear_view_cache()
    return {r.analysis_id: r for r in run_analyze(study, analyses, inspect_trajectories=True,
                                                   parameters=kw.pop("parameters", PARAMS),
                                                   **kw).results}


@requires_gmx
@requires_real_traj
def test_identical_request_reuses_every_observable(tmp_path, monkeypatch):
    study = _study(tmp_path)
    first = _run(study)
    assert all(r.status == AnalysisStatus.SUCCESS for r in first.values()), \
        {k: r.message for k, r in first.items()}
    prov_hashes = {k: sha(r.provenance_path) for k, r in first.items()}
    _guard(monkeypatch)
    second = _run(study)                                  # new invocation, disk state only
    for k, r in second.items():
        assert r.status == AnalysisStatus.CACHED and r.cached, (k, r.message)
        assert r.provenance_path == first[k].provenance_path
        assert [a.storage.path for a in r.arrays] == [a.storage.path for a in first[k].arrays]
        assert r.reused_from["definition_token"] == first[k].definition_token
        chosen = r.compatibility["chosen"]
        assert chosen["status"] == "compatible"
        assert {c["outcome"] for c in chosen["checks"]} == {"match"}
        assert sha(first[k].provenance_path) == prov_hashes[k]     # original untouched
        assert Path(first[k].provenance_path).with_name("reuse.json").is_file()
    reso = next(c for c in second["com-distance"].compatibility["chosen"]["checks"]
                if c["dimension"] == "output_resolution")["evidence"]["measured"]
    assert reso["com_distance"]["resolution"] == 0.001                # gmx distance printing


@requires_gmx
@requires_real_traj
def test_changed_annotation_invalidates_only_its_users(tmp_path):
    study = _study(tmp_path)
    _run(study)
    (study / "simforge_annotations.yaml").write_text(SITE_YAML.format(res="177-191"))
    second = _run(study)
    for k in ("sasa", "com-distance"):
        r = second[k]
        assert r.status == AnalysisStatus.SUCCESS and r.compatibility["decision"] == "recomputed"
        (ev,) = r.compatibility["evaluated"]
        assert ev["status"] == "incompatible" and "atoms_sha256" in ev["reason"]
    for k in ("rg", "min-distance", "hbond-count", "rmsd-receptor"):
        assert second[k].status == AnalysisStatus.CACHED, k


@requires_gmx
@requires_real_traj
@pytest.mark.parametrize("oid,change", [
    ("sasa", {"probe": 0.16}), ("sasa", {"ndots": 48}),
    ("hbond-count", {"distance_cutoff": 0.30}), ("hbond-count", {"angle_cutoff": 25}),
    ("min-distance", {"pbc": False}),
    ("min-distance", {"selection_b": "annotation:site"}),
])
def test_changed_parameters_block_reuse(tmp_path, oid, change):
    study = _study(tmp_path)
    _run(study, [oid])
    p = {**PARAMS, oid: {**PARAMS[oid], **change}}
    r = _run(study, [oid], parameters=p)[oid]
    assert r.status == AnalysisStatus.SUCCESS and r.compatibility["decision"] == "recomputed"
    assert r.compatibility["evaluated"][0]["status"] == "incompatible"


@requires_gmx
@requires_real_traj
def test_component_override_invalidates_ligand_like_selection(tmp_path, monkeypatch):
    """Changing the atoms behind a component group invalidates results using it;
    results that do not use it stay reusable."""
    study = _study(tmp_path)
    first = _run(study, ["min-distance", "rg"])
    # simulate an explicit component re-definition by editing the membrane group
    # atoms in the semantic index the next run builds
    import analysis.campaign.structure.index_builder as ib
    real = ib.build_semantic_index
    def narrower(**kw):
        idx = real(**kw)
        text = Path(idx.path).read_text().split("[ Membrane ]\n", 1)
        head, rest = text[0], text[1].split("\n[", 1)
        atoms = rest[0].split()
        Path(idx.path).write_text(head + "[ Membrane ]\n" + " ".join(atoms[:-1])
                                  + ("\n[" + rest[1] if len(rest) > 1 else "\n"))
        return idx
    monkeypatch.setattr("analysis.campaign.orchestration.study_analyzer.build_semantic_index",
                        narrower)
    second = _run(study, ["min-distance", "rg"])
    assert second["min-distance"].compatibility["decision"] == "recomputed"
    assert second["rg"].status == AnalysisStatus.CACHED


@requires_gmx
@requires_real_traj
def test_modified_output_and_filename_only_never_reused(tmp_path):
    study = _study(tmp_path)
    first = _run(study, ["rg", "com-distance"])
    xvg = Path(first["rg"].arrays[0].storage.path)
    xvg.write_text(xvg.read_text() + "\n")                       # tampered after analysis
    Path(first["com-distance"].provenance_path).unlink()         # only the file name remains
    second = _run(study, ["rg", "com-distance"])
    assert second["rg"].compatibility["evaluated"][0]["status"] == "incompatible"
    assert "changed after the analysis" in second["rg"].compatibility["evaluated"][0]["reason"]
    ev = second["com-distance"].compatibility["evaluated"][0]
    assert ev["status"] == "insufficient_evidence" and ev["candidate"]["tier"] == "tier3_name_only"
    assert second["rg"].status == second["com-distance"].status == AnalysisStatus.SUCCESS


@requires_gmx
@requires_real_traj
def test_inactive_annotation_is_not_resurrected_by_an_old_result(tmp_path):
    study = _study(tmp_path)
    _run(study, ["com-distance"])
    y = study / "simforge_annotations.yaml"
    y.write_text(y.read_text().replace("origin: user_yaml", "origin: derived"))
    r = _run(study, ["com-distance"])["com-distance"]
    assert r.status == AnalysisStatus.REVIEW_REQUIRED and "proposed" in r.message


@requires_gmx
@requires_real_traj
def test_other_location_same_content_reused_different_source_not(tmp_path, monkeypatch):
    a = _study(tmp_path / "a")
    _run(a, ["rg"])
    copy = tmp_path / "copy.xtc"
    shutil.copy(REAL_XTC, copy)                                   # same content, other path
    b = _study(tmp_path / "b", xtc=copy)
    r = _run(b, ["rg"], reuse_roots=[a / "simforge_analysis"])["rg"]
    assert r.status == AnalysisStatus.CACHED
    assert r.provenance_path.startswith(str(a))
    # a different trajectory (first 25 frames) is a different source/view/timeline
    short = tmp_path / "short.xtc"
    assert run_gmx(["trjconv", "-f", str(REAL_XTC), "-o", str(short), "-e", "480"],
                   stdin="0\n").ok
    c = _study(tmp_path / "c", xtc=short)
    r = _run(c, ["rg"], reuse_roots=[a / "simforge_analysis"])["rg"]
    assert r.status == AnalysisStatus.SUCCESS
    ev = r.compatibility["evaluated"][0]
    assert ev["status"] == "incompatible"
    assert {x["dimension"] for x in ev["checks"] if x["outcome"] == "mismatch"} >= {"view", "timeline"}


@requires_gmx
@requires_real_traj
def test_two_compatible_candidates_choose_current_location(tmp_path):
    a = _study(tmp_path / "a")
    _run(a, ["rg"])
    shutil.copytree(a / "simforge_analysis", tmp_path / "mirror", symlinks=True)
    r = _run(a, ["rg"], reuse_roots=[tmp_path / "mirror"])["rg"]
    assert r.status == AnalysisStatus.CACHED
    assert len([e for e in r.compatibility["evaluated"] if e["status"] == "compatible"]) == 2
    assert r.compatibility["chosen"]["candidate"]["in_current_location"] is True


@requires_gmx
@requires_real_traj
def test_explicit_import_is_labelled_not_equivalent(tmp_path):
    from analysis.campaign.results import import_external_series
    study = _study(tmp_path)
    ext = tmp_path / "rg_elsewhere.xvg"
    ext.write_text('@    xaxis  label "Time (ps)"\n@    yaxis  label "Radius (nm)"\n0 2.9\n20 3.0\n')
    arr = import_external_series(ext, quantity="radius_of_gyration", unit="nm", time_unit="ps")
    r = _run(study, ["rg"], use_imported={"rg": arr})["rg"]
    assert r.status == AnalysisStatus.CACHED and "not independently verified" in r.message
    assert r.compatibility["chosen"]["status"] == "explicit_import_only"
    assert r.reused_from["user_selected"] is True
    ext.write_text(ext.read_text() + "40 3.1\n")                    # file changed
    bad = _run(study, ["rg"], use_imported={"rg": arr})["rg"]
    assert bad.status == AnalysisStatus.REVIEW_REQUIRED
