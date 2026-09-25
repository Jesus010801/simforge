"""Phase 7 — generic, annotation-aware observables (rg, sasa, com-distance,
min-distance, hbond-count) on the existing registry / ResultArray / policy /
annotation machinery."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from analysis.campaign.annotations import (
    evaluate_com_to_com_axis, resolve_system_annotations,
)
from analysis.campaign.gmx import run_gmx
from analysis.campaign.models import (
    AnalysisStatus, FrameAlignmentMode, ObservablePurpose as P, SemanticIndex,
    SemanticIndexGroup, SystemRecord, TrajectoryRequirements as TR, TrajectoryView,
)
from analysis.campaign.observables import registry
from analysis.campaign.observables.base import AnalysisContext
from analysis.campaign.observables.selections import SelectionRef
from analysis.campaign.structure.component_detector import detect_components
from core.structural_annotation import AxisAnnotation, AxisDefinition
from tests.analysis.campaign.conftest import (
    REAL_GRO, REAL_TPR, REAL_XTC, requires_gmx, requires_real_traj,
)
from tests.analysis.campaign.test_annotations import aset, atoms_of_dimer, rs, sel, write_gro

GENERIC = ("rg", "sasa", "com-distance", "min-distance", "hbond-count")


# ═══════════════════════════════════════════════════════════════════════════════
# Synthetic context (no gmx execution): homodimer (chains A, B) + LIG
# ═══════════════════════════════════════════════════════════════════════════════

def make_ctx(tmp_path: Path, *decls, params=None, ligand_extra=(), sub="w",
             topology=True) -> AnalysisContext:
    d = tmp_path / sub
    d.mkdir(exist_ok=True)
    gro = write_gro(d / "ref.gro")
    rec = SystemRecord("s", "c", "r", structure_path=str(gro))
    rec.components = detect_components(structure_path=str(gro)).components
    rec.annotations, _ = resolve_system_annotations(rec, extra=[aset(*decls)])
    atoms = atoms_of_dimer()
    groups = {
        "System": list(range(1, len(atoms) + 1)),
        "Receptor": [i for i, a in enumerate(atoms, 1) if a[2] != "LIG"],
        "Ligand": [i for i, a in enumerate(atoms, 1) if a[2] == "LIG"] + list(ligand_extra),
        "Membrane": [1],
    }
    for r in rec.annotations:
        if r.group_name:
            groups[r.group_name] = r.atom_ids
    ndx = d / "semantic.ndx"
    ndx.write_text("".join(f"[ {n} ]\n{' '.join(map(str, a))}\n" for n, a in groups.items()))
    tpr = d / "md.tpr"
    tpr.write_bytes(b"stub-tpr")
    view = TrajectoryView(kind="raw", requirements=TR(), path=str(d / "md.xtc"), cache_key="view-1")
    return AnalysisContext(
        system=rec, semantic_index=SemanticIndex(path=str(ndx), groups=[
            SemanticIndexGroup(n, len(a)) for n, a in groups.items()]),
        trajectory_view=view, topology_path=str(tpr) if topology else str(gro),
        structure_path=str(gro), output_dir=d / "out", parameters=params or {},
        gmx="/nonexistent/gmx")


SITE = rs("site", sel(chain="A", residues="10-12"))
SUB_A = rs("subunit", sel(chain="A"), kind="selected_subunit")
SUB_B = rs("subunit", sel(chain="B"), kind="selected_subunit")
OTHER = rs("other", sel(chain="B", residues="1-3"), kind="custom")


def token(obs_id, ctx):
    return registry.get(obs_id).definition_token(ctx)


# ═══════════════════════════════════════════════════════════════════════════════
# Registry / contract
# ═══════════════════════════════════════════════════════════════════════════════

def test_generic_observables_are_ordinary_registry_entries():
    expected = {"rg": ("radius_of_gyration", "nm", P.INTRAMOLECULAR_SHAPE),
                "sasa": ("sasa", "nm^2", P.INTRAMOLECULAR_SHAPE),
                "com-distance": ("com_distance", "nm", P.INTER_COMPONENT_GEOMETRY),
                "min-distance": ("minimum_distance", "nm", P.INTER_COMPONENT_GEOMETRY),
                "hbond-count": ("hydrogen_bond_count", "count", P.INTER_COMPONENT_GEOMETRY)}
    for oid, (q, unit, purpose) in expected.items():
        spec = registry.get(oid)
        (arr,) = spec.output_schema()
        assert (arr.quantity, arr.unit, arr.axis_signature()) == (q, unit, ("time",))
        assert spec.purpose == purpose and spec.dependencies() == []
        assert spec.to_dict()["parameters_schema"]


def test_requirements_never_request_reconstruction():
    for oid in ("com-distance", "min-distance", "hbond-count"):
        req = registry.get(oid).trajectory_requirements({})
        assert req.minimum_image_distances and req.view_kind() == "raw"
        assert not (req.requires_whole_molecules or req.centering_target or req.fit_selection
                    or req.reconstruction)
    for oid in ("rg", "sasa"):
        assert registry.get(oid).trajectory_requirements({}).view_kind() == "raw"
    nopbc = registry.get("min-distance").trajectory_requirements({"min-distance": {"pbc": False}})
    assert not nopbc.minimum_image_distances


def test_selection_ref_parsing():
    assert str(SelectionRef.parse("component:ligand")) == "component:ligand"
    assert SelectionRef.parse("annotation:site").kind == "annotation"
    for bad in ("ligand", "component:banana", "group:Receptor", "annotation:"):
        with pytest.raises(ValueError):
            SelectionRef.parse(bad)


# ═══════════════════════════════════════════════════════════════════════════════
# Definition identity
# ═══════════════════════════════════════════════════════════════════════════════

PAIR = {"com-distance": {"selection_a": "component:ligand", "selection_b": "annotation:site"}}


def test_annotation_changes_only_its_users(tmp_path):
    base = make_ctx(tmp_path, SITE, OTHER, params=PAIR, sub="base")
    moved = make_ctx(tmp_path, rs("site", sel(chain="A", residues="10-13")), OTHER,
                     params=PAIR, sub="moved")
    unrelated = make_ctx(tmp_path, SITE, rs("other", sel(chain="B", residues="1-4"), kind="custom"),
                         params=PAIR, sub="unrelated")
    described = make_ctx(tmp_path, rs("site", sel(chain="A", residues="10-12"),
                                      description="edited text"), OTHER, params=PAIR, sub="desc")
    t = token("com-distance", base)
    assert t and token("com-distance", moved) != t                      # used annotation changed
    assert token("com-distance", unrelated) == t                         # unrelated annotation
    assert token("com-distance", described) == t                         # description only
    base.output_dir = tmp_path / "elsewhere"
    assert token("com-distance", base) == t                              # output directory
    # a real, unrelated observable is untouched by the site change
    assert token("rg", base) == token("rg", moved)


def test_component_atoms_are_part_of_identity(tmp_path):
    a = make_ctx(tmp_path, SITE, params=PAIR, sub="a")
    b = make_ctx(tmp_path, SITE, params=PAIR, sub="b", ligand_extra=(1,))
    assert token("com-distance", a) != token("com-distance", b)
    ev = registry.get("com-distance").definition_evidence(a)
    assert ev["selections"]["selection_a"]["atoms"]["n_atoms"] == 10
    assert ev["selections"]["selection_b"]["annotation"]["id"] == "site"
    assert ev["coordinate_semantics"]["pbc"].startswith("minimum-image")


def test_subunit_choice_changes_identity(tmp_path):
    p = {"rg": {"selection": "annotation:subunit"}}
    ta = token("rg", make_ctx(tmp_path, SUB_A, params=p, sub="A"))
    tb = token("rg", make_ctx(tmp_path, SUB_B, params=p, sub="B"))
    assert ta and tb and ta != tb


def test_metric_pbc_and_cutoffs_are_part_of_identity(tmp_path):
    sel_pair = {"selection_a": "component:ligand", "selection_b": "annotation:site"}
    ctx = make_ctx(tmp_path, SITE, params={"com-distance": sel_pair, "min-distance": sel_pair,
                                           "hbond-count": sel_pair}, sub="m")
    assert token("com-distance", ctx) != token("min-distance", ctx)        # metric
    t_pbc = token("min-distance", ctx)
    ctx.parameters["min-distance"] = {**sel_pair, "pbc": False}
    assert token("min-distance", ctx) != t_pbc
    t_hb = token("hbond-count", ctx)
    assert token("hbond-count", ctx) == t_hb
    ctx.parameters["hbond-count"] = {**sel_pair, "distance_cutoff": 0.30}
    assert token("hbond-count", ctx) != t_hb
    ctx.parameters["hbond-count"] = {**sel_pair, "angle_cutoff": 25}
    assert token("hbond-count", ctx) != t_hb


def test_rmsd_token_ignores_other_observables_parameter_blocks(tmp_path):
    ctx = make_ctx(tmp_path, SITE, sub="r")
    from analysis.campaign.models import SemanticIndexGroup as G
    ev1 = registry.get("rmsd-receptor").definition_evidence(ctx)
    ctx.parameters = {"com-distance": {"selection_a": "component:ligand"}}
    ev2 = registry.get("rmsd-receptor").definition_evidence(ctx)
    assert ev1 == ev2


# ═══════════════════════════════════════════════════════════════════════════════
# Blocking: inactive / unresolved annotations, missing .tpr, applicability
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize("decl,needle", [
    (rs("site", sel(chain="A", residues="10-12"), origin="derived"), "proposed"),
    (rs("site", sel("author_pdb", chain="A", residues="10")), "unsupported"),
    (rs("site", sel(chain="A", residues="10-40")), "not found"),
    (rs("site", sel(residues="10")), "review_required"),
])
def test_inactive_annotations_block(tmp_path, decl, needle):
    ctx = make_ctx(tmp_path, decl, params=PAIR)
    r = registry.get("com-distance").execute(ctx)
    assert r.status == AnalysisStatus.REVIEW_REQUIRED and r.arrays == []
    assert "site" in r.message and needle in r.message
    assert registry.get("com-distance").definition_token(ctx) is None


def test_missing_tpr_and_missing_params(tmp_path):
    ctx = make_ctx(tmp_path, SITE, params=PAIR, topology=False)
    assert registry.get("com-distance").execute(ctx).status == AnalysisStatus.SKIPPED
    bare = make_ctx(tmp_path, SITE, sub="bare")
    r = registry.get("com-distance").execute(bare)
    assert r.status == AnalysisStatus.REVIEW_REQUIRED and "selection_a" in r.message


def test_applicability_without_execution(tmp_path):
    ctx = make_ctx(tmp_path, rs("site", sel(chain="A", residues="10"), origin="derived"))
    ok = registry.get("rg").applicability(ctx.system, {}, ctx.semantic_index)
    assert ok["applicable"]                                                   # receptor default
    blocked = registry.get("com-distance").applicability(
        ctx.system, PAIR, ctx.semantic_index)
    assert not blocked["applicable"] and blocked["missing_annotations"] == ["site"]
    missing = registry.get("com-distance").applicability(ctx.system, {}, ctx.semantic_index)
    assert not missing["applicable"] and "required" in missing["reasons"][0]


def test_dry_run_command(tmp_path):
    ctx = make_ctx(tmp_path, SITE, params=PAIR)
    ctx.dry_run = True
    r = registry.get("com-distance").execute(ctx)
    assert r.status == AnalysisStatus.PLANNED
    cmd = r.data_summary["command"]
    assert cmd[:2] == ["gmx", "distance"]
    assert 'com of group "Ligand" plus com of group "Ann_site"' in cmd


# ═══════════════════════════════════════════════════════════════════════════════
# Axis evaluation and explicit import
# ═══════════════════════════════════════════════════════════════════════════════

def test_com_to_com_axis_evaluation(tmp_path):
    ctx = make_ctx(tmp_path, rs("up", sel(chain="A", residues="1-2"), kind="gate"),
                   rs("down", sel(chain="A", residues="29-30"), kind="gate"),
                   AxisAnnotation(id="pore_axis", kind="pore_axis", definition=AxisDefinition(
                       type="com_to_com", from_annotation="up", to_annotation="down")))
    axis = next(a for a in ctx.system.annotations if a.annotation_id == "pore_axis")
    up = np.array([[0, 0, 0], [0, 0, 1], [1, 1, 1]], float)
    down = np.array([[0, 0, 2], [0, 0, 1], [1, 1, 4]], float)          # frame 1 degenerate
    vec, length, unit = evaluate_com_to_com_axis(axis, {"up": up, "down": down})
    assert length.tolist() == [2.0, 0.0, 3.0]
    assert unit[0].tolist() == [0, 0, 1] and np.isnan(unit[1]).all()    # not normalised away
    assert unit[2].tolist() == [0, 0, 1]
    site = next(a for a in ctx.system.annotations if a.annotation_id == "up")
    with pytest.raises(ValueError):
        evaluate_com_to_com_axis(site, {"up": up, "down": down})


def test_import_external_series(tmp_path):
    from analysis.campaign.results import import_external_series, read_column
    f = tmp_path / "custom_rmsd.xvg"
    f.write_text('@    xaxis  label "Time (ps)"\n@    yaxis  label "Distance (nm)"\n0 1.0\n20 1.5\n')
    arr = import_external_series(f, quantity="custom_distance", unit="nm", time_unit="ps")
    assert arr.quantity == "custom_distance" and arr.attrs["externally_supplied"]
    assert arr.storage.fingerprint.digest and read_column(arr.storage) == [1.0, 1.5]
    with pytest.raises(ValueError):
        import_external_series(f, quantity="x", unit="kJ/mol", time_unit="ps")   # contradicts header
    with pytest.raises(ValueError):
        import_external_series(f, quantity="x", unit="nm", time_unit="ns")
    c = tmp_path / "series.csv"
    c.write_text("time,value\n0,3\n10,4\n")
    assert import_external_series(c, quantity="q", unit="nm", time_unit="ps").unit == "nm"


# ═══════════════════════════════════════════════════════════════════════════════
# Real GROMACS: bundled GLP-1R membrane system with its real build-spec annotations
# ═══════════════════════════════════════════════════════════════════════════════

def _real_run(tmp_path) -> Path:
    run = tmp_path / "run"
    (run / "metadata").mkdir(parents=True)
    spec = json.loads((REAL_GRO.parents[2] / "metadata" / "run_info.json").read_text())["yaml_source"]
    (run / "metadata" / "run_info.json").write_text(json.dumps({"yaml_source": spec}))
    step = run / "steps" / "11_production_md"
    step.mkdir(parents=True)
    for f in (REAL_XTC, REAL_TPR, REAL_GRO):
        (step / f.name).symlink_to(f)
    return run


def _data_rows(p) -> list[str]:
    return [l for l in Path(p).read_text().splitlines() if l and l[0] not in "#@"]


@requires_gmx
@requires_real_traj
def test_real_generic_observables(tmp_path, monkeypatch):
    import analysis.campaign.diagnostics as dmod
    from analysis.campaign.orchestration.study_analyzer import run_analyze
    calls = []
    real = dmod.diagnose_system
    monkeypatch.setattr(dmod, "diagnose_system", lambda rec, **kw: calls.append(1) or real(rec, **kw))
    run = _real_run(tmp_path)
    params = {
        "rg": {"selection": "component:receptor"},
        "sasa": {"selection": "annotation:tm_1", "surface": "component:receptor"},
        "com-distance": {"selection_a": "annotation:tm_1", "selection_b": "annotation:tm_2"},
        "min-distance": {"selection_a": "component:receptor", "selection_b": "component:membrane"},
        "hbond-count": {"selection_a": "component:receptor", "selection_b": "component:membrane"},
    }
    res = run_analyze(run, list(GENERIC) + ["rmsd-receptor"], inspect_trajectories=True,
                      parameters=params)
    by = {r.analysis_id: r for r in res.results}
    assert len(calls) == 1                                           # diagnostics once
    expected = {"rg": ("radius_of_gyration", "nm"), "sasa": ("sasa", "nm^2"),
                "com-distance": ("com_distance", "nm"),
                "min-distance": ("minimum_distance", "nm"),
                "hbond-count": ("hydrogen_bond_count", "count")}
    for oid, (q, unit) in expected.items():
        r = by[oid]
        assert r.status == AnalysisStatus.SUCCESS, (oid, r.message)
        (arr,) = r.arrays
        prov = json.loads(Path(r.provenance_path).read_text())
        assert (arr.quantity, arr.unit) == (q, unit)
        t = arr.axis("time")
        assert t.unit == "ps" and t.alignment.mode == FrameAlignmentMode.ONE_ROW_PER_FRAME
        assert t.attrs["rows"] == t.attrs["frames"] == 51
        assert arr.view_ref == prov["trajectory_view"]["cache_key"]
        assert arr.definition_token == r.definition_token == prov["result"]["definition_token"]
        assert Path(arr.storage.path).is_file() and arr.storage.fingerprint
        assert prov["trajectory_view"]["kind"] == "raw"                  # no derived view
    assert by["hbond-count"].arrays[0].attrs["unit_source"].startswith("fallback")
    assert by["sasa"].arrays[0].attrs["unit_source"] == "xvg header"
    # inter-component observables were admitted as minimum-image analyses
    prov = json.loads(Path(by["com-distance"].provenance_path).read_text())
    (dec,) = prov["trajectory_view"]["decisions"]
    assert (dec["operation"], dec["classification"]) == ("minimum_image", "valid")
    # generic Rg in a membrane system requested no fit / centring at all
    assert json.loads(Path(by["rg"].provenance_path).read_text())["trajectory_view"]["decisions"] == []
    # (the whole view on disk belongs to rmsd-receptor, as before; generic ones used raw)
    assert not list(run.rglob("fitted.xtc")) and not list(run.rglob("centered.xtc"))
    # the recorded command reproduces the numbers independently
    for oid in ("rg", "com-distance"):
        r = by[oid]
        argv = r.data_summary["command"][1:]
        out = tmp_path / f"indep_{oid}.xvg"
        flag = "-o" if "-o" in argv else "-oall"
        argv[argv.index(flag) + 1] = str(out)
        assert run_gmx(argv, stdin=r.data_summary["stdin"]).ok
        assert _data_rows(out) == _data_rows(r.arrays[0].storage.path)
    assert by["rmsd-receptor"].status == AnalysisStatus.SUCCESS       # RMSD coexists


@requires_gmx
@requires_real_traj
def test_real_unresolved_annotations_block(tmp_path):
    from analysis.campaign.orchestration.study_analyzer import run_analyze
    run = _real_run(tmp_path)
    res = run_analyze(run, ["rg", "com-distance"], inspect_trajectories=False, parameters={
        "rg": {"selection": "annotation:tm_6"},
        "com-distance": {"selection_a": "annotation:intracellular_1",
                         "selection_b": "annotation:tm_1"}})
    by = {r.analysis_id: r for r in res.results}
    rg, dist = by["rg"], by["com-distance"]
    assert rg.status == AnalysisStatus.REVIEW_REQUIRED and "tm_6" in rg.message
    assert "unresolved" in rg.message and "401, 402, 403, 404" in rg.message
    assert dist.status == AnalysisStatus.REVIEW_REQUIRED and "intracellular_1" in dist.message
    assert rg.arrays == [] and dist.arrays == []


def test_gmx_failure_message_is_the_gromacs_explanation():
    from analysis.campaign.gmx import gmx_error_summary
    stderr = ("Reading frame 0 time 0.000\n"
              "-------------------------------------------------------\n"
              "Program:     gmx hbond, version 2025.2\n"
              "Source file: src/gromacs/trajectoryanalysis/modules/hbond.cpp (line 612)\n"
              "Function:    void gmx::analysismodules::(anonymous namespace)::Hbond::linkDA("
              "gmx::analysismodules::(anonymous namespace)::t_info *)\n\n"
              "Inconsistency in user input:\n"
              "Selection Ligand' has no donors AND has no acceptors! Nothing to be done.\n\n"
              "For more information and tips for troubleshooting, please check the GROMACS\n"
              "website at https://manual.gromacs.org/current/user-guide/run-time-errors.html\n"
              "-------------------------------------------------------\n")
    assert gmx_error_summary(stderr) == ("Inconsistency in user input: Selection Ligand' has no "
                                         "donors AND has no acceptors! Nothing to be done.")
    assert gmx_error_summary("plain failure text") == "plain failure text"
