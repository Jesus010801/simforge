"""Labelled N-dimensional result contract (ResultArray / Axis / StorageRef) and
the RMSD migration onto it."""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from analysis.campaign.models import (
    AnalysisResult, AnalysisStatus, Axis, AxisKind, FrameAlignment, FrameAlignmentMode,
    MissingSemantics, ObservablePurpose, ResultArray, SemanticIndex, SemanticIndexGroup,
    StorageFormat, StorageRef, SystemRecord, TrajectoryRequirements, TrajectoryView, ViewKind,
)
from analysis.campaign.observables import registry
from analysis.campaign.observables.base import AnalysisContext, ObservableSpec
from analysis.campaign.results import (
    atom_set_hash, definition_token, index_group_evidence, read_column, xvg_units,
)
from runtime.xvg_parser import parse_xvg
from tests.analysis.campaign.conftest import (
    REAL_GRO, REAL_TPR, REAL_XTC, requires_gmx, requires_real_traj,
)


def _roundtrip(arr: ResultArray) -> ResultArray:
    back = ResultArray.from_dict(json.loads(json.dumps(arr.to_dict())))
    assert back.to_dict() == arr.to_dict()
    return back


# ═══════════════════════════════════════════════════════════════════════════════
# Contract shapes (architecture tests — no observable is implemented here)
# ═══════════════════════════════════════════════════════════════════════════════

def test_rmsd_time_series_roundtrip():
    xvg = StorageRef(path="/data/rmsd.xvg", format=StorageFormat.XVG, column=1)
    arr = ResultArray(
        name="rmsd", quantity="rmsd", unit="nm",
        axes=[Axis("time", AxisKind.TIME, unit="ps",
                   values_ref=StorageRef("/data/rmsd.xvg", StorageFormat.XVG, column=0),
                   alignment=FrameAlignment(FrameAlignmentMode.ONE_ROW_PER_FRAME,
                                            trajectory_digest="abc", view_ref="k1"))],
        storage=xvg, view_ref="k1", definition_token="tok",
        missing=MissingSemantics.NONE_EXPECTED,
    )
    back = _roundtrip(arr)
    assert back.axis_signature() == ("time",) and back.realised and back.validate() == []
    assert back.axes[0].alignment.mode == FrameAlignmentMode.ONE_ROW_PER_FRAME


def test_apl_time_by_leaflet():
    arr = ResultArray(
        name="apl", quantity="area_per_lipid", unit="nm^2",
        axes=[Axis("time", AxisKind.TIME, unit="ps"),
              Axis("leaflet", AxisKind.CATEGORY, values=["upper", "lower"])],
        storage=StorageRef("/d/apl.npy", StorageFormat.NPY),
    )
    back = _roundtrip(arr)
    assert back.axis_signature() == ("time", "category")
    assert back.axis("leaflet").values == ["upper", "lower"] and back.validate() == []


def test_rmsf_static_residue_profile():
    arr = ResultArray(
        name="rmsf", quantity="rmsf", unit="nm",
        axes=[Axis("residue", AxisKind.CATEGORY, values=["A:45", "A:46", "A:47"])],
        attrs={"window_ps": [0.0, 1000.0]},
    )
    back = _roundtrip(arr)
    assert back.axis_signature() == ("category",) and not back.realised
    assert back.attrs["window_ps"] == [0.0, 1000.0]


def test_pore_radius_time_by_z_with_frame():
    z = Axis("z", AxisKind.COORDINATE, unit="nm", frame="pore_axis:test",
             values_ref=StorageRef("/d/pore.npz", StorageFormat.NPZ, member="z"))
    arr = ResultArray(
        name="pore_radius", quantity="pore_radius", unit="nm",
        axes=[Axis("time", AxisKind.TIME, unit="ps"), z],
        storage=StorageRef("/d/pore.npz", StorageFormat.NPZ, member="radius"),
        missing=MissingSemantics.UNDEFINED_STATE,
    )
    back = _roundtrip(arr)
    assert back.axis("z").frame == "pore_axis:test"          # annotation need not exist
    assert back.axis_signature() == ("time", "coordinate")
    assert back.missing == MissingSemantics.UNDEFINED_STATE and back.validate() == []


def test_contact_matrices():
    res = ["A:1", "A:2", "B:7"]
    static = ResultArray(
        name="contacts", quantity="contact_frequency", unit=None,
        axes=[Axis("residue_i", AxisKind.CATEGORY, values=res),
              Axis("residue_j", AxisKind.CATEGORY, values=res)])
    timed = ResultArray(
        name="contacts_t", quantity="contact", unit=None,
        axes=[Axis("time", AxisKind.TIME, unit="ps"),
              Axis("residue_i", AxisKind.CATEGORY, values=res),
              Axis("residue_j", AxisKind.CATEGORY, values=res)])
    assert _roundtrip(static).axis_signature() == ("category", "category")
    assert _roundtrip(timed).axis_signature() == ("time", "category", "category")


def test_membrane_thickness_map_and_order_parameter():
    thick = ResultArray(name="thickness", quantity="membrane_thickness", unit="nm",
                        axes=[Axis("time", AxisKind.TIME, unit="ps"),
                              Axis("x", AxisKind.COORDINATE, unit="nm", frame="box"),
                              Axis("y", AxisKind.COORDINATE, unit="nm", frame="box")])
    order = ResultArray(name="scd", quantity="order_parameter", unit=None,
                        axes=[Axis("time", AxisKind.TIME, unit="ps"),
                              Axis("carbon", AxisKind.INDEX, values=list(range(2, 17)))])
    assert _roundtrip(thick).axis_signature() == ("time", "coordinate", "coordinate")
    assert _roundtrip(order).axis_signature() == ("time", "index")


def test_unknown_axis_kind_is_still_valid():
    arr = ResultArray(name="x", quantity="custom", axes=[Axis("q", "reciprocal_space")])
    assert _roundtrip(arr).validate() == []


def test_duplicate_time_values_never_deduplicated():
    ax = Axis("time", AxisKind.TIME, unit="ps", values=[100.0, 100.0, 100.5])
    back = Axis.from_dict(json.loads(json.dumps(ax.to_dict())))
    assert back.values == [100.0, 100.0, 100.5] and back.validate() == []


def test_validation_problems():
    assert ResultArray(name="", quantity="").validate()
    dup = ResultArray(name="a", quantity="q", axes=[Axis("t", "time"), Axis("t", "time")])
    assert any("duplicate axis" in p for p in dup.validate())
    bad_align = Axis("z", AxisKind.COORDINATE, alignment=FrameAlignment())
    assert any("only applies to time" in p for p in bad_align.validate())
    both = Axis("t", "time", values=[1.0], values_ref=StorageRef("f.xvg", "xvg", column=0))
    assert any("both inline" in p for p in both.validate())


def test_storage_ref_validation():
    assert StorageRef("f.xvg", StorageFormat.XVG, column=1).validate() == []
    assert StorageRef("f.csv", StorageFormat.CSV, column=0).validate() == []
    assert StorageRef("f.npy", StorageFormat.NPY).validate() == []
    assert StorageRef("f.npz", StorageFormat.NPZ, member="x").validate() == []
    assert StorageRef("f.xvg", StorageFormat.XVG).validate()               # needs column
    assert StorageRef("f.npz", StorageFormat.NPZ).validate()               # needs member
    assert StorageRef("f.h5", "hdf5", column=0).validate()                 # unknown format
    assert StorageRef("f.csv", StorageFormat.CSV, column=-1).validate()
    with pytest.raises(NotImplementedError):
        read_column(StorageRef("f.npz", StorageFormat.NPZ, member="x"))
    with pytest.raises(ValueError):
        read_column(StorageRef("f.xvg", StorageFormat.XVG))


_XVG = """# comment
@    title "RMSD"
@    xaxis  label "Time (ps)"
@    yaxis  label "RMSD (nm)"
@TYPE xy
@ subtitle "Backbone after lsq fit to Backbone"
  100.0000000    0.1000000
  100.0000000    0.2000000
  100.5000000    0.3000000
"""


def test_xvg_storage_reads_raw_columns_in_order(tmp_path):
    p = tmp_path / "rmsd.xvg"
    p.write_text(_XVG)
    assert read_column(StorageRef(str(p), StorageFormat.XVG, column=0)) == [100.0, 100.0, 100.5]
    assert read_column(StorageRef(str(p), StorageFormat.XVG, column=1)) == [0.1, 0.2, 0.3]
    with pytest.raises(ValueError):
        read_column(StorageRef(str(p), StorageFormat.XVG, column=2))
    assert xvg_units(p) == ("ps", "nm")


def test_csv_storage(tmp_path):
    p = tmp_path / "apl.csv"
    p.write_text("# produced by test\ntime_ps,upper,lower\n0,0.61,0.62\n20,0.63,0.60\n")
    assert read_column(StorageRef(str(p), StorageFormat.CSV, column=0)) == [0.0, 20.0]
    assert read_column(StorageRef(str(p), StorageFormat.CSV, column=2)) == [0.62, 0.60]
    with pytest.raises(ValueError):
        read_column(StorageRef(str(p), StorageFormat.CSV, column=5))


# ═══════════════════════════════════════════════════════════════════════════════
# AnalysisResult compatibility
# ═══════════════════════════════════════════════════════════════════════════════

def test_legacy_analysis_result_deserialises():
    legacy = {"analysis_id": "rmsd-receptor", "system_id": "s1", "status": "success",
              "output_files": ["/x/rmsd.xvg"], "data_summary": {"mean_nm": 0.2}}
    r = AnalysisResult.from_dict(legacy)
    assert r.arrays == [] and r.definition_token is None
    assert r.output_files == ["/x/rmsd.xvg"] and r.data_summary["mean_nm"] == 0.2
    assert r.trajectory_view_kind == ViewKind.RAW


def test_analysis_result_roundtrip_with_arrays():
    r = AnalysisResult(analysis_id="a", system_id="s", status=AnalysisStatus.SUCCESS,
                       arrays=[ResultArray(name="rmsd", quantity="rmsd", unit="nm",
                                           axes=[Axis("time", AxisKind.TIME, unit="ps")])],
                       definition_token="t", definition_evidence={"k": 1})
    d = json.loads(json.dumps(r.to_dict()))
    back = AnalysisResult.from_dict(d)
    assert back.to_dict() == r.to_dict()
    # legacy keys still present and unchanged in shape
    assert {"output_files", "data_summary", "status"} <= set(d)


# ═══════════════════════════════════════════════════════════════════════════════
# ObservableSpec extensions
# ═══════════════════════════════════════════════════════════════════════════════

def test_observable_spec_defaults_are_inert():
    class Plain(ObservableSpec):
        id = "plain"
    p = Plain()
    assert p.purpose is None and p.output_schema() == [] and p.dependencies() == []
    assert p.definition_evidence(None) is None and p.definition_token(None) is None
    d = p.to_dict()
    assert d["purpose"] is None and d["output_schema"] == []


def test_rmsd_specs_declare_output_schema_and_purpose():
    for spec in registry.all_specs():
        (arr,) = spec.output_schema()
        assert arr.quantity == "rmsd" and arr.unit == "nm" and not arr.realised
        assert arr.axis_signature() == ("time",) and arr.axes[0].unit == "ps"
        assert spec.dependencies() == []
        entry = spec.to_dict()
        assert entry["output_schema"][0]["quantity"] == "rmsd"
    assert registry.get("rmsd-receptor").purpose == ObservablePurpose.INTRAMOLECULAR_SHAPE
    assert (registry.get("rmsd-peptide-receptor-frame").purpose
            == ObservablePurpose.INTER_COMPONENT_GEOMETRY)


# ═══════════════════════════════════════════════════════════════════════════════
# RMSD definition token (no gmx needed)
# ═══════════════════════════════════════════════════════════════════════════════

_NDX = """[ Receptor ]
1 2 3 4 5 6
[ Receptor_Backbone ]
1 2 3
[ Peptide_Backbone ]
7 8
"""


def _ctx(tmp_path, *, ndx=_NDX, ref_bytes=b"TPRDATA", cache_key="view-1",
         out="out1", params=None, dry_run=True, view_path=None):
    ndx_p = tmp_path / "semantic.ndx"
    ndx_p.write_text(ndx)
    tpr = tmp_path / "md.tpr"
    tpr.write_bytes(ref_bytes)
    idx = SemanticIndex(path=str(ndx_p), groups=[
        SemanticIndexGroup(name=n, n_atoms=0) for n in
        ("Receptor", "Receptor_Backbone", "Peptide_Backbone")])
    view = TrajectoryView(kind=ViewKind.WHOLE, requirements=TrajectoryRequirements(
        requires_whole_molecules=True), path=view_path or str(tmp_path / "whole.xtc"),
        cache_key=cache_key)
    return AnalysisContext(
        system=SystemRecord(system_id="s1", condition_id="c", replicate_id="r1"),
        semantic_index=idx, trajectory_view=view, topology_path=str(tpr),
        structure_path=None, output_dir=tmp_path / out, parameters=params or {},
        gmx="/nonexistent/gmx", dry_run=dry_run)


def test_rmsd_token_deterministic_and_path_independent(tmp_path):
    spec = registry.get("rmsd-receptor")
    t1 = spec.definition_token(_ctx(tmp_path))
    t2 = spec.definition_token(_ctx(tmp_path))
    t3 = spec.definition_token(_ctx(tmp_path, out="somewhere/else"))
    assert t1 and t1 == t2 == t3
    ev = spec.definition_evidence(_ctx(tmp_path))
    assert str(tmp_path) not in json.dumps(ev)          # no machine-local paths


@pytest.mark.parametrize("change", [
    {"ndx": _NDX.replace("[ Receptor_Backbone ]\n1 2 3", "[ Receptor_Backbone ]\n1 2 4")},
    {"ref_bytes": b"OTHER-REFERENCE"},
    {"cache_key": "view-2"},
    {"params": {"b": 100}},
])
def test_rmsd_token_changes_with_definition(tmp_path, change):
    spec = registry.get("rmsd-receptor")
    base = spec.definition_token(_ctx(tmp_path))
    other = tmp_path / "other"
    other.mkdir()
    changed = spec.definition_token(_ctx(other, **change))
    same = spec.definition_token(_ctx(other))
    assert base == same != changed


def test_rmsd_token_differs_between_rmsd_definitions(tmp_path):
    ctx = _ctx(tmp_path)
    toks = {registry.get(a).definition_token(ctx)
            for a in ("rmsd-receptor", "rmsd-peptide-intrinsic", "rmsd-peptide-receptor-frame")}
    assert len(toks) == 3


def test_rmsd_token_absent_without_evidence(tmp_path):
    ctx = _ctx(tmp_path, ndx="[ Receptor ]\n1 2\n")
    assert registry.get("rmsd-receptor").definition_token(ctx) is None


def test_group_evidence_is_content_based(tmp_path):
    p = tmp_path / "i.ndx"
    p.write_text("[ G ]\n3 1 2 2\n")
    ev = index_group_evidence(p, "G")
    # repeated atom 2 => multiplicity evidence is added (Phase 3); set hash unchanged
    assert {k: ev[k] for k in ("name", "n_atoms", "atoms_sha256")} == {
        "name": "G", "n_atoms": 3, "atoms_sha256": atom_set_hash([1, 2, 3])}
    assert ev["n_entries"] == 4 and ev["multiset_sha256"]
    p.write_text("[ G ]\n3 1 2\n")
    assert set(index_group_evidence(p, "G")) == {"name", "n_atoms", "atoms_sha256"}
    assert index_group_evidence(p, "missing") is None
    assert definition_token({"a": 1, "b": 2}) == definition_token({"b": 2, "a": 1})


def test_rmsd_dry_run_command_unchanged(tmp_path):
    ctx = _ctx(tmp_path)
    r = registry.get("rmsd-receptor").execute(ctx)
    assert r.status == AnalysisStatus.PLANNED
    assert r.data_summary["planned_command"] == [
        "gmx", "rms", "-s", ctx.topology_path, "-f", ctx.trajectory_view.path,
        "-n", ctx.semantic_index.path, "-o", str(ctx.output_dir / "rmsd.xvg"), "-tu", "ps"]
    assert r.data_summary["stdin"] == "Receptor_Backbone\nReceptor_Backbone\n"
    frame = registry.get("rmsd-peptide-receptor-frame").execute(ctx)
    assert frame.data_summary["planned_command"][-1] == "-nofit"
    assert frame.data_summary["stdin"] == "Peptide_Backbone\n"
    assert r.definition_token and r.arrays == []           # planned: nothing realised


def test_rmsd_result_array_describes_existing_xvg(tmp_path):
    view_file = tmp_path / "whole.xtc"
    view_file.write_bytes(b"\x00" * 32)
    ctx = _ctx(tmp_path, view_path=str(view_file))
    xvg = tmp_path / "rmsd.xvg"
    xvg.write_text(_XVG)
    before = xvg.read_bytes()
    spec = registry.get("rmsd-receptor")
    arr = spec._result_array(ctx, xvg, "tok")
    assert xvg.read_bytes() == before                      # file untouched
    assert arr.quantity == "rmsd" and arr.unit == "nm" and arr.validate() == []
    t = arr.axis("time")
    assert t.kind == AxisKind.TIME and t.unit == "ps"
    assert read_column(t.values_ref) == [100.0, 100.0, 100.5]
    assert read_column(arr.storage) == [0.1, 0.2, 0.3]
    assert arr.storage.path == str(xvg.resolve()) and arr.storage.column == 1
    assert arr.view_ref == "view-1" and t.alignment.view_ref == "view-1"
    assert t.alignment.mode == FrameAlignmentMode.ONE_ROW_PER_FRAME
    assert t.alignment.trajectory_digest
    assert arr.definition_token == "tok" and arr.missing == MissingSemantics.NONE_EXPECTED


def test_rmsd_result_array_unit_undeclared(tmp_path):
    ctx = _ctx(tmp_path)
    xvg = tmp_path / "rmsd.xvg"
    xvg.write_text("0.0 0.1\n20.0 0.2\n")
    arr = registry.get("rmsd-receptor")._result_array(ctx, xvg, None)
    assert arr.unit is None and arr.attrs["unit_source"] == "undeclared"
    assert arr.axis("time").unit == "ps"
    assert "-tu ps" in arr.axis("time").attrs["unit_source"]


# ═══════════════════════════════════════════════════════════════════════════════
# Real RMSD run: identical numbers, typed description, deterministic token
# ═══════════════════════════════════════════════════════════════════════════════

@pytest.fixture
def real_study(tmp_path):
    root = tmp_path / "campaign"
    d = root / "semaglutide" / "rep1"
    d.mkdir(parents=True)
    (d / "production_run.xtc").symlink_to(REAL_XTC)
    (d / "system_topol.tpr").symlink_to(REAL_TPR)
    (d / "reference.gro").symlink_to(REAL_GRO)
    return root


def _data_rows(path: Path) -> list[str]:
    return [l for l in path.read_text().splitlines() if l and l[0] not in "#@"]


def _headers(path: Path) -> list[str]:
    return [l for l in path.read_text().splitlines() if l.startswith("@")]


@requires_gmx
@requires_real_traj
def test_real_rmsd_migration(real_study, tmp_path):
    from analysis.campaign.gmx import run_gmx
    from analysis.campaign.orchestration.study_analyzer import clear_view_cache, run_analyze
    from analysis.campaign.trajectory.time_index import index_trajectory

    clear_view_cache()
    res = run_analyze(real_study, ["rmsd-receptor"], inspect_trajectories=True)
    (r,) = res.results
    assert r.status == AnalysisStatus.SUCCESS, r.message
    xvg = Path(next(f for f in r.output_files if f.endswith("rmsd.xvg")))
    assert xvg.name == "rmsd.xvg"

    # legacy summary intact
    assert r.data_summary["n_frames"] == 51 and "mean_nm" in r.data_summary

    # the XVG is exactly what an independent run of the recorded command produces
    log = (xvg.parent / "gmx_rms.log").read_text().splitlines()
    argv = log[0][2:].split()
    assert argv[1:3] == ["rms", "-s"] and argv[-2:] == ["-tu", "ps"]
    indep = tmp_path / "indep.xvg"
    argv[argv.index("-o") + 1] = str(indep)
    rerun = run_gmx(argv[1:], stdin="Receptor_Backbone\nReceptor_Backbone\n")
    assert rerun.ok
    assert _data_rows(indep) == _data_rows(xvg) and _headers(indep) == _headers(xvg)

    # typed description of that same file
    (arr,) = r.arrays
    data = parse_xvg(xvg)
    assert arr.unit == "nm" and arr.axis("time").unit == "ps"
    assert read_column(arr.storage) == data.series[0].values
    assert read_column(arr.axis("time").values_ref) == data.time_ps
    assert arr.storage.path == str(xvg.resolve())
    prov = json.loads(Path(r.provenance_path).read_text())
    assert arr.view_ref == prov["trajectory_view"]["cache_key"]
    assert prov["result"]["arrays"][0]["quantity"] == "rmsd"
    assert prov["result"]["definition_token"] == r.definition_token
    # Phase 5: same view and numbers, plus the decision that authorised it
    assert prov["trajectory_view"]["policy_planned"] is True
    (dec,) = prov["trajectory_view"]["decisions"]
    assert (dec["operation"], dec["classification"], dec["rule_id"]) == (
        "make_whole", "valid", "make_whole.tpr_connectivity")
    assert dec["purpose"] == "intramolecular_shape" and dec["applied"] is True

    # one row per frame of the analysed view file (Phase 1 time index)
    view_path = prov["trajectory_view"]["path"]
    tidx, _ = index_trajectory(view_path)
    assert tidx.times_ps == data.time_ps
    assert tidx.fingerprint.digest == arr.axis("time").alignment.trajectory_digest

    # deterministic token across a second run that reuses the cached view
    clear_view_cache()
    r2 = run_analyze(real_study, ["rmsd-receptor"], inspect_trajectories=True).results[0]
    assert r2.definition_token == r.definition_token
    assert _data_rows(Path(r2.arrays[0].storage.path)) == _data_rows(indep)
