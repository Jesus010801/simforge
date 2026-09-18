"""Dry-run mechanistic planner: mdfit.xtc is mandatory, nothing is executed."""
from __future__ import annotations

from pathlib import Path

from analysis.campaign import mechanistic as mp
from tests.analysis.campaign.test_legacy import system


def _sys(tmp_path, name="A3-HMG-R", *, with_fit=True):
    d, _ = system(tmp_path, name, with_fit=with_fit)
    (d / "md.tpr").write_bytes(b"tpr")
    return d


def test_planner_never_executes_and_uses_mdfit(tmp_path):
    d = _sys(tmp_path)
    report = mp.plan_mechanistic(tmp_path)
    assert report["execution"] == "disabled; dry-run planning only"
    assert report["analysis_trajectory"] == "mdfit.xtc"
    assert report["raw_md_xtc_is_forbidden_as_analysis_input"] is True
    s = next(s for s in report["systems"] if s["system"] == "A3-HMG-R")
    assert Path(s["analysis_trajectory"]).name == "mdfit.xtc"
    for o in s["observables"]:
        if o["source_trajectory"] is not None:
            assert Path(o["source_trajectory"]).name == "mdfit.xtc"
            assert "md.xtc" not in o["source_trajectory"]


def test_missing_mdfit_blocks_every_non_extension_observable(tmp_path):
    _sys(tmp_path, "A3-HMG-R", with_fit=False)
    report = mp.plan_mechanistic(tmp_path)
    s = next(s for s in report["systems"] if s["system"] == "A3-HMG-R")
    assert s["analysis_trajectory"] is None
    non_ext = [o for o in s["observables"]
               if o["status"] != mp.NEW_EXTENSION]
    assert non_ext and all(o["status"] == mp.REVIEW_REQUIRED for o in non_ext)
    assert all("mdfit.xtc" in o["reason"] for o in non_ext)


def test_clustering_and_state_populations_are_new_extensions(tmp_path):
    _sys(tmp_path)
    report = mp.plan_mechanistic(tmp_path)
    s = report["systems"][0]
    ext = {o["observable"]: o for o in s["observables"]
           if o["status"] == mp.NEW_EXTENSION}
    assert set(ext) == {"clustering", "conformational_state_populations"}
    for o in ext.values():
        assert "not historically present for A1" in o["reason"]


def test_historical_raw_output_is_superseded_pending(tmp_path):
    ref = tmp_path / "Mecanismo_inhibitorio"
    a1 = ref / "A1"
    a1.mkdir(parents=True)
    (a1 / "md.xtc").write_bytes(b"raw 200 ns")
    (a1 / "md.tpr").write_bytes(b"tpr")
    (a1 / "index.ndx").write_text("[ Protein ]\n1\n[ LIG ]\n2\n")
    (a1 / "mdfit.xtc").write_bytes(b"fit")            # reprocessed reference
    (a1 / "rmsd_A1.xvg").write_text("@\n0 0.1\n")     # historical, raw-derived
    (a1 / "rmsf_A1.xvg").write_text("@\n1 0.1\n")

    report = mp.plan_mechanistic(tmp_path / "nonexistent", reference=ref)
    s = next(s for s in report["systems"] if s["system"] == "A1")
    st = {o["observable"]: o["status"] for o in s["observables"]}
    assert st["protein_rmsd"] == mp.SUPERSEDED_PENDING
    assert st["rmsf"] == mp.SUPERSEDED_PENDING
    # an observable with no historical file is simply PLANNED on mdfit.xtc
    assert st["radius_of_gyration"] == mp.PLANNED
    assert st["clustering"] == mp.NEW_EXTENSION


def test_window_is_independent_of_production_length(tmp_path):
    # a 200 ns fixture still plans the same profile; the planner does not
    # shorten or lengthen the profile based on production length.
    d, _ = system(tmp_path, "HMG-R-200ns-A6", duration=200000, with_fit=True)
    (d / "md.tpr").write_bytes(b"tpr")
    report = mp.plan_mechanistic(tmp_path)
    s = next(s for s in report["systems"] if s["system"] == "HMG-R-200ns-A6")
    keys = {o["observable"] for o in s["observables"]}
    assert {"protein_rmsd", "pca", "fel", "dccm", "mmpbsa",
            "clustering", "conformational_state_populations"} <= keys
