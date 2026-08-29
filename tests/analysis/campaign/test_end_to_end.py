"""Full pipeline against the bundled protein-membrane trajectory (needs gmx)."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from analysis.campaign.models import AnalysisStatus
from analysis.campaign.orchestration.study_analyzer import clear_view_cache, run_analyze
from tests.analysis.campaign.conftest import (
    REAL_GRO, REAL_TPR, REAL_XTC, requires_gmx, requires_real_traj,
)

pytestmark = [requires_gmx, requires_real_traj]


@pytest.fixture
def two_replicate_study(tmp_path):
    root = tmp_path / "campaign"
    for rep in ("rep1", "rep2"):
        d = root / "semaglutide" / rep
        d.mkdir(parents=True)
        (d / "production_run.xtc").symlink_to(REAL_XTC)
        (d / "system_topol.tpr").symlink_to(REAL_TPR)
        (d / "reference.gro").symlink_to(REAL_GRO)
    return root


def test_full_rmsd_receptor_run(two_replicate_study):
    res = run_analyze(two_replicate_study, ["rmsd-receptor"], inspect_trajectories=True)
    assert len(res.results) == 2
    for r in res.results:
        assert r.status == AnalysisStatus.SUCCESS, r.message
        assert r.data_summary["n_frames"] == 51
        assert "mean_nm" in r.data_summary
        assert r.trajectory_view_kind == "whole"

        xvg = [f for f in r.output_files if f.endswith("rmsd.xvg")]
        assert xvg and Path(xvg[0]).is_file()

        prov = json.loads(Path(r.provenance_path).read_text())
        assert prov["environment"]["gromacs_version"]
        assert prov["trajectory_view"]["kind"] == "whole"
        ops = [o["operation"] for o in prov["trajectory_view"]["operations"]]
        assert "make_whole" in ops
        assert "Receptor_Backbone" in [g["name"] for g in prov["semantic_index"]["groups"]]

        ndx = Path(prov["semantic_index"]["path"])
        assert ndx.is_file() and ndx.name == "semantic_index.ndx"


def test_second_run_reuses_cached_trajectory_view(two_replicate_study):
    run_analyze(two_replicate_study, ["rmsd-receptor"], inspect_trajectories=True)
    whole = next(two_replicate_study.rglob("whole.xtc"))
    mtime_after_first = whole.stat().st_mtime_ns

    clear_view_cache()          # drop the in-process cache; disk cache remains
    res2 = run_analyze(two_replicate_study, ["rmsd-receptor"], inspect_trajectories=True)

    prov = json.loads(Path(res2.results[0].provenance_path).read_text())
    assert prov["trajectory_view"]["reused"] is True
    assert whole.stat().st_mtime_ns == mtime_after_first


def test_peptide_analysis_on_receptor_only_system_is_skipped(two_replicate_study):
    res = run_analyze(two_replicate_study, ["rmsd-peptide-intrinsic"],
                      inspect_trajectories=True)
    assert res.results
    for r in res.results:
        assert r.status in (AnalysisStatus.SKIPPED, AnalysisStatus.REVIEW_REQUIRED)
        assert r.message
        # no misleading xvg produced
        assert not any(f.endswith("rmsd.xvg") for f in r.output_files)
