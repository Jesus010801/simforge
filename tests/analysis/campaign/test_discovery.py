"""Study discovery: filename-independence, association evidence, traversal order."""
from __future__ import annotations

from pathlib import Path

import pytest

from analysis.campaign.discovery import discover_study
from analysis.campaign.manifest import build_manifest
from tests.analysis.campaign.conftest import write_gro, write_pdb, write_top


def _cand_by_dir(result):
    return {Path(c.sim_dir).name: c for c in result.candidates}


def test_all_trajectories_discovered_regardless_of_name(study_tree):
    res = discover_study(study_tree)
    assert len(res.candidates) == 4
    names = {Path(p).name for c in res.candidates for p in c.production_trajectory_paths}
    assert names == {"final.xtc", "production.xtc", "md_100ns.xtc", "run3.trr"}


def test_filename_never_defines_identity(synthetic_study_tree):
    root = synthetic_study_tree([
        ("GLP1", "TIP3P", "rep1", "semaglutide.xtc", "topol.top", "conf.pdb"),
    ])
    res = discover_study(root)
    cand = res.candidates[0]
    # partner comes from the directory (GLP1), NOT the misleading filename
    assert cand.partner_hint == "GLP1"
    assert "semaglutide" not in (cand.partner_hint or "").lower()


def test_different_topology_names_still_associated(study_tree):
    res = discover_study(study_tree)
    for c in res.candidates:
        assert c.topology_path is not None
        assert any("same_directory" in e for e in c.association_evidence)


def test_stage_matched_tpr_is_unambiguous(synthetic_study_tree):
    """A production trajectory 'md.xtc' next to 'md.tpr' is a clean stage match —
    the presence of other-stage tprs is NOT an ambiguity."""
    root = synthetic_study_tree([("GLP1", "TIP3P", "rep1", "md.xtc", "", "conf.pdb")],
                                make_topology=False)
    d = root / "GLP1" / "TIP3P" / "rep1"
    (d / "md.tpr").write_bytes(b"TPR1")
    (d / "nvt.tpr").write_bytes(b"TPR2")
    (d / "npt.tpr").write_bytes(b"TPR3")

    r1 = discover_study(root)
    r2 = discover_study(root)
    c1, c2 = r1.candidates[0], r2.candidates[0]
    assert c1.topology_path == c2.topology_path
    assert Path(c1.topology_path).name == "md.tpr"
    assert not any(w.code in ("ambiguous_topology", "ambiguous_production_topology")
                   for w in c1.warnings)


def test_multiple_production_tpr_without_stem_match_is_review_required(synthetic_study_tree):
    """Two production-stage tprs, neither stem-matching the trajectory -> review_required,
    deterministic pick, never silent."""
    root = synthetic_study_tree([("GLP1", "TIP3P", "rep1", "traj_final.xtc", "", "conf.pdb")],
                                make_topology=False)
    d = root / "GLP1" / "TIP3P" / "rep1"
    (d / "production_a.tpr").write_bytes(b"TPR1")
    (d / "run_b.tpr").write_bytes(b"TPR2")

    r1 = discover_study(root)
    r2 = discover_study(root)
    assert r1.candidates[0].topology_path == r2.candidates[0].topology_path
    assert any(w.code == "ambiguous_production_topology" and w.severity == "review_required"
               for w in r1.candidates[0].warnings)


def test_missing_topology_warns_but_candidate_created(synthetic_study_tree):
    root = synthetic_study_tree([("GLP1", "TIP3P", "rep1", "md.xtc", "", "conf.pdb")],
                                make_topology=False)
    res = discover_study(root)
    assert len(res.candidates) == 1
    assert res.candidates[0].topology_path is None
    assert any(w.code == "no_topology" for w in res.candidates[0].warnings)


def test_stray_trajectory_in_output_dir_is_own_candidate(synthetic_study_tree):
    root = synthetic_study_tree()
    stray = root / "analysis_output"
    stray.mkdir()
    (stray / "leftover.xtc").write_bytes(b"FAKE")
    res = discover_study(root)
    dirs = {Path(c.sim_dir).name for c in res.candidates}
    assert "analysis_output" in dirs
    assert len(res.candidates) == 5


def test_traversal_order_independence(study_tree, monkeypatch):
    baseline = build_manifest(study_tree, inspect_trajectories=False, gmx="/nonexistent")
    baseline_ids = {s.system_id: sorted(s.association_evidence) for s in baseline.systems}

    real_rglob = Path.rglob
    real_iterdir = Path.iterdir
    monkeypatch.setattr(Path, "rglob",
                        lambda self, pat: list(reversed(list(real_rglob(self, pat)))))
    monkeypatch.setattr(Path, "iterdir",
                        lambda self: list(reversed(list(real_iterdir(self)))))

    shuffled = build_manifest(study_tree, inspect_trajectories=False, gmx="/nonexistent")
    shuffled_ids = {s.system_id: sorted(s.association_evidence) for s in shuffled.systems}
    assert shuffled_ids == baseline_ids


def test_trajectory_segments_single_candidate_ordered(synthetic_study_tree):
    root = synthetic_study_tree([("GLP1", "TIP3P", "rep1", "md.part0001.xtc", "topol.top", "conf.pdb")])
    d = root / "GLP1" / "TIP3P" / "rep1"
    (d / "md.part0002.xtc").write_bytes(b"P2")
    (d / "md.part0003.xtc").write_bytes(b"P3")

    res = discover_study(root)
    assert len(res.candidates) == 1
    paths = [Path(p).name for p in res.candidates[0].production_trajectory_paths]
    assert paths == ["md.part0001.xtc", "md.part0002.xtc", "md.part0003.xtc"]
    assert any("segment" in e for e in res.candidates[0].discovery_evidence)


def test_derived_trajectories_separated_from_primary(synthetic_study_tree):
    root = synthetic_study_tree([("GLP1", "TIP3P", "rep1", "md.xtc", "topol.top", "conf.pdb")])
    d = root / "GLP1" / "TIP3P" / "rep1"
    (d / "md_nojump.xtc").write_bytes(b"NJ")
    (d / "md_center.xtc").write_bytes(b"CN")

    res = discover_study(root)
    cand = res.candidates[0]
    assert [Path(p).name for p in cand.production_trajectory_paths] == ["md.xtc"]
    assert {Path(p).name for p in cand.derived_trajectory_paths} == {"md_nojump.xtc", "md_center.xtc"}
    assert any("pre-processed" in e or "PBC/fit" in e for e in cand.discovery_evidence)


def test_nested_layout(synthetic_study_tree, tmp_path):
    root = tmp_path / "study"
    deep = root / "proteinX" / "conditionY" / "rep1" / "deeper"
    deep.mkdir(parents=True)
    (deep / "traj.xtc").write_bytes(b"FAKE")
    write_top(deep / "topol.top")
    write_pdb(deep / "conf.pdb", {"A": 30})

    res = discover_study(root)
    assert len(res.candidates) == 1
    c = res.candidates[0]
    assert c.replicate_id == "rep01"
    assert c.production_trajectory_paths


@pytest.mark.parametrize("marker,expected", [
    ("rep1", "rep01"), ("replica2", "rep02"), ("r3", "rep03"),
    ("run4", "rep04"), ("3", "rep03"),
])
def test_replicate_marker_detection(synthetic_study_tree, tmp_path, marker, expected):
    root = tmp_path / f"s_{marker}"
    d = root / "GLP1" / marker
    d.mkdir(parents=True)
    (d / "md.xtc").write_bytes(b"FAKE")
    res = discover_study(root)
    assert res.candidates[0].replicate_id == expected


def test_no_replicate_marker_defaults_rep01_with_warning(tmp_path):
    root = tmp_path / "study"
    d = root / "GLP1"
    d.mkdir(parents=True)
    (d / "md.xtc").write_bytes(b"FAKE")
    res = discover_study(root)
    c = res.candidates[0]
    assert c.replicate_id == "rep01"
    assert any(w.code == "no_replicate_marker" for w in c.warnings)


@pytest.mark.parametrize("dirname,model", [
    ("TIP3P", "TIP3P"), ("tip3p", "TIP3P"), ("SPCE", "SPC/E"), ("tip4p", "TIP4P"),
])
def test_water_model_detection_case_insensitive(tmp_path, dirname, model):
    root = tmp_path / "study"
    d = root / "GLP1" / dirname / "rep1"
    d.mkdir(parents=True)
    (d / "md.xtc").write_bytes(b"FAKE")
    res = discover_study(root)
    assert res.candidates[0].water_model == model


def test_bare_integer_dir_is_not_mistaken_for_replicate(tmp_path):
    root = tmp_path / "study"
    for temp in ("300", "310"):
        d = root / "GLP1" / temp / "rep1"
        d.mkdir(parents=True)
        (d / "md.xtc").write_bytes(b"F")
    res = discover_study(root)
    reps = {c.replicate_id for c in res.candidates}
    assert reps == {"rep01"}                       # both are replicate 1
    dims = [c.condition_dimensions for c in res.candidates]
    assert any("300" in str(d.values()) for d in dims)


def test_missing_root_returns_error_warning(tmp_path):
    res = discover_study(tmp_path / "nope")
    assert res.candidates == []
    assert any(w.code == "study_root_missing" for w in res.warnings)


def test_no_trajectories_error(tmp_path):
    (tmp_path / "topol.top").write_text("[ molecules ]\nSOL 1\n")
    res = discover_study(tmp_path)
    assert res.candidates == []
    assert any(w.code == "no_trajectories" for w in res.warnings)
