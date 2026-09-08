"""Adversarial index and discovery checks; all fixtures are disposable."""
import json
from pathlib import Path
import pytest
from analysis.campaign.discovery import discover_study
from analysis.campaign.manifest import build_manifest
from analysis.campaign.legacy import annotate_legacy
from analysis.campaign.structure.index_groups import parse_index_groups, resolve_index_group
from tests.analysis.campaign.test_legacy import system, snapshot

@pytest.mark.parametrize("content", ["", "; comments only", "[ ]\n1", "[[ LIG ]]\n1", "1 2\n[ LIG ]", "[ LIG ]\n0", "[ LIG ]\n-1", "[ LIG ]\none", "[ LIG \n1"])
def test_malformed_index_is_rejected(tmp_path, content):
    d, rec = system(tmp_path)
    (d / "index.ndx").write_text(content)
    with pytest.raises(ValueError):
        parse_index_groups(d / "index.ndx")
    report = annotate_legacy(rec, d)
    assert report["status"] == "review_required"
    assert not report["required_groups_found"]
    assert any("invalid index" in x for x in report["issues"])

@pytest.mark.parametrize("ligand_id,site_id", [(0, 3), (14, 23), (20, 10)])
def test_every_role_uses_actual_header_position(tmp_path, ligand_id, site_id):
    d, rec = system(tmp_path, "HMG-R-25ns-A6", ligand_id, site_id)
    groups = parse_index_groups(d / "index.ndx")
    expected = {g.name: g.group_id for g in groups}
    report = annotate_legacy(rec, d)
    for role in report["selections"].values():
        if role["value"]:
            assert role["group_id"] == expected[role["value"]]
            assert str(d / "index.ndx") in role["provenance"][0]
    assert report["selections"]["ligand"]["group_id"] == ligand_id
    assert report["selections"]["active_site"]["group_id"] == site_id


def test_controlled_aliases_and_duplicates(tmp_path):
    d, rec = system(tmp_path)
    p = d / "index.ndx"
    p.write_text(p.read_text().replace("[ LIG ]", "[ A6 ]"))
    report = annotate_legacy(rec, d)
    assert report["selections"]["ligand"]["value"] == "A6"
    assert report["selections"]["ligand"]["group_id"] == 13
    p.write_text(p.read_text() + "[ A6 ]\n1\n")
    assert annotate_legacy(rec, d)["selections"]["ligand"]["value"] is None
    with pytest.raises(ValueError, match="ambiguous"):
        resolve_index_group(p, "LIG", ("A6",))


def test_canonical_backup_and_multiple_index_discovery(tmp_path):
    d, _ = system(tmp_path)
    (d / "#index.ndx.1#").write_text("[ LIG ]\n1\n")
    assert discover_study(tmp_path).candidates[0].index_path == str(d / "index.ndx")
    (d / "alternate.ndx").write_text("[ LIG ]\n1\n")
    candidate = discover_study(tmp_path).candidates[0]
    assert candidate.index_path is None
    assert any(w.code == "ambiguous_index" for w in candidate.warnings)


def test_missing_trajectory_stays_visible(tmp_path):
    d, _ = system(tmp_path)
    (d / "md.xtc").unlink()
    (d / "md.tpr").write_bytes(b"placeholder")
    rec = build_manifest(tmp_path, inspect_trajectories=False).systems[0]
    assert rec.production_trajectory_paths == []
    assert rec.legacy_study["status"] == "review_required"
    assert any("no trajectory" in x for x in rec.legacy_study["issues"])


def test_multiple_production_requires_review(tmp_path):
    d, _ = system(tmp_path)
    (d / "production.xtc").write_bytes(b"other run")
    rec = build_manifest(tmp_path, inspect_trajectories=False).systems[0]
    assert len(rec.production_trajectory_paths) == 2
    assert rec.legacy_study["status"] == "review_required"
    assert any(w.code == "ambiguous_production_trajectories" for w in rec.warnings)


def test_spaces_determinism_and_source_integrity(tmp_path):
    root = tmp_path / "legacy project with spaces"
    root.mkdir()
    system(root, "AA-A6 sample")
    before = snapshot(root)
    reports = [build_manifest(root, inspect_trajectories=False).systems[0].legacy_study for _ in range(2)]
    assert json.dumps(reports[0], sort_keys=True) == json.dumps(reports[1], sort_keys=True)
    assert snapshot(root) == before


def test_empty_xtc_does_not_silently_validate_alternate_trr(tmp_path):
    d, _ = system(tmp_path)
    (d / "md.xtc").write_bytes(b"")
    (d / "md.trr").write_bytes(b"nonempty alternate")
    rec = build_manifest(tmp_path, inspect_trajectories=False).systems[0]
    assert rec.production_trajectory_paths == [str(d / "md.xtc")]
    assert len(rec.trajectory_artifacts) == 2
    assert any(w.code == "empty_production_trajectory" for w in rec.warnings)
    assert rec.legacy_study["status"] == "review_required"


def test_nearby_index_does_not_borrow_sibling(tmp_path):
    left, _ = system(tmp_path, "AA-A6")
    right, _ = system(tmp_path, "AG-A6")
    (right / "index.ndx").unlink()
    candidates = {Path(c.sim_dir).name: c for c in discover_study(tmp_path).candidates}
    assert candidates[left.name].index_path == str(left / "index.ndx")
    assert candidates[right.name].index_path is None
