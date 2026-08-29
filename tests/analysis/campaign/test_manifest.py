"""StudyManifest: build, group, serialise, reload, honour hand-edits."""
from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from analysis.campaign.manifest import build_manifest, load_manifest, write_manifest
from analysis.campaign.models import ClassificationState, StudyManifest
from tests.analysis.campaign.conftest import write_pdb, write_top


NONGMX = dict(inspect_trajectories=False, gmx="/nonexistent")


def test_build_manifest_systems_and_conditions(study_tree):
    m = build_manifest(study_tree, **NONGMX)
    assert len(m.systems) == 4
    # GLP1+TIP3P has 2 replicates -> one condition with both system ids
    conds = {c.condition_id: c for c in m.conditions}
    two_rep = [c for c in conds.values() if len(c.system_ids) == 2]
    assert len(two_rep) == 1
    assert set(two_rep[0].system_ids) == {
        s.system_id for s in m.systems
        if s.partner == "GLP1" and s.water_model == "TIP3P"
    }
    # every system's condition_id points back to a real condition
    for s in m.systems:
        assert s.condition_id in conds


def test_canonical_ids_unique_and_deterministic(study_tree):
    m1 = build_manifest(study_tree, **NONGMX)
    m2 = build_manifest(study_tree, **NONGMX)
    ids1 = [s.system_id for s in m1.systems]
    ids2 = [s.system_id for s in m2.systems]
    assert ids1 == ids2
    assert len(set(ids1)) == len(ids1)


def test_to_dict_from_dict_roundtrip(study_tree):
    m = build_manifest(study_tree, **NONGMX)
    d = m.to_dict()
    back = StudyManifest.from_dict(d)
    assert [s.system_id for s in back.systems] == [s.system_id for s in m.systems]
    assert [c.condition_id for c in back.conditions] == [c.condition_id for c in m.conditions]
    assert len(back.ambiguities) == len(m.ambiguities)
    assert len(back.warnings) == len(m.warnings)
    for s0, s1 in zip(m.systems, back.systems):
        assert [c.component_type for c in s0.components] == [c.component_type for c in s1.components]


def test_write_manifest_yaml_and_json(study_tree, tmp_path):
    m = build_manifest(study_tree, **NONGMX)
    paths = write_manifest(m, tmp_path / "out")
    assert Path(paths["yaml"]).is_file() and Path(paths["json"]).is_file()
    y = yaml.safe_load(Path(paths["yaml"]).read_text())
    j = json.loads(Path(paths["json"]).read_text())
    assert y["stats"]["n_systems"] == j["stats"]["n_systems"] == 4

    reloaded = load_manifest(paths["yaml"])
    assert [s.system_id for s in reloaded.systems] == [s.system_id for s in m.systems]


def _ambiguous_tree(tmp_path):
    root = tmp_path / "study"
    sd = root / "GLP1" / "TIP3P" / "rep1"
    sd.mkdir(parents=True)
    (sd / "md.xtc").write_bytes(b"FAKE")
    write_top(sd / "topol.top")
    write_pdb(sd / "conf.pdb", {"A": 120, "B": 25})
    return root


def test_ambiguities_populated_when_classification_ambiguous(tmp_path):
    m = build_manifest(_ambiguous_tree(tmp_path), **NONGMX)
    assert m.ambiguities
    assert all(a.kind == "component_classification" for a in m.ambiguities)
    assert m.systems[0].classification_state == ClassificationState.AMBIGUOUS


def test_manifest_is_a_contract_hand_edited_resolution_is_honoured(tmp_path):
    m = build_manifest(_ambiguous_tree(tmp_path), **NONGMX)
    out = tmp_path / "out"
    paths = write_manifest(m, out)

    # hand-edit: resolve the ambiguity
    loaded = load_manifest(paths["yaml"])
    for a in loaded.ambiguities:
        a.resolution = "receptor=A;peptide=B"
    write_manifest(loaded, out)

    final = load_manifest(paths["yaml"])
    rec = final.systems[0]
    assert rec.classification_state == "resolved"
    assert rec.user_overridden is True
    receptor = rec.component("receptor")
    peptide = rec.component("peptide")
    assert receptor.chain_ids == ["A"] and receptor.classification_state == "resolved"
    assert peptide.chain_ids == ["B"] and peptide.classification_state == "resolved"


def test_unbalanced_replicates_is_info_not_invalid(synthetic_study_tree):
    root = synthetic_study_tree([
        ("GLP1", "TIP3P", "rep1", "a.xtc", "topol.top", "c.pdb"),
        ("GLP1", "TIP3P", "rep2", "b.xtc", "topol.top", "c.pdb"),
        ("GLP1", "TIP3P", "rep3", "c.xtc", "topol.top", "c.pdb"),
        ("Semaglutide", "TIP3P", "rep1", "d.xtc", "topol.top", "c.pdb"),
    ])
    m = build_manifest(root, **NONGMX)
    unbalanced = [w for w in m.warnings if w.code == "unbalanced_replicates"]
    assert len(unbalanced) == 1
    assert unbalanced[0].severity == "info"


def test_stats_shape(study_tree):
    m = build_manifest(study_tree, **NONGMX)
    stats = m.stats()
    assert stats["partners"] == sorted(stats["partners"])
    assert set(stats["water_models"]) == {"SPC/E", "TIP3P"}
    assert stats["n_systems"] == 4
    assert isinstance(stats["balanced"], bool)
    assert isinstance(stats["replicate_counts"], dict)
