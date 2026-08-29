"""Multi-level study validation."""
from __future__ import annotations

import pytest

from analysis.campaign.models import (
    ClassificationState, ComponentType, MolecularComponent, StudyManifest,
    SystemRecord, ValidationState,
)
from analysis.campaign.validator import validate_study


def _component(ctype, state=ClassificationState.RESOLVED, chain="A"):
    return MolecularComponent(
        component_type=ctype, label=ctype, chain_ids=[chain],
        classification_state=state, confidence=0.9,
    )


def _system(sid, components, *, topology="t.tpr", structure="c.pdb", traj=("md.xtc",)):
    return SystemRecord(
        system_id=sid, condition_id="cond", replicate_id="rep01",
        receptor="R", partner="P", water_model="TIP3P",
        topology_path=topology, structure_path=structure,
        production_trajectory_paths=list(traj), components=components,
    )


def _manifest(systems):
    return StudyManifest(study_root="/x", generated_utc="", simforge_version="test",
                         systems=systems)


def test_all_resolved_receptor_systems_valid():
    m = _manifest([
        _system("s1", [_component(ComponentType.RECEPTOR)]),
        _system("s2", [_component(ComponentType.RECEPTOR)]),
    ])
    rep = validate_study(m, ["rmsd-receptor"])
    assert all(v == ValidationState.VALID for v in rep.system_states.values())
    assert rep.study_state in (ValidationState.VALID, ValidationState.WARNING)
    for o in rep.observable_statuses:
        assert o.state == ValidationState.VALID
    assert rep.counts["rmsd-receptor"]["ready"] == 2


def test_peptide_problem_does_not_block_receptor_analysis():
    m = _manifest([
        _system("s1", [
            _component(ComponentType.RECEPTOR),
            _component(ComponentType.PEPTIDE, ClassificationState.AMBIGUOUS),
        ]),
    ])
    rep = validate_study(m, ["rmsd-receptor", "rmsd-peptide-intrinsic"])
    by = {(o.system_id, o.analysis_id): o for o in rep.observable_statuses}
    assert by[("s1", "rmsd-receptor")].state == ValidationState.VALID
    # the peptide analysis is flagged, not silently valid
    assert by[("s1", "rmsd-peptide-intrinsic")].state != ValidationState.VALID
    assert rep.counts["rmsd-peptide-intrinsic"]["review_required"] == 1


def test_peptide_absent_is_skipped_not_review():
    m = _manifest([_system("s1", [_component(ComponentType.RECEPTOR)])])
    rep = validate_study(m, ["rmsd-peptide-intrinsic"])
    o = rep.observable_statuses[0]
    assert o.state in (ValidationState.INVALID,)  # "not present" -> skipped bucket
    assert rep.counts["rmsd-peptide-intrinsic"]["skipped"] == 1
    assert rep.counts["rmsd-peptide-intrinsic"]["review_required"] == 0


def test_unknown_analysis_id_error_but_others_validate():
    m = _manifest([_system("s1", [_component(ComponentType.RECEPTOR)])])
    rep = validate_study(m, ["rmsd-receptor", "totally-bogus"])
    assert any(w.code == "unknown_analysis" and w.severity == "error" for w in rep.warnings)
    assert "rmsd-receptor" in rep.counts
    assert rep.counts["rmsd-receptor"]["ready"] == 1


def test_system_without_topology_is_invalid():
    m = _manifest([_system("s1", [_component(ComponentType.RECEPTOR)], topology=None)])
    rep = validate_study(m, ["rmsd-receptor"])
    assert rep.system_states["s1"] == ValidationState.INVALID
    o = rep.observable_statuses[0]
    assert o.state == ValidationState.INVALID
    assert o.reason


def test_report_lines_and_serialisation():
    m = _manifest([_system("s1", [_component(ComponentType.RECEPTOR)])])
    rep = validate_study(m, ["rmsd-receptor"])
    assert isinstance(rep.lines, list) and rep.lines
    assert all(isinstance(x, str) for x in rep.lines)
    d = rep.to_dict()
    assert d["study_state"] == rep.study_state
    assert "counts" in d and "observable_statuses" in d
    assert d["counts"]["rmsd-receptor"].keys() >= {"ready", "skipped", "review_required"}


def test_empty_manifest_is_invalid():
    rep = validate_study(_manifest([]), ["rmsd-receptor"])
    assert rep.study_state == ValidationState.INVALID
