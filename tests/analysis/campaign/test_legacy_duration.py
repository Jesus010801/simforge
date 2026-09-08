"""Adversarial evidence precedence and missing-production checks."""
from pathlib import Path
import pytest
from analysis.campaign.legacy_duration import duration_evidence
from analysis.campaign.models import SystemRecord, TrajectoryInspection, CampaignWarning


def record(tmp_path, ns=200):
    trajectory = tmp_path / "md.xtc"
    trajectory.write_bytes(b"fixture")
    (tmp_path / "md.mdp").write_text(f"dt = .001\nnsteps = {ns * 1000000}\n")
    return SystemRecord("test", "test", "rep01", production_trajectory_paths=[str(trajectory)])


def test_measured_duration_preserves_conflicting_mdp(tmp_path):
    rec = record(tmp_path)
    rec.trajectory_inspection = TrajectoryInspection(inspected=True, total_duration_ps=25000)
    duration, issues = duration_evidence(rec, tmp_path)
    assert duration["value"] == 25
    assert duration["confidence"] == .95
    assert any("200 ns" in source for source in duration["provenance"])
    assert "measured trajectory duration conflicts with intended production MDP duration" in issues


def test_failed_check_never_promotes_partial_timing(tmp_path):
    rec = record(tmp_path)
    rec.trajectory_inspection = TrajectoryInspection(inspected=True, total_duration_ps=25000,
        warnings=[CampaignWarning("gmx_check_failed", "truncated", "warn")])
    duration, issues = duration_evidence(rec, tmp_path)
    assert duration["value"] == 200
    assert duration["confidence"] == .65
    assert any("inspection failed" in issue for issue in issues)


def test_missing_production(tmp_path):
    rec = record(tmp_path)
    rec.production_trajectory_paths = []
    assert "missing production trajectory" in duration_evidence(rec, tmp_path)[1]


@pytest.mark.parametrize("ns", [-1, float("nan"), float("inf"), 0])
def test_invalid_mdp_duration(tmp_path, ns):
    rec = record(tmp_path, ns)
    duration, issues = duration_evidence(rec, tmp_path)
    assert duration["value"] is None
    assert any("invalid or unbounded" in issue for issue in issues)


def test_multiple_production_mdp_conflict(tmp_path):
    rec = record(tmp_path)
    second = tmp_path / "prod2.xtc"
    second.write_bytes(b"fixture")
    second.with_suffix(".mdp").write_text("dt = .001\nnsteps = 25000000\n")
    rec.production_trajectory_paths.append(str(second))
    assert "conflicting production MDP durations" in duration_evidence(rec, tmp_path)[1]


def test_tpr_conflict_preserved(tmp_path):
    rec = record(tmp_path, 25)
    rec.trajectory_inspection = TrajectoryInspection(inspected=True, total_duration_ps=200000,
        tpr_duration_ps=200000, topology_path=str(tmp_path / "md.tpr"))
    duration, issues = duration_evidence(rec, tmp_path)
    assert duration["value"] == 200
    assert "TPR intended duration conflicts with production MDP duration" in issues
    assert "measured trajectory duration conflicts with intended production MDP duration" in issues
    assert any("gmx dump" in p for p in duration["provenance"])


def test_offset_measured_span_not_end_time(tmp_path):
    rec = record(tmp_path, 25)
    rec.trajectory_inspection = TrajectoryInspection(inspected=True, start_time_ps=5000,
        end_time_ps=30000, total_duration_ps=25000, tpr_duration_ps=25000)
    duration, issues = duration_evidence(rec, tmp_path)
    assert duration["value"] == 25
    assert not issues


def test_tpr_timing_roundtrip():
    original = TrajectoryInspection(tpr_duration_ps=200000)
    assert TrajectoryInspection.from_dict(original.to_dict()).tpr_duration_ps == 200000


def test_inspector_failed_check_not_measured(tmp_path, monkeypatch):
    from analysis.campaign.trajectory import inspector
    from analysis.campaign.gmx import GmxResult
    rec = record(tmp_path)
    monkeypatch.setattr(inspector, "gmx_available", lambda _: True)
    monkeypatch.setattr(inspector, "run_gmx", lambda *a, **k: GmxResult(
        ["gmx", "check"], 1, "", "Reading frame 0 time 0\nLast frame 2500 time 25000\nCoords 2501 10"))
    result = inspector.inspect_trajectory(trajectory_paths=rec.production_trajectory_paths,
        fingerprint_segments=False)
    assert result.inspected is False
    assert result.total_duration_ps is None
    assert result.n_frames is None


def test_inspector_tpr_and_measured_timing(tmp_path, monkeypatch):
    from analysis.campaign.trajectory import inspector
    from analysis.campaign.gmx import GmxResult
    rec = record(tmp_path)
    calls = []
    def run(args, **kwargs):
        calls.append(args)
        text = "dt = .001\nnsteps = 200000000\n" if args[0] == "dump" else "Reading frame 0 time 5000\nLast frame 2500 time 30000\nCoords 2501 10"
        return GmxResult(args, 0, text, "")
    monkeypatch.setattr(inspector, "gmx_available", lambda _: True)
    monkeypatch.setattr(inspector, "run_gmx", run)
    result = inspector.inspect_trajectory(trajectory_paths=rec.production_trajectory_paths,
        topology_path=tmp_path / "md.tpr", fingerprint_segments=False)
    assert result.inspected is True
    assert result.tpr_duration_ps == 200000
    assert result.total_duration_ps == 25000
    assert calls[0][0] == "dump"


def test_inspector_success_without_timing_is_incomplete(tmp_path, monkeypatch):
    from analysis.campaign.trajectory import inspector
    from analysis.campaign.gmx import GmxResult
    rec = record(tmp_path)
    monkeypatch.setattr(inspector, "gmx_available", lambda _: True)
    monkeypatch.setattr(inspector, "run_gmx", lambda *a, **k: GmxResult(["gmx", "check"], 0, "", ""))
    result = inspector.inspect_trajectory(trajectory_paths=rec.production_trajectory_paths,
        fingerprint_segments=False)
    assert result.inspected is False
    assert result.total_duration_ps is None
    assert any(w.code == "gmx_check_incomplete" for w in result.warnings)
