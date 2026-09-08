"""Inspection failure paths must preserve read-only behavior."""
import subprocess

import pytest

from analysis.campaign import gmx
from analysis.campaign.orchestration.study_analyzer import run_inspect
from tests.analysis.campaign.test_legacy import system


def test_gmx_timeout_bytes_are_normalized(monkeypatch):
    monkeypatch.setattr(gmx, 'gmx_available', lambda binary: True)
    monkeypatch.setattr(gmx, 'gmx_version', lambda binary: 'test')
    def timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired(args[0], 1, output=b'partial\xff', stderr=b'failed')
    monkeypatch.setattr(subprocess, 'run', timeout)
    result = gmx.run_gmx(['check', '-f', 'missing.xtc'], timeout=1)
    assert result.returncode == 124
    assert isinstance(result.stdout, str)
    assert result.stderr == 'failed'


def test_inspection_export_refuses_dangling_symlink(tmp_path):
    system(tmp_path)
    output = tmp_path / 'export'
    output.mkdir()
    target = tmp_path / 'must_not_be_created.json'
    (output / 'study_manifest.json').symlink_to(target)
    with pytest.raises(FileExistsError):
        run_inspect(tmp_path, output_dir=output, inspect_trajectories=False)
    assert not target.exists()
    assert (output / 'study_manifest.json').is_symlink()
    assert not (output / 'study_manifest.yaml').exists()
