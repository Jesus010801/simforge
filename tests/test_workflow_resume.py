from __future__ import annotations

import json
from pathlib import Path

import pytest

from executors.execution_state import StepStatus
from runtime.executor import RuntimeExecutor
from runtime.journal import JournalWriter
from runtime.resume import ResumeValidationError, prepare_resume


STEPS = [
    ("prepare", "01_prepare"),
    ("build", "02_build"),
    ("finish", "03_finish"),
]


def _workspace(tmp_path: Path) -> Path:
    workspace = tmp_path / "run"
    (workspace / "metadata").mkdir(parents=True)
    manifest_steps = []
    for index, (step_id, dir_name) in enumerate(STEPS):
        step_dir = workspace / "steps" / dir_name
        step_dir.mkdir(parents=True)
        required = [] if index == 0 else [f"../{STEPS[index - 1][1]}/{STEPS[index - 1][0]}.out"]
        metadata = {
            "step_id": step_id,
            "required_inputs": required,
            "expected_outputs": [f"{step_id}.out"],
            "params": {},
            "blocking": True,
        }
        (step_dir / "metadata.json").write_text(json.dumps(metadata))
        (step_dir / "run.sh").write_text(
            f'printf "%s\\n" "{step_id}" >> executed.log\n'
            f'touch "{step_id}.out"\n'
        )
        manifest_steps.append(
            {
                "step_id": step_id,
                "dir_name": dir_name,
                "depends_on": [] if index == 0 else [STEPS[index - 1][0]],
            }
        )
    (workspace / "metadata" / "execution_manifest.json").write_text(
        json.dumps({"system_type": "test", "n_steps": 3, "steps": manifest_steps})
    )
    return workspace


def _output(workspace: Path, index: int) -> Path:
    step_id, dir_name = STEPS[index]
    return workspace / "steps" / dir_name / f"{step_id}.out"


def test_resume_valid_step_skips_prior_and_executes_downstream(tmp_path: Path) -> None:
    workspace = _workspace(tmp_path)
    _output(workspace, 0).touch()

    report = prepare_resume(workspace, "02_build")
    state = RuntimeExecutor(
        workspace,
        dry_run=False,
        resume_prior_steps=set(report["skipped_prior_steps"]),
    ).run()

    prior = state.steps[0]
    assert prior.status == StepStatus.DONE
    assert prior.execution_reason == "skipped_due_to_resume"
    assert not (workspace / "steps" / "01_prepare" / "executed.log").exists()
    assert (workspace / "steps" / "02_build" / "executed.log").exists()
    assert (workspace / "steps" / "03_finish" / "executed.log").exists()


def test_missing_upstream_output_blocks_resume(tmp_path: Path) -> None:
    workspace = _workspace(tmp_path)

    with pytest.raises(ResumeValidationError, match=r"steps/01_prepare/prepare\.out is missing"):
        prepare_resume(workspace, "02_build")

    report = json.loads((workspace / "metadata" / "resume_report.json").read_text())
    assert report["resumed"] is False
    assert report["missing_outputs"] == ["steps/01_prepare/prepare.out"]


def test_invalid_step_id_has_clear_error(tmp_path: Path) -> None:
    workspace = _workspace(tmp_path)
    with pytest.raises(ResumeValidationError, match=r"Invalid --from-step 'unknown'"):
        prepare_resume(workspace, "unknown")


def test_clean_from_removes_only_selected_and_downstream_outputs(tmp_path: Path) -> None:
    workspace = _workspace(tmp_path)
    for index in range(3):
        _output(workspace, index).touch()

    report = prepare_resume(workspace, "02_build", clean_from=True)

    assert _output(workspace, 0).exists()
    assert not _output(workspace, 1).exists()
    assert not _output(workspace, 2).exists()
    assert report["clean_from_enabled"] is True
    assert len(report["cleaned_outputs"]) == 2


def test_resume_report_and_journal_record_provenance(tmp_path: Path) -> None:
    workspace = _workspace(tmp_path)
    _output(workspace, 0).touch()

    prepare_resume(workspace, "02_build")

    report = json.loads((workspace / "metadata" / "resume_report.json").read_text())
    assert report["resumed"] is True
    assert report["from_step"] == "build"
    assert report["skipped_prior_steps"] == ["prepare"]
    assert report["missing_outputs"] == []
    assert report["timestamp"]

    events = JournalWriter(workspace).read_all()
    event_types = [event.event_type.value for event in events]
    assert event_types == [
        "RESUME_STARTED",
        "RESUME_VALIDATED",
        "RESUME_SKIPPED_STEP",
    ]


def test_selected_step_required_input_is_validated(tmp_path: Path) -> None:
    workspace = _workspace(tmp_path)
    # Remove the prior expected output declaration to prove selected-step
    # preflight inputs are also part of early resume validation.
    metadata_path = workspace / "steps" / "01_prepare" / "metadata.json"
    metadata = json.loads(metadata_path.read_text())
    metadata["expected_outputs"] = []
    metadata_path.write_text(json.dumps(metadata))

    with pytest.raises(ResumeValidationError, match=r"steps/01_prepare/prepare\.out is missing"):
        prepare_resume(workspace, "02_build")
