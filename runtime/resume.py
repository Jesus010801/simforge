"""Explicit manifest-driven workflow resume support."""
from __future__ import annotations

import json
import shutil
from datetime import datetime, timezone
from pathlib import Path


class ResumeValidationError(RuntimeError):
    """Raised when an explicit resume request is unsafe."""


def _load_manifest(workspace: Path) -> dict:
    manifest_path = workspace / "metadata" / "execution_manifest.json"
    if not manifest_path.exists():
        raise ResumeValidationError(
            f"Cannot resume: execution manifest is missing: {manifest_path}"
        )
    return json.loads(manifest_path.read_text())


def _resolve_step(steps: list[dict], requested: str) -> int:
    matches = [
        index for index, entry in enumerate(steps)
        if requested in (entry["step_id"], entry["dir_name"])
    ]
    if not matches:
        available = ", ".join(
            f"{entry['step_id']} ({entry['dir_name']})" for entry in steps
        )
        raise ResumeValidationError(
            f"Invalid --from-step '{requested}'. Available steps: {available}"
        )
    return matches[0]


def _step_metadata(workspace: Path, entry: dict) -> tuple[Path, dict]:
    step_dir = workspace / "steps" / entry["dir_name"]
    metadata_path = step_dir / "metadata.json"
    if not step_dir.is_dir():
        raise ResumeValidationError(
            f"Cannot resume: step directory is missing: {step_dir}"
        )
    if not metadata_path.exists():
        raise ResumeValidationError(
            f"Cannot resume: step metadata is missing: {metadata_path}"
        )
    metadata = json.loads(metadata_path.read_text())
    return step_dir, metadata


def _append_journal(workspace: Path, event_type: str, step_id: str, **data) -> None:
    path = workspace / "metadata" / "execution_journal.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    event = {
        "event_type": event_type,
        "workspace_id": workspace.name,
        "step_id": step_id,
        "timestamp": datetime.now(timezone.utc).timestamp(),
        "severity": "INFO",
        "message": data.pop("message", ""),
        "data": data,
    }
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(event, ensure_ascii=False) + "\n")


def _remove_declared_outputs(step_dir: Path, metadata: dict) -> list[str]:
    removed: list[str] = []
    for relative in metadata.get("expected_outputs", []):
        output = (step_dir / relative).resolve()
        try:
            output.relative_to(step_dir.resolve())
        except ValueError as exc:
            raise ResumeValidationError(
                f"Refusing to clean output outside step directory: {output}"
            ) from exc
        if output.is_dir():
            shutil.rmtree(output)
            removed.append(str(output))
        elif output.exists() or output.is_symlink():
            output.unlink()
            removed.append(str(output))
    return removed


def prepare_resume(
    workspace: Path | str,
    from_step: str,
    *,
    clean_from: bool = False,
) -> dict:
    """Validate upstream artifacts, optionally clean downstream, and report."""
    workspace = Path(workspace).resolve()
    manifest = _load_manifest(workspace)
    steps = manifest.get("steps", [])
    selected_index = _resolve_step(steps, from_step)
    selected = steps[selected_index]
    prior = steps[:selected_index]
    downstream = steps[selected_index:]
    canonical_step = selected["step_id"]

    _append_journal(
        workspace,
        "RESUME_STARTED",
        canonical_step,
        message=f"Explicit resume requested from {selected['dir_name']}",
        requested_from_step=from_step,
        clean_from_enabled=clean_from,
    )

    reused_outputs: list[str] = []
    missing_outputs: list[str] = []
    for entry in prior:
        step_dir, metadata = _step_metadata(workspace, entry)
        for relative in metadata.get("expected_outputs", []):
            output = step_dir / relative
            if output.exists():
                reused_outputs.append(str(output.resolve()))
            else:
                missing_outputs.append(str(output.relative_to(workspace)))

    # Validate the selected step inputs before cleanup or execution.
    selected_dir, selected_metadata = _step_metadata(workspace, selected)
    for relative in selected_metadata.get("required_inputs", []):
        required_input = (selected_dir / relative).resolve()
        if required_input.exists():
            resolved = str(required_input)
            if resolved not in reused_outputs:
                reused_outputs.append(resolved)
        else:
            try:
                missing = str(required_input.relative_to(workspace))
            except ValueError:
                missing = str(required_input)
            if missing not in missing_outputs:
                missing_outputs.append(missing)

    report = {
        "resumed": not missing_outputs,
        "from_step": canonical_step,
        "from_step_dir": selected["dir_name"],
        "skipped_prior_steps": [entry["step_id"] for entry in prior],
        "reused_outputs": reused_outputs,
        "missing_outputs": missing_outputs,
        "clean_from_enabled": clean_from,
        "cleaned_outputs": [],
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    report_path = workspace / "metadata" / "resume_report.json"
    if missing_outputs:
        report_path.write_text(json.dumps(report, indent=2))
        details = "\n".join(f"- {path} is missing" for path in missing_outputs)
        raise ResumeValidationError(
            f"Cannot resume from {selected['dir_name']} because:\n{details}"
        )

    if clean_from:
        for entry in downstream:
            step_dir, metadata = _step_metadata(workspace, entry)
            report["cleaned_outputs"].extend(
                _remove_declared_outputs(step_dir, metadata)
            )

    report_path.write_text(json.dumps(report, indent=2))
    _append_journal(
        workspace,
        "RESUME_VALIDATED",
        canonical_step,
        message=f"Resume inputs validated for {selected['dir_name']}",
        reused_outputs=len(reused_outputs),
        clean_from_enabled=clean_from,
    )
    for entry in prior:
        _append_journal(
            workspace,
            "RESUME_SKIPPED_STEP",
            entry["step_id"],
            message="Step reused and skipped due to explicit resume",
            reason="skipped_due_to_resume",
        )
    return report
