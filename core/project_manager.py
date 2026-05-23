"""
core/project_manager.py — Project lifecycle management.

Each compile creates a new immutable timestamped run under the project directory.
Previous runs are never deleted or overwritten.

Structure:
  {output_dir}/{project_name}/
    project.json                      ← run registry
    runs/
      2026-05-23_14-30-00/            ← one per compile
        metadata/execution_manifest.json
        steps/
        inputs/
        ...
"""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path


class ProjectManager:
    """Manages project lifecycle — run creation, provenance, run history."""

    @staticmethod
    def project_dir(output_dir: str | Path, project_name: str) -> Path:
        return Path(output_dir) / project_name

    @staticmethod
    def create_run_dir(project_dir: Path, timestamp: str | None = None) -> Path:
        """
        Create a fresh timestamped run directory under project_dir/runs/.
        Returns the new (empty) run path.
        """
        ts = timestamp or datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        run_dir = project_dir / "runs" / ts
        run_dir.mkdir(parents=True, exist_ok=True)
        return run_dir

    @staticmethod
    def update_project_registry(project_dir: Path, run_info: dict) -> None:
        """Append run_info to project.json, creating the file if absent."""
        registry_path = project_dir / "project.json"
        if registry_path.exists():
            registry = json.loads(registry_path.read_text())
        else:
            registry = {"runs": []}
        registry["runs"].append(run_info)
        registry_path.write_text(json.dumps(registry, indent=4))

    @staticmethod
    def get_run_history(project_dir: Path) -> list[dict]:
        """Return all run entries from project.json, newest first."""
        registry_path = project_dir / "project.json"
        if not registry_path.exists():
            return []
        registry = json.loads(registry_path.read_text())
        return list(reversed(registry.get("runs", [])))
