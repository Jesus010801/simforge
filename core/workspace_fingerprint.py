# core/workspace_fingerprint.py
"""
Workspace fingerprinting and staleness detection.

A workspace is "fresh" when its build_signature matches what would be produced
by compiling the same YAML with the current builder source files.

If either the YAML or any builder changes, the signature diverges and SimForge
refuses to execute the stale workspace — preventing silent scientific errors.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from datetime import datetime

SIMFORGE_VERSION = "0.1.0"

# Builder files that, when changed, invalidate existing workspaces
_BUILDER_FILES: list[str] = [
    "builders/step_builders/_utils.py",
    "builders/step_builders/assembly_builder.py",
    "builders/step_builders/minimization_builder.py",
    "builders/step_builders/equilibration_builder.py",
    "builders/step_builders/production_builder.py",
    "builders/step_builders/preparation_builder.py",
    "builders/step_builders/analysis_builder.py",
    "builders/workspace_builder.py",
    "builders/builder_registry.py",
]

_PROJECT_ROOT = Path(__file__).parent.parent


def compute_builder_signature() -> str:
    """
    SHA256 of all builder source files (sorted, deterministic).

    Changes whenever any builder template or utility is modified.
    Used to detect stale workspaces without needing the original YAML.
    """
    h = hashlib.sha256()
    h.update(SIMFORGE_VERSION.encode())
    for rel in _BUILDER_FILES:
        p = _PROJECT_ROOT / rel
        if p.exists():
            h.update(rel.encode())
            h.update(p.read_bytes())
    return h.hexdigest()[:24]


def compute_build_signature(yaml_source: str) -> str:
    """
    Full build signature: YAML content + builder sources.

    Stored in the manifest at compile time.
    """
    h = hashlib.sha256()
    try:
        h.update(Path(yaml_source).read_bytes())
    except Exception:
        h.update(yaml_source.encode())
    h.update(compute_builder_signature().encode())
    return h.hexdigest()[:24]


def compute_template_hash(scripts: list[Path]) -> str:
    """Hash of generated script files for per-step provenance."""
    h = hashlib.sha256()
    for p in sorted(scripts):
        if p.exists():
            h.update(p.name.encode())
            h.update(p.read_bytes())
    return h.hexdigest()[:16]


def check_workspace_freshness(manifest_path: Path) -> tuple[bool, str]:
    """
    Compare the workspace's stored builder_signature to the current code.

    Returns:
        (is_fresh, message)
        is_fresh=True  → workspace is up-to-date, safe to execute
        is_fresh=False → workspace was compiled with different builders
    """
    if not manifest_path.exists():
        return False, "Manifest not found — workspace may be incomplete. Recompile required."

    try:
        manifest = json.loads(manifest_path.read_text())
    except Exception as exc:
        return False, f"Cannot read manifest: {exc}"

    stored = manifest.get("builder_signature")
    if stored is None:
        # Legacy workspace without fingerprint — warn but allow
        return True, "Legacy workspace (no fingerprint). Cannot verify freshness."

    current = compute_builder_signature()
    if stored != current:
        compiled_at = manifest.get("compiled_at", "unknown time")
        return False, (
            f"Workspace was compiled with different builders (compiled at {compiled_at}).\n"
            f"  stored  builder_signature: {stored}\n"
            f"  current builder_signature: {current}\n"
            f"Run: simforge recompile <yaml>"
        )

    return True, "OK"
