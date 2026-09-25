"""Persistent, identity-addressed intermediate store.

Layout::

    <root>/<artifact_identity[:2]>/<artifact_identity>/
        manifest.json          written LAST — its presence marks a complete entry
        <artifact files>       content-hashed in the manifest
    <root>/.building-<identity>-<pid>-*/   work in progress (never read)
    <root>/.locks/<identity>.lock          one builder per identity (flock)
    <root>/.quarantine/…                   corrupted entries moved aside, never reused

Reuse requires the manifest *and* every file to match: schema, identities,
dependency identities and each file's sha256.  A published entry is
immutable: republishing the same identity with different content is a
reproducibility error (:class:`IntermediateConflict`), never a silent
replacement.  No garbage collection (entries accumulate by design for now).
"""
from __future__ import annotations

import fcntl
import json
import os
import shutil
import tempfile
import time
from contextlib import contextmanager
from dataclasses import replace
from pathlib import Path
from typing import Optional

from analysis.campaign.intermediates.base import (
    INTERMEDIATE_SCHEMA, IntermediateResult, IntermediateStatus,
)
from analysis.campaign.models import ResultArray, StorageRef

MANIFEST = "manifest.json"


class IntermediateConflict(RuntimeError):
    """Same claimed identity, different content — a reproducibility violation."""


def _sha256(path: Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class IntermediateStore:
    def __init__(self, root: str | Path):
        self.root = Path(root)

    def entry_dir(self, identity: str) -> Path:
        return self.root / identity[:2] / identity

    @contextmanager
    def lock(self, identity: str):
        d = self.root / ".locks"
        d.mkdir(parents=True, exist_ok=True)
        with open(d / f"{identity}.lock", "w") as fh:
            fcntl.flock(fh, fcntl.LOCK_EX)
            try:
                yield
            finally:
                fcntl.flock(fh, fcntl.LOCK_UN)

    def new_build_dir(self, identity: str) -> Path:
        self.root.mkdir(parents=True, exist_ok=True)
        return Path(tempfile.mkdtemp(prefix=f".building-{identity[:24]}-{os.getpid()}-",
                                     dir=self.root))

    # ── reading ────────────────────────────────────────────────────────────
    def load(self, identity: str, *, expect_dependencies: Optional[list[str]] = None
             ) -> tuple[Optional[IntermediateResult], str]:
        """``(result, "")`` for a valid entry, else ``(None, why)``."""
        entry = self.entry_dir(identity)
        mf = entry / MANIFEST
        if not mf.is_file():
            return None, "no cached entry"
        try:
            m = json.loads(mf.read_text())
        except (OSError, ValueError):
            return None, "manifest unreadable"
        if m.get("schema") != INTERMEDIATE_SCHEMA or m.get("status") != "complete":
            return None, "manifest schema/status not supported"
        if m.get("artifact_identity") != identity:
            return None, "manifest identity does not match its location"
        deps = [d["artifact_identity"] for d in m.get("dependencies", [])]
        if expect_dependencies is not None and deps != list(expect_dependencies):
            return None, "dependency identities differ"
        for name, sha in (m.get("files") or {}).items():
            p = entry / name
            if Path(name).name != name or not p.is_file():
                return None, f"artifact {name} missing"
            if _sha256(p) != sha:
                return None, f"artifact {name} changed after publication"
        return self._result(m, entry, IntermediateStatus.CACHED), ""

    def _result(self, m: dict, entry: Path, status: str) -> IntermediateResult:
        def absolute(ref: Optional[StorageRef]):
            return replace(ref, path=str(entry / ref.path)) if ref is not None else None
        arrays = []
        for a in m.get("arrays", []):
            arr = ResultArray.from_dict(a)
            arr.storage = absolute(arr.storage)
            for ax in arr.axes:
                ax.values_ref = absolute(ax.values_ref)
            arrays.append(arr)
        return IntermediateResult(
            intermediate_id=m["intermediate_id"], version=m["version"], request=m["request"],
            status=status, definition_identity=m["definition_identity"],
            input_identity=m["input_identity"], artifact_identity=m["artifact_identity"],
            definition_evidence=m.get("definition_evidence", {}), view_ref=m.get("view_ref"),
            arrays=arrays,
            artifacts=[absolute(StorageRef.from_dict(x)) for x in m.get("artifacts", [])],
            dependencies=m.get("dependencies", []), provenance=m.get("provenance", {}),
            entry=str(entry))

    # ── publishing ─────────────────────────────────────────────────────────
    def publish(self, build: Path, result: IntermediateResult) -> IntermediateResult:
        """Write the manifest last in ``build`` and move it into place atomically."""
        identity = result.artifact_identity

        def relative(ref: Optional[StorageRef]):
            if ref is None:
                return None
            p = Path(ref.path)
            if p.parent.resolve() != build.resolve():
                raise ValueError(f"artifact {p} is not inside the build directory")
            return replace(ref, path=p.name, attrs={**(ref.attrs or {}), "sha256": _sha256(p)})
        arrays = []
        for a in result.arrays:
            arr = ResultArray.from_dict(a.to_dict())
            arr.storage = relative(arr.storage)
            for ax in arr.axes:
                ax.values_ref = relative(ax.values_ref)
            arrays.append(arr.to_dict())
        artifacts = [relative(x).to_dict() for x in result.artifacts]
        files = {p.name: _sha256(p) for p in sorted(build.iterdir()) if p.is_file()}
        manifest = {"schema": INTERMEDIATE_SCHEMA, "status": "complete",
                    "intermediate_id": result.intermediate_id, "version": result.version,
                    "request": result.request, "artifact_identity": identity,
                    "definition_identity": result.definition_identity,
                    "input_identity": result.input_identity,
                    "definition_evidence": result.definition_evidence,
                    "view_ref": result.view_ref, "dependencies": result.dependencies,
                    "arrays": arrays, "artifacts": artifacts, "files": files,
                    "provenance": result.provenance}
        tmp = build / (MANIFEST + ".tmp")
        tmp.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        os.replace(tmp, build / MANIFEST)
        final = self.entry_dir(identity)
        final.parent.mkdir(parents=True, exist_ok=True)
        if final.exists():
            existing, why = self.load(identity)
            if existing is not None:
                old = json.loads((final / MANIFEST).read_text()).get("files")
                shutil.rmtree(build, ignore_errors=True)
                if old != files:
                    raise IntermediateConflict(
                        f"{result.intermediate_id} {identity[:12]}: recomputation produced "
                        f"different content under the same identity (non-reproducible)")
                return existing
            q = self.root / ".quarantine"
            q.mkdir(parents=True, exist_ok=True)
            os.rename(final, q / f"{identity}-{time.time_ns()}")    # corrupt: keep for inspection
            result.provenance["replaced_invalid_entry"] = why
        os.rename(build, final)
        out = self._result(manifest, final, result.status)
        out.provenance = result.provenance
        return out
