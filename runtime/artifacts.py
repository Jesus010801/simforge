from __future__ import annotations

import hashlib
import json
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal


SemanticRole = Literal[
    "topology",
    "coordinates",
    "trajectory",
    "checkpoint",
    "analysis",
    "parameter",
    "log",
    "other",
]


@dataclass
class ArtifactRef:
    path:          str          # absolute path
    checksum:      str          # SHA-256 hex
    semantic_role: SemanticRole
    step_id:       str          # which step produced it
    created_at:    float        = field(default_factory=time.time)
    size_bytes:    int          = 0

    def to_dict(self) -> dict:
        return {
            "path":          self.path,
            "checksum":      self.checksum,
            "semantic_role": self.semantic_role,
            "step_id":       self.step_id,
            "created_at":    self.created_at,
            "size_bytes":    self.size_bytes,
        }

    @classmethod
    def from_dict(cls, d: dict) -> ArtifactRef:
        return cls(**d)


@dataclass
class ArtifactModification:
    step_id:        str
    checksum_after: str
    timestamp:      float = field(default_factory=time.time)
    description:    str   = ""

    def to_dict(self) -> dict:
        return {
            "step_id":        self.step_id,
            "checksum_after": self.checksum_after,
            "timestamp":      self.timestamp,
            "description":    self.description,
        }


@dataclass
class ArtifactLineage:
    ref:           ArtifactRef
    modifications: list[ArtifactModification] = field(default_factory=list)

    def record_modification(self, step_id: str, path: Path, description: str = "") -> None:
        self.modifications.append(ArtifactModification(
            step_id        = step_id,
            checksum_after = checksum(path),
            description    = description,
        ))

    def to_dict(self) -> dict:
        return {
            "ref":           self.ref.to_dict(),
            "modifications": [m.to_dict() for m in self.modifications],
        }


def checksum(path: Path) -> str:
    h = hashlib.sha256()
    try:
        with path.open("rb") as fh:
            for chunk in iter(lambda: fh.read(65536), b""):
                h.update(chunk)
    except FileNotFoundError:
        return ""
    return h.hexdigest()


def make_ref(path: Path, semantic_role: SemanticRole, step_id: str) -> ArtifactRef:
    try:
        size = path.stat().st_size
    except FileNotFoundError:
        size = 0
    return ArtifactRef(
        path          = str(path.resolve()),
        checksum      = checksum(path),
        semantic_role = semantic_role,
        step_id       = step_id,
        size_bytes    = size,
    )


class ArtifactRegistry:
    """
    Tracks file provenance across all steps.

    Persisted to <workspace_dir>/metadata/artifact_registry.json.
    Thread-safe writes.
    """

    def __init__(self, workspace_dir: Path) -> None:
        meta = workspace_dir / "metadata"
        meta.mkdir(parents=True, exist_ok=True)
        self._path = meta / "artifact_registry.json"
        self._lock = threading.Lock()
        self._lineages: dict[str, ArtifactLineage] = {}
        self._load()

    def register(
        self,
        path:          Path,
        semantic_role: SemanticRole,
        step_id:       str,
    ) -> ArtifactRef:
        ref = make_ref(path, semantic_role, step_id)
        with self._lock:
            self._lineages[str(path.resolve())] = ArtifactLineage(ref=ref)
            self._persist()
        return ref

    def record_modification(
        self,
        path:        Path,
        step_id:     str,
        description: str = "",
    ) -> None:
        key = str(path.resolve())
        with self._lock:
            if key not in self._lineages:
                return
            self._lineages[key].record_modification(step_id, path, description)
            self._persist()

    def lineage(self, path: Path) -> ArtifactLineage | None:
        return self._lineages.get(str(path.resolve()))

    def all_lineages(self) -> dict[str, ArtifactLineage]:
        return dict(self._lineages)

    # ── Persistence ───────────────────────────────────────────────────────────

    def _persist(self) -> None:
        data = {k: v.to_dict() for k, v in self._lineages.items()}
        self._path.write_text(json.dumps(data, indent=2), encoding="utf-8")

    def _load(self) -> None:
        if not self._path.exists():
            return
        try:
            data = json.loads(self._path.read_text(encoding="utf-8"))
            for k, v in data.items():
                ref  = ArtifactRef.from_dict(v["ref"])
                mods = [ArtifactModification(**m) for m in v.get("modifications", [])]
                self._lineages[k] = ArtifactLineage(ref=ref, modifications=mods)
        except Exception:
            pass
