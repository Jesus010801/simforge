from __future__ import annotations

import hashlib
import json
import threading
from pathlib import Path
from typing import Any

from runtime.artifacts import checksum


class StepCacheManager:
    """
    Content-addressed step cache.

    Fingerprint = SHA-256 of:
      - step.params (JSON-serialized, sorted keys)
      - checksums of all declared input files

    On cache hit the step is skipped entirely; on miss the fingerprint is
    stored so the next run can detect an unchanged step.

    Persisted to <workspace_dir>/metadata/step_cache.json.
    Thread-safe writes.
    """

    def __init__(self, workspace_dir: Path) -> None:
        meta = workspace_dir / "metadata"
        meta.mkdir(parents=True, exist_ok=True)
        self._path  = meta / "step_cache.json"
        self._lock  = threading.Lock()
        self._cache: dict[str, str] = {}   # step_id → fingerprint
        self._load()

    # ── Public API ────────────────────────────────────────────────────────────

    def fingerprint(self, params: dict[str, Any], input_paths: list[Path]) -> str:
        h = hashlib.sha256()
        h.update(json.dumps(params, sort_keys=True).encode())
        for p in sorted(str(ip) for ip in input_paths):
            h.update(p.encode())
            h.update(checksum(Path(p)).encode())
        return h.hexdigest()

    def is_cached(self, step_id: str, fingerprint: str) -> bool:
        with self._lock:
            return self._cache.get(step_id) == fingerprint

    def record(self, step_id: str, fingerprint: str) -> None:
        with self._lock:
            self._cache[step_id] = fingerprint
            self._persist()

    def invalidate(self, step_id: str) -> None:
        with self._lock:
            self._cache.pop(step_id, None)
            self._persist()

    def clear(self) -> None:
        with self._lock:
            self._cache.clear()
            self._persist()

    # ── Persistence ───────────────────────────────────────────────────────────

    def _persist(self) -> None:
        self._path.write_text(json.dumps(self._cache, indent=2), encoding="utf-8")

    def _load(self) -> None:
        if not self._path.exists():
            return
        try:
            self._cache = json.loads(self._path.read_text(encoding="utf-8"))
        except Exception:
            self._cache = {}
