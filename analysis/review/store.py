"""Content-addressed, write-once store for the (small) files a session references.

Campaign output directories are mutable: a later recomputation rewrites
``observables/<id>/…`` in place.  A review session therefore freezes every
result / provenance / time-index file it references here, named by its
sha256, so what a session points to can never change underneath it.
Trajectories are never copied.
"""
from __future__ import annotations

import hashlib
import io
import os
import tempfile
from pathlib import Path

STORE_DIRNAME = "store"


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _place(store: Path, digest: str, suffix: str, write) -> Path:
    dest = Path(store) / digest[:2] / f"{digest}{suffix}"
    if dest.is_file() and sha256_file(dest) == digest:
        return dest                                   # write-once: already present and intact
    dest.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=dest.parent, prefix=".put-")
    try:
        with os.fdopen(fd, "wb") as f:
            write(f)
        if sha256_file(tmp) != digest:
            raise OSError(f"store copy of {digest[:12]} changed while being written")
        os.replace(tmp, dest)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)
    return dest


def put_file(store: str | Path, src: str | Path) -> tuple[Path, str]:
    """Freeze ``src``; returns ``(store path, sha256)``."""
    src = Path(src)
    digest = sha256_file(src)

    def write(f):
        with open(src, "rb") as s:
            for chunk in iter(lambda: s.read(1 << 20), b""):
                f.write(chunk)
    return _place(Path(store), digest, src.suffix, write), digest


def put_bytes(store: str | Path, data: bytes, suffix: str) -> tuple[Path, str]:
    digest = hashlib.sha256(data).hexdigest()
    return _place(Path(store), digest, suffix, lambda f: f.write(data)), digest


def put_array(store: str | Path, values, dtype) -> tuple[Path, str]:
    """A 1-D array as a deterministic ``.npy`` (no timestamps inside)."""
    import numpy as np
    buf = io.BytesIO()
    np.save(buf, np.asarray(values, dtype=dtype), allow_pickle=False)
    return put_bytes(store, buf.getvalue(), ".npy")


def verify(path: str | Path, digest: str) -> bool:
    p = Path(path)
    return p.is_file() and sha256_file(p) == digest
