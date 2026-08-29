"""Source-file fingerprinting for reproducible study discovery.

Two strategies, chosen by size:

* **fast**   – ``size`` + ``mtime_ns`` + sha256 of the first and last
  ``FAST_EDGE_BYTES``.  O(1) in file size; safe for multi-GB trajectories.
  Not collision-proof: two different trajectories with identical size, mtime
  and identical edges would collide (astronomically unlikely for real MD
  output, but the manifest records ``mode="fast"`` so a reader knows).
* **strong** – full-content sha256 (reuses :func:`runtime.artifacts.checksum`).

Rationale follows the existing SimForge convention (``runtime/artifacts.py``,
``core/provenance.py`` both use streaming sha256); we add the fast mode purely
to avoid hashing hundreds of GB of trajectory during discovery.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

from analysis.campaign.models import SourceFingerprint

FAST_EDGE_BYTES = 1 << 20          # 1 MiB from each end
STRONG_MAX_BYTES = 64 << 20        # <= 64 MiB: hash in full even without asking


def _edge_digest(path: Path, edge: int) -> str:
    h = hashlib.sha256()
    size = path.stat().st_size
    with path.open("rb") as fh:
        head = fh.read(min(edge, size))
        h.update(head)
        if size > edge:
            fh.seek(max(size - edge, edge))
            h.update(fh.read(edge))
    return h.hexdigest()


def _full_digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def fingerprint_file(path: Path, *, strong: bool = False) -> SourceFingerprint:
    """Return a :class:`SourceFingerprint` for ``path``.

    ``strong=True`` forces a full-content hash.  Otherwise files at or below
    ``STRONG_MAX_BYTES`` are hashed in full and larger ones use the fast mode.
    """
    p = Path(path)
    st = p.stat()
    size = st.st_size

    use_strong = strong or size <= STRONG_MAX_BYTES
    if use_strong:
        digest = _full_digest(p)
        return SourceFingerprint(
            path=str(p.resolve()), size_bytes=size, mtime_ns=st.st_mtime_ns,
            mode="strong", digest=digest,
            note="full-content sha256",
        )

    edge_digest = _edge_digest(p, FAST_EDGE_BYTES)
    combined = hashlib.sha256(
        f"{size}:{st.st_mtime_ns}:{edge_digest}".encode()
    ).hexdigest()
    return SourceFingerprint(
        path=str(p.resolve()), size_bytes=size, mtime_ns=st.st_mtime_ns,
        mode="fast", digest=combined, edge_bytes=FAST_EDGE_BYTES,
        note=(
            f"fast fingerprint (size + mtime + sha256 of first/last "
            f"{FAST_EDGE_BYTES} bytes); not collision-proof"
        ),
    )


def fingerprint_map(paths: dict[str, Path], *, strong: bool = False) -> dict[str, SourceFingerprint]:
    out: dict[str, SourceFingerprint] = {}
    for label, path in paths.items():
        if path is None:
            continue
        p = Path(path)
        if p.is_file():
            out[label] = fingerprint_file(p, strong=strong)
    return out
