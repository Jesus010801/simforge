"""Canonical scientific identities for discovered systems.

A canonical id is:

* **deterministic**   – a pure function of resolved semantic identity;
* **filesystem-safe** – ``[a-z0-9_]`` only, lowercase;
* **human-readable**  – ``glp1r__semaglutide__tip3p__rep02`` when identity is known;
* **filename-independent** – arbitrary trajectory/topology names never appear;
* **stable across rediscovery** and **independent of directory traversal order**.

Unknown fields use the literal token ``unknown`` (never invented).  Collisions
between two systems that resolve to the same canonical id are broken with a
short suffix derived from a *stable* fingerprint of the system's source files —
again never traversal order.
"""
from __future__ import annotations

import hashlib
import re

UNKNOWN = "unknown"
_DROP_HYPHEN_RE = re.compile(r"(?<=[a-z0-9])-(?=[a-z0-9])")  # join "GLP-1R" -> "glp1r"
_TO_USCORE_RE = re.compile(r"[^a-z0-9]+")


def sanitize_token(value: str | None) -> str:
    """Lowercase; join intra-word hyphens; collapse other runs to ``_``.

    ``"GLP-1R" -> "glp1r"``, ``"Exendin-4" -> "exendin4"``,
    ``"TIP3P" -> "tip3p"``, ``"water model" -> "water_model"``.
    """
    if not value:
        return UNKNOWN
    s = value.strip().lower()
    s = _DROP_HYPHEN_RE.sub("", s)
    s = _TO_USCORE_RE.sub("_", s).strip("_")
    return s or UNKNOWN


def format_replicate(replicate_id: str | int | None) -> str:
    """Normalise a replicate identifier to ``repNN`` when numeric, else a token."""
    if replicate_id is None:
        return f"rep_{UNKNOWN}"
    s = str(replicate_id).strip().lower()
    m = re.search(r"(\d+)", s)
    if m:
        return f"rep{int(m.group(1)):02d}"
    tok = sanitize_token(s)
    return tok if tok.startswith("rep") else f"rep_{tok}"


def canonical_system_id(
    *,
    receptor: str | None,
    partner: str | None,
    water_model: str | None,
    replicate_id: str | int | None,
    extra_dimensions: dict[str, str] | None = None,
) -> str:
    """Build the canonical id from resolved identity fields.

    ``extra_dimensions`` (e.g. ``{"temperature": "310K"}``) are appended in
    sorted-key order so additional study axes never reorder the core fields.
    """
    parts = [
        sanitize_token(receptor),
        sanitize_token(partner),
        sanitize_token(water_model),
    ]
    for key in sorted((extra_dimensions or {})):
        val = sanitize_token(extra_dimensions[key])
        if val != UNKNOWN:
            parts.append(f"{sanitize_token(key)}_{val}")
    parts.append(format_replicate(replicate_id))
    return "__".join(parts)


def stable_disambiguator(tokens: list[str], length: int = 6) -> str:
    """Deterministic short suffix from content tokens (sorted, order-free)."""
    joined = "\x00".join(sorted(t for t in tokens if t))
    return hashlib.sha256(joined.encode()).hexdigest()[:length]


def deduplicate_ids(pairs: list[tuple[str, list[str]]]) -> dict[int, str]:
    """Resolve canonical-id collisions deterministically.

    ``pairs`` is ``[(canonical_id, disambiguation_tokens), ...]`` in *any*
    order.  Returns ``{original_index: final_id}``.  Colliding ids get a
    ``__<hash>`` suffix computed from their own tokens; if even that collides
    (identical tokens), a numeric ``__dupNN`` in fingerprint-sorted order is
    appended so the result is still deterministic.
    """
    # group indices by base id
    groups: dict[str, list[int]] = {}
    for idx, (cid, _tokens) in enumerate(pairs):
        groups.setdefault(cid, []).append(idx)

    final: dict[int, str] = {}
    for cid, idxs in groups.items():
        if len(idxs) == 1:
            final[idxs[0]] = cid
            continue
        # deterministic order within the collision group: by disambiguator then tokens
        annotated = sorted(
            idxs,
            key=lambda i: (stable_disambiguator(pairs[i][1]), tuple(sorted(pairs[i][1]))),
        )
        seen: dict[str, int] = {}
        for i in annotated:
            suffix = stable_disambiguator(pairs[i][1])
            candidate = f"{cid}__{suffix}"
            n = seen.get(candidate, 0)
            seen[candidate] = n + 1
            if n:
                candidate = f"{candidate}__dup{n:02d}"
            final[i] = candidate
    return final
