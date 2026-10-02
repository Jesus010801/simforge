# ADR-012 — Artifact Registry and Content-Addressed Storage

**Status:** Accepted

## Context

Phase 1 defines immutable `ArtifactRef` values containing a supplied internal
artifact ID, SHA-256 digest, and optional media type and byte size. It deliberately
does not store or verify content. ADR-004 requires SHA-256 for physical artifacts
used by a KnowledgeBundle. The Phase 2 entity registry is descriptor-only and
explicitly excludes artifact/CAS persistence. Phase 3 in the implementation plan
calls for `ArtifactRef` semantics, SHA-256 verification, and in-memory/filesystem
registries.

The system needs a byte-storage contract that does not conflate artifact IDs,
content digests, source accessions, filenames, or scientific equivalence, and
that remains independently testable across volatile and filesystem backends.

## Decision

1. **Frozen dependencies.** Phase 1 `ArtifactRef`, `EntityId`, canonical
   serialization, the Phase 2 registry API, and Phase 2 errors remain unchanged.
   The artifact registry composes with Phase 1 values and does not depend on or
   mutate `EntityRegistry`.
2. **Three distinct identities.** Artifact IDs are caller-provided and are not
   derived from content. Content identity is the SHA-256 address of exact bytes.
   ArtifactRef equality is equality of the complete canonical Phase 1 reference.
   Content identity, ArtifactRef equality, and artifact-ID binding identity are
   distinct. Metadata is excluded from byte identity but included in exact
   ArtifactRef/binding equality.
3. **Permanent exact binding.** First successful registration binds one ID in a
   logical registry permanently to one exact complete `ArtifactRef`. Exact replay
   requires the same ID and every reference field (`sha256`, `size_bytes`,
   `media_type`, and any other Phase 1 field), with supplied bytes satisfying
   the reference; it succeeds idempotently and returns the existing binding.
   Any different reference field is `BINDING_CONFLICT`. There is no metadata
   enrichment/update operation or artifact revision. Different IDs may share
   identical bytes and one blob.
4. **Digest and collision contract.** SHA-256 of exact unmodified bytes is the
   only Phase 3 content address. The operational assumption is that equal
   digests represent content identity for addressing. If a digest location is
   occupied and candidate bytes are available, compare them: equal bytes reuse
   the blob; unequal bytes raise `INTEGRITY_FAILURE` with reason `HASH_COLLISION`
   and never overwrite the existing blob. This is not scientific equivalence.
5. **Bytes-only API.** Registration accepts exact `bytes` only; paths, streams,
   iterators/chunks, URLs, and mutable buffers are rejected. `register(reference,
   content)` validates digest and declared size and returns the committed
   reference. `get_ref(artifact_id)` retrieves the binding. `read(artifact_id)`
   returns immutable owned bytes only after validating binding, blob existence,
   declared size when present, and SHA-256. `verify(artifact_id)` performs the
   same checks and returns the bound `ArtifactRef` on success. Future streaming
   must preserve these digest, size, binding, retry, and atomic visibility
   semantics regardless of chunking.
6. **Immutable content without revisions.** There is no update, delete, artifact
   revision, mutable head, registry revision counter, or automatic garbage
   collection. Different reference under one ID conflicts; another immutable
   artifact requires another caller ID. Failed writes may leave invisible orphans.
7. **Filesystem namespaces and metadata.** Under one configured root, immutable
   blobs use deterministic canonical lowercase SHA-256 paths and bindings use a
   safe deterministic encoding/hash of artifact ID; raw IDs are never paths.
   Mappings are root-confined and temporary files are private under controlled
   same-filesystem directories. Binding metadata is exactly
   `ArtifactRef.to_canonical_bytes()` (`sski-json-v1`), persisted unchanged and
   revalidated as Phase 1 `ArtifactRef`; corrupt metadata is an integrity
   failure. Optional bookkeeping is not scientific provenance.
8. **Filesystem transaction.** Validate input, compute/check digest and size,
   stage complete bytes, compare an occupied digest by byte equality, prepare
   canonical binding metadata, acquire cross-process root coordination, re-check
   both the ID binding and digest location, publish a missing blob without
   replacement, then atomically publish the binding. Binding publication is
   the logical commit point. Before
   it, the registration is invisible; temporary/orphan blobs may remain and
   retry is safe. After it, the registration is committed and retries use
   exact-replay/conflict rules. Readers never observe partial content or binding.
9. **Concurrency and supported filesystems.** All handles/processes sharing a
   root coordinate binding publication; an in-process lock alone is insufficient.
   Same-ID identical concurrent writes converge to one binding and succeed;
   same-ID conflicts have at most one winner; different IDs for equal bytes may
   both succeed. Readers see pre-commit absence or complete post-commit state.
   Only local filesystems with required same-filesystem atomic publication and
   cross-process locking are supported. Arbitrary network filesystems and
   distributed stores are outside the guarantee.
10. **Persistence boundary.** Operation atomicity differs from restart
    persistence, process-crash consistency, and power-loss durability. Committed
    files that remain intact can be reconstructed after restart; pre-commit
    temporaries/orphans are not registrations. Phase 3 does not promise survival
    across arbitrary power loss, filesystem corruption, controller-cache loss,
    or hardware failure, nor stronger `fsync` durability.
11. **Digest normalization.** Phase 1 accepts exactly 64 hex characters,
    accepts uppercase and canonicalizes lowercase, and rejects `sha256:`.
    Phase 3 preserves this without changing Phase 1. Hash exact bytes only; no
    text, Unicode, newline, path, or metadata normalization participates.
12. **Structured errors.** A Phase 3-specific base exposes stable category and
    reason, relevant artifact ID, and field path. Binding conflicts expose the
    existing and requested references; integrity errors expose applicable
    expected/actual digest or size; storage I/O errors identify the operation and
    chain the underlying cause. Public categories/classes are
    `INVALID_ARGUMENT` (`InvalidArtifactOperationError`), `NOT_FOUND`
    (`ArtifactNotFoundError`), `BINDING_CONFLICT`
    (`ArtifactBindingConflictError`), `INTEGRITY_FAILURE`
    (`ArtifactIntegrityError`), and `STORAGE_IO_FAILURE`
    (`ArtifactStorageIOError`). Integrity reasons cover digest/size mismatch,
    corrupt binding/blob, missing blob, and hash collision. Malformed syntax is
    invalid input; a well-formed declaration that disagrees with bytes is an
    integrity failure. I/O causes are chained; programming defects are not
    translated, and Python argument-binding misuse may remain `TypeError`.
13. **Layer exclusions.** Provenance, authentication/authorization, licensing,
    federation, import, evidence persistence, resolution, policies, bundle
    compilation, physical KnowledgeBundle, BundleReader, and runtime access
    remain outside this registry.

## Consequences

The in-memory and filesystem backends can share one behavioral conformance suite.
The filesystem backend may leave unreachable temporary/orphan bytes after an
interruption, but must never expose a reference before verified content is
available. Existing references become retrievable only after successful
registration in a registry containing their bytes. Registry presence does not
prove who supplied content or whether it is scientifically valid.

This ADR was accepted following final design review. Acceptance does not authorize implementation; implementation requires separate user authorization against the accepted Phase 3 design and exact allowlist.
