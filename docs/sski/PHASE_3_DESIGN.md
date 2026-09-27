# Phase 3 design — Artifact Registry and content-addressed storage

**Status:** Proposed; design-only, pending adversarial review and ADR acceptance.
**Scope source:** Phase 3 in `09_IMPLEMENTATION_PLAN.md`: “Artifact Registry/content addressing” and “ArtifactRef semantics + SHA256 verification + in-memory/filesystem registry.”
**Implementation authorization:** None. A later explicit implementation authorization is required.

## 1. Purpose

Phase 1 defines an immutable `ArtifactRef` containing a caller-supplied `EntityId`,
a SHA-256 digest, and optional media type and size. It deliberately provides no
byte availability, filesystem access, or registry. Phase 3 designs the missing
artifact storage boundary: register bytes against an immutable reference, retain
and retrieve those bytes by verified content address, and provide volatile-memory
and local-filesystem reference backends.

This phase makes artifact availability and byte integrity testable before later
source federation or KnowledgeBundle compilation. A matching digest establishes
byte equality under SHA-256; it does not establish scientific truth, authorship,
source authenticity, licensing, or scientific equivalence.

## 2. Scope

The proposed artifact registry owns immutable byte content and a registry-local
binding from the existing artifact ID to one `ArtifactRef`. It has two backend
forms: an in-memory reference backend and a local filesystem backend. Both expose
the same observable contract. The filesystem backend is a local content store,
not an object-store service or a general-purpose persistence framework.

The phase specifies registration, reference lookup, verified byte reads, integrity
checks, immutability, retry behavior, backend boundaries, and conformance tests.
It does not compile or validate a complete workflow knowledge context.

## 3. Explicit non-goals

- Source networking, downloads, federation, source caches, adapters, importers,
  release discovery, or source-specific parsing.
- Claims/evidence persistence, scientific provenance graphs, identity resolution,
  equivalence inference, or artifact-to-entity scientific assertions.
- Artifact transformation, decompression, format parsing, MIME sniffing, chemical
  or structural validation, or derived-artifact generation.
- Artifact deletion, replacement, garbage collection, retention policy, or
  multi-location replication.
- PostgreSQL, SQL schemas/migrations, object stores, distributed services, or a
  general repository/unit-of-work framework.
- KnowledgeSnapshot, resolver, policies, compiler, physical KnowledgeBundle,
  BundleReader, runtime/HPC integration, SimulationPlan, or RunDerivationLedger.
- Changes to Phase 1 or Phase 2 production code, tests, exports, or contracts.

## 4. Relationship to Phase 1

Phase 1 remains the only owner of scientific/domain semantics and canonical
serialization. Phase 3 composes with the exact frozen `ArtifactRef`, `EntityId`,
and `sski-json-v1` representation. It does not add provenance or storage fields
to those values. A reference is content-pinned but is not proof that its bytes
exist; Phase 3's registry supplies that separate availability guarantee.

The 15 Phase 1 public objects do not need changes. `ArtifactRef.artifact_id` is
the supplied opaque internal ID, never an external accession and never derived
from content. `sha256` is the digest of exact stored bytes. Optional `size_bytes`,
when present, must equal the actual byte length. Optional `media_type` remains
descriptive metadata; storage does not infer or validate scientific format.

## 5. Relationship to Phase 2

Phase 2 remains a descriptor-only entity registry. It does not own artifact bytes,
and the artifact registry does not add artifact descriptors to, or alter, the
Phase 2 registry. Phase 3 may use the `EntityId` carried by `ArtifactRef`, but it
does not call `EntityRegistry`, modify its revision histories, use its operation
journal, or change its error hierarchy.

The composition is additive: a later compiler can combine already-resolved
Phase 1/2 references with verified artifact bytes. That compiler is not part of
Phase 3.

## 6. Data and semantic model

The proposal separates three concepts:

1. **Artifact ID:** the caller-supplied logical identity from `ArtifactRef`.
2. **Content address:** SHA-256 of the exact byte sequence.
3. **Stored object:** immutable bytes retrievable only after the reference and
   digest are checked.

The registry metadata is a binding `artifact_id -> ArtifactRef`. The byte object
is addressed by digest, so equal byte sequences can be physically deduplicated
even when references use different IDs. An artifact ID is not a path, filename,
source accession, or content hash. A registry does not infer that equal bytes
represent scientifically equivalent artifacts.

**Accepted-for-proposal identity rule:** artifact ID is caller-provided and is
not derived from the digest. The first successful registration binds that ID
permanently to one exact, complete immutable `ArtifactRef`. Exact replay requires
the same ID, SHA-256, `size_bytes`, `media_type`, every other `ArtifactRef`
field, and matching bytes; it succeeds idempotently and returns the existing
binding. Any difference in the reference is a binding conflict. Different IDs
may bind to identical bytes and share one digest-addressed blob. There is no
metadata update/enrichment operation and no artifact revision.

Keep these distinct: **content identity** is the SHA-256 address of exact bytes;
**ArtifactRef equality** compares the complete canonical Phase 1 reference;
**artifact-ID binding identity** is the permanent registry-local association
between the caller ID and that complete reference. Optional metadata is excluded
from byte identity but participates in exact binding equality. A matching digest
is not scientific equivalence.

## 7. Proposed public API

Proposed module ownership and API (not accepted or implemented):

| Module | Public objects |
|---|---|
| `simforge.knowledge.storage.artifact_registry` | `ArtifactRegistry`, `ArtifactRegistryError`, `InvalidArtifactOperationError`, `ArtifactNotFoundError`, `ArtifactBindingConflictError`, `ArtifactIntegrityError`, `ArtifactStorageIOError` |
| `simforge.knowledge.storage.memory_artifact_registry` | `InMemoryArtifactRegistry` |
| `simforge.knowledge.storage.filesystem_artifact_registry` | `FilesystemArtifactRegistry` |

Proposed operations:

```text
register(reference: ArtifactRef, content: bytes) -> ArtifactRef
get_ref(artifact_id: EntityId) -> ArtifactRef
read(artifact_id: EntityId) -> bytes
verify(artifact_id: EntityId) -> ArtifactRef
```

`register` accepts exact `bytes` only, computes SHA-256 over those exact bytes,
and validates the reference digest and optional declared size. It returns the
existing or newly committed binding. `get_ref` looks up the binding by ID.
`read` returns immutable owned bytes only after validating the stored binding,
blob existence, declared size when present, and SHA-256. `verify` performs the
same checks without returning content and returns the bound `ArtifactRef` on
success. Missing identity, binding conflict, invalid input, integrity failure,
and storage I/O failure have distinct structured errors. No public delete,
update, path query, content-hash helper, pagination, or generic query API is
proposed.

Bytes-only input is a Phase 3 decision: mutable buffers, text, paths, open files,
streams, iterators/chunks, and URLs are rejected rather than coerced. Future
streaming APIs must preserve the same digest, byte-size, ID-binding, retry, and
atomic-visibility semantics regardless of chunking. No public API exposes paths.

## 8. Identity rules

Artifact IDs are valid Phase 1 `EntityId` values supplied by the caller. Phase 3
does not generate IDs. External identifiers never participate in artifact lookup
or storage paths. Within one logical registry, `get_ref` is exact by Artifact ID.
Byte storage and integrity lookup use SHA-256. Multiple artifact IDs may bind to
the same digest and share the immutable stored bytes.

The one-ID/one-reference binding is registry-local; independent registries may
reuse an ID. Cross-registry ID coordination belongs to a later Global Knowledge
Plane service contract. A complete canonical `ArtifactRef`, including every
Phase 1 field, participates in binding equality; optional metadata does not
change the byte content address but a metadata difference under an existing ID
is a binding conflict.

## 9. Revision and versioning semantics

Artifact content and registry bindings are immutable. The registry has no
revision counter, mutable head, update, or history API. Any different reference,
including one differing only in metadata, conflicts under the bound ID; a
different immutable artifact requires a different caller-provided ID. Existing
`ArtifactRef` values remain usable as long as the registry retains their binding
and bytes.

Storage format versions, if needed for filesystem metadata, are storage-format
versions, not scientific entity or artifact revisions. They must not alter the
meaning or canonical bytes of Phase 1 `ArtifactRef`.

## 10. Provenance semantics

The registry records byte identity and optional Phase 1 media/size metadata only.
It does not assert who produced the bytes, which source/release supplied them,
which importer parsed them, or whether a claim is true. Phase 1 `ArtifactRef` has
no provenance field and remains unchanged. Later evidence/importer layers must
retain source attribution separately and may refer to the immutable artifact
reference. Hash verification is integrity checking, not authenticity verification.

## 11. Mutation model

The only logical mutation is first binding an artifact ID to its immutable
reference and ensuring the corresponding digest object is present. Existing
bindings and objects are never overwritten. Same-reference registration is a
successful no-op retry. Different-reference registration for an existing ID is
a conflict. There is no administrative withdrawal because deleting a binding
would make previously issued references backend-dependent; retention and
garbage-collection policy are deferred.

## 12. Atomicity and concurrency

Each public write is linearizable for a single logical registry. Concurrent
identical registration of one ID/reference has one logical effect and all
successful callers observe the same reference. Concurrent registration of one
ID to different references has exactly one winner and a structured conflict for
the other. Different IDs with the same digest may both register. Reads concurrent
with publication observe either the prior complete state or the complete new
binding; no reference may become visible before its verified bytes are available.

The in-memory backend prepares detached replacement state and publishes it at
one boundary under a per-instance lock. The filesystem backend uses separate
root-confined namespaces for immutable digest-keyed blobs and artifact-ID
bindings. Blob paths derive deterministically from canonical lowercase SHA-256;
binding paths derive deterministically from a safe encoding or hash of the
artifact ID, never by interpolating an arbitrary ID as a path. Temporary files
are created exclusively inside controlled directories on the same filesystem.
Fanout depth is an implementation detail; deterministic mapping, namespace
separation, and root confinement are contract requirements.

Filesystem `register` has this logical sequence: validate the reference and exact
bytes input; compute and check digest and declared size; stage complete bytes in
a private temporary file; if the digest location exists, compare its bytes and
reuse it only if equal (otherwise raise integrity/collision failure); prepare
the exact canonical reference metadata; acquire cross-process coordination for
the root; re-check the artifact-ID binding and digest location while coordinated;
publish a missing immutable blob without replacing existing content; then
atomically publish the binding. The
atomic binding publication is the logical commit point. Before that point the
registration is invisible through the API. Temporary files and even an
unreferenced published blob may remain after failure and are invisible; retry is
safe. After binding publication the registration is committed; retry follows
normal exact-replay or binding-conflict rules. No rollback of a committed
binding is required.

Publication uses temporary-file plus same-filesystem atomic publication/rename
semantics. Readers see no partial blob or binding and never a binding to a
partially written blob. All handles/processes sharing a supported root coordinate
through a cross-process mechanism (for example, a root-scoped interprocess file
lock); an in-process lock alone is insufficient. Same-ID identical concurrent
registrations converge to one binding; conflicting registrations have at most
one winner; different IDs for equal bytes may both succeed. Reader, verifier,
and writer observations are either pre-commit absence or complete post-commit
state.

The supported filesystem is a local filesystem with the required same-filesystem
atomic publication and locking semantics. Arbitrary network filesystems,
distributed stores, and filesystems lacking those primitives are outside the
contract and must not be represented as supported.

## 13. Idempotency and retry

There is no operation ID in `ArtifactRef` and Phase 3 does not add one. Retry
identity is the complete canonical `ArtifactRef` binding plus exact content
bytes as verified by digest and size. After successful registration, exact
replay returns the existing reference and does not create another binding. Any
different reference for that ID conflicts. A failed registration before binding
publication reserves no ID; retrying the same valid command can succeed. If an
orphan blob remains, retry verifies its bytes and reuses it only when equal.

## 14. Validation and failure model

All registry ingress revalidates caller-supplied Phase 1 references and takes
owned input. Structurally malformed references, non-exact `bytes`, malformed
IDs/digest syntax, invalid size type/range, and invalid metadata are structured
invalid-operation errors with field paths. A well-formed declared digest or size
that disagrees with supplied/stored bytes is an integrity failure, not invalid
input. No coercion or message parsing is required.

Proposed public failures:

- `ArtifactRegistryError`: common structured base, exposing stable `category`,
  `reason`, optional `artifact_id`, and optional `field_path` fields.
- `InvalidArtifactOperationError`: category `INVALID_ARGUMENT`; malformed
  registration/query input, digest syntax, size type/range, or `ArtifactRef`.
- `ArtifactNotFoundError`: category `NOT_FOUND`; no binding exists for the ID,
  with the requested `artifact_id`.
- `ArtifactBindingConflictError`: category `BINDING_CONFLICT`; the ID already
  has a different complete `ArtifactRef` binding, with `artifact_id`,
  `existing_reference`, and `requested_reference`.
- `ArtifactIntegrityError`: category `INTEGRITY_FAILURE`, with stable reasons
  such as `DIGEST_MISMATCH`, `SIZE_MISMATCH`, `CORRUPT_BINDING`, `MISSING_BLOB`,
  `CORRUPT_BLOB`, and `HASH_COLLISION`; include applicable expected/actual
  digest or size fields.
- `ArtifactStorageIOError`: category `STORAGE_IO_FAILURE`; physical I/O or
  coordination failures, with an operation field and underlying cause preserved
  through exception chaining.

These are the minimum stable cross-backend public error categories. Backend
`OSError`/`FileNotFoundError` is not the normal semantic API. Implementations
translate attributable storage/coordination failures while preserving causality;
they must not catch unrelated programming defects as storage failures. Ordinary
Python argument-binding misuse remains `TypeError`.

## 15. Serialization and canonicalization

Raw content is stored byte-for-byte; no text, Unicode, newline, compression, or
scientific normalization is applied before hashing. SHA-256 over the exact
unmodified bytes is the Phase 3 content address; empty, NUL-containing, and
non-UTF-8 byte sequences are valid payloads. No filename or metadata enters the
digest. Phase 1 accepts exactly 64 hexadecimal characters, accepts uppercase
hexadecimal input and canonicalizes it to lowercase, and rejects a `sha256:`
prefix. Phase 3 preserves those rules and does not modify Phase 1.

The operational storage assumption is that equal SHA-256 digests denote the
same content address. When a digest location is occupied and candidate bytes
are locally available, the backend must compare exact bytes: equal bytes reuse
the blob; different bytes raise `ArtifactIntegrityError` with reason
`HASH_COLLISION` and never overwrite the published blob. This is a storage
integrity safeguard, not scientific identity or equivalence inference.

The filesystem binding metadata is exactly the existing
`ArtifactRef.to_canonical_bytes()` (`sski-json-v1`) bytes; these bytes are
persisted as-is and deserialized/revalidated into the frozen Phase 1
`ArtifactRef`. Corrupt or malformed binding data is an integrity failure.
Binding metadata is atomically published in the deterministic, root-confined
binding namespace. Blob and binding namespaces are separate. Optional
filesystem bookkeeping is not scientific provenance; no second scientific
canonicalization or provenance model is introduced.

## 16. Persistence and backend boundary

`InMemoryArtifactRegistry` represents one independent, volatile logical
registry. It has no module-global state or import-time contents. The
`FilesystemArtifactRegistry` is a caller-configured supported local root. This
contract distinguishes operation atomicity, restart persistence, process-crash
consistency, and machine/power-loss durability. Completed committed files that
remain intact can be reconstructed after process restart; pre-commit temporary
files and orphan blobs are not visible as registrations. Phase 3 does not promise
survival across arbitrary power loss, filesystem corruption, controller-cache
loss, or hardware failure, and makes no stronger `fsync`-based power-loss
durability guarantee.

Both implementations satisfy the same public conformance tests. Future
PostgreSQL/object-store or distributed implementations may be added behind the
protocol in later phases; they must preserve ID-binding, digest, retry,
atomicity, and error semantics. No backend-specific capabilities leak through
the protocol.

## 17. Dependency direction

```text
Phase 1 domain (ArtifactRef, EntityId, canonical bytes)
                 ↑
Phase 3 artifact registry protocol/errors
                 ↑
       memory / filesystem backends
```

The domain imports no registry or storage code. Storage imports only the public
Phase 1 domain contracts, its Phase 3 protocol, and approved standard-library
facilities. Phase 3 does not import Phase 2 identity/storage internals, SimForge
application packages, or any future SSKI layer. SimForge remains downstream of
the immutable bundle boundary.

## 18. Deterministic behavior

For the same reference and bytes, all backends return equal `ArtifactRef` values
and equal bytes. SHA-256 is computed from exact bytes; no timestamp, filename,
insertion order, random ID, current working directory, or backend path affects
the digest. Exact repeated writes are idempotent. A missing reference is
distinguished from registered-but-corrupt content. No list API is proposed;
future enumeration must define stable ordering explicitly.

## 19. Security and trust assumptions

The configured filesystem root is trusted configuration and is not writable by
an untrusted actor. Internal blob and binding paths are deterministic and
root-confined; raw IDs are never paths. Implementations reject traversal and
absolute-path injection, create temporary files inside controlled registry
directories, and ensure symlink behavior cannot escape the configured root.
The public API never returns internal writable paths. This registry is not an
authorization service or tenant-isolation boundary. SHA-256 detects changes
relative to a trusted reference but does not prove who supplied the reference.
Media type and source attribution are not trusted without separate policy or
provenance checks.

## 20. Conformance test strategy

One backend-independent behavioral suite runs unchanged against a factory for
each backend and future storage implementations. Shared assertions cover
registration, exact replay, complete-reference binding conflicts, missing
artifacts, verified read/verify, digest and size mismatch, same bytes under
different IDs, different metadata under one ID, corruption detection,
defensive byte ownership, determinism, and absence of revisions/deletion/update.
The factory provisions isolated logical registries; filesystem tests can also
open multiple handles to the same root. Shared assertions use only the public
protocol and returned values. Backend-specific tests may inject I/O failures
without weakening shared semantics.

The in-memory and filesystem backends must pass the same registration, read,
verification, conflict, retry, ownership, and concurrency cases. The filesystem
test setup must verify reopen behavior after successful writes. No tests require
network access, external data, or SimForge execution.

## 21. Adversarial test matrix

- Valid registration and retrieval; absent ID; exact retry; multiple IDs sharing
  one digest; same ID with any differing `ArtifactRef` field conflicts.
- Wrong SHA-256, uppercase/malformed digest handling through Phase 1 validation,
  declared size absent/correct/incorrect, empty content, and large byte payloads.
- `bytes` subclass, `bytearray`, `memoryview`, text, path, stream, and mutable
  caller input ownership; returned bytes cannot mutate stored state.
- Mutated/bypassed `ArtifactRef` fields and nested/caller-owned scalar objects.
- Content changed after registration; truncated, replaced, missing, or corrupted
  blob; malformed reference metadata; verify and read failure classification.
- Occupied digest with equal bytes reuses the blob; occupied digest with unequal
  bytes raises integrity/collision failure and never overwrites it.
- Same-ID/same-reference concurrent writes across processes; same-ID/conflicting
  reference race; different IDs/same digest race; reader and verifier during
  atomic binding publication.
- Failure injection before/during temp write, after complete staging, before
  blob publication, after blob publication but before binding publication, and
  during binding publication; no pre-commit visible binding and safe retry.
- Filesystem restart/reopen, orphan handling, deterministic root confinement,
  symlink/path escape, unwritable root, unsupported atomic/locking semantics,
  corrupt canonical binding metadata, and independent roots.
- Import boundaries, no global mutable state, exact public exports, no SQL,
  network, source federation, importers, compiler, resolver, policy, bundle, or
  runtime dependencies.

## 22. Proposed implementation allowlist

This is a proposal for a later, separately authorized implementation pass.

### New production files

```text
simforge/knowledge/storage/artifact_registry.py
simforge/knowledge/storage/memory_artifact_registry.py
simforge/knowledge/storage/filesystem_artifact_registry.py
```

### New test files

```text
tests/knowledge/artifacts/conftest.py
tests/knowledge/artifacts/registry_contract.py
tests/knowledge/artifacts/test_records.py
tests/knowledge/artifacts/test_memory_registry.py
tests/knowledge/artifacts/test_filesystem_registry.py
tests/knowledge/artifacts/test_concurrency.py
tests/knowledge/artifacts/test_architecture.py
```

### Design documentation in this pass

```text
docs/sski/PHASE_3_DESIGN.md
docs/sski/ADR/ADR-012-artifact-registry-content-addressing.md
```

No Phase 1 or Phase 2 file, parent export, project configuration, unrelated
application code, schema/migration, network client, or runtime code is allowed.
The eventual production allowlist and implementation authority require ADR/design
acceptance and a separate user authorization.

## 23. Explicit deferred functionality

Source and importer integration; licensing decisions; claims/evidence records;
artifact metadata provenance; multiple locations/replication; remote/object
storage; SQL/PostgreSQL; garbage collection/retention; streaming large objects;
content transformation; format validation; immutable bundle compilation;
KnowledgeSnapshot; policies/resolution; BundleReader; runtime access; and
RunDerivationLedger integration remain deferred to their planned phases.

## 24. Migration and compatibility implications

No migration is required for Phase 1 or Phase 2 because they remain unchanged.
Existing `ArtifactRef` values are accepted by composition only when matching
bytes are registered. They do not promise that content is already stored.
Future bundle compilation must fail before the SSKI/SimForge boundary when a
required artifact cannot be retrieved or verified; that compiler behavior is
not implemented or designed as a Phase 3 API here.

Filesystem persistence begins with an empty caller-selected root. No existing
artifact repository or legacy SimForge `ArtifactRef` implementation is migrated
or aliased. If a later persistent backend needs a schema migration, that is a
separate reviewed phase.

## 25. Decision status and implementation gate

This proposal closes the implementation-blocking choices for Phase 3: IDs are
caller-provided and permanently bind to the complete exact `ArtifactRef`; bytes
alone define SHA-256 content addressing; digest collisions with unequal local
bytes fail as integrity errors; registration is bytes-only; filesystem binding
publication is the commit point coordinated across processes on supported local
filesystems; restart reconstruction is qualified by intact committed files and
there is no power-loss durability promise; canonical Phase 1 bytes are the
binding metadata; and the public structured error categories are defined in §14.
`read` always verifies before returning bytes, and `verify` returns the bound
reference on success.

No unresolved architectural decision blocks implementation. ADR-012 remains
Proposed and must be accepted, and a separate explicit user authorization is
still required before implementation. Incidental choices such as fanout depth
and the specific cross-process locking primitive may be selected during
implementation only if they meet the observable contract in §12 and are tested.
