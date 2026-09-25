# Phase 2 design — Universal Entity Registry

**Status:** Proposed; ready for architectural review, not approved for implementation.
**Scope:** Design only. Phase 1 remains frozen.
**Recommendation:** Option A, narrowly limited to identity registration, revisioned registry records, external mapping lookup, and an in-memory reference implementation. No identity resolution or general persistence layer.

## 1. Purpose and architectural basis

Establish a small mutable Global Knowledge Plane component that can remember internal identities, preserve their registration history, and return attributed external mappings without treating those mappings as scientific equivalence.

The execution ledger is decisive: [Phase 2](09_IMPLEMENTATION_PLAN.md) explicitly specifies “Registry protocol + in-memory implementation + aliases/deprecation semantics.” Phase 3 owns artifact registration; Phase 4 owns general storage abstractions; Phase 15 owns PostgreSQL. The storage strategy in [the architecture](01_ARCHITECTURE.md) describes the eventual system, not permission to accelerate those phases.

This design implements the identity registry portion of that architecture. It does not implement the Scientific Knowledge Store. An identity descriptor is not the scientific payload belonging to that identity. The registry cannot answer whether a scientific assertion is true, whether its payload is available, or whether a workflow has complete knowledge.

Alternatives:

| Option | Disposition |
|---|---|
| A: identity + registry foundation | Recommended, bounded as below. |
| B: persistence/storage foundation | Deferred to Phases 4 and 15; no durable backend is needed to specify registry behavior. |
| C: identity-resolution foundation | Defer resolution algorithms and decision records. Phase 2 preserves candidates, evidence references, ambiguity and history needed later. |
| D: smaller combination | Interfaces alone omit the reference implementation explicitly planned for Phase 2. A registry protocol plus in-memory behavior is the smallest coherent planned combination. |

The hard boundary remains `SSKI -> KnowledgeBundle + DecisionContext -> SimForge`. Nothing in Phase 2 crosses it.

## 2. Inspection evidence and baseline

Read both repository instruction files, specifications 00–09, all ten accepted ADRs, the Phase 0 repository audit, the Phase 1 preimplementation baseline, every final domain module, and all eleven domain test files. Code takes precedence over older conceptual examples:

- `EntityId` is supplied, kind-free `sf:UUIDv4`; it does not generate IDs. The kind-prefixed examples in specification 03 are conceptual, not the implemented wire contract.
- `Entity` is an identity descriptor with `sf_id`, `entity_type`, `entity_revision`, and lifecycle. Tests explicitly distinguish a descriptor from a concrete `Protein`.
- The seven descriptor kinds are protein, protein_sequence, chemical_identity, chemical_species, experimental_structure, claim, evidence. There is no artifact kind in `Entity`.
- `ExternalIdentifier` targets an exact internal revision and type. It carries separate accession/source versions, validity, lifecycle, candidate/confirmed status, evidence references or a pinned confirmation rule, and import provenance.
- Confirmation validation checks attribution structure; it does not establish scientific truth or the existence of referenced evidence.
- Scientific models contain references across layers: for example, `ExperimentalStructure` requires a confirmed source mapping, claim/evidence references, and an artifact reference. Storing that object would expand Phase 2 beyond an identity registry.
- Canonical `sski-json-v1` serialization and deeply immutable values already exist. Public methods can be used without importing private domain helper classes.
- Existing architecture tests verify the exact fifteen-object domain API, static dependencies, clean-process isolation, and legacy type separation.

**Actual checkout:** HEAD observed as `b8f93bf12d1af65f1e737fc9470557aeaf687755`, with unrelated campaign/membrane changes. It is not a clean checkout of the supplied `f6bae5b`. `git diff f6bae5b -- simforge/knowledge/domain tests/knowledge/domain docs/sski` was empty before this document was created. The relevant Phase 1 content therefore matches that checkpoint. No reset, cleanup, staging, or commit is authorized. Recheck the SSKI baseline before a later implementation pass.

## 3. Exact scope and exclusions

### In scope

1. A synchronous, backend-independent `EntityRegistry` protocol.
2. Exact `Entity` descriptors only: internal identity, kind, internal revision, lifecycle. Concrete subclasses are rejected rather than silently reduced to descriptors.
3. Immutable registry snapshots containing a descriptor, its explicitly supplied external mappings, and administrative change attribution.
4. Append-only per-entity registry history and a mutable current head.
5. Atomic registration/revision with optimistic concurrency and retry-safe command IDs.
6. Exact namespace/accession lookup returning every matching declaration, including ambiguity.
7. Explicit entity lifecycle changes, attributed alias lifecycle declarations, and administrative alias membership withdrawal as separate operations, with retained historical reads.
8. A process-local, thread-safe in-memory reference implementation and a reusable registry conformance test suite.

Claim/evidence **identity descriptors** may be registered because they are already valid `Entity` kinds. This does not register a `Claim` or `Evidence` payload, its links, or its provenance graph. Callers must construct descriptors explicitly; no generic payload projection helper is provided.

### Out of scope

- Scientific entity payload repositories, including Protein/sequence/chemical/structure payload storage, sequence or chemical indices, and payload consistency across repositories.
- Claim/Evidence persistence, claims-for-subject queries, evidence-for-claim queries, contradiction graphs, confidence assessment, and evidence dereferencing.
- Artifact metadata registry, SHA256-to-location indices, CAS byte storage, multiple artifact locations, verification of bytes/corruption, download, or repair. These belong to Phase 3 and later storage/federation work. Phase 1 digest-bearing provenance remains carried as data; it is not an integrity verification result.
- General repository framework, unit-of-work abstraction, database schema/migrations, SQLAlchemy, PostgreSQL, SQLite, DuckDB, Parquet, filesystem persistence or recovery, object stores, and distributed services.
- UniProt, RCSB/wwPDB, PDB CCD, PubChem, membrane database access, and any other external source access; networking, transport/cache/retry infrastructure, SourceAdapter implementations or contracts, KnowledgeImporter implementations or contracts. Federation contracts begin in Phase 5.
- Scientific identity equivalence, exact-content deduplication of entities, probabilistic matching, sequence similarity, chemical-string heuristics, merging/splitting entities, canonical entity selection, and automatic redirects.
- Policy evaluation, Knowledge Resolver, compiler, physical KnowledgeBundle, KnowledgeSnapshot construction, BundleReader, runtime worker integration, RunDerivationLedger implementation, SimulationPlan redesign, and ADR-010 execution integration.
- CLI/configuration/packaging changes, legacy adapters, compatibility aliases, future-phase scaffolding, and changes to Phase 1 production or tests.

## 4. Package and dependency direction

Proposed additions to existing scaffolding:

```text
simforge/knowledge/
  domain/                       # unchanged
  identity/
    __init__.py                 # exports contract objects only
    records.py                 # registry command/snapshot/lookup values
    registry.py                # EntityRegistry protocol
    errors.py                  # typed registry failures
  storage/
    __init__.py                 # no eager backend export
    memory_registry.py         # InMemoryEntityRegistry

tests/knowledge/identity/
  conftest.py
  registry_contract.py          # shared behavioral assertions
  test_records.py
  test_registry.py
  test_aliases.py
  test_concurrency.py
  test_conformance.py
  test_architecture.py
```

Import arrows point toward dependencies:

```text
storage.memory_registry -> identity.registry / records / errors -> domain public API
```

`identity` must not import `storage`. `domain` must not import either. Neither layer imports existing SimForge application packages. Putting the single memory backend in `storage/` respects package ownership; it does not authorize Phase 4's general repositories or a backend plugin framework.

New immutable record types use frozen dataclasses with explicit validation, immutable tuples, owned scalar values, and validated copies of nested public domain objects. Do not subclass or import `_DomainModel`, `_Ref`, `_ImportProvenance`, or other private domain helpers. Domain values retain their existing validation and canonical serialization. New registry objects deliberately have no persistent wire format in Phase 2; no second canonical scientific serializer is introduced.

## 5. Exact proposed API

All method arguments below are typed; integer revisions reject booleans/coercion. Revision bounds are:

- `expected_registry_revision`: integer >= 0; zero means the caller expects no entity.
- Stored or explicitly requested `registry_revision`: integer >= 1.
- Explicit `entity_revision`: integer >= 1.

`None` is permitted only for the optional `get` registry revision and selects the current head. Invalid types or numeric ranges at registry-operation ingress raise `InvalidRegistryOperationError` with `code=INVALID_ARGUMENT`; no coercion is allowed. A valid but absent explicit revision on a known entity raises its corresponding typed not-found error. Section 9 defines constructor and nested-validation boundaries.

Public registry values reject unknown fields and retain no caller-owned mutable containers. Identifiers/strings are owned ordinary scalar values. UTC timestamps are supplied timezone-aware instants normalized to UTC; no registry clock is read.

| Object | Module under `simforge.knowledge` | Responsibility and fields | Mutability | Serialization / persistence | Visibility |
|---|---|---|---|---|---|
| `RegistryChange` | `identity.records` | `operation_id: UUID` (v4), `actor: str`, `reason: str`, `recorded_at: datetime`; actor/reason nonempty. Administrative attribution, not scientific evidence. | Frozen | No wire contract; retained in memory history and retry journal. | Public |
| `RegistryWrite` | `identity.records` | `entity: Entity`, `aliases: tuple[ExternalIdentifier, ...]`, `expected_registry_revision: int`, `change: RegistryChange`. Complete proposed state for one identity. | Frozen | No wire contract; normalized content retained for retry comparison. | Public |
| `RegistryRecord` | `identity.records` | `registry_revision: int`, `entity: Entity`, `aliases: tuple[ExternalIdentifier, ...]`, `change: RegistryChange`. One committed snapshot. | Frozen | Nested domain values keep existing serialization; enclosing record is process-local, not a persisted schema. | Public |
| `ExternalIdentifierKey` | `identity.records` | `namespace: str`, `accession: str`; exact validated lookup key, with Phase 1 NFC/lowercase-namespace rules. | Frozen | No wire contract; in-memory index key. | Public |
| `ExternalIdentifierMatch` | `identity.records` | `entity: Entity` (current descriptor), `registry_revision: int`, `alias: ExternalIdentifier` (exact target revision). Makes current lifecycle and pinned mapping target distinct. | Frozen | No wire contract; derived read result. | Public |
| `EntityRegistry` | `identity.registry` | Protocol with the five methods below. No implementation/default backend. | Interface | Not serializable or persistent. | Public |
| `InMemoryEntityRegistry` | `storage.memory_registry` | Implements the protocol for one independent logical registry. No constructor configuration. | Mutable, internally synchronized | Volatile; no load/save/export API. | Public, imported from this module only |
| `RegistryError` | `identity.errors` | Base exception for classified registry failures. | Exception object; structured attributes are read-only | No wire/persistence contract. | Public |
| `EntityNotFoundError` | `identity.errors` | Unknown internal ID. | As above | As above | Public |
| `RegistryRevisionNotFoundError` | `identity.errors` | Known ID, absent registry revision. | As above | As above | Public |
| `EntityRevisionNotFoundError` | `identity.errors` | Known ID, absent internal entity revision. | As above | As above | Public |
| `RevisionConflictError` | `identity.errors` | Expected registry head differs from actual head. | As above | As above | Public |
| `OperationConflictError` | `identity.errors` | Operation ID reused for a different normalized command. | As above | As above | Public |
| `InvalidRegistryOperationError` | `identity.errors` | Closed reason code and structured context for invalid state/transition. | As above | As above | Public |

Internal implementation details only: a lock; dictionaries for `(EntityId, registry_revision)` snapshots, current heads, entity revision descriptors, and lookup membership; a successful-operation journal. No additional public classes, enums, aliases, ID factory, or generic query objects. Private helpers may validate/copy/sort these values. Literal annotations define closed reason codes without extra public enum types.

Conceptual signatures (design notation, not executable implementation):

```text
write(command: RegistryWrite) -> RegistryRecord
get(entity_id: EntityId, *, registry_revision: int | None = None) -> RegistryRecord
get_entity_revision(entity_id: EntityId, entity_revision: int) -> Entity
history(entity_id: EntityId) -> tuple[RegistryRecord, ...]
lookup_external(key: ExternalIdentifierKey) -> tuple[ExternalIdentifierMatch, ...]
```

`get` without a revision returns the head, including deprecated/retracted identities. Explicit revisions never fall forward. `history` is ascending by registry revision. An unknown lookup key returns an empty tuple; an unknown internal ID raises a typed failure. No pagination, predicate language, latest-source-version selector, claims/evidence/artifact queries, or multi-query snapshot API is justified yet.

## 6. Data lifecycle and revision model

### Identity creation

The caller allocates an opaque UUIDv4 independent of source accessions and supplies a validated `EntityId`. The registry does not derive, recycle, generate, or change it. Registration is `write` with expected registry revision zero, an unused ID, and entity revision one. Any valid Phase 1 lifecycle can be registered: a first observation may already be deprecated. The registry supplies registry revision one. A descriptor reserves a kind for that ID permanently within this registry.

There is no globally coordinated ID authority in this phase. UUID collisions with an existing ID produce a conflict, never overwrite. Future artifact/policy registries must define cross-registry role validation separately; Phase 2 does not claim whole-plane role uniqueness.

### Two revision axes

- **Entity revision** identifies a caller-declared domain descriptor revision. It is not a source release or registry transaction number.
- **Registry revision** identifies the complete registration state (descriptor, aliases, change attribution) for one ID. Each successful new write advances it by exactly one.

A write against head `r` must declare expected registry revision `r`. It may retain the identical entity descriptor at its current entity revision, or supply precisely the next entity revision with the same ID and kind. A lifecycle change requires the next entity revision. Reusing an entity revision with a different descriptor is forbidden. No gaps, backward steps, or initial revision greater than one. A new entity revision may have otherwise identical fields: a later payload owner might have changed scientific content that the registry does not store. Registering that revision does **not** attest that such a payload exists.

Alias-only changes advance registry revision without inventing a scientific revision. Internal entity revisions and source versions never track each other implicitly. Old aliases continue to point at the exact old entity revision until the caller explicitly supplies a new declaration.

The snapshot store is append-only; only heads and derived lookup indices are mutable. There is no event replay engine, timestamp-as-version ordering, history pruning, physical delete, tombstone erasure, or automatic compaction. Full snapshots are adequate for the reference backend.

### Alias set and changes

Every command supplies the full alias tuple for this entity's current registry state. Administrative withdrawal of an `ExternalIdentifier` from the CURRENT registry index is a registry membership operation, not a source/scientific retraction. There is no implicit alias union or patch behavior.

If an ACTIVE declaration was attributed to source release R, withdrawing it means omitting it from the new current registry snapshot and index. The original declaration remains unchanged in historical revisions. `RegistryChange` records the administrative membership change. The registry MUST NOT synthesize DEPRECATED or RETRACTED lifecycle state merely because an administrator removed the mapping.

When changing a previously ACTIVE alias, DEPRECATED or RETRACTED may appear only in a changed domain declaration with provenance appropriate to that changed assertion. A first observation that is already DEPRECATED or RETRACTED likewise requires attribution appropriate to that declaration; registry membership changes never manufacture either state. The provenance of the earlier ACTIVE assertion is NOT automatically sufficient provenance for a later retraction. Phase 2 preserves and structurally validates supplied attribution; it does not authenticate source content or adjudicate scientific evidence. No Phase 1 model changes are required.

```text
administrative withdrawal != source-attributed lifecycle retraction
```

Withdrawal of alias X from entity A and assignment of X to entity B are TWO independent atomic writes. Phase 2 has no atomic "move alias" primitive. Withdrawal followed by assignment can expose no current mapping between writes; assignment followed by withdrawal can expose temporary overlap, which the many-to-many index permits. Neither order erases either entity's history or asserts entity equivalence.

Every alias must target this same internal ID and kind, and an entity revision registered previously or introduced by this command. Validate all aliases before publishing any state. Their evidence references are retained but not dereferenced. Their provenance must pass the unchanged Phase 1 contract. Exact duplicate canonical alias declarations within one write are rejected. Distinct declarations, including contradictory status, provenance, target revision, source release, or validity interval, are preserved.

New snapshots and returned values are deeply immutable. Domain values are revalidated and copied at ingress; concrete Entity subclasses are rejected before any conversion. No permissive dump-to-Entity path may discard scientific payload fields.

### Lifecycle and corrections

Allowed transitions are ACTIVE -> DEPRECATED/RETRACTED, DEPRECATED -> ACTIVE/RETRACTED, and RETRACTED -> ACTIVE/DEPRECATED, plus unchanged lifecycle. Every lifecycle change requires a new entity revision, expected-head match, actor, reason and time. Retaining an unchanged lifecycle does not by itself require a new entity revision. Restoration is an explicit correction; it does not undo history. Phase 2 does not assign scientific meaning to the reason text or authorize callers by actor name.

No replacement entity ID or automatic redirect is attached to deprecation. One deprecated identity may later relate to multiple successors. Merge/split/supersession relationships require a separately reviewed identity-resolution contract; using a single `replaced_by` field now would make split semantics incorrect. Source aliases can be retired and newly declared without asserting a merge or split.

### Historical reconstruction

Pin `(EntityId, registry_revision)` to recover the exact registry state consulted; pin `(EntityId, entity_revision)` for the descriptor. A registry record includes mappings and administrative attribution as they were known in that write. Source valid-time fields remain separate from transaction order. This provides per-entity history only while the in-memory instance lives. Durable recovery, whole-plane consistent snapshots, and reconstruction of scientific payload/evidence/bytes remain later requirements; no Phase 2 output is a KnowledgeSnapshot or a complete workflow reproducibility record.

## 7. External identity behavior and lookup

Normalize only according to existing domain rules: namespace must already be a lowercase valid token; accession is NFC text with its case and whitespace preserved. Do not strip accession suffixes, case-fold accession text, infer versions, or introduce provider-specific parsers.

The secondary lookup key is `(namespace, accession)`. Neither component is an internal primary key. Accession version, source record revision, source release/provenance, target entity revision, validity, status and lifecycle remain on each returned declaration. No source-specific uniqueness constraint is assumed.

`lookup_external` returns matching aliases from **current registry snapshots**, including candidate, confirmed, deprecated, retracted, future-dated and expired declarations. It does not consult a clock or filter them silently. The result includes the current descriptor so a caller can distinguish a retired entity from a retired alias. Historical aliases omitted from current state remain available through pinned snapshots/history, not a fabricated global time-travel query.

Ordering is deterministic: sort by internal ID string, then canonical alias bytes. Alias tuples in records use canonical alias-byte order. Duplicate exact declarations within one record fail; identical lookup keys on different IDs are allowed and return both. A confirmed mapping does not reserve an accession globally. Two confirmed mappings to different entities are visible ambiguity, not a database integrity error or an instruction to choose whichever arrived first.

A subsequent resolver can read candidates, inspect pinned rule/evidence references, and propose a new explicitly attributed write. Phase 2 provides no `resolve`, `confirm`, `merge`, `split`, or `canonicalize` method, and no empty resolver interface scaffolding. The read/write protocol is the available seam. Equal sequences, hashes, SMILES, InChI representations, or labels never merge IDs. Predictions remain predictions; indexing does not promote them into fact.

## 8. Atomicity, concurrency, and retry contract

A successful write atomically commits one snapshot, its head pointer, the entity-revision descriptor if new, the external lookup projection, and the successful-operation journal entry. Failure commits none of those. No entity becomes visible without its submitted aliases, and no failed alias validation advances a revision.

Write atomicity and per-call linearizability are PUBLIC backend-independent contracts over logical registry state. Each read obtains a consistent immutable result for that call. These guarantees do not create a transaction spanning separate calls. A process-local lock is only an implementation detail of `InMemoryEntityRegistry`, not a protocol requirement. Cross-entity writes, multi-record unit-of-work transactions, distributed coordination implementations, and cross-store atomicity are deferred. Phase 2 implements process-local access only; future persistent handles must obey the same public state semantics.

`operation_id` belongs to LOGICAL REGISTRY STATE, not an individual Python object or backend handle. Successful operation IDs are unique across all entities in that logical registry, and the successful-operation journal is part of that state.

- For `InMemoryEntityRegistry`, one instance represents one independent logical registry; its operation journal belongs to that registry state. Discarding the instance discards that volatile state.
- For future persistent backends, multiple handles to the same logical registry share its operation journal. Retrying an already successful operation ID through another handle returns the original committed result. Reopening the same persistent logical registry preserves operation-ID semantics and its journal.
- A fresh empty registry is a different logical registry with no previous operation history. Opening another handle to existing state does not create a fresh registry.

Compare a retry's normalized full command, including expected revision and attribution, with the originally accepted command. Domain comparison uses canonical bytes; commands normalize alias order and UTC instants. Do not compare strings from exception messages or object identity. `recorded_at` is caller-provided attribution included in command comparison, not a commit clock. Equal or backdated timestamps are permitted; registry revision order remains authoritative.

- Identical successful retry returns the original `RegistryRecord`, even if the head has since advanced. It does not reapply the change.
- Reusing that operation ID with different content raises `OperationConflictError` before head checks.
- Failed operations do not reserve an operation ID.
- A new operation ID with an unchanged proposed state is a new attributed snapshot and advances registry revision. Idempotence is by operation ID, not an inferred scientific equivalence.
- Two concurrent registrations of the same internal ID with different operation IDs: exactly one succeeds, the other gets `RevisionConflictError(expected=0, actual=1)`.
- Two identical commands with the same operation ID: both return the same committed record.
- The same operation ID concurrently submitted with different normalized, otherwise valid commands: one commits and the other receives `OperationConflictError`, including when commands target different EntityIds. The conflict check spans the logical registry, not only one entity.
- Two commands against the same expected head: one succeeds, one conflicts; no last-writer-wins.
- Two different IDs sharing an external key: both registrations succeed; lookup exposes both. This is not a uniqueness violation.

A revision conflict requires rereading and a deliberate new command; there is no automatic overwrite/retry loop. The in-memory backend promises no durability or crash recovery. A later database backend must meet this observable contract with atomic transactions/constraints; no particular SQL isolation setting or schema is mandated now.

## 9. Failure model

Direct Phase 1 domain construction keeps its existing validation behavior unchanged. Registry semantic/input validation uses the public structured registry exceptions; messages are explanatory only and never the machine-readable contract.

Ordinary Python signature/constructor misuse, such as missing required arguments, unexpected keyword arguments, or duplicate argument binding, raises `TypeError` where Python rejects the call before validation. Once arguments are bound, invalid registry field values and operation inputs raise `InvalidRegistryOperationError` with structured `code` and `field_path`. Invalid types or numeric ranges use `code=INVALID_ARGUMENT`, including booleans, floats, strings, negative expected revisions, and zero/negative explicit requested revisions. No coercion occurs. The specialized codes below distinguish otherwise well-typed but invalid registry transitions and targets.

Nested Phase 1 validation failures encountered during registry construction/revalidation at ingress are translated to `InvalidRegistryOperationError(code=INVALID_ARGUMENT)`. Preserve structured validation locations by prefixing the domain location with the registry field path (for example, `('aliases', 0, 'provenance', 'source_record', 'version_status')`). The exception additionally retains an immutable `validation_issues` tuple of `(field_path, validation_code)` pairs for all nested issues, using the domain validator's structured error codes. Other registry failures use an empty tuple. No client or implementation may parse exception message strings to determine meaning. This translation belongs to the registry boundary and does not modify Phase 1 errors.

The following ingress outcomes are normative for every backend:

| Input or operation | Outcome |
|---|---|
| expected_registry_revision < 0, or any revision supplied as bool/float/string | InvalidRegistryOperationError with code=INVALID_ARGUMENT and the input field path; no coercion. |
| Explicit registry_revision or entity_revision <= 0 | InvalidRegistryOperationError with code=INVALID_ARGUMENT; do not classify it as a missing revision. |
| Positive registry_revision absent for a known EntityId | RegistryRevisionNotFoundError. |
| Positive entity_revision absent for a known EntityId | EntityRevisionNotFoundError. |
| Valid expected_registry_revision mismatches current head | RevisionConflictError; the expected head is a precondition, not a historical query. |
| Python rejects argument binding before validation | TypeError; semantic registry failures do not use this route. |
| Nested domain validation fails at registry ingress | InvalidRegistryOperationError with code=INVALID_ARGUMENT and structured prefixed locations/validation codes. |

For reads, a valid explicit revision that is absent on a known entity raises `RegistryRevisionNotFoundError` or `EntityRevisionNotFoundError`. An unknown identity raises `EntityNotFoundError` first, after argument validation. For writes, a valid expected revision that differs from the current head raises `RevisionConflictError`, including actual zero for an absent entity; expected revision is a compare-and-swap precondition, not a historical-read request.

| Failure | Required structured fields | Condition |
|---|---|---|
| `EntityNotFoundError` | `entity_id` | Read of an unregistered identity. |
| `RegistryRevisionNotFoundError` | `entity_id`, `registry_revision` | Known identity, missing registry snapshot. |
| `EntityRevisionNotFoundError` | `entity_id`, `entity_revision` | Known identity, missing descriptor revision. |
| `RevisionConflictError` | `entity_id`, `expected`, `actual` | Head mismatch; actual zero denotes absence. |
| `OperationConflictError` | `operation_id` | Successful command ID reused with changed normalized command. |
| `InvalidRegistryOperationError` | `code`, `entity_id` when available, `field_path: tuple[str | int, ...]`, `validation_issues: tuple[tuple[tuple[str | int, ...], str], ...]` | Invalid typed request or transition; nested issues preserve locations and validation codes. |

Closed `InvalidRegistryOperationError.code` values: `INVALID_ARGUMENT`, `UNSUPPORTED_ENTITY_MODEL`, `ENTITY_TYPE_CHANGE`, `ENTITY_REVISION_SEQUENCE`, `ENTITY_REVISION_REWRITE`, `ALIAS_TARGET_MISMATCH`, `ALIAS_TARGET_REVISION_MISSING`, `DUPLICATE_ALIAS_DECLARATION`.

Validation order: structural command validation; successful-operation replay/conflict check; expected-head check; entity transition validation; alias target/duplicate validation; atomic commit. For reads, validate arguments, then identity existence, then requested revision existence. No partial writes for any failure. Validation details identify field paths without requiring clients to parse strings.

Do **not** add `DuplicateExternalIdentifierError`: multiple attributed mappings are legal. Do **not** add `IdentityConflictError`: scientific conflict adjudication is deferred; lookup returns the conflicting declarations. Do **not** add `IntegrityFailureError`: this phase does not verify bytes or operate a CAS. Internal programming defects remain defects, not misleading scientific error categories.

## 10. Test design (future implementation only)

No tests are implemented or executed in this design pass. The future tests must check behavior independently of private index layout.

| Area | Required cases |
|---|---|
| Typed records | Strict revisions/UUIDs/timestamps; immutable nested values; caller container/scalar mutation; unknown fields; invalid namespace; preserved accession case/whitespace; exact Entity accepted and concrete subclasses rejected. |
| Registration | Initial revision one, unknown-ID reads, first observation deprecated/retracted, duplicate same ID, independent registries, claim/evidence descriptors accepted but payloads rejected. |
| Domain revision invariants | Same kind forever; gaps/backward steps/rewrite rejected; alias-only writes preserve entity revision; lifecycle changes require new entity revision; source version unrelated to both internal counters. |
| Lifecycle/history | Every allowed transition including correction; current and pinned reads; retained withdrawn aliases; no redirect, deletion or historical mutation; old reads unchanged after later writes. |
| Alias attribution | Missing target/revision/type mismatch; exact duplicate declaration rejected; distinct evidence/provenance declarations retained; evidence references not dereferenced; canonical order. |
| External collisions | Same external key across IDs, source versions, candidates and confirmations; contradictory confirmations remain visible; no canonical winner; equal chemical representations or sequence fingerprints cannot deduplicate descriptors. |
| Validity | Open/closed validity metadata preserved; expired/future aliases returned; no ambient time filtering; source supersession represented as declarations, not entity merge. |
| Atomicity | Invalid last alias leaves no descriptor/head/index/journal change; replacement of current alias set is all-or-nothing; failed operation ID can be reused. |
| Concurrency | Barriers, not sleeps: same-ID competing creation; identical concurrent retry; competing revisions; different IDs with same external key; reader sees complete pre- or post-write state. Exactly specified success/conflict counts. |
| Retry semantics | Same operation after head advances returns original result; normalized alias order/timezone equivalence; changed actor/reason/expected revision causes operation conflict; no-op with new operation advances history. |
| Determinism | Fixed IDs/timestamps produce identical descriptor bytes, alias order, history and errors across processes/hash seeds; insertion order does not choose an identity. Thread scheduling may choose a write winner, never corrupt state. |
| Backend conformance | The same shared behavioral suite runs unchanged against memory and future persistent backends through a test-only fixture supplying a fresh logical registry and handles to that same state. Memory may supply the same instance for both handles; future persistent fixtures supply distinct handles. Fixtures own setup/cleanup, never a public reset/clear API. Only memory fixtures are implemented in Phase 2; future backend provisioning is deferred. |
| Fresh interpreter | Import identity contracts and memory backend independently in child interpreters; block scientific network operations, database connections, application/engine imports and filesystem persistence. Exercise writes/reads under guards. Account for import-time Python file loading separately. |
| Dependency direction | AST/import-graph checks: domain unchanged, identity does not import storage/private domain helpers, memory imports contracts, parent imports do not eagerly initialize backend. No new public domain names or compatibility aliases. |
| Compatibility | Existing fifteen-domain-object tests unchanged and passing; legacy Evidence, ArtifactRef, ChemicalIdentity remain distinct. |
| Concurrent operation-ID collision | Submit the same operation ID with different normalized commands, both for the same EntityId and for different EntityIds. With otherwise valid commands, exactly one commits; the other raises OperationConflictError, with no partial second write. |
| Administrative withdrawal versus retraction | Start with ACTIVE alias X attributed to release R. Omit X in a new snapshot; current lookup loses X, historical X remains ACTIVE with unchanged provenance, and RegistryChange attributes membership withdrawal. No DEPRECATED/RETRACTED declaration is synthesized. A separately supplied lifecycle declaration retains attribution appropriate to the changed assertion; earlier ACTIVE provenance is not treated as automatic support. |
| Alias reassignment | Withdraw X from A, then assign X to B as two independent atomic writes. Observe an empty intermediate lookup. Reverse the write ordering and observe permitted temporary overlap. Both histories remain intact; no atomic move or identity equivalence is inferred. |
| Source release/semantics changes | Submit distinct attributed declarations under the same namespace/accession key across source release and source-semantics changes. Preserve source metadata and all supplied current declarations without deduplication or inferred equivalence; prior snapshots remain unchanged. |
| Equal timestamps | Commit two distinct operations with equal recorded_at values; registry revisions advance and determine history order. |
| Backdated timestamps | Commit a later operation with an earlier recorded_at value; accept it and retain ascending registry revision order. Time never determines revision validity or identity. |
| Shared logical-registry retry | Commit through handle A and retry through handle B for the same logical state, including after head advancement: return the original result without a new revision. Future persistent fixtures must use different handles; reopening the same persistent state must preserve the journal. A fresh empty logical registry has no previous operation history. |

The conformance suite's behavioral assertions must be reusable unchanged against `InMemoryEntityRegistry` and a future persistent backend. Backend-specific fixtures may supply handles and arrange state, but must not weaken expected outcomes or replace assertions. The shared-state retry case must compare the original committed result, history length, and current head to prove that retry did not append a revision. In the future persistent fixture, the two handles must be distinct and connected to the SAME logical registry; reopening tests retain that same state. These fixture provisions do not authorize persistent-backend implementation in Phase 2.

Future verification sequence: new registry suite plus unchanged `tests/knowledge/domain/`; focused existing regressions enumerated in the Phase 0 audit; existing non-GROMACS/non-RDKit gate `pytest -m "not gmx and not rdkit" -q -p no:cacheprovider`. Preserve optional CI tests and report unavailable environments. Capture an implementation-time baseline rather than assuming the dirty current tree matches past counts.

The Phase 1 ledger reports 548 domain passes and no new isolated regressions; it also documents five existing isolated failures (four ligand hydrogenation failures plus the membrane one-scale topology test). The preimplementation baseline lists only the original four. These are historical reports, not results rerun or independently reproduced in this design pass.

## 11. Compatibility and migration

No migration from existing SimForge models is needed or permitted. `core.md_knowledge.evidence.Evidence`, `runtime.artifacts.ArtifactRef`, and `ligand.campaign.ChemicalIdentity` retain their separate contracts. No imports, conversions, adapters, changes to consumers, or compatibility aliases are introduced.

No durable SSKI registry exists to migrate. Memory contents are intentionally disposable. Freezing Phase 1 means no new entity kinds, public domain references, serialization versions, or revision fields are added there. A future durable schema and payload repositories require their own review, including consistency between registered descriptors and scientific payloads.

## 12. Implementation sequence after approval

1. Recheck relevant Phase 1 tree against the approved checkpoint and capture current regression baseline without altering unrelated work.
2. Review and accept the bounded registry ADR described below. Confirm this exact scope and allowlist; design readiness is not implementation authorization.
3. Add typed records/errors/protocol and their contract validation tests.
4. Add the single in-memory backend; cover registration, reads, revisions and retry behavior with the conformance suite.
5. Add alias/lifecycle, atomicity and concurrency cases; no later-layer collaborators.
6. Run the planned domain, registry, architecture and regression checks. Inspect implementation against the SSKI invariants and this design.
7. Fix all Phase 2 BLOCKER/MAJOR findings within scope; broader changes require a revised design.
8. Record changed files, actual test results, baseline failures, limitations and the architecture audit. Update only the Phase 2 ledger section. Commit only if the later implementation instruction authorizes it; this pass explicitly forbids a commit.

## 13. Review concerns and decisions

Severity means: BLOCKER prevents coherent implementation; MAJOR threatens an invariant unless the specified disposition is adopted; MINOR is a bounded limitation or documentation issue.

| Severity | Concern | Proposed disposition / review question |
|---|---|---|
| BLOCKER if rejected | The descriptor/payload boundary must be explicit. A generic `Entity` repository accepting subclasses would silently become a scientific store or lose fields. | Accept exact descriptors only, rejecting subclasses. If full payload storage is desired, redesign scope before implementation rather than expanding implicitly. |
| MAJOR | Mutable alias metadata must not rewrite pinned scientific revisions. | Accept independent per-entity registry revisions and immutable historical snapshots. |
| MAJOR | Assuming one accession maps to exactly one entity hides scientific disagreement. | Accept many-to-many declaration lookup; no uniqueness error or automatic canonical selection. |
| MAJOR | Confirmed mappings can reference evidence not present locally. | Preserve Phase 1 structural checks; explicitly disclaim evidentiary closure. Claim/evidence storage and resolution remain deferred. |
| MAJOR | Whole-plane snapshots and atomic descriptor/payload writes are not provided. | Accept per-call/per-entity guarantees now; revisit cross-store consistency in Phase 4 and pinned snapshots before compilation. |
| MAJOR | Supersession support in specification 03 includes merged/split records. | Retain aliases and history now; defer scientific relationship/merge/split semantics. Do not pretend deprecation alone models scientific supersession. |
| MAJOR | Process-local history is not durable historical reconstruction. | Reference backend only; persistent retention and recovery are Phase 15 concerns. No compiler/runtime may depend on this backend as a reproducibility artifact. |
| MAJOR, deferred | ADR-010 integration gap and legacy artifact/evidence/chemical semantics remain. | No Phase 2 bridge or redesign. |
| MINOR | Requested clean checkpoint differs from checkout HEAD and unrelated dirty work exists. | Relevant SSKI files match `f6bae5b`; preserve all unrelated work and recheck before implementation. |
| MINOR | Snapshot history and operation journal are unbounded in memory; lookup/history have no pagination. | Accept for reference scale. Production retention/performance requirements are not inferred now. |
| MINOR | Registry change actor is attribution, not authenticated authorization. | No security/service interface in this phase; later service design owns authorization. |
| MINOR | Older kind-prefixed ID examples differ from final Phase 1. | Use frozen kind-free UUIDv4 contract; document clarification in the new ADR without changing Phase 1. |

**Review recommendation:** adopt the dispositions above. There is no unresolved architectural BLOCKER in this proposed scope. Its tradeoffs require ordinary design approval, not missing scientific inputs or an unchosen implementation alternative. If reviewers reject the descriptor-only, ambiguous-lookup, or independent-registry-revision decisions, return to design before coding.

Future questions intentionally deferred: scientific merge/split decision provenance and reversal; cross-store transactional closure; global snapshots; source-specific namespace/version interpretation; payload repository schema; artifact ID/digest and multi-location semantics; durable retention/migration and production scale. No placeholder interfaces for these are authorized.

## 14. ADR required before implementation

[ADR-011 — Revisioned Entity Registry and Non-Resolving Alias Lookup](ADR/ADR-011-revisioned-entity-registry.md) records the durable decisions specified in this remediation. Its status is Proposed pending final review; creation does not authorize implementation.

The ADR records descriptor-only scope with Phase 1 ownership of scientific semantics; distinct entity/registry revisions; immutable history and mutable heads/indices; declarations versus current membership; administrative withdrawal versus source-attributed lifecycle changes; non-resolving many-to-many lookup; logical-registry operation IDs and original-result retries; optimistic compare-and-swap; backend-independent atomicity and linearizability; the memory reference implementation and deferred persistence.

One ADR is sufficient. Lock implementation, Python containers, PostgreSQL schema, exception wording, test organization, and other incidental mechanics are deliberately not durable ADR decisions. Detailed validation boundaries and conformance scenarios remain in this design. ADRs 001–010 and Phase 1 remain unchanged. Final acceptance of ADR-011 is required before the separately authorized implementation pass.

## 15. Exact future implementation allowlist

This is a proposed allowlist for a subsequent explicitly authorized implementation pass, not permission to edit these files now.

### New production files

```text
simforge/knowledge/identity/__init__.py
simforge/knowledge/identity/records.py
simforge/knowledge/identity/registry.py
simforge/knowledge/identity/errors.py
simforge/knowledge/storage/__init__.py
simforge/knowledge/storage/memory_registry.py
```

### New test files

```text
tests/knowledge/identity/conftest.py
tests/knowledge/identity/registry_contract.py
tests/knowledge/identity/test_records.py
tests/knowledge/identity/test_registry.py
tests/knowledge/identity/test_aliases.py
tests/knowledge/identity/test_concurrency.py
tests/knowledge/identity/test_conformance.py
tests/knowledge/identity/test_architecture.py
```

### Documentation

```text
docs/sski/ADR/ADR-011-revisioned-entity-registry.md   # review/approval prerequisite
docs/sski/PHASE_2_DESIGN.md                        # approved clarifications only
docs/sski/audits/PHASE_2_IMPLEMENTATION_AUDIT.md    # actual findings/results
docs/sski/09_IMPLEMENTATION_PLAN.md                # Phase 2 status/results only
```

No other files may be created, edited, moved or deleted. Existing `.gitkeep` files may remain. In particular, `simforge/knowledge/__init__.py`, all domain production/tests, parent package exports, root conftest, project dependency/configuration files, application packages, and unrelated dirty files remain unchanged.

Allowed production dependencies: public `simforge.knowledge.domain` objects/methods, the new identity modules, standard-library `typing`, `dataclasses`, `datetime`, `uuid`, `unicodedata`, `re`, and `threading` (backend only). No dependency additions. Validation helpers stay private to identity and cannot change scientific normalization rules. Tests may additionally use pytest, AST/import inspection, subprocesses, pathlib, JSON, environment control, and standard-library concurrency utilities. Network, database, filesystem persistence, chemistry execution, engines and existing application imports are forbidden in new production code.

## 16. Design-pass validation record

The original design pass created only `docs/sski/PHASE_2_DESIGN.md`. This remediation edits that document and creates `docs/sski/ADR/ADR-011-revisioned-entity-registry.md`, resolving the two MAJOR and two MINOR review findings without changing scope. No production/test code or existing application packages were modified. No tests were implemented or run, and no commit was created. Validation consists of required document/code/test inspection, baseline comparison, and documentation diff/whitespace and scope checks. No test-pass claims are made for the proposed implementation. The normal sandbox initially failed with `bwrap: setting up uid map: Permission denied`; repository inspection proceeded through approved escalated command execution.
