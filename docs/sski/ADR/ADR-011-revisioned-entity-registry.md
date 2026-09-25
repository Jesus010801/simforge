# ADR-011 — Revisioned Entity Registry and Non-Resolving Alias Lookup

**Status:** Accepted

## Context
The mutable Global Knowledge Plane needs a registry foundation that preserves identity and provenance without introducing scientific resolution or persistent infrastructure in Phase 2. Phase 1 provides frozen scientific domain contracts.

## Decision
Phase 2 is a descriptor-only entity registry. Phase 1 remains the owner of scientific domain semantics. The registry accepts explicit Entity descriptors and does not flatten concrete scientific payloads into generic records or modify domain objects.

Entity revision and registry revision are distinct and scoped per EntityId. Historical snapshots are immutable; current heads and external identifier indices are mutable. Changes preserve prior revisions rather than overwriting history.

ExternalIdentifier values are attributed declarations, never internal primary keys. Current alias membership is distinct from historical/source assertions. Administrative alias withdrawal omits the declaration from a new current snapshot and index, preserves the original declaration unchanged in history, and records the membership action in RegistryChange. It is not scientific/source retraction and must not synthesize DEPRECATED or RETRACTED lifecycle state. DEPRECATED or RETRACTED for a previously ACTIVE alias requires a changed domain declaration with provenance appropriate to that changed assertion; the earlier ACTIVE declaration's provenance is not automatically sufficient. Administrative withdrawal != source-attributed lifecycle retraction.

Alias lookup is non-resolving and may return multiple EntityIds. External identifier lookup != scientific identity resolution != entity equivalence. Matching accessions or chemical representation strings never establish scientific equivalence. Merge/split implementation remains deferred.

`operation_id` belongs to logical registry state, not a Python object or backend handle. Its successful-operation journal belongs to that state. An identical successful retry returns the original committed result, even after the head advances; reuse with a different normalized command fails with a structured operation conflict.

Writes use optimistic compare-and-swap against the expected registry revision. Expected zero requires absence; an existing head requires an exact match. A competing write cannot silently overwrite the winner. Each successful write atomically commits its snapshot, head/index changes, and operation journal entry. Atomicity and per-call linearizability are backend-independent public contracts; separate writes do not form an atomic alias move or cross-entity transaction.

InMemoryEntityRegistry is the Phase 2 semantic reference implementation. Each instance represents one independent logical registry with its own volatile state and journal. For future persistent backends, handles to the same logical registry share the journal; retry through another handle and reopening the same persistent registry preserve operation-ID semantics. A fresh empty logical registry has no previous operation history. Persistent backend implementation remains deferred.

## Consequences
Registry history preserves attribution without treating administrative withdrawal as a source assertion or lookup as resolution. Backend implementations must satisfy the same observable behavioral contract; the detailed validation and conformance specification is in [the Phase 2 design](../PHASE_2_DESIGN.md). The same behavioral conformance assertions apply to the memory reference implementation and future persistent backends.

Phase 2 does not implement claims/evidence persistence, artifacts/CAS, federation, SourceAdapters, KnowledgeImporters, networking, PostgreSQL, policies, resolver, compiler, physical KnowledgeBundle, BundleReader, runtime/HPC integration, SimulationPlan redesign, or ADR-010 integration. Phase 1 production and tests remain unchanged.

Lock implementation, Python container types, PostgreSQL schema, exception wording, test filenames/organization, and incidental implementation mechanics are not fixed by this ADR. This ADR is accepted following final design review; implementation remains separately authorized.
