# Phase 2 implementation audit — Revisioned Entity Registry

## Scope and checkpoint

Implementation was explicitly authorized on `feat/sski-phase2-registry`, starting
from clean `ee4e543486aff7274ae7a5d831dc7f038bc3e0b6`. ADR-011 is Accepted.
The user's explicit authorization and accepted substantive design govern this
implementation; historical Proposed/design-only wording in PHASE_2_DESIGN.md
was not treated as a new approval requirement. Neither that design nor ADR-011
was modified. No files were staged or committed.

The implementation uses only the Phase 2 allowlist. Phase 1 production/tests,
parent package exports, root test configuration, application code, and all
later-phase scaffolding remain unchanged. Both counters are per EntityId,
following the accepted design and ADR, rather than a registry-global counter.

## Files created

Production:

- `simforge/knowledge/identity/__init__.py`
- `simforge/knowledge/identity/records.py`
- `simforge/knowledge/identity/registry.py`
- `simforge/knowledge/identity/errors.py`
- `simforge/knowledge/storage/__init__.py`
- `simforge/knowledge/storage/memory_registry.py`

Tests:

- `tests/knowledge/identity/conftest.py`
- `tests/knowledge/identity/registry_contract.py`
- `tests/knowledge/identity/test_records.py`
- `tests/knowledge/identity/test_registry.py`
- `tests/knowledge/identity/test_aliases.py`
- `tests/knowledge/identity/test_concurrency.py`
- `tests/knowledge/identity/test_conformance.py`
- `tests/knowledge/identity/test_architecture.py`

Documentation:

- `docs/sski/audits/PHASE_2_IMPLEMENTATION_AUDIT.md`

Existing files modified: only the Phase 2 section of
`docs/sski/09_IMPLEMENTATION_PLAN.md`.

## Exact public API

| Module | Public objects |
|---|---|
| `identity.records` | `RegistryChange`, `RegistryWrite`, `RegistryRecord`, `ExternalIdentifierKey`, `ExternalIdentifierMatch` |
| `identity.registry` | `EntityRegistry` |
| `identity.errors` | `RegistryError`, `EntityNotFoundError`, `RegistryRevisionNotFoundError`, `EntityRevisionNotFoundError`, `RevisionConflictError`, `OperationConflictError`, `InvalidRegistryOperationError` |
| `storage.memory_registry` | `InMemoryEntityRegistry` |

The identity package exports exactly its thirteen contract objects. The storage
package exports no backend eagerly. The registry protocol and implementation
have exactly these five public operations:

```text
write(command: RegistryWrite) -> RegistryRecord
get(entity_id: EntityId, *, registry_revision: int | None = None) -> RegistryRecord
get_entity_revision(entity_id: EntityId, entity_revision: int) -> Entity
history(entity_id: EntityId) -> tuple[RegistryRecord, ...]
lookup_external(key: ExternalIdentifierKey) -> tuple[ExternalIdentifierMatch, ...]
```

No extra public classes, enums, aliases, factories, wire formats, or convenience
operations were introduced. The memory constructor takes no configuration.

## Contract and architecture findings

No unresolved BLOCKER or MAJOR finding remains in the reviewed Phase 2 code.
The audit was performed by the implementing agent, using source inspection and
the independent behavioral assertions; it is not a claim of external approval.

- **Descriptor boundary:** Exact Phase 1 Entity instances only. All seven rich
  scientific subclasses are explicitly tested for rejection before projection.
  Claim/evidence descriptors are accepted; payloads are not stored.
- **Dependency direction:** Public domain API -> identity contracts -> memory
  implementation. Static tests prohibit reverse dependencies, private domain
  helpers, application imports, networking, database/persistence, engines, and
  future-phase dependencies. Frozen Phase 1 remains byte-for-byte unchanged.
- **Revision semantics:** Entity and registry revisions start at one per ID;
  expected zero requires absence. Entity descriptors retain current revision
  unchanged or introduce precisely the next revision. Kind changes, gaps,
  backward revisions, and descriptor rewrites fail with their specified codes.
- **Atomicity and linearizability:** A private per-instance lock covers each
  public operation. Validation, snapshot construction, copied descriptor map,
  new complete index projection, and detached return value are prepared before
  publishing history, head, descriptors, index, and journal. Readers hold the
  same lock. No dependence on dictionary atomicity or cross-call transactions.
- **Idempotency:** The instance-owned successful-operation journal spans all
  entities. Normalized full-command comparison includes expected revision and
  attribution; domain values use Phase 1 canonical bytes, alias order is sorted,
  and times are UTC. Replay/conflict precedes CAS. Failed commands reserve no
  operation ID. Exact retries return detached copies of the original result,
  including after same-entity and unrelated-entity writes.
- **Alias semantics:** Lookup keys are lowercase namespace tokens and NFC
  accessions preserving case/whitespace. Every current declaration is retained,
  sorted by ID then canonical alias bytes, without lifecycle/time filtering,
  uniqueness selection, scientific equivalence, or evidence dereferencing.
- **Withdrawal:** Omission creates a new membership snapshot. Historical ACTIVE
  declaration/provenance remain ACTIVE and unchanged. No scientific lifecycle
  assertion is synthesized. Separately supplied retraction declarations retain
  supplied source attribution. Phase 2 structurally validates that attribution;
  it does not authenticate provenance or judge whether a source truly retracted
  its assertion. Two independent reassignment writes expose the permitted empty
  or overlapping intermediate mapping states.
- **Ownership:** Frozen slotted records, owned scalars/UTC offsets, revalidated
  nested domain copies, immutable tuples, and detached returned records prevent
  caller inputs or returned objects from changing internal history. Tests also
  deliberately bypass normal frozen assignment to verify copy isolation.
- **Failures:** Exact seven-class hierarchy, read-only structured fields, closed
  invalid-operation Literal codes, entity context when available, and all nested
  validation locations/codes. Ordinary signature misuse stays TypeError.
  Implementation defects are not translated into scientific conflicts.
- **Time:** No clock reads, generated IDs, timestamp-based identity, or revision
  ordering. Equal/backdated timestamps and timezone-equivalent retries pass.
- **Logical state:** Each memory instance has independent heads, history,
  descriptors, index, and journal. No module-global mutable registry state.
- **Serialization:** No registry object is added to sski-json-v1. The private
  command-comparison tuple is neither a persistence schema nor a content ID.

Audit corrections made within scope: retain available EntityId context in
invalid-input errors; classify malformed timezone offsets as INVALID_ARGUMENT;
retain EntityId context for invalid lookup-result aliases. Regression tests
cover each correction.

## Test organization and implementation findings

`registry_contract.py` contains backend-independent behavioral classes.
`test_conformance.py` runs those exact assertions using test-only provisioning:
a fresh registry factory and two handles plus a reopen callback for the same
logical state. Memory supplies the same instance. Future persistent fixtures
must supply distinct handles/reopening of the same durable state; no persistent
fixture or backend is implemented here.

Coverage includes registration for every descriptor kind/lifecycle, exact and
head reads, both histories, strict revision ingress, every lifecycle transition,
CAS, duplicate ID behavior, all retry outcomes, global operation-ID collisions,
barrier-synchronized conflicting/identical writes, complete lookup/write
boundaries, alias ambiguity/order/validity, withdrawal/reassignment, source
release/semantics changes, copy isolation, structured nested errors, independent
instances, exact API, cold imports, and deterministic results across hash seeds.

The test-first run failed as expected because the production modules did not yet
exist. Incremental Phase 2 runs passed as coverage was added. Combined collection
then exposed a collision between the approved `test_architecture.py` filename
and the frozen domain file of the same name. The allowlisted local conftest now
namespaces only Phase 2 test modules through pytest's collector hook. No test
rename, extra package file, root configuration change, or Phase 1 change was
needed. A subsequent combined run passed (740 tests at that point).

## Executed validation

Commands used the existing environment and `-q -p no:cacheprovider`.

| Gate | Before implementation | Final |
|---|---|---|
| `pytest tests/knowledge/identity` | Expected missing-module failure | 205 passed |
| `pytest tests/knowledge/domain` | 548 passed | 548 passed |
| Established focused SimForge gate | 770 passed, 2 skipped | 770 passed, 2 skipped |
| `pytest -m "not gmx and not rdkit"` | 3250 passed, 16 skipped, 4 deselected, 8 xfailed, 1 xpassed, 5 failed | 3455 passed, 16 skipped, 4 deselected, 8 xfailed, 1 xpassed, 5 failed |

The exact focused regression paths (unchanged from the Phase 0 audit) were:

```text
core/test_parser.py core/test_decision_engine.py core/test_semantic_inference.py
core/test_structural_annotation.py core/test_md_knowledge.py core/test_compiler.py
workflows/test_workflow_graph.py builders/test_workspace_builder.py
builders/test_workspace_selfcontained.py core/test_system_spec.py
core/test_system_build.py core/test_system_validate.py runtime/test_runtime.py
runtime/test_checkpoint_recovery.py executors/test_executor.py
executors/test_dag_blocking.py tests/test_workflow_resume.py ligand/test_campaign.py
ligand/test_normalization.py ligand/test_chemical_perception.py
tests/test_pdb2gmx_provenance.py tests/test_topology_assembly.py test_cli_smoke.py
test_cli_platform.py test_cli_ligand.py test_cli_integrate.py
tests/analysis/campaign/test_manifest.py tests/analysis/campaign/test_canonical_ids.py
tests/analysis/campaign/test_fingerprint.py
```

The broad baseline reproduced these five existing failures:

1. `ligand/test_integrate.py::TestPrepareCLI::test_prepare_with_ambiguous_h_ratio_blocks_without_rdkit`
2. `ligand/test_integrate.py::TestHydrogenationCLIOption::test_hydrogenation_none_unknown_status_exits_nonzero`
3. `ligand/test_integrate.py::TestHydrationStatus::test_status_unknown_for_plausible_ratio`
4. `ligand/test_integrate.py::TestHydrationStatus::test_complex_with_h_status_is_unknown`
5. `tests/test_membrane_water_v2_stability.py::test_one_scale_connection_is_reported_as_unresolved_topology`

The final failure set was compared programmatically with the captured baseline
and is identical: no new failures. All 205 additional passing tests are Phase 2
tests. The broad gate is not a zero-failure result. These application failures
were not repaired in this phase. A broad run started
before the last audit corrections was stopped and restarted against final code;
only the completed final run is reported in the final column.

Existing optional GROMACS/RDKit jobs were preserved and not separately run.
Both tools are present in this environment; unavailable-tool claims are not
made. The prescribed broad gate deliberately deselects their marked tests.
Existing Pydantic deprecation warnings remain unrelated.

The default command sandbox failed to start (`bwrap: setting up uid map:
Permission denied`). Approved escalated command execution was used for the original
implementation validation. During subsequent targeted remediation, an escalated
tool action was rejected by automatic approval after its usage limit; work resumed
with permitted workspace operations. No source changes were made outside the
approved Phase 2 paths.

## Targeted adversarial remediation validation

The adversarial implementation review reported 0 BLOCKER, 2 MAJOR, and 2 MINOR
findings. The two MAJOR findings were exception atomicity in `write()` and nested
TypeError leakage at the registry validation boundary. The two MINOR findings
were weak concurrency assertions and shadowable structured exception fields.
All four findings were addressed without changing Phase 1 or the accepted design.

`InMemoryEntityRegistry.write()` now prepares a complete replacement `_State`
containing history, heads, entity revisions, the current alias index, and the
successful-operation journal while holding the instance lock. One assignment
publishes that state after all preparation succeeds. Failure injection raises
`MemoryError` at every executed Python line occurrence in the preparation/write
path for both registration and revision writes. Each injected failure preserved
the full public pre-write observation; exact retry then committed once and
subsequent retry returned the committed result. This exercises the old
history-before-head and alias-before-journal failure windows as well as the other
pre-publication boundaries.

The identity boundary now snapshots and normalizes caller-provided nested
ExternalIdentifier timestamps narrowly before Phase 1 revalidation. Invalid
`tzinfo.utcoffset()` results are reported as structured `INVALID_ARGUMENT`
failures with their nested field path. Ordinary signature misuse remains
`TypeError`, and unrelated model-validation programming errors still propagate.
Registry exception fields use read-only data descriptors backed by owned
immutable values, so ordinary assignment, deletion, and `__dict__` shadowing do
not alter public structured data.

The reusable concurrency assertions now compare lookup reads with exact complete
pre-write and post-write snapshots, including revision, match ordering, aliases,
and current descriptor. Deterministically coordinated races cover same-entity
CAS, identical and conflicting operation IDs, independent entity targets, shared
keys, and lookup during publication. Disposable disabled-synchronization and
permanently stale-lookup controls fail those shared assertions. Backend-specific
scheduling remains in the memory test fixture; the reusable contract itself does
not inspect backend internals.

Final remediation validation, using `-q -p no:cacheprovider`:

| Gate | Result |
|---|---|
| `pytest tests/knowledge/identity/test_records.py` | 105 passed |
| `pytest tests/knowledge/identity/test_registry.py` | 20 passed |
| `pytest tests/knowledge/identity/test_concurrency.py` | 16 passed |
| `pytest tests/knowledge/identity` | 234 passed |
| `pytest tests/knowledge/domain` | 548 passed |
| `pytest tests/knowledge/domain tests/knowledge/identity` | 782 passed |
| Established focused SimForge gate | 770 passed, 2 skipped |
| `pytest -m "not gmx and not rdkit"` | 3484 passed, 16 skipped, 4 deselected, 8 xfailed, 1 xpassed, 5 failed |

The broad run's five failures are exactly the four known ligand hydrogenation
failures and the known membrane one-scale topology failure listed above; no
additional failure appeared. It is not a repository zero-failure result. The
Phase 1 production and test files remain unchanged. No Phase 3 functionality was
introduced. The review findings are considered remediated based on the targeted
and full validation above; no design deviation was made.

Final scope/whitespace checks passed: all 16 changed paths are allowlisted;
`git diff --check -- <all Phase 2 allowlisted paths>` is clean; new untracked
source and documentation files were separately scanned for trailing whitespace.
The only tracked diff is the Phase 2 ledger section. Phase 1 has
no diff against ee4e543. The staging area is empty, and final `git status --short`
contains only the fifteen new files and that one ledger modification.

## Deviations, assumptions, and deferred work

No design/API/scope deviation. Local test-module namespacing is an implementation
detail needed by the exact allowlist and existing pytest configuration.

The reference backend is volatile and keeps unbounded history/journal in memory;
index updates favor explicit all-or-nothing semantics over production-scale
performance. Actor text is attribution, not authentication. Evidence/source
truth, cross-store consistency, durable recovery, whole-plane snapshots,
persistence, scientific identity resolution, merges/splits, claims/evidence
storage, artifact registry/CAS, federation, networking, policies, compiler,
physical bundles, runtime integration, and ADR-010 integration remain deferred.
No Phase 3 functionality was introduced.
