# SSKI Implementation Plan

This is the execution ledger for Codex-assisted development. Do not skip phases merely because later work appears straightforward.

## PHASE 0 — Architecture contracts and repository audit
**Status:** DONE

Goal: establish repository-aware constraints.

Allowed: `AGENTS.md`, `docs/sski/`, `simforge/knowledge/AGENTS.md`, scaffolding only.

Forbidden: production SSKI implementation, API code, DB schemas, CLI changes, DecisionEngine refactors.

Exit: no unresolved BLOCKER about package placement/integration boundaries.

Phase 0 validation: reviewed and accepted with no BLOCKER findings. `simforge/knowledge/` is approved as the canonical SSKI package location. Existing execution/artifact/evidence/chemical models remain isolated from the SSKI domain layer. The ADR-010 integration gap is deferred beyond Phase 1.

---

## PHASE 1 — Pure scientific domain model
**Status:** DONE

Objects: EntityId, ExternalIdentifier, Entity, Protein, ProteinSequence, ChemicalIdentity, ChemicalSpecies, ExperimentalStructure, ArtifactRef, Claim, Evidence, PolicyRef, Decision, DecisionContext, KnowledgeBundleManifest.

Allowed: `simforge/knowledge/domain/`, `tests/knowledge/domain/`, minimal exports.

Forbidden: HTTP, SQLAlchemy/PostgreSQL, SQLite/DuckDB, simulation engines, CLI, source-specific adapters.

Acceptance: external IDs are not PKs; Claim != Evidence; deterministic serialization; SHA256 artifact semantics; no persistence/network/engine imports.

Implementation: all 15 approved public contracts are implemented in the 12 allowed domain files, with 11 new domain test files. The approved contract uses kind-free UUIDv4 identities, separate internal/source revisions, deeply immutable Pydantic 2 values, typed SHA256 digests, canonical NFC/Decimal `sski-json-v1`, versioned sequence fingerprints, distinct claims/evidence and chemical identities/species, extensible experimental methods, four-state decisions, pinned deferred operations, and value-level manifest consistency. Parent-package exports remain unchanged.

Final validation:

- **Final domain validation:** 548 passed, no failures, existing warnings only.
- **Final focused regression gate:** 770 passed, 2 skipped, 0 failed; matches the original focused baseline.
- **Final isolated repository validation:** 3250 passed, 16 skipped, 4 deselected, 8 xfailed, 1 xpassed, 5 failed.

The isolated repository validation was performed in a detached worktree based on `ad212fc` — `docs(knowledge): establish SSKI architecture contracts`. Only the final Phase 1 overlay was copied into that worktree:

- `simforge/knowledge/domain/`
- `tests/knowledge/domain/`
- `docs/sski/09_IMPLEMENTATION_PLAN.md`

Four isolated failures are the previously documented ligand hydrogenation baseline failures:

1. `ligand/test_integrate.py::TestPrepareCLI::test_prepare_with_ambiguous_h_ratio_blocks_without_rdkit`
2. `ligand/test_integrate.py::TestHydrogenationCLIOption::test_hydrogenation_none_unknown_status_exits_nonzero`
3. `ligand/test_integrate.py::TestHydrationStatus::test_status_unknown_for_plausible_ratio`
4. `ligand/test_integrate.py::TestHydrationStatus::test_complex_with_h_status_is_unknown`

The fifth isolated failure is `tests/test_membrane_water_v2_stability.py::test_one_scale_connection_is_reported_as_unresolved_topology`. It was independently reproduced against clean `ad212fc` with no Phase 1 overlay and against the same baseline with the final Phase 1 overlay. Ten repeated executions on the clean baseline failed consistently; ten repeated executions with the overlay also failed consistently. This membrane-water failure is pre-existing relative to Phase 1, not an SSKI regression.

**Phase 1 introduced zero new failures relative to the isolated baseline.** The repository validation still has the five pre-existing failures listed above; it is not a zero-failure result.

A separate repository-wide run in the user's active dirty working tree produced additional `analysis/campaign` failures caused by concurrent user-owned campaign development. These failures disappear in the isolated `ad212fc` + Phase 1 worktree and are explicitly excluded from Phase 1 regression attribution.

Static dependency checks, fresh-child-interpreter isolation, cross-process canonical serialization, golden vectors, immutability/bypass checks, and manifest consistency tests passed.

Final architecture review: **Phase 1 final review passed.** No unresolved BLOCKER findings and no unresolved MAJOR findings remain. The scalar ownership, DEFERRED scientific-only schema, and ExperimentalStructure classification issues were remediated and adversarially revalidated. Exactly 15 approved public objects remain. No Phase 2 implementation was introduced. Unrelated user-owned work was preserved.

Deferred: physical bundle/compiler/runtime guarantees, persistence, federation, registry access, scientific authenticity/equivalence resolution, and the accepted SimulationPlan/ADR-010 integration gap. Phase 2 has not started.

---

## PHASE 2 — Universal Entity Registry
**Status:** IMPLEMENTED — READY FOR REVIEW (UNCOMMITTED)

Implemented the accepted descriptor-only registry foundation and ADR-011:
immutable per-entity descriptor/registry histories, mutable current heads and
non-resolving alias index, typed records/errors and five-method EntityRegistry
protocol, and synchronized InMemoryEntityRegistry. Successful-operation journals
belong to logical registry state; exact retries return original results before
CAS checks. Administrative withdrawal preserves historical source declarations.
Phase 1 production and tests remain unchanged. No Phase 3 functionality,
persistence, scientific resolution, or runtime integration was introduced.

Initial implementation validation: 205 Phase 2 tests passed; the subsequent
targeted remediation validation passed 234 Phase 2 tests, 548 frozen domain
tests, and 782 combined domain/identity tests. The established focused regression
gate passed 770 tests with 2 skipped. The broad gate passed 3484 tests, with
16 skipped, 4 deselected, 8 xfailed, 1 xpassed, and 5 failed. The five failures
exactly match the recorded pre-existing baseline (3250 passed before Phase 1):
the four previously documented ligand hydrogenation failures and the membrane
one-scale topology failure. No new failures were introduced. The adversarial
implementation review reported 0 BLOCKER, 2 MAJOR, and 2 MINOR findings. All
four were remediated and the targeted failure-injection, nested-validation,
structured-error, stale-lookup, and concurrency negative controls passed. The
architecture audit has no unresolved Phase 2 BLOCKER/MAJOR findings. Allowlist,
whitespace, unchanged Phase 1, and unstaged-work checks passed.

See [Phase 2 implementation audit](audits/PHASE_2_IMPLEMENTATION_AUDIT.md) for the
exact public API, changed files, reusable backend-conformance suite, architectural
findings, commands, limitations, and actual baseline/final results. No staging
or commit is authorized or performed.

## PHASE 3 — Artifact Registry/content addressing
**Status:** TODO

ArtifactRef semantics + SHA256 verification + in-memory/filesystem registry.

## PHASE 4 — Storage abstractions
**Status:** TODO

Repository protocols independent of PostgreSQL. In-memory first.

## PHASE 5 — Federation contracts
**Status:** TODO

SourceAdapter + KnowledgeImporter interfaces with enforced responsibility separation.

## PHASE 6 — UniProt vertical slice
**Status:** TODO

UniProt accession -> Protein + ProteinSequence + provenance. Fixture-based Internet-free tests.

## PHASE 7 — RCSB/wwPDB vertical slice
**Status:** TODO

Structures, assemblies, chains, revision-aware artifacts.

## PHASE 8 — PDB CCD chemical vertical slice
**Status:** TODO

Non-polymer chemical identity preserving ChemicalIdentity/ChemicalSpecies distinction.

## PHASE 9 — Policy Engine + Knowledge Resolver
**Status:** TODO

Typed deterministic policies and four decision states. No generic DSL yet.

## PHASE 10 — Knowledge Compiler
**Status:** TODO

Resolved workflow context -> deterministic compilation model.

## PHASE 11 — Immutable KnowledgeBundle
**Status:** TODO

Filesystem bundle + manifest + contracts + provenance + hashed artifacts.

## PHASE 12 — SimForge BundleReader boundary
**Status:** TODO

Consume bundles without external federation. Do not refactor the whole DecisionEngine in same patch.

## PHASE 13 — Protein + ligand end-to-end MVP
**Status:** TODO

Target:

```bash
simforge knowledge prepare --protein P00533 --ligand gefitinib --out egfr_gefitinib_bundle/
simforge run --knowledge-bundle egfr_gefitinib_bundle/
```

Second command must work without Internet.

## PHASE 14 — MembraneTemplate and orientation
**Status:** TODO

Evidence-aware transmembrane workflows and explicit fallback templates.

## PHASE 15 — PostgreSQL production backend
**Status:** TODO

Mutable Global Knowledge Plane behind stable repository contracts.

## PHASE 16 — DuckDB/Parquet Data Plane
**Status:** TODO

Optimize large workflow-local collections and bundle export.

## PHASE 17 — HPC distribution/integration tests
**Status:** TODO

Validate staging, no-network runtime, fail-fast behavior and concurrent bundle consumption.

# Per-phase Codex protocol

1. READ specs/ADRs/AGENTS
2. INSPECT repository
3. REPORT affected files, risks, proposed changes
4. DO NOT EDIT until scope is coherent
5. IMPLEMENT only current phase
6. ADD tests
7. RUN focused tests
8. RUN relevant regressions
9. REVIEW diff as architecture auditor
10. FIX BLOCKER/MAJOR only
11. UPDATE this plan
12. COMMIT one conceptual change
