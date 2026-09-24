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
**Status:** READY FOR DESIGN

Objects: EntityId, ExternalIdentifier, Entity, Protein, ProteinSequence, ChemicalIdentity, ChemicalSpecies, ExperimentalStructure, ArtifactRef, Claim, Evidence, PolicyRef, Decision, DecisionContext, KnowledgeBundleManifest.

Allowed: `simforge/knowledge/domain/`, `tests/knowledge/domain/`, minimal exports.

Forbidden: HTTP, SQLAlchemy/PostgreSQL, SQLite/DuckDB, simulation engines, CLI, source-specific adapters.

Acceptance: external IDs are not PKs; Claim != Evidence; deterministic serialization; SHA256 artifact semantics; no persistence/network/engine imports.

---

## PHASE 2 — Universal Entity Registry
**Status:** TODO

Registry protocol + in-memory implementation + aliases/deprecation semantics.

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
