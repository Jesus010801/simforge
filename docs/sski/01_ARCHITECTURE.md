# SSKI Architecture Specification

## Canonical architecture

```text
                     MUTABLE SCIENTIFIC WORLD
                              │
             ┌────────────────┴─────────────────┐
             │                                  │
       Scientific sources                 Local knowledge
             │                                  │
             └────────────────┬─────────────────┘
                              ▼
                    SOURCE FEDERATION
                              │
                              ▼
                   KNOWLEDGE IMPORTERS
                              │
                              ▼
┌──────────────────────────────────────────────────────────────┐
│                 GLOBAL KNOWLEDGE PLANE                       │
│ Universal Entity Registry                                    │
│ Evidence Graph / Claims / Provenance                          │
│ Identity Resolution                                          │
│ Sequence / Chemical / Structure indices                      │
│ Artifact Registry                                            │
│ Knowledge Resolver                                           │
│ PostgreSQL + object store + Parquet                          │
└─────────────────────────────┬────────────────────────────────┘
                              │
                       Knowledge Compiler
                              │
                              ▼
┌──────────────────────────────────────────────────────────────┐
│                       CONTROL PLANE                          │
│ Workflow request                                             │
│ Policy Engine                                                │
│ Decision Engine                                              │
│ Resolver                                                     │
│ Artifact Fetcher                                             │
│ Bundle Builder                                               │
└─────────────────────────────┬────────────────────────────────┘
                              │
                    immutable KnowledgeBundle
                              │
                              ▼
┌──────────────────────────────────────────────────────────────┐
│                         DATA PLANE                           │
│ SQLite / DuckDB                                              │
│ Parquet                                                      │
│ mmCIF / SDF / FASTA / auxiliary artifacts                   │
│ DecisionContext                                              │
│ policy snapshot                                              │
│ hashes                                                       │
└─────────────────────────────┬────────────────────────────────┘
                              │
══════════════════════════════╪══════════════════════════════════
                 HARD SSKI / SIMFORGE BOUNDARY
══════════════════════════════╪══════════════════════════════════
                              ▼
┌──────────────────────────────────────────────────────────────┐
│                   SIMFORGE COMPUTE WORKER                    │
│ System Preparation                                           │
│ Parameterization                                             │
│ Assembly                                                     │
│ Engine-independent SimulationPlan                            │
│ Execution backend                                            │
└─────────────────────────────┬────────────────────────────────┘
                              │
                              ▼
                     RunDerivationLedger
                              │
                              ▼
                           Results
```

## Responsibilities

### Global Knowledge Plane

Responsible for canonical internal identities, external mappings, claims, evidence, source releases, provenance, artifact indexing and cross-source mappings. It is mutable because scientific knowledge changes.

### Control Plane

May query the Global Plane, trigger source federation before job submission, apply policies, select artifacts, detect insufficient knowledge and build an immutable bundle.

### Data Plane

Workflow-local representation distributed to compute environments.

- SQLite: identity lookup, mappings, compact relational metadata
- DuckDB/Parquet: analytical collections and larger local tables
- filesystem/object artifacts: mmCIF, SDF, FASTA, transforms

### Compute Worker

Consumes the bundle. It does not federate, access scientific APIs, reinterpret missing source knowledge or mutate the bundle. It may execute transformations explicitly authorized by a preparation contract.

## Knowledge vs derivation

```text
KnowledgeBundle
    ↓
protonation / reconstruction / topology / membrane assembly
    ↓
RunDerivationLedger
    ↓
SimulationPlan
```

The KnowledgeBundle records the resolved scientific input. The RunDerivationLedger records transformations produced by SimForge.

## Storage strategy

- **PostgreSQL:** mutable normalized knowledge and transactional state.
- **Object store/filesystem CAS:** large raw/normalized artifacts keyed by content hash.
- **Parquet:** analytical datasets, snapshots and candidate collections.
- **SQLite/DuckDB:** bundle-local access.
- **Graph database:** optional projection later; not canonical source of truth for MVP.

## Core invariants

- Domain semantics are independent of persistence.
- External IDs are aliases, not canonical SSKI primary keys.
- Every scientific assertion is distinguishable from its supporting evidence.
- Bundle compilation is deterministic for the same pinned inputs, policies and source versions.
- Runtime execution is independent of network availability.
