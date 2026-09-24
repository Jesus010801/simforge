# SSKI Package Agent Rules

These rules apply to all files under `simforge/knowledge/`.

Read the relevant SSKI specifications and ADRs before editing.

## Dependency boundary

`simforge/knowledge/` must remain simulation-engine independent.

Forbidden architectural dependencies include direct imports from:

- GROMACS-specific builders/executors
- OpenMM-specific runtime code
- AMBER-specific runtime code

The package may expose engine-neutral scientific context only.

## Layer responsibilities

### `domain/`

Pure scientific/domain models.

Must not depend on:

- HTTP/network libraries
- database libraries
- runtime executors
- source-specific API clients

### `identity/`

Internal/external identity mapping and resolution contracts.

Similarity may create candidate mappings but must not silently assert identity.

### `federation/`

Transport, caching, integrity, source release metadata.

No scientific decision logic.

### `importers/`

Scientific parsing and normalization.

No network access.

### `storage/`

Persistence protocols/backends.

Must not define scientific meaning.

### `policies/`

Versioned deterministic policy models/evaluation.

Must preserve:

- RESOLVED
- RESOLVED_WITH_WARNING
- DEFERRED
- ABORT

### `resolver/`

Combines knowledge candidates and policies into decisions.

### `compiler/`

Builds immutable KnowledgeBundles from resolved contexts.

### `runtime/`

Read-only bundle validation/access helpers only.

Must not perform federation.

## Invariants

- External IDs are never internal primary keys.
- Claim != Evidence.
- ChemicalIdentity != ChemicalSpecies != SimulationMicrostate.
- Protein != ProteinSequence != structure chain.
- PredictedStructure != ExperimentalStructure.
- Every bundle artifact is hash-addressable/verifiable.
- KnowledgeBundle is immutable after compilation.
- Derived runtime artifacts belong outside the original bundle.
