# Codex Prompt — Phase 1 Domain Design

Read:

- `AGENTS.md`
- `simforge/knowledge/AGENTS.md`
- `docs/sski/02_SCIENTIFIC_KNOWLEDGE_MODEL.md`
- `docs/sski/03_IDENTITY_AND_MAPPING.md`
- `docs/sski/04_EVIDENCE_AND_PROVENANCE.md`
- `docs/sski/05_KNOWLEDGE_BUNDLE_SPEC.md`
- relevant ADRs
- Phase 0 audit

Inspect existing SimForge model/test conventions.

## Scope

Design ONLY the pure SSKI domain model for:

- EntityId
- ExternalIdentifier
- Entity
- Protein
- ProteinSequence
- ChemicalIdentity
- ChemicalSpecies
- ExperimentalStructure
- ArtifactRef
- Claim
- Evidence
- PolicyRef
- Decision
- DecisionContext
- KnowledgeBundleManifest

## Constraints

- no HTTP
- no database dependency
- no SQLite/DuckDB
- no GROMACS/OpenMM/AMBER
- no CLI changes
- no SourceAdapter implementation
- no source-specific types
- external identifiers are not internal PKs
- Claim and Evidence remain distinct
- serialization deterministic
- ArtifactRef enforces SHA256 semantics

## First response

Do not edit files yet.

Return proposed files, model responsibilities, invariants, reuse/new decisions, test matrix and risks/collisions. Wait for an implementation instruction.
