# SSKI — Scientific Knowledge Infrastructure for SimForge

**Specification status:** Draft normative specification  
**Architecture line:** SSKI 0.1

## Purpose

SSKI lets SimForge consume broad, current scientific knowledge without treating external databases as infallible truth and without requiring compute workers to query the Internet during execution.

> SSKI does not need to physically store all scientific knowledge. It must be able to identify it, locate it, normalize it, preserve its provenance, resolve it under explicit policies, remember what is necessary, and reconstruct exactly what was used.

```text
Mutable external knowledge
        ↓
SSKI
        ↓
KnowledgeBundle + DecisionContext
        ↓
SimForge
```

## Problem A — universal scientific representation

Represent and resolve:

- proteins and sequences
- isoforms and variants
- experimental and predicted structures
- chains, assemblies, residues and mappings
- chemical identities, species and simulation microstates
- ligands, cofactors, ions and solvents
- lipids and membrane contexts
- claims, evidence, provenance and contradictions
- artifacts and derived scientific artifacts

## Problem B — efficient large-scale distribution

Allow hundreds or thousands of workflows to execute without:

- direct Internet access from workers
- repeated API calls
- PostgreSQL bottlenecks
- mutable runtime knowledge
- hidden provenance
- queue stalls caused by federation latency

## Hard boundary

SSKI resolves **what the scientific context is**.

SimForge resolves **how to prepare, parameterize, assemble, plan and execute the simulation**.

SSKI MUST NOT decide engine-specific runtime parameters.

## Core components

- Universal Entity Registry
- Source Federation
- Knowledge Importers
- Scientific Knowledge Store
- Evidence/Provenance model
- Identity Resolver
- Knowledge Resolver
- Policy Engine
- Artifact Registry
- Knowledge Compiler
- KnowledgeBundle
- DecisionContext
- RunDerivationLedger

## Mutability model

| Component | Mutability |
|---|---|
| Global Knowledge Plane | Mutable/versioned |
| KnowledgeBundle | Immutable |
| RunDerivationLedger | Append-only |
| Scientific artifacts | Content-addressed |
| Simulation outputs | Content-addressed/versioned |

## Normative documents

1. `01_ARCHITECTURE.md`
2. `02_SCIENTIFIC_KNOWLEDGE_MODEL.md`
3. `03_IDENTITY_AND_MAPPING.md`
4. `04_EVIDENCE_AND_PROVENANCE.md`
5. `05_KNOWLEDGE_BUNDLE_SPEC.md`
6. `06_POLICY_ENGINE_SPEC.md`
7. `07_FEDERATION_AND_SOURCES.md`
8. `08_RUNTIME_AND_HPC_CONTRACT.md`
9. `09_IMPLEMENTATION_PLAN.md`
10. `ADR/`

Long-form design history is retained under `reference/` and is subordinate to these normative specifications.
