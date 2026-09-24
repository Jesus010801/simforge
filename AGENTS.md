# SimForge Agent Instructions

## General repository rules

1. Inspect relevant existing code and tests before editing.
2. Preserve backward compatibility unless the task explicitly authorizes a breaking change.
3. Do not modify unrelated files.
4. One conceptual change per patch/commit.
5. Every behavior change requires regression tests.
6. Do not silently redesign architecture while implementing a feature.
7. Prefer explicit typed contracts over implicit conventions.
8. Report files changed, tests executed, failures, assumptions, and intentionally deferred work.

# SSKI — mandatory architectural rules

Before modifying SSKI, read:

- `docs/sski/00_OVERVIEW.md`
- `docs/sski/01_ARCHITECTURE.md`
- the specification relevant to the active phase
- `docs/sski/09_IMPLEMENTATION_PLAN.md`
- all accepted ADRs under `docs/sski/ADR/`

## SSKI invariants

- SSKI transforms mutable external scientific knowledge into reproducible, immutable scientific contexts consumable by SimForge.
- External database identifiers MUST NOT be internal primary keys.
- Scientific claims MUST remain distinct from evidence supporting them.
- Predictions/inferences MUST NOT be silently promoted to factual assertions.
- External information included in a decision MUST retain provenance.
- Artifacts participating in a KnowledgeBundle MUST be cryptographically hashed.
- KnowledgeBundle is immutable after compilation.
- Compute workers MUST NOT contact external scientific APIs.
- Compute workers MUST NOT query the Control Plane or Global Knowledge Plane.
- Compute workers MUST NOT mutate a KnowledgeBundle.
- Missing required bundle data MUST fail fast.
- Source adapters handle transport/cache/integrity, not scientific interpretation.
- Knowledge importers handle parsing/normalization/claim emission, not networking.
- `simforge/knowledge/` MUST remain simulation-engine independent.
- SSKI MUST NOT encode GROMACS/OpenMM/AMBER-specific execution behavior.

Hard boundary:

```text
SSKI -> KnowledgeBundle + DecisionContext -> SimForge
```

## Development process

For each phase:

1. inspect
2. propose
3. review
4. implement only approved scope
5. test
6. architecture audit
7. fix BLOCKER/MAJOR findings
8. commit
9. update `docs/sski/09_IMPLEMENTATION_PLAN.md`

Never implement multiple future phases merely because they appear straightforward.
