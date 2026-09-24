# ADR-002 — Immutable KnowledgeBundle

**Status:** Accepted

## Decision
KnowledgeBundle becomes immutable after successful compilation. Runtime-derived outputs are written outside the bundle and tracked by RunDerivationLedger.

## Consequences
Updated knowledge requires a new bundle, never in-place mutation.
