# ADR-008 — Policy-Driven Deterministic Resolution

**Status:** Accepted

## Decision
Workflow-scoped scientific conflicts are resolved through explicit versioned policies.

Allowed states: RESOLVED, RESOLVED_WITH_WARNING, DEFERRED, ABORT.

## Consequences
A decision must be reconstructible from its knowledge inputs and policy version.
