# ADR-003 — No Network Access from Compute Workers

**Status:** Accepted

## Decision
Compute workers may not call scientific APIs or query SSKI services. Missing mandatory knowledge causes a fail-fast runtime error.

## Consequences
All federation occurs before job submission.
