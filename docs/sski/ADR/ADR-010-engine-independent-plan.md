# ADR-010 — Engine-Independent SimulationPlan

**Status:** Accepted

## Decision
Scientific knowledge resolution must not directly emit engine-specific configuration. SimForge consumes resolved context into an engine-independent SimulationPlan before backend translation.

## Consequences
Knowledge infrastructure evolves independently of individual simulation engines.
