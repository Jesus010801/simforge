# ADR-001 — Global Knowledge Plane vs Runtime Data Plane

**Status:** Accepted

## Context
Scientific knowledge is mutable and may require federation; HPC workers must execute reproducibly without external service dependencies.

## Decision
Separate a mutable Global Knowledge Plane from a workflow-scoped immutable runtime Data Plane. Knowledge is compiled before compute execution.

## Consequences
Workers cannot use the Global Plane directly. Required knowledge must be staged in a KnowledgeBundle.
