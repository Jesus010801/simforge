# ADR-009 — Hard SSKI / SimForge Boundary

**Status:** Accepted

## Decision
SSKI terminates at KnowledgeBundle + DecisionContext.

SimForge owns preparation, parameterization, assembly, engine-independent SimulationPlan and execution.

## Consequences
SSKI cannot contain GROMACS/OpenMM/AMBER execution decisions.
