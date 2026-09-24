# ADR-007 — KnowledgeSnapshot + RunDerivationLedger

**Status:** Accepted

## Decision
Use two provenance domains: KnowledgeSnapshot records what SSKI supplied/decided; RunDerivationLedger records what SimForge calculated/transformed.

## Consequences
Generated topologies and protonated structures are derived artifacts, not mutations of the original KnowledgeBundle.
