# ADR-005 — SourceAdapter / KnowledgeImporter Separation

**Status:** Accepted

## Decision
SourceAdapter owns transport/cache/integrity only. KnowledgeImporter owns scientific parsing/normalization/entity/claim emission only.

## Consequences
KnowledgeImporters are fixture-testable without network access.
