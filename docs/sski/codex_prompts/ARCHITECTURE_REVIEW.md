# Codex Prompt — SSKI Architecture Review

Review the current diff only. Do NOT modify files.

Audit against root `AGENTS.md`, `simforge/knowledge/AGENTS.md`, relevant specs, ADRs and active phase scope.

Look for:

- domain/persistence coupling
- scientific assumptions encoded as facts
- external IDs used as primary keys
- missing provenance
- Claim/Evidence collapse
- mutable KnowledgeBundle behavior
- network access outside SourceAdapters
- scientific parsing inside SourceAdapters
- network access inside KnowledgeImporters
- GROMACS/OpenMM/AMBER coupling in SSKI
- silent fallback behavior
- missing regression tests
- premature abstractions
- unrelated edits

Classify BLOCKER / MAJOR / MINOR. Do not recommend unrelated cleanup.
