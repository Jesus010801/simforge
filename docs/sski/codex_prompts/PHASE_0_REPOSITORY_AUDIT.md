# Codex Prompt — Phase 0 Repository Audit

Read first:

- `AGENTS.md`
- `simforge/knowledge/AGENTS.md`
- `docs/sski/00_OVERVIEW.md`
- `docs/sski/01_ARCHITECTURE.md`
- `docs/sski/09_IMPLEMENTATION_PLAN.md`
- all accepted ADRs in `docs/sski/ADR/`

Then inspect the current repository, especially existing core models, decision engine, structural annotation, semantic inference, MD knowledge modules, manifest/provenance implementations, CLI organization and test conventions.

## Task

Do NOT implement SSKI production code.

Produce a repository-aware audit answering:

1. Which current modules overlap conceptually with SSKI?
2. Which current models should be reused, adapted or isolated?
3. Which dependencies could violate the hard SSKI/SimForge boundary?
4. What is the safest package placement?
5. What exact files should Phase 1 be allowed to create/edit?
6. Which regression suites must Phase 1 run?
7. Are there naming collisions with `Entity`, `Artifact`, `Decision`, `Manifest`, etc.?
8. Which assumptions in the SSKI docs conflict with the actual repository?

Classify findings as BLOCKER, MAJOR or MINOR.

Do not make unrelated refactor suggestions. End with a proposed Phase 1 patch scope only.
