# Codex Prompt — Phase 1 Implementation

Implement only the approved Phase 1 design.

Allowed paths:

- `simforge/knowledge/domain/`
- `tests/knowledge/domain/`
- minimal `simforge/knowledge/__init__.py` exports if required

Forbidden:

- external API access
- persistence implementation
- SQLAlchemy/PostgreSQL
- SQLite/DuckDB
- simulation engine imports
- CLI modifications
- unrelated refactors

Requirements:

1. implement every approved invariant
2. add regression/architecture tests
3. run focused tests
4. run relevant pre-existing regression subset identified in Phase 0
5. inspect final diff for architecture violations

Report files changed, tests added/run, deviations and intentionally deferred work.

Do not begin Phase 2.
