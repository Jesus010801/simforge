# simforge/cli.py — package entry point
# Re-exports the root cli so all commands are available under `simforge`.
# With editable install (pip install -e .) the project root is on sys.path,
# so `from cli import cli` resolves to the top-level cli.py.
from cli import cli as app  # noqa: E402

__all__ = ["app"]
