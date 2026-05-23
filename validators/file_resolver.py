# validators/file_resolver.py
"""
Smart structure file resolver.

When a component path is a directory, searches for candidate structures
and either auto-selects a single match (with warning) or raises a
StructureFileError with actionable suggestions.

Supported formats: .pdb  .gro  .cif
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional


STRUCTURE_EXTENSIONS = (".pdb", ".gro", ".cif")


class StructureFileError(Exception):
    """Raised when a component path cannot be resolved to a structure file."""

    def __init__(
        self,
        message: str,
        component_id: str = "",
        candidates: list[Path] | None = None,
        directory: Path | None = None,
    ):
        super().__init__(message)
        self.component_id = component_id
        self.candidates   = candidates or []
        self.directory    = directory

    def rich_lines(self) -> list[str]:
        """Return Rich-formatted lines for display in the CLI."""
        lines = [
            f"[red]✗[/red] [bold]{self.component_id}[/bold] expects a structure file, "
            "but received a [yellow]directory[/yellow].",
        ]
        if self.directory:
            lines.append(f"  Directory: [dim]{self.directory}[/dim]")
        if self.candidates:
            lines.append("")
            lines.append("  [bold]Detected candidate structures:[/bold]")
            for c in self.candidates:
                lines.append(f"    · [cyan]{c.name}[/cyan]")
            lines.append("")
            lines.append(f"  [dim]Did you mean:[/dim]  [green]{self.candidates[0]}[/green]")
        elif self.directory:
            lines.append("")
            lines.append(
                f"  [dim]No structure files (*{', *'.join(STRUCTURE_EXTENSIONS)}) "
                "found in that directory.[/dim]"
            )
        return lines


def _find_candidates(directory: Path) -> list[Path]:
    candidates: list[Path] = []
    for ext in STRUCTURE_EXTENSIONS:
        candidates.extend(sorted(directory.glob(f"*{ext}")))
    return candidates


def resolve_structure_file(
    raw_path: str | Path,
    component_id: str = "component",
) -> tuple[Path, Optional[str]]:
    """
    Resolve a raw path from a YAML component to a concrete structure file.

    Returns:
        (resolved_path, warning_message | None)
        warning_message is non-None only when the path was auto-resolved
        from a directory containing a single candidate.

    Raises:
        StructureFileError: if path is a directory with 0 or 2+ candidates,
                            or if path does not exist at all.
    """
    path = Path(raw_path)

    if not path.exists():
        raise FileNotFoundError(f"{component_id}: file not found: {path}")

    if not path.is_dir():
        return path, None

    candidates = _find_candidates(path)

    if not candidates:
        raise StructureFileError(
            f"{component_id}: no structure files in directory {path}",
            component_id=component_id,
            candidates=[],
            directory=path,
        )

    if len(candidates) == 1:
        warning = (
            f"{component_id}: path is a directory — "
            f"auto-selected sole candidate: {candidates[0].name}"
        )
        return candidates[0], warning

    raise StructureFileError(
        f"{component_id}: ambiguous — {len(candidates)} structure files found in {path}",
        component_id=component_id,
        candidates=candidates,
        directory=path,
    )
