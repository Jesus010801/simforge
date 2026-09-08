"""Read original GROMACS indices; IDs are zero-based header positions."""
from dataclasses import dataclass
from pathlib import Path
import re


@dataclass(frozen=True)
class IndexGroup:
    name: str
    group_id: int
    atom_ids: tuple[int, ...]


def parse_index_groups(index_file: str | Path) -> list[IndexGroup]:
    groups = []
    name = None
    atoms = []
    for line in Path(index_file).read_text().splitlines():
        line = line.split(';', 1)[0].strip()
        if not line:
            continue
        match = re.fullmatch(r'\[\s*(.*?)\s*\]', line)
        if match:
            if not match[1] or "[" in match[1] or "]" in match[1]:
                raise ValueError("index headers require a nonempty group name without brackets")
            if name is not None:
                groups.append(IndexGroup(name, len(groups), tuple(atoms)))
            name, atoms = match[1], []
        else:
            if name is None:
                raise ValueError('atom list before first index header')
            values = [int(x) for x in line.split()]
            if any(x <= 0 for x in values):
                raise ValueError('index atom IDs must be positive')
            atoms.extend(values)
    if name is not None:
        groups.append(IndexGroup(name, len(groups), tuple(atoms)))
    if not groups:
        raise ValueError("index contains no group headers")
    return groups


def resolve_index_group(index_file: str | Path, group_name: str,
                        aliases: tuple[str, ...] = ()) -> IndexGroup | None:
    """Exact match first; explicit aliases only. Duplicate matches are errors."""
    groups = parse_index_groups(index_file)
    hits = [g for g in groups if g.name == group_name]
    if not hits:
        hits = [g for g in groups if g.name in aliases]
    if len(hits) > 1:
        raise ValueError(f'ambiguous index group {group_name}: {[g.name for g in hits]}')
    return hits[0] if hits else None
