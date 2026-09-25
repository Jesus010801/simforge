"""Topological molecule partition — which atoms form one molecule (Phase 13.5).

Authoritative source: the run input (.tpr) molecule blocks, read with
``gmx dump -s``.  Molecules are contiguous atom ranges in topology order:
molblocks in order, ``#molecules`` copies of their moltype, each moltype with
``atom (N)`` atoms.  Nothing is inferred from chain ids or residue numbering.

This is topology evidence, not a biological annotation: a Receptor component
may be one molecule or four (an oligomer); a Membrane holds thousands.
"""
from __future__ import annotations

import bisect
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

_MOLBLOCK = re.compile(r"^\s+molblock \((\d+)\):")
_MOLTYPE_REF = re.compile(r"^\s+moltype\s+=\s+(\d+)\s+\"(.*)\"")
_NMOL = re.compile(r"^\s+#molecules\s+=\s+(\d+)")
_MOLTYPE = re.compile(r"^\s+moltype \((\d+)\):")
_ATOMS = re.compile(r"^\s+atom \((\d+)\):")
_NAME = re.compile(r'^\s+name="(.*)"')

_CACHE: dict[str, "MoleculePartition"] = {}


class MoleculePartitionError(RuntimeError):
    """The molecule structure could not be established from the topology."""


@dataclass
class MoleculePartition:
    """Molecule k spans 1-based atoms ``starts[k] .. starts[k+1]-1``."""
    starts: list[int]                          # len = n_molecules + 1 (1-based, exclusive end)
    moltypes: list[str]                        # moltype name per molecule
    source: dict = field(default_factory=dict)  # topology path + digest + method

    @property
    def n_molecules(self) -> int:
        return len(self.starts) - 1

    @property
    def n_atoms(self) -> int:
        return self.starts[-1] - 1

    def molecule_of(self, atom_id: int) -> int:
        if not 1 <= atom_id <= self.n_atoms:
            raise ValueError(f"atom {atom_id} outside 1..{self.n_atoms}")
        return bisect.bisect_right(self.starts, atom_id) - 1

    def atoms_of(self, molecule: int) -> range:
        return range(self.starts[molecule], self.starts[molecule + 1])

    def select(self, atom_ids: Iterable[int]) -> "SelectionMoleculePartition":
        ids = sorted(set(int(a) for a in atom_ids))
        per: dict[int, int] = {}
        for a in ids:
            m = self.molecule_of(a)
            per[m] = per.get(m, 0) + 1
        mols = sorted(per)
        complete = all(per[m] == len(self.atoms_of(m)) for m in mols)
        return SelectionMoleculePartition(
            n_atoms=len(ids), molecules=mols, atoms_per_molecule=[per[m] for m in mols],
            moltypes=sorted({self.moltypes[m] for m in mols}), complete_molecules=complete,
            source_identity=self.source.get("digest"))

    def closure(self, atom_ids: Iterable[int]) -> list[int]:
        """Every atom of every molecule the selection touches (clusterable unit)."""
        out: list[int] = []
        for m in sorted({self.molecule_of(int(a)) for a in atom_ids}):
            out.extend(self.atoms_of(m))
        return out


@dataclass
class SelectionMoleculePartition:
    n_atoms: int
    molecules: list[int]                       # 0-based molecule indices touched
    atoms_per_molecule: list[int]
    moltypes: list[str]
    complete_molecules: bool                   # every touched molecule fully selected
    source_identity: Optional[str] = None

    @property
    def n_molecules(self) -> int:
        return len(self.molecules)

    @property
    def multi_molecule(self) -> bool:
        return self.n_molecules > 1

    def evidence(self) -> dict:
        return {"n_molecules": self.n_molecules, "moltypes": self.moltypes,
                "complete_molecules": self.complete_molecules,
                "topology": self.source_identity}


def parse_dump(lines: Iterable[str]) -> tuple[list[tuple[int, int]], dict[int, tuple[str, int]]]:
    """``gmx dump -s`` text → ``([(moltype, n_molecules) per molblock], {moltype: (name, n_atoms)})``."""
    blocks: list[list[int]] = []
    moltypes: dict[int, list] = {}
    cur_block = cur_type = None
    for line in lines:
        m = _MOLBLOCK.match(line)
        if m:
            blocks.append([None, None]); cur_block = len(blocks) - 1; cur_type = None
            continue
        m = _MOLTYPE.match(line)
        if m:
            cur_type = int(m.group(1)); moltypes[cur_type] = [None, None]; cur_block = None
            continue
        if cur_block is not None:
            m = _MOLTYPE_REF.match(line)
            if m:
                blocks[cur_block][0] = int(m.group(1))
                continue
            m = _NMOL.match(line)
            if m:
                blocks[cur_block][1] = int(m.group(1))
                continue
        if cur_type is not None:
            m = _NAME.match(line)
            if m and moltypes[cur_type][0] is None:
                moltypes[cur_type][0] = m.group(1)
                continue
            m = _ATOMS.match(line)
            if m and moltypes[cur_type][1] is None:
                moltypes[cur_type][1] = int(m.group(1))
    if not blocks or any(b[0] is None or b[1] is None for b in blocks):
        raise MoleculePartitionError("no complete molblock records in the topology dump")
    for t, (name, n) in moltypes.items():
        if n is None:
            raise MoleculePartitionError(f"moltype {t} ({name}) has no atom count")
    return [(b[0], b[1]) for b in blocks], {t: (v[0], v[1]) for t, v in moltypes.items()}


def partition_from_dump(lines: Iterable[str], source: Optional[dict] = None) -> MoleculePartition:
    blocks, types = parse_dump(lines)
    starts, names, pos = [1], [], 1
    for t, n in blocks:
        if t not in types:
            raise MoleculePartitionError(f"molblock refers to unknown moltype {t}")
        name, size = types[t]
        for _ in range(n):
            pos += size
            starts.append(pos)
            names.append(name)
    return MoleculePartition(starts, names, dict(source or {}))


def molecule_partition(topology_path: Optional[str], gmx: str = "gmx"
                       ) -> Optional[MoleculePartition]:
    """The molecule partition of a .tpr (cached by content digest); None when
    no .tpr is available (a .gro/.pdb carries no molecule blocks — no guess)."""
    if not topology_path or not str(topology_path).lower().endswith(".tpr"):
        return None
    p = Path(topology_path)
    if not p.is_file():
        return None
    from analysis.campaign.fingerprint import fingerprint_file
    digest = fingerprint_file(p).digest
    if digest in _CACHE:
        return _CACHE[digest]
    proc = subprocess.Popen([gmx, "dump", "-s", str(p)], stdout=subprocess.PIPE,
                            stderr=subprocess.DEVNULL, text=True)
    try:
        part = partition_from_dump(proc.stdout, {"path": str(p), "digest": digest,
                                                 "method": "gmx dump -s (molblocks)"})
    finally:
        proc.stdout.close()
        proc.wait()
    if proc.returncode != 0:
        raise MoleculePartitionError(f"gmx dump -s {p.name} failed (rc={proc.returncode})")
    _CACHE[digest] = part
    return part
