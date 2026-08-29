"""Lightweight structure/topology parsers for component inference.

Reuses ``utils.gro_parser`` where possible.  ``utils.gro_parser`` has no chain
concept (``.gro`` carries none) so a minimal PDB reader is added here purely to
recover chain IDs and per-chain residue sequences.  No Biopython/MDAnalysis
dependency — consistent with the rest of SimForge.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from analysis.campaign.structure.residue_classes import classify_resname


@dataclass
class ChainInfo:
    chain_id: str
    residues: list[tuple[int, str]] = field(default_factory=list)  # (resSeq, resName) ordered, unique
    atom_count: int = 0
    has_ter: bool = False

    @property
    def residue_count(self) -> int:
        return len(self.residues)

    @property
    def resname_set(self) -> set[str]:
        return {r for _, r in self.residues}

    def class_histogram(self) -> dict[str, int]:
        hist: dict[str, int] = {}
        for _, rn in self.residues:
            k = classify_resname(rn)
            hist[k] = hist.get(k, 0) + 1
        return hist

    def dominant_class(self) -> str:
        hist = self.class_histogram()
        return max(hist, key=hist.get) if hist else "other"

    def polymer_fraction(self) -> float:
        if not self.residues:
            return 0.0
        hist = self.class_histogram()
        poly = hist.get("amino_acid", 0) + hist.get("nucleotide", 0)
        return poly / len(self.residues)


@dataclass
class StructureModel:
    path: str
    source_format: str                       # "pdb" | "gro"
    chains: list[ChainInfo] = field(default_factory=list)
    atom_count: int = 0
    residue_class_histogram: dict[str, int] = field(default_factory=dict)
    has_chain_ids: bool = False
    warnings: list[str] = field(default_factory=list)

    def resname_set(self) -> set[str]:
        out: set[str] = set()
        for c in self.chains:
            out |= c.resname_set
        return out


_PDB_ATOM_RE = re.compile(r"^(ATOM|HETATM)")


def parse_pdb_structure(path: Path) -> StructureModel:
    """Parse chain IDs and per-chain residue sequences from a PDB file."""
    p = Path(path)
    model = StructureModel(path=str(p.resolve()), source_format="pdb")
    chains: dict[str, ChainInfo] = {}
    seen_res: dict[tuple[str, int, str], bool] = {}
    stop_after_model_1 = False

    for line in p.read_text(errors="replace").splitlines():
        rec = line[:6].strip()
        if rec == "ENDMDL":
            stop_after_model_1 = True
            break
        if rec == "TER":
            cid = line[21:22].strip() or " "
            if cid in chains:
                chains[cid].has_ter = True
            continue
        if not _PDB_ATOM_RE.match(line):
            continue
        try:
            chain_id = line[21:22].strip() or "_"
            resseq = int(line[22:26])
            resname = line[17:20].strip()
        except (ValueError, IndexError):
            continue
        ci = chains.setdefault(chain_id, ChainInfo(chain_id=chain_id))
        ci.atom_count += 1
        model.atom_count += 1
        key = (chain_id, resseq, resname)
        if key not in seen_res:
            seen_res[key] = True
            ci.residues.append((resseq, resname))

    model.chains = [chains[k] for k in sorted(chains)]
    model.has_chain_ids = any(c.chain_id not in ("_", " ") for c in model.chains)
    for c in model.chains:
        for k, v in c.class_histogram().items():
            model.residue_class_histogram[k] = model.residue_class_histogram.get(k, 0) + v
    if stop_after_model_1:
        model.warnings.append("multi-model PDB: only MODEL 1 parsed")
    if not model.chains:
        model.warnings.append("no ATOM/HETATM records parsed")
    return model


def parse_gro_structure(path: Path) -> StructureModel:
    """Parse residues from a GRO file (no chain IDs available)."""
    from utils.gro_parser import parse_gro

    p = Path(path)
    model = StructureModel(path=str(p.resolve()), source_format="gro")
    try:
        gro = parse_gro(p)
    except (ValueError, OSError) as exc:
        model.warnings.append(f"gro parse failed: {exc}")
        return model

    model.atom_count = gro.atom_count
    # GRO carries no chain IDs.  Reconstruct *molecule-like* segments: start a new
    # segment when the residue class flips between polymer and non-polymer, or
    # when the residue number resets/decreases (a new moleculetype in GROMACS).
    segments: list[ChainInfo] = []
    cur: Optional[ChainInfo] = None
    prev_resnum: Optional[int] = None
    prev_is_poly: Optional[bool] = None
    seen_res: set[int] = set()
    poly_letters = iter("ABCDEFGHIJKLMNOPQRSTUVWXYZ")

    for atom in gro.atoms:
        is_poly = classify_resname(atom.residue_name) in ("amino_acid", "nucleotide")
        reset = prev_resnum is not None and atom.residue_number < prev_resnum - 1
        flip = prev_is_poly is not None and is_poly != prev_is_poly
        if cur is None or reset or flip:
            cid = next(poly_letters, "_") if is_poly else "_"
            cur = ChainInfo(chain_id=cid)
            segments.append(cur)
            seen_res = set()
        cur.atom_count += 1
        if atom.residue_number not in seen_res:
            seen_res.add(atom.residue_number)
            cur.residues.append((atom.residue_number, atom.residue_name))
        prev_resnum = atom.residue_number
        prev_is_poly = is_poly

    # merge all non-polymer segments back into class buckets on one "_" chain
    poly = [s for s in segments if s.chain_id != "_"]
    nonpoly = ChainInfo(chain_id="_")
    for s in segments:
        if s.chain_id == "_":
            nonpoly.residues.extend(s.residues)
            nonpoly.atom_count += s.atom_count
    model.chains = poly + ([nonpoly] if nonpoly.residues else [])
    model.has_chain_ids = False
    for c in model.chains:
        for k, v in c.class_histogram().items():
            model.residue_class_histogram[k] = model.residue_class_histogram.get(k, 0) + v
    model.warnings.append(
        f"GRO has no chain IDs; reconstructed {len(poly)} polymer segment(s) "
        f"from residue-number/class boundaries"
    )
    return model


def parse_structure(path: Path) -> StructureModel:
    p = Path(path)
    if p.suffix.lower() == ".pdb":
        return parse_pdb_structure(p)
    if p.suffix.lower() == ".gro":
        return parse_gro_structure(p)
    raise ValueError(f"unsupported structure format: {p.suffix}")


# ── Topology [ molecules ] ────────────────────────────────────────────────────

@dataclass
class TopologyMolecules:
    path: str
    molecules: list[tuple[str, int]] = field(default_factory=list)  # (name, count) in order
    warnings: list[str] = field(default_factory=list)

    def names(self) -> list[str]:
        return [n for n, _ in self.molecules]


def parse_top_molecules(path: Path) -> TopologyMolecules:
    """Parse the ``[ molecules ]`` section of a ``.top`` file.

    Section-scoped (mirrors ``core.system_validate._parse_molecules``); does not
    resolve ``#include`` — only the literal ``[ molecules ]`` block matters for
    component naming/counts.
    """
    p = Path(path)
    tm = TopologyMolecules(path=str(p.resolve()))
    try:
        text = p.read_text(errors="replace")
    except OSError as exc:
        tm.warnings.append(f"cannot read topology: {exc}")
        return tm

    in_block = False
    for raw in text.splitlines():
        line = raw.split(";", 1)[0].strip()
        if not line:
            continue
        if line.startswith("["):
            in_block = line.strip("[] ").lower() == "molecules"
            continue
        if not in_block:
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                tm.molecules.append((parts[0], int(parts[1])))
            except ValueError:
                tm.warnings.append(f"unparseable [ molecules ] line: {line!r}")
    if not tm.molecules:
        tm.warnings.append("no [ molecules ] section found")
    return tm
