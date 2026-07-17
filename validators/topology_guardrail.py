# validators/topology_guardrail.py
"""
Guardrail preventing pdb2gmx from being run on mixed protein-lipid GRO files.

pdb2gmx must only run on:
  - Original protein PDB files  (correct chain/TER/terminus semantics)
  - Pure-protein GRO files      (single-chain, no lipid/solvent contamination)
  - Pure-lipid GRO files        (no protein terminus issues)

pdb2gmx must NEVER run on:
  - embed_in_bilayer/system.gro     (mixed protein + lipids, no TER records)
  - membrane_embedding/converged.gro (same)
  - Any GRO containing both protein residues and lipid/solvent residues
"""
from __future__ import annotations

from pathlib import Path

from core.topology_models import TopologyBuildError, is_mixed_system_gro

_EMBED_PATH_KEYWORDS = (
    "embed_in_bilayer",
    "membrane_embedding",
    "converged",
)


def check_pdb2gmx_input_safety(path: "str | Path", context: str = "") -> None:
    """
    Raise TopologyBuildError if pdb2gmx is about to run on a mixed/embedded file.

    Checks:
    1. Path contains known embed-step keywords (fast, no I/O)
    2. GRO file structurally contains both protein and membrane/solvent residues

    Safe inputs: .pdb files, pure-protein GRO, pure-lipid GRO.
    """
    path = Path(path)
    path_lower = str(path).lower()

    for kw in _EMBED_PATH_KEYWORDS:
        if kw in path_lower:
            raise TopologyBuildError(
                f"pdb2gmx must not run on an embedded protein-lipid coordinate file.\n"
                f"  Input: '{path}'\n"
                f"  Path contains '{kw}' — this file is from an embed/convergence step.\n"
                f"  Fix: run pdb2gmx on inputs/protein_1.pdb (generate_protein_topology),\n"
                f"  then use assemble_system_topology for the final combined topology."
            )

    if path.suffix == ".gro" and path.exists() and is_mixed_system_gro(path):
        raise TopologyBuildError(
            f"pdb2gmx must not run on a mixed protein-lipid GRO file.\n"
            f"  Input: '{path}'\n"
            f"  File contains both protein and membrane/solvent residues.\n"
            f"  Fix: run pdb2gmx on inputs/protein_1.pdb (generate_protein_topology)\n"
            f"  and assemble the combined topology with assemble_system_topology."
        )
