"""
Optional hydrogenation backends for ligand preparation.

Auto-selection priority: RDKit → Open Babel → ManualRequired.
All imports of rdkit and obabel are lazy/optional — this module never hard-depends
on either.

Public API:
  HydrogenationResult     — dataclass returned by every backend
  HydrogenationBackend    — base class (not a Protocol to keep Python 3.9 compat)
  RDKitHydrogenationBackend
  OpenBabelHydrogenationBackend
  ManualRequiredBackend
  select_backend(mode)    → HydrogenationBackend | None   (None iff mode=="none")
"""

from __future__ import annotations

import importlib.util
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


# ── Result dataclass ──────────────────────────────────────────────────────────

@dataclass
class HydrogenationResult:
    success: bool
    output_path: Optional[Path] = None
    backend_name: str = "none"      # "rdkit" | "openbabel" | "manual_required"
    n_hydrogen_atoms: int = 0
    confidence: float = 0.0
    warning: Optional[str] = None
    # Chemical-perception confidence metadata (see ligand.chemical_perception).
    # "high" | "medium" | "low" | "unknown". A successful hydrogenation call
    # does NOT by itself imply high confidence -- see decide_ligpargen_readiness.
    connectivity_confidence: str = "unknown"
    bond_order_confidence: str = "unknown"
    hydrogenation_confidence: str = "unknown"
    method: str = "none"
    formal_charge: Optional[int] = None


# ── Availability probes (module-level for easy mocking in tests) ──────────────

def _rdkit_available() -> bool:
    """Return True if rdkit is importable."""
    return importlib.util.find_spec("rdkit") is not None


def _obabel_available() -> bool:
    """Return True if the obabel executable is in PATH."""
    return shutil.which("obabel") is not None


# ── Shared H-counting helper ──────────────────────────────────────────────────

def _count_h_atoms(pdb_lines: list[str]) -> int:
    """Count hydrogen atoms in ATOM/HETATM lines using element + name heuristics."""
    count = 0
    for line in pdb_lines:
        record = line[:6].rstrip()
        if record not in ("ATOM", "HETATM"):
            continue
        # Element column (cols 76-78, 0-indexed 76:78) — most reliable
        if len(line) >= 78:
            element = line[76:78].strip().upper()
            if element in ("H", "D"):
                count += 1
                continue
        # Atom name (cols 12-16) fallback
        if len(line) >= 14:
            atom_name = line[12:16].strip()
            if atom_name and atom_name[0].upper() in ("H", "D"):
                count += 1
    return count


# ── Backend base class ────────────────────────────────────────────────────────

class HydrogenationBackend:
    """Abstract base for hydrogenation backends."""
    name: str = "base"

    def is_available(self) -> bool:
        return False

    def add_hydrogens(self, input_pdb: Path, output_pdb: Path) -> HydrogenationResult:
        raise NotImplementedError


# ── RDKit backend ─────────────────────────────────────────────────────────────

class RDKitHydrogenationBackend(HydrogenationBackend):
    """
    Hydrogenates via geometry-based bond-order perception
    (``ligand.chemical_perception.add_hydrogens_with_confidence``), NOT the
    legacy ``Chem.MolFromPDBFile()`` + ``Chem.AddHs()`` combination.

    That legacy combination trusts the PDB element column verbatim and does
    not perceive aromaticity/bond order from geometry, so it silently
    over-protonates molecules whenever the input has an element-column
    defect or non-trivial ring systems -- this was the root cause of a real
    mis-hydrogenation regression (an 18-carbon, 1-chlorine ligand with a
    "CL" atom mistyped to element "C" was hydrogenated to ~C19H34O3 instead
    of the correct C18H15ClO3). See ``ligand.chemical_perception`` for the
    full analysis. This backend never falls back to that legacy path --
    when geometry-based perception cannot converge confidently (typically:
    a heavy-atom-only pose with no explicit hydrogens), it fails explicitly
    so the ``auto`` cascade can try Open Babel or report manual_required,
    rather than silently emitting a wrong molecule.
    """
    name = "rdkit"

    def is_available(self) -> bool:
        return _rdkit_available()

    def add_hydrogens(self, input_pdb: Path, output_pdb: Path) -> HydrogenationResult:
        try:
            from ligand.chemical_perception import add_hydrogens_with_confidence
        except ImportError:
            return HydrogenationResult(
                success=False,
                backend_name=self.name,
                warning="RDKit is not installed.",
            )

        result = add_hydrogens_with_confidence(input_pdb, output_pdb)
        if not result.success:
            return HydrogenationResult(
                success=False,
                backend_name=self.name,
                warning=result.warning,
                connectivity_confidence=result.connectivity_confidence,
                bond_order_confidence=result.bond_order_confidence,
                method=result.method,
            )

        return HydrogenationResult(
            success=True,
            output_path=result.output_path,
            backend_name=self.name,
            n_hydrogen_atoms=result.n_hydrogen_atoms,
            confidence=1.0,
            connectivity_confidence=result.connectivity_confidence,
            bond_order_confidence=result.bond_order_confidence,
            hydrogenation_confidence=result.hydrogenation_confidence,
            method=result.method,
            formal_charge=result.formal_charge,
        )


# ── Open Babel backend ────────────────────────────────────────────────────────

class OpenBabelHydrogenationBackend(HydrogenationBackend):
    name = "openbabel"

    def is_available(self) -> bool:
        return _obabel_available()

    def add_hydrogens(self, input_pdb: Path, output_pdb: Path) -> HydrogenationResult:
        try:
            proc = subprocess.run(
                [
                    "obabel",
                    "-ipdb", str(input_pdb),
                    "-opdb",
                    "-O", str(output_pdb),
                    "-h",
                ],
                capture_output=True, text=True, timeout=30,
            )
        except subprocess.TimeoutExpired:
            return HydrogenationResult(
                success=False,
                backend_name=self.name,
                warning="Open Babel timed out while adding hydrogens (30 s limit).",
            )
        except Exception as exc:
            return HydrogenationResult(
                success=False,
                backend_name=self.name,
                warning=f"Open Babel raised an unexpected error: {exc}",
            )

        if proc.returncode != 0:
            return HydrogenationResult(
                success=False,
                backend_name=self.name,
                warning=(
                    f"Open Babel exited with code {proc.returncode}: "
                    f"{proc.stderr.strip() or '(no stderr output)'}"
                ),
            )

        if not output_pdb.exists():
            return HydrogenationResult(
                success=False,
                backend_name=self.name,
                warning=(
                    "Open Babel ran successfully but did not create the output file. "
                    "The ligand may have structural issues."
                ),
            )

        n_h = _count_h_atoms(output_pdb.read_text().splitlines())
        if n_h == 0:
            output_pdb.unlink(missing_ok=True)
            return HydrogenationResult(
                success=False,
                backend_name=self.name,
                warning=(
                    "Open Babel ran but the output still contains no hydrogen atoms. "
                    "The ligand may have an unusual structure — add hydrogens manually."
                ),
            )

        return HydrogenationResult(
            success=True,
            output_path=output_pdb,
            backend_name=self.name,
            n_hydrogen_atoms=n_h,
            confidence=0.9,
            # Open Babel performs its own valence/aromaticity perception, which
            # is generally more robust than the retired legacy RDKit PDB-flavor
            # path -- but it is not externally validated and has itself
            # produced an incorrect result for at least one real regression
            # ligand. "medium", not "high": on its own it cannot satisfy
            # decide_ligpargen_readiness()'s bar for ready_for_ligpargen.
            connectivity_confidence="medium",
            bond_order_confidence="medium",
            hydrogenation_confidence="medium",
            method="openbabel_perception",
        )


# ── Manual-required fallback ──────────────────────────────────────────────────

class ManualRequiredBackend(HydrogenationBackend):
    name = "manual_required"

    def is_available(self) -> bool:
        return True

    def add_hydrogens(self, input_pdb: Path, output_pdb: Path) -> HydrogenationResult:
        return HydrogenationResult(
            success=False,
            backend_name=self.name,
            warning=(
                "No automatic hydrogenation backend is available. "
                "Install RDKit (conda install -c conda-forge rdkit) or "
                "Open Babel (conda install -c conda-forge openbabel), or "
                "add hydrogens manually in Avogadro, PyMOL, or Discovery Studio."
            ),
        )


# ── Backend selection ─────────────────────────────────────────────────────────

def select_backend(mode: str) -> Optional[HydrogenationBackend]:
    """
    Return the best available backend for the given mode.

    Args:
        mode: "auto" | "rdkit" | "openbabel" | "none"

    Returns:
        A HydrogenationBackend, or None when mode is "none" (skip hydrogenation).

    Priority for "auto": RDKit → Open Babel → ManualRequired.
    """
    if mode == "none":
        return None

    if mode == "auto":
        candidates: list[HydrogenationBackend] = [
            RDKitHydrogenationBackend(),
            OpenBabelHydrogenationBackend(),
            ManualRequiredBackend(),
        ]
    elif mode == "rdkit":
        candidates = [RDKitHydrogenationBackend(), ManualRequiredBackend()]
    elif mode == "openbabel":
        candidates = [OpenBabelHydrogenationBackend(), ManualRequiredBackend()]
    else:
        candidates = [ManualRequiredBackend()]

    for backend in candidates:
        if backend.is_available():
            return backend

    return ManualRequiredBackend()
