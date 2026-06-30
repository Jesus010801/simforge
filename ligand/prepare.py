"""
Ligand extraction and preparation from a docked protein-ligand complex PDB.

Steps:
  1. Detect HETATM residue names (auto or user-supplied).
  2. Extract ligand HETATM records to ligand_for_ligpargen.pdb.
  3. Extract ATOM records (protein) to protein_only.pdb.
  4. Assess hydrogen completeness:
       has_hydrogens        = n_hydrogen_atoms > 0
       hydrogenation_status = "complete" | "incomplete" | "missing" | "unknown"
       hydrogenation_complete = True only when status is "complete"
  5. Trigger hydrogenation when status is "missing" or "incomplete":
       auto (default): RDKit -> Open Babel -> manual_required
       rdkit / openbabel: force that backend
       none:  skip; block if status is not "complete"
  6. Estimate formal charge via RDKit (warns if unavailable).
  7. Write ligand_report.yaml.

Public API:
  detect_ligand_residues(pdb_path) -> list[str]
  extract_ligand_from_complex(
      complex_pdb, ligand_resname, out_dir,
      hydrogenation_mode="auto",
  ) -> LigandPrepareResult
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import yaml

from ligand.hydrogenation import select_backend

_STANDARD_RESIDUES: frozenset[str] = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL",
    "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "CYX", "CYM",
    "ACE", "NME", "NHE",
    "HOH", "WAT", "SOL", "TIP", "SPC",
    "NA", "CL", "K", "MG", "CA", "ZN", "FE", "MN", "CU",
    "NA+", "CL-", "ION",
})

# If H/heavy_atom ratio falls below this, the molecule is classified "incomplete"
# without needing RDKit.  Catches cases like 28 heavy + 1 H (ratio 0.036) while
# leaving ambiguous small-molecule fixtures (ratio 0.5) as "unknown".
_INCOMPLETE_RATIO_THRESHOLD = 0.15


@dataclass
class LigandPrepareResult:
    success: bool
    ligand_pdb: Optional[Path] = None
    protein_pdb: Optional[Path] = None
    report_path: Optional[Path] = None
    ligand_resname: str = ""
    ligand_atom_count: int = 0
    protein_atom_count: int = 0
    formal_charge: Optional[int] = None
    charge_estimated: bool = False
    # Hydrogen-handling fields
    has_hydrogens: bool = False                  # n_hydrogen_atoms > 0
    n_hydrogen_atoms: int = 0
    n_heavy_atoms: int = 0
    hydrogenation_status: str = "unknown"        # "complete"|"incomplete"|"missing"|"unknown"
    hydrogenation_complete: bool = False         # True iff status == "complete"
    hydrogenation_required: bool = False         # True iff status in (missing,incomplete) and not performed
    hydrogenation_backend: str = "none"          # "none"|"rdkit"|"openbabel"|"manual_required"|"skipped"
    hydrogenation_performed: bool = False
    hydrogenated_ligand_pdb: Optional[Path] = None
    hydrogenation_confidence: float = 0.0
    recommended_ligpargen_input: Optional[Path] = None
    warnings: list[str] = field(default_factory=list)
    error: Optional[str] = None


# ── Public API ────────────────────────────────────────────────────────────────

def detect_ligand_residues(pdb_path: str | Path) -> list[str]:
    """Return sorted unique ligand residue names from HETATM records."""
    path = Path(pdb_path)
    hetatm_resnames: set[str] = set()
    atom_nonstd: set[str] = set()

    for line in path.read_text().splitlines():
        record = line[:6].rstrip()
        resname = line[17:20].strip() if len(line) >= 20 else ""
        if not resname or resname in _STANDARD_RESIDUES:
            continue
        if record == "HETATM":
            hetatm_resnames.add(resname)
        elif record == "ATOM":
            atom_nonstd.add(resname)

    return sorted(hetatm_resnames or atom_nonstd)


def extract_ligand_from_complex(
    complex_pdb: str | Path,
    ligand_resname: str,
    out_dir: str | Path,
    hydrogenation_mode: str = "auto",
) -> LigandPrepareResult:
    """
    Split a docked complex PDB into protein-only and ligand-only files.

    Hydrogen completeness is assessed and auto-hydrogenation is triggered when
    the status is "missing" or "incomplete".  Use ``hydrogenation_mode="none"``
    to skip hydrogenation; the call will still block downstream if the ligand is
    not provably complete.

    Args:
        complex_pdb:        Path to the docked complex PDB.
        ligand_resname:     Residue name of the ligand (e.g. "E20").
        out_dir:            Directory for output files.
        hydrogenation_mode: "auto" | "rdkit" | "openbabel" | "none"
    """
    complex_pdb = Path(complex_pdb)
    out_dir = Path(out_dir)
    ligand_resname = ligand_resname.strip().upper()

    if not complex_pdb.exists():
        return LigandPrepareResult(
            success=False,
            error=f"Complex PDB not found: {complex_pdb}",
        )

    out_dir.mkdir(parents=True, exist_ok=True)

    protein_lines: list[str] = []  # ATOM, non-ligand HETATM, and TER records
    protein_atom_count: int = 0   # only ATOM/HETATM (not TER)
    ligand_lines: list[str] = []
    all_conect: list[str] = []
    warnings: list[str] = []

    for line in complex_pdb.read_text().splitlines():
        record = line[:6].rstrip()
        resname = line[17:20].strip() if len(line) >= 20 else ""

        if record == "ATOM":
            protein_lines.append(line)
            protein_atom_count += 1
        elif record == "HETATM":
            if resname == ligand_resname:
                ligand_lines.append(line)
            else:
                protein_lines.append(line)
                protein_atom_count += 1
        elif record == "TER":
            # Preserve TER records so pdb2gmx can recognise chain boundaries
            # in multichain proteins/dimers.
            protein_lines.append(line)
        elif record == "CONECT":
            all_conect.append(line)

    if not ligand_lines:
        candidates = detect_ligand_residues(complex_pdb)
        return LigandPrepareResult(
            success=False,
            ligand_resname=ligand_resname,
            error=(
                f"No records found for residue '{ligand_resname}' in {complex_pdb.name}. "
                f"Detected ligand residues: {candidates or ['(none)']}"
            ),
        )

    ligand_serials = _parse_serials(ligand_lines)
    conect_lines = [
        ln for ln in all_conect if _conect_first_serial(ln) in ligand_serials
    ]

    # ── Write raw ligand PDB ──────────────────────────────────────────────────
    lig_pdb = out_dir / "ligand_for_ligpargen.pdb"
    parts = ligand_lines[:]
    if conect_lines:
        parts += conect_lines
    parts.append("END")
    lig_pdb.write_text("\n".join(parts) + "\n")

    # ── Write protein PDB ─────────────────────────────────────────────────────
    prot_pdb = out_dir / "protein_only.pdb"
    # If the source PDB already has TER records (multichain), preserve them as-is.
    # Otherwise append a TER so pdb2gmx correctly terminates the single chain.
    has_ter = any(l[:6].rstrip() == "TER" for l in protein_lines)
    end_suffix = "\nEND\n" if has_ter else "\nTER\nEND\n"
    prot_pdb.write_text("\n".join(protein_lines) + end_suffix)

    # ── Hydrogen completeness assessment ──────────────────────────────────────
    n_h_atoms = _count_h_atoms(ligand_lines)
    n_heavy_atoms = _count_heavy_atoms(ligand_lines)
    has_h = n_h_atoms > 0

    hydrogenation_status, hydrogenation_complete = _determine_hydrogenation_status(
        ligand_lines, n_h_atoms, n_heavy_atoms
    )

    # ── Hydrogenation ─────────────────────────────────────────────────────────
    needs_hydrogenation = hydrogenation_status in ("missing", "incomplete")
    hydrogenated_pdb: Optional[Path] = None
    hydrogenation_performed = False
    hydrogenation_backend = "none"
    hydrogenation_confidence = 1.0 if hydrogenation_complete else 0.0

    if not hydrogenation_complete:
        if hydrogenation_mode == "none":
            # User explicitly chose not to hydrogenate
            hydrogenation_backend = "skipped"
            hydrogenation_confidence = 0.0
            if hydrogenation_status == "missing":
                warnings.append(
                    "Hydrogenation was skipped (--hydrogenation none). "
                    "Ligand has no explicit hydrogen atoms. "
                    "Do not submit ligand_for_ligpargen.pdb to LigParGen without hydrogens -- "
                    "the force-field parameters will be incorrect."
                )
            elif hydrogenation_status == "incomplete":
                warnings.append(
                    f"Hydrogenation was skipped (--hydrogenation none). "
                    f"Ligand has {n_h_atoms} H atom(s) for {n_heavy_atoms} heavy atoms "
                    f"(ratio {n_h_atoms/n_heavy_atoms:.2f}), which is insufficient. "
                    f"Do not submit this file to LigParGen without complete protonation."
                )
            else:  # "unknown"
                warnings.append(
                    "Hydrogenation was skipped (--hydrogenation none). "
                    "Ligand hydrogen completeness could not be verified from the PDB alone. "
                    "The file may not be LigParGen-ready."
                )

        elif needs_hydrogenation:
            backend = select_backend(hydrogenation_mode)
            h_pdb = out_dir / "ligand_for_ligpargen_H.pdb"
            h_result = backend.add_hydrogens(lig_pdb, h_pdb)
            hydrogenation_backend = h_result.backend_name

            if h_result.success:
                hydrogenated_pdb = h_result.output_path
                hydrogenation_performed = True
                hydrogenation_confidence = h_result.confidence
                hydrogenation_status = "complete"
                hydrogenation_complete = True
                if h_result.n_hydrogen_atoms:
                    n_h_atoms = h_result.n_hydrogen_atoms
            else:
                hydrogenation_confidence = 0.0
                if h_result.warning:
                    warnings.append(h_result.warning)
        # status == "unknown" and mode != "none": leave as-is, no action

    hydrogenation_required = (
        hydrogenation_status in ("missing", "incomplete") and not hydrogenation_performed
    )

    # Recommended LigParGen input
    if hydrogenated_pdb is not None:
        recommended_input: Optional[Path] = hydrogenated_pdb
    elif hydrogenation_complete:
        recommended_input = lig_pdb
    elif hydrogenation_status == "unknown":
        recommended_input = lig_pdb  # best effort; user should verify
    else:
        recommended_input = None  # not safe to submit

    # ── Estimate formal charge ────────────────────────────────────────────────
    charge_source = hydrogenated_pdb if hydrogenated_pdb is not None else lig_pdb
    formal_charge, charge_estimated = _estimate_formal_charge(charge_source, warnings)

    # ── Write report ──────────────────────────────────────────────────────────
    report_data = _build_report(
        ligand_resname=ligand_resname,
        ligand_atom_count=len(ligand_lines),
        protein_atom_count=protein_atom_count,
        formal_charge=formal_charge,
        charge_estimated=charge_estimated,
        has_hydrogens=has_h,
        n_hydrogen_atoms=n_h_atoms,
        n_heavy_atoms=n_heavy_atoms,
        hydrogenation_status=hydrogenation_status,
        hydrogenation_complete=hydrogenation_complete,
        hydrogenation_required=hydrogenation_required,
        hydrogenation_backend=hydrogenation_backend,
        hydrogenation_performed=hydrogenation_performed,
        hydrogenated_ligand_pdb=hydrogenated_pdb,
        hydrogenation_confidence=hydrogenation_confidence,
        recommended_ligpargen_input=recommended_input,
        lig_pdb=lig_pdb,
        prot_pdb=prot_pdb,
        warnings=warnings,
    )
    report_path = out_dir / "ligand_report.yaml"
    report_path.write_text(yaml.dump(report_data, default_flow_style=False, allow_unicode=True))

    return LigandPrepareResult(
        success=True,
        ligand_pdb=lig_pdb,
        protein_pdb=prot_pdb,
        report_path=report_path,
        ligand_resname=ligand_resname,
        ligand_atom_count=len(ligand_lines),
        protein_atom_count=protein_atom_count,
        formal_charge=formal_charge,
        charge_estimated=charge_estimated,
        has_hydrogens=has_h,
        n_hydrogen_atoms=n_h_atoms,
        n_heavy_atoms=n_heavy_atoms,
        hydrogenation_status=hydrogenation_status,
        hydrogenation_complete=hydrogenation_complete,
        hydrogenation_required=hydrogenation_required,
        hydrogenation_backend=hydrogenation_backend,
        hydrogenation_performed=hydrogenation_performed,
        hydrogenated_ligand_pdb=hydrogenated_pdb,
        hydrogenation_confidence=hydrogenation_confidence,
        recommended_ligpargen_input=recommended_input,
        warnings=warnings,
    )


# ── Hydrogen detection and counting ──────────────────────────────────────────

def _has_hydrogens(pdb_lines: list[str]) -> bool:
    """True if any ATOM/HETATM record is a hydrogen (element column or atom name)."""
    for line in pdb_lines:
        record = line[:6].rstrip()
        if record not in ("ATOM", "HETATM"):
            continue
        if len(line) >= 78:
            element = line[76:78].strip().upper()
            if element in ("H", "D"):
                return True
        if len(line) >= 14:
            atom_name = line[12:16].strip()
            if atom_name and atom_name[0].upper() in ("H", "D"):
                return True
    return False


def _count_h_atoms(pdb_lines: list[str]) -> int:
    """Count H/D atoms in ATOM/HETATM lines (element column, then atom name fallback)."""
    count = 0
    for line in pdb_lines:
        record = line[:6].rstrip()
        if record not in ("ATOM", "HETATM"):
            continue
        if len(line) >= 78:
            element = line[76:78].strip().upper()
            if element in ("H", "D"):
                count += 1
                continue
        if len(line) >= 14:
            atom_name = line[12:16].strip()
            if atom_name and atom_name[0].upper() in ("H", "D"):
                count += 1
    return count


def _count_heavy_atoms(pdb_lines: list[str]) -> int:
    """Count non-hydrogen atoms in ATOM/HETATM lines."""
    count = 0
    for line in pdb_lines:
        record = line[:6].rstrip()
        if record not in ("ATOM", "HETATM"):
            continue
        if len(line) >= 78:
            element = line[76:78].strip().upper()
            if element:
                if element not in ("H", "D"):
                    count += 1
                continue  # element column present — trust it
        if len(line) >= 14:
            atom_name = line[12:16].strip()
            if atom_name and atom_name[0].upper() not in ("H", "D"):
                count += 1
    return count


def _determine_hydrogenation_status(
    pdb_lines: list[str],
    n_h_atoms: int,
    n_heavy_atoms: int,
) -> tuple[str, bool]:
    """
    Return (hydrogenation_status, hydrogenation_complete).

    Status values:
      "missing"    -- zero H atoms; completeness=False
      "incomplete" -- H present but count is chemically insufficient; completeness=False
      "complete"   -- completeness proven via RDKit; completeness=True
      "unknown"    -- H present, ratio plausible, but can't prove without bond orders

    For PDB input without bond orders, completeness can only be proven by RDKit.
    The heuristic flag for "incomplete" uses H/heavy ratio < _INCOMPLETE_RATIO_THRESHOLD.
    This catches severe under-protonation (e.g. 28 heavy + 1 H = 0.036) while
    leaving ambiguous cases (e.g. 2 heavy + 1 H = 0.5) as "unknown".
    """
    if n_h_atoms == 0:
        return "missing", False

    # Try RDKit for exact completeness check
    try:
        from rdkit import Chem

        pdb_block = "\n".join(pdb_lines) + "\nEND\n"
        mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False, sanitize=True)
        if mol is not None:
            mol_no_h = Chem.RemoveHs(mol)
            mol_all_h = Chem.AddHs(mol_no_h)
            expected_h = mol_all_h.GetNumAtoms() - mol_no_h.GetNumAtoms()
            if n_h_atoms >= expected_h:
                return "complete", True
            else:
                return "incomplete", False
    except (ImportError, Exception):
        pass

    # Without RDKit: apply ratio heuristic
    if n_heavy_atoms > 0 and (n_h_atoms / n_heavy_atoms) < _INCOMPLETE_RATIO_THRESHOLD:
        return "incomplete", False

    # H present, ratio within plausible range, but can't prove — be conservative
    return "unknown", False


# ── Internal helpers ──────────────────────────────────────────────────────────

def _parse_serials(pdb_lines: list[str]) -> set[int]:
    serials: set[int] = set()
    for line in pdb_lines:
        try:
            serials.add(int(line[6:11]))
        except (ValueError, IndexError):
            pass
    return serials


def _conect_first_serial(line: str) -> Optional[int]:
    try:
        return int(line[6:11])
    except (ValueError, IndexError):
        return None


def _estimate_formal_charge(
    pdb_path: Path, warnings: list[str]
) -> tuple[Optional[int], bool]:
    try:
        from rdkit import Chem

        mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=False, sanitize=True)
        if mol is None:
            warnings.append(
                "RDKit could not parse the ligand PDB for charge estimation. "
                "Verify the charge manually before submitting to LigParGen."
            )
            return None, False
        return Chem.GetFormalCharge(mol), True
    except ImportError:
        warnings.append(
            "RDKit not available -- formal charge could not be estimated. "
            "Determine the charge manually before submitting to LigParGen."
        )
        return None, False
    except Exception as exc:
        warnings.append(f"Charge estimation failed ({exc}). Verify charge manually.")
        return None, False


def _build_report(
    ligand_resname: str,
    ligand_atom_count: int,
    protein_atom_count: int,
    formal_charge: Optional[int],
    charge_estimated: bool,
    has_hydrogens: bool,
    n_hydrogen_atoms: int,
    n_heavy_atoms: int,
    hydrogenation_status: str,
    hydrogenation_complete: bool,
    hydrogenation_required: bool,
    hydrogenation_backend: str,
    hydrogenation_performed: bool,
    hydrogenated_ligand_pdb: Optional[Path],
    hydrogenation_confidence: float,
    recommended_ligpargen_input: Optional[Path],
    lig_pdb: Path,
    prot_pdb: Path,
    warnings: list[str],
) -> dict:
    charge_str = str(formal_charge) if formal_charge is not None else "UNKNOWN"

    not_ready = (
        hydrogenation_status in ("missing", "incomplete")
        or hydrogenation_backend == "skipped"
    )

    if not_ready and not hydrogenation_performed:
        next_steps: list[str] = [
            "WARNING: The ligand is not LigParGen-ready. "
            f"Hydrogenation status: {hydrogenation_status}. "
            "Do NOT submit ligand_for_ligpargen.pdb to LigParGen without complete protonation.",
            "Option A -- install RDKit (conda install -c conda-forge rdkit) and re-run",
            "Option B -- install Open Babel (conda install -c conda-forge openbabel) and re-run",
            "Option C -- add hydrogens manually in Avogadro, PyMOL, or Discovery Studio",
        ]
    else:
        submit_file = recommended_ligpargen_input or lig_pdb
        next_steps = [
            f"Submit {submit_file.name} to LigParGen with charge={charge_str}",
            "Download the resulting .itp and .gro files from LigParGen",
            (
                f"Run: simforge ligand integrate "
                f"--protein {prot_pdb} "
                f"--ligand-itp <LIGAND>.itp "
                f"--ligand-gro <LIGAND>.gro "
                f"--out system/"
            ),
        ]

    return {
        "ligand_resname": ligand_resname,
        "ligand_atom_count": ligand_atom_count,
        "protein_atom_count": protein_atom_count,
        "formal_charge": formal_charge,
        "charge_estimated_by_rdkit": charge_estimated,
        "has_hydrogens": has_hydrogens,
        "n_hydrogen_atoms": n_hydrogen_atoms,
        "n_heavy_atoms": n_heavy_atoms,
        "hydrogenation_status": hydrogenation_status,
        "hydrogenation_complete": hydrogenation_complete,
        "hydrogenation_required": hydrogenation_required,
        "hydrogenation_backend": hydrogenation_backend,
        "hydrogenation_performed": hydrogenation_performed,
        "hydrogenated_ligand_pdb": str(hydrogenated_ligand_pdb) if hydrogenated_ligand_pdb else None,
        "hydrogenation_confidence": hydrogenation_confidence,
        "recommended_ligpargen_input": str(recommended_ligpargen_input) if recommended_ligpargen_input else None,
        "outputs": {
            "ligand_pdb": str(lig_pdb),
            "protein_pdb": str(prot_pdb),
        },
        "warnings": warnings,
        "next_steps": next_steps,
    }
