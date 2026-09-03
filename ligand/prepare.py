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
    # Chemical-perception fields (see ligand.chemical_perception). A successful
    # hydrogenation call does NOT by itself imply these are "high" -- see
    # parameterization_status / readiness_block_reasons, which are the
    # authoritative LigParGen-readiness signal.
    element_corrections: list[dict] = field(default_factory=list)
    element_assignment_confidence: str = "unknown"      # high|medium|low|unknown
    connectivity_confidence: str = "unknown"
    bond_order_confidence: str = "unknown"
    hydrogenation_chemical_confidence: str = "unknown"
    chemical_identity_method: str = "none"               # chemical_reference|pdb_geometry_valence_optimization|none
    chemical_identity_confidence: str = "unknown"
    parameterization_status: str = "chemical_validation_required"  # ready_for_ligpargen|chemical_validation_required
    readiness_block_reasons: list[str] = field(default_factory=list)


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
    chemical_reference: Optional[str | Path] = None,
) -> LigandPrepareResult:
    """
    Split a docked complex PDB into protein-only and ligand-only files.

    Ligand atom elements are normalized first (see
    ``ligand.chemical_perception.normalize_ligand_elements``), then either:

      - ``chemical_reference`` is given: bonds/orders/formal-charge/
        stereochemistry come from that authoritative source, mapped onto
        this pose's heavy-atom coordinates
        (``ligand.chemical_perception.map_reference_onto_pose``), or
      - PDB-only: hydrogen completeness/hydrogenation goes through
        geometry-based bond-order perception, never the legacy
        ``MolFromPDBFile()+AddHs()`` combination (see
        ``ligand.chemical_perception`` module docstring for why).

    ``LigandPrepareResult.parameterization_status`` is the authoritative
    LigParGen-readiness signal -- it is ``"ready_for_ligpargen"`` only when
    every chemical-confidence dimension is high, NOT merely because
    hydrogenation raised no exception. Use ``hydrogenation_mode="none"`` to
    skip automatic hydrogenation entirely; the call still reports readiness
    (which will not be ready without existing complete hydrogens).

    Args:
        complex_pdb:        Path to the docked complex PDB.
        ligand_resname:     Residue name of the ligand (e.g. "E20").
        out_dir:            Directory for output files.
        hydrogenation_mode: "auto" | "rdkit" | "openbabel" | "none"
        chemical_reference: Optional path to a chemically explicit reference
                             (.sdf/.mol/.pdb) treated as the authoritative
                             chemistry source for this ligand.
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

    # ── Element normalization (ligand-only, before any chemistry) ──────────────
    from ligand.chemical_perception import normalize_ligand_elements

    elem_norm = normalize_ligand_elements(lig_pdb.read_text())
    element_corrections = [
        {
            "atom_index": c.atom_index,
            "atom_name": c.atom_name,
            "original_element": c.original_element,
            "corrected_element": c.corrected_element,
            "reason": c.reason,
        }
        for c in elem_norm.corrections
    ]
    element_assignment_confidence = elem_norm.confidence
    if elem_norm.corrections:
        lig_pdb.write_text(elem_norm.pdb_text)
        warnings.append(
            f"Corrected {len(elem_norm.corrections)} ligand atom element assignment(s): "
            + "; ".join(
                f"atom {c.atom_index} ({c.atom_name}): {c.original_element} -> {c.corrected_element}"
                for c in elem_norm.corrections
            )
        )
    # Re-derive from the (possibly corrected) file for downstream counting.
    ligand_lines = [
        l for l in lig_pdb.read_text().splitlines()
        if l[:6].rstrip() in ("ATOM", "HETATM")
    ]

    # ── Hydrogen completeness assessment ──────────────────────────────────────
    n_h_atoms = _count_h_atoms(ligand_lines)
    n_heavy_atoms = _count_heavy_atoms(ligand_lines)
    has_h = n_h_atoms > 0

    hydrogenated_pdb: Optional[Path] = None
    hydrogenation_performed = False
    hydrogenation_backend = "none"
    hydrogenation_confidence = 0.0
    hydrogenation_chemical_confidence = "unknown"
    connectivity_confidence = "unknown"
    bond_order_confidence = "unknown"
    chemical_identity_method = "none"
    chemical_identity_confidence = "unknown"
    formal_charge: Optional[int] = None
    charge_estimated = False

    if chemical_reference is not None:
        # ── Reference-guided path: authoritative chemistry ─────────────────────
        from ligand.chemical_perception import map_reference_onto_pose, normalize_ligand_residue_metadata

        chemical_identity_method = "chemical_reference"
        mapping_result = map_reference_onto_pose(chemical_reference, lig_pdb)

        if not mapping_result.success:
            hydrogenation_status, hydrogenation_complete = "unknown", False
            connectivity_confidence = "low"
            bond_order_confidence = "low"
            chemical_identity_confidence = "low"
            warnings.append(f"Reference-guided preparation failed: {mapping_result.error}")
        else:
            from rdkit import Chem

            hydrogenated_pdb = out_dir / "ligand_for_ligpargen_H.pdb"
            pdb_block = Chem.MolToPDBBlock(mapping_result.mol)
            fixed_block = normalize_ligand_residue_metadata(pdb_block, resname=ligand_resname)
            hydrogenated_pdb.write_text(fixed_block)

            hydrogenation_performed = True
            hydrogenation_backend = "rdkit_reference_guided"
            hydrogenation_confidence = 1.0
            hydrogenation_chemical_confidence = "high"
            connectivity_confidence = "high"
            bond_order_confidence = "high"
            chemical_identity_confidence = "high"
            hydrogenation_status, hydrogenation_complete = "complete", True
            n_h_atoms = mapping_result.n_hydrogen_atoms
            formal_charge = mapping_result.formal_charge
            charge_estimated = True
            warnings.extend(mapping_result.warnings)

        needs_hydrogenation = False  # never fall through to backend cascade below

    else:
        # ── PDB-only path ────────────────────────────────────────────────────
        chemical_identity_method = "pdb_geometry_valence_optimization"

        # Probe geometry-based perception once up front so connectivity/
        # bond-order confidence is captured even when the ligand already has
        # enough hydrogens and no backend call ends up being needed below
        # (_determine_hydrogenation_status performs the same probe
        # internally to decide status/complete, but does not expose
        # confidence through its existing tuple return).
        if n_h_atoms > 0:
            from ligand.chemical_perception import perceive_chemistry

            _perception_probe = perceive_chemistry(lig_pdb)
            connectivity_confidence = _perception_probe.connectivity_confidence
            bond_order_confidence = _perception_probe.bond_order_confidence

        hydrogenation_status, hydrogenation_complete = _determine_hydrogenation_status(
            ligand_lines, n_h_atoms, n_heavy_atoms
        )
        needs_hydrogenation = hydrogenation_status in ("missing", "incomplete")
        hydrogenation_confidence = 1.0 if hydrogenation_complete else 0.0
        if hydrogenation_complete:
            if bond_order_confidence == "high":
                hydrogenation_chemical_confidence = "high"
            elif n_h_atoms > 0 and _perception_probe.method in ("unavailable", "legacy_pdb_flavor"):
                # RDKit isn't installed/too old to independently verify bond
                # orders -- but the input already had a complete hydrogen
                # set and SimForge computed/fabricated nothing, so this
                # isn't the failure mode behind the historical
                # mis-hydrogenation bug. Accept at "medium": sufficient for
                # decide_ligpargen_readiness's relaxed bar when
                # hydrogenation_performed=False, never for the "high" bar
                # required when SimForge itself adds hydrogens.
                connectivity_confidence = "medium"
                bond_order_confidence = "medium"
                hydrogenation_chemical_confidence = "medium"
            else:
                hydrogenation_chemical_confidence = "low"
            if formal_charge is None and n_h_atoms > 0 and _perception_probe.used_charge is not None:
                formal_charge = _perception_probe.used_charge
                charge_estimated = True

    if chemical_reference is None and not hydrogenation_complete:
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
            connectivity_confidence = h_result.connectivity_confidence
            bond_order_confidence = h_result.bond_order_confidence

            if h_result.success:
                hydrogenated_pdb = h_result.output_path
                hydrogenation_performed = True
                hydrogenation_confidence = h_result.confidence
                hydrogenation_chemical_confidence = h_result.hydrogenation_confidence
                hydrogenation_status = "complete"
                hydrogenation_complete = True
                if h_result.n_hydrogen_atoms:
                    n_h_atoms = h_result.n_hydrogen_atoms
                if h_result.formal_charge is not None:
                    formal_charge = h_result.formal_charge
                    charge_estimated = True
            else:
                hydrogenation_confidence = 0.0
                hydrogenation_chemical_confidence = "low"
                if h_result.warning:
                    warnings.append(h_result.warning)
        # status == "unknown" and mode != "none": leave as-is, no action

    if chemical_reference is None:
        # PDB-only chemical identity is never more than "medium" confidence
        # (no external reference), even when geometry-based perception fully
        # converged -- see decide_ligpargen_readiness for the hard gate.
        chemical_identity_confidence = "medium" if bond_order_confidence == "high" else "low"

    # ── Normalize residue metadata on whatever backend produced hydrogens ──────
    # RDKit's AddHs()/some Open Babel output writes new atoms into a generic
    # "UNL" residue with a blank chain; force every atom back to the ligand's
    # own resname/chain/resid.
    if hydrogenated_pdb is not None:
        from ligand.chemical_perception import normalize_ligand_residue_metadata

        hydrogenated_pdb.write_text(
            normalize_ligand_residue_metadata(hydrogenated_pdb.read_text(), resname=ligand_resname)
        )

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

    # ── Estimate formal charge (fallback only -- prefer the value already
    # derived from geometry-based/reference-guided chemical perception) ────────
    if formal_charge is None:
        charge_source = hydrogenated_pdb if hydrogenated_pdb is not None else lig_pdb
        formal_charge, charge_estimated = _estimate_formal_charge(charge_source, warnings)

    # ── LigParGen readiness gate ────────────────────────────────────────────────
    from ligand.chemical_perception import decide_ligpargen_readiness

    parameterization_status, readiness_block_reasons = decide_ligpargen_readiness(
        element_assignment_confidence=element_assignment_confidence,
        connectivity_confidence=connectivity_confidence,
        bond_order_confidence=bond_order_confidence,
        hydrogenation_status=hydrogenation_status,
        hydrogenation_confidence=hydrogenation_chemical_confidence,
        hydrogenation_performed=hydrogenation_performed,
    )

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
        element_corrections=element_corrections,
        element_assignment_confidence=element_assignment_confidence,
        connectivity_confidence=connectivity_confidence,
        bond_order_confidence=bond_order_confidence,
        hydrogenation_chemical_confidence=hydrogenation_chemical_confidence,
        chemical_identity_method=chemical_identity_method,
        chemical_identity_confidence=chemical_identity_confidence,
        parameterization_status=parameterization_status,
        readiness_block_reasons=readiness_block_reasons,
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
        element_corrections=element_corrections,
        element_assignment_confidence=element_assignment_confidence,
        connectivity_confidence=connectivity_confidence,
        bond_order_confidence=bond_order_confidence,
        hydrogenation_chemical_confidence=hydrogenation_chemical_confidence,
        chemical_identity_method=chemical_identity_method,
        chemical_identity_confidence=chemical_identity_confidence,
        parameterization_status=parameterization_status,
        readiness_block_reasons=readiness_block_reasons,
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

    # Try geometry-based bond-order perception for an exact completeness
    # check (ligand.chemical_perception.perceive_chemistry). This replaced
    # the legacy Chem.MolFromPDBBlock()+RemoveHs()+AddHs() round-trip, which
    # does not perceive aromaticity/bond order from geometry and therefore
    # mis-classifies (and, if used for hydrogenation itself, over-protonates)
    # ring-containing and element-column-defective ligands -- see
    # ligand.chemical_perception module docstring for the real regression
    # this caused. When perception does not converge (e.g. RDKit absent, or
    # too old for rdDetermineBonds), fall through to the ratio heuristic
    # exactly as before.
    try:
        import tempfile

        from ligand.chemical_perception import perceive_chemistry

        with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as f:
            f.write("\n".join(pdb_lines) + "\nEND\n")
            tmp_path = f.name
        try:
            perception = perceive_chemistry(tmp_path)
        finally:
            Path(tmp_path).unlink(missing_ok=True)

        if perception.success:
            from rdkit import Chem

            mol_no_h = Chem.RemoveHs(perception.mol)
            mol_all_h = Chem.AddHs(mol_no_h)
            expected_h = mol_all_h.GetNumAtoms() - mol_no_h.GetNumAtoms()
            if n_h_atoms >= expected_h:
                return "complete", True
            else:
                return "incomplete", False
    except Exception:
        pass

    # Without confident geometry-based perception: apply ratio heuristic
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
    element_corrections: list[dict],
    element_assignment_confidence: str,
    connectivity_confidence: str,
    bond_order_confidence: str,
    hydrogenation_chemical_confidence: str,
    chemical_identity_method: str,
    chemical_identity_confidence: str,
    parameterization_status: str,
    readiness_block_reasons: list[str],
) -> dict:
    charge_str = str(formal_charge) if formal_charge is not None else "UNKNOWN"
    ready = parameterization_status == "ready_for_ligpargen"

    if not ready:
        next_steps: list[str] = [
            "WARNING: chemical_validation_required -- this ligand is NOT "
            "LigParGen-ready. A successful hydrogenation call is not, by "
            "itself, sufficient (see readiness_block_reasons below). Do NOT "
            "submit ligand_for_ligpargen.pdb to LigParGen.",
        ] + [f"  - {reason}" for reason in readiness_block_reasons] + [
            "Option A -- provide --ligand-reference <SDF/MOL/PDB> with the "
            "correct connectivity/bond orders/charge for this ligand",
            "Option B -- supply a pose PDB that already retains explicit, "
            "chemically correct hydrogens",
            "Option C -- resolve chemistry manually (Avogadro, PyMOL, "
            "Discovery Studio) and re-run with --hydrogenation none",
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
        # ── Chemical perception provenance ──────────────────────────────────
        "element_corrections": element_corrections,
        "element_assignment_confidence": element_assignment_confidence,
        "connectivity_confidence": connectivity_confidence,
        "bond_order_confidence": bond_order_confidence,
        "hydrogenation_chemical_confidence": hydrogenation_chemical_confidence,
        "chemical_identity_method": chemical_identity_method,
        "chemical_identity_confidence": chemical_identity_confidence,
        "parameterization_status": parameterization_status,
        "readiness_block_reasons": readiness_block_reasons,
        "outputs": {
            "ligand_pdb": str(lig_pdb),
            "protein_pdb": str(prot_pdb),
        },
        "warnings": warnings,
        "next_steps": next_steps,
    }
