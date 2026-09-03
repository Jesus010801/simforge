"""
Ligand chemical-perception layer: element normalization, geometry-based
connectivity/bond-order determination, reference-guided chemistry mapping,
and the LigParGen-readiness gate.

Why this module exists
-----------------------
SimForge's historical PDB-only hydrogenation path called
``Chem.MolFromPDBFile(..., sanitize=True)`` followed by ``Chem.AddHs()``.
That combination has two independent failure modes that compound:

  1. It trusts the PDB element column verbatim. A common docking-tool
     defect writes atom name "CL" with element column "C" (chlorine
     truncated to its first character). RDKit then treats the atom as an
     under-bonded carbon and fills it with extra implicit hydrogens.

  2. ``MolFromPDBFile``'s legacy connectivity perception does not assign
     bond orders or aromaticity from geometry alone -- every aromatic ring
     atom looks "unsaturated" by default sp3-carbon valence rules, so
     ``AddHs`` over-protonates the whole molecule, not just the mis-typed
     atom.

For a real regression ligand (A6, C18H15ClO3, 22 heavy atoms) this produced
an approximately C19H34O3 molecule that was reported as
"Hydrogenation: complete / Parameterization: ready_for_ligpargen" -- wrong,
and confidently so.

This module fixes both failure modes and, critically, separates
"hydrogenation completed without raising an exception" from "the resulting
chemistry is trustworthy enough to submit to LigParGen":

  - ``normalize_ligand_elements``   -- context-aware, ligand-only element
    correction for the small set of genuinely unambiguous PDB writer
    defects (CL/BR/F/I truncated to their first character). Never touches
    unrelated names (e.g. protein "CA" stays carbon-alpha).

  - ``perceive_chemistry``          -- geometry-based bond-order
    determination via ``rdkit.Chem.rdDetermineBonds.DetermineBonds``
    (RDKit's xyz2mol-derived valence-optimization algorithm), searching a
    small set of plausible net charges. This replaces the legacy
    ``MolFromPDBFile`` connectivity path entirely for automatic
    hydrogenation -- that legacy path is not used anywhere in this module.
    Empirically (see ligand/test_chemical_perception.py), this converges
    reliably when explicit hydrogens are present in the input pose, and
    reliably *fails to converge* -- safely, loudly -- for heavy-atom-only
    poses, which is the honest answer: a heavy skeleton's geometry alone
    does not determine where hydrogens belong.

  - ``map_reference_onto_pose``     -- when a chemically explicit reference
    (SDF/MOL/PDB) is available, its bonds/orders/charge/stereochemistry are
    authoritative. Heavy atoms are mapped reference<->docked-pose (by atom
    name when both sides share a naming scheme, falling back to geometric
    nearest-neighbour matching), the reference's *coordinates* are replaced
    by the docked pose's, and hydrogens are added against the trusted bond
    graph.

  - ``decide_ligpargen_readiness``  -- the single gate: every confidence
    dimension (element assignment, connectivity, bond order, hydrogenation)
    must independently be "high" for ``ready_for_ligpargen``. Anything less
    is ``chemical_validation_required`` -- explicit and actionable, never
    silent.

  - ``normalize_ligand_residue_metadata`` -- RDKit's ``AddHs`` writes newly
    created hydrogens into a generic "UNL" residue with a blank chain.
    This rewrites every atom's resname/chain/resid to match the ligand's
    original identity.

RDKit is imported lazily throughout (consistent with ``ligand.rdkit_reader``)
-- this module is safely importable without RDKit installed; callers get a
clear failure result, never an ImportError from module load.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

# ── Confidence vocabulary ────────────────────────────────────────────────────

ConfidenceLevel = str  # "high" | "medium" | "low" | "unknown"

_CONFIDENCE_RANK: dict[str, int] = {"high": 3, "medium": 2, "low": 1, "unknown": 0}

# Net charges tried, in order, when determining bond orders from geometry.
# Neutral first (the common case), then the most common small ionic charges.
_CHARGE_SEARCH_ORDER: tuple[int, ...] = (0, -1, 1, -2, 2)


# ═══════════════════════════════════════════════════════════════════════════
# 1. Element normalization
# ═══════════════════════════════════════════════════════════════════════════

# Ligand atom NAMES that unambiguously identify a halogen even when the PDB
# element column is wrong or blank. Deliberately a short, explicit allowlist
# -- NOT a general "guess the element from the name" rule. In particular,
# "CA" is intentionally absent: within a ligand HETATM block it is
# overwhelmingly a generic ligand carbon name, and calcium *ions* are
# filtered out upstream (ligand.prepare._STANDARD_RESIDUES) before this
# function ever sees them, so no context-dependent guessing is needed here.
_LIGAND_HALOGEN_ELEMENT_BY_NAME: dict[str, str] = {
    "CL": "Cl",
    "BR": "Br",
    "F": "F",
    "I": "I",
}


@dataclass(frozen=True)
class ElementCorrection:
    atom_index: int
    atom_name: str
    original_element: str
    corrected_element: str
    reason: str = "ligand_atom_name_element_conflict"


@dataclass
class ElementNormalizationResult:
    pdb_text: str
    corrections: list[ElementCorrection] = field(default_factory=list)
    confidence: ConfidenceLevel = "high"
    warnings: list[str] = field(default_factory=list)


def _proper_case_element(symbol: str) -> str:
    """RDKit's canonical symbol casing (e.g. "Cl", not "CL"/"cl")."""
    symbol = symbol.strip()
    if len(symbol) > 1:
        return symbol[0].upper() + symbol[1:].lower()
    return symbol.upper()


def _read_original_element(line: str, atom_name: str) -> tuple[str, bool]:
    """
    Best-effort element extraction, tolerant of ligand PDB output that does
    not follow strict 78-column PDB widths (common for docking-tool output
    that stops after occupancy/temperature-factor with no element/segID
    columns at all).

    Returns (element, confidently_read). Strategy, in order:
      1. Strict columns 77-78, when the line is long enough to have them.
      2. The last whitespace-separated token, when it plausibly looks like
         a bare element symbol (1-2 alphabetic characters) -- catches
         shorter, loosely-formatted lines that still end with the element.
      3. The atom name's first character, as an absolute last resort (not
         "confidently read" -- callers should not treat this as proof the
         element was actually present in the source file).
    """
    if len(line) >= 78:
        candidate = line[76:78].strip()
        if candidate:
            return candidate.upper(), True
    tokens = line.split()
    if tokens:
        last = tokens[-1]
        if 1 <= len(last) <= 2 and last.isalpha():
            return last.upper(), True
    return (atom_name[0].upper() if atom_name else "X"), False


def _rebuild_pdb_atom_line(line: str, element: str) -> str:
    """
    Rewrite an ATOM/HETATM line with ``element`` correctly serialized into
    PDB columns 77-78, guaranteed regardless of the input line's original
    width.

    Record type, serial, atom name, altloc, residue name, chain, residue
    number, insertion code, and x/y/z coordinates (columns 1-54) are
    preserved byte-for-byte -- this function never touches geometry.
    Occupancy/temperature-factor (columns 55-66) are preserved when the
    input line has them, defaulted to "1.00"/"0.00" otherwise.
    """
    if len(line) < 54:
        return line  # not enough of a record to safely reconstruct

    head = line[:54]
    occ_temp = line[54:66] if len(line) >= 66 else f"{1.00:6.2f}{0.00:6.2f}"
    return head + occ_temp + " " * 10 + element.rjust(2)


def normalize_ligand_elements(pdb_text: str) -> ElementNormalizationResult:
    """
    Correct unambiguous atom-name/element-column contradictions in a
    ligand-only PDB text (HETATM/ATOM records for a single non-standard
    residue), and unconditionally re-serialize every atom's element into
    proper PDB columns 77-78 -- regardless of whether the *input* line
    followed strict PDB column widths.

    This is not merely a correction pass: it is the sole point where a
    normalized element is written back into the PDB representation that
    everything downstream (bond-order perception, hydrogenation backends,
    the final ligand_for_ligpargen.pdb) reads from. A hydrogenation
    backend's own writer (Open Babel, RDKit) generally serializes elements
    correctly *given a well-formed input* -- the historical gap was that a
    loosely-formatted input (common for real docking-tool PDB output, which
    often stops after occupancy/temperature-factor with no element/segID
    columns at all) was previously skipped outright (``len(line) < 78``),
    so neither the correction nor a valid element field was ever produced
    for such atoms.

    Only fires corrections for names in ``_LIGAND_HALOGEN_ELEMENT_BY_NAME``
    -- every correction is therefore unambiguous by construction, which is
    why ``confidence`` stays "high" whenever this function returns
    normally. Columns 1-54 (record/serial/name/resname/chain/resid/
    coordinates) are always preserved byte-for-byte; no other field is
    ever altered.
    """
    corrections: list[ElementCorrection] = []
    out_lines: list[str] = []

    for line in pdb_text.splitlines():
        record = line[:6].rstrip()
        if record not in ("ATOM", "HETATM") or len(line) < 54:
            out_lines.append(line)
            continue

        atom_name = line[12:16].strip().upper()
        try:
            serial = int(line[6:11])
        except ValueError:
            serial = -1

        original_element, confidently_read = _read_original_element(line, atom_name)
        expected = _LIGAND_HALOGEN_ELEMENT_BY_NAME.get(atom_name)

        if expected is not None and original_element != expected.upper():
            corrections.append(ElementCorrection(
                atom_index=serial,
                atom_name=atom_name,
                original_element=original_element if confidently_read else "(blank)",
                corrected_element=expected,
                reason="ligand_atom_name_element_conflict",
            ))
            final_element = expected
        elif expected is not None:
            final_element = expected  # already correct; keep canonical casing
        else:
            final_element = original_element

        # Always written in RDKit's canonical proper-case symbol (e.g. "Cl",
        # not "CL"/"cl") -- matches the convention used everywhere else in
        # this codebase (ligand.rdkit_reader.heavy_atom_elements,
        # ligand.pose_rewriter._element_from_mass), so downstream
        # element-list comparisons (e.g. campaign ligand-identity grouping
        # against a --ligand-reference) don't silently fail on a case
        # mismatch, and the element field is guaranteed present and
        # correctly positioned regardless of the input line's own width.
        out_lines.append(_rebuild_pdb_atom_line(line, _proper_case_element(final_element)))

    trailing_newline = "\n" if pdb_text.endswith("\n") else ""
    return ElementNormalizationResult(
        pdb_text="\n".join(out_lines) + trailing_newline,
        corrections=corrections,
        confidence="high",
    )


# ═══════════════════════════════════════════════════════════════════════════
# 2. Geometry-based connectivity / bond-order perception
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class ChemicalPerceptionResult:
    success: bool
    mol: Optional[object] = None  # rdkit.Chem.Mol when success is True
    used_charge: Optional[int] = None
    method: str = "none"  # "geometry_valence_optimization" | "legacy_pdb_flavor" | "none"
    connectivity_confidence: ConfidenceLevel = "low"
    bond_order_confidence: ConfidenceLevel = "low"
    warning: Optional[str] = None


def perceive_chemistry(
    pdb_path: str | Path,
    charge_candidates: tuple[int, ...] = _CHARGE_SEARCH_ORDER,
) -> ChemicalPerceptionResult:
    """
    Determine connectivity and bond orders for a ligand PDB purely from 3D
    geometry, using RDKit's valence-optimization algorithm
    (``rdDetermineBonds.DetermineBonds``, RDKit >= 2022.09).

    Always applies ``normalize_ligand_elements`` first, internally, so this
    function is safe to call directly regardless of what the caller has
    already done. This is not merely defense in depth: an un-normalized
    element (e.g. a "CL"-named atom whose PDB element column says "C") can
    make ``DetermineBonds`` *converge* to a self-consistent but wrong
    structure at a shifted net charge (empirically: the real A6 regression
    ligand converges to a plausible-looking C19H15O3(-1) if fed uncorrected)
    -- silent convergence to the wrong answer is exactly the failure mode
    this module exists to prevent, so element correctness cannot be left to
    caller discipline.

    Tries each of ``charge_candidates`` in order and accepts the first net
    charge for which a self-consistent, sanitizable bond-order assignment
    is found. Returns success=False (never a guessed/low-confidence mol)
    when none converge -- this is the expected, correct outcome for a
    heavy-atom-only pose with no explicit hydrogens, where valence
    completion is fundamentally underdetermined by geometry alone.
    """
    try:
        from rdkit import Chem
    except ImportError:
        return ChemicalPerceptionResult(success=False, method="unavailable", warning="RDKit is not installed.")

    try:
        from rdkit.Chem import rdDetermineBonds
    except ImportError:
        return ChemicalPerceptionResult(
            success=False,
            method="legacy_pdb_flavor",
            warning=(
                "RDKit is installed but lacks rdkit.Chem.rdDetermineBonds "
                "(requires RDKit >= 2022.09) for geometry-based bond-order "
                "determination. Upgrade RDKit, or provide --ligand-reference."
            ),
        )

    if not Path(pdb_path).exists():
        return ChemicalPerceptionResult(success=False, warning=f"Ligand PDB not found: {pdb_path}")

    from validators.ligand_parsers.pdb_parser import parse_pdb_ligand

    norm = normalize_ligand_elements(Path(pdb_path).read_text())
    normalized_tmp: Optional[Path] = None
    if norm.corrections:
        import tempfile

        with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as f:
            f.write(norm.pdb_text)
            normalized_tmp = Path(f.name)
        pdb_path = normalized_tmp

    parsed = parse_pdb_ligand(Path(pdb_path))
    if normalized_tmp is not None:
        normalized_tmp.unlink(missing_ok=True)
    if parsed.get("error"):
        return ChemicalPerceptionResult(success=False, warning=f"Cannot parse ligand PDB: {parsed['error']}")

    atoms = parsed["atoms"]
    if not atoms:
        return ChemicalPerceptionResult(success=False, warning="No atoms found in ligand PDB.")

    from rdkit.Geometry import Point3D

    def _build_mol():
        rw = Chem.RWMol()
        for a in atoms:
            rw.AddAtom(Chem.Atom(a["element"]))
        conf = Chem.Conformer(rw.GetNumAtoms())
        for i, a in enumerate(atoms):
            conf.SetAtomPosition(i, Point3D(float(a["x"]), float(a["y"]), float(a["z"])))
        rw.AddConformer(conf)
        return rw.GetMol()

    attempts: list[str] = []
    for charge in charge_candidates:
        mol = _build_mol()
        try:
            rdDetermineBonds.DetermineBonds(mol, charge=charge)
            Chem.SanitizeMol(mol)
        except Exception as exc:
            attempts.append(f"charge={charge}: {exc}")
            continue

        return ChemicalPerceptionResult(
            success=True,
            mol=mol,
            used_charge=charge,
            method="geometry_valence_optimization",
            connectivity_confidence="high",
            bond_order_confidence="high",
        )

    return ChemicalPerceptionResult(
        success=False,
        method="geometry_valence_optimization_failed",
        connectivity_confidence="low",
        bond_order_confidence="low",
        warning=(
            "Geometry-based bond-order determination did not converge for any "
            f"plausible net charge ({', '.join(str(c) for c in charge_candidates)}). "
            "This is expected for a heavy-atom-only pose with no explicit "
            "hydrogens -- valence completion cannot be resolved from geometry "
            "alone. Provide --ligand-reference (SDF/MOL) or a pose retaining "
            "explicit hydrogens. Attempts: " + "; ".join(attempts)
        ),
    )


@dataclass
class ChemicalPreparationResult:
    success: bool
    output_path: Optional[Path] = None
    n_hydrogen_atoms: int = 0
    formal_charge: Optional[int] = None
    method: str = "none"
    connectivity_confidence: ConfidenceLevel = "low"
    bond_order_confidence: ConfidenceLevel = "low"
    hydrogenation_confidence: ConfidenceLevel = "low"
    warning: Optional[str] = None


def add_hydrogens_with_confidence(
    input_pdb: str | Path,
    output_pdb: str | Path,
) -> ChemicalPreparationResult:
    """
    Add hydrogens ONLY when geometry-based bond-order perception succeeds
    with high confidence. Never falls back to naive
    ``MolFromPDBFile()+AddHs()`` -- that path is the exact mechanism behind
    SimForge's historical mis-hydrogenation bug (see module docstring).
    """
    perception = perceive_chemistry(input_pdb)
    if not perception.success:
        return ChemicalPreparationResult(
            success=False,
            method=perception.method,
            connectivity_confidence=perception.connectivity_confidence,
            bond_order_confidence=perception.bond_order_confidence,
            warning=perception.warning,
        )

    from rdkit import Chem

    mol = perception.mol
    mol_h = Chem.AddHs(mol, addCoords=True)
    pdb_block = Chem.MolToPDBBlock(mol_h)
    if not pdb_block:
        return ChemicalPreparationResult(
            success=False,
            method=perception.method,
            connectivity_confidence=perception.connectivity_confidence,
            bond_order_confidence=perception.bond_order_confidence,
            warning="RDKit could not write the hydrogenated PDB.",
        )

    output_pdb = Path(output_pdb)
    output_pdb.write_text(pdb_block)

    n_h = sum(1 for a in mol_h.GetAtoms() if a.GetSymbol() == "H")
    return ChemicalPreparationResult(
        success=True,
        output_path=output_pdb,
        n_hydrogen_atoms=n_h,
        formal_charge=Chem.GetFormalCharge(mol),
        method=perception.method,
        connectivity_confidence=perception.connectivity_confidence,
        bond_order_confidence=perception.bond_order_confidence,
        hydrogenation_confidence="high",
    )


# ═══════════════════════════════════════════════════════════════════════════
# 3. Reference-guided chemistry mapping
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class ReferenceMappingResult:
    success: bool
    mol: Optional[object] = None  # hydrogenated rdkit.Chem.Mol with receptor coordinates
    formal_charge: Optional[int] = None
    n_hydrogen_atoms: int = 0
    mapping_method: str = "none"  # "atom_name" | "geometric_nearest_neighbor"
    warnings: list = field(default_factory=list)
    error: Optional[str] = None


def _read_heavy_atom_names(path: Path) -> list[str]:
    names: list[str] = []
    for line in Path(path).read_text().splitlines():
        record = line[:6].rstrip()
        if record not in ("ATOM", "HETATM"):
            continue
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if element in ("H", "D"):
            continue
        names.append(line[12:16].strip().upper())
    return names


def map_reference_onto_pose(
    reference_path: str | Path,
    docked_heavy_pdb: str | Path,
) -> ReferenceMappingResult:
    """
    Build a chemically correct molecule by combining a trusted chemical
    reference's bonds/orders/formal-charge/stereochemistry with a docked
    pose's heavy-atom coordinates.

    Heavy-atom correspondence is established by atom name when both the
    reference and the docked pose share a consistent PDB naming scheme
    (cross-validated against element compatibility so a name match can
    never silently override a genuine element conflict); otherwise by
    centroid-aligned geometric nearest-neighbour matching per element
    (reusing ``ligand.pose_rewriter._match_atoms``). Fails explicitly on a
    heavy-atom count mismatch or an unresolvable mapping -- never guesses.

    Only heavy-atom coordinates are taken from the docked pose; hydrogens
    are (re)placed by RDKit against the reference's trusted bond graph.
    """
    try:
        from rdkit import Chem
    except ImportError:
        return ReferenceMappingResult(success=False, error="RDKit is required for reference-guided preparation.")

    import numpy as np
    from rdkit.Geometry import Point3D

    from ligand.pose_rewriter import _match_atoms
    from ligand.rdkit_reader import load_mol
    from validators.ligand_parsers.pdb_parser import parse_pdb_ligand

    reference_path = Path(reference_path)
    docked_heavy_pdb = Path(docked_heavy_pdb)
    ref_names: list[str] = []

    if reference_path.suffix.lower() == ".pdb":
        # A PDB reference carries no explicit bond orders, so -- exactly like
        # the docked pose -- it must go through geometry-based perception
        # rather than the legacy MolFromPDBFile()+sanitize path (which
        # under-perceives aromatic/ring bonds; see module docstring). This
        # also lets us keep the PDB atom names (needed for Step 1 below),
        # which a from-scratch RWMol would not otherwise carry.
        perception = perceive_chemistry(reference_path)
        if not perception.success:
            return ReferenceMappingResult(
                success=False,
                error=f"Chemical reference could not be resolved from geometry: {perception.warning}",
            )
        ref_mol = perception.mol
        ref_names_all = _read_heavy_atom_names(reference_path)  # heavy-only, matches RemoveHs order below
    else:
        try:
            ref_mol = load_mol(reference_path)
        except (FileNotFoundError, ValueError) as exc:
            return ReferenceMappingResult(success=False, error=f"Cannot load chemical reference: {exc}")
        has_conformer = ref_mol.GetNumConformers() > 0
        ref_mol = Chem.AddHs(ref_mol, addCoords=has_conformer)
        try:
            Chem.SanitizeMol(ref_mol)
        except Exception as exc:
            return ReferenceMappingResult(success=False, error=f"Chemical reference failed sanitization: {exc}")
        ref_names_all = []

    ref_noh = Chem.RemoveHs(ref_mol)
    if ref_noh.GetNumConformers() == 0:
        return ReferenceMappingResult(
            success=False,
            error="Chemical reference has no 3D conformer; cannot establish geometric atom correspondence.",
        )
    ref_conf = ref_noh.GetConformer()
    ref_elements = [a.GetSymbol() for a in ref_noh.GetAtoms()]
    ref_names = ref_names_all if len(ref_names_all) == ref_noh.GetNumAtoms() else []
    ref_coords = np.array(
        [[ref_conf.GetAtomPosition(i).x, ref_conf.GetAtomPosition(i).y, ref_conf.GetAtomPosition(i).z]
         for i in range(ref_noh.GetNumAtoms())]
    )

    docked = parse_pdb_ligand(docked_heavy_pdb)
    if docked.get("error"):
        return ReferenceMappingResult(success=False, error=f"Cannot parse docked pose: {docked['error']}")
    docked_heavy = [a for a in docked["atoms"] if a["element"].upper() != "H"]
    docked_elements = [a["element"] for a in docked_heavy]
    docked_coords = np.array([[a["x"], a["y"], a["z"]] for a in docked_heavy], dtype=float)
    docked_names = [n for n in _read_heavy_atom_names(docked_heavy_pdb)]

    if len(ref_elements) != len(docked_elements):
        return ReferenceMappingResult(
            success=False,
            error=(
                f"Heavy-atom count mismatch: reference has {len(ref_elements)} heavy "
                f"atoms, docked pose has {len(docked_elements)}."
            ),
        )

    n_atoms = len(ref_elements)
    mapping: list[Optional[int]] = [None] * n_atoms
    warnings: list[str] = []
    name_seeded = 0

    # Step 1 -- pre-seed correspondence from atom names that are unique on
    # BOTH sides (e.g. a single "CL" among many generically-named "C"
    # atoms). Trusted regardless of the docked side's element column --
    # the whole point is to stay correct even when that column is wrong
    # (a "CL"-named atom mistyped as element "C"). Disagreements are
    # recorded as warnings, not treated as fatal: the reference's element
    # is authoritative once mapped. Ambiguous names (appearing more than
    # once, e.g. generic "C") are deliberately left for Step 2.
    if docked_names and len(docked_names) == len(docked_elements):
        from collections import Counter
        ref_name_counts = Counter(ref_names)
        docked_name_counts = Counter(docked_names)
        docked_idx_by_name: dict[str, int] = {}
        for i, n in enumerate(docked_names):
            if docked_name_counts[n] == 1:
                docked_idx_by_name[n] = i

        for ref_idx, name in enumerate(ref_names):
            if not name or ref_name_counts[name] != 1:
                continue
            docked_idx = docked_idx_by_name.get(name)
            if docked_idx is None:
                continue
            mapping[ref_idx] = docked_idx
            name_seeded += 1
            if docked_elements[docked_idx].upper() != ref_elements[ref_idx].upper():
                warnings.append(
                    f"Atom name '{name}' mapped by name despite a docked-pose "
                    f"element disagreement (reference={ref_elements[ref_idx]}, "
                    f"docked={docked_elements[docked_idx]}); trusting the "
                    "reference element."
                )

    # Step 2 -- element + centroid-aligned nearest-neighbour matching
    # (ligand.pose_rewriter._match_atoms) for whatever remains unmapped.
    remaining_ref_idx = [i for i in range(n_atoms) if mapping[i] is None]
    used_docked_idx = {j for j in mapping if j is not None}
    remaining_docked_idx = [j for j in range(n_atoms) if j not in used_docked_idx]

    if remaining_ref_idx:
        sub_ref_elements = [ref_elements[i] for i in remaining_ref_idx]
        sub_ref_coords = ref_coords[remaining_ref_idx]
        sub_docked_elements = [docked_elements[j] for j in remaining_docked_idx]
        sub_docked_coords = docked_coords[remaining_docked_idx]
        try:
            sub_mapping = _match_atoms(sub_ref_elements, sub_ref_coords, sub_docked_elements, sub_docked_coords)
        except ValueError as exc:
            return ReferenceMappingResult(success=False, error=f"Heavy-atom mapping failed: {exc}")
        for local_ref, local_docked in zip(remaining_ref_idx, sub_mapping):
            mapping[local_ref] = remaining_docked_idx[local_docked]

    mapping_method = (
        "atom_name" if name_seeded == n_atoms
        else "atom_name+geometric_nearest_neighbor" if name_seeded > 0
        else "geometric_nearest_neighbor"
    )

    rw = Chem.RWMol(ref_noh)
    conf = rw.GetConformer()
    for ref_idx, docked_idx in enumerate(mapping):
        x, y, z = docked_coords[docked_idx]
        conf.SetAtomPosition(ref_idx, Point3D(float(x), float(y), float(z)))

    mol_repositioned = rw.GetMol()
    try:
        Chem.SanitizeMol(mol_repositioned)
    except Exception as exc:
        return ReferenceMappingResult(success=False, error=f"Repositioned molecule failed sanitization: {exc}")

    mol_h = Chem.AddHs(mol_repositioned, addCoords=True)

    return ReferenceMappingResult(
        success=True,
        mol=mol_h,
        formal_charge=Chem.GetFormalCharge(mol_repositioned),
        n_hydrogen_atoms=sum(1 for a in mol_h.GetAtoms() if a.GetSymbol() == "H"),
        mapping_method=mapping_method,
        warnings=warnings,
    )


# ═══════════════════════════════════════════════════════════════════════════
# 4. Residue metadata normalization
# ═══════════════════════════════════════════════════════════════════════════

def normalize_ligand_residue_metadata(
    pdb_text: str,
    resname: str,
    chain: str = "A",
    resid: int = 1,
) -> str:
    """
    Rewrite resname/chain/resid to a single consistent identity across every
    ATOM/HETATM line.

    RDKit's ``AddHs`` (and some Open Babel output modes) write newly created
    hydrogens into a generic "UNL" residue with a blank chain, leaving a PDB
    where heavy atoms say e.g. "LIG A 1" and their own hydrogens say
    "UNL   1" -- invalid for any downstream tool that keys off residue
    identity. This does not touch coordinates, element, or atom name fields.
    """
    resname_field = resname.strip().upper()[:3].rjust(3)
    chain_field = (chain.strip() or "A")[0]
    resid_field = f"{resid:>4}"

    out_lines: list[str] = []
    for line in pdb_text.splitlines():
        record = line[:6].rstrip()
        if record in ("ATOM", "HETATM") and len(line) >= 26:
            line = line[:17] + resname_field + " " + chain_field + resid_field + line[26:]
        out_lines.append(line)

    trailing_newline = "\n" if pdb_text.endswith("\n") else ""
    return "\n".join(out_lines) + trailing_newline


# ═══════════════════════════════════════════════════════════════════════════
# 5. LigParGen readiness gate
# ═══════════════════════════════════════════════════════════════════════════

def decide_ligpargen_readiness(
    element_assignment_confidence: ConfidenceLevel,
    connectivity_confidence: ConfidenceLevel,
    bond_order_confidence: ConfidenceLevel,
    hydrogenation_status: str,
    hydrogenation_confidence: ConfidenceLevel,
    hydrogenation_performed: bool = True,
) -> tuple[str, list[str]]:
    """
    The single LigParGen-readiness decision point.

    Returns ``(parameterization_status, block_reasons)`` where status is
    ``"ready_for_ligpargen"`` only when hydrogenation is complete AND every
    confidence dimension clears its required bar; otherwise
    ``"chemical_validation_required"`` with an explicit, actionable list of
    which dimensions fell short.

    A successful ``AddHs()`` call is not, by itself, sufficient --
    ``hydrogenation_confidence`` reflects whether the bond graph AddHs
    operated on was itself trustworthy, not merely whether the RDKit call
    raised no exception.

    ``hydrogenation_performed`` controls how strict the connectivity/
    bond-order/hydrogenation bar is:

      - True  (SimForge itself added hydrogens): requires "high" -- this is
        exactly the mechanism behind the historical mis-hydrogenation bug
        (RDKit confidently fabricating the wrong hydrogens), so it must be
        airtight.
      - False (the input already had a complete hydrogen set and SimForge
        computed/added nothing): requires only "medium" -- SimForge is not
        the source of any potential error here, so a heuristic-level
        assessment of the pre-existing input is an acceptable bar, matching
        the tool's pre-existing trust model for already-complete ligands.

    ``element_assignment_confidence`` always requires "high" regardless --
    element normalization needs no RDKit and is unambiguous by construction
    (see ``normalize_ligand_elements``), so there is no reason to relax it.
    """
    reasons: list[str] = []
    required = "high" if hydrogenation_performed else "medium"
    required_rank = _CONFIDENCE_RANK[required]

    if hydrogenation_status != "complete":
        reasons.append(f"hydrogenation_status is '{hydrogenation_status}' (must be 'complete')")
    if _CONFIDENCE_RANK.get(element_assignment_confidence, 0) < _CONFIDENCE_RANK["high"]:
        reasons.append(f"element_assignment_confidence is '{element_assignment_confidence}' (must be 'high')")
    if _CONFIDENCE_RANK.get(connectivity_confidence, 0) < required_rank:
        reasons.append(f"connectivity_confidence is '{connectivity_confidence}' (must be '{required}')")
    if _CONFIDENCE_RANK.get(bond_order_confidence, 0) < required_rank:
        reasons.append(f"bond_order_confidence is '{bond_order_confidence}' (must be '{required}')")
    if _CONFIDENCE_RANK.get(hydrogenation_confidence, 0) < required_rank:
        reasons.append(f"hydrogenation_confidence is '{hydrogenation_confidence}' (must be '{required}')")

    status = "ready_for_ligpargen" if not reasons else "chemical_validation_required"
    return status, reasons
