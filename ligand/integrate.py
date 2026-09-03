"""
Protein–ligand GROMACS system assembly from LigParGen outputs.

Given a prepared protein PDB and LigParGen-generated .itp/.gro files, this
module builds a complete GROMACS-ready system without manual topology editing.
Multichain proteins (dimers, etc.) are fully supported in two modes:

  Embedded topology mode (monomer / inline multichain):
    pdb2gmx embeds all [ moleculetype ] blocks directly in topol.top.
    SimForge extracts them into protein.itp and generates a new master top.

  Master topology injection mode (pdb2gmx split-chain output):
    pdb2gmx writes per-chain topologies to topol_Protein_chain_A.itp, etc.
    topol.top itself has no inline [ moleculetype ].  SimForge detects this
    and patches the existing topol.top in-place with ligand includes and
    molecule entries, leaving per-chain files untouched.

Steps performed by assemble_system():
  1. gmx pdb2gmx on protein → protein.gro, raw topol.top, posre.itp
     1a. First attempt without -ignh.
     1b. On failure, inspect output: if GROMACS diagnoses a hydrogen/protonation
         mismatch (e.g. "Option -ignh will ignore all hydrogens in the input"),
         retry with -ignh. Any other failure is fatal immediately.
  2. Detect topology mode (embedded vs master_injection).
  3-4. Extract [ atomtypes ] block from ligand .itp → ligand_atomtypes.itp
  5. Write cleaned ligand .itp (atomtypes removed) to out_dir.
  6. Assemble final topol.top (mode-dependent):
       embedded:          split protein into protein.itp → generate new topol.top
       master_injection:  inject ligand includes/molecules into existing topol.top
  7. Merge protein.gro + ligand.gro → complex.gro
  8. Validate the assembled system
  9. Write assembly_report.yaml

Public API:
  Pdb2gmxMeta                                  — dataclass with both-attempt metadata
  extract_atomtypes_section(itp_text)          -> tuple[str, str]
  remove_atomtypes_from_itp(itp_path, out)     -> None
  parse_molecules_section(topol_path)          -> list[tuple[str, int]]
  detect_chain_itp_includes(topol_path)        -> list[str]
  inject_ligand_into_master_topol(...)         -> None
  split_topol_top(...)                         -> tuple[list[tuple[str, int]], str]
  generate_topol_top(...)                      -> Path
  merge_gro_files(...)                         -> int
  validate_system(...)                         -> list[str]
  run_pdb2gmx(protein_pdb, out_dir, ...)       -> tuple[Path, Path, Path, Pdb2gmxMeta]
  assemble_system(...)                         -> AssemblyResult

Multi-component (protein + N parameterized non-protein molecules — ligand,
cofactor, substrate, inhibitor, ...) extensions. ``assemble_system`` above is
a single-component wrapper around ``assemble_system_multi``; all other
single-component functions above are unmodified and remain independently
usable/tested.

  ParameterizedMolecule                        — one non-protein component (itp+gro+id+role)
  resolve_and_merge_atomtypes(...)             -> AtomtypeMergeResult
  inject_components_into_master_topol(...)     -> None
  generate_topol_top_multi(...)                -> Path
  generate_topol_top_preparameterized(...)     -> Path
  merge_gro_files_multi(...)                   -> int
  validate_system_multi(...)                   -> list[str]
  assemble_system_multi(...)                   -> AssemblyResult

Protein input modes (assemble_system_multi, keyword-only)
-----------------------------------------------------------
  MODE A (protein_pdb=...):  raw PDB -> gmx pdb2gmx -> protein topology.
  MODE B (protein_gro=..., protein_topology_itps=[...],
          protein_posre_itps=[...] optional):  an already-parameterized
  protein (pdb2gmx already run upstream) -> pdb2gmx is never invoked.

  A pre-parameterized protein is a first-class input, not a workaround:
  reconstructing a raw PDB from an already-processed GRO does not reliably
  round-trip through pdb2gmx (post-processing can rename atoms in ways the
  force field's .rtp templates no longer recognize, e.g. "Atom O1 in residue
  ASN ... was not found in rtp entry ASN"). When the caller already has
  protein coordinates plus their chain topology .itp file(s) from a prior
  pdb2gmx run, Mode B consumes them directly instead of forcing a second,
  unsafe pdb2gmx pass. The two modes are mutually exclusive.

Cross-component atomtype handling
----------------------------------
Independently-run LigParGen submissions each restart their local OPLS type
numbering from scratch (opls_800, opls_801, ...). Two components parameterized
in separate submissions can therefore reuse the same type *name* for two
chemically different atom types. ``resolve_and_merge_atomtypes`` detects this:
identical same-name definitions are deduplicated; conflicting same-name
definitions are namespaced per component (``opls_800`` -> ``A1_opls_800`` /
``COA_opls_800``) and every reference to the renamed name inside that
component's own topology is rewritten consistently. Charges, masses, sigma,
epsilon, bonded parameters, atom numbering, and coordinates are never touched
-- only the symbolic atomtype identifier changes. A conflict that cannot be
safely namespaced raises ``ValueError`` rather than silently picking one
definition.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import yaml

from utils.gro_parser import GroAtom, GroFile, parse_gro, write_gro
from utils.itp_parser import ItpFile, parse_itp

# Regex for extracting the filename from an #include "..." or #include '...' line.
_INCLUDE_RE = re.compile(r'#include\s+["\'](.+?)["\']')

# Molecule names that indicate solvent or ions in the [ molecules ] section.
_SOLVENT_ION_NAMES: frozenset[str] = frozenset({
    "SOL", "HOH", "WAT", "TIP3", "SPC",
    "NA", "CL", "K", "MG", "CA", "ZN", "FE", "MN", "CU",
    "NA+", "CL-", "ION",
})

# Map common shorthand force-field names to their GROMACS ff directory names.
_FF_DIR: dict[str, str] = {
    "oplsaa": "oplsaa.ff",
    "opls-aa": "oplsaa.ff",
    "charmm36": "charmm36-jul2022.ff",
    "amber99sb": "amber99sb.ff",
    "amber99sb-ildn": "amber99sb-ildn.ff",
    "gromos54a7": "gromos54a7.ff",
}


# ── pdb2gmx metadata ─────────────────────────────────────────────────────────

@dataclass
class Pdb2gmxMeta:
    """Records both the initial pdb2gmx attempt and any -ignh fallback."""
    initial_exit_code: int = -1
    initial_command: list[str] = field(default_factory=list)
    initial_error_summary: str = ""
    fallback_used: bool = False
    fallback_reason: str = ""
    fallback_command: list[str] = field(default_factory=list)
    final_exit_code: int = -1


# ── Multi-component molecule descriptor ────────────────────────────────────────

@dataclass
class ParameterizedMolecule:
    """
    One already-parameterized non-protein molecule to fold into the system:
    a ligand, cofactor, substrate, inhibitor, or any other small molecule with
    its own LigParGen (or equivalent) .itp/.gro pair.

    ``component_id`` identifies the component for reporting and for namespacing
    atomtype names on conflict (e.g. "A1", "COA") — it is independent of, and
    not assumed to equal, the GROMACS molecule name inside the .itp's
    [ moleculetype ] section (that name is parsed from the file, never
    inferred from the filename or from ``component_id``).

    ``role`` is metadata only ("ligand" | "cofactor" | "substrate" |
    "inhibitor" | ... or None) — the assembler treats every component
    identically regardless of role.
    """
    component_id: str
    topology_path: Path
    coordinate_path: Path
    role: Optional[str] = None


# ── Result type ───────────────────────────────────────────────────────────────

@dataclass
class AssemblyResult:
    success: bool
    protein_gro: Optional[Path] = None
    protein_itp: Optional[Path] = None  # None in master_injection/preparameterized modes
    ligand_atomtypes_itp: Optional[Path] = None
    ligand_itp: Optional[Path] = None
    topol_top: Optional[Path] = None
    complex_gro: Optional[Path] = None
    report_path: Optional[Path] = None
    topology_mode: str = "embedded"  # "embedded" | "master_injection"
    protein_topology_includes: list[str] = field(default_factory=list)
    protein_molecules: list[dict] = field(default_factory=list)  # [{name, count}, ...]
    protein_mol_name: str = ""  # first chain name, for backwards compatibility
    ligand_mol_name: str = ""
    protein_atom_count: int = 0
    ligand_atom_count: int = 0
    total_atom_count: int = 0
    pdb2gmx_fallback_used: bool = False
    pdb2gmx_meta: Optional[Pdb2gmxMeta] = None
    warnings: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    error: Optional[str] = None
    # Populated by assemble_system_multi (and, for the single-component case,
    # by the assemble_system wrapper) — one entry per ParameterizedMolecule:
    # {id, molecule_name, atom_count, count, role, topology, coordinates,
    #  atomtype_renames}. ligand_* fields above stay populated only when
    # exactly one component is present (back-compat with assemble_system).
    components: list[dict] = field(default_factory=list)


# ── ITP atomtypes extraction ──────────────────────────────────────────────────

def extract_atomtypes_section(itp_text: str) -> tuple[str, str]:
    """
    Extract the [ atomtypes ] section from raw ITP text.

    Returns:
        (atomtypes_block, itp_without_atomtypes)

    If no [ atomtypes ] section exists, returns ("", itp_text) unchanged.
    """
    lines = itp_text.splitlines(keepends=True)

    start: Optional[int] = None
    end: Optional[int] = None
    in_atomtypes = False

    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("[") and not stripped.startswith(";"):
            section = stripped.strip("[]").strip().lower()
            if section == "atomtypes":
                start = i
                in_atomtypes = True
            elif in_atomtypes:
                end = i
                in_atomtypes = False

    if start is None:
        return "", itp_text

    end_idx = end if end is not None else len(lines)
    atomtypes_block = "".join(lines[start:end_idx])
    remaining = "".join(lines[:start] + lines[end_idx:])
    return atomtypes_block, remaining


def remove_atomtypes_from_itp(
    itp_path: str | Path, out_path: str | Path
) -> None:
    """Write itp_path to out_path with the [ atomtypes ] section stripped."""
    _, cleaned = extract_atomtypes_section(Path(itp_path).read_text())
    Path(out_path).write_text(cleaned)


# ── Cross-component atomtype namespace resolution ──────────────────────────────

_VALID_ID_RE = re.compile(r"[^A-Za-z0-9_]")


def _sanitize_component_id(raw: str) -> str:
    """Return ``raw`` made safe for use as a GROMACS atomtype-name prefix."""
    s = _VALID_ID_RE.sub("_", raw.strip())
    if not s:
        s = "MOL"
    if s[0].isdigit():
        s = f"C_{s}"
    return s


@dataclass
class _AtomtypeEntry:
    name: str
    params: str  # everything after the name, whitespace-joined, comment stripped


def _parse_atomtypes_block(block_text: str) -> list[_AtomtypeEntry]:
    """
    Parse an [ atomtypes ] block's data lines into (name, params) entries.

    Raises ValueError on a line that cannot be split into a name plus at
    least one parameter — such a line cannot be safely deduplicated or
    conflict-checked, so it is treated as unsupported rather than guessed at.
    """
    entries: list[_AtomtypeEntry] = []
    for line in block_text.splitlines():
        s = line.split(";")[0].strip()
        if not s or s.startswith("["):
            continue
        parts = s.split()
        if len(parts) < 2:
            raise ValueError(
                f"Malformed [ atomtypes ] line (cannot isolate a name and "
                f"parameters): {line!r}"
            )
        entries.append(_AtomtypeEntry(name=parts[0], params=" ".join(parts[1:])))
    return entries


def _params_equal(a: str, b: str, tol: float = 1e-6) -> bool:
    """True if two atomtype parameter strings are equivalent (numeric tolerance)."""
    ta, tb = a.split(), b.split()
    if len(ta) != len(tb):
        return False
    for x, y in zip(ta, tb):
        try:
            fx, fy = float(x), float(y)
        except ValueError:
            if x != y:
                return False
            continue
        if abs(fx - fy) > tol * max(1.0, abs(fx), abs(fy)):
            return False
    return True


@dataclass
class AtomtypeMergeResult:
    combined_block: str                        # final "[ atomtypes ]\n..." text
    renames: dict[str, dict[str, str]]          # component_id -> {old_name: new_name}
    warnings: list[str] = field(default_factory=list)


def resolve_and_merge_atomtypes(
    components: list[tuple[str, str]],
) -> AtomtypeMergeResult:
    """
    Merge N components' [ atomtypes ] blocks into one, namespacing conflicts.

    ``components`` is an ordered list of (component_id, atomtypes_block_text)
    — block_text may be "" for a component with no [ atomtypes ] section.

    Resolution rules, applied per atomtype name shared by 2+ components:
      - identical parameters everywhere it's defined -> deduplicated, kept
        under its original name (only one definition is emitted).
      - different parameters somewhere -> a genuine conflict. The type is
        renamed in *every* component that defines it (not just the
        "losing" one) to ``<sanitized_component_id>_<original_name>``, and
        each renamed component's own definition is kept under its new,
        now-unique name. This is deterministic regardless of component order.

    A generated namespaced name that itself collides with an unrelated,
    differently-parameterized definition (only possible if a component's own
    id text is adversarially chosen to collide) raises ValueError rather than
    silently overwriting either definition.

    Names that appear in only one component are passed through unchanged.
    """
    per_component_entries: dict[str, list[_AtomtypeEntry]] = {}
    for cid, block in components:
        if not block.strip():
            per_component_entries[cid] = []
            continue
        try:
            per_component_entries[cid] = _parse_atomtypes_block(block)
        except ValueError as exc:
            raise ValueError(f"[{cid}] {exc}") from exc

    by_name: dict[str, list[tuple[str, _AtomtypeEntry]]] = {}
    for cid, entries in per_component_entries.items():
        for e in entries:
            by_name.setdefault(e.name, []).append((cid, e))

    conflicting_names: set[str] = set()
    for name, occs in by_name.items():
        if len(occs) < 2:
            continue
        first_params = occs[0][1].params
        if not all(_params_equal(first_params, e.params) for _, e in occs[1:]):
            conflicting_names.add(name)

    renames: dict[str, dict[str, str]] = {cid: {} for cid, _ in components}
    seen_params: dict[str, str] = {}  # final_name -> params already emitted under it
    final_entries: list[_AtomtypeEntry] = []

    for cid, _block in components:
        for e in per_component_entries[cid]:
            if e.name in conflicting_names:
                new_name = f"{_sanitize_component_id(cid)}_{e.name}"
                if new_name in seen_params and not _params_equal(seen_params[new_name], e.params):
                    raise ValueError(
                        f"Atomtype namespace collision: generated name "
                        f"'{new_name}' for component '{cid}' atomtype "
                        f"'{e.name}' collides with an existing, differently "
                        f"parameterized definition. Cannot safely resolve "
                        f"automatically — rename component '{cid}' or its "
                        f"conflicting atomtype and re-run."
                    )
                renames[cid][e.name] = new_name
                final_name = new_name
            else:
                final_name = e.name

            if final_name in seen_params:
                if not _params_equal(seen_params[final_name], e.params):
                    raise ValueError(
                        f"Unresolvable atomtype conflict for '{final_name}' "
                        f"(component '{cid}')."
                    )
                continue  # identical duplicate — already emitted
            seen_params[final_name] = e.params
            final_entries.append(_AtomtypeEntry(name=final_name, params=e.params))

    if final_entries:
        lines = ["[ atomtypes ]"]
        for e in final_entries:
            lines.append(f"  {e.name:<14s} {e.params}")
        combined_block = "\n".join(lines) + "\n"
    else:
        combined_block = ""

    return AtomtypeMergeResult(combined_block=combined_block, renames=renames)


def _apply_atomtype_renames(itp_text: str, renames: dict[str, str]) -> str:
    """
    Rewrite every whole-word reference to a renamed atomtype in ``itp_text``.

    Word-boundary matching means only the exact symbolic identifier is
    touched — numeric fields, atom indices, and unrelated identifiers that
    merely contain the old name as a substring are left untouched. In
    practice this affects the atom-type column of [ atoms ]; it is applied
    to the whole remaining text (not just [ atoms ]) so any bonded/pair type
    table that references atomtype names by name (e.g. [ pairtypes ],
    [ nonbond_params ]) would also be caught, though standard per-molecule
    LigParGen output does not use those sections.
    """
    if not renames:
        return itp_text
    for old_name in sorted(renames, key=len, reverse=True):
        pattern = re.compile(r"\b" + re.escape(old_name) + r"\b")
        itp_text = pattern.sub(renames[old_name], itp_text)
    return itp_text


# ── topol.top splitting ───────────────────────────────────────────────────────

def parse_molecules_section(topol_path: str | Path) -> list[tuple[str, int]]:
    """
    Parse the [ molecules ] section from a pdb2gmx-generated topol.top.

    Returns list of (molecule_name, count) tuples in order.  For a fresh
    pdb2gmx run (no solvation), this contains only the protein chain entries,
    e.g. [("Protein_chain_A", 1), ("Protein_chain_B", 1)] for a dimer.
    Returns an empty list when no [ molecules ] section is found.
    """
    lines = Path(topol_path).read_text().splitlines()
    in_molecules = False
    entries: list[tuple[str, int]] = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("[") and not stripped.startswith(";"):
            section = stripped.strip("[]").strip().lower()
            in_molecules = section == "molecules"
            continue
        if in_molecules:
            if not stripped or stripped.startswith(";"):
                continue
            parts = stripped.split()
            if len(parts) >= 2:
                try:
                    entries.append((parts[0], int(parts[1])))
                except ValueError:
                    pass
    return entries


def detect_chain_itp_includes(topol_path: str | Path) -> list[str]:
    """
    Find per-chain ITP filenames referenced in a pdb2gmx master topol.top.

    Returns an ordered list of basenames such as
    ``['topol_Protein_chain_A.itp', 'topol_Protein_chain_B.itp']``.
    Returns an empty list for embedded (inline moleculetype) topologies.

    Detection criterion: an ``#include`` whose quoted filename starts with
    ``topol_``, ends with ``.itp``, and contains no directory separator.
    """
    found: list[str] = []
    for line in Path(topol_path).read_text().splitlines():
        s = line.strip()
        if not s.startswith("#include"):
            continue
        m = _INCLUDE_RE.match(s)
        if m:
            fname = m.group(1)
            if "/" not in fname and fname.startswith("topol_") and fname.endswith(".itp"):
                found.append(fname)
    return found


def _scan_master_topol_insertion_points(orig_lines: list[str]) -> tuple[int, int, int]:
    """
    Locate the three positions needed to inject includes/molecules into a
    pdb2gmx master topol.top: the forcefield #include line, the last
    per-chain topol_*.itp #include line, and the line index before which new
    [ molecules ] entries should be inserted (before the first solvent/ion
    entry, else after the last protein entry, else at end of file).

    Raises:
        ValueError if the forcefield include or a per-chain include are not
        found (indicates this is not a valid master-mode topol.top).
    """
    ff_idx: Optional[int] = None
    last_chain_idx: Optional[int] = None
    last_protein_mol_idx: Optional[int] = None
    first_solvent_mol_idx: Optional[int] = None
    in_molecules = False

    for i, line in enumerate(orig_lines):
        s = line.strip()

        if s.startswith("[") and not s.startswith(";"):
            in_molecules = s.strip("[]").strip().lower() == "molecules"

        if in_molecules and s and not s.startswith(";") and not s.startswith("["):
            parts = s.split()
            if len(parts) >= 2:
                try:
                    int(parts[1])
                    name = parts[0]
                    if name.upper() in _SOLVENT_ION_NAMES:
                        if first_solvent_mol_idx is None:
                            first_solvent_mol_idx = i
                    else:
                        last_protein_mol_idx = i
                except ValueError:
                    pass

        if not s.startswith("#include"):
            continue
        m = _INCLUDE_RE.match(s)
        if not m:
            continue
        fname = m.group(1)
        if "forcefield.itp" in fname:
            ff_idx = i
        elif "/" not in fname and fname.startswith("topol_") and fname.endswith(".itp"):
            last_chain_idx = i

    if ff_idx is None:
        raise ValueError(
            "No forcefield #include found in topol.top — "
            "cannot determine where to insert ligand atomtypes."
        )
    if last_chain_idx is None:
        raise ValueError(
            "No topol_*.itp per-chain #include found — "
            "this does not appear to be a pdb2gmx master topology."
        )

    if first_solvent_mol_idx is not None:
        mol_insert_before = first_solvent_mol_idx   # insert before first solvent
    elif last_protein_mol_idx is not None:
        mol_insert_before = last_protein_mol_idx + 1  # insert after last protein
    else:
        mol_insert_before = len(orig_lines)            # append

    return ff_idx, last_chain_idx, mol_insert_before


def inject_ligand_into_master_topol(
    topol_path: str | Path,
    ligand_itp_name: str,
    ligand_atomtypes_itp_name: str,
    ligand_mol_name: str,
) -> None:
    """
    Patch a pdb2gmx master topol.top in-place to include the ligand.

    Three modifications are applied:

    1. Insert ``#include "<ligand_atomtypes_itp_name>"`` on the line
       immediately after the forcefield ``#include``.
    2. Insert ``#include "<ligand_itp_name>"`` on the line immediately after
       the last per-chain topology include (files matching ``topol_*.itp``).
    3. Append the ligand molecule entry to the ``[ molecules ]`` section,
       after all protein entries and before any solvent/ion entries.

    Raises:
        ValueError if the forcefield include or a per-chain include are not
        found (indicates this is not a valid master-mode topol.top).
    """
    topol_path = Path(topol_path)
    orig_lines = topol_path.read_text().splitlines()

    ff_idx, last_chain_idx, mol_insert_before = _scan_master_topol_insertion_points(orig_lines)

    # Build the modified file in a single streaming pass.
    new_lines: list[str] = []
    for i, line in enumerate(orig_lines):
        new_lines.append(line)
        if i == ff_idx:
            new_lines.append(f'#include "{ligand_atomtypes_itp_name}"')
        if i == last_chain_idx:
            new_lines.append(f'#include "{ligand_itp_name}"')

    # Adjust mol_insert_before for lines inserted before it.
    shifts = int(ff_idx < mol_insert_before) + int(last_chain_idx < mol_insert_before)
    adjusted_mol_idx = mol_insert_before + shifts

    lig_entry = f"{ligand_mol_name}            1"
    new_lines = (
        new_lines[:adjusted_mol_idx] + [lig_entry] + new_lines[adjusted_mol_idx:]
    )

    topol_path.write_text("\n".join(new_lines) + "\n")


def inject_components_into_master_topol(
    topol_path: str | Path,
    atomtypes_itp_name: str,
    components: list[tuple[str, str, int]],
) -> None:
    """
    Patch a pdb2gmx master topol.top in-place to include N components.

    ``components`` is an ordered list of (itp_filename, molecule_name, count).
    Modifications mirror ``inject_ligand_into_master_topol`` but generalized:

    1. Insert a single ``#include "<atomtypes_itp_name>"`` after the
       forcefield #include (one combined atomtypes file for all components).
    2. Insert one ``#include "<itp_filename>"`` per component, in order,
       immediately after the last per-chain topology include.
    3. Append one molecule entry per component to [ molecules ], in order,
       after all protein entries and before any solvent/ion entries.

    Raises:
        ValueError if the forcefield include or a per-chain include are not
        found (indicates this is not a valid master-mode topol.top).
    """
    topol_path = Path(topol_path)
    orig_lines = topol_path.read_text().splitlines()

    ff_idx, last_chain_idx, mol_insert_before = _scan_master_topol_insertion_points(orig_lines)

    new_lines: list[str] = []
    for i, line in enumerate(orig_lines):
        new_lines.append(line)
        if i == ff_idx:
            new_lines.append(f'#include "{atomtypes_itp_name}"')
        if i == last_chain_idx:
            for itp_name, _mol_name, _count in components:
                new_lines.append(f'#include "{itp_name}"')

    n_inserted_before_mol = int(ff_idx < mol_insert_before)
    if last_chain_idx < mol_insert_before:
        n_inserted_before_mol += len(components)
    adjusted_mol_idx = mol_insert_before + n_inserted_before_mol

    mol_entries = [f"{mol_name}            {count}" for _itp_name, mol_name, count in components]
    new_lines = new_lines[:adjusted_mol_idx] + mol_entries + new_lines[adjusted_mol_idx:]

    topol_path.write_text("\n".join(new_lines) + "\n")


def split_topol_top(
    topol_path: str | Path,
    out_dir: str | Path,
    protein_itp_name: str = "protein.itp",
) -> tuple[list[tuple[str, int]], str]:
    """
    Read a pdb2gmx-generated topol.top and extract the protein topology.

    Writes ``protein_itp_name`` to ``out_dir`` containing all [ section ]
    blocks that belong to the protein (from the first [ moleculetype ] up to,
    but not including, the first water/ions .ff/ include).  For multichain
    proteins, all chain topologies are included in a single protein.itp.

    Returns:
        (protein_molecules, water_ions_block)
        - protein_molecules: list of (mol_name, count) from pdb2gmx [ molecules ]
          — for a dimer: [("Protein_chain_A", 1), ("Protein_chain_B", 1)]
        - water_ions_block: raw text from the water include to just before
          [ system ] — inserted verbatim into the new topol.top

    Raises:
        ValueError if required sections are not found.
    """
    topol_path = Path(topol_path)
    out_dir = Path(out_dir)
    lines = topol_path.read_text().splitlines(keepends=True)

    ff_include_idx: Optional[int] = None
    first_mol_idx: Optional[int] = None
    water_include_idx: Optional[int] = None
    system_idx: Optional[int] = None

    for i, line in enumerate(lines):
        s = line.strip()

        if ff_include_idx is None and s.startswith("#include") and "forcefield.itp" in s:
            ff_include_idx = i
            continue

        if ff_include_idx is not None and first_mol_idx is None:
            if _is_section(line) and _section_name(line) == "moleculetype":
                first_mol_idx = i
                continue

        if first_mol_idx is not None and water_include_idx is None:
            # Water/ions include: #include "*.ff/*.itp" (not the forcefield one)
            if s.startswith("#include") and ".ff/" in s and "forcefield.itp" not in s:
                water_include_idx = i
                continue

        if water_include_idx is not None and system_idx is None:
            if _is_section(line) and _section_name(line) == "system":
                system_idx = i

    if first_mol_idx is None:
        raise ValueError(
            f"No [ moleculetype ] section found in {topol_path}. "
            "Is this a valid pdb2gmx-generated topol.top?"
        )
    if water_include_idx is None:
        raise ValueError(
            f"No water/ions .itp #include found in {topol_path}. "
            "Cannot determine where the protein topology ends."
        )

    # Protein topology body: from first [ moleculetype ] to just before water include
    # For multichain proteins this spans all chain topologies in order.
    protein_lines = lines[first_mol_idx:water_include_idx]

    # Water+ions block: from water include to just before [ system ]
    water_end = system_idx if system_idx is not None else len(lines)
    water_ions_block = "".join(lines[water_include_idx:water_end]).rstrip()

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / protein_itp_name).write_text("".join(protein_lines))

    # Use [ molecules ] as the authoritative source for which protein molecules exist.
    protein_molecules = parse_molecules_section(topol_path)
    if not protein_molecules:
        # Fallback when [ molecules ] is absent: derive from first [ moleculetype ]
        mol_name = _mol_name_from_lines(protein_lines)
        protein_molecules = [(mol_name, 1)]

    return protein_molecules, water_ions_block


def _is_section(line: str) -> bool:
    s = line.strip()
    return s.startswith("[") and not s.startswith(";")


def _section_name(line: str) -> str:
    return line.strip().strip("[]").strip().lower()


def _mol_name_from_lines(lines: list[str]) -> str:
    """Return the molecule name from a [ moleculetype ] block."""
    in_mol = False
    for line in lines:
        if _is_section(line):
            in_mol = _section_name(line) == "moleculetype"
            continue
        if in_mol:
            s = line.strip()
            if not s or s.startswith(";"):
                continue
            return s.split()[0]
    return "Protein"


# ── topol.top generation ──────────────────────────────────────────────────────

def generate_topol_top(
    out_dir: str | Path,
    forcefield: str,
    water_model: str,
    protein_itp_name: str,
    ligand_itp_name: str,
    ligand_atomtypes_itp_name: str,
    protein_molecules: list[tuple[str, int]],
    ligand_mol_name: str,
    water_ions_block: str = "",
) -> Path:
    """
    Write a clean master topol.top with the correct GROMACS include order.

    ``protein_molecules`` is a list of (molecule_name, count) tuples exactly
    as produced by pdb2gmx's [ molecules ] section.  For a monomer pass
    [("Protein_chain_A", 1)]; for a dimer pass
    [("Protein_chain_A", 1), ("Protein_chain_B", 1)].

    Include order:
      1. #include "ff/forcefield.itp"
      2. #include ligand_atomtypes_itp  ← must precede ligand itp
      3. #include protein_itp
      4. #include ligand_itp
      5. water + ions block (verbatim from pdb2gmx output, or sensible default)
      6. [ system ] / [ molecules ]
    """
    out_dir = Path(out_dir)
    ff = forcefield.lower()
    ff_dir = _FF_DIR.get(ff, f"{ff}.ff")

    if not water_ions_block:
        wm = water_model.lower()
        water_ions_block = (
            f'#include "{ff_dir}/{wm}.itp"\n'
            "\n"
            "#ifdef POSRES_WATER\n"
            "; Position restraint for each water molecule\n"
            "[ position_restraints ]\n"
            "; i funct       fcx        fcy        fcz\n"
            "  1    1       1000       1000       1000\n"
            "#endif\n"
            "\n"
            f'#include "{ff_dir}/ions.itp"'
        )

    mol_entries = [f"{name}           {count}" for name, count in protein_molecules]
    mol_entries.append(f"{ligand_mol_name}            1")

    sections = [
        "; Generated by SimForge ligand integrate",
        "",
        "; Force-field parameters",
        f'#include "{ff_dir}/forcefield.itp"',
        "",
        "; Ligand atomtypes (must come before ligand topology)",
        f'#include "{ligand_atomtypes_itp_name}"',
        "",
        "; Protein topology",
        f'#include "{protein_itp_name}"',
        "",
        "; Ligand topology",
        f'#include "{ligand_itp_name}"',
        "",
        water_ions_block,
        "",
        "[ system ]",
        "Protein-Ligand Complex",
        "",
        "[ molecules ]",
        "; Compound        #mols",
    ] + mol_entries

    top_path = out_dir / "topol.top"
    top_path.write_text("\n".join(sections) + "\n")
    return top_path


def _default_water_ions_block(ff_dir: str, water_model: str) -> str:
    wm = water_model.lower()
    return (
        f'#include "{ff_dir}/{wm}.itp"\n'
        "\n"
        "#ifdef POSRES_WATER\n"
        "; Position restraint for each water molecule\n"
        "[ position_restraints ]\n"
        "; i funct       fcx        fcy        fcz\n"
        "  1    1       1000       1000       1000\n"
        "#endif\n"
        "\n"
        f'#include "{ff_dir}/ions.itp"'
    )


def _render_multi_component_topol_top(
    out_dir: Path,
    forcefield: str,
    water_model: str,
    protein_include_names: list[str],
    atomtypes_itp_name: str,
    protein_molecules: list[tuple[str, int]],
    components: list[tuple[str, str, int]],
    water_ions_block: str = "",
    protein_section_comment: str = "; Protein topology",
) -> Path:
    """
    Shared topol.top renderer for both generate_topol_top_multi (single
    protein.itp, post-pdb2gmx-split) and generate_topol_top_preparameterized
    (N pre-existing chain topology files, pdb2gmx never run). Only the
    protein include list and its section comment differ between the two.
    """
    ff = forcefield.lower()
    ff_dir = _FF_DIR.get(ff, f"{ff}.ff")

    if not water_ions_block:
        water_ions_block = _default_water_ions_block(ff_dir, water_model)

    mol_entries = [f"{name}           {count}" for name, count in protein_molecules]
    mol_entries += [f"{mol_name}            {count}" for _itp, mol_name, count in components]

    sections = [
        "; Generated by SimForge ligand integrate",
        "",
        "; Force-field parameters",
        f'#include "{ff_dir}/forcefield.itp"',
        "",
        "; Component atomtypes (must come before component topologies)",
        f'#include "{atomtypes_itp_name}"',
        "",
        protein_section_comment,
    ]
    for name in protein_include_names:
        sections.append(f'#include "{name}"')
    sections += ["", "; Non-protein component topologies"]
    for itp_name, _mol_name, _count in components:
        sections.append(f'#include "{itp_name}"')
    sections += [
        "",
        water_ions_block,
        "",
        "[ system ]",
        "Protein-Component Complex",
        "",
        "[ molecules ]",
        "; Compound        #mols",
    ] + mol_entries

    top_path = out_dir / "topol.top"
    top_path.write_text("\n".join(sections) + "\n")
    return top_path


def generate_topol_top_multi(
    out_dir: str | Path,
    forcefield: str,
    water_model: str,
    protein_itp_name: str,
    atomtypes_itp_name: str,
    protein_molecules: list[tuple[str, int]],
    components: list[tuple[str, str, int]],
    water_ions_block: str = "",
) -> Path:
    """
    Write a clean master topol.top for a protein + N-component system.

    ``components`` is an ordered list of (itp_filename, molecule_name, count).
    Include order mirrors ``generate_topol_top``, generalized to N components:

      1. #include "ff/forcefield.itp"
      2. #include atomtypes_itp_name       ← one combined include for all components
      3. #include protein_itp
      4. #include <component itp> for each component, in order
      5. water + ions block
      6. [ system ] / [ molecules ]
    """
    return _render_multi_component_topol_top(
        Path(out_dir), forcefield, water_model, [protein_itp_name],
        atomtypes_itp_name, protein_molecules, components, water_ions_block,
    )


def generate_topol_top_preparameterized(
    out_dir: str | Path,
    forcefield: str,
    water_model: str,
    protein_topology_names: list[str],
    atomtypes_itp_name: str,
    protein_molecules: list[tuple[str, int]],
    components: list[tuple[str, str, int]],
    water_ions_block: str = "",
) -> Path:
    """
    Write a master topol.top for a pre-parameterized protein (pdb2gmx never
    run) + N non-protein components.

    ``protein_topology_names`` is the ordered list of already-existing chain
    topology .itp filenames (copied verbatim into out_dir by the caller) —
    one per protein chain, e.g. ``["topol_Protein_chain_A.itp", ...]``. Each
    file's own position-restraint #include (already embedded by pdb2gmx when
    it originally generated these files) is preserved as-is; no new
    position-restraint include is generated here.

    Include order:
      1. #include "ff/forcefield.itp"
      2. #include atomtypes_itp_name
      3. #include <protein chain topology> for each chain, in order
      4. #include <component itp> for each component, in order
      5. water + ions block
      6. [ system ] / [ molecules ]
    """
    return _render_multi_component_topol_top(
        Path(out_dir), forcefield, water_model, list(protein_topology_names),
        atomtypes_itp_name, protein_molecules, components, water_ions_block,
        protein_section_comment="; Protein topology (pre-parameterized — pdb2gmx not run)",
    )


# ── GRO merging ───────────────────────────────────────────────────────────────

def merge_gro_files(
    protein_gro: str | Path,
    ligand_gro: str | Path,
    out_path: str | Path,
    title: str = "Protein-Ligand Complex",
) -> int:
    """
    Merge protein and ligand .gro files into a single complex.gro.

    - Box vectors come from the protein .gro.
    - Atom serial numbers are renumbered sequentially.
    - Residue numbers are preserved as-is from each source file.

    Returns:
        Total atom count written (protein + ligand).
    """
    return merge_gro_files_multi(protein_gro, [ligand_gro], out_path, title=title)


def merge_gro_files_multi(
    protein_gro: str | Path,
    component_gros: list[str | Path],
    out_path: str | Path,
    title: str = "Protein-Ligand Complex",
) -> int:
    """
    Merge a protein .gro and N component .gro files into a single complex.gro.

    - Box vectors come from the protein .gro.
    - Atoms are concatenated in order: protein, then each component in the
      order given in ``component_gros``.
    - Atom serial numbers are renumbered sequentially across the whole file.
    - Residue numbers are preserved as-is from each source file.
    - Coordinates are copied verbatim — no transformation, translation, or
      re-centering is applied to any component.

    Returns:
        Total atom count written (protein + all components).
    """
    prot = parse_gro(protein_gro)

    merged: list[GroAtom] = list(prot.atoms)
    for gro_path in component_gros:
        merged.extend(parse_gro(gro_path).atoms)

    # Renumber atom indices sequentially (GROMACS wraps at 99999)
    for idx, atom in enumerate(merged, start=1):
        atom.atom_number = idx

    merged_gro = GroFile(title=title, atoms=merged, box=prot.box)
    write_gro(merged_gro, out_path)
    return len(merged)


# ── Validation ────────────────────────────────────────────────────────────────

def validate_system(
    out_dir: str | Path,
    protein_mol_name: str,
    ligand_mol_name: str,
    ligand_itp_path: str | Path,
    complex_gro_path: str | Path,
    topol_top_path: str | Path,
    ligand_atomtypes_itp_path: str | Path,
    run_grompp: bool = True,
) -> list[str]:
    """
    Validate the assembled system.

    Returns a list of error strings; an empty list means the system is valid.

    Checks:
      1. Ligand molecule name in .itp matches the name used in [ molecules ].
      2. Atom count header in complex.gro matches actual atom lines.
      3. ligand_atomtypes.itp exists on disk.
      4. topol.top includes ligand_atomtypes BEFORE the ligand .itp.
      5. (Optional) gmx grompp dry run if GROMACS is in PATH.
    """
    errors: list[str] = []

    # 1. Molecule name consistency
    try:
        itp = parse_itp(ligand_itp_path)
        itp_name = itp.moleculetype.name if itp.moleculetype else None
        if itp_name != ligand_mol_name:
            errors.append(
                f"Molecule name mismatch: ligand .itp has '{itp_name}' but "
                f"[ molecules ] expects '{ligand_mol_name}'."
            )
    except Exception as exc:
        errors.append(f"Cannot parse ligand .itp for name check: {exc}")

    # 2. Atom count in complex.gro
    try:
        gro = parse_gro(complex_gro_path)
        if gro.atom_count != len(gro.atoms):
            errors.append(
                f"complex.gro atom count mismatch: header declares {gro.atom_count} "
                f"but {len(gro.atoms)} atom lines are present."
            )
        if gro.atom_count == 0:
            errors.append("complex.gro contains no atoms.")
    except Exception as exc:
        errors.append(f"Cannot parse complex.gro: {exc}")

    # 3. ligand_atomtypes.itp present
    if not Path(ligand_atomtypes_itp_path).exists():
        errors.append(
            f"ligand_atomtypes.itp not found: {ligand_atomtypes_itp_path}"
        )

    # 4. Include order in topol.top
    try:
        top_text = Path(topol_top_path).read_text()
        at_name = Path(ligand_atomtypes_itp_path).name
        lig_name = Path(ligand_itp_path).name
        at_pos = top_text.find(at_name)
        lig_pos = top_text.find(lig_name)
        if at_pos == -1:
            errors.append(
                f"topol.top does not include '{at_name}'. "
                "Ligand atomtypes must be defined before the ligand topology."
            )
        elif lig_pos != -1 and at_pos > lig_pos:
            errors.append(
                f"topol.top include order error: '{at_name}' must appear "
                f"BEFORE '{lig_name}'."
            )
    except Exception as exc:
        errors.append(f"Cannot read topol.top for include-order check: {exc}")

    # 5. Optional gmx grompp dry run
    if run_grompp and shutil.which("gmx"):
        errors.extend(
            _grompp_check(Path(out_dir), Path(complex_gro_path), Path(topol_top_path))
        )

    return errors


def validate_system_multi(
    out_dir: str | Path,
    protein_mol_name: str,
    components: list[dict],
    complex_gro_path: str | Path,
    topol_top_path: str | Path,
    atomtypes_itp_path: str | Path,
    run_grompp: bool = True,
) -> list[str]:
    """
    Validate an assembled protein + N-component system.

    ``components`` is a list of dicts, one per component, each with keys
    ``component_id``, ``molecule_name``, ``itp_path`` (the cleaned/renamed
    component .itp actually written to out_dir).

    Returns a list of error strings; an empty list means the system is valid.

    Checks (generalizing ``validate_system`` to N components):
      1. Each component's molecule name in its (cleaned) .itp matches the
         name used in [ molecules ].
      2. Atom count header in complex.gro matches actual atom lines.
      3. component_atomtypes.itp exists on disk.
      4. topol.top includes the atomtypes file BEFORE every component .itp.
      5. (Optional) gmx grompp dry run if GROMACS is in PATH.
    """
    errors: list[str] = []

    # 1. Molecule name consistency, per component
    for comp in components:
        cid = comp["component_id"]
        try:
            itp = parse_itp(comp["itp_path"])
            itp_name = itp.moleculetype.name if itp.moleculetype else None
            if itp_name != comp["molecule_name"]:
                errors.append(
                    f"[{cid}] Molecule name mismatch: component .itp has "
                    f"'{itp_name}' but [ molecules ] expects "
                    f"'{comp['molecule_name']}'."
                )
        except Exception as exc:
            errors.append(f"[{cid}] Cannot parse component .itp for name check: {exc}")

    # 2. Atom count in complex.gro
    try:
        gro = parse_gro(complex_gro_path)
        if gro.atom_count != len(gro.atoms):
            errors.append(
                f"complex.gro atom count mismatch: header declares {gro.atom_count} "
                f"but {len(gro.atoms)} atom lines are present."
            )
        if gro.atom_count == 0:
            errors.append("complex.gro contains no atoms.")
    except Exception as exc:
        errors.append(f"Cannot parse complex.gro: {exc}")

    # 3. component_atomtypes.itp present
    if not Path(atomtypes_itp_path).exists():
        errors.append(f"Component atomtypes file not found: {atomtypes_itp_path}")

    # 4. Include order in topol.top
    try:
        top_text = Path(topol_top_path).read_text()
        at_name = Path(atomtypes_itp_path).name
        at_pos = top_text.find(at_name)
        if at_pos == -1:
            errors.append(
                f"topol.top does not include '{at_name}'. "
                "Component atomtypes must be defined before component topologies."
            )
        else:
            for comp in components:
                comp_itp_name = Path(comp["itp_path"]).name
                comp_pos = top_text.find(comp_itp_name)
                if comp_pos != -1 and at_pos > comp_pos:
                    errors.append(
                        f"topol.top include order error: '{at_name}' must "
                        f"appear BEFORE '{comp_itp_name}'."
                    )
    except Exception as exc:
        errors.append(f"Cannot read topol.top for include-order check: {exc}")

    # 5. Optional gmx grompp dry run
    if run_grompp and shutil.which("gmx"):
        errors.extend(
            _grompp_check(Path(out_dir), Path(complex_gro_path), Path(topol_top_path))
        )

    return errors


def _min_box_vector(gro: Path) -> "float | None":
    """Shortest box edge (nm) from a .gro file's last line, or None."""
    try:
        lines = [ln for ln in gro.read_text().splitlines() if ln.strip()]
        vals = [float(x) for x in lines[-1].split()[:3]] if len(lines) >= 3 else []
        pos = [v for v in vals if v > 0]
        return min(pos) if pos else None
    except Exception:
        return None


def _grompp_check(out_dir: Path, gro: Path, top: Path) -> list[str]:
    """gmx grompp as a strict topology syntax check (``-maxwarn 0``).

    Returns a list of strings. Fatal problems are returned verbatim (the caller
    treats them as errors). Non-fatal grompp warnings are returned prefixed
    ``[grompp-warning]`` so the caller can surface them without failing the
    build — SimForge never silently swallows a grompp warning.
    """
    mdp = out_dir / "_gmx_check.mdp"
    tpr = out_dir / "_gmx_check.tpr"
    # Scale the check cut-off to the (often tight, pre-solvation) box so grompp
    # does not reject "cut-off > half box" before it can check the topology.
    min_box = _min_box_vector(gro)
    rc = 1.0 if not min_box else max(0.30, min(1.0, round(0.45 * min_box, 3)))
    mdp.write_text(
        "integrator   = steep\n"
        "nsteps       = 0\n"
        "emtol        = 100\n"
        "nstlist      = 1\n"
        "cutoff-scheme = Verlet\n"
        f"rlist        = {rc}\n"
        f"rcoulomb     = {rc}\n"
        f"rvdw         = {rc}\n"
    )
    try:
        result = subprocess.run(
            [
                "gmx", "grompp",
                "-f", str(mdp),
                "-c", str(gro),
                "-p", str(top),
                "-o", str(tpr),
                "-maxwarn", "0",
            ],
            capture_output=True, text=True, cwd=out_dir, timeout=30,
        )
        combined = result.stdout + result.stderr
        warns = [ln.strip() for ln in combined.splitlines()
                 if ln.strip().lower().startswith("warning")]
        if result.returncode == 0:
            return [f"[grompp-warning] {w}" for w in warns]
        # Non-zero exit: re-run tolerantly to tell "only warnings" from "fatal".
        tol = subprocess.run(
            ["gmx", "grompp", "-f", str(mdp), "-c", str(gro), "-p", str(top),
             "-o", str(tpr), "-maxwarn", "10"],
            capture_output=True, text=True, cwd=out_dir, timeout=30,
        )
        if tol.returncode == 0:
            tol_all = tol.stdout + tol.stderr
            tw = [ln.strip() for ln in tol_all.splitlines()
                  if ln.strip().lower().startswith("warning")] or warns
            return [f"[grompp-warning] {w}" for w in tw]
        return [f"gmx grompp dry run failed:\n{result.stderr.strip()}"]
    except subprocess.TimeoutExpired:
        return ["gmx grompp dry run timed out (30 s)."]
    except Exception as exc:
        return [f"gmx grompp dry run raised an exception: {exc}"]
    finally:
        mdp.unlink(missing_ok=True)
        tpr.unlink(missing_ok=True)
    return []


# ── pdb2gmx H-mismatch detection ─────────────────────────────────────────────

# GROMACS prints these when the failure is caused by a hydrogen/protonation
# mismatch in the input PDB.  Either pattern is sufficient to trigger -ignh.
_IGNH_TRIGGER_PATTERNS: tuple[str, ...] = (
    "option -ignh will ignore all hydrogens in the input",
    "for a hydrogen, this can be a different protonation state",
)


def _should_retry_with_ignh(combined_output: str) -> tuple[bool, str]:
    """
    Inspect combined pdb2gmx stderr+stdout.

    Returns:
        (should_retry, reason_line)

    Returns True only when GROMACS explicitly diagnoses a hydrogen or
    protonation-state mismatch.  Missing residues, unknown residue names,
    missing force-field files, atom-type errors, and malformed PDB input are
    all excluded — they produce no -ignh suggestion from GROMACS.
    """
    lower = combined_output.lower()
    for pattern in _IGNH_TRIGGER_PATTERNS:
        if pattern in lower:
            for line in combined_output.splitlines():
                if pattern in line.lower():
                    return True, line.strip()
            return True, pattern
    return False, ""


# ── GROMACS pdb2gmx wrapper ───────────────────────────────────────────────────

def _build_pdb2gmx_cmd(
    protein_pdb: Path,
    protein_gro: Path,
    topol_top: Path,
    posre_itp: Path,
    forcefield: str,
    water_model: str,
    ignh: bool = False,
) -> list[str]:
    cmd = [
        "gmx", "pdb2gmx",
        "-f", str(protein_pdb),
        "-o", str(protein_gro),
        "-p", str(topol_top),
        "-i", str(posre_itp),
        "-ff", forcefield,
        "-water", water_model,
    ]
    if ignh:
        cmd.append("-ignh")
    return cmd


def run_pdb2gmx(
    protein_pdb: str | Path,
    out_dir: str | Path,
    forcefield: str = "oplsaa",
    water_model: str = "spce",
) -> tuple[Path, Path, Path, Pdb2gmxMeta]:
    """
    Run ``gmx pdb2gmx`` on a protein PDB with adaptive -ignh fallback.

    Attempt 1: run without ``-ignh``.
    If it fails and GROMACS output contains a hydrogen/protonation-mismatch
    diagnosis (see ``_should_retry_with_ignh``), attempt 2 is run with
    ``-ignh``.  All other failures are propagated immediately.

    Returns:
        (protein_gro, topol_top, posre_itp, Pdb2gmxMeta)

    Raises:
        RuntimeError if gmx is absent or both attempts fail.
    """
    if not shutil.which("gmx"):
        raise RuntimeError(
            "GROMACS (gmx) not found in PATH. "
            "Install GROMACS or activate its environment before running integrate."
        )

    protein_pdb = Path(protein_pdb)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    protein_gro = out_dir / "protein.gro"
    topol_top = out_dir / "topol.top"
    posre_itp = out_dir / "posre.itp"

    # ── Attempt 1: without -ignh ──────────────────────────────────────────────
    cmd1 = _build_pdb2gmx_cmd(
        protein_pdb, protein_gro, topol_top, posre_itp,
        forcefield, water_model, ignh=False,
    )
    r1 = subprocess.run(cmd1, capture_output=True, text=True, cwd=out_dir)
    combined1 = r1.stdout + "\n" + r1.stderr

    meta = Pdb2gmxMeta(
        initial_exit_code=r1.returncode,
        initial_command=cmd1,
        initial_error_summary=r1.stderr.strip()[:500] if r1.returncode != 0 else "",
        final_exit_code=r1.returncode,
    )

    if r1.returncode == 0:
        return protein_gro, topol_top, posre_itp, meta

    # ── Attempt 1 failed — check for H-mismatch ───────────────────────────────
    should_retry, reason = _should_retry_with_ignh(combined1)

    if not should_retry:
        raise RuntimeError(
            f"gmx pdb2gmx failed (exit {r1.returncode}):\n"
            f"{r1.stderr.strip()}"
        )

    # ── Attempt 2: with -ignh ─────────────────────────────────────────────────
    cmd2 = _build_pdb2gmx_cmd(
        protein_pdb, protein_gro, topol_top, posre_itp,
        forcefield, water_model, ignh=True,
    )
    r2 = subprocess.run(cmd2, capture_output=True, text=True, cwd=out_dir)

    meta.fallback_used = True
    meta.fallback_reason = reason
    meta.fallback_command = cmd2
    meta.final_exit_code = r2.returncode

    if r2.returncode != 0:
        raise RuntimeError(
            f"gmx pdb2gmx failed even with -ignh (exit {r2.returncode}):\n"
            f"{r2.stderr.strip()}"
        )

    return protein_gro, topol_top, posre_itp, meta


# ── Top-level pipeline ────────────────────────────────────────────────────────

def _referenced_posre_includes(itp_text: str) -> list[str]:
    """Return posre_*.itp filenames referenced by #include in itp_text."""
    found = []
    for line in itp_text.splitlines():
        s = line.strip()
        if not s.startswith("#include"):
            continue
        m = _INCLUDE_RE.match(s)
        if m and "posre" in m.group(1).lower():
            found.append(Path(m.group(1)).name)
    return found


def assemble_system_multi(
    *,
    protein_pdb: Optional[str | Path] = None,
    protein_gro: Optional[str | Path] = None,
    protein_topology_itps: Optional[list[str | Path]] = None,
    protein_posre_itps: Optional[list[str | Path]] = None,
    components: list[ParameterizedMolecule],
    out_dir: str | Path,
    forcefield: str = "oplsaa",
    water_model: str = "spce",
    protein_itp_name: str = "protein.itp",
    run_grompp: bool = True,
    atomtypes_filename: str = "component_atomtypes.itp",
) -> AssemblyResult:
    """
    Full protein + N-parameterized-non-protein-component GROMACS assembly.

    Generalizes ``assemble_system`` to any number of independently
    LigParGen-parameterized molecules (ligand, cofactor, substrate,
    inhibitor, ...) coexisting in the same system. ``assemble_system`` is a
    single-component wrapper around this function.

    Two mutually exclusive protein input modes:

      MODE A — raw protein PDB (pass ``protein_pdb``):
        gmx pdb2gmx runs on it, exactly as before. Supports embedded topology
        mode (inline moleculetype) and master topology injection mode
        (per-chain topol_*.itp files) depending on what pdb2gmx produces.

      MODE B — pre-parameterized protein (pass ``protein_gro`` +
      ``protein_topology_itps``, optionally ``protein_posre_itps``):
        pdb2gmx is never invoked. The caller supplies protein coordinates
        that have *already* been through pdb2gmx (or an equivalent
        parameterization) plus its chain topology .itp file(s) directly.
        This exists because re-deriving a raw PDB from an already-processed
        GRO (e.g. by naively converting coordinates back to PDB) does not
        reliably round-trip through pdb2gmx: post-processing can rename
        atoms in ways the force field's .rtp templates no longer recognize
        (pdb2gmx then fails with e.g. "Atom O1 in residue ASN ... was not
        found in rtp entry ASN"). A pre-parameterized protein is a
        first-class input, not a workaround — do not attempt to sanitize or
        re-run such a protein through pdb2gmx.

    Each protein topology file in ``protein_topology_itps`` must contain
    exactly one [ moleculetype ] (one file per chain, as pdb2gmx itself
    writes them) — its molecule name and atom count are parsed from the
    file, not assumed from its filename. Each file's own position-restraint
    #include (already embedded when pdb2gmx originally wrote it) is
    preserved verbatim; no new restraint include is generated for it.

    Each component's coordinates are copied verbatim into complex.gro in the
    order given (protein, then each component) — this function assumes the
    caller has already resolved system-specific pose vs. parameterization
    geometry (see ``ligand.pose_rewriter``) and passes the correct .gro for
    each component; no relocation or transformation is applied here.

    Cross-component [ atomtypes ] are merged via
    ``resolve_and_merge_atomtypes`` — identical definitions are deduplicated,
    conflicting same-name definitions are namespaced per component, and an
    unresolvable conflict raises before anything is written past that point.
    This works identically regardless of protein input mode.

    All intermediate and final files are written to ``out_dir``.
    GROMACS must be in PATH for Mode A (required for pdb2gmx); Mode B never
    calls GROMACS except for the optional grompp dry-run validation.
    """
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    preparameterized = protein_gro is not None or bool(protein_topology_itps)
    if preparameterized and protein_pdb is not None:
        return AssemblyResult(
            success=False,
            error="protein_pdb (raw PDB, runs pdb2gmx) and protein_gro/"
            "protein_topology_itps (pre-parameterized, skips pdb2gmx) are "
            "mutually exclusive.",
        )
    if not preparameterized and protein_pdb is None:
        return AssemblyResult(
            success=False,
            error="Provide either protein_pdb (raw PDB) or protein_gro + "
            "protein_topology_itps (pre-parameterized protein).",
        )
    if preparameterized and (protein_gro is None or not protein_topology_itps):
        return AssemblyResult(
            success=False,
            error="Pre-parameterized protein input requires both protein_gro "
            "and at least one entry in protein_topology_itps.",
        )

    if not components:
        return AssemblyResult(
            success=False,
            error="No components provided; at least one parameterized "
            "non-protein molecule is required.",
        )

    components = [
        ParameterizedMolecule(
            component_id=c.component_id,
            topology_path=Path(c.topology_path).resolve(),
            coordinate_path=Path(c.coordinate_path).resolve(),
            role=c.role,
        )
        for c in components
    ]

    warnings: list[str] = []
    errors: list[str] = []

    # 0. Pre-flight: each component's ITP atom count must match its GRO. ──────
    component_itps: dict[str, ItpFile] = {}
    component_gros: dict[str, GroFile] = {}
    for c in components:
        try:
            itp = parse_itp(c.topology_path)
        except Exception as exc:
            return AssemblyResult(
                success=False,
                error=f"[{c.component_id}] Cannot parse topology {c.topology_path}: {exc}",
            )
        try:
            gro = parse_gro(c.coordinate_path)
        except Exception as exc:
            return AssemblyResult(
                success=False,
                error=f"[{c.component_id}] Cannot parse coordinates {c.coordinate_path}: {exc}",
            )
        if len(itp.atoms) != gro.atom_count:
            return AssemblyResult(
                success=False,
                error=(
                    f"[{c.component_id}] Atom count mismatch: "
                    f"{c.topology_path.name} declares {len(itp.atoms)} atoms "
                    f"but {c.coordinate_path.name} has {gro.atom_count}."
                ),
            )
        component_itps[c.component_id] = itp
        component_gros[c.component_id] = gro

    # 1-2. Protein input (mode-dependent) ─────────────────────────────────────
    pdb2gmx_meta: Optional[Pdb2gmxMeta] = None

    if preparameterized:
        # MODE B — pre-parameterized protein; pdb2gmx is never invoked.
        protein_gro_src = Path(protein_gro).resolve()
        protein_topology_srcs = [Path(p).resolve() for p in protein_topology_itps]
        protein_posre_srcs = [Path(p).resolve() for p in (protein_posre_itps or [])]

        protein_molecules_list: list[tuple[str, int]] = []
        protein_include_names: list[str] = []
        total_protein_itp_atoms = 0
        for itp_path in protein_topology_srcs:
            try:
                p_itp = parse_itp(itp_path)
            except Exception as exc:
                return AssemblyResult(
                    success=False,
                    error=f"Cannot parse protein topology {itp_path}: {exc}",
                )
            if not p_itp.moleculetype:
                return AssemblyResult(
                    success=False,
                    error=f"No [ moleculetype ] found in protein topology {itp_path}.",
                )
            protein_molecules_list.append((p_itp.moleculetype.name, 1))
            protein_include_names.append(itp_path.name)
            total_protein_itp_atoms += len(p_itp.atoms)

        try:
            protein_gro_parsed = parse_gro(protein_gro_src)
        except Exception as exc:
            return AssemblyResult(
                success=False, error=f"Cannot parse protein GRO {protein_gro_src}: {exc}"
            )

        if total_protein_itp_atoms != protein_gro_parsed.atom_count:
            return AssemblyResult(
                success=False,
                error=(
                    f"Pre-parameterized protein atom count mismatch: protein "
                    f"topology declares {total_protein_itp_atoms} atoms across "
                    f"{len(protein_topology_srcs)} file(s) but "
                    f"{protein_gro_src.name} has {protein_gro_parsed.atom_count}."
                ),
            )

        # Copy protein GRO, chain topology files, and posre files verbatim.
        protein_gro = out_dir / "protein.gro"
        protein_gro.write_text(protein_gro_src.read_text())
        for itp_path in protein_topology_srcs:
            (out_dir / itp_path.name).write_text(itp_path.read_text())
        for posre_path in protein_posre_srcs:
            (out_dir / posre_path.name).write_text(posre_path.read_text())

        # Position restraints are already wired via #include inside each
        # chain topology file — verify what's referenced actually made it
        # into out_dir, warn (don't block) if it did not, since a missing
        # posre file only matters to a later -DPOSRES equilibration run, not
        # to this assembly step.
        provided_posre_names = {p.name for p in protein_posre_srcs}
        for itp_path in protein_topology_srcs:
            for ref in _referenced_posre_includes(itp_path.read_text()):
                if ref not in provided_posre_names and not (out_dir / ref).exists():
                    warnings.append(
                        f"{itp_path.name} references position-restraint include "
                        f"'{ref}' which was not supplied via protein_posre_itps "
                        f"and is not present in {out_dir}. Assembly will still "
                        f"succeed, but a later -DPOSRES run will fail until it "
                        f"is added."
                    )

        chain_itp_includes = protein_include_names
        topology_mode = "preparameterized"
        raw_top: Optional[Path] = None

    else:
        # MODE A — raw protein PDB; run pdb2gmx as before.
        protein_pdb = Path(protein_pdb).resolve()
        try:
            protein_gro, raw_top, _, pdb2gmx_meta = run_pdb2gmx(
                protein_pdb, out_dir, forcefield=forcefield, water_model=water_model
            )
        except RuntimeError as exc:
            return AssemblyResult(success=False, error=str(exc))

        chain_itp_includes = detect_chain_itp_includes(raw_top)
        topology_mode = "master_injection" if chain_itp_includes else "embedded"

    # 3. Extract each component's [ atomtypes ] and resolve namespace conflicts ─
    extracted: dict[str, tuple[str, str]] = {}  # id -> (atomtypes_block, remaining_text)
    for c in components:
        raw_text = c.topology_path.read_text()
        block, remaining = extract_atomtypes_section(raw_text)
        if not block:
            warnings.append(
                f"[{c.component_id}] No [ atomtypes ] section found in "
                f"{c.topology_path.name}."
            )
        extracted[c.component_id] = (block, remaining)

    try:
        merge_result = resolve_and_merge_atomtypes(
            [(c.component_id, extracted[c.component_id][0]) for c in components]
        )
    except ValueError as exc:
        return AssemblyResult(
            success=False, error=f"Atomtype conflict resolution failed: {exc}"
        )

    atomtypes_path = out_dir / atomtypes_filename
    atomtypes_path.write_text(merge_result.combined_block or "; No [ atomtypes ] found\n")

    # 4. Write cleaned + renamed component .itp files ─────────────────────────
    component_molecule_names: dict[str, str] = {}
    for c in components:
        _, remaining = extracted[c.component_id]
        renamed_text = _apply_atomtype_renames(
            remaining, merge_result.renames.get(c.component_id, {})
        )
        cleaned_path = out_dir / c.topology_path.name
        cleaned_path.write_text(renamed_text)

        mtype = component_itps[c.component_id].moleculetype
        component_molecule_names[c.component_id] = mtype.name if mtype else c.component_id

    # 5. Assemble final topol.top (mode-dependent) ────────────────────────────
    protein_itp_path: Optional[Path] = None
    component_specs = [
        (c.topology_path.name, component_molecule_names[c.component_id], 1)
        for c in components
    ]

    if topology_mode == "master_injection":
        protein_molecules_list = parse_molecules_section(raw_top)
        try:
            inject_components_into_master_topol(
                raw_top,
                atomtypes_itp_name=atomtypes_filename,
                components=component_specs,
            )
        except ValueError as exc:
            return AssemblyResult(success=False, error=str(exc))
        topol_top_path = raw_top

    elif topology_mode == "preparameterized":
        topol_top_path = generate_topol_top_preparameterized(
            out_dir=out_dir,
            forcefield=forcefield,
            water_model=water_model,
            protein_topology_names=protein_include_names,
            atomtypes_itp_name=atomtypes_filename,
            protein_molecules=protein_molecules_list,
            components=component_specs,
        )

    else:  # embedded
        try:
            protein_molecules_list, water_ions_block = split_topol_top(
                raw_top, out_dir, protein_itp_name
            )
        except ValueError as exc:
            return AssemblyResult(success=False, error=str(exc))

        protein_itp_path = out_dir / protein_itp_name

        topol_top_path = generate_topol_top_multi(
            out_dir=out_dir,
            forcefield=forcefield,
            water_model=water_model,
            protein_itp_name=protein_itp_name,
            atomtypes_itp_name=atomtypes_filename,
            protein_molecules=protein_molecules_list,
            components=component_specs,
            water_ions_block=water_ions_block,
        )

    protein_mol_name = protein_molecules_list[0][0] if protein_molecules_list else "Protein"

    # 6. Merge GRO files (protein + each component, in order; no relocation) ──
    complex_gro_path = out_dir / "complex.gro"
    try:
        total_atoms = merge_gro_files_multi(
            protein_gro, [c.coordinate_path for c in components], complex_gro_path
        )
    except Exception as exc:
        return AssemblyResult(success=False, error=f"Failed to merge .gro files: {exc}")

    prot_gro = parse_gro(protein_gro)

    # 7. Validate ─────────────────────────────────────────────────────────────
    component_validation = [
        {
            "component_id": c.component_id,
            "molecule_name": component_molecule_names[c.component_id],
            "itp_path": out_dir / c.topology_path.name,
        }
        for c in components
    ]
    validation_errs = validate_system_multi(
        out_dir=out_dir,
        protein_mol_name=protein_mol_name,
        components=component_validation,
        complex_gro_path=complex_gro_path,
        topol_top_path=topol_top_path,
        atomtypes_itp_path=atomtypes_path,
        run_grompp=run_grompp,
    )
    # grompp warnings are surfaced, not swallowed, but do not fail assembly.
    _gw_prefix = "[grompp-warning]"
    warnings.extend(e[len(_gw_prefix):].strip() for e in validation_errs
                    if e.startswith(_gw_prefix))
    errors.extend(e for e in validation_errs if not e.startswith(_gw_prefix))

    # 8. Report ───────────────────────────────────────────────────────────────
    _meta = pdb2gmx_meta
    prot_mols_dicts = [{"name": n, "count": c} for n, c in protein_molecules_list]

    components_report: list[dict] = []
    for c in components:
        gro = component_gros[c.component_id]
        components_report.append({
            "id": c.component_id,
            "molecule_name": component_molecule_names[c.component_id],
            "atom_count": gro.atom_count,
            "count": 1,
            "role": c.role,
            "topology": c.topology_path.name,
            "coordinates": c.coordinate_path.name,
            "atomtype_renames": dict(merge_result.renames.get(c.component_id, {})),
        })

    report_data = {
        "success": not errors,
        "topology_mode": topology_mode,
        "protein_input_mode": "preparameterized" if preparameterized else "raw_pdb",
        "pdb2gmx_used": not preparameterized,
        "protein_topology_includes": chain_itp_includes,
        "protein": {"molecules": prot_mols_dicts},
        "protein_molecules": prot_mols_dicts,
        "components": components_report,
        "protein_mol_name": protein_mol_name,
        "protein_atom_count": prot_gro.atom_count,
        "total_atom_count": total_atoms,
        "pdb2gmx_initial_attempt_exit_code": _meta.initial_exit_code if _meta else None,
        "pdb2gmx_initial_command": " ".join(_meta.initial_command) if _meta else None,
        "pdb2gmx_initial_error_summary": _meta.initial_error_summary if _meta else None,
        "pdb2gmx_fallback_used": _meta.fallback_used if _meta else False,
        "pdb2gmx_fallback_reason": _meta.fallback_reason if _meta else "",
        "pdb2gmx_fallback_command": " ".join(_meta.fallback_command) if (_meta and _meta.fallback_command) else None,
        "pdb2gmx_final_exit_code": _meta.final_exit_code if _meta else None,
        "outputs": {
            "protein_gro": str(protein_gro),
            "protein_itp": str(protein_itp_path) if protein_itp_path else None,
            "component_atomtypes_itp": str(atomtypes_path),
            "topol_top": str(topol_top_path),
            "complex_gro": str(complex_gro_path),
        },
        "warnings": warnings,
        "errors": errors,
    }

    # Back-compat report fields for the single-component (legacy ligand) case.
    if len(components) == 1:
        only = components[0]
        report_data["ligand_molecule"] = {
            "name": component_molecule_names[only.component_id], "count": 1,
        }
        report_data["ligand_mol_name"] = component_molecule_names[only.component_id]
        report_data["ligand_atom_count"] = component_gros[only.component_id].atom_count
        report_data["outputs"]["ligand_atomtypes_itp"] = str(atomtypes_path)
        report_data["outputs"]["ligand_itp"] = str(out_dir / only.topology_path.name)

    report_path = out_dir / "assembly_report.yaml"
    report_path.write_text(
        yaml.dump(report_data, default_flow_style=False, allow_unicode=True)
    )

    fallback_used = _meta.fallback_used if _meta else False
    result = AssemblyResult(
        success=not errors,
        protein_gro=protein_gro,
        protein_itp=protein_itp_path,
        topol_top=topol_top_path,
        complex_gro=complex_gro_path,
        report_path=report_path,
        topology_mode=topology_mode,
        protein_topology_includes=chain_itp_includes,
        protein_molecules=prot_mols_dicts,
        protein_mol_name=protein_mol_name,
        protein_atom_count=prot_gro.atom_count,
        total_atom_count=total_atoms,
        pdb2gmx_fallback_used=fallback_used,
        pdb2gmx_meta=_meta,
        warnings=warnings,
        errors=errors,
        components=components_report,
    )

    # Back-compat AssemblyResult fields for the single-component case.
    if len(components) == 1:
        only = components[0]
        result.ligand_atomtypes_itp = atomtypes_path
        result.ligand_itp = out_dir / only.topology_path.name
        result.ligand_mol_name = component_molecule_names[only.component_id]
        result.ligand_atom_count = component_gros[only.component_id].atom_count

    return result


def assemble_system(
    protein_pdb: str | Path,
    ligand_itp: str | Path,
    ligand_gro: str | Path,
    out_dir: str | Path,
    forcefield: str = "oplsaa",
    water_model: str = "spce",
    protein_itp_name: str = "protein.itp",
    run_grompp: bool = True,
) -> AssemblyResult:
    """
    Full protein–ligand GROMACS system assembly.

    Thin single-component wrapper around ``assemble_system_multi`` — kept for
    backward compatibility with the original ``simforge ligand integrate``
    API. Supports embedded topology mode (inline moleculetype) and master
    topology injection mode (per-chain topol_*.itp files). All paths are
    resolved to absolute before subprocess invocation so relative paths work
    correctly regardless of the current working directory.

    All intermediate and final files are written to ``out_dir``.
    GROMACS must be in PATH (required for pdb2gmx in step 1).
    """
    ligand_itp = Path(ligand_itp)
    ligand_gro = Path(ligand_gro)
    component = ParameterizedMolecule(
        component_id=ligand_itp.stem,
        topology_path=ligand_itp,
        coordinate_path=ligand_gro,
    )
    return assemble_system_multi(
        protein_pdb=protein_pdb,
        components=[component],
        out_dir=out_dir,
        forcefield=forcefield,
        water_model=water_model,
        protein_itp_name=protein_itp_name,
        run_grompp=run_grompp,
        atomtypes_filename="ligand_atomtypes.itp",
    )


# ── Internal helpers ──────────────────────────────────────────────────────────

def _mol_name_from_itp(itp_path: str | Path) -> str:
    """Return the molecule name from [ moleculetype ] or 'LIG' as fallback."""
    try:
        itp = parse_itp(itp_path)
        if itp.moleculetype:
            return itp.moleculetype.name
    except Exception:
        pass
    return "LIG"
