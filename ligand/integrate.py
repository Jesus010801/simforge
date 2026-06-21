"""
Protein–ligand GROMACS system assembly from LigParGen outputs.

Given a prepared protein PDB and LigParGen-generated .itp/.gro files, this
module builds a complete GROMACS-ready system without manual topology editing.

Steps performed by assemble_system():
  1. gmx pdb2gmx on protein → protein.gro, raw topol.top, posre.itp
     1a. First attempt without -ignh.
     1b. On failure, inspect output: if GROMACS diagnoses a hydrogen/protonation
         mismatch (e.g. "Option -ignh will ignore all hydrogens in the input"),
         retry with -ignh. Any other failure is fatal immediately.
  2. Split topol.top → protein.itp (protein topology body extracted)
  3. Extract [ atomtypes ] block from ligand .itp → ligand_atomtypes.itp
  4. Write cleaned ligand .itp (atomtypes removed) to out_dir
  5. Generate master topol.top with correct include order
  6. Merge protein.gro + ligand.gro → complex.gro
  7. Validate the assembled system
  8. Write assembly_report.yaml

Public API:
  Pdb2gmxMeta                              — dataclass with both-attempt metadata
  extract_atomtypes_section(itp_text)      -> tuple[str, str]
  remove_atomtypes_from_itp(itp_path, out_path) -> None
  split_topol_top(...)                     -> tuple[str, str]
  generate_topol_top(...)                  -> Path
  merge_gro_files(...)                     -> int
  validate_system(...)                     -> list[str]
  run_pdb2gmx(protein_pdb, out_dir, ...)   -> tuple[Path, Path, Path, Pdb2gmxMeta]
  assemble_system(...)                     -> AssemblyResult
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import yaml

from utils.gro_parser import GroAtom, GroFile, parse_gro, write_gro
from utils.itp_parser import parse_itp


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


# ── Result type ───────────────────────────────────────────────────────────────

@dataclass
class AssemblyResult:
    success: bool
    protein_gro: Optional[Path] = None
    protein_itp: Optional[Path] = None
    ligand_atomtypes_itp: Optional[Path] = None
    ligand_itp: Optional[Path] = None
    topol_top: Optional[Path] = None
    complex_gro: Optional[Path] = None
    report_path: Optional[Path] = None
    protein_mol_name: str = ""
    ligand_mol_name: str = ""
    protein_atom_count: int = 0
    ligand_atom_count: int = 0
    total_atom_count: int = 0
    pdb2gmx_fallback_used: bool = False
    pdb2gmx_meta: Optional[Pdb2gmxMeta] = None
    warnings: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    error: Optional[str] = None


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


# ── topol.top splitting ───────────────────────────────────────────────────────

def split_topol_top(
    topol_path: str | Path,
    out_dir: str | Path,
    protein_itp_name: str = "protein.itp",
) -> tuple[str, str]:
    """
    Read a pdb2gmx-generated topol.top and extract the protein topology.

    Writes ``protein_itp_name`` to ``out_dir`` containing all [ section ]
    blocks that belong to the protein (from [ moleculetype ] up to, but not
    including, the first water/ions .ff/ include).

    Returns:
        (protein_mol_name, water_ions_block)
        - protein_mol_name: molecule name from [ moleculetype ]
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

    # Protein topology body: from [ moleculetype ] to just before water include
    protein_lines = lines[first_mol_idx:water_include_idx]
    protein_mol_name = _mol_name_from_lines(protein_lines)

    # Water+ions block: from water include to just before [ system ]
    water_end = system_idx if system_idx is not None else len(lines)
    water_ions_block = "".join(lines[water_include_idx:water_end]).rstrip()

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / protein_itp_name).write_text("".join(protein_lines))

    return protein_mol_name, water_ions_block


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
    protein_mol_name: str,
    ligand_mol_name: str,
    water_ions_block: str = "",
) -> Path:
    """
    Write a clean master topol.top with the correct GROMACS include order.

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

    # Map common shorthand to GROMACS ff directory names
    _FF_DIR = {
        "oplsaa": "oplsaa.ff",
        "opls-aa": "oplsaa.ff",
        "charmm36": "charmm36-jul2022.ff",
        "amber99sb": "amber99sb.ff",
        "amber99sb-ildn": "amber99sb-ildn.ff",
        "gromos54a7": "gromos54a7.ff",
    }
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
        f"{protein_mol_name}           1",
        f"{ligand_mol_name}            1",
    ]

    top_path = out_dir / "topol.top"
    top_path.write_text("\n".join(sections) + "\n")
    return top_path


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
    prot = parse_gro(protein_gro)
    lig = parse_gro(ligand_gro)

    merged: list[GroAtom] = list(prot.atoms) + list(lig.atoms)

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


def _grompp_check(out_dir: Path, gro: Path, top: Path) -> list[str]:
    """Run gmx grompp with a minimal MDP as a topology syntax check."""
    mdp = out_dir / "_gmx_check.mdp"
    tpr = out_dir / "_gmx_check.tpr"
    mdp.write_text(
        "integrator   = steep\n"
        "nsteps       = 0\n"
        "emtol        = 100\n"
        "nstlist      = 1\n"
        "cutoff-scheme = Verlet\n"
    )
    try:
        result = subprocess.run(
            [
                "gmx", "grompp",
                "-f", str(mdp),
                "-c", str(gro),
                "-p", str(top),
                "-o", str(tpr),
                "-maxwarn", "5",
            ],
            capture_output=True, text=True, cwd=out_dir, timeout=30,
        )
        if result.returncode != 0:
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

    All intermediate and final files are written to ``out_dir``.
    GROMACS must be in PATH (required for pdb2gmx in step 1).
    """
    out_dir = Path(out_dir)
    ligand_itp = Path(ligand_itp)
    ligand_gro = Path(ligand_gro)
    out_dir.mkdir(parents=True, exist_ok=True)

    warnings: list[str] = []
    errors: list[str] = []

    # 1. pdb2gmx ──────────────────────────────────────────────────────────────
    pdb2gmx_meta: Optional[Pdb2gmxMeta] = None
    try:
        protein_gro, raw_top, _, pdb2gmx_meta = run_pdb2gmx(
            protein_pdb, out_dir, forcefield=forcefield, water_model=water_model
        )
    except RuntimeError as exc:
        return AssemblyResult(success=False, error=str(exc))

    # 2. Split topol.top → protein.itp ────────────────────────────────────────
    try:
        protein_mol_name, water_ions_block = split_topol_top(
            raw_top, out_dir, protein_itp_name
        )
    except ValueError as exc:
        return AssemblyResult(success=False, error=str(exc))

    protein_itp_path = out_dir / protein_itp_name

    # 3-4. Extract ligand [ atomtypes ] ───────────────────────────────────────
    itp_text = ligand_itp.read_text()
    atomtypes_block, _ = extract_atomtypes_section(itp_text)

    if not atomtypes_block:
        warnings.append(
            f"No [ atomtypes ] section found in {ligand_itp.name}. "
            "ligand_atomtypes.itp will be empty — "
            "ensure atomtypes are defined elsewhere."
        )

    atomtypes_path = out_dir / "ligand_atomtypes.itp"
    atomtypes_path.write_text(atomtypes_block or "; No [ atomtypes ] found\n")

    # 5. Write cleaned ligand .itp ─────────────────────────────────────────────
    ligand_mol_name = _mol_name_from_itp(ligand_itp)
    ligand_itp_out = out_dir / ligand_itp.name
    remove_atomtypes_from_itp(ligand_itp, ligand_itp_out)

    # 6. Generate master topol.top ────────────────────────────────────────────
    topol_top_path = generate_topol_top(
        out_dir=out_dir,
        forcefield=forcefield,
        water_model=water_model,
        protein_itp_name=protein_itp_name,
        ligand_itp_name=ligand_itp.name,
        ligand_atomtypes_itp_name="ligand_atomtypes.itp",
        protein_mol_name=protein_mol_name,
        ligand_mol_name=ligand_mol_name,
        water_ions_block=water_ions_block,
    )

    # 7. Merge GRO files ──────────────────────────────────────────────────────
    complex_gro_path = out_dir / "complex.gro"
    try:
        total_atoms = merge_gro_files(protein_gro, ligand_gro, complex_gro_path)
    except Exception as exc:
        return AssemblyResult(
            success=False, error=f"Failed to merge .gro files: {exc}"
        )

    prot_gro = parse_gro(protein_gro)
    lig_gro_parsed = parse_gro(ligand_gro)

    # 8. Validate ─────────────────────────────────────────────────────────────
    validation_errs = validate_system(
        out_dir=out_dir,
        protein_mol_name=protein_mol_name,
        ligand_mol_name=ligand_mol_name,
        ligand_itp_path=ligand_itp_out,
        complex_gro_path=complex_gro_path,
        topol_top_path=topol_top_path,
        ligand_atomtypes_itp_path=atomtypes_path,
        run_grompp=run_grompp,
    )
    errors.extend(validation_errs)

    # 9. Report ───────────────────────────────────────────────────────────────
    _meta = pdb2gmx_meta  # may be None if pdb2gmx was never reached (shouldn't happen here)
    report_data = {
        "success": not errors,
        "protein_mol_name": protein_mol_name,
        "ligand_mol_name": ligand_mol_name,
        "protein_atom_count": prot_gro.atom_count,
        "ligand_atom_count": lig_gro_parsed.atom_count,
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
            "protein_itp": str(protein_itp_path),
            "ligand_atomtypes_itp": str(atomtypes_path),
            "ligand_itp": str(ligand_itp_out),
            "topol_top": str(topol_top_path),
            "complex_gro": str(complex_gro_path),
        },
        "warnings": warnings,
        "errors": errors,
    }
    report_path = out_dir / "assembly_report.yaml"
    report_path.write_text(
        yaml.dump(report_data, default_flow_style=False, allow_unicode=True)
    )

    fallback_used = _meta.fallback_used if _meta else False
    return AssemblyResult(
        success=not errors,
        protein_gro=protein_gro,
        protein_itp=protein_itp_path,
        ligand_atomtypes_itp=atomtypes_path,
        ligand_itp=ligand_itp_out,
        topol_top=topol_top_path,
        complex_gro=complex_gro_path,
        report_path=report_path,
        protein_mol_name=protein_mol_name,
        ligand_mol_name=ligand_mol_name,
        protein_atom_count=prot_gro.atom_count,
        ligand_atom_count=lig_gro_parsed.atom_count,
        total_atom_count=total_atoms,
        pdb2gmx_fallback_used=fallback_used,
        pdb2gmx_meta=_meta,
        warnings=warnings,
        errors=errors,
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
