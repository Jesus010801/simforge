"""
Ligand campaign layer — discovery, pose extraction, ligand identity/dedup,
and manifest generation for multi-complex, same-ligand workflows.

Built entirely on top of the existing single-complex ligand pipeline
(``ligand.prepare``, ``ligand.hydrogenation``). This module does not
duplicate hydrogenation, charge estimation, or PDB-splitting logic — it
reuses ``ligand.prepare.extract_ligand_from_complex`` for chemical
preparation and adds a thinner, non-chemical extraction step
(``extract_complex_components``) for receptor-specific pose isolation.

Pipeline (see docs referenced in the campaign CLI help):

  1. discover_complexes()            — find one complex PDB per system dir
  2. extract_complex_components()    — split into protein_only.pdb / ligand_pose.pdb
                                        (no hydrogenation, no optimization)
  3. compute_ligand_signature() /
     compute_chemical_identity()     — determine ligand chemical identity
  4. group_ligand_occurrences()      — deduplicate ligand occurrences
  5. prepare_campaign()              — orchestrates 1-4, hydrogenates each
                                        unique ligand exactly once, writes
                                        the campaign manifest + report

Public API:
  ComplexEntry, CampaignDiscoveryError, discover_complexes
  ExtractedComplex, extract_complex_components
  LigandSignature, compute_ligand_signature
  ChemicalIdentity, compute_chemical_identity
  LigandGroup, group_ligand_occurrences
  SystemManifestEntry, LigandManifestEntry, CampaignManifest
  CampaignPrepareResult, prepare_campaign
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from pydantic import BaseModel, Field

from ligand.prepare import (
    _conect_first_serial,
    _count_h_atoms,
    _count_heavy_atoms,
    _parse_serials,
    detect_ligand_residues,
    extract_ligand_from_complex,
)

_CAMPAIGN_DIR_NAME = "simforge_campaign"
_SYSTEM_SUBDIR_NAME = "simforge"
_IGNORED_TOP_DIR_NAMES = frozenset({_CAMPAIGN_DIR_NAME, _SYSTEM_SUBDIR_NAME})


# ═══════════════════════════════════════════════════════════════════════════
# Phase 1 — discovery
# ═══════════════════════════════════════════════════════════════════════════

class CampaignDiscoveryError(Exception):
    """Raised when campaign discovery cannot resolve a system unambiguously."""


@dataclass(frozen=True)
class ComplexEntry:
    system_id: str
    directory: Path
    complex_pdb: Path


def discover_complexes(root: str | Path) -> list[ComplexEntry]:
    """
    Recursively discover one candidate complex PDB per immediate subdirectory
    of ``root``. The subdirectory name becomes the system id.

    SimForge-generated output directories (``simforge/``, ``simforge_campaign/``)
    are ignored. A directory containing more than one candidate PDB is
    ambiguous and raises ``CampaignDiscoveryError`` rather than guessing.

    Returns a deterministically ordered list (sorted by system_id).
    """
    root = Path(root)
    if not root.exists():
        raise CampaignDiscoveryError(f"Campaign root not found: {root}")
    if not root.is_dir():
        raise CampaignDiscoveryError(f"Campaign root is not a directory: {root}")

    entries: list[ComplexEntry] = []
    ambiguous: list[str] = []

    for directory in sorted(p for p in root.iterdir() if p.is_dir()):
        if directory.name in _IGNORED_TOP_DIR_NAMES or directory.name.startswith("."):
            continue

        candidates = sorted(
            p for p in directory.rglob("*.pdb")
            if _SYSTEM_SUBDIR_NAME not in p.relative_to(directory).parts
        )
        if not candidates:
            continue
        if len(candidates) > 1:
            ambiguous.append(
                f"  {directory.name}/: {', '.join(p.relative_to(directory).as_posix() for p in candidates)}"
            )
            continue

        entries.append(ComplexEntry(
            system_id=directory.name,
            directory=directory,
            complex_pdb=candidates[0],
        ))

    if ambiguous:
        raise CampaignDiscoveryError(
            "Ambiguous complex directories — expected exactly one candidate PDB "
            "per system directory, found multiple:\n" + "\n".join(ambiguous) +
            "\nMove or remove the extra files, or split them into separate directories."
        )

    return sorted(entries, key=lambda e: e.system_id)


# ═══════════════════════════════════════════════════════════════════════════
# Phase 2 — pose extraction (no chemical preparation)
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class ExtractedComplex:
    success: bool
    protein_pdb: Optional[Path] = None
    ligand_pose_pdb: Optional[Path] = None
    ligand_resname: str = ""
    ligand_atom_count: int = 0
    protein_atom_count: int = 0
    n_hydrogen_atoms: int = 0
    n_heavy_atoms: int = 0
    error: Optional[str] = None


def extract_complex_components(
    complex_pdb: str | Path,
    ligand_resname: str,
    out_dir: str | Path,
) -> ExtractedComplex:
    """
    Split a docked complex PDB into ``protein_only.pdb`` and
    ``ligand_pose.pdb`` without any chemical modification.

    Unlike ``ligand.prepare.extract_ligand_from_complex``, this performs no
    hydrogenation, no charge estimation, and no optimization — the ligand
    pose file preserves the exact original heavy-atom (and any existing
    hydrogen) coordinates from the complex. Use this for campaign pose
    extraction, where chemical preparation happens once per unique ligand
    rather than once per complex.
    """
    complex_pdb = Path(complex_pdb)
    out_dir = Path(out_dir)
    ligand_resname = ligand_resname.strip().upper()

    if not complex_pdb.exists():
        return ExtractedComplex(success=False, error=f"Complex PDB not found: {complex_pdb}")

    out_dir.mkdir(parents=True, exist_ok=True)

    protein_lines: list[str] = []
    protein_atom_count = 0
    ligand_lines: list[str] = []
    all_conect: list[str] = []

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
            protein_lines.append(line)
        elif record == "CONECT":
            all_conect.append(line)

    if not ligand_lines:
        candidates = detect_ligand_residues(complex_pdb)
        return ExtractedComplex(
            success=False,
            ligand_resname=ligand_resname,
            error=(
                f"No records found for residue '{ligand_resname}' in {complex_pdb.name}. "
                f"Detected ligand residues: {candidates or ['(none)']}"
            ),
        )

    ligand_serials = _parse_serials(ligand_lines)
    conect_lines = [ln for ln in all_conect if _conect_first_serial(ln) in ligand_serials]

    lig_pdb = out_dir / "ligand_pose.pdb"
    parts = ligand_lines[:] + conect_lines + ["END"]
    lig_pdb.write_text("\n".join(parts) + "\n")

    prot_pdb = out_dir / "protein_only.pdb"
    has_ter = any(l[:6].rstrip() == "TER" for l in protein_lines)
    end_suffix = "\nEND\n" if has_ter else "\nTER\nEND\n"
    prot_pdb.write_text("\n".join(protein_lines) + end_suffix)

    return ExtractedComplex(
        success=True,
        protein_pdb=prot_pdb,
        ligand_pose_pdb=lig_pdb,
        ligand_resname=ligand_resname,
        ligand_atom_count=len(ligand_lines),
        protein_atom_count=protein_atom_count,
        n_hydrogen_atoms=_count_h_atoms(ligand_lines),
        n_heavy_atoms=_count_heavy_atoms(ligand_lines),
    )


# ═══════════════════════════════════════════════════════════════════════════
# Phase 3 — ligand identity
# ═══════════════════════════════════════════════════════════════════════════

@dataclass(frozen=True)
class LigandSignature:
    """
    Structural identity descriptor derived from a PDB ligand pose.

    NOT proof of absolute chemical identity (protonation state, tautomer,
    and bond orders are not recoverable from coordinates alone). Use
    ``ChemicalIdentity`` (from an explicit SDF/MOL reference) for a
    higher-confidence identity source.
    """
    resname: str
    heavy_atom_count: int
    atom_names: tuple[str, ...]
    elements: tuple[str, ...]

    def matches(self, other: "LigandSignature") -> bool:
        """
        Equivalence used for grouping: same residue name, same heavy-atom
        count, same element composition (order-independent). Atom naming is
        intentionally excluded — cosmetic naming differences between
        extraction tools should not fracture an otherwise identical ligand
        into separate groups.
        """
        return (
            self.resname == other.resname
            and self.heavy_atom_count == other.heavy_atom_count
            and tuple(sorted(self.elements)) == tuple(sorted(other.elements))
        )


def compute_ligand_signature(ligand_pose_pdb: str | Path, ligand_resname: str) -> LigandSignature:
    """Compute a structural signature from a ligand pose PDB (heavy atoms only).

    Applies the same context-aware element normalization used before
    chemical preparation (``ligand.chemical_perception.normalize_ligand_elements``)
    so a docking-tool element-column defect (e.g. a "CL"-named atom whose
    element column says "C") does not corrupt identity/dedup grouping --
    without this, such a ligand's signature would show one carbon too many
    and one chlorine too few, causing it to be spuriously excluded from an
    otherwise-matching --ligand-reference identity group.
    """
    from ligand.chemical_perception import normalize_ligand_elements

    path = Path(ligand_pose_pdb)
    normalized_text = normalize_ligand_elements(path.read_text()).pdb_text
    names: list[str] = []
    elements: list[str] = []

    for line in normalized_text.splitlines():
        record = line[:6].rstrip()
        if record not in ("ATOM", "HETATM"):
            continue
        raw_element = line[76:78].strip() if len(line) >= 78 else ""
        atom_name = line[12:16].strip() if len(line) >= 16 else ""
        if raw_element:
            # Proper-case symbol (e.g. "Cl", not "CL"/"cl") -- matches
            # RDKit's canonical convention (ligand.rdkit_reader,
            # ligand.pose_rewriter) so element-list comparisons elsewhere
            # (e.g. against a --ligand-reference's chemical identity) are
            # never defeated by a case mismatch alone.
            element = (
                raw_element[0].upper() + raw_element[1:].lower()
                if len(raw_element) > 1 else raw_element.upper()
            )
        else:
            element = atom_name[0].upper() if atom_name else "X"
        if element.upper() in ("H", "D"):
            continue
        names.append(atom_name)
        elements.append(element)

    order = sorted(range(len(names)), key=lambda i: (elements[i], names[i]))
    return LigandSignature(
        resname=ligand_resname.strip().upper(),
        heavy_atom_count=len(names),
        atom_names=tuple(names[i] for i in order),
        elements=tuple(elements[i] for i in order),
    )


@dataclass(frozen=True)
class ChemicalIdentity:
    """Authoritative chemical identity derived from an explicit reference file."""
    source_path: Path
    canonical_smiles: Optional[str] = None
    inchikey: Optional[str] = None
    heavy_atom_count: int = 0
    elements: tuple[str, ...] = ()


def compute_chemical_identity(reference_path: str | Path) -> ChemicalIdentity:
    """
    Derive a canonical chemical identity (SMILES / InChIKey) from an explicit
    SDF/MOL/PDB reference using RDKit.

    RDKit is imported lazily; callers should catch ``ImportError`` and fall
    back to structural-signature identity with reduced confidence, per
    SimForge's existing optional-RDKit convention (see ``ligand.rdkit_reader``).

    Raises:
        ImportError:  RDKit is not installed.
        FileNotFoundError / ValueError:  reference file missing or unparsable.
    """
    from ligand.rdkit_reader import heavy_atom_elements, load_mol

    path = Path(reference_path)
    mol = load_mol(path)

    from rdkit import Chem

    mol_noh = Chem.RemoveHs(mol)
    smiles = Chem.MolToSmiles(mol_noh, canonical=True) or None

    inchikey: Optional[str] = None
    try:
        inchi = Chem.MolToInchi(mol_noh)
        if inchi:
            inchikey = Chem.InchiToInchiKey(inchi) or None
    except Exception:
        inchikey = None

    return ChemicalIdentity(
        source_path=path,
        canonical_smiles=smiles,
        inchikey=inchikey,
        heavy_atom_count=mol_noh.GetNumAtoms(),
        elements=tuple(sorted(heavy_atom_elements(mol_noh))),
    )


@dataclass
class LigandGroup:
    ligand_id: str
    resname: str
    system_ids: list[str] = field(default_factory=list)
    signature: Optional[LigandSignature] = None
    chemical_identity: Optional[ChemicalIdentity] = None
    identity_method: str = "pdb_atom_signature"
    identity_confidence: str = "structural_signature"
    representative_system_id: str = ""


def group_ligand_occurrences(
    occurrences: list[tuple[str, LigandSignature]],
    ligand_resname: str,
    chemical_identity: Optional[ChemicalIdentity] = None,
) -> tuple[dict[str, LigandGroup], list[str]]:
    """
    Group per-system ligand signatures into reusable-ligand groups.

    Without an explicit chemical reference, occurrences are grouped by
    structural-signature equality. Occurrences sharing a residue name but
    with a divergent signature are never silently merged — they are split
    into a separate group (``<resname>_2``, ``<resname>_3``, ...) with a
    warning explaining why.

    With an explicit chemical reference, all occurrences are assigned to a
    single high-confidence group, but any occurrence whose structural
    signature is inconsistent with the reference (different heavy-atom
    count or element composition) is excluded from the group with a
    warning rather than forced in.

    Returns (groups keyed by ligand_id, warnings).
    """
    warnings: list[str] = []
    groups: dict[str, LigandGroup] = {}
    resname = ligand_resname.strip().upper()

    if chemical_identity is not None:
        group = LigandGroup(
            ligand_id=resname,
            resname=resname,
            identity_method="chemical_reference",
            identity_confidence="high",
            chemical_identity=chemical_identity,
        )
        for system_id, sig in occurrences:
            if (
                chemical_identity.heavy_atom_count
                and sig.heavy_atom_count != chemical_identity.heavy_atom_count
            ):
                warnings.append(
                    f"{system_id}: heavy-atom count {sig.heavy_atom_count} does not match "
                    f"chemical reference ({chemical_identity.heavy_atom_count}); excluded from '{resname}'."
                )
                continue
            if (
                chemical_identity.elements
                and tuple(sorted(sig.elements)) != chemical_identity.elements
            ):
                warnings.append(
                    f"{system_id}: element composition does not match the chemical reference; "
                    f"excluded from '{resname}'."
                )
                continue
            group.system_ids.append(system_id)
            if group.signature is None:
                group.signature = sig

        if not group.system_ids:
            warnings.append(f"No occurrences matched the chemical reference for '{resname}'.")
            return {}, warnings

        group.representative_system_id = group.system_ids[0]
        return {resname: group}, warnings

    # ── Structural-signature grouping ────────────────────────────────────
    known: list[tuple[LigandSignature, str]] = []
    suffix_counters: dict[str, int] = {}

    for system_id, sig in occurrences:
        ligand_id: Optional[str] = None
        for existing_sig, existing_id in known:
            if existing_sig.matches(sig):
                ligand_id = existing_id
                break

        if ligand_id is None:
            base = sig.resname
            if base not in suffix_counters:
                suffix_counters[base] = 1
                ligand_id = base
            else:
                suffix_counters[base] += 1
                ligand_id = f"{base}_{suffix_counters[base]}"
                warnings.append(
                    f"Residue '{base}' has more than one distinct structural signature; "
                    f"'{system_id}' assigned to separate group '{ligand_id}' rather than merged."
                )
            known.append((sig, ligand_id))
            groups[ligand_id] = LigandGroup(
                ligand_id=ligand_id,
                resname=base,
                signature=sig,
                identity_method="pdb_atom_signature",
                identity_confidence="structural_signature",
            )

        groups[ligand_id].system_ids.append(system_id)

    for group in groups.values():
        group.representative_system_id = group.system_ids[0]

    return groups, warnings


# ═══════════════════════════════════════════════════════════════════════════
# Phase 5/6 — manifest
# ═══════════════════════════════════════════════════════════════════════════

class SystemManifestEntry(BaseModel):
    id: str
    directory: str
    complex_pdb: str
    protein: str
    ligand_pose: str
    ligand_id: str
    ligand_resname: str
    protein_atom_count: int = 0
    ligand_pose_atom_count: int = 0


class LigandManifestEntry(BaseModel):
    resname: str
    occurrences: int
    identity_method: str
    identity_confidence: str
    hydrogenation_status: str = "unknown"
    hydrogenation_performed: bool = False
    hydrogenation_backend: str = "none"
    parameterization_status: str = "chemical_validation_required"  # "chemical_validation_required" | "ready_for_ligpargen"
    # Chemical-perception provenance (see ligand.chemical_perception). These
    # are the authoritative signal behind parameterization_status --
    # hydrogenation_performed=True does not by itself mean the chemistry is
    # trustworthy enough to submit to LigParGen.
    element_corrections: list[dict] = Field(default_factory=list)
    element_assignment_confidence: str = "unknown"
    connectivity_confidence: str = "unknown"
    bond_order_confidence: str = "unknown"
    chemical_identity_method: str = "none"
    chemical_identity_confidence: str = "unknown"
    readiness_block_reasons: list[str] = Field(default_factory=list)
    ligand_reference: Optional[str] = None
    ligand_for_ligpargen: Optional[str] = None
    ligand_report: Optional[str] = None
    representative_system_id: str = ""
    warnings: list[str] = Field(default_factory=list)


class CampaignManifest(BaseModel):
    schema_version: int = 1
    campaign_id: str
    root: str
    systems: list[SystemManifestEntry] = Field(default_factory=list)
    ligands: dict[str, LigandManifestEntry] = Field(default_factory=dict)
    warnings: list[str] = Field(default_factory=list)

    def write(self, path: str | Path) -> None:
        Path(path).write_text(self.model_dump_json(indent=2) + "\n")

    @classmethod
    def load(cls, path: str | Path) -> "CampaignManifest":
        return cls.model_validate_json(Path(path).read_text())

    def validate_paths(self, root: str | Path) -> list[str]:
        """Check every referenced relative path still exists under ``root``,
        and that each ligand's parameterization input still reflects current
        element-normalization logic.

        The second check exists because resuming from an existing manifest
        (see ``prepare_campaign``) otherwise trusts on-disk artifacts purely
        by presence: a ``ligand_for_ligpargen`` file written before this
        module's element-normalization fix existed (or by any other process
        that reintroduced an un-normalized element) would be served forever,
        silently, on every subsequent resume -- exactly the same class of
        "confidently wrong" failure this module exists to prevent (see
        ``ligand.chemical_perception`` module docstring). A campaign whose
        stored artifacts disagree with current normalization is therefore
        treated as invalid, the same as a missing file: the caller is asked
        to remove/regenerate it rather than have stale chemistry served
        unnoticed.
        """
        root = Path(root)
        errors: list[str] = []
        for s in self.systems:
            for label, rel in (
                ("complex_pdb", s.complex_pdb),
                ("protein", s.protein),
                ("ligand_pose", s.ligand_pose),
            ):
                if not (root / rel).exists():
                    errors.append(f"System '{s.id}': missing {label} at {rel}")

        from ligand.chemical_perception import normalize_ligand_elements

        for ligand_id, lig in self.ligands.items():
            if not lig.ligand_for_ligpargen:
                continue
            path = root / lig.ligand_for_ligpargen
            if not path.exists():
                errors.append(f"Ligand '{ligand_id}': missing ligand_for_ligpargen at {lig.ligand_for_ligpargen}")
                continue
            stale_corrections = normalize_ligand_elements(path.read_text()).corrections
            if stale_corrections:
                atoms = ", ".join(f"{c.atom_name} ({c.original_element}->{c.corrected_element})" for c in stale_corrections)
                errors.append(
                    f"Ligand '{ligand_id}': {lig.ligand_for_ligpargen} contains un-normalized "
                    f"element(s) not reflected in stored provenance ({atoms}); this campaign was "
                    "prepared before the current element-normalization logic or its artifacts were "
                    "modified externally. Remove simforge_campaign/ and rerun batch-prepare."
                )
        return errors


# ═══════════════════════════════════════════════════════════════════════════
# Phase 4-7 — campaign preparation orchestration
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class CampaignPrepareResult:
    success: bool
    manifest: Optional[CampaignManifest] = None
    manifest_path: Optional[Path] = None
    report_path: Optional[Path] = None
    systems_discovered: int = 0
    unique_ligands: int = 0
    system_errors: dict[str, str] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    error: Optional[str] = None
    resumed: bool = False


def prepare_campaign(
    root: str | Path,
    ligand_resname: str,
    ligand_reference: Optional[str | Path] = None,
    hydrogenation_mode: str = "auto",
) -> CampaignPrepareResult:
    """
    Discover complexes, extract per-system poses, deduplicate ligand
    identity, hydrogenate each unique ligand exactly once, and write
    ``simforge_campaign/campaign_manifest.json`` + ``report.md``.

    Safe to rerun: if a valid manifest already exists it is validated
    against disk and returned unchanged rather than regenerated.
    """
    root = Path(root)
    ligand_resname = ligand_resname.strip().upper()
    campaign_dir = root / _CAMPAIGN_DIR_NAME
    manifest_path = campaign_dir / "campaign_manifest.json"
    report_path = campaign_dir / "report.md"

    # ── Resume: validate rather than regenerate ──────────────────────────
    if manifest_path.exists():
        try:
            manifest = CampaignManifest.load(manifest_path)
        except Exception as exc:
            return CampaignPrepareResult(
                success=False,
                error=(
                    f"Existing campaign manifest is unreadable ({manifest_path}): {exc}. "
                    "Remove or fix it before re-running batch-prepare."
                ),
            )
        path_errors = manifest.validate_paths(root)
        if path_errors:
            return CampaignPrepareResult(
                success=False,
                error=(
                    "Existing campaign manifest references missing files:\n"
                    + "\n".join(f"  - {e}" for e in path_errors)
                    + f"\nRemove {manifest_path} to regenerate, or restore the missing files."
                ),
            )
        return CampaignPrepareResult(
            success=True,
            manifest=manifest,
            manifest_path=manifest_path,
            report_path=report_path if report_path.exists() else None,
            systems_discovered=len(manifest.systems),
            unique_ligands=len(manifest.ligands),
            warnings=manifest.warnings,
            resumed=True,
        )

    # ── Phase 1: discovery ────────────────────────────────────────────────
    try:
        complexes = discover_complexes(root)
    except CampaignDiscoveryError as exc:
        return CampaignPrepareResult(success=False, error=str(exc))

    if not complexes:
        return CampaignPrepareResult(
            success=False,
            error=f"No protein-ligand complex PDBs discovered under {root}.",
        )

    # ── Phase 2: pose extraction ────────────────────────────────────────
    extracted: dict[str, ExtractedComplex] = {}
    signatures: list[tuple[str, LigandSignature]] = []
    system_errors: dict[str, str] = {}

    for entry in complexes:
        out_dir = entry.directory / _SYSTEM_SUBDIR_NAME
        result = extract_complex_components(entry.complex_pdb, ligand_resname, out_dir)
        if not result.success:
            system_errors[entry.system_id] = result.error or "extraction failed"
            continue
        extracted[entry.system_id] = result
        signatures.append((entry.system_id, compute_ligand_signature(result.ligand_pose_pdb, ligand_resname)))

    if not signatures:
        return CampaignPrepareResult(
            success=False,
            error=(
                f"No ligand poses could be extracted for residue '{ligand_resname}'. "
                "Per-system errors:\n"
                + "\n".join(f"  - {sid}: {err}" for sid, err in system_errors.items())
            ),
            system_errors=system_errors,
            systems_discovered=len(complexes),
        )

    # ── Phase 3: chemical identity ──────────────────────────────────────
    warnings: list[str] = []
    chemical_identity: Optional[ChemicalIdentity] = None
    ref_path: Optional[Path] = None
    if ligand_reference is not None:
        ref_path = Path(ligand_reference)
        try:
            chemical_identity = compute_chemical_identity(ref_path)
        except ImportError as exc:
            warnings.append(
                f"--ligand-reference {ref_path} was provided but RDKit is unavailable ({exc}); "
                "falling back to structural-signature identity."
            )
        except Exception as exc:
            warnings.append(
                f"Could not derive chemical identity from --ligand-reference {ref_path} ({exc}); "
                "falling back to structural-signature identity."
            )

    groups, group_warnings = group_ligand_occurrences(signatures, ligand_resname, chemical_identity)
    warnings.extend(group_warnings)

    if not groups:
        return CampaignPrepareResult(
            success=False,
            error="Ligand grouping produced no valid groups.",
            warnings=warnings,
            system_errors=system_errors,
            systems_discovered=len(complexes),
        )

    # ── Phase 4: hydrogenate/prepare each unique ligand once ────────────
    campaign_dir.mkdir(parents=True, exist_ok=True)
    ligands_dir = campaign_dir / "ligands"
    manifest_ligands: dict[str, LigandManifestEntry] = {}

    for ligand_id, group in groups.items():
        lig_dir = ligands_dir / ligand_id
        param_dir = lig_dir / "parameterization"
        param_dir.mkdir(parents=True, exist_ok=True)
        (lig_dir / "ligpargen").mkdir(parents=True, exist_ok=True)

        rep_system = group.representative_system_id
        rep_pose = extracted[rep_system].ligand_pose_pdb

        scratch_dir = param_dir / "_scratch"
        group_has_reference = (
            ref_path is not None and chemical_identity is not None and group.chemical_identity is chemical_identity
        )
        prep_result = extract_ligand_from_complex(
            rep_pose, group.resname, scratch_dir, hydrogenation_mode=hydrogenation_mode,
            chemical_reference=ref_path if group_has_reference else None,
        )

        ligand_for_ligpargen: Optional[Path] = None
        ligand_report: Optional[Path] = None
        if prep_result.success:
            source = prep_result.recommended_ligpargen_input or prep_result.ligand_pdb
            dest = param_dir / "ligand_for_ligpargen.pdb"
            dest.write_text(Path(source).read_text())
            ligand_for_ligpargen = dest

            report_dest = param_dir / "ligand_report.yaml"
            report_dest.write_text(Path(prep_result.report_path).read_text())
            ligand_report = report_dest

            warnings.extend(f"[{ligand_id}] {w}" for w in prep_result.warnings)
        else:
            warnings.append(f"Ligand '{ligand_id}': preparation failed — {prep_result.error}")

        shutil.rmtree(scratch_dir, ignore_errors=True)

        ligand_reference_rel: Optional[str] = None
        if group_has_reference:
            ref_dir = lig_dir / "reference"
            ref_dir.mkdir(parents=True, exist_ok=True)
            ref_dest = ref_dir / f"ligand_reference{ref_path.suffix}"
            ref_dest.write_bytes(ref_path.read_bytes())
            ligand_reference_rel = str(ref_dest.relative_to(root))

        param_status = (
            prep_result.parameterization_status
            if prep_result.success
            else "chemical_validation_required"
        )

        manifest_ligands[ligand_id] = LigandManifestEntry(
            resname=group.resname,
            occurrences=len(group.system_ids),
            identity_method=group.identity_method,
            identity_confidence=group.identity_confidence,
            hydrogenation_status=prep_result.hydrogenation_status if prep_result.success else "unknown",
            hydrogenation_performed=prep_result.hydrogenation_performed if prep_result.success else False,
            hydrogenation_backend=prep_result.hydrogenation_backend if prep_result.success else "none",
            parameterization_status=param_status,
            element_corrections=prep_result.element_corrections if prep_result.success else [],
            element_assignment_confidence=prep_result.element_assignment_confidence if prep_result.success else "unknown",
            connectivity_confidence=prep_result.connectivity_confidence if prep_result.success else "unknown",
            bond_order_confidence=prep_result.bond_order_confidence if prep_result.success else "unknown",
            chemical_identity_method=prep_result.chemical_identity_method if prep_result.success else "none",
            chemical_identity_confidence=prep_result.chemical_identity_confidence if prep_result.success else "unknown",
            readiness_block_reasons=prep_result.readiness_block_reasons if prep_result.success else [],
            ligand_reference=ligand_reference_rel,
            ligand_for_ligpargen=str(ligand_for_ligpargen.relative_to(root)) if ligand_for_ligpargen else None,
            ligand_report=str(ligand_report.relative_to(root)) if ligand_report else None,
            representative_system_id=rep_system,
        )

    # ── Phase 5/6: systems + manifest ────────────────────────────────────
    system_to_ligand = {sid: lid for lid, g in groups.items() for sid in g.system_ids}

    manifest_systems: list[SystemManifestEntry] = []
    for entry in complexes:
        if entry.system_id not in extracted or entry.system_id not in system_to_ligand:
            continue
        result = extracted[entry.system_id]
        manifest_systems.append(SystemManifestEntry(
            id=entry.system_id,
            directory=str(entry.directory.relative_to(root)),
            complex_pdb=str(entry.complex_pdb.relative_to(root)),
            protein=str(result.protein_pdb.relative_to(root)),
            ligand_pose=str(result.ligand_pose_pdb.relative_to(root)),
            ligand_id=system_to_ligand[entry.system_id],
            ligand_resname=ligand_resname,
            protein_atom_count=result.protein_atom_count,
            ligand_pose_atom_count=result.ligand_atom_count,
        ))

    manifest = CampaignManifest(
        campaign_id=ligand_resname,
        root=str(root.resolve()),
        systems=manifest_systems,
        ligands=manifest_ligands,
        warnings=warnings,
    )
    manifest.write(manifest_path)
    report_path.write_text(_render_prepare_report(manifest, system_errors))

    return CampaignPrepareResult(
        success=True,
        manifest=manifest,
        manifest_path=manifest_path,
        report_path=report_path,
        systems_discovered=len(complexes),
        unique_ligands=len(groups),
        system_errors=system_errors,
        warnings=warnings,
    )


def _render_prepare_report(manifest: CampaignManifest, system_errors: dict[str, str]) -> str:
    lines = [
        "# SimForge Ligand Campaign — Prepare Report",
        "",
        f"- Campaign: {manifest.campaign_id}",
        f"- Systems discovered: {len(manifest.systems)}",
        f"- Unique ligands: {len(manifest.ligands)}",
        "",
    ]
    for ligand_id, lig in manifest.ligands.items():
        lines += [
            f"## {ligand_id}",
            "",
            f"- Occurrences: {lig.occurrences}",
            f"- Identity method: {lig.identity_method} ({lig.identity_confidence})",
            f"- Elements: {lig.element_assignment_confidence}"
            + (f" ({len(lig.element_corrections)} correction(s))" if lig.element_corrections else ""),
            f"- Connectivity: {lig.connectivity_confidence}",
            f"- Bond orders: {lig.bond_order_confidence}",
            f"- Hydrogens: {lig.hydrogenation_status}",
            f"- Hydrogenation performed: {lig.hydrogenation_performed} ({lig.hydrogenation_backend})",
            f"- Chemical validation: {lig.chemical_identity_confidence} ({lig.chemical_identity_method})",
            f"- Parameterization: {lig.parameterization_status}",
        ]
        if lig.readiness_block_reasons:
            lines.append("- Blocked because:")
            lines += [f"    - {r}" for r in lig.readiness_block_reasons]
        if lig.ligand_for_ligpargen:
            lines.append(f"- Upload to LigParGen: `{lig.ligand_for_ligpargen}`")
        lines.append("")

    if system_errors:
        lines.append("## System errors")
        lines.append("")
        for sid, err in system_errors.items():
            lines.append(f"- {sid}: {err}")
        lines.append("")

    if manifest.warnings:
        lines.append("## Warnings")
        lines.append("")
        for w in manifest.warnings:
            lines.append(f"- {w}")
        lines.append("")

    return "\n".join(lines)
