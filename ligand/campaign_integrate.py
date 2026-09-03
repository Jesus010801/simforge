"""
Batch integration for ligand campaigns.

Validates each unique ligand's LigParGen output exactly once, then
assembles a GROMACS system per receptor by reusing the existing
single-complex assembly pipeline. Does not duplicate pdb2gmx invocation,
topology splitting, atomtypes extraction, GRO merging, or grompp
validation — all of that lives in ``ligand.integrate.assemble_system`` and
is called once per system here.

Public API:
  SystemIntegrationOutcome
  LigandIntegrationOutcome
  CampaignIntegrateResult
  integrate_campaign(root, ...) -> CampaignIntegrateResult
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from ligand.campaign import _CAMPAIGN_DIR_NAME, _SYSTEM_SUBDIR_NAME, CampaignManifest
from ligand.integrate import assemble_system
from ligand.ligpargen_import_validator import LigParGenImportValidator
from ligand.normalization import _GROMACS_NAME_RE, LigandIdentity
from ligand.pose_rewriter import LigandPoseRewriter


@dataclass
class SystemIntegrationOutcome:
    system_id: str
    success: bool
    pose_gro: Optional[Path] = None
    system_dir: Optional[Path] = None
    error: Optional[str] = None
    warnings: list[str] = field(default_factory=list)


@dataclass
class LigandIntegrationOutcome:
    ligand_id: str
    success: bool
    normalized_itp: Optional[Path] = None
    normalized_gro: Optional[Path] = None
    error: Optional[str] = None
    warnings: list[str] = field(default_factory=list)
    systems: dict[str, SystemIntegrationOutcome] = field(default_factory=dict)


@dataclass
class CampaignIntegrateResult:
    success: bool
    manifest: Optional[CampaignManifest] = None
    ligands: dict[str, LigandIntegrationOutcome] = field(default_factory=dict)
    systems_total: int = 0
    systems_succeeded: int = 0
    error: Optional[str] = None


def _find_ligpargen_files(
    ligpargen_dir: Path, ligand_id: str
) -> tuple[Optional[Path], Optional[Path], Optional[str]]:
    """
    Locate the .itp/.gro pair the user placed in a ligand's ligpargen/ dir.

    Prefers files literally named ``<ligand_id>.itp`` / ``<ligand_id>.gro``;
    falls back to accepting a single .itp/.gro pair with any name.
    """
    if not ligpargen_dir.exists():
        return None, None, f"LigParGen output directory not found: {ligpargen_dir}"

    preferred_itp = ligpargen_dir / f"{ligand_id}.itp"
    preferred_gro = ligpargen_dir / f"{ligand_id}.gro"
    if preferred_itp.exists() and preferred_gro.exists():
        return preferred_itp, preferred_gro, None

    itps = sorted(p for p in ligpargen_dir.glob("*.itp") if p.is_file())
    gros = sorted(p for p in ligpargen_dir.glob("*.gro") if p.is_file())
    if len(itps) == 1 and len(gros) == 1:
        return itps[0], gros[0], None

    return None, None, (
        f"Could not find a single .itp/.gro pair in {ligpargen_dir}. "
        f"Expected '{ligand_id}.itp' and '{ligand_id}.gro' (or exactly one of "
        f"each). Found itp={[p.name for p in itps]}, gro={[p.name for p in gros]}."
    )


def integrate_campaign(
    root: str | Path,
    forcefield: str = "oplsaa",
    water_model: str = "spce",
    run_grompp: bool = True,
) -> CampaignIntegrateResult:
    """
    Validate each unique ligand's LigParGen output once, reconstruct a
    receptor-specific ligand pose for every system in its group, and call
    ``assemble_system`` independently per system.

    Failures in one system (or one ligand's validation) are recorded against
    that system/ligand id and do not abort the rest of the campaign.
    """
    root = Path(root)
    manifest_path = root / _CAMPAIGN_DIR_NAME / "campaign_manifest.json"

    if not manifest_path.exists():
        return CampaignIntegrateResult(
            success=False,
            error=(
                f"No campaign manifest found at {manifest_path}. "
                "Run 'simforge ligand batch-prepare' first."
            ),
        )

    try:
        manifest = CampaignManifest.load(manifest_path)
    except Exception as exc:
        return CampaignIntegrateResult(success=False, error=f"Cannot read campaign manifest: {exc}")

    path_errors = manifest.validate_paths(root)
    if path_errors:
        return CampaignIntegrateResult(
            success=False,
            manifest=manifest,
            error=(
                "Campaign manifest references missing files:\n"
                + "\n".join(f"  - {e}" for e in path_errors)
            ),
        )

    systems_by_ligand: dict[str, list] = {}
    for s in manifest.systems:
        systems_by_ligand.setdefault(s.ligand_id, []).append(s)

    ligand_outcomes: dict[str, LigandIntegrationOutcome] = {}
    systems_total = 0
    systems_succeeded = 0

    for ligand_id, systems in systems_by_ligand.items():
        outcome = LigandIntegrationOutcome(ligand_id=ligand_id, success=False)
        ligand_outcomes[ligand_id] = outcome

        if not _GROMACS_NAME_RE.match(ligand_id):
            outcome.error = (
                f"Ligand id '{ligand_id}' is not GROMACS-safe (must be 1-5 "
                "uppercase alphanumeric characters) and cannot be used as the "
                "residue/moleculetype name."
            )
            systems_total += len(systems)
            continue

        ligpargen_dir = root / _CAMPAIGN_DIR_NAME / "ligands" / ligand_id / "ligpargen"
        itp_path, gro_path, err = _find_ligpargen_files(ligpargen_dir, ligand_id)
        if err:
            outcome.error = err
            systems_total += len(systems)
            continue

        identity = LigandIdentity(
            component_id=ligand_id.lower(),
            display_name=ligand_id,
            source_filename=itp_path.name,
            internal_id=ligand_id,
            residue_name=ligand_id,
            moleculetype=ligand_id,
        )
        validation = LigParGenImportValidator().validate(
            ligand_gro=gro_path,
            ligand_itp=itp_path,
            identity=identity,
            work_dir=ligpargen_dir,
        )
        if not validation.valid:
            outcome.error = "LigParGen output validation failed: " + "; ".join(validation.errors)
            outcome.warnings = validation.warnings
            systems_total += len(systems)
            continue

        outcome.normalized_itp = validation.itp_path
        outcome.normalized_gro = validation.gro_path
        outcome.warnings = list(validation.warnings)

        rewriter = LigandPoseRewriter()
        for s in systems:
            systems_total += 1
            sys_outcome = SystemIntegrationOutcome(system_id=s.id, success=False)
            outcome.systems[s.id] = sys_outcome

            system_dir = root / s.directory / _SYSTEM_SUBDIR_NAME
            pose_result = rewriter.rewrite_heavy_atom_transfer(
                reference_gro=validation.gro_path,
                reference_itp=validation.itp_path,
                target_pose_pdb=root / s.ligand_pose,
                output_dir=system_dir,
            )
            if not pose_result.success:
                sys_outcome.error = f"Pose reconstruction failed: {pose_result.error}"
                continue
            sys_outcome.pose_gro = pose_result.output_path
            sys_outcome.warnings.extend(pose_result.warnings)

            assembly_out_dir = system_dir / "system"
            assembly = assemble_system(
                protein_pdb=root / s.protein,
                ligand_itp=validation.itp_path,
                ligand_gro=pose_result.output_path,
                out_dir=assembly_out_dir,
                forcefield=forcefield,
                water_model=water_model,
                run_grompp=run_grompp,
            )
            sys_outcome.system_dir = assembly_out_dir
            if assembly.error:
                sys_outcome.error = assembly.error
                continue
            if not assembly.success:
                sys_outcome.error = "; ".join(assembly.errors) or "Assembly validation failed"
                continue

            sys_outcome.success = True
            sys_outcome.warnings.extend(assembly.warnings)
            systems_succeeded += 1

        outcome.success = bool(outcome.systems) and all(so.success for so in outcome.systems.values())

    overall_success = systems_total > 0 and systems_succeeded == systems_total

    return CampaignIntegrateResult(
        success=overall_success,
        manifest=manifest,
        ligands=ligand_outcomes,
        systems_total=systems_total,
        systems_succeeded=systems_succeeded,
    )
