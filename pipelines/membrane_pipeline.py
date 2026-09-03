# pipelines/membrane_pipeline.py
"""
MembraneWorkflowOPLSAA — pipeline de membrana DPPC + OPLS-AA.

Captura fielmente el workflow experto de tutorial_membrana.txt sin
generalización prematura.  Cada step refleja una decisión científica
real; los parámetros físicos vienen de core/membrane_knowledge.py.

DAG generado (Phase 11 topology architecture):
    generate_protein_topology  PREPARATION   AUTO   (pdb2gmx on original protein PDB only)
    orient_protein             PREPARATION   AUTOMATED (with structural_annotation) / GUIDED (without)
    match_box_to_bilayer       ASSEMBLY      AUTOMATED (MatchBoxBuilder — box_match_report gate)
    embed_in_bilayer           ASSEMBLY      AUTOMATED (MoveMembAdapter Python — no gfortran)
    generate_topology          PREPARATION   AUTO   (assemble embed-time topology from components)
    membrane_embedding         MEMBRANE_EMBEDDING AUTO  (shrink loop)
    assemble_system_topology   ASSEMBLY      AUTO   (final topol.top, counts from converged.gro)
    solvate_membrane           ASSEMBLY      AUTO
    clean_water                ASSEMBLY      AUTOMATED (WaterDeletorAdapter Python)
    add_ions                   ASSEMBLY      AUTO
    energy_minimization        MINIMIZATION  AUTO   (-DPOSRES -DSTRONG_POSRES)
    equilibration              EQUILIBRATION AUTO   (semiisotropic, Berendsen NPT)
    production_md              PRODUCTION    AUTO   (NH+PR, dt=0.001, semiisotropic)
    analysis_*                 ANALYSIS      AUTO

NO generalizar antes de que este pipeline funcione end-to-end con DPPC.
"""

from __future__ import annotations

from pipelines.base_pipeline import BasePipeline
from core.models import SystemState
from core.execution_models import (
    AutomationLevel,
    SimulationPlan,
    SimulationStep,
    PlanStatus,
    StepStage,
    StepType,
    WorkflowPolicy,
)
from core import membrane_knowledge as mk


class MembraneWorkflowOPLSAA(BasePipeline):

    pipeline_type = "protein-membrane"

    def build_plan(self, state: SystemState) -> SimulationPlan:

        mem    = state.environment.membrane
        lipid  = (mem.type or "DPPC").upper()
        ff     = state.forcefields.protein or "opls-aa"
        T      = state.environment.temperature_K
        wm     = state.environment.solvent.water_model

        # ── Physical constants from knowledge layer ────────────────────────────
        try:
            apl_target  = mk.apl_target(lipid, ff, T)
        except KeyError:
            apl_target  = 62.0

        try:
            lipid_resname = mk.lipid_residue_name(lipid, ff)
        except KeyError:
            lipid_resname = "DPP"

        try:
            atom_names = mk.lipid_atom_names(lipid, ff)
        except KeyError:
            atom_names = mk.LipidAtomNames(headgroup_ref="O33", tail_middle="C50")

        inflate_f = mk.inflation_factor("single_pass_tm")
        defaults  = mk.MEMBRANE_EQUILIBRATION_DEFAULTS

        # ── Embedding config (user-overridable) ───────────────────────────────
        _mem_cfg = (getattr(state, "config", None) or {}).get("membrane", {})
        embed_cfg = _mem_cfg.get("embedding", {})
        # Default pinned to the original protmemfiles tutorial's one-shot overlap-
        # removal cutoff (docs/Prot-Memb_FILES/tutorial_membrana.txt:36 — cutoff
        # arg "14" -> 1.4 nm). A previous default of 0.14 nm was 10x too small and
        # left protein-clashing lipids unremoved — see
        # docs/audits/simforge_vs_protmemfiles_audit.md.
        lipid_exclusion_cutoff_nm  = float(embed_cfg.get("lipid_exclusion_cutoff_nm", 1.4))
        # Shrink-loop (per-iteration deflate) overlap cutoff. The original
        # protmemfiles method applies the overlap-removal cutoff exactly once
        # (docs/Prot-Memb_FILES/tutorial_membrana.txt:36), then runs every
        # subsequent shrink-loop deflate call with cutoff=0
        # (docs/Prot-Memb_FILES/ScriptCamilo-Jorge.sh:18,38 and
        # run_inflategro.sh:41,48 — literal cutoff argument "0"). Reusing the
        # nonzero one-shot cutoff on every deflate iteration causes the loop to
        # keep deleting newly-compressed annular lipids as the box shrinks,
        # producing a lipid-free void around the receptor even though the
        # global APL converges — see docs/audits/ for the GLP-1R evidence
        # (478 -> 400 DPP over 29 iterations, entirely inside the loop).
        shrink_loop_cutoff_nm      = float(embed_cfg.get("shrink_loop_cutoff_nm", 0.0))
        trapped_lipid_policy       = embed_cfg.get("trapped_lipid_policy", "warn")
        annular_repair_enabled     = bool(embed_cfg.get("annular_repair_enabled", False))
        annular_repair_allow_insertion = bool(embed_cfg.get("annular_repair_allow_insertion", False))
        pre_exclude_trapped_lipids = bool(embed_cfg.get("pre_exclude_trapped_lipids", True))
        max_pre_excluded_lipids    = int(embed_cfg.get("max_pre_excluded_lipids", 50))
        save_embedding_iterations  = embed_cfg.get("save_embedding_iterations", "final_only")
        save_embedding_iteration_stride = int(
            embed_cfg.get("save_embedding_iteration_stride", 5)
        )
        quality_diagnostics        = bool(embed_cfg.get("quality_diagnostics", True))
        quality_analysis_region    = embed_cfg.get("quality_analysis_region", "tm_footprint")
        quality_footprint_margin   = float(embed_cfg.get("quality_footprint_margin_nm", 0.5))
        tm_mask_validation         = bool(embed_cfg.get("tm_mask_validation", True))
        tm_mask_policy             = str(embed_cfg.get("tm_mask_policy", "warn"))
        tm_aware_exclusion         = bool(embed_cfg.get("tm_aware_exclusion", False))
        tm_aware_exclusion_policy  = str(embed_cfg.get("tm_aware_exclusion_policy", "apply" if tm_aware_exclusion else "off"))
        if not tm_aware_exclusion:
            tm_aware_exclusion_policy = "off"
        tm_aware_max_removed       = int(embed_cfg.get("tm_aware_max_removed_lipids", 50))

        # ── Phase 9B: Optimizer config (passed through from embed_cfg) ─────────
        optimizer_enabled          = embed_cfg.get("optimizer_enabled")   # None → builder auto-detects
        optimizer_policy           = embed_cfg.get("optimizer_policy")
        optimizer_max_candidates   = embed_cfg.get("optimizer_max_candidates")
        optimizer_z_shifts         = embed_cfg.get("optimizer_z_shift_offsets_nm")
        optimizer_mask_paddings    = embed_cfg.get("optimizer_mask_padding_nm")
        optimizer_min_improvement  = embed_cfg.get("optimizer_min_score_improvement")
        optimizer_run_minimization = embed_cfg.get("optimizer_run_candidate_minimization")
        optimizer_excl_safety      = embed_cfg.get("optimizer_exclusion_safety_limit")

        # ── Backend & Quality gates config ────────────────────────────────────
        backend = embed_cfg.get("backend", "inflategro")
        if backend not in ("inflategro", "tm_aware"):
            raise ValueError(
                f"Invalid membrane.embedding.backend: '{backend}'. "
                "Allowed values: 'inflategro', 'tm_aware'."
            )

        quality_gates = embed_cfg.get("quality_gates", {})
        max_void_fraction_local = float(quality_gates.get("max_void_fraction_local", 0.25))
        max_void_area_local_nm2 = float(quality_gates.get("max_void_area_local_nm2", 1.0))
        min_tm_burial_score = float(quality_gates.get("min_tm_burial_score", 0.80))
        max_trapped_lipids = int(quality_gates.get("max_trapped_lipids", 0))
        max_soluble_domain_core_atoms = int(quality_gates.get("max_soluble_domain_core_atoms", 10))

        if backend == "tm_aware":
            tm_aware_exclusion = True
            tm_aware_exclusion_policy = "strict_apply"
            tm_mask_validation = True
            tm_mask_policy = "strict"

        # ── Phase 10A: Interface builder (lipid refill) config ────────────────
        _iface_b_cfg = embed_cfg.get("interface_builder", {})
        iface_builder_enabled      = bool(_iface_b_cfg.get("enabled", False))
        iface_builder_policy       = str( _iface_b_cfg.get("policy", "warn"))
        iface_builder_max_inserted = int( _iface_b_cfg.get("max_inserted_lipids", 40))
        iface_builder_max_per_cl   = int( _iface_b_cfg.get("max_lipids_per_gap_cluster", 4))
        iface_builder_target_dist  = float(_iface_b_cfg.get("target_contact_distance_nm", 0.45))
        iface_builder_prot_clash   = float(_iface_b_cfg.get("protein_clash_cutoff_nm", 0.20))
        iface_builder_lip_clash    = float(_iface_b_cfg.get("lipid_clash_cutoff_nm", 0.18))
        iface_builder_seed         = int(  _iface_b_cfg.get("deterministic_seed", 17))

        # ── Production duration ────────────────────────────────────────────────
        if state.environment.duration_ns is not None:
            prod_ns = state.environment.duration_ns
        else:
            prod_ns = defaults.prod_nsteps * defaults.prod_dt / 1000.0  # 500ns

        prod_nsteps = int(prod_ns * 1000.0 / defaults.prod_dt)

        # ── Build plan ────────────────────────────────────────────────────────
        plan = SimulationPlan(
            status=PlanStatus.READY,
            inferred_system_type=state.inferred_system_type,
            workflow_policy=WorkflowPolicy(
                temperature_K=T,
                production_time_ns=prod_ns,
            ),
        )

        # ── Protein component id ──────────────────────────────────────────────
        protein = next(
            (c for c in state.components if c.role in ("protein", "peptide")),
            None,
        )
        prot_id   = protein.id   if protein else "protein_1"
        prot_file = protein.file if protein else "protein.pdb"

        # ─────────────────────────────────────────────────────────────────────
        # Step 0: generate protein topology (pdb2gmx on original protein PDB)
        # Phase 11: pdb2gmx MUST run on the original PDB (chain/TER/terminus
        # semantics preserved).  It must NEVER run on embed_in_bilayer/system.gro
        # (mixed protein-lipid GRO loses those semantics → OXT terminus error).
        # ─────────────────────────────────────────────────────────────────────
        plan.steps.append(SimulationStep(
            step_id="generate_protein_topology",
            title="Generar topología proteína (pdb2gmx sobre PDB original)",
            stage=StepStage.PREPARATION,
            step_type=StepType.AUTOMATIC,
            automation_level=AutomationLevel.AUTOMATED,
            engine="gromacs:pdb2gmx_protein",
            target_components=[prot_id],
            params={
                "source_file": prot_file,
                "forcefield":  "opls-aa-membrane",
                "water_model": "none",
                "note": (
                    "pdb2gmx runs ONLY on the original protein PDB.\n"
                    "NEVER on embed_in_bilayer/system.gro — that file is mixed "
                    "protein+lipid and loses chain/TER/terminus semantics."
                ),
            },
            notes=[
                "Phase 11: pdb2gmx on original protein PDB only",
                "Outputs: protein_processed.gro, topol.top, posre.itp, protein_topology_manifest.json",
            ],
        ))

        # ─────────────────────────────────────────────────────────────────────
        # Step 1: orient protein
        #   AUTOMATED when structural_annotation has EC + IC topology (new format)
        #             OR environment.membrane.orientation has the legacy fields.
        #   GUIDED    otherwise (user must rotate manually and provide GRO).
        # ─────────────────────────────────────────────────────────────────────
        sa = state.structural_annotation
        ec_residues: str | None = None
        ic_residues: str | None = None
        tm_residues: str | None = None
        ec_target_side: str = "+z"  # GROMACS convention default

        if sa is not None and sa.membrane_topology is not None:
            mt = sa.membrane_topology
            # Join multiple EC/IC regions with commas — structural_annotation
            # supports multi-segment topology, builder uses flat range string.
            if mt.extracellular_regions:
                ec_residues = ",".join(mt.extracellular_regions)
            if mt.intracellular_regions:
                ic_residues = ",".join(mt.intracellular_regions)
            if mt.transmembrane_segments:
                from core.structural_annotation import TransmembraneSegment
                tm_parts = [
                    (s.residues if isinstance(s, TransmembraneSegment) else s)
                    for s in mt.transmembrane_segments
                ]
                tm_residues = ",".join(tm_parts)
            # Explicit orientation axis from structural_annotation
            if sa.orientation is not None:
                ec_target_side = sa.orientation.extracellular_side
        else:
            # Backward compat: read legacy environment.membrane.orientation
            legacy = getattr(mem, "orientation", None)
            if legacy is not None:
                ec_residues = getattr(legacy, "extracellular_residues", None)
                ic_residues = getattr(legacy, "intracellular_residues", None)
                tm_residues = getattr(legacy, "tm_segments", None)
                # Legacy format had no axis declaration → assume +z

        has_orient = bool(ec_residues and ic_residues)

        # ── Registration config (Phase 8A) ────────────────────────────────────
        # Must come after tm_residues/ec_residues/ic_residues are resolved above.
        reg_cfg = _mem_cfg.get("registration", {})
        has_tm_annotation = bool(tm_residues)
        adaptive_registration_default = has_tm_annotation and backend == "tm_aware"
        adaptive_registration   = bool(reg_cfg.get("adaptive", adaptive_registration_default))
        reg_search_min          = float(reg_cfg.get("search_min_nm",       -1.5))
        reg_search_max          = float(reg_cfg.get("search_max_nm",        1.5))
        reg_search_step         = float(reg_cfg.get("search_step_nm",       0.1))
        reg_policy              = str(reg_cfg.get("policy", "warn"))
        reg_max_allowed_shift   = float(reg_cfg.get("max_allowed_shift_nm", 1.5))
        reg_hc_half_thickness   = float(reg_cfg.get("hydrophobic_half_thickness_nm", 1.25))

        orient_step_type = StepType.AUTOMATIC if has_orient else StepType.MANUAL
        orient_params: dict = {"source_file": prot_file}

        if has_orient:
            orient_params["_orientation"] = {
                "extracellular_residues": ec_residues,
                "intracellular_residues": ic_residues,
                "tm_segments":            tm_residues,
                "extracellular_side":     ec_target_side,
            }
        else:
            orient_params["note"] = (
                "1. gmx editconf -f protein.pdb -o protein_princ.gro -c -d 1.5 -bt triclinic -princ\n"
                "2. Inspect in VMD/PyMOL — identify rotation needed to align TM helix with Z\n"
                "3. gmx editconf -f protein_princ.gro -o protein_oriented.gro -rotate 0 ROT 0 -c -d 1.5 -bt triclinic -princ\n"
                "   Tip: add structural_annotation.membrane_topology to YAML "
                "to make this step automatic."
            )

        plan.steps.append(SimulationStep(
            step_id="orient_protein",
            title="Orientar proteína (eje TM alineado con Z)",
            stage=StepStage.PREPARATION,
            step_type=orient_step_type,
            automation_level=(
                AutomationLevel.AUTOMATED if has_orient else AutomationLevel.GUIDED
            ),
            engine="gromacs:editconf+orient",
            target_components=[prot_id],
            params=orient_params,
            notes=(
                ["Rotation computed from EC/IC Cα centres of mass via structural_annotation"]
                if has_orient else
                ["Add structural_annotation.membrane_topology to YAML to enable automatic rotation"]
            ),
        ))

        # ─────────────────────────────────────────────────────────────────────
        # Bilayer selection — authoritative for XY box dimensions (Phase 3).
        # Computed here so both match_box_to_bilayer and embed_in_bilayer
        # receive the same bilayer_file.
        # ─────────────────────────────────────────────────────────────────────
        bilayer      = mk.bilayer_for_box(12.84, 12.89, lipid)
        bilayer_file = bilayer.filename if bilayer else "dppc512_whole.gro"

        # ─────────────────────────────────────────────────────────────────────
        # Step 2: match box to bilayer (automated — MatchBoxBuilder computes
        # bounding box from protein_oriented.gro; selected bilayer GRO is
        # authoritative for XY box dimensions; Z remains protein-aware)
        # ─────────────────────────────────────────────────────────────────────
        match_box_params: dict = {
            "lipid":        lipid,
            "bilayer_file": bilayer_file,   # baked into helper for XY authority
            # Phase 8A: adaptive membrane registration
            "adaptive_registration":          adaptive_registration,
            "reg_search_min_nm":              reg_search_min,
            "reg_search_max_nm":              reg_search_max,
            "reg_search_step_nm":             reg_search_step,
            "reg_policy":                     reg_policy,
            "reg_max_allowed_shift_nm":       reg_max_allowed_shift,
            "reg_hydrophobic_half_thickness_nm": reg_hc_half_thickness,
        }
        # Bake structural annotation into registration helper for TM/EC/IC scoring
        if tm_residues:
            match_box_params["tm_residues"] = tm_residues
        if ec_residues:
            match_box_params["ec_residues"] = ec_residues
        if ic_residues:
            match_box_params["ic_residues"] = ic_residues

        reg_note = (
            f"Phase 8A: adaptive membrane registration enabled "
            f"(search {reg_search_min:+.1f}→{reg_search_max:+.1f} nm, step {reg_search_step} nm)"
            if adaptive_registration else
            "Phase 8A: adaptive membrane registration disabled (no TM annotation or user override)"
        )
        plan.steps.append(SimulationStep(
            step_id="match_box_to_bilayer",
            title="Calcular y validar caja de bicapa",
            stage=StepStage.ASSEMBLY,
            step_type=StepType.AUTOMATIC,
            automation_level=AutomationLevel.AUTOMATED,
            engine="gromacs:box_match",
            depends_on=["orient_protein"],
            params=match_box_params,
            notes=[
                "Reads protein_oriented.gro → computes protein bounding box",
                f"XY = selected bilayer ({bilayer_file}) box XY — bilayer is authoritative",
                "Z = bilayer_thickness + protein_Z + 2×solvent_padding (protein-aware)",
                f"Lipid: {lipid} — parameters from core/bilayer_geometry.py",
                "Validates protein fits in bilayer with ≥1 nm periodic-image margin",
                "Produces box_match_report.json (gate) + protein_boxed.gro",
                reg_note,
            ],
        ))

        # ─────────────────────────────────────────────────────────────────────
        # Step 3: embed in bilayer
        # TM residues (if annotated) drive the bilayer Z alignment so that
        # large EC/IC domains do not mis-center the hydrophobic core.
        # When absent, the adapter falls back to protein Z-centre and emits a
        # warning in the output report.
        # ─────────────────────────────────────────────────────────────────────
        embed_params: dict = {
            "bilayer_file":              bilayer_file,
            "lipid":                     lipid,
            "lipid_residue_name":        lipid_resname,
            "pre_exclude_trapped_lipids": pre_exclude_trapped_lipids,
            "max_pre_excluded_lipids":    max_pre_excluded_lipids,
            "tm_mask_validation":         tm_mask_validation,
            "tm_mask_policy":             tm_mask_policy,
            "tm_aware_exclusion":         tm_aware_exclusion,
            "tm_aware_exclusion_policy":  tm_aware_exclusion_policy,
            "tm_aware_max_removed_lipids": tm_aware_max_removed,
        }
        if tm_residues:
            embed_params["tm_residues"] = tm_residues  # range string, e.g. "51-75"

        # Phase 10A: always forward interface_builder params so the builder can
        # bake them into run_interface_refill.py
        embed_params["iface_builder_enabled"]    = iface_builder_enabled
        embed_params["iface_builder_policy"]     = iface_builder_policy
        embed_params["iface_builder_max_inserted"]  = iface_builder_max_inserted
        embed_params["iface_builder_max_per_cluster"] = iface_builder_max_per_cl
        embed_params["iface_builder_target_dist"]    = iface_builder_target_dist
        embed_params["iface_builder_prot_clash"]     = iface_builder_prot_clash
        embed_params["iface_builder_lip_clash"]      = iface_builder_lip_clash
        embed_params["iface_builder_seed"]           = iface_builder_seed

        embed_notes = [
            "MoveMembAdapter (Python) aligns bilayer midplane to TM-region Z-centre — no gfortran required"
            if tm_residues else
            "MoveMembAdapter: no TM annotation — falling back to full protein Z-centre; "
            "add transmembrane_segments to structural_annotation for accurate placement"
        ]
        plan.steps.append(SimulationStep(
            step_id="embed_in_bilayer",
            title="Embutir proteína en bicapa lipídica",
            stage=StepStage.ASSEMBLY,
            step_type=StepType.AUTOMATIC,
            automation_level=AutomationLevel.AUTOMATED,
            engine="movememb+genrestr",
            depends_on=["match_box_to_bilayer"],
            params=embed_params,
            notes=embed_notes,
        ))

        # ─────────────────────────────────────────────────────────────────────
        # Step 4: membrane embedding shrink loop / TM-aware embedding (meta-step)
        #
        # Depends on embed_in_bilayer and generate_protein_topology.
        # assemble_system_topology is NOT a dependency — it runs AFTER embedding.
        # At runtime the embedding step calls run_assemble_system.py with
        # embed_in_bilayer/system.gro to build a bootstrap topology for grompp.
        # ─────────────────────────────────────────────────────────────────────
        embedding_title = (
            f"TM-aware embedding & relaxation ({lipid})"
            if backend == "tm_aware"
            else f"Shrink loop — convergencia APL ({lipid}, target ≤ {apl_target + mk.APL_CONVERGENCE_TOLERANCE:.0f} Å²)"
        )
        embedding_engine = (
            "gromacs:tm_aware_embedding"
            if backend == "tm_aware"
            else "gromacs+perl:inflategro"
        )
        embedding_params = {
            "lipid":               lipid,
            "lipid_residue_name":  lipid_resname,
            "forcefield":          ff,
            "water_model":         wm,
            "temperature_K":       T,
            "apl_target_ang2":     apl_target,
            "apl_tolerance_ang2":  mk.APL_CONVERGENCE_TOLERANCE,
            "inflate_factor":      inflate_f,
            "deflate_factor":      mk.SHRINK_DEFLATION_FACTOR,
            "max_iterations":      mk.SHRINK_MAX_ITERATIONS,
            "gridsize":            5,
            "cutoff":              lipid_exclusion_cutoff_nm,
            "shrink_loop_cutoff_nm": shrink_loop_cutoff_nm,
            "trapped_lipid_policy": trapped_lipid_policy,
            "annular_repair_enabled": annular_repair_enabled,
            "annular_repair_allow_insertion": annular_repair_allow_insertion,
            "inflategro_script":   "inflategro-Jorge.pl",
            "input_gro":           "../embed_in_bilayer/system.gro",
            "topol_top":           "bootstrap_topol.top",  # generated at runtime by run_assemble_system.py
            # Iteration snapshot config
            "save_embedding_iterations":       save_embedding_iterations,
            "save_embedding_iteration_stride": save_embedding_iteration_stride,
            # Post-convergence quality diagnostics
            "quality_diagnostics":       quality_diagnostics,
            "quality_analysis_region":   quality_analysis_region,
            "quality_footprint_margin_nm": quality_footprint_margin,
            # Backend & Quality gates config
            "backend":                       backend,
            "max_void_fraction_local":       max_void_fraction_local,
            "max_void_area_local_nm2":       max_void_area_local_nm2,
            "min_tm_burial_score":           min_tm_burial_score,
            "max_trapped_lipids":            max_trapped_lipids,
            "max_soluble_domain_core_atoms": max_soluble_domain_core_atoms,
            "tm_aware_max_removed_lipids":   tm_aware_max_removed,
        }
        if tm_residues:
            embedding_params["tm_residues"] = tm_residues

        # Phase 9B: only forward optimizer params when explicitly set in YAML
        if optimizer_enabled is not None:
            embedding_params["optimizer_enabled"] = optimizer_enabled
        if optimizer_policy is not None:
            embedding_params["optimizer_policy"] = optimizer_policy
        if optimizer_max_candidates is not None:
            embedding_params["optimizer_max_candidates"] = optimizer_max_candidates
        if optimizer_z_shifts is not None:
            embedding_params["optimizer_z_shift_offsets_nm"] = optimizer_z_shifts
        if optimizer_mask_paddings is not None:
            embedding_params["optimizer_mask_padding_nm"] = optimizer_mask_paddings
        if optimizer_min_improvement is not None:
            embedding_params["optimizer_min_score_improvement"] = optimizer_min_improvement
        if optimizer_run_minimization is not None:
            embedding_params["optimizer_run_candidate_minimization"] = optimizer_run_minimization
        if optimizer_excl_safety is not None:
            embedding_params["optimizer_exclusion_safety_limit"] = optimizer_excl_safety

        plan.steps.append(SimulationStep(
            step_id="membrane_embedding",
            title=embedding_title,
            stage=StepStage.MEMBRANE_EMBEDDING,
            step_type=StepType.AUTOMATIC,
            engine=embedding_engine,
            blocking=True,
            depends_on=["embed_in_bilayer", "generate_protein_topology"],
            params=embedding_params,
        ))

        # ─────────────────────────────────────────────────────────────────────
        # Step 5: assemble_system_topology (post-embedding)
        #
        # Runs AFTER membrane_embedding so final molecule counts derive from
        # membrane_embedding/converged.gro — the authoritative post-shrink
        # coordinate file evaluated by Phase 9A.
        #
        # Assembles topol.top from:
        #   - generate_protein_topology/ (protein, pdb2gmx on original PDB)
        #   - pdb2gmx on lipids_only.gro (lipid moleculetype; safe)
        #   - membrane_embedding/converged.gro (final molecule counts)
        # ─────────────────────────────────────────────────────────────────────
        plan.steps.append(SimulationStep(
            step_id="assemble_system_topology",
            title="Ensamblar topología final del sistema (post-embedding)",
            stage=StepStage.ASSEMBLY,
            step_type=StepType.AUTOMATIC,
            automation_level=AutomationLevel.AUTOMATED,
            engine="topology:assemble_system",
            depends_on=["membrane_embedding", "generate_protein_topology"],
            target_components=[prot_id],
            params={
                "forcefield":   "opls-aa-membrane",
                "water_model":  wm,
            },
            notes=[
                "Runs AFTER membrane_embedding — final counts from converged.gro",
                "pdb2gmx only on lipids_only.gro (safe: lipid-only input)",
                "Outputs: topol.top (authoritative) + topology_assembly_report.json",
            ],
        ))

        # ─────────────────────────────────────────────────────────────────────
        # Step 6: solvate membrane system
        # ─────────────────────────────────────────────────────────────────────
        plan.steps.append(SimulationStep(
            step_id="solvate_membrane",
            title="Solvatación del sistema membrana",
            stage=StepStage.ASSEMBLY,
            step_type=StepType.AUTOMATIC,
            engine="gromacs:solvate",
            depends_on=["assemble_system_topology"],
            params={
                "water_model": wm,
                "water_gro":   "spc216.gro",
                "input_gro":   "../membrane_embedding/converged.gro",
                "topol_top":   "../assemble_system_topology/topol.top",
            },
        ))

        # ─────────────────────────────────────────────────────────────────────
        # Step 7: remove water inside bilayer (water_deletor.pl)
        # ─────────────────────────────────────────────────────────────────────
        plan.steps.append(SimulationStep(
            step_id="clean_water",
            title=f"Eliminar agua interior de bicapa ({atom_names.headgroup_ref}/{atom_names.tail_middle})",
            stage=StepStage.ASSEMBLY,
            step_type=StepType.AUTOMATIC,
            automation_level=AutomationLevel.AUTOMATED,
            engine="python:water_deletor",
            depends_on=["solvate_membrane"],
            params={
                "ref_atom":    atom_names.headgroup_ref,
                "middle_atom": atom_names.tail_middle,
                "nwater":      3,
                "tm_residues": tm_residues,
            },
        ))

        # ─────────────────────────────────────────────────────────────────────
        # Step 8: add ions
        # ─────────────────────────────────────────────────────────────────────
        plan.steps.append(SimulationStep(
            step_id="add_ions",
            title="Neutralización e iones fisiológicos",
            stage=StepStage.ASSEMBLY,
            step_type=StepType.AUTOMATIC,
            engine="gromacs:genion",
            depends_on=["clean_water"],
            params={
                "concentration": state.environment.ions.concentration,
                "positive_ion":  state.environment.ions.positive,
                "negative_ion":  state.environment.ions.negative,
            },
        ))

        # ─────────────────────────────────────────────────────────────────────
        # Step 9: energy minimization (with strong position restraints)
        # ─────────────────────────────────────────────────────────────────────
        plan.steps.append(SimulationStep(
            step_id="energy_minimization",
            title="Minimización energética (POSRES + STRONG_POSRES)",
            stage=StepStage.MINIMIZATION,
            step_type=StepType.AUTOMATIC,
            engine="gromacs",
            blocking=True,
            depends_on=["add_ions"],
            params={
                "integrator": "steep",
                "emtol":      mk.SHRINK_LOOP_EMTOL,
                "emstep":     0.01,
                "nsteps":     50_000,
                "define":     "-DPOSRES -DSTRONG_POSRES",
                "rcoulomb":   1.2,
                "rvdw":       1.2,
                "disp_corr":  "EnerPres",
            },
        ))

        # ─────────────────────────────────────────────────────────────────────
        # Step 10: equilibration (NVT 100ps + NPT 1ns, semiisotropic)
        # ─────────────────────────────────────────────────────────────────────
        plan.steps.append(SimulationStep(
            step_id="equilibration",
            title="Equilibración NVT (100 ps) + NPT semiisotropic (1 ns)",
            stage=StepStage.EQUILIBRATION,
            step_type=StepType.AUTOMATIC,
            engine="gromacs",
            depends_on=["energy_minimization"],
            params={
                "dt":          defaults.nvt_dt,
                "nvt_nsteps":  defaults.nvt_nsteps,   # 25000 = 100ps
                "npt_nsteps":  defaults.npt_nsteps,   # 150000 = 1ns (Berendsen equilibration)
                "temperature": T,
                "tc_grps":     defaults.nvt_tc_grps,  # "system"
                "tau_t":       str(defaults.nvt_tau_t),
                "ref_t":       str(T),
                "constraints": "all-bonds",
                # Membrane-specific pressure coupling
                "pcoupltype":  defaults.npt_pcoupltype,    # "semiisotropic"
                "pcoupl_npt":  defaults.npt_pcoupl,        # "Berendsen"
                "ref_p_xy":    defaults.npt_ref_p_xy,      # 0.5
                "ref_p_z":     defaults.npt_ref_p_z,       # 0.5
                "tau_p":       defaults.npt_tau_p,         # 5.0
                "rcoulomb":    1.2,
                "rvdw":        1.2,
                "disp_corr":   "EnerPres",
            },
        ))

        # ─────────────────────────────────────────────────────────────────────
        # Step 11: production MD
        # ─────────────────────────────────────────────────────────────────────
        prod_tau_t = " ".join(["0.5"] * len(defaults.nvt_tc_grps.split()))
        plan.steps.append(SimulationStep(
            step_id="production_md",
            title=f"Producción MD ({prod_ns:.0f} ns, NH+PR, semiisotropic, dt={defaults.prod_dt})",
            stage=StepStage.PRODUCTION,
            step_type=StepType.AUTOMATIC,
            engine="gromacs",
            depends_on=["equilibration"],
            params={
                "dt":                defaults.prod_dt,      # 0.001 — OPLS-AA lipid constraint
                "nsteps":            prod_nsteps,
                "temperature":       T,
                "tc_grps":           defaults.nvt_tc_grps,  # "system" (single group)
                "tau_t":             prod_tau_t,
                "ref_t":             str(T),
                "constraints":       defaults.prod_constraints,  # "h-bonds"
                "tcoupl":            defaults.prod_tcoupl,       # "Nose-Hoover"
                # Semiisotropic Parrinello-Rahman
                "pcoupltype":        defaults.prod_pcoupltype,   # "semiisotropic"
                "ref_p_xy":          defaults.prod_ref_p_xy,     # 1.0
                "ref_p_z":           defaults.prod_ref_p_z,      # 1.0
                "tau_p":             defaults.prod_tau_p,        # 2.0
                "rcoulomb":          1.2,
                "rvdw":              1.2,
                "disp_corr":         "EnerPres",
                "nstxout_compressed": 20_000,
                "nstenergy":          20_000,
                "nstlog":             20_000,
            },
        ))

        # ─────────────────────────────────────────────────────────────────────
        # Steps 12+: analysis
        # ─────────────────────────────────────────────────────────────────────
        for analysis in state.analysis:
            plan.steps.append(SimulationStep(
                step_id=f"analysis_{analysis.type}",
                title=f"Análisis {analysis.type}",
                stage=StepStage.ANALYSIS,
                step_type=StepType.AUTOMATIC,
                engine="analysis_pipeline",
                depends_on=["production_md"],
                params={"analysis_type": analysis.type},
            ))

        plan.notes.append(f"Pipeline: MembraneWorkflowOPLSAA — {lipid} + {ff} + {wm}")
        plan.notes.append(f"APL target: ≤ {apl_target + mk.APL_CONVERGENCE_TOLERANCE:.0f} Å²  |  Production: {prod_ns:.0f} ns  |  T: {T} K")
        plan.notes.append(f"Steps generados: {len(plan.steps)}")

        return plan
