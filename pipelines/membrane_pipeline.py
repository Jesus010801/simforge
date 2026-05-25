# pipelines/membrane_pipeline.py
"""
MembraneWorkflowOPLSAA — pipeline de membrana DPPC + OPLS-AA.

Captura fielmente el workflow experto de tutorial_membrana.txt sin
generalización prematura.  Cada step refleja una decisión científica
real; los parámetros físicos vienen de core/membrane_knowledge.py.

DAG generado:
    orient_protein          PREPARATION   AUTOMATED (with structural_annotation) / GUIDED (without)
    match_box_to_bilayer    ASSEMBLY      AUTOMATED (MatchBoxBuilder — box_match_report gate)
    embed_in_bilayer        ASSEMBLY      AUTOMATED (MoveMembAdapter Python — no gfortran)
    generate_topology       PREPARATION   AUTO   (pdb2gmx oplsaa_membrane.ff)
    membrane_embedding      MEMBRANE_EMBEDDING AUTO  (shrink loop)
    solvate_membrane        ASSEMBLY      AUTO
    clean_water             ASSEMBLY      AUTOMATED (WaterDeletorAdapter Python)
    add_ions                ASSEMBLY      AUTO
    energy_minimization     MINIMIZATION  AUTO   (-DPOSRES -DSTRONG_POSRES)
    equilibration           EQUILIBRATION AUTO   (semiisotropic, Berendsen NPT)
    production_md           PRODUCTION    AUTO   (NH+PR, dt=0.001, semiisotropic)
    analysis_*              ANALYSIS      AUTO

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
            params=orient_params,
            notes=(
                ["Rotation computed from EC/IC Cα centres of mass via structural_annotation"]
                if has_orient else
                ["Add structural_annotation.membrane_topology to YAML to enable automatic rotation"]
            ),
        ))

        # ─────────────────────────────────────────────────────────────────────
        # Step 2: match box to bilayer (automated — MatchBoxBuilder computes
        # bounding box from protein_oriented.gro and recommends dimensions)
        # ─────────────────────────────────────────────────────────────────────
        plan.steps.append(SimulationStep(
            step_id="match_box_to_bilayer",
            title="Calcular y validar caja de bicapa",
            stage=StepStage.ASSEMBLY,
            step_type=StepType.AUTOMATIC,
            automation_level=AutomationLevel.AUTOMATED,
            engine="gromacs:box_match",
            depends_on=["orient_protein"],
            params={
                "lipid": lipid,
            },
            notes=[
                "Reads protein_oriented.gro → computes protein bounding box",
                "Recommends XY = footprint + 2×padding, Z = bilayer + protein + 2×solvent",
                f"Lipid: {lipid} — parameters from core/bilayer_geometry.py",
                "Produces box_match_report.json (gate) + protein_boxed.gro",
            ],
        ))

        # ─────────────────────────────────────────────────────────────────────
        # Step 3: embed in bilayer (manual — requires MoveMemb Z-displacement)
        # ─────────────────────────────────────────────────────────────────────
        bilayer     = mk.bilayer_for_box(12.84, 12.89, lipid)
        bilayer_file = bilayer.filename if bilayer else "dppc512_whole.gro"
        plan.steps.append(SimulationStep(
            step_id="embed_in_bilayer",
            title="Embutir proteína en bicapa lipídica",
            stage=StepStage.ASSEMBLY,
            step_type=StepType.AUTOMATIC,
            automation_level=AutomationLevel.AUTOMATED,
            engine="movememb+genrestr",
            depends_on=["match_box_to_bilayer"],
            params={
                "bilayer_file":       bilayer_file,
                "lipid":              lipid,
                "lipid_residue_name": lipid_resname,   # GRO resname (e.g. "DPP" for DPPC OPLS-AA)
            },
            notes=["MoveMembAdapter (Python) aligns bilayer Z-midplane to protein Z-centre — no gfortran required"],
        ))

        # ─────────────────────────────────────────────────────────────────────
        # Step 4: generate topology (pdb2gmx with oplsaa_membrane.ff)
        # ─────────────────────────────────────────────────────────────────────
        plan.steps.append(SimulationStep(
            step_id="generate_topology",
            title="Generar topología (OPLS-AA membrana)",
            stage=StepStage.PREPARATION,
            step_type=StepType.AUTOMATIC,
            engine="gromacs:pdb2gmx",
            depends_on=["embed_in_bilayer"],
            target_components=[prot_id],
            params={
                "source_file": "system.gro",
                "source_step": "embed_in_bilayer",   # input = embed_in_bilayer/system.gro
                "forcefield":  "opls-aa-membrane",   # → _FF_GROMACS_NAME → oplsaa_membrane
                "water_model": wm,
                "note": (
                    "Requires oplsaa_membrane.ff in working directory or GMXLIB path.\n"
                    "Copy Prot-Memb_FILES/oplsaa_membrane.ff to the step directory."
                ),
            },
        ))

        # ─────────────────────────────────────────────────────────────────────
        # Step 5: membrane embedding shrink loop (meta-step — fully automatic)
        # ─────────────────────────────────────────────────────────────────────
        plan.steps.append(SimulationStep(
            step_id="membrane_embedding",
            title=f"Shrink loop — convergencia APL ({lipid}, target ≤ {apl_target + mk.APL_CONVERGENCE_TOLERANCE:.0f} Å²)",
            stage=StepStage.MEMBRANE_EMBEDDING,
            step_type=StepType.AUTOMATIC,
            engine="gromacs+perl:inflategro",
            blocking=True,
            depends_on=["embed_in_bilayer", "generate_topology"],
            params={
                "lipid":               lipid,
                "lipid_residue_name":  lipid_resname,
                "forcefield":          ff,
                "temperature_K":       T,
                "apl_target_ang2":     apl_target,
                "apl_tolerance_ang2":  mk.APL_CONVERGENCE_TOLERANCE,
                "inflate_factor":      inflate_f,
                "deflate_factor":      mk.SHRINK_DEFLATION_FACTOR,
                "max_iterations":      mk.SHRINK_MAX_ITERATIONS,
                "gridsize":            5,
                "cutoff":              0.0,
                "inflategro_script":   "inflategro-Jorge.pl",
                "input_gro":           "../embed_in_bilayer/system.gro",
                "topol_top":           "../generate_topology/topol.top",
            },
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
            depends_on=["membrane_embedding"],
            params={
                "water_model": wm,
                "water_gro":   "spc216.gro",
                "input_gro":   "../membrane_embedding/converged.gro",
                "topol_top":   "../generate_topology/topol.top",
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
