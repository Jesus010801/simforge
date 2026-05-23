# builders/step_builders/production_builder.py

from __future__ import annotations

from pathlib import Path
import json

from core.execution_models import SimulationStep
from builders.step_builders._utils import rel as _rel, mdrun_block as _mdrun_block, mdrun_resume_block as _mdrun_resume_block


# ═══════════════════════════════════════════════════════════════════════════════
# Production Builder
# ═══════════════════════════════════════════════════════════════════════════════

class ProductionBuilder:
    """
    Genera archivos para dinámica molecular
    de producción en GROMACS.
    """

    def build(
        self,
        step:         SimulationStep,
        step_dir:     Path,
        step_dir_map: dict = {},
    ) -> None:

        p                   = step.params
        dt                  = p.get("dt",                   0.002)
        nsteps              = p.get("nsteps",               5_000_000)
        temperature         = p.get("temperature",          300.0)
        pressure            = p.get("pressure",             1.0)
        tc_grps             = p.get("tc_grps",              "Protein Non-Protein")
        tau_t               = p.get("tau_t",                " ".join(["0.1"] * len(tc_grps.split())))
        ref_t               = p.get("ref_t",                " ".join([str(temperature)] * len(tc_grps.split())))
        constraints         = p.get("constraints",          "h-bonds")
        nstxout_compressed  = p.get("nstxout_compressed",   5_000)
        nstenergy           = p.get("nstenergy",             1_000)
        nstlog              = p.get("nstlog",                1_000)
        # Extended params — membrane overrides
        tcoupl              = p.get("tcoupl",               "V-rescale")
        pcoupltype          = p.get("pcoupltype",           "isotropic")
        tau_p               = p.get("tau_p",                2.0)
        ref_p_xy            = p.get("ref_p_xy",             pressure)
        ref_p_z             = p.get("ref_p_z",             pressure)
        rcoulomb            = p.get("rcoulomb",              1.0)
        rvdw                = p.get("rvdw",                  1.0)
        disp_corr           = p.get("disp_corr",            "")
        hardware            = p.get("hardware",              "auto")

        # ── Inter-step paths ─────────────────────────────────────────────────
        # npt.gro / npt.cpt: output de equilibration (direct dep)
        eq_dir = next(
            (step_dir_map[d] for d in step.depends_on if "equilibration" in d and d in step_dir_map),
            None,
        )
        eq_ref = _rel(step_dir, eq_dir) if eq_dir else "../equilibration"

        # topol.top chain (final topology after ion addition):
        #   add_ions > generate_topology (membrane, no ion step) > assemble_system
        topol_dir = (
            step_dir_map.get("add_ions")
            or step_dir_map.get("generate_topology")
            or step_dir_map.get("assemble_system")
        )
        topol_ref = _rel(step_dir, topol_dir) if topol_dir else "../add_ions"

        # ────────────────────────────────────────────────────────────────────
        # md.mdp
        # ────────────────────────────────────────────────────────────────────

        if pcoupltype == "semiisotropic":
            pressure_block = (
                f"pcoupl                  = Parrinello-Rahman\n"
                f"pcoupltype              = semiisotropic\n"
                f"tau_p                   = {tau_p}\n"
                f"ref_p                   = {ref_p_xy}  {ref_p_z}\n"
                f"compressibility         = 4.5e-5  4.5e-5"
            )
        else:
            pressure_block = (
                f"pcoupl                  = Parrinello-Rahman\n"
                f"pcoupltype              = isotropic\n"
                f"tau_p                   = {tau_p}\n"
                f"ref_p                   = {pressure}\n"
                f"compressibility         = 4.5e-5"
            )

        disp_block = f"\nDispCorr                = {disp_corr}" if disp_corr else ""

        mdp_text = (
            f"title                   = Production MD\n\n"
            f"integrator              = md\n"
            f"dt                      = {dt}\n"
            f"nsteps                  = {nsteps}\n\n"
            f"tcoupl                  = {tcoupl}\n"
            f"tc-grps                 = {tc_grps}\n"
            f"tau_t                   = {tau_t}\n"
            f"ref_t                   = {ref_t}\n\n"
            f"{pressure_block}\n\n"
            f"constraints             = {constraints}\n\n"
            f"cutoff-scheme           = Verlet\n"
            f"coulombtype             = PME\n"
            f"rcoulomb                = {rcoulomb}\n"
            f"rvdw                    = {rvdw}\n\n"
            f"nstxout-compressed      = {nstxout_compressed}\n"
            f"nstenergy               = {nstenergy}\n"
            f"nstlog                  = {nstlog}\n\n"
            f"pbc                     = xyz{disp_block}\n"
        )

        mdp_path = (
            step_dir / "md.mdp"
        )

        mdp_path.write_text(
            mdp_text.strip()
        )

        # ────────────────────────────────────────────────────────────────────
        # Run script
        # ────────────────────────────────────────────────────────────────────

        run_script = f"""#!/bin/bash
# ─── Production MD ───────────────────────────────────────────────────────────
# Paths resueltos desde DAG

EQ_DIR="{eq_ref}"
TOPOL_DIR="{topol_ref}"

gmx grompp \\
    -f md.mdp \\
    -c "$EQ_DIR/npt.gro" \\
    -t "$EQ_DIR/npt.cpt" \\
    -p "$TOPOL_DIR/topol.top" \\
    -o md.tpr \\
    -maxwarn 1

{_mdrun_block("md", hardware, stage="md")}"""

        run_path = (
            step_dir / "run_md.sh"
        )

        run_path.write_text(
            run_script.strip()
        )

        # ────────────────────────────────────────────────────────────────────
        # Resume script (checkpoint recovery — no grompp, uses existing md.tpr)
        # ────────────────────────────────────────────────────────────────────

        resume_script = f"""#!/bin/bash
# ─── Production MD (resume from checkpoint) ──────────────────────────────────
# Resumes from md.cpt — md.tpr must already exist in this directory.
# Do NOT call grompp; the .tpr is reused as-is.

{_mdrun_resume_block("md", hardware, stage="md")}"""

        (step_dir / "run_md_resume.sh").write_text(resume_script.strip())

        # ────────────────────────────────────────────────────────────────────
        # Metadata
        # ────────────────────────────────────────────────────────────────────

        metadata = {
            "step_id":          step.step_id,
            "stage":            step.stage.value,
            "engine":           step.engine,
            "step_type":        step.step_type.value,
            "blocking":         step.blocking,
            "generated_by":     "ProductionBuilder",
            "simulation_type":  "production_md",
            "expected_outputs": ["md.xtc", "md.edr", "md.log", "md.gro", "md.cpt"],
            "required_inputs":  [f"{eq_ref}/npt.gro", f"{eq_ref}/npt.cpt"],
            "params": {
                "dt": dt, "nsteps": nsteps, "temperature": temperature,
                "pressure": pressure, "tc_grps": tc_grps, "constraints": constraints,
                "nstxout_compressed": nstxout_compressed,
            },
        }

        metadata_path = (
            step_dir / "metadata.json"
        )

        metadata_path.write_text(
            json.dumps(
                metadata,
                indent=4,
            )
        )