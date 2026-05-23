# builders/step_builders/minimization_builder.py

from __future__ import annotations

from pathlib import Path
import json

from core.execution_models import SimulationStep
from builders.step_builders._utils import rel as _rel, mdrun_block as _mdrun_block


# ═══════════════════════════════════════════════════════════════════════════════
# Minimization Builder
# ═══════════════════════════════════════════════════════════════════════════════

class MinimizationBuilder:
    """
    Genera archivos para minimización energética en GROMACS.

    Soporta sistemas proteína-en-agua y sistemas de membrana:
        - define: "-DPOSRES -DSTRONG_POSRES" para EM con restraints de membrana
        - rcoulomb/rvdw: 1.2 nm para membrana (vs 1.0 nm estándar)
        - disp_corr: "EnerPres" para sistemas de membrana
        - topol_ref: prefiere generate_topology > assemble_system
    """

    def build(
        self,
        step:         SimulationStep,
        step_dir:     Path,
        step_dir_map: dict = {},
    ) -> None:

        p          = step.params
        integrator = p.get("integrator", "steep")
        emtol      = p.get("emtol",      1000.0)
        emstep     = p.get("emstep",     0.01)
        nsteps     = p.get("nsteps",     50_000)
        hardware   = p.get("hardware",   "auto")
        define     = p.get("define",     "")
        rcoulomb   = p.get("rcoulomb",   1.0)
        rvdw       = p.get("rvdw",       1.0)
        disp_corr  = p.get("disp_corr",  "")

        define_block = f"define                  = {define}\n\n" if define else ""
        disp_block   = f"\nDispCorr                = {disp_corr}" if disp_corr else ""

        mdp_text = (
            f"{define_block}"
            f"integrator              = {integrator}\n"
            f"emtol                   = {emtol}\n"
            f"emstep                  = {emstep}\n"
            f"nsteps                  = {nsteps}\n\n"
            f"cutoff-scheme           = Verlet\n\n"
            f"nstlist                 = 10\n"
            f"coulombtype             = PME\n"
            f"rcoulomb                = {rcoulomb}\n"
            f"rvdw                    = {rvdw}\n\n"
            f"pbc                     = xyz{disp_block}\n"
        )

        (step_dir / "em.mdp").write_text(mdp_text)

        # ── Inter-step paths ─────────────────────────────────────────────────
        ions_dir = next(
            (step_dir_map[d] for d in step.depends_on if "ions" in d and d in step_dir_map),
            None,
        )
        ions_ref = _rel(step_dir, ions_dir) if ions_dir else "../add_ions"

        # topol.top chain (final topology after ion addition):
        #   add_ions > generate_topology (membrane, no ion step) > assemble_system
        topol_dir = (
            step_dir_map.get("add_ions")
            or step_dir_map.get("generate_topology")
            or step_dir_map.get("assemble_system")
        )
        topol_ref = _rel(step_dir, topol_dir) if topol_dir else "../add_ions"

        run_script = f"""#!/bin/bash
# ─── Energy minimization ─────────────────────────────────────────────────────
IONS_DIR="{ions_ref}"
TOPOL_DIR="{topol_ref}"

gmx grompp \\
    -f em.mdp \\
    -c "$IONS_DIR/aaions.gro" \\
    -p "$TOPOL_DIR/topol.top" \\
    -o em.tpr \\
    -maxwarn 1

{_mdrun_block("em", hardware, stage="minimization")}"""

        (step_dir / "run.sh").write_text(run_script.strip())
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":          step.step_id,
            "stage":            step.stage.value,
            "engine":           step.engine,
            "step_type":        step.step_type.value,
            "blocking":         step.blocking,
            "generated_by":     "MinimizationBuilder",
            "expected_outputs": ["em.gro", "em.edr", "em.log"],
            "required_inputs":  [f"{ions_ref}/aaions.gro", f"{topol_ref}/topol.top"],
            "params": {
                "integrator": integrator, "emtol": emtol,
                "emstep": emstep, "nsteps": nsteps,
                "define": define, "rcoulomb": rcoulomb, "rvdw": rvdw,
            },
        }, indent=4))