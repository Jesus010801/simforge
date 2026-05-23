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
    Genera archivos necesarios para
    minimización energética en GROMACS.
    """

    def build(
        self,
        step:         SimulationStep,
        step_dir:     Path,
        step_dir_map: dict = {},
    ) -> None:

        # ────────────────────────────────────────────────────────────────────
        # em.mdp
        # ────────────────────────────────────────────────────────────────────

        p = step.params
        integrator = p.get("integrator", "steep")
        emtol      = p.get("emtol",      1000.0)
        emstep     = p.get("emstep",     0.01)
        nsteps     = p.get("nsteps",     50_000)
        hardware   = p.get("hardware",   "auto")

        mdp_text = f"""
integrator  = {integrator}
emtol       = {emtol}
emstep      = {emstep}
nsteps      = {nsteps}

cutoff-scheme = Verlet

nstlist     = 10
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0

pbc         = xyz
"""

        mdp_path = step_dir / "em.mdp"

        mdp_path.write_text(
            mdp_text.strip()
        )

        # ────────────────────────────────────────────────────────────────────
        # Inter-step paths (resolved from DAG)
        # ────────────────────────────────────────────────────────────────────

        ions_dir = next(
            (step_dir_map[d] for d in step.depends_on if "ions" in d and d in step_dir_map),
            None,
        )
        ions_ref = _rel(step_dir, ions_dir) if ions_dir else "../add_ions"

        assemble_dir = step_dir_map.get("assemble_system")
        topol_ref    = _rel(step_dir, assemble_dir) if assemble_dir else "../assemble_system"

        # ────────────────────────────────────────────────────────────────────
        # run.sh
        # ────────────────────────────────────────────────────────────────────

        run_script = f"""#!/bin/bash
# ─── Energy minimization ─────────────────────────────────────────────────────
# Paths resueltos desde DAG

IONS_DIR="{ions_ref}"
TOPOL_DIR="{topol_ref}"

gmx grompp \\
    -f em.mdp \\
    -c "$IONS_DIR/aaions.gro" \\
    -p "$TOPOL_DIR/topol.top" \\
    -o em.tpr \\
    -maxwarn 1

{_mdrun_block("em", hardware)}"""

        run_path = step_dir / "run.sh"

        run_path.write_text(
            run_script.strip()
        )

        # ────────────────────────────────────────────────────────────────────
        # metadata.json
        # ────────────────────────────────────────────────────────────────────

        metadata = {
            "step_id":          step.step_id,
            "stage":            step.stage.value,
            "engine":           step.engine,
            "step_type":        step.step_type.value,
            "blocking":         step.blocking,
            "generated_by":     "MinimizationBuilder",
            "expected_outputs": ["em.gro", "em.edr", "em.log"],
            "params":           {"integrator": integrator, "emtol": emtol, "emstep": emstep, "nsteps": nsteps},
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