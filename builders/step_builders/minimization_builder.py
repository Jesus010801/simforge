# builders/step_builders/minimization_builder.py

from __future__ import annotations

from pathlib import Path
import json

from core.execution_models import (
    SimulationStep,
)


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
        # run.sh
        # ────────────────────────────────────────────────────────────────────

        run_script = """
gmx grompp \
    -f em.mdp \
    -c aaions.gro \
    -p topol.top \
    -o em.tpr \
    -maxwarn 1

gmx mdrun \
    -v \
    -deffnm em \
    -nb gpu
"""

        run_path = step_dir / "run.sh"

        run_path.write_text(
            run_script.strip()
        )

        # ────────────────────────────────────────────────────────────────────
        # metadata.json
        # ────────────────────────────────────────────────────────────────────

        metadata = {
            "step_id":      step.step_id,
            "stage":        step.stage.value,
            "engine":       step.engine,
            "step_type":    step.step_type.value,
            "blocking":     step.blocking,
            "generated_by": "MinimizationBuilder",
            "params":       {"integrator": integrator, "emtol": emtol, "emstep": emstep, "nsteps": nsteps},
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