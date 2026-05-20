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
        step: SimulationStep,
        step_dir: Path,
    ) -> None:

        # ────────────────────────────────────────────────────────────────────
        # em.mdp
        # ────────────────────────────────────────────────────────────────────

        mdp_text = """
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 50000

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
            "step_id": step.step_id,
            "stage": step.stage.value,
            "engine": step.engine,
            "generated_by": (
                "MinimizationBuilder"
            ),
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