# builders/step_builders/production_builder.py

from __future__ import annotations

from pathlib import Path
import json

from core.execution_models import (
    SimulationStep,
)


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
        step: SimulationStep,
        step_dir: Path,
    ) -> None:

        # ────────────────────────────────────────────────────────────────────
        # md.mdp
        # ────────────────────────────────────────────────────────────────────

        mdp_text = """
title                   = Production MD

integrator              = md
dt                      = 0.002
nsteps                  = 5000000

tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1 0.1
ref_t                   = 300 300

pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5

constraints             = h-bonds

cutoff-scheme           = Verlet
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0

nstxout-compressed      = 5000
nstenergy               = 1000
nstlog                  = 1000

pbc                     = xyz
"""

        mdp_path = (
            step_dir / "md.mdp"
        )

        mdp_path.write_text(
            mdp_text.strip()
        )

        # ────────────────────────────────────────────────────────────────────
        # Run script
        # ────────────────────────────────────────────────────────────────────

        run_script = """
gmx grompp \
    -f md.mdp \
    -c npt.gro \
    -t npt.cpt \
    -p topol.top \
    -o md.tpr

gmx mdrun \
    -v \
    -deffnm md \
    -nb gpu
"""

        run_path = (
            step_dir / "run_md.sh"
        )

        run_path.write_text(
            run_script.strip()
        )

        # ────────────────────────────────────────────────────────────────────
        # Metadata
        # ────────────────────────────────────────────────────────────────────

        metadata = {
            "step_id": step.step_id,
            "stage": step.stage.value,
            "engine": step.engine,
            "simulation_type": (
                "production_md"
            ),
            "expected_outputs": [
                "md.xtc",
                "md.edr",
                "md.log",
                "md.gro",
                "md.cpt",
            ],
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