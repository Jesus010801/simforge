# builders/step_builders/equilibration_builder.py

from __future__ import annotations

from pathlib import Path
import json

from core.execution_models import (
    SimulationStep,
)


# ═══════════════════════════════════════════════════════════════════════════════
# Equilibration Builder
# ═══════════════════════════════════════════════════════════════════════════════

class EquilibrationBuilder:
    """
    Genera archivos para equilibrio NVT/NPT.
    """

    def build(
        self,
        step: SimulationStep,
        step_dir: Path,
    ) -> None:

        # ────────────────────────────────────────────────────────────────────
        # NVT MDP
        # ────────────────────────────────────────────────────────────────────

        nvt_mdp = """
title                   = NVT equilibration

integrator              = md
dt                      = 0.002
nsteps                  = 50000

tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1 0.1
ref_t                   = 300 300

pcoupl                  = no

constraints             = h-bonds

cutoff-scheme           = Verlet
coulombtype             = PME
rcoulomb                = 1.0
rvdw                    = 1.0

pbc                     = xyz
"""

        nvt_path = (
            step_dir / "nvt.mdp"
        )

        nvt_path.write_text(
            nvt_mdp.strip()
        )

        # ────────────────────────────────────────────────────────────────────
        # NPT MDP
        # ────────────────────────────────────────────────────────────────────

        npt_mdp = """
title                   = NPT equilibration

integrator              = md
dt                      = 0.002
nsteps                  = 50000

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

pbc                     = xyz
"""

        npt_path = (
            step_dir / "npt.mdp"
        )

        npt_path.write_text(
            npt_mdp.strip()
        )

        # ────────────────────────────────────────────────────────────────────
        # Run scripts
        # ────────────────────────────────────────────────────────────────────

        nvt_script = """
gmx grompp \
    -f nvt.mdp \
    -c em.gro \
    -r em.gro \
    -p topol.top \
    -o nvt.tpr

gmx mdrun \
    -v \
    -deffnm nvt \
    -nb gpu
"""

        nvt_run = (
            step_dir / "run_nvt.sh"
        )

        nvt_run.write_text(
            nvt_script.strip()
        )

        npt_script = """
gmx grompp \
    -f npt.mdp \
    -c nvt.gro \
    -r nvt.gro \
    -t nvt.cpt \
    -p topol.top \
    -o npt.tpr

gmx mdrun \
    -v \
    -deffnm npt \
    -nb gpu
"""

        npt_run = (
            step_dir / "run_npt.sh"
        )

        npt_run.write_text(
            npt_script.strip()
        )

        # ────────────────────────────────────────────────────────────────────
        # Metadata
        # ────────────────────────────────────────────────────────────────────

        metadata = {
            "step_id": step.step_id,
            "stage": step.stage.value,
            "engine": step.engine,
            "equilibration_type": [
                "NVT",
                "NPT",
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