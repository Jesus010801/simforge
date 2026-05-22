# builders/step_builders/equilibration_builder.py

from __future__ import annotations

from pathlib import Path
import json

from core.execution_models import SimulationStep
from builders.step_builders._utils import rel as _rel


# ═══════════════════════════════════════════════════════════════════════════════
# Equilibration Builder
# ═══════════════════════════════════════════════════════════════════════════════

class EquilibrationBuilder:
    """
    Genera archivos para equilibrio NVT/NPT.
    """

    def build(
        self,
        step:         SimulationStep,
        step_dir:     Path,
        step_dir_map: dict = {},
    ) -> None:

        p           = step.params
        dt          = p.get("dt",          0.002)
        nvt_nsteps  = p.get("nvt_nsteps",  50_000)
        npt_nsteps  = p.get("npt_nsteps",  50_000)
        temperature = p.get("temperature", 300.0)
        pressure    = p.get("pressure",    1.0)
        tc_grps     = p.get("tc_grps",     "Protein Non-Protein")
        tau_t       = p.get("tau_t",       " ".join(["0.1"] * len(tc_grps.split())))
        ref_t       = p.get("ref_t",       " ".join([str(temperature)] * len(tc_grps.split())))
        constraints = p.get("constraints", "h-bonds")

        # ── Inter-step paths ─────────────────────────────────────────────────
        # em.gro: output de energy_minimization (direct dep)
        em_dir = next(
            (step_dir_map[d] for d in step.depends_on if "minimization" in d and d in step_dir_map),
            None,
        )
        em_ref = _rel(step_dir, em_dir) if em_dir else "../energy_minimization"

        # topol.top: vive en assemble_system (no es dep directo, pero siempre presente)
        assemble_dir = step_dir_map.get("assemble_system")
        topol_ref    = _rel(step_dir, assemble_dir) if assemble_dir else "../assemble_system"

        # ────────────────────────────────────────────────────────────────────
        # NVT MDP
        # ────────────────────────────────────────────────────────────────────

        nvt_mdp = f"""
title                   = NVT equilibration

integrator              = md
dt                      = {dt}
nsteps                  = {nvt_nsteps}

tcoupl                  = V-rescale
tc-grps                 = {tc_grps}
tau_t                   = {tau_t}
ref_t                   = {ref_t}

pcoupl                  = no

constraints             = {constraints}

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

        npt_mdp = f"""
title                   = NPT equilibration

integrator              = md
dt                      = {dt}
nsteps                  = {npt_nsteps}

tcoupl                  = V-rescale
tc-grps                 = {tc_grps}
tau_t                   = {tau_t}
ref_t                   = {ref_t}

pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = {pressure}
compressibility         = 4.5e-5

constraints             = {constraints}

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

        nvt_script = f"""#!/bin/bash
# ─── NVT equilibration ───────────────────────────────────────────────────────
# Paths resueltos desde DAG

EM_DIR="{em_ref}"
TOPOL_DIR="{topol_ref}"

gmx grompp \\
    -f nvt.mdp \\
    -c "$EM_DIR/em.gro" \\
    -r "$EM_DIR/em.gro" \\
    -p "$TOPOL_DIR/topol.top" \\
    -o nvt.tpr

gmx mdrun \\
    -v \\
    -deffnm nvt \\
    -nb gpu
"""

        nvt_run = (
            step_dir / "run_nvt.sh"
        )

        nvt_run.write_text(
            nvt_script.strip()
        )

        npt_script = f"""#!/bin/bash
# ─── NPT equilibration ───────────────────────────────────────────────────────
# nvt.gro y nvt.cpt son outputs locales del step NVT anterior

TOPOL_DIR="{topol_ref}"

gmx grompp \\
    -f npt.mdp \\
    -c nvt.gro \\
    -r nvt.gro \\
    -t nvt.cpt \\
    -p "$TOPOL_DIR/topol.top" \\
    -o npt.tpr

gmx mdrun \\
    -v \\
    -deffnm npt \\
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
            "step_id":           step.step_id,
            "stage":             step.stage.value,
            "engine":            step.engine,
            "step_type":         step.step_type.value,
            "blocking":          step.blocking,
            "generated_by":      "EquilibrationBuilder",
            "equilibration_type": ["NVT", "NPT"],
            "params": {
                "dt": dt, "nvt_nsteps": nvt_nsteps, "npt_nsteps": npt_nsteps,
                "temperature": temperature, "pressure": pressure,
                "tc_grps": tc_grps, "constraints": constraints,
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