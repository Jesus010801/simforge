# builders/step_builders/production_builder.py

from __future__ import annotations

from pathlib import Path
import json

from core.execution_models import SimulationStep
from builders.step_builders._utils import rel as _rel


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

        # ── Inter-step paths ─────────────────────────────────────────────────
        # npt.gro / npt.cpt: output de equilibration (direct dep)
        eq_dir = next(
            (step_dir_map[d] for d in step.depends_on if "equilibration" in d and d in step_dir_map),
            None,
        )
        eq_ref = _rel(step_dir, eq_dir) if eq_dir else "../equilibration"

        # topol.top: vive en assemble_system (no es dep directo, pero siempre presente)
        assemble_dir = step_dir_map.get("assemble_system")
        topol_ref    = _rel(step_dir, assemble_dir) if assemble_dir else "../assemble_system"

        # ────────────────────────────────────────────────────────────────────
        # md.mdp
        # ────────────────────────────────────────────────────────────────────

        mdp_text = f"""
title                   = Production MD

integrator              = md
dt                      = {dt}
nsteps                  = {nsteps}

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

nstxout-compressed      = {nstxout_compressed}
nstenergy               = {nstenergy}
nstlog                  = {nstlog}

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
    -o md.tpr

gmx mdrun \\
    -v \\
    -deffnm md \\
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
            "step_id":          step.step_id,
            "stage":            step.stage.value,
            "engine":           step.engine,
            "step_type":        step.step_type.value,
            "blocking":         step.blocking,
            "generated_by":     "ProductionBuilder",
            "simulation_type":  "production_md",
            "expected_outputs": ["md.xtc", "md.edr", "md.log", "md.gro", "md.cpt"],
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