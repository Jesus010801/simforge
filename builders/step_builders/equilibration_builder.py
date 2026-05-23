# builders/step_builders/equilibration_builder.py

from __future__ import annotations

from pathlib import Path
import json

from core.execution_models import SimulationStep
from builders.step_builders._utils import rel as _rel, mdrun_block as _mdrun_block


class EquilibrationBuilder:
    """
    Genera archivos para equilibrio NVT/NPT.

    Soporta tanto sistemas proteína-en-agua (isotropic, Parrinello-Rahman)
    como sistemas de membrana (semiisotropic, Berendsen NPT) — controlado
    completamente por step.params sin lógica condicional en el builder.
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
        hardware    = p.get("hardware",    "auto")
        # Pressure coupling — membrane uses semiisotropic + Berendsen
        pcoupltype  = p.get("pcoupltype",  "isotropic")
        pcoupl_npt  = p.get("pcoupl_npt",  "Parrinello-Rahman")
        ref_p_xy    = p.get("ref_p_xy",    pressure)
        ref_p_z     = p.get("ref_p_z",     pressure)
        tau_p       = p.get("tau_p",       2.0)
        # Cutoffs — membrane uses 1.2 nm; protein-in-water uses 1.0 nm
        rcoulomb    = p.get("rcoulomb",    1.0)
        rvdw        = p.get("rvdw",        1.0)
        disp_corr   = p.get("disp_corr",   "")

        # Pressure block differs for isotropic vs semiisotropic
        if pcoupltype == "semiisotropic":
            npt_pressure_block = (
                f"pcoupl                  = {pcoupl_npt}\n"
                f"pcoupltype              = semiisotropic\n"
                f"tau_p                   = {tau_p}\n"
                f"ref_p                   = {ref_p_xy}  {ref_p_z}\n"
                f"compressibility         = 4.5e-5  4.5e-5"
            )
        else:
            npt_pressure_block = (
                f"pcoupl                  = {pcoupl_npt}\n"
                f"pcoupltype              = isotropic\n"
                f"tau_p                   = {tau_p}\n"
                f"ref_p                   = {pressure}\n"
                f"compressibility         = 4.5e-5"
            )

        disp_block = f"\nDispCorr                = {disp_corr}" if disp_corr else ""

        # ── Inter-step paths ─────────────────────────────────────────────────
        em_dir = next(
            (step_dir_map[d] for d in step.depends_on if "minimization" in d and d in step_dir_map),
            None,
        )
        em_ref = _rel(step_dir, em_dir) if em_dir else "../energy_minimization"

        # topol.top — membrane: generate_topology; protein: assemble_system
        topol_dir = (
            step_dir_map.get("generate_topology")
            or step_dir_map.get("assemble_system")
        )
        topol_ref = _rel(step_dir, topol_dir) if topol_dir else "../assemble_system"

        # ── NVT MDP ──────────────────────────────────────────────────────────

        nvt_mdp = (
            f"title                   = NVT equilibration\n\n"
            f"integrator              = md\n"
            f"dt                      = {dt}\n"
            f"nsteps                  = {nvt_nsteps}\n\n"
            f"tcoupl                  = V-rescale\n"
            f"tc-grps                 = {tc_grps}\n"
            f"tau_t                   = {tau_t}\n"
            f"ref_t                   = {ref_t}\n\n"
            f"pcoupl                  = no\n\n"
            f"gen_vel                 = yes\n"
            f"gen_temp                = {temperature}\n"
            f"gen_seed                = -1\n\n"
            f"constraints             = {constraints}\n\n"
            f"cutoff-scheme           = Verlet\n"
            f"coulombtype             = PME\n"
            f"rcoulomb                = {rcoulomb}\n"
            f"rvdw                    = {rvdw}\n\n"
            f"pbc                     = xyz{disp_block}\n"
        )
        (step_dir / "nvt.mdp").write_text(nvt_mdp)

        # ── NPT MDP ──────────────────────────────────────────────────────────

        npt_mdp = (
            f"title                   = NPT equilibration\n\n"
            f"integrator              = md\n"
            f"dt                      = {dt}\n"
            f"nsteps                  = {npt_nsteps}\n\n"
            f"continuation            = yes\n\n"
            f"tcoupl                  = V-rescale\n"
            f"tc-grps                 = {tc_grps}\n"
            f"tau_t                   = {tau_t}\n"
            f"ref_t                   = {ref_t}\n\n"
            f"{npt_pressure_block}\n\n"
            f"constraints             = {constraints}\n\n"
            f"cutoff-scheme           = Verlet\n"
            f"coulombtype             = PME\n"
            f"rcoulomb                = {rcoulomb}\n"
            f"rvdw                    = {rvdw}\n\n"
            f"pbc                     = xyz{disp_block}\n"
        )
        (step_dir / "npt.mdp").write_text(npt_mdp)

        # ── Run scripts ───────────────────────────────────────────────────────

        (step_dir / "run_nvt.sh").write_text(
            f"#!/bin/bash\n"
            f"# ─── NVT equilibration ───────────────────────────────────────────────────────\n"
            f"EM_DIR=\"{em_ref}\"\n"
            f"TOPOL_DIR=\"{topol_ref}\"\n\n"
            f"gmx grompp \\\n"
            f"    -f nvt.mdp \\\n"
            f"    -c \"$EM_DIR/em.gro\" \\\n"
            f"    -r \"$EM_DIR/em.gro\" \\\n"
            f"    -p \"$TOPOL_DIR/topol.top\" \\\n"
            f"    -o nvt.tpr \\\n"
            f"    -maxwarn 1\n\n"
            f"{_mdrun_block('nvt', hardware)}"
        )

        (step_dir / "run_npt.sh").write_text(
            f"#!/bin/bash\n"
            f"# ─── NPT equilibration ───────────────────────────────────────────────────────\n"
            f"TOPOL_DIR=\"{topol_ref}\"\n\n"
            f"gmx grompp \\\n"
            f"    -f npt.mdp \\\n"
            f"    -c nvt.gro \\\n"
            f"    -r nvt.gro \\\n"
            f"    -t nvt.cpt \\\n"
            f"    -p \"$TOPOL_DIR/topol.top\" \\\n"
            f"    -o npt.tpr \\\n"
            f"    -maxwarn 1\n\n"
            f"{_mdrun_block('npt', hardware)}"
        )

        (step_dir / "run.sh").write_text(
            "#!/bin/bash\nset -e\n"
            "bash \"$(dirname \"$0\")/run_nvt.sh\"\n"
            "bash \"$(dirname \"$0\")/run_npt.sh\"\n"
        )

        # ── Metadata ──────────────────────────────────────────────────────────

        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":            step.step_id,
            "stage":              step.stage.value,
            "engine":             step.engine,
            "step_type":          step.step_type.value,
            "blocking":           step.blocking,
            "generated_by":       "EquilibrationBuilder",
            "equilibration_type": ["NVT", "NPT"],
            "expected_outputs":   ["npt.gro", "npt.cpt"],
            "params": {
                "dt": dt, "nvt_nsteps": nvt_nsteps, "npt_nsteps": npt_nsteps,
                "temperature": temperature, "pcoupltype": pcoupltype,
                "pcoupl_npt": pcoupl_npt, "tc_grps": tc_grps,
            },
        }, indent=4))
