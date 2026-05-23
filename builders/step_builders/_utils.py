# builders/step_builders/_utils.py

from __future__ import annotations

import os
from pathlib import Path


def rel(from_dir: Path, to_dir: Path) -> str:
    """Relative path from from_dir to to_dir for use in shell scripts."""
    return os.path.relpath(to_dir, from_dir)


def mdrun_block(
    deffnm: str,
    hardware: str = "auto",
    stage: str = "md",
    extra_gpu_flags: str = "",
    extra_cpu_flags: str = "",
) -> str:
    """
    Genera el bloque bash para gmx mdrun con soporte GPU/CPU/auto.

    hardware:
        "gpu"  — siempre usa flags GPU compatibles con el stage
        "cpu"  — siempre usa CPU
        "auto" — detecta nvidia-smi en runtime y elige automáticamente

    stage:
        "minimization" — non-dynamical integrators (steep, cg, l-bfgs).
                         GROMACS >=2024 no soporta -pme gpu con estos integrators.
                         Solo -nb gpu es válido.
        "md"           — dynamical integrators (md, md-vv).
                         Soporta todos los flags GPU (-pme gpu, -bonded gpu, etc.)
    """
    if stage == "minimization":
        # Non-dynamical integrators: only neighbor-list on GPU is supported.
        # -pme gpu, -bonded gpu, -pmefft gpu require a dynamical integrator.
        gpu_cmd = (
            f"gmx mdrun -v -deffnm {deffnm}"
            f" -gpu_id 0 -nb gpu -pme cpu"
            f" -ntmpi 1 -ntomp $(nproc) -pin on"
            + (f" {extra_gpu_flags}" if extra_gpu_flags else "")
        )
    else:
        # Dynamical integrators (md, md-vv): full GPU offloading.
        gpu_cmd = (
            f"gmx mdrun -v -deffnm {deffnm}"
            f" -gpu_id 0 -pme gpu -bonded gpu -nb gpu -update cpu"
            f" -ntmpi 1 -ntomp $(nproc) -nstlist 150 -pin on -tunepme no -pmefft gpu -dlb auto"
            + (f" {extra_gpu_flags}" if extra_gpu_flags else "")
        )

    cpu_cmd = (
        f"gmx mdrun -v -deffnm {deffnm}"
        f" -nb cpu -pme cpu -ntmpi 1 -ntomp $(nproc) -pin on"
        + (f" {extra_cpu_flags}" if extra_cpu_flags else "")
    )

    if hardware == "gpu":
        return gpu_cmd + "\n"
    if hardware == "cpu":
        return cpu_cmd + "\n"

    # auto: detección en runtime
    return f"""if command -v nvidia-smi &>/dev/null && nvidia-smi &>/dev/null 2>&1; then
    {gpu_cmd}
else
    {cpu_cmd}
fi
"""


def mdrun_resume_block(
    deffnm: str,
    hardware: str = "auto",
    stage: str = "md",
) -> str:
    """
    Generates a bash block for `gmx mdrun` that resumes from a checkpoint.

    Equivalent to mdrun_block() but adds `-cpi {deffnm}.cpt -append`.
    No grompp call — md.tpr must already exist in the working directory.
    """
    return mdrun_block(
        deffnm,
        hardware=hardware,
        stage=stage,
        extra_gpu_flags=f"-cpi {deffnm}.cpt -append",
        extra_cpu_flags=f"-cpi {deffnm}.cpt -append",
    )
