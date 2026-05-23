# builders/step_builders/_utils.py

from __future__ import annotations

import os
from pathlib import Path


def rel(from_dir: Path, to_dir: Path) -> str:
    """Relative path from from_dir to to_dir for use in shell scripts."""
    return os.path.relpath(to_dir, from_dir)


def mdrun_block(deffnm: str, hardware: str = "auto", extra_gpu_flags: str = "", extra_cpu_flags: str = "") -> str:
    """
    Genera el bloque bash para gmx mdrun con soporte GPU/CPU/auto.

    hardware:
        "gpu"  — siempre usa flags GPU optimizados
        "cpu"  — siempre usa CPU
        "auto" — detecta nvidia-smi en runtime y elige automáticamente
    """
    gpu_cmd = (
        f"gmx mdrun -v -deffnm {deffnm}"
        f" -gpu_id 0 -pme gpu -bonded gpu -nb gpu -update cpu"
        f" -ntomp 10 -nstlist 150 -pin on -tunepme no -pmefft gpu -dlb auto"
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
