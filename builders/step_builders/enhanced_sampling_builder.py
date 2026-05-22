# builders/step_builders/enhanced_sampling_builder.py

from __future__ import annotations

from pathlib import Path
import json

from core.execution_models import SimulationStep


class EnhancedSamplingBuilder:
    """
    Genera instrucciones para enhanced sampling.

    REST2: Replica Exchange with Solute Tempering
    Metadinámica: via PLUMED

    Estos steps son externos — requieren PLUMED compilado
    con GROMACS y configuración específica del sistema.
    """

    def build(
        self,
        step:         SimulationStep,
        step_dir:     Path,
        step_dir_map: dict = {},
    ) -> None:

        sid = step.step_id

        if "rest2" in sid:
            self._build_rest2(step, step_dir)
        elif "metadyn" in sid or "metad" in sid:
            self._build_metadynamics(step, step_dir)
        else:
            self._build_generic(step, step_dir)

    # ── REST2 ─────────────────────────────────────────────────────────────────

    def _build_rest2(
        self,
        step:     SimulationStep,
        step_dir: Path,
    ) -> None:

        readme = f"""# REST2 Enhanced Sampling

## Engine
`plumed:gromacs` — Replica Exchange with Solute Tempering (REST2)

## Notas
{chr(10).join(f'- {n}' for n in step.notes)}

## Requisitos
- GROMACS compilado con soporte MPI
- PLUMED instalado y patcheado en GROMACS

## Configuración básica REST2
```bash
# 1. Definir réplicas (temperatura efectiva del soluto)
#    Típico: 4-8 réplicas entre 300K y 450K

# 2. Preparar MDP para cada réplica
#    Ver rest2_template.mdp

# 3. Configurar PLUMED
#    Ver plumed.dat

# 4. Lanzar
mpirun -np 4 gmx_mpi mdrun \\
    -v \\
    -deffnm rest2 \\
    -multidir rep0 rep1 rep2 rep3 \\
    -replex 1000 \\
    -hrex \\
    -plumed plumed.dat
```

## Referencias
- REST2: Terakawa et al. JCTC 2011
- PLUMED: https://www.plumed.org
"""

        plumed_template = """# plumed.dat — template REST2
# Adaptar según el sistema

# Definir grupo del soluto (ligando)
MOLINFO STRUCTURE=complex.pdb

lig: GROUP ATOMS=1-35  # ajustar índices del ligando

# REST2: escalar interacciones del soluto
REST2 ...
  SOLUTE=lig
  TEMP=300
... REST2

PRINT ARG=* FILE=colvar STRIDE=1000
"""

        (step_dir / "README.md").write_text(readme)
        (step_dir / "plumed_template.dat").write_text(plumed_template)
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":   step.step_id,
            "stage":     step.stage.value,
            "engine":    step.engine,
            "blocking":  step.blocking,
            "step_type": step.step_type.value,
            "notes":     step.notes,
            "expected_outputs": [
                "rest2.xtc",
                "rest2.edr",
                "rest2.log",
            ],
        }, indent=4))

    # ── Metadinámica ──────────────────────────────────────────────────────────

    def _build_metadynamics(
        self,
        step:     SimulationStep,
        step_dir: Path,
    ) -> None:

        (step_dir / "README.md").write_text(
            "# Metadinámica\n\n"
            "Requiere PLUMED + definición de collective variables (CVs).\n\n"
            "Ver documentación: https://www.plumed.org/doc-v2.9/user-doc/html/\n"
        )
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":   step.step_id,
            "stage":     step.stage.value,
            "engine":    step.engine,
            "step_type": step.step_type.value,
        }, indent=4))

    # ── Genérico ──────────────────────────────────────────────────────────────

    def _build_generic(
        self,
        step:     SimulationStep,
        step_dir: Path,
    ) -> None:

        (step_dir / "README.md").write_text(
            f"# Enhanced Sampling: {step.step_id}\n\n"
            f"Engine: {step.engine}\n\n"
            f"Step externo — configurar manualmente.\n"
        )
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":   step.step_id,
            "stage":     step.stage.value,
            "engine":    step.engine,
            "step_type": step.step_type.value,
        }, indent=4))