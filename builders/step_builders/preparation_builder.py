# builders/step_builders/preparation_builder.py

from __future__ import annotations

from pathlib import Path
import json

from core.execution_models import SimulationStep


class PreparationBuilder:
    """
    Genera instrucciones para el stage de preparación.

    Para proteínas: pdb2gmx
    Para ligandos: conversión a formato parametrizable

    No genera scripts ejecutables automáticamente porque
    la preparación requiere decisiones del usuario
    (protonación, forcefield, opciones de terminales).
    """

    def build(
        self,
        step:         SimulationStep,
        step_dir:     Path,
        step_dir_map: dict = {},
    ) -> None:

        targets = step.target_components
        engine  = step.engine

        if engine == "gromacs:pdb2gmx":
            self._build_protein_prep(step, step_dir)
        else:
            self._build_ligand_prep(step, step_dir)

    # ── Proteína ──────────────────────────────────────────────────────────────

    def _build_protein_prep(
        self,
        step:     SimulationStep,
        step_dir: Path,
    ) -> None:

        target = (
            step.target_components[0]
            if step.target_components
            else "protein"
        )

        commands = f"""#!/bin/bash
# ─── Preparación de proteína: {target} ───────────────────────────────────────
# Engine: gromacs:pdb2gmx
# Ejecutar manualmente — requiere selección interactiva de forcefield

# 1. Verificar archivo de entrada
#    El archivo PDB debe estar limpio (sin HETATM inesperados, sin cadenas rotas)

# 2. Generar topología con pdb2gmx
gmx pdb2gmx \\
    -f {target}.pdb \\
    -o {target}_processed.gro \\
    -p topol.top \\
    -ignh \\
    -ter

# Flags importantes:
#   -ignh     → ignorar hidrógenos existentes, agregar nuevos
#   -ter      → modo interactivo para terminales (N y C)
#
# Seleccionar en modo interactivo:
#   Forcefield → según configs/hmg_competition.yaml (charmm36)
#   Water      → tip3p

# 3. Verificar outputs
#   {target}_processed.gro  → estructura procesada
#   topol.top               → topología
#   posre.itp               → restraints de posición (generado automáticamente)
"""

        readme = f"""# Preparación: {target}

## Qué hace este step
Convierte el PDB de la proteína a formato GROMACS con topología completa.

## Engine
`gromacs:pdb2gmx`

## Notas
{chr(10).join(f'- {n}' for n in step.notes)}

## Outputs esperados
- `{target}_processed.gro`
- `topol.top`
- `posre.itp`

## Cómo ejecutar
```bash
bash commands.sh
```
"""

        (step_dir / "commands.sh").write_text(commands)
        (step_dir / "README.md").write_text(readme)
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":          step.step_id,
            "stage":            step.stage.value,
            "engine":           step.engine,
            "target":           target,
            "step_type":        step.step_type.value,
            "expected_outputs": [
                f"{target}_processed.gro",
                "topol.top",
                "posre.itp",
            ],
        }, indent=4))

    # ── Ligando ───────────────────────────────────────────────────────────────

    def _build_ligand_prep(
        self,
        step:     SimulationStep,
        step_dir: Path,
    ) -> None:

        target = (
            step.target_components[0]
            if step.target_components
            else "ligand"
        )

        commands = f"""#!/bin/bash
# ─── Preparación de ligando: {target} ────────────────────────────────────────
# Engine: ligand_preparation
# Ejecutar manualmente

# Opción A: si tienes SDF limpio
#   obabel {target}.pdb -O {target}.sdf --gen3d

# Opción B: si ya tienes SDF
#   cp {target}.sdf .

# Verificar estructura en Avogadro o PyMOL antes de parametrizar
echo "Verificar {target}.sdf antes de continuar a parametrización"
"""

        readme = f"""# Preparación: {target}

## Qué hace este step
Prepara el ligando para parametrización.
Convierte de PDB a SDF con conectividad explícita.

## Recomendación
Usar OpenBabel o RDKit para conversión limpia.
Verificar visualmente en Avogadro o PyMOL.

## Notas
{chr(10).join(f'- {n}' for n in step.notes) if step.notes else '- Sin notas adicionales'}

## Outputs esperados
- `{target}.sdf`
"""

        (step_dir / "commands.sh").write_text(commands)
        (step_dir / "README.md").write_text(readme)
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":          step.step_id,
            "stage":            step.stage.value,
            "engine":           step.engine,
            "target":           target,
            "step_type":        step.step_type.value,
            "expected_outputs": [f"{target}.sdf"],
        }, indent=4))
