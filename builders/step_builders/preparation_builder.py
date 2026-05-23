# builders/step_builders/preparation_builder.py

from __future__ import annotations

from pathlib import Path
import json
import shutil

from core.execution_models import SimulationStep

# GROMACS uses "oplsaa" (no dash) as the -ff argument; the ontology uses "opls-aa"
_FF_GROMACS_NAME: dict[str, str] = {
    "opls-aa":          "oplsaa",
    "opls-aa-membrane": "oplsaa_membrane",
}


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

        target      = step.target_components[0] if step.target_components else "protein"
        source_file = step.params.get("source_file")
        forcefield  = step.params.get("forcefield", "charmm36")
        water_model = step.params.get("water_model", "tip3p")

        # Translate ontology FF name → GROMACS -ff argument
        ff_gmx = _FF_GROMACS_NAME.get(forcefield, forcefield)

        # Copy source PDB into the step directory at compile time
        if source_file:
            src = Path(source_file)
            dst = step_dir / f"{target}.pdb"
            if src.exists() and src.resolve() != dst.resolve():
                shutil.copy2(src, dst)

        script = f"""#!/bin/bash
# ─── Preparación de proteína: {target} ───────────────────────────────────────

gmx pdb2gmx \\
    -f {target}.pdb \\
    -o {target}_processed.gro \\
    -p topol.top \\
    -ff {ff_gmx} \\
    -water {water_model} \\
    -ignh
"""

        (step_dir / "run.sh").write_text(script)
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
