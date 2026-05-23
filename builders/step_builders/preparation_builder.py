# builders/step_builders/preparation_builder.py

from __future__ import annotations

from pathlib import Path
import json

from core.execution_models import SimulationStep
from builders.step_builders._utils import rel as _rel

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

        workspace_root = step_dir_map.get("__workspace_root__")
        inputs_dir     = (Path(workspace_root) / "inputs") if workspace_root else None
        inputs_ref     = _rel(step_dir, inputs_dir) if inputs_dir else "../../inputs"

        engine = step.engine

        if engine == "gromacs:pdb2gmx":
            self._build_protein_prep(step, step_dir, inputs_ref)
        else:
            self._build_ligand_prep(step, step_dir, inputs_ref)

    # ── Proteína ──────────────────────────────────────────────────────────────

    def _build_protein_prep(
        self,
        step:       SimulationStep,
        step_dir:   Path,
        inputs_ref: str,
    ) -> None:

        target      = step.target_components[0] if step.target_components else "protein"
        source_file = step.params.get("source_file")
        forcefield  = step.params.get("forcefield", "charmm36")
        water_model = step.params.get("water_model", "tip3p")

        ff_gmx = _FF_GROMACS_NAME.get(forcefield, forcefield)

        # Derive the staged filename (component_id + original extension)
        src_ext  = Path(source_file).suffix if source_file else ".pdb"
        pdb_name = f"{target}{src_ext}"

        script = f"""#!/bin/bash
# ─── Preparación de proteína: {target} ───────────────────────────────────────
# El PDB de entrada vive en workspace/inputs/ — workspace auto-contenido.

INPUTS_DIR="{inputs_ref}"

gmx pdb2gmx \\
    -f "$INPUTS_DIR/{pdb_name}" \\
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
            "required_inputs":  [f"{inputs_ref}/{pdb_name}"],
        }, indent=4))

    # ── Ligando ───────────────────────────────────────────────────────────────

    def _build_ligand_prep(
        self,
        step:       SimulationStep,
        step_dir:   Path,
        inputs_ref: str,
    ) -> None:

        target = (
            step.target_components[0]
            if step.target_components
            else "ligand"
        )
        source_file = step.params.get("source_file")
        src_ext     = Path(source_file).suffix if source_file else ".pdb"
        pdb_name    = f"{target}{src_ext}"

        commands = f"""#!/bin/bash
# ─── Preparación de ligando: {target} ────────────────────────────────────────
# El PDB de entrada vive en workspace/inputs/ — workspace auto-contenido.
# Ejecutar manualmente.

INPUTS_DIR="{inputs_ref}"

# Opción A: convertir desde inputs (recomendado)
#   obabel "$INPUTS_DIR/{pdb_name}" -O {target}.sdf --gen3d

# Opción B: si ya tienes SDF listo
#   cp /path/to/{target}.sdf .

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
        meta: dict = {
            "step_id":          step.step_id,
            "stage":            step.stage.value,
            "engine":           step.engine,
            "target":           target,
            "step_type":        step.step_type.value,
            "expected_outputs": [f"{target}.sdf"],
        }
        if source_file:
            meta["required_inputs"] = [f"{inputs_ref}/{pdb_name}"]
        (step_dir / "metadata.json").write_text(json.dumps(meta, indent=4))
