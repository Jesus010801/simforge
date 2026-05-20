# builders/step_builders/parametrization_builder.py

from __future__ import annotations

from pathlib import Path
import json

from core.execution_models import SimulationStep


class ParametrizationBuilder:
    """
    Genera instrucciones para parametrización de ligandos.

    CGenFF: via ParamChem (web) o CHARMM
    GAFF:   via antechamber (AmberTools)
    OpenFF: via openff-toolkit
    """

    def build(
        self,
        step:     SimulationStep,
        step_dir: Path,
    ) -> None:

        target = (
            step.target_components[0]
            if step.target_components
            else "ligand"
        )

        engine = step.engine.lower()

        if "cgenff" in engine:
            self._build_cgenff(step, step_dir, target)
        elif "gaff" in engine:
            self._build_gaff(step, step_dir, target)
        else:
            self._build_generic(step, step_dir, target)

    # ── CGenFF ────────────────────────────────────────────────────────────────

    def _build_cgenff(
        self,
        step:     SimulationStep,
        step_dir: Path,
        target:   str,
    ) -> None:

        is_blocking = step.blocking

        commands = f"""#!/bin/bash
# ─── Parametrización CGenFF: {target} ────────────────────────────────────────
{'# ⚠  BLOCKING: revisión manual requerida antes de continuar' if is_blocking else ''}

# Opción A: ParamChem online (recomendado)
#   1. Ir a https://cgenff.umaryland.edu
#   2. Subir {target}.mol2 o {target}.sdf
#   3. Descargar {target}.str
#   4. Revisar penalizaciones (penalty score)
#      - Score < 10  → parámetros confiables
#      - Score 10-50 → revisar manualmente
#      - Score > 50  → requiere QM

# Opción B: CHARMM local (si disponible)
#   cgenff {target}.mol2 > {target}.str

# Post-parametrización: convertir a formato GROMACS
python cgenff_charmm2gmx.py {target} {target}.mol2 {target}.str charmm36.ff
"""

        readme = f"""# Parametrización CGenFF: {target}

## Engine
`cgenff` — CHARMM General Force Field

## Estado
{'⚠ **BLOCKING** — revisión manual requerida antes de producción' if is_blocking else 'Automático con revisión recomendada'}

## Outputs esperados
- `{target}.str`      — parámetros CGenFF
- `{target}.itp`      — topología GROMACS
- `{target}.prm`      — parámetros en formato CHARMM

## Criterio de aceptación
Penalty score < 10 en ParamChem para todos los átomos.
"""

        (step_dir / "commands.sh").write_text(commands)
        (step_dir / "README.md").write_text(readme)
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":          step.step_id,
            "stage":            step.stage.value,
            "engine":           step.engine,
            "target":           target,
            "blocking":         is_blocking,
            "step_type":        step.step_type.value,
            "expected_outputs": [
                f"{target}.str",
                f"{target}.itp",
            ],
        }, indent=4))

    # ── GAFF ──────────────────────────────────────────────────────────────────

    def _build_gaff(
        self,
        step:     SimulationStep,
        step_dir: Path,
        target:   str,
    ) -> None:

        commands = f"""#!/bin/bash
# ─── Parametrización GAFF: {target} ──────────────────────────────────────────

# Requiere: AmberTools (antechamber, parmchk2, tleap)

# 1. Asignar cargas AM1-BCC
antechamber \\
    -i {target}.sdf \\
    -fi sdf \\
    -o {target}.mol2 \\
    -fo mol2 \\
    -c bcc \\
    -s 2

# 2. Verificar parámetros faltantes
parmchk2 \\
    -i {target}.mol2 \\
    -f mol2 \\
    -o {target}.frcmod

# 3. Generar topología con tleap
#    (ver template.tleap en este directorio)
"""

        (step_dir / "commands.sh").write_text(commands)
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":          step.step_id,
            "stage":            step.stage.value,
            "engine":           step.engine,
            "target":           target,
            "blocking":         step.blocking,
            "step_type":        step.step_type.value,
            "expected_outputs": [
                f"{target}.mol2",
                f"{target}.frcmod",
            ],
        }, indent=4))

    # ── Genérico ──────────────────────────────────────────────────────────────

    def _build_generic(
        self,
        step:     SimulationStep,
        step_dir: Path,
        target:   str,
    ) -> None:

        (step_dir / "README.md").write_text(
            f"# Parametrización: {target}\n\n"
            f"Engine: {step.engine}\n\n"
            "Consultar documentación del forcefield para instrucciones.\n"
        )
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":   step.step_id,
            "stage":     step.stage.value,
            "engine":    step.engine,
            "target":    target,
            "blocking":  step.blocking,
            "step_type": step.step_type.value,
        }, indent=4))
