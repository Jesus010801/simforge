# builders/step_builders/validation_builder.py

from __future__ import annotations

from pathlib import Path
import json

from core.execution_models import SimulationStep


class ValidationBuilder:
    """
    Genera instrucciones para steps de validación manual.

    Estos steps no tienen script ejecutable — requieren
    intervención humana (revisar parametrización, validar pose, etc.)
    El builder escribe metadata.json con step_type correcto
    y README.md con instrucciones específicas.
    """

    def build(
        self,
        step:     SimulationStep,
        step_dir: Path,
    ) -> None:

        sid = step.step_id

        if "review_parametrization" in sid:
            self._build_param_review(step, step_dir)
        elif "validate_pose" in sid:
            self._build_pose_validation(step, step_dir)
        else:
            self._build_generic(step, step_dir)

    # ── Revisión de parametrización ───────────────────────────────────────────

    def _build_param_review(
        self,
        step:     SimulationStep,
        step_dir: Path,
    ) -> None:

        target = (
            step.target_components[0]
            if step.target_components
            else "ligand"
        )

        readme = f"""# Revisión Manual de Parametrización: {target}

## ⚠ Este step requiere intervención manual

## Qué revisar
{chr(10).join(f'- {n}' for n in step.notes)}

## Criterio de aceptación
- Penalty score < 10 en ParamChem para todos los átomos
- Geometría post-minimización similar al input (RMSD < 0.5Å)
- Cargas parciales razonables (sin valores > ±1.5e para átomos orgánicos típicos)

## Cuando esté listo
Marcar como completado y continuar con el siguiente step.

## Herramientas sugeridas
- ParamChem: https://cgenff.umaryland.edu
- VMD para inspección visual de cargas
- CHARMM o GROMACS para minimización de prueba
"""

        (step_dir / "README.md").write_text(readme)
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":   step.step_id,
            "stage":     step.stage.value,
            "engine":    step.engine,
            "target":    target,
            "blocking":  step.blocking,
            "step_type": step.step_type.value,
            "notes":     step.notes,
        }, indent=4))

    # ── Validación de pose ────────────────────────────────────────────────────

    def _build_pose_validation(
        self,
        step:     SimulationStep,
        step_dir: Path,
    ) -> None:

        target = (
            step.target_components[0]
            if step.target_components
            else "ligand"
        )

        readme = f"""# Validación de Pose Inicial: {target}

## Qué verificar
{chr(10).join(f'- {n}' for n in step.notes)}

## Cómo validar
1. Abrir el complejo en PyMOL o VMD
2. Verificar que el ligando está dentro del sitio activo
3. Revisar clashes (distancias < 2Å entre átomos no enlazados)
4. Confirmar contactos esperados con residuos clave

## Comandos PyMOL útiles
```
# Ver contactos dentro de 4Å
select contacts, byres ({target} around 4)
show sticks, contacts

# Verificar clashes
find_clashes {target}, cutoff=2.0
```
"""

        (step_dir / "README.md").write_text(readme)
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":   step.step_id,
            "stage":     step.stage.value,
            "engine":    step.engine,
            "target":    target,
            "blocking":  step.blocking,
            "step_type": step.step_type.value,
            "notes":     step.notes,
        }, indent=4))

    # ── Genérico ──────────────────────────────────────────────────────────────

    def _build_generic(
        self,
        step:     SimulationStep,
        step_dir: Path,
    ) -> None:

        (step_dir / "README.md").write_text(
            f"# Validación: {step.step_id}\n\n"
            f"Engine: {step.engine}\n\n"
            f"Notas:\n"
            + "\n".join(f"- {n}" for n in step.notes)
        )
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":   step.step_id,
            "stage":     step.stage.value,
            "engine":    step.engine,
            "blocking":  step.blocking,
            "step_type": step.step_type.value,
        }, indent=4))