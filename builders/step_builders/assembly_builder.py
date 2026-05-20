# builders/step_builders/assembly_builder.py

from __future__ import annotations

from pathlib import Path
import json

from core.execution_models import SimulationStep


class AssemblyBuilder:
    """
    Genera instrucciones y scripts para el stage de assembly.

    Steps posibles:
        assemble_system  → combinar proteína + ligandos
        solvate_system   → gmx solvate
        add_ions         → gmx genion
        build_membrane   → CHARMM-GUI (externo)
    """

    def build(
        self,
        step:     SimulationStep,
        step_dir: Path,
    ) -> None:

        sid = step.step_id

        if sid == "assemble_system":
            self._build_assemble(step, step_dir)
        elif sid == "solvate_system":
            self._build_solvate(step, step_dir)
        elif sid == "add_ions":
            self._build_ions(step, step_dir)
        elif sid == "build_membrane":
            self._build_membrane(step, step_dir)
        else:
            self._build_generic(step, step_dir)

    # ── assemble_system ───────────────────────────────────────────────────────

    def _build_assemble(
        self,
        step:     SimulationStep,
        step_dir: Path,
    ) -> None:

        script = """#!/bin/bash
# ─── Assembly: combinar proteína y ligandos ───────────────────────────────────

# Combinar GRO de proteína + ligandos
# Ajustar según número de ligandos del sistema

# Proteína procesada (del step de preparation)
PROTEIN="../../01_prepare_protein_1/protein_1_processed.gro"

# Ligandos parametrizados
LIGAND_1="../../03_parametrize_substrate_1/substrate_1.gro"
LIGAND_2="../../04_parametrize_ligand_1/ligand_1.gro"

# Combinar estructuras
cat $PROTEIN $LIGAND_1 $LIGAND_2 > complex_raw.gro

# Actualizar número de átomos en la primera línea (suma total)
# Editar manualmente o usar script de python:
python3 -c "
import sys
lines = open('complex_raw.gro').readlines()
n_atoms = sum(1 for l in lines[2:-1] if l.strip())
lines[1] = f'{n_atoms}\\n'
open('complex.gro', 'w').writelines(lines)
print(f'Complex: {n_atoms} atoms')
"

# Combinar topologías (editar topol.top manualmente para incluir ligandos)
echo "Editar topol.top para incluir ligand ITP files"
echo "Agregar al final de topol.top:"
echo '; Ligandos'
echo '#include "substrate_1.itp"'
echo '#include "ligand_1.itp"'
"""

        (step_dir / "run.sh").write_text(script)
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":          step.step_id,
            "stage":            step.stage.value,
            "engine":           step.engine,
            "step_type":        step.step_type.value,
            "expected_outputs": ["complex.gro", "topol.top"],
        }, indent=4))

    # ── solvate_system ────────────────────────────────────────────────────────

    def _build_solvate(
        self,
        step:     SimulationStep,
        step_dir: Path,
    ) -> None:

        script = """#!/bin/bash
# ─── Solvatación ─────────────────────────────────────────────────────────────

# Definir caja de simulación
gmx editconf \\
    -f ../assemble_system/complex.gro \\
    -o box.gro \\
    -c \\
    -d 1.2 \\
    -bt dodecahedron

# Agregar agua TIP3P
gmx solvate \\
    -cp box.gro \\
    -cs spc216.gro \\
    -o solvated.gro \\
    -p ../assemble_system/topol.top
"""

        (step_dir / "run.sh").write_text(script)
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":          step.step_id,
            "stage":            step.stage.value,
            "engine":           step.engine,
            "step_type":        step.step_type.value,
            "expected_outputs": ["solvated.gro"],
        }, indent=4))

    # ── add_ions ──────────────────────────────────────────────────────────────

    def _build_ions(
        self,
        step:     SimulationStep,
        step_dir: Path,
    ) -> None:

        script = """#!/bin/bash
# ─── Adición de iones ────────────────────────────────────────────────────────

# Preparar archivo de entrada para genion
gmx grompp \\
    -f ions.mdp \\
    -c ../solvate_system/solvated.gro \\
    -p ../assemble_system/topol.top \\
    -o ions.tpr \\
    -maxwarn 2

# Agregar iones (neutralizar + 0.15M NaCl)
echo "SOL" | gmx genion \\
    -s ions.tpr \\
    -o aaions.gro \\
    -p ../assemble_system/topol.top \\
    -pname NA \\
    -nname CL \\
    -neutral \\
    -conc 0.15
"""

        ions_mdp = """; ions.mdp — mínimo para genion
integrator  = steep
nsteps      = 0
pbc         = xyz
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
"""

        (step_dir / "run.sh").write_text(script)
        (step_dir / "ions.mdp").write_text(ions_mdp)
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":          step.step_id,
            "stage":            step.stage.value,
            "engine":           step.engine,
            "step_type":        step.step_type.value,
            "expected_outputs": ["aaions.gro"],
        }, indent=4))

    # ── build_membrane ────────────────────────────────────────────────────────

    def _build_membrane(
        self,
        step:     SimulationStep,
        step_dir: Path,
    ) -> None:

        (step_dir / "README.md").write_text(
            "# Build Membrane\n\n"
            "Este step requiere CHARMM-GUI (externo).\n\n"
            "1. Ir a https://charmm-gui.org → Membrane Builder\n"
            "2. Subir la proteína procesada\n"
            "3. Configurar la bicapa lipídica\n"
            "4. Descargar el sistema y continuar desde aquí\n"
        )
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":   step.step_id,
            "stage":     step.stage.value,
            "engine":    step.engine,
            "step_type": "external",
        }, indent=4))

    # ── genérico ──────────────────────────────────────────────────────────────

    def _build_generic(
        self,
        step:     SimulationStep,
        step_dir: Path,
    ) -> None:

        (step_dir / "README.md").write_text(
            f"# Assembly: {step.step_id}\n\n"
            f"Engine: {step.engine}\n\n"
            "Instrucciones específicas pendientes.\n"
        )
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":   step.step_id,
            "stage":     step.stage.value,
            "engine":    step.engine,
            "step_type": step.step_type.value,
        }, indent=4))
