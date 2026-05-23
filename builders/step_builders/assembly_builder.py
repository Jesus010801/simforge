# builders/step_builders/assembly_builder.py

from __future__ import annotations

from pathlib import Path
import json

from core.execution_models import SimulationStep
from builders.step_builders._utils import rel as _rel


def _component_gro(step_id: str) -> str | None:
    """
    Expected GRO output filename for a dep step used in assembly.
    Returns None for steps that don't produce a structure for assembly
    (e.g. prepare_substrate_*, prepare_ligand_*, review_*, validate_*).
    """
    if step_id.startswith("prepare_protein_"):
        component = step_id[len("prepare_"):]          # "protein_1"
        return f"{component}_processed.gro"
    if step_id.startswith("parametrize_"):
        component = step_id[len("parametrize_"):]      # "substrate_1", "ligand_1"
        return f"{component}.gro"
    return None


class AssemblyBuilder:
    """
    Genera instrucciones y scripts para el stage de assembly.

    Steps posibles:
        assemble_system  → combinar proteína + ligandos
        solvate_system   → gmx solvate
        add_ions         → gmx genion
        build_membrane   → CHARMM-GUI (externo)

    Reads all scientific params from step.params (populated by decision_engine).
    Reads inter-step paths from step_dir_map (built by WorkspaceBuilder).
    """

    def build(
        self,
        step:         SimulationStep,
        step_dir:     Path,
        step_dir_map: dict[str, Path] = {},
    ) -> None:

        sid = step.step_id

        if sid == "assemble_system":
            self._build_assemble(step, step_dir, step_dir_map)
        elif sid == "solvate_system":
            self._build_solvate(step, step_dir, step_dir_map)
        elif sid == "add_ions":
            self._build_ions(step, step_dir, step_dir_map)
        elif sid == "build_membrane":
            self._build_membrane(step, step_dir, step_dir_map)
        else:
            self._build_generic(step, step_dir, step_dir_map)

    # ── assemble_system ───────────────────────────────────────────────────────

    def _build_assemble(
        self,
        step:         SimulationStep,
        step_dir:     Path,
        step_dir_map: dict[str, Path],
    ) -> None:

        # Derive input GRO paths from DAG depends_on
        inputs: list[tuple[str, str]] = []   # (VAR_NAME, rel/path/to/file.gro)
        for dep_id in step.depends_on:
            dep_dir = step_dir_map.get(dep_id)
            if dep_dir is None:
                continue
            gro = _component_gro(dep_id)
            if gro is None:
                continue
            component = dep_id.split("_", 1)[1]          # "protein_1", "ligand_1"
            var_name  = component.upper()                 # "PROTEIN_1", "LIGAND_1"
            inputs.append((var_name, f"{_rel(step_dir, dep_dir)}/{gro}"))

        # Locate protein prep dir to hand off topol.top
        protein_prep_dep = next(
            (d for d in step.depends_on if d.startswith("prepare_protein_") and d in step_dir_map),
            None,
        )
        protein_prep_ref = (
            _rel(step_dir, step_dir_map[protein_prep_dep])
            if protein_prep_dep else None
        )

        var_lines = "\n".join(f'{name}="{path}"' for name, path in inputs)
        cat_args  = " ".join(f"${name}" for name, _ in inputs)
        itp_lines = "\n".join(
            f'echo \'#include "{dep_id.split("_", 1)[1]}.itp"\''
            for dep_id in step.depends_on
            if dep_id.startswith("parametrize_")
        )

        topol_block = ""
        if protein_prep_ref:
            topol_block = f"""
# Handoff topology from protein prep
PROTEIN_PREP_DIR="{protein_prep_ref}"
cp "$PROTEIN_PREP_DIR/topol.top" topol.top
cp "$PROTEIN_PREP_DIR/posre.itp" posre.itp
"""

        script = f"""#!/bin/bash
# ─── Assembly: combinar proteína y ligandos ───────────────────────────────────
# Paths resueltos desde DAG — no editar manualmente

{var_lines}
{topol_block}
# Combinar estructuras
cat {cat_args} > complex_raw.gro

# Actualizar número de átomos
python3 -c "
lines = open('complex_raw.gro').readlines()
n_atoms = sum(1 for l in lines[2:-1] if l.strip())
lines[1] = f'{{n_atoms}}\\n'
open('complex.gro', 'w').writelines(lines)
print(f'Complex: {{n_atoms}} atoms')
"

# Combinar topologías
echo "Editar topol.top para incluir ligand ITP files:"
{itp_lines}
"""

        params_effective = {
            "inputs": [{"var": n, "path": p} for n, p in inputs],
        }

        (step_dir / "run.sh").write_text(script)
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":          step.step_id,
            "stage":            step.stage.value,
            "engine":           step.engine,
            "step_type":        step.step_type.value,
            "blocking":         step.blocking,
            "generated_by":     "AssemblyBuilder",
            "expected_outputs": ["complex.gro", "topol.top"],
            "params":           params_effective,
        }, indent=4))

    # ── solvate_system ────────────────────────────────────────────────────────

    def _build_solvate(
        self,
        step:         SimulationStep,
        step_dir:     Path,
        step_dir_map: dict[str, Path],
    ) -> None:

        p            = step.params
        box_type     = p.get("box_type",     "dodecahedron")
        box_distance = p.get("box_distance",  1.2)
        water_gro    = p.get("water_gro",    "spc216.gro")
        water_model  = p.get("water_model",  "tip3p")

        # Resolve assemble_system path (direct dep)
        assemble_dir = next(
            (step_dir_map[d] for d in step.depends_on if "assemble" in d and d in step_dir_map),
            None,
        )
        assemble_ref = _rel(step_dir, assemble_dir) if assemble_dir else "../assemble_system"

        script = f"""#!/bin/bash
# ─── Solvatación ─────────────────────────────────────────────────────────────
# water_model={water_model}  box_type={box_type}  d={box_distance}nm
# Paths resueltos desde DAG

ASSEMBLE_DIR="{assemble_ref}"

# Definir caja de simulación
gmx editconf \\
    -f "$ASSEMBLE_DIR/complex.gro" \\
    -o box.gro \\
    -c \\
    -d {box_distance} \\
    -bt {box_type} \\
    -princ

# Agregar agua
gmx solvate \\
    -cp box.gro \\
    -cs {water_gro} \\
    -o solvated.gro \\
    -p "$ASSEMBLE_DIR/topol.top"
"""

        (step_dir / "run.sh").write_text(script)
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":          step.step_id,
            "stage":            step.stage.value,
            "engine":           step.engine,
            "step_type":        step.step_type.value,
            "blocking":         step.blocking,
            "generated_by":     "AssemblyBuilder",
            "expected_outputs": ["solvated.gro"],
            "params": {
                "box_type": box_type, "box_distance": box_distance,
                "water_model": water_model, "water_gro": water_gro,
            },
        }, indent=4))

    # ── add_ions ──────────────────────────────────────────────────────────────

    def _build_ions(
        self,
        step:         SimulationStep,
        step_dir:     Path,
        step_dir_map: dict[str, Path],
    ) -> None:

        p             = step.params
        concentration = p.get("concentration",  0.15)
        positive_ion  = p.get("positive_ion",   "NA")
        negative_ion  = p.get("negative_ion",   "CL")

        # solvate_system is direct dep; assemble_system holds topol.top
        solvate_dir  = next(
            (step_dir_map[d] for d in step.depends_on if "solvate" in d and d in step_dir_map),
            None,
        )
        assemble_dir = step_dir_map.get("assemble_system")

        solvate_ref  = _rel(step_dir, solvate_dir)  if solvate_dir  else "../solvate_system"
        assemble_ref = _rel(step_dir, assemble_dir) if assemble_dir else "../assemble_system"

        ions_mdp = """; ions.mdp — mínimo para genion
integrator    = steep
nsteps        = 0
pbc           = xyz
cutoff-scheme = Verlet
coulombtype   = PME
rcoulomb      = 1.0
rvdw          = 1.0
"""

        script = f"""#!/bin/bash
# ─── Adición de iones ────────────────────────────────────────────────────────
# concentration={concentration}M  +={positive_ion}  -={negative_ion}
# Paths resueltos desde DAG

SOLVATE_DIR="{solvate_ref}"
ASSEMBLE_DIR="{assemble_ref}"

gmx grompp \\
    -f ions.mdp \\
    -c "$SOLVATE_DIR/solvated.gro" \\
    -p "$ASSEMBLE_DIR/topol.top" \\
    -o ions.tpr \\
    -maxwarn 2

echo "SOL" | gmx genion \\
    -s ions.tpr \\
    -o aaions.gro \\
    -p "$ASSEMBLE_DIR/topol.top" \\
    -pname {positive_ion} \\
    -nname {negative_ion} \\
    -neutral \\
    -conc {concentration}
"""

        (step_dir / "ions.mdp").write_text(ions_mdp)
        (step_dir / "run.sh").write_text(script)
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":          step.step_id,
            "stage":            step.stage.value,
            "engine":           step.engine,
            "step_type":        step.step_type.value,
            "blocking":         step.blocking,
            "generated_by":     "AssemblyBuilder",
            "expected_outputs": ["aaions.gro"],
            "params": {
                "concentration": concentration,
                "positive_ion":  positive_ion,
                "negative_ion":  negative_ion,
            },
        }, indent=4))

    # ── build_membrane ────────────────────────────────────────────────────────

    def _build_membrane(
        self,
        step:         SimulationStep,
        step_dir:     Path,
        step_dir_map: dict[str, Path],
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
        step:         SimulationStep,
        step_dir:     Path,
        step_dir_map: dict[str, Path],
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
