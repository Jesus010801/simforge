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
        elif sid == "solvate_membrane":
            self._build_solvate_membrane(step, step_dir, step_dir_map)
        elif sid == "clean_water":
            self._build_clean_water(step, step_dir, step_dir_map)
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

        p             = step.params
        box_type      = p.get("box_type",      "triclinic")
        box_distance  = p.get("box_distance",   1.2)
        water_gro     = p.get("water_gro",     "spc216.gro")
        water_model   = p.get("water_model",   "tip3p")
        box_alignment = p.get("box_alignment", "")  # "principal_axes" → adds -princ

        # -princ aligns the molecule along principal axes before box calculation.
        # Required for elongated peptides to avoid oversized boxes, but BLOCKS
        # with an interactive "Select group" TTY prompt. Use only as explicit opt-in.
        princ_flag = " \\\n    -princ" if box_alignment == "principal_axes" else ""

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

# Copiar topología — gmx solvate la modifica in-place (añade SOL).
# La copia local asegura que assemble_system/topol.top quede intacta
# y que los pasos siguientes lean la topología actualizada desde aquí.
cp "$ASSEMBLE_DIR/topol.top" topol.top

# Definir caja de simulación
gmx editconf \\
    -f "$ASSEMBLE_DIR/complex.gro" \\
    -o box.gro \\
    -c \\
    -d {box_distance} \\
    -bt {box_type}{princ_flag}

# Agregar agua
gmx solvate \\
    -cp box.gro \\
    -cs {water_gro} \\
    -o solvated.gro \\
    -p topol.top
"""

        (step_dir / "run.sh").write_text(script)
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":          step.step_id,
            "stage":            step.stage.value,
            "engine":           step.engine,
            "step_type":        step.step_type.value,
            "blocking":         step.blocking,
            "generated_by":     "AssemblyBuilder",
            "expected_outputs": ["solvated.gro", "topol.top"],
            "params": {
                "box_type": box_type, "box_distance": box_distance,
                "water_model": water_model, "water_gro": water_gro,
                "box_alignment": box_alignment,
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

        # Input GRO: clean_water (system_clean.gro) > solvate step (solvated.gro)
        clean_dir   = step_dir_map.get("clean_water")
        solvate_dir = next(
            (step_dir_map[d] for d in step.depends_on if "solvate" in d and d in step_dir_map),
            None,
        )
        if clean_dir:
            input_ref = _rel(step_dir, clean_dir)
            input_gro = "system_clean.gro"
        elif solvate_dir:
            input_ref = _rel(step_dir, solvate_dir)
            input_gro = "solvated.gro"
        else:
            input_ref = "../solvate_system"
            input_gro = "solvated.gro"

        # topol.top chain: clean_water > solvate_* > assemble_system
        # Both solvate_membrane and solvate_system now keep a local topol.top.
        topol_src_dir = (
            clean_dir
            or solvate_dir
            or step_dir_map.get("assemble_system")
        )
        topol_src_ref = _rel(step_dir, topol_src_dir) if topol_src_dir else "../assemble_system"

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
INPUT_DIR="{input_ref}"
TOPOL_SRC="{topol_src_ref}"

# Copiar topología — gmx genion la modifica in-place (reemplaza SOL por iones).
# La copia local garantiza que el paso anterior quede sin modificar.
cp "$TOPOL_SRC/topol.top" topol.top

gmx grompp \\
    -f ions.mdp \\
    -c "$INPUT_DIR/{input_gro}" \\
    -p topol.top \\
    -o ions.tpr \\
    -maxwarn 2

echo "SOL" | gmx genion \\
    -s ions.tpr \\
    -o aaions.gro \\
    -p topol.top \\
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
            "expected_outputs": ["aaions.gro", "topol.top"],
            "params": {
                "concentration": concentration,
                "positive_ion":  positive_ion,
                "negative_ion":  negative_ion,
            },
        }, indent=4))

    # ── solvate_membrane ──────────────────────────────────────────────────────

    def _build_solvate_membrane(
        self,
        step:         SimulationStep,
        step_dir:     Path,
        step_dir_map: dict[str, Path],
    ) -> None:
        """
        Solvatación de sistema membrana.

        El input es converged.gro del shrink loop (membrane_embedding).
        gmx solvate añade agua SPC/E o TIP3P respetando la forma de la caja.
        La topología se copia localmente para que genion/clean_water la modifiquen
        sin alterar el directorio generate_topology.
        """
        p         = step.params
        water_gro = p.get("water_gro", "spc216.gro")

        # Input GRO viene del shrink loop
        embed_dir = next(
            (step_dir_map[d] for d in step.depends_on if "embedding" in d and d in step_dir_map),
            None,
        )
        embed_ref = _rel(step_dir, embed_dir) if embed_dir else "../membrane_embedding"

        # Topología base desde generate_topology
        topol_dir = step_dir_map.get("generate_topology")
        topol_ref = _rel(step_dir, topol_dir) if topol_dir else "../generate_topology"

        script = f"""#!/bin/bash
# ─── Solvatación (sistema membrana) ─────────────────────────────────────────
# gmx solvate añade agua respetando la caja ya definida por el shrink loop.
# La topología se copia aquí para que genion/clean_water la modifiquen
# sin alterar el directorio generate_topology.
EMBED_DIR="{embed_ref}"
TOPOL_DIR="{topol_ref}"

# Copiar topología — los pasos siguientes la modifican en local
cp "$TOPOL_DIR/topol.top" topol.top

gmx solvate \\
    -cp "$EMBED_DIR/converged.gro" \\
    -cs {water_gro} \\
    -o solvated.gro \\
    -p topol.top
"""

        (step_dir / "run.sh").write_text(script)
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":          step.step_id,
            "stage":            step.stage.value,
            "engine":           step.engine,
            "step_type":        step.step_type.value,
            "blocking":         step.blocking,
            "generated_by":     "AssemblyBuilder",
            "expected_outputs": ["solvated.gro", "topol.top"],
            "params":           {"water_gro": water_gro},
        }, indent=4))

    # ── clean_water ───────────────────────────────────────────────────────────

    def _build_clean_water(
        self,
        step:         SimulationStep,
        step_dir:     Path,
        step_dir_map: dict[str, Path],
    ) -> None:
        """
        Elimina moléculas de agua del interior de la bicapa lipídica.

        Usa WaterDeletorAdapter (reimplementación Python de water_deletor.pl).
        El script actualiza automáticamente el conteo SOL en topol.top.
        """
        p           = step.params
        ref_atom    = p.get("ref_atom",    "O33")
        middle_atom = p.get("middle_atom", "C50")
        nwater      = p.get("nwater",      3)

        solvate_dir = next(
            (step_dir_map[d] for d in step.depends_on if "solvate" in d and d in step_dir_map),
            None,
        )
        solvate_ref = _rel(step_dir, solvate_dir) if solvate_dir else "../solvate_membrane"

        script = f"""#!/usr/bin/env python3
# ─── Eliminar agua interior de bicapa ────────────────────────────────────────
# Reimplementación Python de water_deletor.pl (Lemkul 2017).
# Ejecutar: python3 run_clean_water.py
import sys
from pathlib import Path

# Resolución de paths relativa a la ubicación del script
SCRIPT_DIR  = Path(__file__).parent.resolve()
PROJECT_ROOT = SCRIPT_DIR.parent.parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from adapters.water_deletor_adapter import WaterDeletorAdapter

SOLVATE_DIR = (SCRIPT_DIR / "{solvate_ref}").resolve()
gro_in      = SOLVATE_DIR / "solvated.gro"
gro_out     = SCRIPT_DIR / "system_clean.gro"
topol_src   = SOLVATE_DIR / "topol.top"
topol_local = SCRIPT_DIR / "topol.top"

adapter = WaterDeletorAdapter()
result  = adapter.run(
    gro_in=gro_in,
    gro_out=gro_out,
    ref_atom="{ref_atom}",
    middle_atom="{middle_atom}",
    nwater={nwater},
    verbose=True,
)

if not result.success:
    print(f"ERROR: {{result.error_message}}", file=sys.stderr)
    sys.exit(1)

print(result.stdout)
waters_removed = result.metadata["waters_removed"]

# Actualizar conteo SOL en topol.top ─────────────────────────────────────────
import re, shutil
shutil.copy2(topol_src, topol_local)
text = topol_local.read_text()

# Encuentra la línea SOL en la sección [ molecules ]
def update_sol_count(text, n_removed):
    lines = text.splitlines()
    out = []
    for line in lines:
        m = re.match(r'^(SOL)\\s+(\\d+)', line)
        if m:
            old = int(m.group(2))
            new = old - n_removed
            print(f"  topol.top SOL: {{old}} → {{new}}")
            line = f"SOL              {{new}}"
        out.append(line)
    return "\\n".join(out)

topol_local.write_text(update_sol_count(text, waters_removed) + "\\n")
print(f"Output: {{gro_out}}")
print(f"topol.top updated: {{topol_local}}")
"""

        (step_dir / "run_clean_water.py").write_text(script)
        # Wrapper bash para compatibilidad con run.sh convention
        (step_dir / "run.sh").write_text(
            "#!/bin/bash\n"
            "# ─── Clean water (bilayer interior) ─────────────────────────────\n"
            'python3 "$(dirname "$0")/run_clean_water.py"\n'
        )
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":          step.step_id,
            "stage":            step.stage.value,
            "engine":           step.engine,
            "step_type":        step.step_type.value,
            "blocking":         step.blocking,
            "generated_by":     "AssemblyBuilder",
            "expected_outputs": ["system_clean.gro", "topol.top"],
            "params":           {"ref_atom": ref_atom, "middle_atom": middle_atom, "nwater": nwater},
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
