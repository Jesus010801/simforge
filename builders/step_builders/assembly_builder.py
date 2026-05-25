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

        if sid == "match_box_to_bilayer":
            from builders.step_builders.match_box_builder import MatchBoxBuilder
            MatchBoxBuilder().build(step, step_dir, step_dir_map)
        elif sid == "embed_in_bilayer":
            self._build_embed_in_bilayer(step, step_dir, step_dir_map)
        elif sid == "assemble_system":
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
# Outputs: system_clean.gro, topol.top, water_report.json
import sys, re, shutil, json
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.resolve()

def _find_root(start):
    p = start
    for _ in range(12):
        if (p / "adapters").is_dir() and (p / "core").is_dir():
            return p
        p = p.parent
    raise RuntimeError(f"SimForge project root not found from {{start}}")

PROJECT_ROOT = _find_root(SCRIPT_DIR)
sys.path.insert(0, str(PROJECT_ROOT))

from adapters.water_deletor_adapter import WaterDeletorAdapter
from validators.membrane_validators import validate_no_water_in_bilayer

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
shutil.copy2(topol_src, topol_local)
text = topol_local.read_text()

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

# ── Water gate: verify no waters remain inside bilayer core ───────────────────
# Threshold: >5 waters → error (gate blocks); 1-5 waters → warning; 0 → pass
_WATER_BLOCK_THRESHOLD = 5
wv = validate_no_water_in_bilayer(gro_out, headgroup_atom="{ref_atom}", tail_atom="{middle_atom}")
n_remain = wv.n_waters_in_bilayer
_w_errors   = [wv.message] if n_remain > _WATER_BLOCK_THRESHOLD else []
_w_warnings = [wv.message] if 0 < n_remain <= _WATER_BLOCK_THRESHOLD else []
w_report = {{
    "passed":             n_remain == 0,
    "waters_removed":     waters_removed,
    "n_waters_remaining": n_remain,
    "bilayer_z_min_nm":   wv.bilayer_z_min_nm,
    "bilayer_z_max_nm":   wv.bilayer_z_max_nm,
    "message":            wv.message,
    "errors":             _w_errors,
    "warnings":           _w_warnings,
    "confidence":         1.0,
}}
(SCRIPT_DIR / "water_report.json").write_text(json.dumps(w_report, indent=2))
print(f"[water_gate] {{wv.message}}")
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
            "automation_level": "automated",
            "blocking":         step.blocking,
            "generated_by":     "AssemblyBuilder",
            "gate":             {"type": "water_report"},
            "expected_outputs": ["system_clean.gro", "topol.top", "water_report.json"],
            "params":           {"ref_atom": ref_atom, "middle_atom": middle_atom, "nwater": nwater},
        }, indent=4))

    # ── embed_in_bilayer ──────────────────────────────────────────────────────

    def _build_embed_in_bilayer(
        self,
        step:         SimulationStep,
        step_dir:     Path,
        step_dir_map: dict[str, Path],
    ) -> None:
        """
        Embeds protein_boxed.gro into a pre-built bilayer using MoveMembAdapter
        (Python reimplementation of MoveMemb.f — no gfortran required).

        Inputs (resolved from DAG):
            match_box_to_bilayer/protein_boxed.gro
            <bilayer_file>  (looked up in CWD then Prot-Memb_FILES/)

        Outputs:
            system.gro       — protein + shifted bilayer (foundational artifact)
            strong_posre.itp — gmx genrestr FC=100000 on Protein group
        """
        p            = step.params
        bilayer_file  = p.get("bilayer_file", "dppc512_whole.gro")
        lipid         = p.get("lipid", "DPPC")
        # GRO residue name differs from the common lipid name (e.g. "DPP" vs "DPPC")
        lipid_resname = p.get("lipid_residue_name", lipid)

        match_box_dir = step_dir_map.get("match_box_to_bilayer")
        match_box_ref = _rel(step_dir, match_box_dir) if match_box_dir else "../match_box_to_bilayer"

        script = f"""#!/usr/bin/env python3
# ─── Embed protein in bilayer ─────────────────────────────────────────────────
# Uses MoveMembAdapter (Python reimpl of MoveMemb.f) to align bilayer midplane
# with protein Z-centre, then generates strong position restraints.
# Inputs:  <match_box_to_bilayer>/protein_boxed.gro  +  {bilayer_file}
# Outputs: system.gro, strong_posre.itp, overlap_report.json
import sys, subprocess, json
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.resolve()

def _find_root(start):
    p = start
    for _ in range(12):
        if (p / "adapters").is_dir() and (p / "core").is_dir():
            return p
        p = p.parent
    raise RuntimeError(f"SimForge project root not found from {{start}}")

PROJECT_ROOT = _find_root(SCRIPT_DIR)
sys.path.insert(0, str(PROJECT_ROOT))

from adapters.movememb_adapter import MoveMembAdapter
from validators.membrane_validators import validate_no_overlap

MATCH_BOX_DIR = (SCRIPT_DIR / "{match_box_ref}").resolve()
BILAYER_FILE  = "{bilayer_file}"
LIPID_RESNAME = "{lipid_resname}"   # GRO residue name (e.g. "DPP" for DPPC OPLS-AA)

# ── Resolve bilayer GRO ───────────────────────────────────────────────────────
bilayer_path = Path(BILAYER_FILE)
if not bilayer_path.exists():
    candidate = PROJECT_ROOT / "Prot-Memb_FILES" / BILAYER_FILE
    if candidate.exists():
        bilayer_path = candidate
    else:
        print(f"ERROR: bilayer '{{BILAYER_FILE}}' not found in CWD or Prot-Memb_FILES/", file=sys.stderr)
        sys.exit(1)

protein_gro = MATCH_BOX_DIR / "protein_boxed.gro"
gro_out     = SCRIPT_DIR / "system.gro"

if not protein_gro.exists():
    print(f"ERROR: protein_boxed.gro not found at {{protein_gro}}", file=sys.stderr)
    sys.exit(1)

# ── MoveMemb: align bilayer Z-midplane with protein Z-centre ──────────────────
adapter = MoveMembAdapter()
result  = adapter.run(
    protein_gro=protein_gro,
    bilayer_gro=bilayer_path,
    gro_out=gro_out,
)

if not result.success:
    print(f"ERROR: MoveMembAdapter failed: {{result.error_message}}", file=sys.stderr)
    sys.exit(1)

print(result.stdout)
m = result.metadata
print(f"Z-shift:          {{m['z_shift_nm']:+.4f}} nm")
print(f"Protein Z:        [{{m['protein_z_min']:.3f}}, {{m['protein_z_max']:.3f}}] nm")
print(f"Bilayer Z (orig): [{{m['bilayer_z_min']:.3f}}, {{m['bilayer_z_max']:.3f}}] nm")
print(f"Combined atoms:   {{m['atoms_total']}}")

# ── gmx genrestr — strong position restraints on protein heavy atoms ──────────
posre_out = SCRIPT_DIR / "strong_posre.itp"
ret = subprocess.run(
    ["gmx", "genrestr",
     "-f", str(gro_out),
     "-o", str(posre_out),
     "-fc", "100000", "100000", "100000"],
    input="Protein\\n", text=True, capture_output=True,
)
if ret.returncode != 0:
    print(f"ERROR: gmx genrestr failed:\\n{{ret.stderr}}", file=sys.stderr)
    sys.exit(1)
print(f"Position restraints: {{posre_out}}")

# ── gmx editconf — renumber residues from 1 ──────────────────────────────────
ret = subprocess.run(
    ["gmx", "editconf",
     "-f", str(gro_out),
     "-o", str(gro_out),
     "-resnr", "1"],
    capture_output=True, text=True,
)
if ret.returncode != 0:
    print(f"ERROR: gmx editconf -resnr 1 failed:\\n{{ret.stderr}}", file=sys.stderr)
    sys.exit(1)
print(f"Output: {{gro_out}}")

# ── Overlap gate: check for protein–lipid clashes ────────────────────────────
ov = validate_no_overlap(gro_out, lipid_residue_name=LIPID_RESNAME)
ov_report = {{
    "passed":          ov.n_clashes == 0,
    "n_clashes":       ov.n_clashes,
    "n_protein_atoms": ov.n_protein_atoms,
    "n_lipid_atoms":   ov.n_lipid_atoms,
    "message":         ov.message,
    "errors":          [] if ov.n_clashes == 0 else [ov.message],
    "warnings":        [],
    "confidence":      1.0,
}}
(SCRIPT_DIR / "overlap_report.json").write_text(json.dumps(ov_report, indent=2))
print(f"[overlap_gate] {{ov.message}}")
"""

        (step_dir / "run_embed.py").write_text(script)
        (step_dir / "run.sh").write_text(
            "#!/bin/bash\n"
            "# ─── embed_in_bilayer (automatic) ──────────────────────────────────\n"
            'python3 "$(dirname "$0")/run_embed.py"\n'
        )
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":          step.step_id,
            "stage":            step.stage.value,
            "engine":           step.engine,
            "step_type":        "automatic",
            "automation_level": "automated",
            "generated_by":     "AssemblyBuilder",
            "gate":             {"type": "overlap_report"},
            "expected_outputs": ["system.gro", "strong_posre.itp", "overlap_report.json"],
            "params": {
                "bilayer_file":       bilayer_file,
                "lipid":              lipid,
                "lipid_residue_name": lipid_resname,
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
