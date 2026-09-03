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
        elif sid == "assemble_system_topology":
            self._build_assemble_system_topology(step, step_dir, step_dir_map)
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

        # Input GRO from membrane_embedding (resolved directly from step_dir_map)
        embed_dir = step_dir_map.get("membrane_embedding")
        embed_ref = _rel(step_dir, embed_dir) if embed_dir else "../membrane_embedding"

        # Topología base: Phase 11 uses assemble_system_topology (fallback: generate_topology)
        topol_dir = (
            step_dir_map.get("assemble_system_topology")
            or step_dir_map.get("generate_topology")
        )
        topol_ref = _rel(step_dir, topol_dir) if topol_dir else "../assemble_system_topology"

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
        tm_residues_str = p.get("tm_residues")
        if tm_residues_str:
            from core.structural_annotation import residues_in_range
            tm_residues_literal = repr(sorted(residues_in_range(tm_residues_str)))
        else:
            tm_residues_literal = "[]"

        solvate_dir = next(
            (step_dir_map[d] for d in step.depends_on if "solvate" in d and d in step_dir_map),
            None,
        )
        solvate_ref = _rel(step_dir, solvate_dir) if solvate_dir else "../solvate_membrane"

        script = f"""#!/usr/bin/env python3
# ─── Eliminar agua interior de bicapa ────────────────────────────────────────
# Reimplementación Python de water_deletor.pl (Lemkul 2017).
# Outputs: system_clean.gro, topol.top, clean_water_report.json, water_report.json
import sys, re, shutil, json
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.resolve()

from adapters.water_deletor_adapter import WaterDeletorAdapter, _parse_gro
from validators.membrane_validators import validate_no_water_in_bilayer

SOLVATE_DIR = (SCRIPT_DIR / "{solvate_ref}").resolve()
gro_in      = SOLVATE_DIR / "solvated.gro"
gro_out     = SCRIPT_DIR / "system_clean.gro"
topol_src   = SOLVATE_DIR / "topol.top"
topol_local = SCRIPT_DIR / "topol.top"

TM_RESIDUES = set({tm_residues_literal})
if TM_RESIDUES:
    from validators.pore_hydration import clean_water_channel_aware
    report = clean_water_channel_aware(
        gro_in, gro_out, tm_residues=TM_RESIDUES, output_dir=SCRIPT_DIR,
        topol_in=topol_src, topol_out=topol_local,
    )
    print(json.dumps(report, indent=2))
    sys.exit(0)

# Count input SOL water molecules (one OW per molecule)
def _count_sol_ow(path, resname="SOL", ow="OW"):
    _, atoms, _ = _parse_gro(path)
    return sum(1 for l in atoms if l[5:10].strip() == resname and l[10:15].strip() == ow)

input_water_molecules = _count_sol_ow(gro_in)

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
n_molecules_removed = result.metadata["waters_removed"]
n_atoms_removed     = n_molecules_removed * {nwater}  # atoms per water molecule

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

topology_updated = False
try:
    topol_local.write_text(update_sol_count(text, n_molecules_removed) + "\\n")
    topology_updated = True
    print(f"topol.top updated: {{topol_local}}")
except Exception as _te:
    print(f"WARNING: topology update failed: {{_te}}", file=sys.stderr)

# ── Water gate: verify no OW remain inside the hydrophobic core ───────────────
# validate_no_water_in_bilayer uses mean-per-leaflet headgroup Z (same as
# WaterDeletorAdapter) — NOT min/max of all headgroups.  After deletion,
# n_remain should be 0 for a correct run.
# Threshold: >5 OW in core → error (gate blocks); 1-5 → warning; 0 → pass.
_WATER_BLOCK_THRESHOLD = 5
wv = validate_no_water_in_bilayer(gro_out, headgroup_atom="{ref_atom}", tail_atom="{middle_atom}")
n_remain = wv.n_waters_in_bilayer
_w_errors   = [wv.message] if n_remain > _WATER_BLOCK_THRESHOLD else []
_w_warnings = [wv.message] if 0 < n_remain <= _WATER_BLOCK_THRESHOLD else []

# ── Write clean_water_report.json (primary structured report) ─────────────────
clean_report = {{
    "input_water_molecules":          input_water_molecules,
    "n_water_molecules_removed":      n_molecules_removed,
    "n_water_atoms_removed":          n_atoms_removed,
    "n_water_oxygens_remaining_in_core": n_remain,
    "bilayer_midplane_z":             wv.bilayer_midplane_z,
    "core_z_min":                     wv.core_z_min,
    "core_z_max":                     wv.core_z_max,
    "z_bot_nm":                       result.metadata["z_bot_nm"],
    "z_top_nm":                       result.metadata["z_top_nm"],
    "output_gro_path":                str(gro_out),
    "topology_updated":               topology_updated,
    "cleanup_passed":                 n_remain == 0,
}}
(SCRIPT_DIR / "clean_water_report.json").write_text(json.dumps(clean_report, indent=2))
print(f"Output: {{gro_out}}")

# ── Write water_report.json (gate-compatible legacy format) ───────────────────
w_report = {{
    "passed":                            n_remain == 0,
    "n_water_molecules_removed":         n_molecules_removed,
    "n_water_atoms_removed":             n_atoms_removed,
    "n_waters_remaining":                n_remain,
    "n_water_oxygens_remaining_in_core": n_remain,
    "bilayer_midplane_z":                wv.bilayer_midplane_z,
    "core_z_min":                        wv.core_z_min,
    "core_z_max":                        wv.core_z_max,
    "bilayer_z_min_nm":                  wv.bilayer_z_min_nm,
    "bilayer_z_max_nm":                  wv.bilayer_z_max_nm,
    "cleanup_passed":                    n_remain == 0,
    "message":                           wv.message,
    "errors":                            _w_errors,
    "warnings":                          _w_warnings,
    "confidence":                        1.0,
}}
(SCRIPT_DIR / "water_report.json").write_text(json.dumps(w_report, indent=2))
print(f"[water_gate] {{wv.message}}")
if n_remain == 0:
    print("[water_gate] PASS — no water oxygens in bilayer hydrophobic core")
elif n_remain <= _WATER_BLOCK_THRESHOLD:
    print(f"[water_gate] WARNING — {{n_remain}} OW in core (boundary; may be acceptable)")
else:
    print(f"[water_gate] FAIL — {{n_remain}} OW remain in core; inspect system_clean.gro", file=sys.stderr)
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
            "expected_outputs": [
                "system_clean.gro",
                "topol.top",
                "clean_water_report.json",
                "water_report.json",
                "pore_hydration_report.json",
            ],
            "params": {"ref_atom": ref_atom, "middle_atom": middle_atom, "nwater": nwater},
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
            <bilayer_file>  (loaded from workspace membrane_assets/)

        Outputs:
            system.gro       — protein + shifted bilayer (foundational artifact)
            strong_posre.itp — gmx genrestr FC=100000 on Protein group
        """
        p            = step.params
        bilayer_file  = p.get("bilayer_file", "dppc512_whole.gro")
        lipid         = p.get("lipid", "DPPC")
        # GRO residue name differs from the common lipid name (e.g. "DPP" vs "DPPC")
        lipid_resname = p.get("lipid_residue_name", lipid)
        pre_exclude_trapped = bool(p.get("pre_exclude_trapped_lipids", True))
        max_pre_excluded    = int(p.get("max_pre_excluded_lipids", 50))
        tm_mask_validation  = bool(p.get("tm_mask_validation", True))
        tm_mask_policy      = str(p.get("tm_mask_policy", "warn"))
        tm_aware_exclusion  = bool(p.get("tm_aware_exclusion", False))
        tm_aware_exclusion_policy = str(p.get("tm_aware_exclusion_policy", "off"))
        tm_aware_max_removed = int(p.get("tm_aware_max_removed_lipids", 50))

        # ── Phase 10A: interface builder (lipid refill) params ─────────────────
        iface_b_enabled     = bool( p.get("iface_builder_enabled",     False))
        iface_b_policy      = str(  p.get("iface_builder_policy",      "warn"))
        iface_b_max_ins     = int(  p.get("iface_builder_max_inserted", 40))
        iface_b_max_per_cl  = int(  p.get("iface_builder_max_per_cluster", 4))
        iface_b_target_dist = float(p.get("iface_builder_target_dist", 0.45))
        iface_b_prot_clash  = float(p.get("iface_builder_prot_clash",  0.20))
        iface_b_lip_clash   = float(p.get("iface_builder_lip_clash",   0.18))
        iface_b_seed        = int(  p.get("iface_builder_seed",        17))

        match_box_dir = step_dir_map.get("match_box_to_bilayer")
        match_box_ref = _rel(step_dir, match_box_dir) if match_box_dir else "../match_box_to_bilayer"

        membrane_assets_dir = step_dir_map.get("__membrane_assets__")
        assets_ref = (
            _rel(step_dir, membrane_assets_dir)
            if membrane_assets_dir
            else "../../membrane_assets"
        )

        prot_top_dir = step_dir_map.get("generate_protein_topology")
        prot_top_ref = _rel(step_dir, prot_top_dir) if prot_top_dir else "../generate_protein_topology"

        # ── TM residue set (baked into the generated script at compile time) ──
        tm_residues_str = p.get("tm_residues")  # e.g. "51-75" or None
        if tm_residues_str:
            from core.structural_annotation import residues_in_range
            tm_set = residues_in_range(tm_residues_str)
            tm_value = f"set({sorted(tm_set)!r})"
        else:
            tm_value = "None"

        script = f"""#!/usr/bin/env python3
# ─── Embed protein in bilayer ─────────────────────────────────────────────────
# Uses MoveMembAdapter (Python reimpl of MoveMemb.f) to align bilayer midplane
# with TM-region Z-centre (or full protein Z-centre as fallback).
# Phase 8A: reads membrane_registration_report.json from match_box step and uses
#           the optimised shift when adaptive registration was enabled.
# Inputs:  <match_box_to_bilayer>/protein_boxed.gro  +  membrane_assets/{bilayer_file}
# Outputs: system.gro, strong_posre.itp, overlap_report.json,
#          pre_shrink_lipid_exclusion_report.json
import sys, subprocess, json
from pathlib import Path

SCRIPT_DIR   = Path(__file__).parent.resolve()
PROT_TOP_DIR = (SCRIPT_DIR / "{prot_top_ref}").resolve()

from adapters.movememb_adapter import MoveMembAdapter
from validators.membrane_validators import validate_no_overlap, detect_lipid_in_cavity
from validators.pre_shrink_exclusion import apply_pre_shrink_exclusion

MATCH_BOX_DIR = (SCRIPT_DIR / "{match_box_ref}").resolve()
BILAYER_FILE  = "{bilayer_file}"
LIPID_RESNAME = "{lipid_resname}"   # GRO residue name (e.g. "DPP" for DPPC OPLS-AA)

# TM residues baked in at compile time from structural_annotation.
# None  → adapter falls back to full protein Z-centre (emits warning in report).
TM_RESIDUES = {tm_value}

# Pre-shrink trapped-lipid exclusion config (baked in from membrane.embedding config).
PRE_EXCLUDE_TRAPPED = {pre_exclude_trapped}
MAX_PRE_EXCLUDED    = {max_pre_excluded}

# TM-aware mask validation config (baked in from membrane.embedding config).
TM_MASK_VALIDATION = {tm_mask_validation}
TM_MASK_POLICY     = "{tm_mask_policy}"

# TM-aware lipid exclusion config (baked in from membrane.embedding config).
TM_AWARE_EXCLUSION = {tm_aware_exclusion}
TM_AWARE_EXCLUSION_POLICY = "{tm_aware_exclusion_policy}"
TM_AWARE_MAX_REMOVED = {tm_aware_max_removed}

# ── Resolve bilayer GRO from workspace membrane_assets ────────────────────────
ASSETS_DIR   = (SCRIPT_DIR / "{assets_ref}").resolve()
bilayer_path = ASSETS_DIR / BILAYER_FILE
if not bilayer_path.exists():
    print(f"ERROR: bilayer '{{BILAYER_FILE}}' not found in membrane_assets: {{bilayer_path}}", file=sys.stderr)
    sys.exit(1)

protein_gro = MATCH_BOX_DIR / "protein_boxed.gro"
gro_out     = SCRIPT_DIR / "system.gro"

if not protein_gro.exists():
    print(f"ERROR: protein_boxed.gro not found at {{protein_gro}}", file=sys.stderr)
    sys.exit(1)

# ── Phase 8A/8B: read adaptive registration shift ────────────────────────────
# Phase 8B (surface-interference) takes precedence when available and enabled.
_explicit_z_shift = None
_reg_source = None

_surf_reg_path = MATCH_BOX_DIR / "membrane_surface_registration_report.json"
if _surf_reg_path.exists():
    try:
        _sreg = json.loads(_surf_reg_path.read_text())
        if _sreg.get("enabled", False) and _sreg.get("actual_z_shift_nm") is not None:
            _explicit_z_shift = float(_sreg["actual_z_shift_nm"])
            _sel_delta = _sreg.get("selected_shift_z_nm", 0.0)
            _reg_source = "Phase 8B (surface-interference)"
            print(f"[embed] Using Phase 8B surface-interference shift: {{_explicit_z_shift:+.3f}} nm "
                  f"(delta={{_sel_delta:+.3f}} nm)")
    except Exception as _e:
        print(f"WARNING: could not read membrane_surface_registration_report.json: {{_e}}", file=sys.stderr)

if _explicit_z_shift is None:
    _reg_report_path = MATCH_BOX_DIR / "membrane_registration_report.json"
    if _reg_report_path.exists():
        try:
            _reg = json.loads(_reg_report_path.read_text())
            if _reg.get("enabled", False) and _reg.get("actual_z_shift_nm") is not None:
                _explicit_z_shift = float(_reg["actual_z_shift_nm"])
                _sel_delta = _reg.get("selected_shift_z_nm", 0.0)
                _reg_source = "Phase 8A (TM-burial)"
                print(f"[embed] Using Phase 8A adaptive registration shift: {{_explicit_z_shift:+.3f}} nm "
                      f"(delta={{_sel_delta:+.3f}} nm from TM-center alignment)")
        except Exception as _e:
            print(f"WARNING: could not read membrane_registration_report.json: {{_e}}", file=sys.stderr)

# ── MoveMemb: align bilayer midplane with TM-region Z-centre ──────────────────
# When adaptive registration produced a shift, pass it explicitly so the adapter
# uses the optimised Z position rather than re-computing from TM-centre.
adapter = MoveMembAdapter()
result  = adapter.run(
    protein_gro=protein_gro,
    bilayer_gro=bilayer_path,
    gro_out=gro_out,
    z_shift_nm=_explicit_z_shift,   # None → auto-compute from TM centre
    tm_residues=TM_RESIDUES if _explicit_z_shift is None else None,
)

if not result.success:
    print(f"ERROR: MoveMembAdapter failed: {{result.error_message}}", file=sys.stderr)
    sys.exit(1)

print(result.stdout)
m = result.metadata

# ── Per-chain strong position restraints (local-index-correct) ───────────────
# Generate one strong_posre_<MolName>.itp per chain using LOCAL atom indices
# parsed from the chain ITP's [ atoms ] section.  Global atom indices from
# gmx genrestr on the full system GRO are wrong for multichain proteins because
# they exceed the local atom count of individual moleculetypes.
import re as _re_posre

def _parse_heavy_local_indices(top_path):
    # Returns 1-based local heavy-atom indices from the first [ atoms ] section.
    _in_atoms = False
    _idx_list = []
    for _l in top_path.read_text().splitlines():
        _ls = _l.strip()
        if _re_posre.match(r'^\[\s*atoms\s*\]', _ls, _re_posre.IGNORECASE):
            _in_atoms = True; continue
        if _in_atoms:
            if _ls.startswith('['):
                break
            if not _ls or _ls.startswith(';'):
                continue
            _pp = _ls.split()
            if len(_pp) >= 5:
                try:
                    _nr    = int(_pp[0])
                    _aname = _pp[4]
                    if not _aname.upper().startswith('H'):
                        _idx_list.append(_nr)
                except (ValueError, IndexError):
                    pass
    return _idx_list

def _write_posre_itp(out_path, indices, mol_name):
    _lines = [
        f"; Strong position restraints for {{mol_name}} — FC=100000 kJ/mol/nm²",
        "[ position_restraints ]",
        "; atom  funct   fcx     fcy     fcz",
    ] + [f"  {{idx:6d}}    1  100000  100000  100000" for idx in indices]
    out_path.write_text("\\n".join(_lines) + "\\n")

_posre_manifest = {{}}
_manifest_f_pm  = PROT_TOP_DIR / "protein_topology_manifest.json"
if _manifest_f_pm.exists():
    _pm = json.loads(_manifest_f_pm.read_text())
    _mol_names_pm = _pm.get("molecule_names", [])
    _mol_itps_pm  = _pm.get("molecule_itps",  [])
    if _mol_itps_pm:
        for _mn_pm, _itp_rel_pm in zip(_mol_names_pm, _mol_itps_pm):
            _itp_abs_pm = (Path(_itp_rel_pm) if Path(_itp_rel_pm).is_absolute()
                           else (PROT_TOP_DIR / _itp_rel_pm)).resolve()
            if not _itp_abs_pm.exists():
                print(f"WARNING: chain ITP not found: {{_itp_abs_pm}}", file=sys.stderr)
                continue
            _hidx_pm = _parse_heavy_local_indices(_itp_abs_pm)
            if not _hidx_pm:
                print(f"WARNING: no heavy atoms in {{_itp_abs_pm.name}}", file=sys.stderr)
                continue
            _sp_name_pm = f"strong_posre_{{_mn_pm}}.itp"
            _write_posre_itp(SCRIPT_DIR / _sp_name_pm, _hidx_pm, _mn_pm)
            _posre_manifest[_mn_pm] = _sp_name_pm
            print(f"[embed] {{_mn_pm}}: {{len(_hidx_pm)}} heavy atoms → {{_sp_name_pm}}")
    else:
        _pt_f_pm = PROT_TOP_DIR / "topol.top"
        _mn0_pm  = _mol_names_pm[0] if _mol_names_pm else "Protein"
        _hidx_pm = _parse_heavy_local_indices(_pt_f_pm) if _pt_f_pm.exists() else []
        if _hidx_pm:
            _write_posre_itp(SCRIPT_DIR / "strong_posre.itp", _hidx_pm, _mn0_pm)
            _posre_manifest[_mn0_pm] = "strong_posre.itp"
            print(f"[embed] {{_mn0_pm}}: {{len(_hidx_pm)}} heavy atoms → strong_posre.itp")
        else:
            print("WARNING: no heavy atoms found for strong restraints", file=sys.stderr)
else:
    print(f"WARNING: protein_topology_manifest.json not found at {{_manifest_f_pm}}", file=sys.stderr)
(SCRIPT_DIR / "strong_posre_manifest.json").write_text(json.dumps(_posre_manifest, indent=2))
print(f"Strong posre files: {{list(_posre_manifest.values())}}")

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

# ── Pre-shrink cavity-trapped lipid exclusion ─────────────────────────────────
# Removes lipids geometrically trapped inside the TM cavity before topology
# generation.  Only lipids angularly surrounded by protein CA atoms are removed
# (is_surrounded_by_protein=True) — surface-adjacent overlaps are left to inflategro.
# Topology is NOT modified here; generate_topology reads the filtered system.gro.
_pre_excl_report = {{
    "enabled":                               PRE_EXCLUDE_TRAPPED,
    "n_lipids_checked":                      0,
    "n_lipids_removed":                      0,
    "removed_resids":                        [],
    "removed_resnames":                      [],
    "removed_reasons":                       [],
    "n_hard_overlap_removed":                0,
    "n_low_overlap_geometrically_trapped_removed": 0,
    "n_candidate_lipids":                    0,
    "n_removed_by_membrane_mask":            0,
    "n_rejected_out_of_slab":                0,
    "n_rejected_soluble_only_contact":       0,
    "membrane_mask_atom_count":              0,
    "membrane_mask_residue_count":           0,
    "slab_z_min":                            None,
    "slab_z_max":                            None,
    "dpp_before":                            0,
    "dpp_after":                             0,
    "annular_occupancy_before":              None,
    "annular_occupancy_after":               None,
    "safety_limit":                          MAX_PRE_EXCLUDED,
    "safety_status":                         "ok",
}}
if PRE_EXCLUDE_TRAPPED:
    try:
        _pre_excl_report = apply_pre_shrink_exclusion(
            gro_path=gro_out,
            lipid_resname=LIPID_RESNAME,
            tm_residues=TM_RESIDUES,
            safety_limit=MAX_PRE_EXCLUDED,
            enabled=True,
        )
        n_excl = _pre_excl_report["n_lipids_removed"]
        if n_excl > 0:
            print(f"[pre_shrink] Removed {{n_excl}} cavity-trapped lipid(s) from system.gro")
            print(f"[pre_shrink] Residue IDs: {{_pre_excl_report['removed_resids'][:10]}}"
                  + (" (+ more)" if n_excl > 10 else ""))
            # Re-renumber residues after removal
            ret2 = subprocess.run(
                ["gmx", "editconf",
                 "-f", str(gro_out), "-o", str(gro_out), "-resnr", "1"],
                capture_output=True, text=True,
            )
            if ret2.returncode != 0:
                print(f"WARNING: gmx editconf -resnr 1 after pre-exclusion failed:\\n{{ret2.stderr}}")
        else:
            print("[pre_shrink] No cavity-trapped lipids found — system.gro unchanged")
    except ValueError as _exc:
        _pre_excl_report["safety_status"] = "exceeded"
        (SCRIPT_DIR / "pre_shrink_lipid_exclusion_report.json").write_text(
            json.dumps(_pre_excl_report, indent=2)
        )
        print(f"ERROR: {{_exc}}", file=sys.stderr)
        sys.exit(1)
(SCRIPT_DIR / "pre_shrink_lipid_exclusion_report.json").write_text(
    json.dumps(_pre_excl_report, indent=2)
)

# ── TM-aware controlled lipid exclusion ────────────────────────────────────────
_excl_report = {{
    "enabled":                   TM_AWARE_EXCLUSION,
    "policy":                    TM_AWARE_EXCLUSION_POLICY,
    "input_gro":                 str(gro_out.name),
    "output_gro":                str(gro_out.name),
    "n_total_lipids":            0,
    "n_bulk_lipids":             0,
    "n_soluble_domain_adjacent": 0,
    "n_tm_surface_overlap":      0,
    "n_tm_cavity_trapped":       0,
    "n_lipids_removed":          0,
    "removed_resids":            [],
    "removed_resnames":          [],
    "removed_classes":           [],
    "safety_limit":              TM_AWARE_MAX_REMOVED,
    "safety_status":             "ok",
    "post_exclusion_n_tm_surface_overlap": 0,
    "post_exclusion_n_tm_cavity_trapped": 0,
    "coordinate_file_modified":  False,
    "topology_modified":         False
}}

if TM_AWARE_EXCLUSION:
    if not TM_RESIDUES:
        msg = ("TM-aware controlled exclusion is enabled, but no TM residues are annotated. "
               "Add tm_segments to structural_annotation or disable exclusion.")
        if TM_AWARE_EXCLUSION_POLICY in ("apply", "strict_apply"):
            print(f"ERROR: {{msg}}", file=sys.stderr)
            sys.exit(1)
        else:
            print(f"WARNING: {{msg}}")
    else:
        try:
            from validators.tm_exclusion_mask import exclude_lipids_tm_aware
            _excl_report = exclude_lipids_tm_aware(
                gro_path=gro_out,
                tm_residues=TM_RESIDUES,
                lipid_resname=LIPID_RESNAME,
                policy=TM_AWARE_EXCLUSION_POLICY,
                safety_limit=TM_AWARE_MAX_REMOVED,
                dry_run=False,
            )
            if _excl_report.get("n_lipids_removed", 0) > 0:
                print(f"[tm_aware_exclusion] Removed {{_excl_report['n_lipids_removed']}} lipid(s) from system.gro")
        except ValueError as _exc:
            print(f"ERROR: {{_exc}}", file=sys.stderr)
            sys.exit(1)
else:
    # If off, write a dummy disabled report
    (SCRIPT_DIR / "tm_aware_lipid_exclusion_report.json").write_text(json.dumps(_excl_report, indent=2))

# ── Overlap + alignment gate ──────────────────────────────────────────────────
ov = validate_no_overlap(gro_out, lipid_residue_name=LIPID_RESNAME)

_align_warnings = []
if not m["tm_annotation_used"]:
    _align_warnings.append(
        "TM annotation absent — bilayer midplane aligned to full protein Z centre "
        "(suboptimal for proteins with large EC/IC domains). "
        "Provide transmembrane_segments in structural_annotation for accurate TM placement."
    )

# Initial lipid–protein clashes are expected at embedding stage — the
# membrane_embedding shrink loop is responsible for resolving them.
_overlap_warnings = (
    [] if ov.n_clashes == 0
    else [f"Initial embed: {{ov.message}} — expected; membrane_embedding shrink loop will resolve"]
)

# ── Lipid-in-cavity detection (reporting only, no remediation) ────────────────
cav = detect_lipid_in_cavity(gro_out, lipid_resname=LIPID_RESNAME, tm_residues=TM_RESIDUES)
_cavity_warnings = []
if cav.n_lipid_residues_in_cavity > 0:
    _cavity_warnings.append(
        f"[lipid_cavity] {{cav.message}} — "
        "visual inspection recommended; membrane_embedding shrink loop may not "
        "fully resolve lipids trapped inside the protein TM channel. "
        f"Affected lipid residue IDs: {{cav.lipid_residue_ids[:10]}}"
        + (" (+ more)" if cav.n_lipid_residues_in_cavity > 10 else "")
    )
    print(f"[cavity_advisory] {{cav.message}}")
else:
    print(f"[cavity_check] {{cav.message}}")

# ── TM-aware mask validation gate ─────────────────────────────────────────────
_tm_errors = []
_tm_warnings = []
if TM_MASK_VALIDATION:
    if not TM_RESIDUES:
        msg = ("TM-aware mask validation is enabled, but no TM residues are annotated. "
               "Add tm_segments to structural_annotation or set membrane.embedding.tm_mask_validation to false.")
        if TM_MASK_POLICY == "strict":
            _tm_errors.append(msg)
        else:
            _tm_warnings.append(msg)
    else:
        try:
            from validators.tm_exclusion_mask import classify_lipids_tm_aware
            res = classify_lipids_tm_aware(
                gro_path=gro_out,
                tm_residues=TM_RESIDUES,
                lipid_resname=LIPID_RESNAME,
            )
            report_data = res["classification_report"]
            counts = report_data["counts"]
            n_bulk = counts.get("bulk_lipid", 0)
            n_surface = counts.get("tm_surface_overlap", 0)
            n_cavity = counts.get("tm_cavity_trapped", 0)
            n_soluble = counts.get("soluble_domain_adjacent", 0)
            n_total = report_data["total_lipids_checked"]
            validation_passed = (n_surface == 0 and n_cavity == 0)

            # Update the JSON report with the summary fields
            report_file = SCRIPT_DIR / "tm_lipid_classification_report.json"
            report_data.update({{
                "n_bulk_lipids": n_bulk,
                "n_tm_surface_overlap": n_surface,
                "n_tm_cavity_trapped": n_cavity,
                "n_soluble_domain_adjacent": n_soluble,
                "n_total_lipids": n_total,
                "tm_mask_policy": TM_MASK_POLICY,
                "validation_passed": validation_passed,
            }})
            report_file.write_text(json.dumps(report_data, indent=2))

            if not validation_passed:
                msg = (f"TM-aware mask validation failed: found {{n_surface}} surface overlaps "
                       f"and {{n_cavity}} cavity-trapped lipids.")
                if TM_MASK_POLICY == "strict":
                    _tm_errors.append(msg)
                else:
                    _tm_warnings.append(msg)
                print(f"[tm_mask_validation] WARNING: {{msg}}")
            else:
                print("[tm_mask_validation] Passed: no suspicious lipids found.")
        except Exception as e:
            msg = f"TM-aware mask validation failed with exception: {{e}}"
            if TM_MASK_POLICY == "strict":
                _tm_errors.append(msg)
            else:
                _tm_warnings.append(msg)
            print(f"[tm_mask_validation] ERROR: {{msg}}")

ov_report = {{
    "passed":                   ov.n_clashes == 0 and (len(_tm_errors) == 0 if TM_MASK_VALIDATION else True),
    "initial_clashes_expected": True,
    "n_clashes":                ov.n_clashes,
    "n_protein_atoms":          ov.n_protein_atoms,
    "n_lipid_atoms":            ov.n_lipid_atoms,
    "message":                  ov.message,
    "tm_center_z":              m["tm_center_z"],
    "protein_center_z":         m["protein_center_z"],
    "bilayer_midplane_z":       m["bilayer_midplane_z"],
    "alignment_method":         m["alignment_method"],
    "tm_annotation_used":       m["tm_annotation_used"],
    "lipid_in_cavity": {{
        "n_lipid_residues_in_cavity": cav.n_lipid_residues_in_cavity,
        "lipid_residue_ids":          cav.lipid_residue_ids,
        "cavity_center_x":            cav.cavity_center_x,
        "cavity_center_y":            cav.cavity_center_y,
        "cavity_z_min":               cav.cavity_z_min,
        "cavity_z_max":               cav.cavity_z_max,
        "cavity_radius_nm":           cav.cavity_radius_nm,
        "message":                    cav.message,
    }},
    "errors":                   _tm_errors,
    "warnings":                 _overlap_warnings + _align_warnings + _cavity_warnings + _tm_warnings,
    "confidence":               1.0,
}}
(SCRIPT_DIR / "overlap_report.json").write_text(json.dumps(ov_report, indent=2))
print(f"[overlap_gate] {{ov.message}}")
for _w in _overlap_warnings + _align_warnings + _tm_warnings:
    print(f"[overlap_advisory] {{_w[:120]}}")
"""

        (step_dir / "run_embed.py").write_text(script)

        # ── Phase 10A: surface-guided lipid refill script ─────────────────────
        simforge_root_embed = str(Path(__file__).resolve().parent.parent.parent)
        self._write_refill_script(
            step_dir           = step_dir,
            assets_ref         = assets_ref,
            bilayer_file       = bilayer_file,
            lipid_resname      = lipid_resname,
            tm_residues_str    = tm_residues_str,
            enabled            = iface_b_enabled,
            policy             = iface_b_policy,
            max_inserted       = iface_b_max_ins,
            max_per_cluster    = iface_b_max_per_cl,
            target_dist        = iface_b_target_dist,
            prot_clash         = iface_b_prot_clash,
            lip_clash          = iface_b_lip_clash,
            seed               = iface_b_seed,
            simforge_root      = simforge_root_embed,
        )

        (step_dir / "run.sh").write_text(
            "#!/bin/bash\n"
            "# ─── embed_in_bilayer (automatic) ──────────────────────────────────\n"
            'python3 "$(dirname "$0")/run_embed.py"\n'
            "# ── Phase 10A: surface-guided lipid refill ──────────────────────────\n"
            'python3 "$(dirname "$0")/run_interface_refill.py"\n'
        )
        meta_params = {
            "bilayer_file":              bilayer_file,
            "lipid":                     lipid,
            "lipid_residue_name":        lipid_resname,
            "tm_residues":               tm_residues_str,
            "pre_exclude_trapped_lipids": pre_exclude_trapped,
            "max_pre_excluded_lipids":    max_pre_excluded,
            "tm_mask_validation":         tm_mask_validation,
            "tm_mask_policy":             tm_mask_policy,
            "tm_aware_exclusion":         tm_aware_exclusion,
            "tm_aware_exclusion_policy":  tm_aware_exclusion_policy,
            "tm_aware_max_removed_lipids": tm_aware_max_removed,
            "iface_builder_enabled":      iface_b_enabled,
            "iface_builder_policy":       iface_b_policy,
        }
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":          step.step_id,
            "stage":            step.stage.value,
            "engine":           step.engine,
            "step_type":        "automatic",
            "automation_level": "automated",
            "generated_by":     "AssemblyBuilder",
            "gate":             {"type": "overlap_report"},
            "expected_outputs": [
                "system.gro",
                "strong_posre_manifest.json",
                "overlap_report.json",
                "pre_shrink_lipid_exclusion_report.json",
                "tm_mask_report.json",
                "tm_lipid_classification_report.json",
                "tm_aware_lipid_exclusion_report.json",
                "interface_lipid_refill_report.json",
            ],
            "params":           meta_params,
        }, indent=4))

    # ── Phase 10A: lipid refill script ───────────────────────────────────────

    def _write_refill_script(
        self,
        step_dir:        Path,
        assets_ref:      str,
        bilayer_file:    str,
        lipid_resname:   str,
        tm_residues_str: str | None,
        enabled:         bool,
        policy:          str,
        max_inserted:    int,
        max_per_cluster: int,
        target_dist:     float,
        prot_clash:      float,
        lip_clash:       float,
        seed:            int,
        simforge_root:   str,
    ) -> None:
        """Write run_interface_refill.py into the embed_in_bilayer step directory."""
        if tm_residues_str:
            from core.structural_annotation import residues_in_range
            tm_set   = sorted(residues_in_range(tm_residues_str))
            tm_value = f"set({tm_set!r})"
        else:
            tm_value = "set()"

        script = f'''#!/usr/bin/env python3
"""Phase 10A Surface-Guided Lipid Refill — generated by AssemblyBuilder."""
import sys
from pathlib import Path

sys.path.insert(0, {simforge_root!r})

from validators.lipid_refill import run_lipid_refill

SCRIPT_DIR  = Path(__file__).parent.resolve()
ASSETS_DIR  = (SCRIPT_DIR / {assets_ref!r}).resolve()

SYSTEM_GRO   = SCRIPT_DIR / "system.gro"
TEMPLATE_GRO = ASSETS_DIR / {bilayer_file!r}
OUTPUT_DIR   = SCRIPT_DIR

TM_RESIDUES     = {tm_value}
LIPID_RESNAME   = {lipid_resname!r}

ENABLED                   = {enabled!r}
POLICY                    = {policy!r}
MAX_INSERTED_LIPIDS       = {max_inserted!r}
MAX_LIPIDS_PER_CLUSTER    = {max_per_cluster!r}
TARGET_CONTACT_DIST_NM    = {target_dist!r}
PROTEIN_CLASH_CUTOFF_NM   = {prot_clash!r}
LIPID_CLASH_CUTOFF_NM     = {lip_clash!r}
DETERMINISTIC_SEED        = {seed!r}

def main():
    if not SYSTEM_GRO.exists():
        print("[refill] system.gro not found — skipping lipid refill", file=sys.stderr)
        return

    try:
        report = run_lipid_refill(
            system_gro                 = SYSTEM_GRO,
            template_gro               = TEMPLATE_GRO,
            tm_residues                = TM_RESIDUES,
            output_gro                 = SYSTEM_GRO,   # in-place update
            output_dir                 = OUTPUT_DIR,
            lipid_resname              = LIPID_RESNAME,
            enabled                    = ENABLED,
            policy                     = POLICY,
            max_inserted_lipids        = MAX_INSERTED_LIPIDS,
            max_lipids_per_gap_cluster = MAX_LIPIDS_PER_CLUSTER,
            target_contact_distance_nm = TARGET_CONTACT_DIST_NM,
            protein_clash_cutoff_nm    = PROTEIN_CLASH_CUTOFF_NM,
            lipid_clash_cutoff_nm      = LIPID_CLASH_CUTOFF_NM,
            deterministic_seed         = DETERMINISTIC_SEED,
            simforge_root              = {simforge_root!r},
        )
        n = report.get("n_inserted_lipids", 0)
        fc_b = report.get("fraction_covered_before", "N/A")
        fc_a = report.get("fraction_covered_after",  "N/A")
        if not ENABLED:
            print("[refill] Phase 10A lipid refill disabled.")
        elif n > 0:
            print(f"[refill] Inserted {{n}} lipid(s) | "
                  f"fraction_covered: {{fc_b}} → {{fc_a}}")
        else:
            print(f"[refill] No lipids inserted | fraction_covered: {{fc_b}}")
        for w in report.get("warning_messages", []):
            print(f"[refill] WARNING: {{w}}")
    except RuntimeError as _e:
        print(f"[refill] STRICT BLOCK: {{_e}}", file=sys.stderr)
        sys.exit(1)
    except Exception as _e:
        print(f"[refill] ERROR: {{_e}}", file=sys.stderr)
        import json as _json
        (OUTPUT_DIR / "interface_lipid_refill_report.json").write_text(
            _json.dumps({{"enabled": ENABLED, "error": str(_e), "n_inserted_lipids": 0}})
        )

if __name__ == "__main__":
    main()
'''
        path = step_dir / "run_interface_refill.py"
        path.write_text(script)
        path.chmod(0o755)

    # ── assemble_system_topology ──────────────────────────────────────────────

    def _build_assemble_system_topology(
        self,
        step:         SimulationStep,
        step_dir:     Path,
        step_dir_map: dict[str, Path],
    ) -> None:
        """
        Assembles the system topol.top from:
        - generate_protein_topology/topol.top  (protein, from pdb2gmx on original PDB)
        - oplsaa_membrane.ff/<lipid>.itp        (lipid moleculetype, static pre-built asset)
        - membrane_embedding/converged.gro      (post-shrink molecule counts)

        Runs AFTER membrane_embedding. Bootstrap mode generates a temporary topology
        at the start of the shrink loop; the final run uses converged.gro counts.
        """
        from builders.step_builders._utils import rel as _rel

        # EMBED_DIR: coordinate source (membrane_embedding/converged.gro after shrink loop)
        # EMBED_BILAYER_DIR: structural assets from embed_in_bilayer (strong_posre_manifest.json)
        embed_dir         = (
            step_dir_map.get("membrane_embedding")
            or step_dir_map.get("embed_in_bilayer")
        )
        embed_bilayer_dir = step_dir_map.get("embed_in_bilayer")
        prot_top_dir      = step_dir_map.get("generate_protein_topology")
        membrane_assets_dir = step_dir_map.get("__membrane_assets__")

        embed_ref        = _rel(step_dir, embed_dir)          if embed_dir          else "../membrane_embedding"
        embed_bilayer_ref = _rel(step_dir, embed_bilayer_dir) if embed_bilayer_dir  else "../embed_in_bilayer"
        prot_top_ref     = _rel(step_dir, prot_top_dir)       if prot_top_dir       else "../generate_protein_topology"
        assets_ref       = _rel(step_dir, membrane_assets_dir) if membrane_assets_dir else "../../membrane_assets"

        forcefield  = step.params.get("forcefield", "opls-aa-membrane")
        from builders.step_builders.preparation_builder import _FF_GROMACS_NAME
        ff_gmx = _FF_GROMACS_NAME.get(forcefield, forcefield)

        script_lines = [
            "#!/usr/bin/env python3",
            "# ─── assemble_system_topology: builds topol.top from converged.gro ─────────",
            "# Runs AFTER membrane_embedding; counts come from membrane_embedding/converged.gro.",
            "# Bootstrap mode (argv[1]=GRO, argv[2]=OUT_TOP): called by membrane_embedding",
            "# at the start of the shrink loop to generate a temporary topology.",
            "import sys, subprocess, shutil, json, re",
            "from pathlib import Path",
            "",
            "SCRIPT_DIR       = Path(__file__).parent.resolve()",
            f'EMBED_DIR        = (SCRIPT_DIR / "{embed_ref}").resolve()',
            f'EMBED_BILAYER_DIR = (SCRIPT_DIR / "{embed_bilayer_ref}").resolve()',
            f'PROT_TOP_DIR     = (SCRIPT_DIR / "{prot_top_ref}").resolve()',
            f'ASSETS_DIR       = (SCRIPT_DIR / "{assets_ref}").resolve()',
            "",
            "# Bootstrap support: argv[1]=input GRO override, argv[2]=output topology path",
            "_gro_override = sys.argv[1] if len(sys.argv) > 1 else None",
            "_top_override = sys.argv[2] if len(sys.argv) > 2 else None",
            "_out_top = Path(_top_override).resolve() if _top_override else SCRIPT_DIR / 'topol.top'",
            "_out_dir = _out_top.parent",
            "",
            "import os as _os",
            "_gmx_env = dict(_os.environ)",
            '_gmx_env["GMXLIB"] = str(ASSETS_DIR) + (_os.pathsep + _gmx_env["GMXLIB"] if "GMXLIB" in _gmx_env else "")',
            "",
            "_LIPID   = {'DPP', 'DPPC', 'POPC', 'POPE', 'POPG', 'POPS', 'CHOL',",
            "            'DLPC', 'DMPC', 'DOPC', 'PLPC', 'LYPC'}",
            "_SOLVENT = {'SOL', 'HOH', 'WAT', 'TIP3', 'TIP4', 'TIP5'}",
            "_IONS    = {'NA', 'CL', 'SOD', 'MG', 'K', 'CA', 'ZN'}",
            "",
            "# ── Step 1: Resolve input coordinate file ────────────────────────────────",
            "if _gro_override:",
            "    FINAL_GRO = Path(_gro_override).resolve()",
            "else:",
            "    FINAL_GRO = EMBED_DIR / 'converged.gro'",
            "    if not FINAL_GRO.exists():",
            "        FINAL_GRO = EMBED_DIR / 'system.gro'",
            "if not FINAL_GRO.exists():",
            "    print(f'[assemble_system_topology] ERROR: coordinate file not found at {FINAL_GRO}', file=sys.stderr)",
            "    sys.exit(1)",
            "",
            "# system_processed.gro: written in final mode (non-bootstrap) for backward compat",
            "if not _top_override:",
            "    shutil.copy2(FINAL_GRO, SCRIPT_DIR / 'system_processed.gro')",
            "",
            "gro_lines = FINAL_GRO.read_text().splitlines()",
            "n_atoms   = int(gro_lines[1].strip())",
            "atom_lines = gro_lines[2: 2 + n_atoms]",
            "",
            "# ── Step 2: Count residues ────────────────────────────────────────────────",
            "residue_counts: dict = {}",
            "prev_key = None",
            "for ln in atom_lines:",
            "    if len(ln) < 10: continue",
            "    resnum  = int(ln[0:5])",
            "    resname = ln[5:10].strip()",
            "    key = (resnum % 100000, resname)",
            "    if key != prev_key:",
            "        residue_counts[resname] = residue_counts.get(resname, 0) + 1",
            "        prev_key = key",
            "",
            "lipid_counts   = {k: v for k, v in residue_counts.items() if k in _LIPID}",
            "solvent_counts = {k: v for k, v in residue_counts.items() if k in _SOLVENT}",
            "ion_counts     = {k: v for k, v in residue_counts.items() if k in _IONS}",
            "",
            "# ── Step 3: Read protein topology from manifest ──────────────────────────",
            "_manifest_f = PROT_TOP_DIR / 'protein_topology_manifest.json'",
            "if _manifest_f.exists():",
            "    _prot_manifest     = json.loads(_manifest_f.read_text())",
            "    _molecule_itps_p   = _prot_manifest.get('molecule_itps', [])",
            "    _ff_incs_p         = _prot_manifest.get('forcefield_includes', [])",
            "    _mol_names_p       = _prot_manifest.get('molecule_names', [])",
            "    # Priority 1: molecule_entries (names + authoritative counts from pdb2gmx [ molecules ])",
            "    _mol_entries_p = _prot_manifest.get('molecule_entries', [])",
            "    if _mol_entries_p:",
            "        protein_mol_lines = [f\"{me['name']} {me.get('count', 1)}\" for me in _mol_entries_p]",
            "    else:",
            "        # Priority 2: molecule_names (just names, count=1 each)",
            "        protein_mol_lines = [f'{n} 1' for n in _mol_names_p]",
            "else:",
            "    # Fallback: parse topol.top directly (single-chain or legacy workspace)",
            "    _molecule_itps_p, _ff_incs_p, _mol_names_p = [], [], []",
            "    _pt_text = (PROT_TOP_DIR / 'topol.top').read_text()",
            "    for _ln in _pt_text.splitlines():",
            "        _s = _ln.lstrip()",
            "        if _s.startswith('#include'):",
            "            _pts = _s.split('\"')",
            "            _inc = _pts[1] if len(_pts) >= 2 else ''",
            "            if '.ff/' in _inc or '.ff\"' in _inc or _inc.endswith('.ff'):",
            "                _ff_incs_p.append(_inc)",
            "            elif _inc.endswith('.itp') and 'posre' not in _inc.lower():",
            "                _molecule_itps_p.append(_inc)",
            "    _mol_m = re.search(r'^\\s*\\[\\s*molecules\\s*\\](.*)',",
            "                       _pt_text, re.MULTILINE | re.IGNORECASE | re.DOTALL)",
            "    protein_mol_lines = []",
            "    if _mol_m:",
            "        for _ml in _mol_m.group(1).splitlines():",
            "            _ml = _ml.strip()",
            "            if _ml and not _ml.startswith(';'):",
            "                _ps = _ml.split()",
            "                if _ps:",
            "                    protein_mol_lines.append(f'{_ps[0]} {_ps[1] if len(_ps) > 1 else 1}')",
            "",
            "# Priority 3: parse [ moleculetype ] from each included ITP (robust fallback)",
            "# Used when manifest has no molecule_entries/molecule_names, or as a cross-check.",
            "if not protein_mol_lines and _molecule_itps_p:",
            "    _in_mt = False",
            "    for _itp_fb in _molecule_itps_p:",
            "        _itp_fb_abs = (Path(_itp_fb) if Path(_itp_fb).is_absolute() else (PROT_TOP_DIR / _itp_fb)).resolve()",
            "        if not _itp_fb_abs.exists(): continue",
            "        _in_mt = False",
            "        for _l in _itp_fb_abs.read_text().splitlines():",
            "            _ls = _l.lstrip()",
            "            if re.match(r'^\\s*\\[\\s*moleculetype\\s*\\]', _ls, re.IGNORECASE):",
            "                _in_mt = True; continue",
            "            if _in_mt and _ls and not _ls.startswith(';') and not _ls.startswith('['):",
            "                protein_mol_lines.append(f'{_ls.split()[0]} 1'); _in_mt = False",
            "            elif _in_mt and _ls.startswith('['):",
            "                _in_mt = False",
            "if protein_mol_lines:",
            "    print(f'[assemble_system_topology] Protein molecules: {protein_mol_lines}')",
            "else:",
            "    print('[assemble_system_topology] WARNING: no protein molecule entries found', file=sys.stderr)",
            "",
            "# ── Step 4: Build protein include block with correct relative paths ────────",
            "# _out_dir: directory that will contain the output topology (may differ from",
            "# SCRIPT_DIR in bootstrap mode). All #include paths are relative to _out_dir",
            "# so grompp finds them when parsing the topology from that location.",
            "_prot_dir_rel  = _os.path.relpath(str(PROT_TOP_DIR), str(_out_dir))",
            "_assets_ff_rel = _os.path.relpath(str(ASSETS_DIR / 'oplsaa_membrane.ff'), str(_out_dir))",
            "_prot_block = []",
            "# Always use a filesystem-relative path so grompp finds the FF without GMXLIB",
            "_prot_block.append(f'#include \"{_assets_ff_rel}/forcefield.itp\"')",
            "_prot_block.append('')",
            "# Read strong_posre_manifest from embed_in_bilayer (maps mol_name → per-chain file)",
            "_sp_mf = EMBED_BILAYER_DIR / 'strong_posre_manifest.json'",
            "_sp_manifest_asm = json.loads(_sp_mf.read_text()) if _sp_mf.exists() else {}",
            "_sp_per_chain_done = set()  # mol names that got per-chain includes injected",
            "for _i_itp, _itp in enumerate(_molecule_itps_p):",
            "    # Skip any FF file that leaked into molecule_itps (e.g. absolute-path forcefield.itp",
            "    # from pdb2gmx that _walk_includes resolved locally). We already have the canonical FF include.",
            "    if '.ff/' in _itp or _itp.endswith('.ff'):",
            "        continue",
            "    # Normalize: manifest may contain absolute or relative paths — always rewrite",
            "    _itp_abs = (Path(_itp) if Path(_itp).is_absolute() else (PROT_TOP_DIR / _itp)).resolve()",
            "    _itp_rel = _os.path.relpath(str(_itp_abs), str(_out_dir))",
            "    _prot_block.append(f'#include \"{_itp_rel}\"')",
            "    # Inject per-chain strong posre immediately after this chain ITP include",
            "    _mn_i = _mol_names_p[_i_itp] if _i_itp < len(_mol_names_p) else None",
            "    if _mn_i and _mn_i in _sp_manifest_asm:",
            "        _sp_src_i = EMBED_BILAYER_DIR / _sp_manifest_asm[_mn_i]",
            "        if _sp_src_i.exists():",
            "            shutil.copy2(_sp_src_i, _out_dir / _sp_manifest_asm[_mn_i])",
            "            _prot_block += ['#ifdef STRONG_POSRES',",
            "                            f'#include \"{_sp_manifest_asm[_mn_i]}\"',",
            "                            '#endif']",
            "            _sp_per_chain_done.add(_mn_i)",
            "# When molecule_itps is empty (single-chain: pdb2gmx inlines topology into topol.top),",
            "# embed cleaned protein topology content so grompp finds the [ moleculetype ] definition.",
            "if not _molecule_itps_p:",
            "    _prot_top_f = PROT_TOP_DIR / 'topol.top'",
            "    if _prot_top_f.exists():",
            "        _in_skip_sec = False",
            "        _pt_out = []",
            "        for _ptl in _prot_top_f.read_text().splitlines():",
            "            _sec_m = re.match(r'^\\s*\\[\\s*(\\w+)\\s*\\]', _ptl)",
            "            if _sec_m:",
            "                _in_skip_sec = _sec_m.group(1).lower() in ('defaults', 'system', 'molecules')",
            "            if _in_skip_sec: continue",
            "            _ptls = _ptl.lstrip()",
            "            if _ptls.startswith('#include') and ('.ff/' in _ptls or _ptls.rstrip().endswith('.ff\"')):",
            "                continue  # skip FF includes — canonical FF include already in prot_block",
            "            if _ptls.startswith('#include'):",
            "                _inc_q = _ptls.split('\"')",
            "                if len(_inc_q) >= 2 and not _os.path.isabs(_inc_q[1]):",
            "                    _inc_abs = (PROT_TOP_DIR / _inc_q[1]).resolve()",
            "                    _inc_rel = _os.path.relpath(str(_inc_abs), str(SCRIPT_DIR))",
            "                    _ptl = f'#include \"{_inc_rel}\"'",
            "            _pt_out.append(_ptl)",
            "        _pt_clean_str = '\\n'.join(_pt_out).strip()",
            "        if _pt_clean_str:",
            "            _prot_block.append('')",
            "            _prot_block.append(_pt_clean_str)",
            "protein_section = '\\n'.join(_prot_block)",
            "",
            "_strong_posre_injected = False",
            "# single-chain fallback: copy strong_posre.itp only when no per-chain files were injected",
            "if not _sp_per_chain_done:",
            "    for _spsrc in [EMBED_BILAYER_DIR / 'strong_posre.itp', _out_dir / 'strong_posre.itp']:",
            "        if _spsrc.exists() and _spsrc != _out_dir / 'strong_posre.itp':",
            "            shutil.copy2(_spsrc, _out_dir / 'strong_posre.itp')",
            "            _strong_posre_injected = True",
            "            break",
            "        elif (_spsrc == _out_dir / 'strong_posre.itp') and _spsrc.exists():",
            "            _strong_posre_injected = True",
            "            break",
            "",
            "# ── Step 5: Resolve lipid topology from pre-built ITP asset ────────────────",
            "# DPP (DPPC) and other OPLS-AA membrane lipids are defined only in aminoacids.rtp",
            "# with no standalone pdb2gmx-compatible ITP. Pre-built ITP files are shipped in",
            "# oplsaa_membrane.ff/ (generated from RTP by core/lipid_itp_builder.py) and",
            "# auto-staged to membrane_assets/ at workspace build time.",
            "_LIPID_ITP_MAP = {'DPP': 'dpp.itp', 'DPPC': 'dpp.itp'}",
            "_lipid_section = ''",
            "if lipid_counts:",
            "    _lip_itp_parts = []",
            "    for _res in lipid_counts:",
            "        _itp_name = _LIPID_ITP_MAP.get(_res)",
            "        if _itp_name:",
            "            _lip_itp_path = ASSETS_DIR / 'oplsaa_membrane.ff' / _itp_name",
            "            if _lip_itp_path.exists():",
            "                _lip_itp_rel = _os.path.relpath(str(_lip_itp_path), str(_out_dir))",
            "                _lip_itp_parts.append(f'#include \"{_lip_itp_rel}\"')",
            "            else:",
            "                print(f'[assemble_system_topology] WARNING: lipid ITP not found: {_lip_itp_path}',",
            "                      file=sys.stderr)",
            "        else:",
            "            print(f'[assemble_system_topology] WARNING: no ITP mapping for lipid {_res!r}',",
            "                  file=sys.stderr)",
            "    if _lip_itp_parts:",
            "        _lipid_section = '\\n'.join(_lip_itp_parts)",
            "",
            "# ── Step 6: Assemble topology (counts from final coordinate file) ─────────",
            "out_parts = [protein_section]",
            "if _lipid_section:",
            "    out_parts += ['', '; ─── Lipid topology ────────────────────────────────────────────────', _lipid_section]",
            "out_parts += [",
            "    '',",
            "    '[ system ]',",
            "    'Protein-membrane system',",
            "    '',",
            "    '[ molecules ]',",
            "    '; Counts derived from final coordinate file (membrane_embedding/converged.gro)',",
            "]",
            "for _ml in protein_mol_lines:",
            "    out_parts.append(_ml)",
            "for _res, _cnt in lipid_counts.items():",
            "    out_parts.append(f'{_res:<20s} {_cnt}')",
            "for _res, _cnt in solvent_counts.items():",
            "    out_parts.append(f'{_res:<20s} {_cnt}')",
            "for _res, _cnt in ion_counts.items():",
            "    out_parts.append(f'{_res:<20s} {_cnt}')",
            "topol_text = '\\n'.join(out_parts) + '\\n'",
            "",
            "if not _sp_per_chain_done and _strong_posre_injected and 'strong_posre.itp' not in topol_text:",
            "    _tlines  = topol_text.splitlines()",
            "    _last_pi = -1",
            "    for _ti, _tl in enumerate(_tlines):",
            "        if '#include' in _tl and _prot_dir_rel in _tl:",
            "            _last_pi = _ti",
            "    _posre_block = ['; Strong position restraints',",
            "                    '#ifdef STRONG_POSRES', '#include \"strong_posre.itp\"', '#endif']",
            "    if _last_pi >= 0:",
            "        _tlines = _tlines[:_last_pi+1] + _posre_block + _tlines[_last_pi+1:]",
            "    else:",
            "        _tlines += _posre_block",
            "    topol_text = '\\n'.join(_tlines) + '\\n'",
            "",
            "_out_top.write_text(topol_text)",
            "print(f'[assemble_system_topology] Topology written to {_out_top}')",
            "",
            "# ── Moleculetype consistency helpers (used in bootstrap and normal mode) ───",
            "def _collect_moleculetypes(f, visited=None):",
            "    if visited is None: visited = set()",
            "    key = str(Path(f).resolve())",
            "    if key in visited: return []",
            "    visited.add(key)",
            "    result = []",
            "    try: txt = Path(f).read_text()",
            "    except OSError: return result",
            "    _in_mt = False",
            "    for _l in txt.splitlines():",
            "        _ls = _l.lstrip()",
            "        if re.match(r'^\\[\\s*moleculetype\\s*\\]', _ls, re.IGNORECASE):",
            "            _in_mt = True; continue",
            "        if _in_mt:",
            "            if _ls and not _ls.startswith(';') and not _ls.startswith('['):",
            "                result.append(_ls.split()[0]); _in_mt = False",
            "            elif _ls.startswith('['):",
            "                _in_mt = False",
            "        if _ls.startswith('#include'):",
            "            _pts = _ls.split('\"')",
            "            if len(_pts) >= 2:",
            "                _cand = (Path(f).parent / _pts[1]).resolve()",
            "                if _cand.exists(): result += _collect_moleculetypes(str(_cand), visited)",
            "    return result",
            "",
            "def _check_moltype_consistency(top_path):",
            "    _defined = set(_collect_moleculetypes(str(top_path)))",
            "    _txt = top_path.read_text()",
            "    _mol_m = re.search(r'^\\s*\\[\\s*molecules\\s*\\](.*)', _txt,",
            "                       re.MULTILINE | re.IGNORECASE | re.DOTALL)",
            "    _mol_names = []",
            "    if _mol_m:",
            "        for _mbl in _mol_m.group(1).splitlines():",
            "            _mbl = _mbl.strip()",
            "            if _mbl and not _mbl.startswith(';') and not _mbl.startswith('#'):",
            "                _mbs = _mbl.split()",
            "                if _mbs: _mol_names.append(_mbs[0])",
            "    _undef = [n for n in _mol_names if n not in _defined]",
            "    return _defined, _mol_names, _undef",
            "",
            "# Bootstrap mode: validate moleculetype consistency before grompp runs",
            "if _top_override:",
            "    _def_bs, _names_bs, _undef_bs = _check_moltype_consistency(_out_top)",
            "    if _undef_bs:",
            "        for _ub in _undef_bs:",
            "            print(f'TOPOLOGY_ASSEMBLY_ERROR: Molecule {_ub!r} in [ molecules ] '",
            "                  f'but no [ moleculetype ] {_ub!r} defined.', file=sys.stderr)",
            "        sys.exit(1)",
            "    print(f'[assemble_system_topology] Bootstrap moleculetype check OK '",
            "          f'({len(_def_bs)} types, {len(_names_bs)} molecules)')",
            "    sys.exit(0)",
            "",
            "# ── Step 6b: Validate topology include graph ─────────────────────────────",
            "def _chk_inc(f, seen=None):",
            "    seen = seen or set()",
            "    key = str(Path(f).resolve())",
            "    if key in seen: return []",
            "    seen.add(key)",
            "    bad = []",
            "    try: txt = Path(f).read_text()",
            "    except OSError: return [f'Cannot read: {f}']",
            "    for _l in txt.splitlines():",
            "        _s = _l.lstrip()",
            "        if not _s.startswith('#include'): continue",
            "        _pts = _s.split('\"')",
            "        if len(_pts) < 2: continue",
            "        _inc = _pts[1]",
            "        _cand = (Path(f).parent / _inc).resolve()",
            "        if _cand.exists(): bad += _chk_inc(str(_cand), seen)",
            "        else: bad.append(f'Cannot resolve {_inc!r} (in {Path(f).name})')",
            "    return bad",
            "_inc_errors = _chk_inc(str(_out_top))",
            "if _inc_errors:",
            "    print('[assemble_system_topology] ERROR: broken topology includes:', file=sys.stderr)",
            "    for _e in _inc_errors: print(f'  {_e}', file=sys.stderr)",
            "    sys.exit(1)",
            "print('[assemble_system_topology] Topology include graph OK')",
            "",
            "# ── Step 6c: Validate exactly one [ defaults ] in expanded topology ─────────",
            "# Count like GROMACS (no dedup across siblings — only cycle-protect ancestors).",
            "def _count_defaults(f, _anc=None):",
            "    if _anc is None: _anc = set()",
            "    key = str(Path(f).resolve())",
            "    if key in _anc: return 0  # cycle guard only",
            "    _child_anc = _anc | {key}",
            "    n = 0",
            "    try: txt = Path(f).read_text()",
            "    except OSError: return 0",
            "    for _l in txt.splitlines():",
            "        if re.match(r'^\\s*\\[\\s*defaults\\s*\\]', _l, re.IGNORECASE): n += 1",
            "        _ls = _l.lstrip()",
            "        if _ls.startswith('#include'):",
            "            _pts = _ls.split('\"')",
            "            if len(_pts) >= 2:",
            "                _cand = (Path(f).parent / _pts[1]).resolve()",
            "                if _cand.exists(): n += _count_defaults(str(_cand), _child_anc)",
            "    return n",
            "_def_count = _count_defaults(str(_out_top))",
            "if _def_count != 1:",
            "    print(f'[assemble_system_topology] TOPOLOGY_DUPLICATE_DEFAULTS_ERROR: '",
            "          f'expanded topology has {_def_count} [ defaults ] section(s) (expected exactly 1).',",
            "          file=sys.stderr)",
            "    print('  Cause: forcefield.itp is included more than once, or a molecule .itp',",
            "          file=sys.stderr)",
            "    print('  includes forcefield.itp itself (standalone topology included wholesale).',",
            "          file=sys.stderr)",
            "    sys.exit(1)",
            "print(f'[assemble_system_topology] [ defaults ] count: {_def_count} (OK)')",
            "",
            "# ── Step 6d: Validate [ moleculetype ] consistency ─────────────────────────",
            "# Every name in [ molecules ] must have a [ moleculetype ] definition.",
            "# Catches 'No such moleculetype DPP' and name mismatches before grompp.",
            "_defined_moltypes, _mol_names_v, _undefined_mols = _check_moltype_consistency(_out_top)",
            "if _undefined_mols:",
            "    _avail = ', '.join(sorted(_defined_moltypes)) or '(none found)'",
            "    for _ut in _undefined_mols:",
            "        print(f'TOPOLOGY_UNDEFINED_MOLECULETYPE_ERROR: [ molecules ] contains \\'{_ut}\\', '",
            "              f'but no included [ moleculetype ] named \\'{_ut}\\' exists.', file=sys.stderr)",
            "    print(f'  Available molecule types include: {_avail}', file=sys.stderr)",
            "    sys.exit(1)",
            "print(f'[assemble_system_topology] [ moleculetype ] check: {len(_defined_moltypes)} types defined, {len(_mol_names_v)} molecule entries — OK')",
            "",
            "# ── Step 6e: Validate [ position_restraints ] index bounds ──────────────",
            "# Every atom index in every included [ position_restraints ] file must be",
            "# <= the local atom count of its enclosing moleculetype.",
            "import sys as _sys_posre",
            "try:",
            "    import importlib.util as _ilu",
            "    _spv_spec = _ilu.spec_from_file_location(",
            "        'strong_posre_validator',",
            "        str(Path(sys.argv[0]).parent.parent.parent / 'validators' / 'strong_posre_validator.py'),",
            "    )",
            "    if _spv_spec is None:",
            "        raise ImportError('strong_posre_validator not found')",
            "    _spv_mod = _ilu.module_from_spec(_spv_spec)",
            "    _spv_spec.loader.exec_module(_spv_mod)",
            "    _posre_errors = _spv_mod.validate_posre_bounds(_out_top)",
            "    if _posre_errors:",
            "        for _pe in _posre_errors:",
            "            print(_pe, file=sys.stderr)",
            "        sys.exit(1)",
            "    print('[assemble_system_topology] [ position_restraints ] bounds check OK')",
            "except (ImportError, OSError) as _e:",
            "    print(f'[assemble_system_topology] WARNING: could not load strong_posre_validator: {_e}', file=sys.stderr)",
            "",
            "# ── Step 7: topology_assembly_report.json ────────────────────────────────",
            "_errors, _warnings = [], []",
            "_top_exists = (SCRIPT_DIR / 'topol.top').exists()",
            "if not _top_exists: _errors.append('topol.top not assembled')",
            "if not _lipid_section and lipid_counts: _warnings.append('lipid topology section is empty')",
            "_report = {",
            "    'passed':         len(_errors) == 0,",
            "    'final_gro':      str(FINAL_GRO),",
            "    'n_atoms':        n_atoms,",
            "    'residue_counts': residue_counts,",
            "    'errors':         _errors,",
            "    'warnings':       _warnings,",
            "}",
            "(SCRIPT_DIR / 'topology_assembly_report.json').write_text(json.dumps(_report, indent=2))",
            "print(f'[assemble_system_topology] OK: {sum(residue_counts.values())} residues, errors={_errors}')",
        ]

        (step_dir / "run_assemble_system.py").write_text("\n".join(script_lines) + "\n")
        (step_dir / "run.sh").write_text(
            "#!/bin/bash\n"
            "# ─── assemble_system_topology: builds topol.top after membrane embedding ──\n"
            "# Reads membrane_embedding/converged.gro for final molecule counts.\n"
            'python3 "$(dirname "$0")/run_assemble_system.py"\n'
        )
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":          step.step_id,
            "stage":            step.stage.value,
            "engine":           step.engine,
            "step_type":        step.step_type.value,
            "expected_outputs": ["topol.top", "system_processed.gro", "topology_assembly_report.json"],
            "required_inputs":  [
                f"{embed_ref}/converged.gro",
                f"{embed_bilayer_ref}/strong_posre_manifest.json",
                f"{prot_top_ref}/topol.top",
                f"{prot_top_ref}/protein_topology_manifest.json",
            ],
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
