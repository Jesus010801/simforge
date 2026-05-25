# builders/step_builders/preparation_builder.py

from __future__ import annotations

from pathlib import Path
import json

from core.execution_models import SimulationStep, StepType
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

        if engine == "gromacs:editconf+orient":
            from builders.step_builders.membrane_orient_builder import MembraneOrientBuilder
            MembraneOrientBuilder().build(step, step_dir, step_dir_map)
        elif step.step_type == StepType.MANUAL:
            self._build_manual_readme(step, step_dir)
        elif engine == "gromacs:pdb2gmx":
            self._build_protein_prep(step, step_dir, inputs_ref, step_dir_map)
        else:
            self._build_ligand_prep(step, step_dir, inputs_ref)

    # ── Proteína ──────────────────────────────────────────────────────────────

    def _build_protein_prep(
        self,
        step:         SimulationStep,
        step_dir:     Path,
        inputs_ref:   str,
        step_dir_map: dict = {},
    ) -> None:

        target      = step.target_components[0] if step.target_components else "protein"
        source_file = step.params.get("source_file")
        forcefield  = step.params.get("forcefield", "charmm36")
        water_model = step.params.get("water_model", "tip3p")
        source_step = step.params.get("source_step")  # e.g. "embed_in_bilayer"

        ff_gmx = _FF_GROMACS_NAME.get(forcefield, forcefield)

        if source_step and step_dir_map and source_step in step_dir_map:
            # Input comes from a prior step's output directory (membrane workflow)
            src_dir        = step_dir_map[source_step]
            src_ref        = _rel(step_dir, src_dir)
            src_filename   = source_file or "system.gro"
            required_input = f"{src_ref}/{src_filename}"
            output_gro     = "system_processed.gro"
            expected_outputs = [output_gro, "topol.top", "posre.itp"]

            embed_ref = (
                _rel(step_dir, step_dir_map["embed_in_bilayer"])
                if "embed_in_bilayer" in step_dir_map
                else "../embed_in_bilayer"
            )
            inject_strong_posre = (source_step == "embed_in_bilayer")
            if inject_strong_posre:
                expected_outputs.append("strong_posre.itp")
            expected_outputs.append("topology_consistency_report.json")

            script_lines = [
                "#!/usr/bin/env python3",
                f"# ─── generate_topology (membrane, input from {source_step}) ──────────────────",
                "import sys, subprocess, shutil, json",
                "from pathlib import Path",
                "",
                "SCRIPT_DIR = Path(__file__).parent.resolve()",
                "",
                "def _find_root(start):",
                "    p = start",
                "    for _ in range(12):",
                "        if (p / 'adapters').is_dir() and (p / 'core').is_dir():",
                "            return p",
                "        p = p.parent",
                "    raise RuntimeError(f'SimForge project root not found from {start}')",
                "",
                "PROJECT_ROOT = _find_root(SCRIPT_DIR)",
                "sys.path.insert(0, str(PROJECT_ROOT))",
                "",
                "ret = subprocess.run(",
                '    ["gmx", "pdb2gmx",',
                f'     "-f", str((SCRIPT_DIR / "{src_ref}" / "{src_filename}").resolve()),',
                f'     "-o", "{output_gro}",',
                '     "-p", "topol.top",',
                f'     "-ff", "{ff_gmx}",',
                f'     "-water", "{water_model}",',
                '     "-ignh"],',
                "    cwd=str(SCRIPT_DIR),",
                "    capture_output=False,",
                ")",
                "if ret.returncode != 0:",
                "    sys.exit(ret.returncode)",
            ]

            if inject_strong_posre:
                script_lines += [
                    "",
                    "# ── Copy strong_posre.itp and inject #ifdef STRONG_POSRES into topol.top ──────",
                    f'posre_src = (SCRIPT_DIR / "{embed_ref}" / "strong_posre.itp").resolve()',
                    "_strong_posre_injected = False",
                    "if posre_src.exists():",
                    '    shutil.copy2(posre_src, SCRIPT_DIR / "strong_posre.itp")',
                    '    top_path = SCRIPT_DIR / "topol.top"',
                    "    lines = top_path.read_text().splitlines()",
                    "    out, inserted = [], False",
                    "    for line in lines:",
                    "        out.append(line)",
                    '        if not inserted and \'#include "posre.itp"\' in line:',
                    "            out.append('; Strong position restraints (shrink loop + EM)')",
                    "            out.append('#ifdef STRONG_POSRES')",
                    "            out.append('#include \"strong_posre.itp\"')",
                    "            out.append('#endif')",
                    "            inserted = True",
                    "    if not inserted:",
                    "        out += ['; Strong position restraints', '#ifdef STRONG_POSRES',",
                    "                '#include \"strong_posre.itp\"', '#endif']",
                    "    top_path.write_text('\\n'.join(out) + '\\n')",
                    "    print('strong_posre.itp copied and topol.top updated (#ifdef STRONG_POSRES)')",
                    "    _strong_posre_injected = True",
                    "else:",
                    "    print(f'WARNING: strong_posre.itp not found at {posre_src} — skipping')",
                ]

            # Topology consistency gate report
            script_lines += [
                "",
                "# ── Topology consistency gate report ─────────────────────────────────────",
                f'_gro_out = SCRIPT_DIR / "{output_gro}"',
                "_gro_atoms = None",
                "if _gro_out.exists():",
                "    try:",
                "        _gro_atoms = int(_gro_out.read_text().splitlines()[1].strip())",
                "    except Exception:",
                "        pass",
                '_top_path = SCRIPT_DIR / "topol.top"',
                "_top_exists = _top_path.exists()",
                '_posre_ok = _top_exists and \'#include "posre.itp"\' in _top_path.read_text() if _top_exists else False',
                "_errors, _warnings = [], []",
                "if not _top_exists:",
                '    _errors.append("topol.top not generated by pdb2gmx")',
                "elif not _posre_ok:",
                '    _errors.append("posre.itp not referenced in topol.top — pdb2gmx may have failed silently")',
            ]
            if inject_strong_posre:
                script_lines += [
                    "if not _strong_posre_injected:",
                    '    _warnings.append("strong_posre.itp not injected into topol.top — check embed_in_bilayer output")',
                ]
            script_lines += [
                "_tc_report = {",
                '    "passed":                 len(_errors) == 0,',
                '    "gro_atoms":              _gro_atoms,',
                '    "top_exists":             _top_exists,',
                '    "posre_included":         _posre_ok,',
            ]
            if inject_strong_posre:
                script_lines.append('    "strong_posre_included":  _strong_posre_injected,')
            script_lines += [
                '    "errors":                 _errors,',
                '    "warnings":               _warnings,',
                '    "confidence":             1.0,',
                "}",
                '(SCRIPT_DIR / "topology_consistency_report.json").write_text(json.dumps(_tc_report, indent=2))',
                'print(f"[topology_gate] errors={_errors}, warnings={_warnings}")',
            ]

            (step_dir / "run_topology.py").write_text("\n".join(script_lines) + "\n")
            (step_dir / "run.sh").write_text(
                "#!/bin/bash\n"
                f"# ─── generate_topology (membrane, source: {source_step}) ────────────\n"
                'python3 "$(dirname "$0")/run_topology.py"\n'
            )

        else:
            src_ext        = Path(source_file).suffix if source_file else ".pdb"
            pdb_name       = f"{target}{src_ext}"
            required_input = f"{inputs_ref}/{pdb_name}"
            output_gro     = f"{target}_processed.gro"
            expected_outputs = [output_gro, "topol.top", "posre.itp"]

            script = (
                f"#!/bin/bash\n"
                f"# ─── Preparación de proteína: {target} ────────────────────────────────────\n"
                f"# El PDB de entrada vive en workspace/inputs/ — workspace auto-contenido.\n\n"
                f'INPUTS_DIR="{inputs_ref}"\n\n'
                f"gmx pdb2gmx \\\\\n"
                f'    -f "$INPUTS_DIR/{pdb_name}" \\\\\n'
                f"    -o {output_gro} \\\\\n"
                f"    -p topol.top \\\\\n"
                f"    -ff {ff_gmx} \\\\\n"
                f"    -water {water_model} \\\\\n"
                f"    -ignh\n"
            )
            (step_dir / "run.sh").write_text(script)

        gate = {"type": "topology_consistency"} if source_step else None
        meta_dict: dict = {
            "step_id":          step.step_id,
            "stage":            step.stage.value,
            "engine":           step.engine,
            "target":           target,
            "step_type":        step.step_type.value,
            "expected_outputs": expected_outputs,
            "required_inputs":  [required_input],
        }
        if gate:
            meta_dict["gate"] = gate
        (step_dir / "metadata.json").write_text(json.dumps(meta_dict, indent=4))

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

    # ── Manual (instrucciones para el usuario) ────────────────────────────────

    def _build_manual_readme(
        self,
        step:     SimulationStep,
        step_dir: Path,
    ) -> None:
        note = step.params.get("note", "Ver documentación del pipeline.")
        notes_lines = "\n".join(f"- {n}" for n in step.notes) if step.notes else "- Sin notas adicionales"

        readme = f"""# {step.title}

## Acción requerida (manual)
Este step requiere intervención del usuario antes de continuar.

## Instrucciones
```
{note}
```

## Notas
{notes_lines}

## Cómo continuar
Una vez completado este step manualmente, ejecuta:
    simforge run <workspace>

El executor detectará los outputs y continuará desde el siguiente step automático.
"""
        (step_dir / "README.md").write_text(readme)
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":          step.step_id,
            "stage":            step.stage.value,
            "engine":           step.engine,
            "step_type":        step.step_type.value,
            "automation_level": "guided",
            "params":           step.params,
        }, indent=4))
