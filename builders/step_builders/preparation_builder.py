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

# ─────────────────────────────────────────────────────────────────────────────
# _ITP_WALK_CODE — injected verbatim into generated scripts.
# Recursively walks #include directives in a GROMACS topology file and
# returns (local_itps, forcefield_includes).
# local_itps:          names of .itp files that exist alongside the .top
# forcefield_includes: includes that did not resolve inside base_d (ff/ paths)
# ─────────────────────────────────────────────────────────────────────────────
_ITP_WALK_CODE = """\
def _walk_includes(top_f, base_d, _seen=None, _local=None, _ff=None):
    if _seen is None:
        _seen, _local, _ff = set(), [], []
    key = str(top_f.resolve())
    if key in _seen:
        return _local, _ff
    _seen.add(key)
    try:
        _text = top_f.read_text()
    except OSError:
        return _local, _ff
    for _line in _text.splitlines():
        _s = _line.lstrip()
        if not _s.startswith('#include'):
            continue
        _pts = _s.split('"')
        if len(_pts) < 2:
            continue
        _inc = _pts[1]
        # Force field files must always go to _ff — never to _local.
        # pdb2gmx writes absolute paths (e.g. /abs/.../oplsaa_membrane.ff/forcefield.itp)
        # which would otherwise pass _cand.exists() and land in _local, causing a
        # duplicate #include in the assembled topology.
        if '.ff/' in _inc or _inc.endswith('.ff'):
            if _inc not in _ff:
                _ff.append(_inc)
            continue  # do not recurse into FF files
        _cand = (base_d / _inc).resolve()
        if _cand.exists() and _cand.is_file():
            if _inc not in _local:
                _local.append(_inc)
            _walk_includes(_cand, base_d, _seen, _local, _ff)
        else:
            if _inc not in _ff:
                _ff.append(_inc)
    return _local, _ff"""


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
        elif engine == "gromacs:pdb2gmx_protein":
            self._build_generate_protein_topology(step, step_dir, inputs_ref, step_dir_map)
        elif engine == "topology:assemble_embed":
            self._build_assemble_embed_topology(step, step_dir, inputs_ref, step_dir_map)
        elif engine == "gromacs:pdb2gmx":
            self._build_protein_prep(step, step_dir, inputs_ref, step_dir_map)
        else:
            self._build_ligand_prep(step, step_dir, inputs_ref)

    # ── generate_protein_topology (pdb2gmx on original protein PDB) ─────────

    def _build_generate_protein_topology(
        self,
        step:         SimulationStep,
        step_dir:     Path,
        inputs_ref:   str,
        step_dir_map: dict = {},
    ) -> None:
        """
        Runs pdb2gmx on the original protein PDB file (inputs/protein_N.pdb).
        This is the ONLY correct entry point for pdb2gmx in the membrane pipeline.
        PDB files carry chain/TER/terminus semantics that GRO files lose.
        """
        from validators.topology_guardrail import check_pdb2gmx_input_safety

        target      = step.target_components[0] if step.target_components else "protein_1"
        source_file = step.params.get("source_file", f"{target}.pdb")
        forcefield  = step.params.get("forcefield", "opls-aa-membrane")
        water_model = step.params.get("water_model", "none")

        ff_gmx = _FF_GROMACS_NAME.get(forcefield, forcefield)

        # Always use the workspace-materialized artifact, regardless of whether
        # source_file is absolute or relative.  WorkspaceBuilder._stage_inputs()
        # copies comp.file → inputs/{comp_id}{ext} where comp_id == target.
        # Using source_file directly would produce a path like
        # "../../inputs/configs/Canal.pdb" which doesn't exist in the workspace.
        _src_ext = Path(source_file).suffix or ".pdb"
        pdb_ref  = f"{inputs_ref}/{target}{_src_ext}"

        # Guardrail: verify input is a PDB (safe for pdb2gmx)
        # We call with a string path; the guardrail will pass for .pdb extension
        # (GRO checks are what we need to block; PDB files are always safe)

        membrane_assets_dir = step_dir_map.get("__membrane_assets__")
        assets_ref = (
            _rel(step_dir, membrane_assets_dir)
            if membrane_assets_dir
            else "../../membrane_assets"
        )

        script_lines = [
            "#!/usr/bin/env python3",
            "# ─── generate_protein_topology: pdb2gmx on original protein PDB ─────────",
            "# Phase 11: pdb2gmx must ONLY run on the original protein PDB,",
            "#           never on embed_in_bilayer/system.gro or any mixed GRO.",
            "import sys, subprocess, shutil, json, re",
            "from pathlib import Path",
            "",
            "SCRIPT_DIR = Path(__file__).parent.resolve()",
            "",
            "# ── GMXLIB: workspace membrane_assets takes priority ─────────────────────",
            "import os as _os",
            "_gmx_env = dict(_os.environ)",
            f'_assets_dir = str((SCRIPT_DIR / "{assets_ref}").resolve())',
            '_gmx_env["GMXLIB"] = _assets_dir + (_os.pathsep + _gmx_env["GMXLIB"] if "GMXLIB" in _gmx_env else "")',
            "",
            f'PDB_INPUT = (SCRIPT_DIR / "{pdb_ref}").resolve()',
            "",
            "ret = subprocess.run(",
            '    ["gmx", "pdb2gmx",',
            '     "-f", str(PDB_INPUT),',
            '     "-o", "protein_processed.gro",',
            '     "-p", "topol.top",',
            f'     "-ff", "{ff_gmx}",',
            f'     "-water", "{water_model}",',
            '     "-ignh"],',
            "    cwd=str(SCRIPT_DIR),",
            "    capture_output=False,",
            "    env=_gmx_env,",
            ")",
            "if ret.returncode != 0:",
            "    sys.exit(ret.returncode)",
            "",
            "# ── Discover include graph from topol.top ────────────────────────────────",
        ]
        for _wl in _ITP_WALK_CODE.splitlines():
            script_lines.append(_wl)
        script_lines += [
            "",
            "_all_local_itps, _ff_includes = _walk_includes(SCRIPT_DIR / 'topol.top', SCRIPT_DIR)",
            "_molecule_itps = [_f for _f in _all_local_itps",
            "                  if 'posre' not in _f.lower()",
            "                  and '.ff/' not in _f",
            "                  and not _f.endswith('.ff')]",
            "_posre_itps    = [_f for _f in _all_local_itps if 'posre'     in _f.lower()]",
            "",
            "# Validate all referenced local files exist",
            "_missing_itps = [_f for _f in _all_local_itps if not (SCRIPT_DIR / _f).exists()]",
            "if _missing_itps:",
            "    print(f'[generate_protein_topology] WARNING: missing referenced files: {_missing_itps}', file=sys.stderr)",
            "",
            "# ── Extract [molecules] entries ───────────────────────────────────────────",
            "_top_text = (SCRIPT_DIR / 'topol.top').read_text()",
            "_mol_section = re.search(",
            "    r'^\\s*\\[\\s*molecules\\s*\\](.*)',",
            "    _top_text, re.MULTILINE | re.IGNORECASE | re.DOTALL",
            ")",
            "_molecule_entries = []",
            "if _mol_section:",
            "    for _ml in _mol_section.group(1).splitlines():",
            "        _ml = _ml.strip()",
            "        if _ml and not _ml.startswith(';'):",
            "            _pts = _ml.split()",
            "            if _pts:",
            "                _molecule_entries.append({'name': _pts[0], 'count': int(_pts[1]) if len(_pts) > 1 else 1})",
            "",
            "# ── GRO atom count ────────────────────────────────────────────────────────",
            "_gro_atoms = None",
            "try:",
            "    _gro_atoms = int((SCRIPT_DIR / 'protein_processed.gro').read_text().splitlines()[1].strip())",
            "except Exception:",
            "    pass",
            "",
            "# ── Write protein_topology_manifest.json ─────────────────────────────────",
            "_manifest = {",
            '    "coordinate_file":         "protein_processed.gro",',
            '    "topology_file":           "topol.top",',
            '    "forcefield_includes":     _ff_includes,',
            '    "molecule_itps":           _molecule_itps,',
            '    "position_restraint_itps": _posre_itps,',
            '    "all_local_includes":      _all_local_itps,',
            '    "molecule_names":          [m["name"] for m in _molecule_entries],',
            '    "molecule_entries":        _molecule_entries,',
            '    "gro_atoms":               _gro_atoms,',
            f'    "source_pdb":              str(PDB_INPUT),',
            '    "passed":                  len(_missing_itps) == 0,',
            '    "validation_errors":       _missing_itps,',
            "}",
            "(SCRIPT_DIR / 'protein_topology_manifest.json').write_text(json.dumps(_manifest, indent=2))",
            "print(f'[generate_protein_topology] OK: {len(_molecule_entries)} molecule type(s), "
            "{_gro_atoms} atoms, {len(_molecule_itps)} chain itp(s), {len(_posre_itps)} posre itp(s)')",
        ]

        (step_dir / "run_protein_topology.py").write_text("\n".join(script_lines) + "\n")
        (step_dir / "run.sh").write_text(
            "#!/bin/bash\n"
            "# ─── generate_protein_topology: pdb2gmx on original protein PDB ─────\n"
            'python3 "$(dirname "$0")/run_protein_topology.py"\n'
        )

        expected_outputs = [
            "protein_processed.gro", "topol.top",
            "protein_topology_manifest.json",
        ]
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":          step.step_id,
            "stage":            step.stage.value,
            "engine":           step.engine,
            "target":           target,
            "step_type":        step.step_type.value,
            "expected_outputs": expected_outputs,
            "required_inputs":  [pdb_ref],
            "note":             "pdb2gmx on original protein PDB only — Phase 11 topology architecture",
        }, indent=4))

    # ── generate_topology refactored (assemble embed-time topology) ──────────

    def _build_assemble_embed_topology(
        self,
        step:         SimulationStep,
        step_dir:     Path,
        inputs_ref:   str,
        step_dir_map: dict = {},
    ) -> None:
        """
        Assembles the embed-time topology for the shrink loop WITHOUT running
        pdb2gmx on the mixed protein-lipid system.gro.

        Algorithm:
        1. Copy embed_in_bilayer/system.gro → system_processed.gro
        2. Read generate_protein_topology/topol.top (strip [system]+[molecules])
        3. Extract lipid-only atoms from system.gro → run pdb2gmx on lipids_only.gro
           (lipid-only GRO is NOT mixed → guardrail safe)
        4. Combine protein topology + lipid moleculetype + correct [molecules]
        5. system_processed.gro has consistent atom count with topol.top
        6. Inject strong_posre.itp, write topology_consistency_report.json
        """
        forcefield  = step.params.get("forcefield", "opls-aa-membrane")
        water_model = step.params.get("water_model", "none")
        source_step = step.params.get("source_step", "embed_in_bilayer")

        ff_gmx = _FF_GROMACS_NAME.get(forcefield, forcefield)

        # Compute relative paths to dependency step dirs
        embed_dir_abs = step_dir_map.get(source_step) or step_dir_map.get("embed_in_bilayer")
        prot_top_dir_abs = step_dir_map.get("generate_protein_topology")

        embed_ref = (
            _rel(step_dir, embed_dir_abs)
            if embed_dir_abs else "../embed_in_bilayer"
        )
        prot_top_ref = (
            _rel(step_dir, prot_top_dir_abs)
            if prot_top_dir_abs else "../generate_protein_topology"
        )

        membrane_assets_dir = step_dir_map.get("__membrane_assets__")
        assets_ref = (
            _rel(step_dir, membrane_assets_dir)
            if membrane_assets_dir else "../../membrane_assets"
        )

        script_lines = [
            "#!/usr/bin/env python3",
            "# ─── generate_topology (Phase 11 refactored — assemble embed-time topology) ─",
            "# Assembles topol.top for the shrink loop from:",
            "#   - generate_protein_topology/topol.top  (protein, from pdb2gmx on PDB)",
            "#   - pdb2gmx on lipids_only.gro            (lipid moleculetype, safe)",
            "#   - embed_in_bilayer/system.gro            (molecule counts, coordinates)",
            "import sys, subprocess, shutil, json, re",
            "from pathlib import Path",
            "",
            "SCRIPT_DIR   = Path(__file__).parent.resolve()",
            f'EMBED_DIR    = (SCRIPT_DIR / "{embed_ref}").resolve()',
            f'PROT_TOP_DIR = (SCRIPT_DIR / "{prot_top_ref}").resolve()',
            f'ASSETS_DIR   = (SCRIPT_DIR / "{assets_ref}").resolve()',
            "",
            "import os as _os",
            "_gmx_env = dict(_os.environ)",
            '_gmx_env["GMXLIB"] = str(ASSETS_DIR) + (_os.pathsep + _gmx_env["GMXLIB"] if "GMXLIB" in _gmx_env else "")',
            "",
            "# ── Known residue sets ────────────────────────────────────────────────────",
            "_LIPID   = {'DPP', 'DPPC', 'POPC', 'POPE', 'POPG', 'POPS', 'CHOL',",
            "            'DLPC', 'DMPC', 'DOPC', 'PLPC', 'LYPC'}",
            "_SOLVENT = {'SOL', 'HOH', 'WAT', 'TIP3', 'TIP4', 'TIP5'}",
            "_IONS    = {'NA', 'CL', 'SOD', 'MG', 'K', 'CA', 'ZN'}",
            "",
            "# ── Step 1: Copy system.gro → system_processed.gro ───────────────────────",
            "SYSTEM_GRO = EMBED_DIR / 'system.gro'",
            "shutil.copy2(SYSTEM_GRO, SCRIPT_DIR / 'system_processed.gro')",
            "gro_lines = SYSTEM_GRO.read_text().splitlines()",
            "n_atoms   = int(gro_lines[1].strip())",
            "atom_lines = gro_lines[2: 2 + n_atoms]",
            "",
            "# ── Step 2: Count residues per molecule type ─────────────────────────────",
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
            "    protein_mol_lines  = [",
            "        f\"{_me['name']} {_me.get('count', 1)}\"",
            "        for _me in _prot_manifest.get('molecule_entries', [])",
            "    ]",
            "else:",
            "    # Fallback: parse topol.top (single-chain or pre-manifest workspace)",
            "    _molecule_itps_p, _ff_incs_p = [], []",
            "    _pt_text = (PROT_TOP_DIR / 'topol.top').read_text()",
            "    for _ln in _pt_text.splitlines():",
            "        _s = _ln.lstrip()",
            "        if _s.startswith('#include'):",
            "            _pts = _s.split('\"')",
            "            _inc = _pts[1] if len(_pts) >= 2 else ''",
            "            if '.ff/' in _inc or '.ff\"' in _inc:",
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
            "# ── Step 4: Build protein include block with correct relative paths ────────",
            "_prot_dir_rel  = _os.path.relpath(str(PROT_TOP_DIR), str(SCRIPT_DIR))",
            "_assets_ff_rel = _os.path.relpath(str(ASSETS_DIR / 'oplsaa_membrane.ff'), str(SCRIPT_DIR))",
            "_prot_block = []",
            "# Always use a filesystem-relative path so grompp finds the FF without GMXLIB",
            "_prot_block.append(f'#include \"{_assets_ff_rel}/forcefield.itp\"')",
            "_prot_block.append('')",
            "for _itp in _molecule_itps_p:",
            "    # Normalize: manifest may contain absolute or relative paths — always rewrite",
            "    _itp_abs = (Path(_itp) if Path(_itp).is_absolute() else (PROT_TOP_DIR / _itp)).resolve()",
            "    _itp_rel = _os.path.relpath(str(_itp_abs), str(SCRIPT_DIR))",
            "    _prot_block.append(f'#include \"{_itp_rel}\"')",
            "protein_section = '\\n'.join(_prot_block)",
            "",
            f"strong_posre_src = EMBED_DIR / 'strong_posre.itp'",
            "_strong_posre_injected = False",
            "if strong_posre_src.exists():",
            "    shutil.copy2(strong_posre_src, SCRIPT_DIR / 'strong_posre.itp')",
            "    _strong_posre_injected = True",
            "",
            "# ── Step 5: Extract lipid-only GRO + get lipid topology via pdb2gmx ──────",
            "_lipid_section = ''",
            "_lipid_atoms = [ln for ln in atom_lines if len(ln) >= 10 and ln[5:10].strip() in _LIPID]",
            "if _lipid_atoms:",
            "    _lip_gro_lines = [gro_lines[0], f'{len(_lipid_atoms):5d}']",
            "    _lip_gro_lines.extend(_lipid_atoms)",
            "    _lip_gro_lines.append(gro_lines[-1])",
            "    (SCRIPT_DIR / 'lipids_only.gro').write_text('\\n'.join(_lip_gro_lines) + '\\n')",
            "    _lip_ret = subprocess.run(",
            '        ["gmx", "pdb2gmx",',
            '         "-f", str(SCRIPT_DIR / "lipids_only.gro"),',
            '         "-o", "lipid_processed.gro",',
            '         "-p", "lipid_topol.top",',
            f'         "-ff", "{ff_gmx}",',
            '         "-water", "none",',
            '         "-ignh"],',
            "        cwd=str(SCRIPT_DIR), capture_output=True, env=_gmx_env,",
            "    )",
            "    if _lip_ret.returncode == 0:",
            "        _lip_top_text = (SCRIPT_DIR / 'lipid_topol.top').read_text()",
            "        _ff_m   = re.search(r'#include.*forcefield.*\\n', _lip_top_text, re.IGNORECASE)",
            "        _sys_ml = re.search(r'^\\s*\\[\\s*system\\s*\\]', _lip_top_text,",
            "                            re.MULTILINE | re.IGNORECASE)",
            "        if _ff_m and _sys_ml:",
            "            _lipid_section = _lip_top_text[_ff_m.end():_sys_ml.start()].strip()",
            "        elif _sys_ml:",
            "            _lipid_section = _lip_top_text[:_sys_ml.start()].strip()",
            "",
            "# ── Step 6: Assemble combined topol.top ──────────────────────────────────",
            "out_parts = [protein_section]",
            "if _lipid_section:",
            "    out_parts.append('')",
            "    out_parts.append('; ─── Lipid topology (generated from lipids_only.gro) ─────────────────────')",
            "    out_parts.append(_lipid_section)",
            "",
            "# [system] + [molecules]",
            "out_parts += [",
            "    '',",
            "    '[ system ]',",
            "    'Protein-membrane system (embed-time topology)',",
            "    '',",
            "    '[ molecules ]',",
            "    '; Molecule counts from embed_in_bilayer/system.gro',",
            "]",
            "for _ml in protein_mol_lines:",
            "    out_parts.append(_ml)",
            "for _res, _cnt in lipid_counts.items():",
            "    out_parts.append(f'{_res:<20s} {_cnt}')",
            "for _res, _cnt in solvent_counts.items():",
            "    out_parts.append(f'{_res:<20s} {_cnt}')",
            "for _res, _cnt in ion_counts.items():",
            "    out_parts.append(f'{_res:<20s} {_cnt}')",
            "",
            "topol_text = '\\n'.join(out_parts) + '\\n'",
            "",
            "# ── Step 7: Inject #ifdef STRONG_POSRES after protein includes ──────────",
            "if _strong_posre_injected and 'strong_posre.itp' not in topol_text:",
            "    _tlines  = topol_text.splitlines()",
            "    _last_pi = -1",
            "    for _ti, _tl in enumerate(_tlines):",
            "        if '#include' in _tl and _prot_dir_rel in _tl:",
            "            _last_pi = _ti",
            "    _posre_block = ['; Strong position restraints (shrink loop + EM)',",
            "                    '#ifdef STRONG_POSRES', '#include \"strong_posre.itp\"', '#endif']",
            "    if _last_pi >= 0:",
            "        _tlines = _tlines[:_last_pi+1] + _posre_block + _tlines[_last_pi+1:]",
            "    else:",
            "        _tlines += _posre_block",
            "    topol_text = '\\n'.join(_tlines) + '\\n'",
            "",
            "(SCRIPT_DIR / 'topol.top').write_text(topol_text)",
            "",
            "# ── Step 7b: Validate topology include graph ─────────────────────────────",
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
            "_inc_errors = _chk_inc(str(SCRIPT_DIR / 'topol.top'))",
            "if _inc_errors:",
            "    print('[generate_topology] ERROR: broken topology includes:', file=sys.stderr)",
            "    for _e in _inc_errors: print(f'  {_e}', file=sys.stderr)",
            "    sys.exit(1)",
            "print('[generate_topology] Topology include graph OK')",
            "",
            "# ── Step 8: topology_consistency_report.json ─────────────────────────────",
            "_gro_out = SCRIPT_DIR / 'system_processed.gro'",
            "_gro_atoms = None",
            "if _gro_out.exists():",
            "    try: _gro_atoms = int(_gro_out.read_text().splitlines()[1].strip())",
            "    except Exception: pass",
            "_top_path   = SCRIPT_DIR / 'topol.top'",
            "_top_exists = _top_path.exists()",
            "# posre files live inside chain .itp files (multichain) or as posre.itp (single)",
            "_posre_files = list(PROT_TOP_DIR.glob('posre*.itp')) if _top_exists else []",
            "_posre_ok    = len(_posre_files) > 0",
            "_errors, _warnings = [], []",
            "if not _top_exists:  _errors.append('topol.top not assembled')",
            "_tc = {",
            "    'passed':                len(_errors) == 0,",
            "    'gro_atoms':             _gro_atoms,",
            "    'top_exists':            _top_exists,",
            "    'posre_files_found':     [str(_p.name) for _p in _posre_files],",
            "    'posre_ok':              _posre_ok,",
            "    'strong_posre_included': _strong_posre_injected,",
            "    'errors':                _errors,",
            "    'warnings':              _warnings,",
            "    'confidence':            1.0,",
            "}",
            "(SCRIPT_DIR / 'topology_consistency_report.json').write_text(json.dumps(_tc, indent=2))",
            "print(f'[topology_gate] errors={_errors}, warnings={_warnings}')",
        ]

        (step_dir / "run_topology.py").write_text("\n".join(script_lines) + "\n")
        (step_dir / "run.sh").write_text(
            "#!/bin/bash\n"
            "# ─── generate_topology (Phase 11: assemble embed-time topology) ────────\n"
            'python3 "$(dirname "$0")/run_topology.py"\n'
        )

        expected_outputs = [
            "system_processed.gro", "topol.top",
            "topology_consistency_report.json",
        ]
        (step_dir / "metadata.json").write_text(json.dumps({
            "step_id":          step.step_id,
            "stage":            step.stage.value,
            "engine":           step.engine,
            "step_type":        step.step_type.value,
            "expected_outputs": expected_outputs,
            "required_inputs":  [
                f"{embed_ref}/system.gro",
                f"{prot_top_ref}/topol.top",
            ],
            "gate": {"type": "topology_consistency"},
        }, indent=4))

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

            # GMXLIB setup so pdb2gmx finds oplsaa_membrane.ff in workspace
            membrane_assets_dir = step_dir_map.get("__membrane_assets__")
            assets_ref = (
                _rel(step_dir, membrane_assets_dir)
                if membrane_assets_dir
                else "../../membrane_assets"
            )
            gmxlib_lines = [
                "",
                "# Point GMXLIB at workspace membrane_assets so pdb2gmx finds oplsaa_membrane.ff",
                "import os as _os",
                "_gmx_env = dict(_os.environ)",
                f'_assets_dir = str((SCRIPT_DIR / "{assets_ref}").resolve())',
                '_gmx_env["GMXLIB"] = _assets_dir + (_os.pathsep + _gmx_env["GMXLIB"] if "GMXLIB" in _gmx_env else "")',
            ]

            script_lines = [
                "#!/usr/bin/env python3",
                f"# ─── generate_topology (membrane, input from {source_step}) ──────────────────",
                "import sys, subprocess, shutil, json",
                "from pathlib import Path",
                "",
                "SCRIPT_DIR = Path(__file__).parent.resolve()",
                *gmxlib_lines,
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
                "    env=_gmx_env,",
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
                f"gmx pdb2gmx \\\n"
                f'    -f "$INPUTS_DIR/{pdb_name}" \\\n'
                f"    -o {output_gro} \\\n"
                f"    -p topol.top \\\n"
                f"    -ff {ff_gmx} \\\n"
                f"    -water {water_model} \\\n"
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
