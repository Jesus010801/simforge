"""
builders/step_builders/embedding_builder.py

Builds the membrane embedding meta-step: a self-contained shrink loop that:
  1. Inflates the bilayer lipids radially (creates space around protein)
  2. Iteratively deflates (scale 0.95) + minimizes until APL converges
  3. Writes per-iteration telemetry to shrink_telemetry.json
  4. Produces a single deterministic output: converged.gro

The entire loop lives in shrink_loop.sh — the executor treats it as one opaque
step.  No adaptive DAG.  Convergence = APL <= apl_target + tolerance.

Physical constants come from core/membrane_knowledge.py.
APL is read directly from area_2.dat in Python — no Fortran AperR needed.
"""

from __future__ import annotations

import json
from pathlib import Path

from core.execution_models import SimulationStep
from builders.step_builders._utils import rel as _rel


class EmbeddingBuilder:

    def build(
        self,
        step:         SimulationStep,
        step_dir:     Path,
        step_dir_map: dict = {},
    ) -> None:

        p = step.params

        # ── Physical parameters (from membrane_knowledge via decision_engine) ──
        lipid               = p.get("lipid",               "DPPC")
        lipid_residue_name  = p.get("lipid_residue_name",  "DPP")
        forcefield          = p.get("forcefield",           "opls-aa")
        temperature_K       = float(p.get("temperature_K", 298.0))
        apl_target_ang2     = float(p.get("apl_target_ang2", 62.0))
        apl_tolerance_ang2  = float(p.get("apl_tolerance_ang2", 2.0))
        inflate_factor      = float(p.get("inflate_factor",  4.0))
        deflate_factor      = float(p.get("deflate_factor",  0.95))
        max_iterations      = int(p.get("max_iterations",    200))
        gridsize            = int(p.get("gridsize",          5))
        cutoff              = float(p.get("cutoff",          0.0))

        # ── Tool paths (set by pipeline or user config) ────────────────────────
        inflategro_script = p.get("inflategro_script", "inflategro-Jorge.pl")
        topol_top         = p.get("topol_top",         "topol.top")
        input_gro         = p.get("input_gro",         "system.gro")

        # Override with DAG-resolved paths if available
        embed_dir = next(
            (step_dir_map[d] for d in step.depends_on
             if "embed" in d and d in step_dir_map),
            None,
        )
        topology_dir = next(
            (step_dir_map[d] for d in step.depends_on
             if "topology" in d or "topol" in d and d in step_dir_map),
            None,
        )
        if embed_dir:
            input_gro = str(Path(_rel(step_dir, embed_dir)) / "system.gro")
        if topology_dir:
            topol_top = str(Path(_rel(step_dir, topology_dir)) / "topol.top")

        # ── Shrink-loop MDP (position-restrained minimization) ─────────────────
        self._write_minim_mdp(step_dir, temperature_K)

        # ── Main script ────────────────────────────────────────────────────────
        self._write_shrink_loop(
            step_dir         = step_dir,
            inflategro_script= inflategro_script,
            input_gro        = input_gro,
            topol_top        = topol_top,
            lipid_name       = lipid_residue_name,
            inflate_factor   = inflate_factor,
            deflate_factor   = deflate_factor,
            apl_target_ang2  = apl_target_ang2,
            apl_tolerance_ang2 = apl_tolerance_ang2,
            max_iterations   = max_iterations,
            gridsize         = gridsize,
            cutoff           = cutoff,
        )

        # ── Metadata ───────────────────────────────────────────────────────────
        metadata = {
            "step_id":        step.step_id,
            "stage":          step.stage.value,
            "engine":         step.engine,
            "step_type":      step.step_type.value,
            "blocking":       step.blocking,
            "generated_by":   "EmbeddingBuilder",
            "gate":           {"type": "apl_report"},
            "expected_outputs": ["converged.gro", "shrink_telemetry.json"],
            "params": {
                "lipid":              lipid,
                "lipid_residue_name": lipid_residue_name,
                "forcefield":         forcefield,
                "temperature_K":      temperature_K,
                "apl_target_ang2":    apl_target_ang2,
                "apl_tolerance_ang2": apl_tolerance_ang2,
                "inflate_factor":     inflate_factor,
                "deflate_factor":     deflate_factor,
                "max_iterations":     max_iterations,
            },
        }
        (step_dir / "metadata.json").write_text(json.dumps(metadata, indent=4))

    # ── MDP for shrink-loop minimization ──────────────────────────────────────

    def _write_minim_mdp(self, step_dir: Path, temperature_K: float) -> None:
        # Tight emtol + both POSRES and STRONG_POSRES active during embedding
        mdp = f"""; Shrink-loop minimization — position restraints active
define          = -DPOSRES -DSTRONG_POSRES
integrator      = steep
emtol           = 1000.0
emstep          = 0.01
nsteps          = 50000

cutoff-scheme   = Verlet
nstlist         = 1
coulombtype     = PME
rcoulomb        = 1.2
rvdw            = 1.2
pbc             = xyz
DispCorr        = EnerPres
"""
        (step_dir / "minim_shrink.mdp").write_text(mdp.strip())

    # ── Main shrink-loop script ───────────────────────────────────────────────

    def _write_shrink_loop(
        self,
        step_dir:           Path,
        inflategro_script:  str,
        input_gro:          str,
        topol_top:          str,
        lipid_name:         str,
        inflate_factor:     float,
        deflate_factor:     float,
        apl_target_ang2:    float,
        apl_tolerance_ang2: float,
        max_iterations:     int,
        gridsize:           int,
        cutoff:             float,
    ) -> None:

        apl_cutoff = apl_target_ang2 + apl_tolerance_ang2

        script = f"""#!/bin/bash
# ─── Membrane embedding — InflateGRO shrink loop ──────────────────────────────
# Generated by EmbeddingBuilder. Do not edit manually.
#
# Physical target: APL ≤ {apl_cutoff:.1f} Å² ({lipid_name}, {apl_target_ang2:.0f} + {apl_tolerance_ang2:.0f} Å² tolerance)
# Convergence: iterative deflation (factor {deflate_factor}) until APL converges.
# Telemetry: shrink_telemetry.json — one entry per iteration.
#
set -e
cd "$(dirname "$0")"

INFLATEGRO="{inflategro_script}"
INPUT_GRO="{input_gro}"
TOPOL="{topol_top}"
LIPID="{lipid_name}"
INFLATE={inflate_factor}
DEFLATE={deflate_factor}
APL_TARGET={apl_cutoff}
MAX_ITER={max_iterations}
GRIDSIZE={gridsize}
CUTOFF={cutoff}

# ── Python helper: read APL in Å² from area_2.dat ────────────────────────────
read_apl() {{
    python3 -c "
import sys
try:
    val = float(open('area_2.dat').read().strip())
    print(int(val * 100))
except Exception as e:
    print(0)
    sys.exit(1)
"
}}

# ── Python helper: append one iteration to shrink_telemetry.json ─────────────
log_iter() {{
    local iter=$1 apl=$2 converged=$3
    python3 -c "
import json, pathlib, datetime
p = pathlib.Path('shrink_telemetry.json')
data = json.loads(p.read_text()) if p.exists() else {{'iterations': [], 'converged': False, 'started_at': '$( date -Iseconds )'}}
data['iterations'].append({{'iter': $iter, 'apl_ang2': $apl, 'converged': $converged == 'true'}})
data['last_apl_ang2'] = $apl
p.write_text(json.dumps(data, indent=2))
"
}}

finalize_telemetry() {{
    local apl=$1 converged=$2 n_iter=$3
    python3 -c "
import json, pathlib, datetime
p = pathlib.Path('shrink_telemetry.json')
data = json.loads(p.read_text()) if p.exists() else {{'iterations': []}}
data['converged'] = $converged == 'true'
data['final_apl_ang2'] = $apl
data['n_iterations'] = $n_iter
data['finished_at'] = datetime.datetime.now().isoformat(timespec='seconds')
p.write_text(json.dumps(data, indent=2))
"
}}

# ── Validate tools ────────────────────────────────────────────────────────────
if [ ! -f "$INFLATEGRO" ]; then
    echo "[embedding] ERROR: inflategro script not found: $INFLATEGRO"
    echo "[embedding] Set inflategro_script in your YAML config or copy the script here."
    exit 1
fi
if ! command -v perl &>/dev/null; then
    echo "[embedding] ERROR: perl not found on PATH"
    exit 1
fi
if ! command -v gmx &>/dev/null; then
    echo "[embedding] ERROR: gmx not found on PATH"
    exit 1
fi

echo "[embedding] Starting membrane embedding shrink loop"
echo "[embedding] Lipid: $LIPID  |  APL target: ≤ ${{APL_TARGET}} Å²  |  Max iterations: $MAX_ITER"

# ── Step 1: initial inflation ─────────────────────────────────────────────────
echo "[embedding] Step 1: inflating bilayer (factor $INFLATE)"
cp "$INPUT_GRO" work.gro

perl "$INFLATEGRO" work.gro $INFLATE $LIPID $CUTOFF inflated.gro $GRIDSIZE area_2.dat
APL=$(read_apl)
echo "[embedding] Post-inflation APL = ${{APL}} Å²"
log_iter 0 $APL false
cp inflated.gro work.gro

# ── Step 2: first minimization (inflated system) ──────────────────────────────
echo "[embedding] Step 2: initial minimization (inflated system)"
gmx grompp -f minim_shrink.mdp -c work.gro -r work.gro -p "$TOPOL" -o work.tpr -maxwarn 2 -quiet
gmx mdrun -s work.tpr -deffnm work -nb gpu -quiet
cp work.gro work_prev.gro

# ── Step 3: shrink loop ───────────────────────────────────────────────────────
echo "[embedding] Step 3: shrink loop (deflate $DEFLATE per iteration)"
ITER=1
APL=$( python3 -c "print(int(float(open('area_2.dat').read().strip()) * 100))" )

while [ "$APL" -gt "$( python3 -c "print(int($APL_TARGET))" )" ] && [ "$ITER" -le "$MAX_ITER" ]; do

    echo "[embedding] Iter $ITER: APL = ${{APL}} Å² > ${{APL_TARGET}} Å² — deflating"

    perl "$INFLATEGRO" work.gro $DEFLATE $LIPID 0 work.gro $GRIDSIZE area_2.dat
    APL=$(read_apl)
    log_iter $ITER $APL false

    gmx grompp -f minim_shrink.mdp -c work.gro -r work.gro -p "$TOPOL" -o work.tpr -maxwarn 2 -quiet
    gmx mdrun -s work.tpr -deffnm work -nb gpu -quiet

    # Clean GROMACS backup files
    rm -f \\#*

    ITER=$(( ITER + 1 ))
done

# ── Step 4: check convergence ─────────────────────────────────────────────────
APL=$(read_apl)
log_iter $ITER $APL true

if [ "$APL" -le "$( python3 -c "print(int($APL_TARGET))" )" ]; then
    echo "[embedding] CONVERGED at iter $ITER: APL = ${{APL}} Å² ≤ ${{APL_TARGET}} Å²"
    CONVERGED=true
else
    echo "[embedding] WARNING: did not converge after $MAX_ITER iterations (APL = ${{APL}} Å²)"
    CONVERGED=false
fi

# ── Step 5: final minimization on converged system ────────────────────────────
echo "[embedding] Step 5: final minimization"
gmx grompp -f minim_shrink.mdp -c work.gro -r work.gro -p "$TOPOL" -o final.tpr -maxwarn 2 -quiet
gmx mdrun -s final.tpr -deffnm final_shrink -nb gpu -quiet
rm -f \\#*

# ── Step 6: write canonical output ───────────────────────────────────────────
cp final_shrink.gro converged.gro
finalize_telemetry $APL $CONVERGED $ITER

echo "[embedding] Done. Output: converged.gro  |  APL = ${{APL}} Å²  |  Iterations: $ITER"

if [ "$CONVERGED" = false ]; then
    echo "[embedding] HINT: APL did not reach target. Try increasing inflate_factor (current: $INFLATE)."
    exit 1
fi
"""
        script_path = step_dir / "run.sh"
        script_path.write_text(script.lstrip())
        script_path.chmod(0o755)
