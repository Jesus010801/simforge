# builders/step_builders/analysis_builder.py

from __future__ import annotations

from pathlib import Path
import json

from core.execution_models import SimulationStep
from builders.step_builders._utils import rel as _rel


# ═══════════════════════════════════════════════════════════════════════════════
# Analysis Builder
# ═══════════════════════════════════════════════════════════════════════════════

# GROMACS built-in group numbers (standard pdb2gmx output)
_GROMACS_DEFAULT_GROUPS = """\
# GROMACS default groups (pdb2gmx standard output):
#   0  System          — all atoms
#   1  Protein         — protein atoms
#   2  Protein-H       — protein heavy atoms
#   3  C-alpha         — Cα atoms only
#   4  Backbone        — N, Cα, C
#   5  MainChain       — backbone + Cβ
#   6  MainChain+Cb    — MainChain + Cβ
#   7  MainChain+H     — MainChain + backbone H
#   8  SideChain       — side chain atoms
#   9  SideChain-H     — side chain heavy atoms
#  10  Prot-Masses     — protein with masses
#  11  non-Protein     — water + ions + ligands
#  12  Other           — non-protein non-water
#  13  SOL             — water molecules
#  14  non-Water       — everything except water
#  15+ NA, CL ...      — individual ion species (if present)
"""


class AnalysisBuilder:
    """
    Genera scripts GROMACS reales para análisis post-simulación.

    Dispatch por analysis_type:
        rmsd             → gmx rms  (Backbone ref, selección configurable)
        rmsf             → gmx rmsf (por residuo, Backbone/C-alpha)
        hydrogen_bonds   → gmx hbond (Protein donors/acceptors)
        distance         → gmx distance (selección vía make_ndx)
        energy           → gmx energy (ETot, Temp, Pressure, Pot)
        radius_of_gyration → gmx gyrate (Protein)

    Grupos: GROMACS built-ins exclusivamente — no se necesita make_ndx
    para rmsd/rmsf/hbond/energy/rgyr. Para distance se emite make_ndx.sh
    con comentarios de guía.
    """

    def build(
        self,
        step:         SimulationStep,
        step_dir:     Path,
        step_dir_map: dict = {},
    ) -> None:

        (step_dir / "tables").mkdir(exist_ok=True)
        (step_dir / "plots").mkdir(exist_ok=True)

        analysis_type = step.params.get(
            "analysis_type",
            step.step_id.removeprefix("analysis_"),
        )

        # Resolver path a producción (md.xtc, md.tpr)
        prod_dir = (
            step_dir_map.get("production_md")
            or step_dir_map.get("production")
        )
        prod_ref = _rel(step_dir, prod_dir) if prod_dir else "../production_md"

        dispatch = {
            "rmsd":               self._build_rmsd,
            "rmsf":               self._build_rmsf,
            "hydrogen_bonds":     self._build_hbond,
            "hbond":              self._build_hbond,
            "distance":           self._build_distance,
            "distance_analysis":  self._build_distance,
            "energy":             self._build_energy,
            "radius_of_gyration": self._build_rgyr,
            "rgyr":               self._build_rgyr,
        }

        builder_fn = dispatch.get(analysis_type, self._build_generic)
        builder_fn(step, step_dir, prod_ref, analysis_type)

    # ── RMSD ─────────────────────────────────────────────────────────────────

    def _build_rmsd(
        self,
        step:          SimulationStep,
        step_dir:      Path,
        prod_ref:      str,
        analysis_type: str,
    ) -> None:

        selection = step.params.get("selection", "Backbone")

        script = f"""#!/bin/bash
# ─── RMSD analysis ───────────────────────────────────────────────────────────
# Reference: Backbone (group 4)
# Selection: {selection}
# Input groups use GROMACS built-in defaults — no make_ndx required.

PROD_DIR="{prod_ref}"

mkdir -p tables plots

# RMSD vs. initial structure (backbone reference + selection)
# Prompt: group 4 (Backbone) for least-squares fit, then selection for RMSD
echo "4 4" | gmx rms \\
    -s "$PROD_DIR/md.tpr" \\
    -f "$PROD_DIR/md.xtc" \\
    -o tables/rmsd_backbone.xvg \\
    -tu ns

# RMSD of selection vs. backbone reference
echo "4 1" | gmx rms \\
    -s "$PROD_DIR/md.tpr" \\
    -f "$PROD_DIR/md.xtc" \\
    -o tables/rmsd_protein.xvg \\
    -tu ns

echo "RMSD analysis complete → tables/rmsd_backbone.xvg  tables/rmsd_protein.xvg"
"""

        (step_dir / "run_analysis.sh").write_text(script.strip())
        self._write_metadata(step, step_dir, analysis_type, {
            "reference_group": "Backbone (4)",
            "selection": selection,
            "outputs": ["tables/rmsd_backbone.xvg", "tables/rmsd_protein.xvg"],
        })

    # ── RMSF ─────────────────────────────────────────────────────────────────

    def _build_rmsf(
        self,
        step:          SimulationStep,
        step_dir:      Path,
        prod_ref:      str,
        analysis_type: str,
    ) -> None:

        script = f"""#!/bin/bash
# ─── RMSF per-residue ────────────────────────────────────────────────────────
# Group 4 (Backbone) — per-residue fluctuation
# Input groups use GROMACS built-in defaults — no make_ndx required.

PROD_DIR="{prod_ref}"

mkdir -p tables plots

echo "4" | gmx rmsf \\
    -s "$PROD_DIR/md.tpr" \\
    -f "$PROD_DIR/md.xtc" \\
    -o tables/rmsf_backbone.xvg \\
    -res

echo "3" | gmx rmsf \\
    -s "$PROD_DIR/md.tpr" \\
    -f "$PROD_DIR/md.xtc" \\
    -o tables/rmsf_ca.xvg \\
    -res

echo "RMSF analysis complete → tables/rmsf_backbone.xvg  tables/rmsf_ca.xvg"
"""

        (step_dir / "run_analysis.sh").write_text(script.strip())
        self._write_metadata(step, step_dir, analysis_type, {
            "groups": ["Backbone (4)", "C-alpha (3)"],
            "outputs": ["tables/rmsf_backbone.xvg", "tables/rmsf_ca.xvg"],
        })

    # ── Hydrogen bonds ────────────────────────────────────────────────────────

    def _build_hbond(
        self,
        step:          SimulationStep,
        step_dir:      Path,
        prod_ref:      str,
        analysis_type: str,
    ) -> None:

        script = f"""#!/bin/bash
# ─── Hydrogen bond analysis ───────────────────────────────────────────────────
# Protein intra-molecular H-bonds (group 1 → group 1)
# Input groups use GROMACS built-in defaults — no make_ndx required.

PROD_DIR="{prod_ref}"

mkdir -p tables plots

# Intra-protein H-bonds
echo "1 1" | gmx hbond \\
    -s "$PROD_DIR/md.tpr" \\
    -f "$PROD_DIR/md.xtc" \\
    -num tables/hbnum_protein.xvg \\
    -dist tables/hbdist_protein.xvg

# Protein–solvent H-bonds
echo "1 13" | gmx hbond \\
    -s "$PROD_DIR/md.tpr" \\
    -f "$PROD_DIR/md.xtc" \\
    -num tables/hbnum_protein_sol.xvg

echo "H-bond analysis complete → tables/hbnum_protein.xvg  tables/hbnum_protein_sol.xvg"
"""

        (step_dir / "run_analysis.sh").write_text(script.strip())
        self._write_metadata(step, step_dir, analysis_type, {
            "analyses": ["Protein-Protein (1 1)", "Protein-SOL (1 13)"],
            "outputs": [
                "tables/hbnum_protein.xvg",
                "tables/hbdist_protein.xvg",
                "tables/hbnum_protein_sol.xvg",
            ],
        })

    # ── Distance ──────────────────────────────────────────────────────────────

    def _build_distance(
        self,
        step:          SimulationStep,
        step_dir:      Path,
        prod_ref:      str,
        analysis_type: str,
    ) -> None:

        selection = step.params.get("selection", "")
        sel_comment = f"# Selection from config: {selection}" if selection else \
                      "# No selection specified — define custom groups in make_ndx.sh"

        make_ndx_script = f"""#!/bin/bash
# ─── Custom index groups for distance analysis ────────────────────────────────
# Edit the 'name' and atom-selection expressions below to match your system.
# Run this BEFORE run_analysis.sh.
{_GROMACS_DEFAULT_GROUPS}
{sel_comment}

PROD_DIR="{prod_ref}"

# Interactive make_ndx — type selections at the prompt, then 'q' to save.
# Example selections:
#   "r 100 & a CA"    — Cα of residue 100
#   "r 200 & a CA"    — Cα of residue 200
# Then name them:
#   "name 16 ResA_CA"
#   "name 17 ResB_CA"
#   "q"
gmx make_ndx \\
    -f "$PROD_DIR/md.tpr" \\
    -o index.ndx

echo "Index file created → index.ndx"
echo "Edit run_analysis.sh with the correct group numbers before running."
"""

        script = f"""#!/bin/bash
# ─── Distance analysis ────────────────────────────────────────────────────────
# Requires index.ndx — run make_ndx.sh first.
{sel_comment}

PROD_DIR="{prod_ref}"

mkdir -p tables plots

# Replace GROUP_A and GROUP_B with the group names from index.ndx
echo "GROUP_A GROUP_B" | gmx distance \\
    -s "$PROD_DIR/md.tpr" \\
    -f "$PROD_DIR/md.xtc" \\
    -n index.ndx \\
    -oav tables/distance_avg.xvg \\
    -tu ns

echo "Distance analysis complete → tables/distance_avg.xvg"
"""

        (step_dir / "make_ndx.sh").write_text(make_ndx_script.strip())
        (step_dir / "run_analysis.sh").write_text(script.strip())
        self._write_metadata(step, step_dir, analysis_type, {
            "selection": selection,
            "requires_make_ndx": True,
            "outputs": ["tables/distance_avg.xvg"],
        })

    # ── Energy ────────────────────────────────────────────────────────────────

    def _build_energy(
        self,
        step:          SimulationStep,
        step_dir:      Path,
        prod_ref:      str,
        analysis_type: str,
    ) -> None:

        script = f"""#!/bin/bash
# ─── Energy terms ────────────────────────────────────────────────────────────
# Extracts: Potential, Kinetic En., Total Energy, Temperature, Pressure.
# No group selection needed — reads directly from .edr.

PROD_DIR="{prod_ref}"

mkdir -p tables

echo "Potential Kinetic-En. Total-Energy Temperature Pressure" | gmx energy \\
    -f "$PROD_DIR/md.edr" \\
    -o tables/energy.xvg

echo "Energy analysis complete → tables/energy.xvg"
"""

        (step_dir / "run_analysis.sh").write_text(script.strip())
        self._write_metadata(step, step_dir, analysis_type, {
            "terms": ["Potential", "Kinetic-En.", "Total-Energy", "Temperature", "Pressure"],
            "outputs": ["tables/energy.xvg"],
        })

    # ── Radius of gyration ────────────────────────────────────────────────────

    def _build_rgyr(
        self,
        step:          SimulationStep,
        step_dir:      Path,
        prod_ref:      str,
        analysis_type: str,
    ) -> None:

        script = f"""#!/bin/bash
# ─── Radius of gyration ───────────────────────────────────────────────────────
# Group 1 (Protein) — overall compactness over trajectory.

PROD_DIR="{prod_ref}"

mkdir -p tables plots

echo "1" | gmx gyrate \\
    -s "$PROD_DIR/md.tpr" \\
    -f "$PROD_DIR/md.xtc" \\
    -o tables/gyrate.xvg

echo "Radius of gyration complete → tables/gyrate.xvg"
"""

        (step_dir / "run_analysis.sh").write_text(script.strip())
        self._write_metadata(step, step_dir, analysis_type, {
            "group": "Protein (1)",
            "outputs": ["tables/gyrate.xvg"],
        })

    # ── Generic fallback ──────────────────────────────────────────────────────

    def _build_generic(
        self,
        step:          SimulationStep,
        step_dir:      Path,
        prod_ref:      str,
        analysis_type: str,
    ) -> None:

        script = f"""#!/bin/bash
# ─── Analysis: {analysis_type} ────────────────────────────────────────────────
# No GROMACS command template available for this analysis type.
# Implement manually using the trajectory files below.

PROD_DIR="{prod_ref}"

mkdir -p tables plots

# Trajectory:  $PROD_DIR/md.xtc
# Structure:   $PROD_DIR/md.tpr
# Energy:      $PROD_DIR/md.edr
# Log:         $PROD_DIR/md.log

echo "Implement analysis: {analysis_type}"
"""

        (step_dir / "run_analysis.sh").write_text(script.strip())
        self._write_metadata(step, step_dir, analysis_type, {
            "note": "no template available — implement manually",
        })

    # ── Shared helpers ────────────────────────────────────────────────────────

    def _write_metadata(
        self,
        step:          SimulationStep,
        step_dir:      Path,
        analysis_type: str,
        extra:         dict,
    ) -> None:

        metadata = {
            "step_id":       step.step_id,
            "stage":         step.stage.value,
            "engine":        step.engine,
            "step_type":     step.step_type.value,
            "blocking":      step.blocking,
            "generated_by":  "AnalysisBuilder",
            "analysis_type": analysis_type,
            "gromacs_groups": "built-in defaults (pdb2gmx standard)",
            "params":        step.params,
        }
        metadata.update(extra)

        (step_dir / "metadata.json").write_text(
            json.dumps(metadata, indent=4)
        )
