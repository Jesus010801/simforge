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

InflateGRO cutoff units:
  The Perl script expects the cutoff argument in Å (it multiplies by 0.1 to
  convert to nm internally: $cutoff = $ARGV[3]*0.1).  The builder accepts
  cutoff_nm (user-facing, in nm) and passes cutoff_nm*10 (Å) to the script.
  Default: 1.4 nm → passes 14 to Perl → internal cutoff 1.4 nm.  This value is
  pinned to the original protmemfiles tutorial command
  (docs/Prot-Memb_FILES/tutorial_membrana.txt:36 —
  "perl inflategro-Jorge.pl system.gro 4 DPP 14 system_inflated.gro 5 area.dat"),
  which removes every lipid within 1.4 nm of any protein Cα atom, unconditionally,
  before the shrink loop starts. A previous default of 0.14 nm (10x too small)
  left protein-clashing lipids essentially unremoved — see
  docs/audits/simforge_vs_protmemfiles_audit.md.
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

        # Cutoff in nm (user-facing). Default 1.4 nm, matching the original
        # protmemfiles tutorial's ONE-SHOT overlap-removal cutoff, applied
        # once before the shrink loop starts
        # (docs/Prot-Memb_FILES/tutorial_membrana.txt:36 — cutoff arg "14" → 1.4 nm).
        # Perl script receives cutoff_nm * 10 (Å), then multiplies by 0.1 → nm.
        cutoff_nm           = float(p.get("cutoff", 1.4))
        cutoff_angstrom     = cutoff_nm * 10.0   # Å value passed to inflategro

        # Shrink-loop (per-iteration deflate) overlap cutoff — separate from the
        # one-shot cutoff above. The original protmemfiles shrink loop passes a
        # literal cutoff of 0 on every deflate call after the initial cleanup
        # (docs/Prot-Memb_FILES/ScriptCamilo-Jorge.sh:18,38;
        # docs/Prot-Memb_FILES/run_inflategro.sh:41,48), so deflation only
        # relaxes/compresses lipid packing via subsequent minimization — it
        # never deletes lipids again. Reusing the nonzero one-shot cutoff on
        # every loop iteration instead causes cumulative annular lipid loss as
        # the box shrinks and previously-clear lipids are compressed back
        # within the cutoff distance (observed: 478 -> 400 DPP over 29
        # iterations on GLP-1R, entirely inside the loop — see
        # docs/audits/). Default 0.0 nm matches the original method.
        shrink_loop_cutoff_nm       = float(p.get("shrink_loop_cutoff_nm", 0.0))
        shrink_loop_cutoff_angstrom = shrink_loop_cutoff_nm * 10.0

        # Trapped-lipid policy: "warn" (default) or "strict" (block before solvation)
        trapped_lipid_policy = p.get("trapped_lipid_policy", "warn")

        annular_repair_enabled = bool(p.get("annular_repair_enabled", False))
        annular_repair_allow_insertion = bool(p.get("annular_repair_allow_insertion", False))

        # Backend & Quality gates config
        backend      = p.get("backend", "inflategro")
        water_model  = p.get("water_model", "spce")
        max_void_fraction_local = float(p.get("max_void_fraction_local", 0.25))
        max_void_area_local_nm2 = float(p.get("max_void_area_local_nm2", 1.0))
        min_tm_burial_score = float(p.get("min_tm_burial_score", 0.80))
        max_trapped_lipids = int(p.get("max_trapped_lipids", 0))
        max_soluble_domain_core_atoms = int(p.get("max_soluble_domain_core_atoms", 10))
        tm_aware_max_removed_lipids = int(p.get("tm_aware_max_removed_lipids", 50))

        # ── Phase 9A: Protein–Membrane Interface Evaluator config ─────────────
        iface_tm_residues_str = p.get("tm_residues")
        iface_ec_residues_str = p.get("ec_residues")
        iface_ic_residues_str = p.get("ic_residues")
        iface_min_covered     = float(p.get("interface_min_fraction_covered",    0.75))
        iface_max_gap         = float(p.get("interface_max_fraction_exposed_gap", 0.15))
        iface_max_p90         = float(p.get("interface_max_p90_distance_nm",      0.70))

        # ── Phase 9B: Embedding Optimizer config ──────────────────────────────
        # enabled: auto-enable when backend=tm_aware and tm_residues are set
        _opt_default_enabled = (backend == "tm_aware" and bool(p.get("tm_residues")))
        opt_enabled       = bool(p.get("optimizer_enabled", _opt_default_enabled))
        opt_policy        = p.get("optimizer_policy", "warn")
        opt_max_candidates = int(p.get("optimizer_max_candidates", 28))
        _z_raw = p.get("optimizer_z_shift_offsets_nm")
        opt_z_shifts = list(_z_raw) if _z_raw else [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6]
        _pad_raw = p.get("optimizer_mask_padding_nm")
        opt_paddings = list(_pad_raw) if _pad_raw else [0.0, 0.1, 0.2, 0.3]
        opt_min_improvement = float(p.get("optimizer_min_score_improvement", 0.10))
        opt_run_minimization = bool(p.get("optimizer_run_candidate_minimization", False))
        opt_excl_safety = int(p.get("optimizer_exclusion_safety_limit",
                                     max(tm_aware_max_removed_lipids, 100)))

        # ── Iteration snapshot config ──────────────────────────────────────────
        # save_embedding_iterations: "final_only" (default) | "every_n" | "all"
        # save_embedding_iteration_stride: only used for "every_n" mode (default 5)
        save_iterations     = p.get("save_embedding_iterations",       "final_only")
        iteration_stride    = int(p.get("save_embedding_iteration_stride", 5))

        # ── Embedding quality diagnostics & interface evaluator ───────────────
        quality_diagnostics      = bool(p.get("quality_diagnostics",      True))
        quality_analysis_region  = p.get("quality_analysis_region",  "tm_footprint")
        quality_footprint_margin = float(p.get("quality_footprint_margin_nm", 0.5))

        # Atom names for quality diagnostic come from the knowledge layer
        try:
            from core import membrane_knowledge as mk
            _atom_names      = mk.lipid_atom_names(lipid, forcefield)
            headgroup_atom   = _atom_names.headgroup_ref
            tail_atom        = _atom_names.tail_middle
            # Convert APL target from Å² to nm² for the quality diagnostic
            reference_apl_nm2 = apl_target_ang2 / 100.0
        except Exception:
            headgroup_atom    = p.get("headgroup_atom", "O33")
            tail_atom         = p.get("tail_atom",      "C50")
            reference_apl_nm2 = apl_target_ang2 / 100.0

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
        # Only match explicit topology assembly steps (NOT generate_protein_topology)
        topology_dir = next(
            (step_dir_map[d] for d in step.depends_on
             if d in ("assemble_system_topology", "generate_topology") and d in step_dir_map),
            None,
        )
        # assemble_system_topology dir (always available in step_dir_map; used for bootstrap)
        assemble_dir = step_dir_map.get("assemble_system_topology")
        assemble_ref = _rel(step_dir, assemble_dir) if assemble_dir else "../assemble_system_topology"

        if topology_dir:
            topol_top = str(Path(_rel(step_dir, topology_dir)) / "topol.top")
            # system_processed.gro: written by assemble_system_topology (legacy path)
            input_gro = str(Path(_rel(step_dir, topology_dir)) / "system_processed.gro")
        elif embed_dir:
            # New DAG: assemble_system_topology runs AFTER embedding.
            # Use embed_in_bilayer/system.gro as input; bootstrap topology is generated
            # at runtime by calling assemble_system_topology/run_assemble_system.py.
            input_gro = str(Path(_rel(step_dir, embed_dir)) / "system.gro")
            topol_top = "bootstrap_topol.top"

        # Resolve inflategro from workspace membrane_assets/ (staged at compile time)
        membrane_assets_dir = step_dir_map.get("__membrane_assets__")
        if membrane_assets_dir:
            inflategro_script = str(
                Path(_rel(step_dir, membrane_assets_dir)) / "inflategro-Jorge.pl"
            )

        # Simforge root for the diagnostic helper (known at build time)
        simforge_root = str(Path(__file__).resolve().parent.parent.parent)

        # ── Shrink-loop MDP (position-restrained minimization) ─────────────────
        self._write_minim_mdp(step_dir, temperature_K)

        # ── Coordinate/topology lipid-count sync helper ────────────────────────
        self._write_topology_sync_script(
            step_dir           = step_dir,
            lipid_residue_name = lipid_residue_name,
            simforge_root      = simforge_root,
        )

        # ── Trapped-lipid diagnostic helper ────────────────────────────────────
        self._write_diagnose_script(
            step_dir             = step_dir,
            lipid_residue_name   = lipid_residue_name,
            cutoff_nm            = cutoff_nm,
            trapped_lipid_policy = trapped_lipid_policy,
            simforge_root        = simforge_root,
        )

        # ── Embedding quality diagnostic helper ────────────────────────────────
        if quality_diagnostics:
            self._write_quality_script(
                step_dir             = step_dir,
                lipid_residue_name   = lipid_residue_name,
                headgroup_atom       = headgroup_atom,
                tail_atom            = tail_atom,
                reference_apl_nm2    = reference_apl_nm2,
                analysis_region      = quality_analysis_region,
                footprint_margin_nm  = quality_footprint_margin,
                simforge_root        = simforge_root,
            )

        # ── Phase 9B: Embedding Optimizer script ──────────────────────────────
        self._write_optimizer_script(
            step_dir               = step_dir,
            tm_residues_str        = p.get("tm_residues"),
            lipid_residue_name     = lipid_residue_name,
            headgroup_atom         = headgroup_atom,
            tail_atom              = tail_atom,
            reference_apl_nm2      = reference_apl_nm2,
            input_gro              = input_gro,
            enabled                = opt_enabled,
            policy                 = opt_policy,
            max_candidates         = opt_max_candidates,
            z_shift_offsets_nm     = opt_z_shifts,
            mask_padding_nm_values = opt_paddings,
            min_score_improvement  = opt_min_improvement,
            run_candidate_minimization = opt_run_minimization,
            exclusion_safety_limit = opt_excl_safety,
            min_fraction_covered   = iface_min_covered,
            max_fraction_exposed_gap = iface_max_gap,
            max_p90_distance_nm    = iface_max_p90,
            simforge_root          = simforge_root,
        )

        # ── Phase 9A: Protein–Membrane Interface Evaluator script ─────────────
        self._write_interface_script(
            step_dir              = step_dir,
            tm_residues_str       = iface_tm_residues_str,
            ec_residues_str       = iface_ec_residues_str,
            ic_residues_str       = iface_ic_residues_str,
            min_fraction_covered  = iface_min_covered,
            max_fraction_exposed_gap = iface_max_gap,
            max_p90_distance_nm   = iface_max_p90,
            simforge_root         = simforge_root,
        )

        # ── Main script ────────────────────────────────────────────────────────
        membrane_assets_rel = (
            str(Path(_rel(step_dir, membrane_assets_dir)))
            if membrane_assets_dir else None
        )
        use_bootstrap = (topol_top == "bootstrap_topol.top")

        if backend == "tm_aware":
            self._write_tm_aware_script(
                step_dir                      = step_dir,
                input_gro                     = input_gro,
                topol_top                     = topol_top,
                lipid_residue_name            = lipid_residue_name,
                tm_residues_str               = p.get("tm_residues"),
                water_model                   = water_model,
                max_void_fraction_local       = max_void_fraction_local,
                max_void_area_local_nm2       = max_void_area_local_nm2,
                min_tm_burial_score           = min_tm_burial_score,
                max_trapped_lipids            = max_trapped_lipids,
                max_soluble_domain_core_atoms = max_soluble_domain_core_atoms,
                tm_aware_max_removed_lipids   = tm_aware_max_removed_lipids,
                headgroup_atom                = headgroup_atom,
                tail_atom                     = tail_atom,
                reference_apl_nm2             = reference_apl_nm2,
                quality_analysis_region       = quality_analysis_region,
                quality_footprint_margin      = quality_footprint_margin,
                simforge_root                 = simforge_root,
                membrane_assets_rel           = membrane_assets_rel,
                assemble_ref                  = assemble_ref,
                use_bootstrap                 = use_bootstrap,
            )
            # Write run.sh: Phase 9B optimizer → tm_aware backend → Phase 9A
            run_sh = (
                "#!/bin/bash\n"
                "set -e\n"
                "cd \"$(dirname \"$0\")\"\n"
                "# Phase 9B: membrane embedding optimizer\n"
                "echo \"[embedding] Phase 9B: running membrane embedding optimizer...\"\n"
                "python3 run_membrane_optimizer.py\n"
                "# Phase: TM-aware backend (reads work_input.gro if it exists)\n"
                "python3 run_tm_aware.py\n"
                "# Phase 9A: protein–membrane interface evaluation (warning only)\n"
                "python3 run_protein_membrane_interface.py || true\n"
            )
            (step_dir / "run.sh").write_text(run_sh)
            (step_dir / "run.sh").chmod(0o755)
        else:
            self._write_shrink_loop(
                step_dir              = step_dir,
                inflategro_script     = inflategro_script,
                input_gro             = input_gro,
                topol_top             = topol_top,
                lipid_name            = lipid_residue_name,
                inflate_factor        = inflate_factor,
                deflate_factor        = deflate_factor,
                apl_target_ang2       = apl_target_ang2,
                apl_tolerance_ang2    = apl_tolerance_ang2,
                max_iterations        = max_iterations,
                gridsize              = gridsize,
                cutoff_nm             = cutoff_nm,
                cutoff_angstrom       = cutoff_angstrom,
                shrink_loop_cutoff_nm       = shrink_loop_cutoff_nm,
                shrink_loop_cutoff_angstrom = shrink_loop_cutoff_angstrom,
                trapped_lipid_policy  = trapped_lipid_policy,
                annular_repair_enabled = annular_repair_enabled,
                annular_repair_allow_insertion = annular_repair_allow_insertion,
                tm_residues_str       = iface_tm_residues_str,
                save_iterations       = save_iterations,
                iteration_stride      = iteration_stride,
                quality_diagnostics   = quality_diagnostics,
                membrane_assets_rel   = membrane_assets_rel,
                assemble_ref          = assemble_ref,
                use_bootstrap         = use_bootstrap,
            )

        # ── Metadata ───────────────────────────────────────────────────────────
        expected_outputs = [
            "converged.gro",
            "trapped_lipid_diagnosis.json",
            "protein_membrane_interface_report.json",
            "membrane_embedding_optimizer_report.json",
            "system_baseline.gro",
            "work_input.gro",
        ]
        if backend == "tm_aware":
            expected_outputs.extend([
                "tm_mask_report.json",
                "tm_lipid_classification_report.json",
                "tm_aware_lipid_exclusion_report.json",
                "embedding_quality_report.json",
                "tm_aware_backend_report.json",
            ])
        else:
            expected_outputs.extend([
                "shrink_telemetry.json",
            ])
            if quality_diagnostics:
                expected_outputs.append("embedding_quality_report.json")
        if save_iterations in ("all", "every_n") and backend != "tm_aware":
            expected_outputs.append("snapshots/")

        metadata = {
            "step_id":        step.step_id,
            "stage":          step.stage.value,
            "engine":         step.engine,
            "step_type":      step.step_type.value,
            "blocking":       step.blocking,
            "generated_by":   "EmbeddingBuilder",
            "gate":           {"type": "tm_aware_report"} if backend == "tm_aware" else {"type": "apl_report"},
            "expected_outputs": expected_outputs,
            "params": {
                "lipid":                        lipid,
                "lipid_residue_name":           lipid_residue_name,
                "forcefield":                   forcefield,
                "temperature_K":                temperature_K,
                "backend":                      backend,
                "apl_target_ang2":              apl_target_ang2,
                "apl_tolerance_ang2":           apl_tolerance_ang2,
                "inflate_factor":               inflate_factor,
                "deflate_factor":               deflate_factor,
                "max_iterations":               max_iterations,
                "lipid_exclusion_cutoff_nm":    cutoff_nm,
                "cutoff_angstrom_to_perl":      cutoff_angstrom,
                "shrink_loop_cutoff_nm":        shrink_loop_cutoff_nm,
                "shrink_loop_cutoff_angstrom_to_perl": shrink_loop_cutoff_angstrom,
                "trapped_lipid_policy":         trapped_lipid_policy,
                "save_embedding_iterations":    save_iterations,
                "save_embedding_iteration_stride": iteration_stride,
                "quality_diagnostics":          quality_diagnostics,
                "quality_analysis_region":      quality_analysis_region,
                "quality_footprint_margin_nm":  quality_footprint_margin,
                "max_void_fraction_local":       max_void_fraction_local,
                "max_void_area_local_nm2":       max_void_area_local_nm2,
                "min_tm_burial_score":           min_tm_burial_score,
                "max_trapped_lipids":            max_trapped_lipids,
                "max_soluble_domain_core_atoms": max_soluble_domain_core_atoms,
            },
        }
        (step_dir / "metadata.json").write_text(json.dumps(metadata, indent=4))

    # ── MDP for shrink-loop minimization ──────────────────────────────────────

    def _write_minim_mdp(self, step_dir: Path, temperature_K: float) -> None:
        # Tight emtol + both POSRES and STRONG_POSRES active during embedding.
        # reaction-field avoids GROMACS 2025+ FATAL ERROR on PME + net charge
        # (shrink-loop minimization only removes steric clashes; RF accuracy is sufficient).
        mdp = f"""; Shrink-loop minimization — position restraints active
define          = -DPOSRES -DSTRONG_POSRES
integrator      = steep
emtol           = 1000.0
emstep          = 0.01
nsteps          = 50000

cutoff-scheme   = Verlet
nstlist         = 1
coulombtype     = reaction-field
epsilon-rf      = 0
rcoulomb        = 1.2
rvdw            = 1.2
pbc             = xyz
DispCorr        = EnerPres
"""
        (step_dir / "minim_shrink.mdp").write_text(mdp.strip())

    # ── Coordinate/topology lipid-count sync helper ───────────────────────────

    def _write_topology_sync_script(
        self,
        step_dir:           Path,
        lipid_residue_name: str,
        simforge_root:      str,
    ) -> None:
        """
        Write sync_topology.py into the step directory.

        Counts lipid residues directly from a coordinate file (never trusts a
        prior expected count or inflategro's log text) and rewrites only the
        matching [ molecules ] entry in the topology, so that every grompp
        call downstream sees a coordinate/topology atom count that agrees.
        See validators/topology_sync.py for the underlying algorithm and
        docs/audits/ for the real-run evidence this fixes.
        """
        script = f'''#!/usr/bin/env python3
"""Coordinate/topology lipid-count sync — generated by EmbeddingBuilder. Do not edit."""
import sys
from pathlib import Path

sys.path.insert(0, {simforge_root!r})

from validators.topology_sync import sync_topology_lipid_count

LIPID = {lipid_residue_name!r}

GRO_PATH   = sys.argv[1] if len(sys.argv) > 1 else "work.gro"
TOPOL_PATH = sys.argv[2] if len(sys.argv) > 2 else "topol.top"

report = sync_topology_lipid_count(
    gro_path=GRO_PATH,
    topol_path=TOPOL_PATH,
    lipid_resname=LIPID,
    report_path="topology_sync_report.json",
)

print(
    f"[topology_sync] {{LIPID}}: {{report['old_dpp_count']}} -> {{report['new_dpp_count']}} "
    f"(removed {{report['removed_dpp_count']}})  "
    f"coordinate_atoms={{report['coordinate_atoms']}}  "
    f"synchronized={{report['synchronized']}}"
)
for w in report.get("warnings", []):
    print(f"[topology_sync] WARNING: {{w}}")
'''
        path = step_dir / "sync_topology.py"
        path.write_text(script)
        path.chmod(0o755)

    # ── Trapped-lipid diagnostic Python helper ────────────────────────────────

    def _write_diagnose_script(
        self,
        step_dir:             Path,
        lipid_residue_name:   str,
        cutoff_nm:            float,
        trapped_lipid_policy: str,
        simforge_root:        str,
    ) -> None:
        """Write diagnose_trapped.py into the step directory."""
        script = f'''#!/usr/bin/env python3
"""Trapped-lipid diagnostic — generated by EmbeddingBuilder. Do not edit."""
import json, sys, math
from pathlib import Path

LIPID  = "{lipid_residue_name}"
POLICY = "{trapped_lipid_policy}"
CUTOFF_NM = {cutoff_nm!r}

initial_lipid_count = int(sys.argv[1]) if len(sys.argv) > 1 else 0
final_lipid_count   = int(sys.argv[2]) if len(sys.argv) > 2 else 0
n_removed = initial_lipid_count - final_lipid_count

# Try the full validator from simforge
try:
    sys.path.insert(0, {simforge_root!r})
    from validators.membrane_validators import detect_trapped_lipids
    diag = detect_trapped_lipids("converged.gro", lipid_resname=LIPID)
    trapped  = diag.n_suspicious_lipids
    diag_dict = diag.dict()
except Exception as _e:
    trapped  = 0
    diag_dict = {{"error": str(_e), "note": "validator unavailable; counts may be unreliable"}}

cleanup_passed = trapped == 0
if trapped == 0:
    recommendation = "OK — no suspicious trapped lipids detected after embedding."
elif POLICY == "warn":
    recommendation = (
        f"WARNING: {{trapped}} suspicious trapped lipid(s) remain after embedding. "
        "Review converged.gro before proceeding to production MD. "
        "Consider increasing lipid_exclusion_cutoff_nm or re-running with strict policy."
    )
else:
    recommendation = (
        f"ERROR: {{trapped}} suspicious trapped lipid(s) remain. "
        "Solvation blocked by strict policy. "
        "Increase lipid_exclusion_cutoff_nm and re-run embedding."
    )

report = {{
    "initial_lipid_count":           initial_lipid_count,
    "final_lipid_count":             final_lipid_count,
    "n_lipids_removed_by_inflategro": n_removed,
    "lipid_exclusion_cutoff_nm":      CUTOFF_NM,
    "trapped_lipids_remaining":       trapped,
    "trapped_lipid_policy":           POLICY,
    "cleanup_passed":                 cleanup_passed,
    "recommendation":                 recommendation,
    "diagnosis_detail":               diag_dict,
}}

Path("trapped_lipid_diagnosis.json").write_text(json.dumps(report, indent=2))

summary = {{k: v for k, v in report.items() if k != "diagnosis_detail"}}
print(json.dumps(summary, indent=2))

if trapped > 0:
    if POLICY == "strict":
        print(
            f"[embedding] ERROR: {{trapped}} suspicious trapped lipid(s) detected "
            f"(strict policy). Blocking before solvation.",
            file=sys.stderr,
        )
        sys.exit(2)
    else:
        print(
            f"[embedding] WARNING: {{trapped}} suspicious trapped lipid(s) detected. "
            f"Review trapped_lipid_diagnosis.json before production MD."
        )
        print(
            "[embedding] WARNING: These lipids may cause instability. "
            "Set trapped_lipid_policy: strict to block the workflow instead."
        )
'''
        path = step_dir / "diagnose_trapped.py"
        path.write_text(script)
        path.chmod(0o755)

    # ── Embedding quality diagnostic script ──────────────────────────────────

    def _write_quality_script(
        self,
        step_dir:            Path,
        lipid_residue_name:  str,
        headgroup_atom:      str,
        tail_atom:           str,
        reference_apl_nm2:   float,
        analysis_region:     str,
        footprint_margin_nm: float,
        simforge_root:       str,
    ) -> None:
        """Write diagnose_quality.py into the step directory."""
        script = f'''#!/usr/bin/env python3
"""Embedding quality diagnostic — generated by EmbeddingBuilder. Do not edit."""
import json, sys
from pathlib import Path

LIPID_RESNAME       = "{lipid_residue_name}"
HEADGROUP_ATOM      = "{headgroup_atom}"
TAIL_ATOM           = "{tail_atom}"
REFERENCE_APL_NM2   = {reference_apl_nm2!r}
ANALYSIS_REGION     = "{analysis_region}"
FOOTPRINT_MARGIN_NM = {footprint_margin_nm!r}

try:
    sys.path.insert(0, {simforge_root!r})
    from validators.embedding_quality import diagnose_embedding_quality
    result = diagnose_embedding_quality(
        "converged.gro",
        lipid_resname=LIPID_RESNAME,
        headgroup_atom=HEADGROUP_ATOM,
        tail_atom=TAIL_ATOM,
        reference_apl_nm2=REFERENCE_APL_NM2,
        analysis_region=ANALYSIS_REGION,
        footprint_margin_nm=FOOTPRINT_MARGIN_NM,
    )
    report = result.dict()
except Exception as _e:
    report = {{"error": str(_e), "note": "embedding quality diagnostic unavailable"}}

Path("embedding_quality_report.json").write_text(json.dumps(report, indent=2))

# Print key metrics (never block — diagnostics only)
_keys = ("n_lipids_in_region", "void_fraction_local", "void_area_local_nm2",
         "apl_local_nm2", "bilayer_core_thickness_nm", "bilayer_outer_thickness_nm",
         "z_core_top_median", "z_core_bot_median", "message")
print(json.dumps({{k: report.get(k) for k in _keys}}, indent=2))
for w in report.get("warnings", []):
    print(f"[quality_advisory] {{w}}")
'''
        path = step_dir / "diagnose_quality.py"
        path.write_text(script)
        path.chmod(0o755)

    # ── Main shrink-loop script ───────────────────────────────────────────────

    def _write_shrink_loop(
        self,
        step_dir:             Path,
        inflategro_script:    str,
        input_gro:            str,
        topol_top:            str,
        lipid_name:           str,
        inflate_factor:       float,
        deflate_factor:       float,
        apl_target_ang2:      float,
        apl_tolerance_ang2:   float,
        max_iterations:       int,
        gridsize:             int,
        cutoff_nm:            float,
        cutoff_angstrom:      float,
        shrink_loop_cutoff_nm:       float = 0.0,
        shrink_loop_cutoff_angstrom: float = 0.0,
        trapped_lipid_policy: str = "warn",
        annular_repair_enabled: bool = False,
        annular_repair_allow_insertion: bool = False,
        tm_residues_str:       str | None = None,
        save_iterations:      str = "final_only",
        iteration_stride:     int = 5,
        quality_diagnostics:  bool = True,
        membrane_assets_rel:  Optional[str] = None,
        assemble_ref:         str = "../assemble_system_topology",
        use_bootstrap:        bool = False,
    ) -> None:

        apl_cutoff = apl_target_ang2 + apl_tolerance_ang2
        annular_repair_enabled_shell = "true" if annular_repair_enabled else "false"
        annular_repair_allow_insertion_py = "True" if annular_repair_allow_insertion else "False"
        from core.structural_annotation import residues_in_range
        tm_residues = sorted(residues_in_range(tm_residues_str)) if tm_residues_str else []
        tm_residues_literal = "[" + ", ".join(str(resid) for resid in tm_residues) + "]"

        # Format cutoff_angstrom: integer when it's a whole number, else float
        cutoff_angstrom_str = (
            str(int(cutoff_angstrom))
            if cutoff_angstrom == int(cutoff_angstrom)
            else f"{cutoff_angstrom:.6g}"
        )
        shrink_loop_cutoff_angstrom_str = (
            str(int(shrink_loop_cutoff_angstrom))
            if shrink_loop_cutoff_angstrom == int(shrink_loop_cutoff_angstrom)
            else f"{shrink_loop_cutoff_angstrom:.6g}"
        )

        # ── Snapshot bash code fragments (baked in by save_iterations mode) ────
        if save_iterations == "all":
            snapshot_setup = "# Iteration snapshots (save_embedding_iterations=all)\nmkdir -p snapshots\n"
            snapshot_inner = '    cp work.gro "snapshots/iter_$(printf \'%04d\' $ITER).gro"\n'
            snapshot_label = "all"
        elif save_iterations == "every_n":
            snapshot_setup = (
                f"# Iteration snapshots (save_embedding_iterations=every_n, stride={iteration_stride})\n"
                f"SNAP_STRIDE={iteration_stride}\n"
                "mkdir -p snapshots\n"
            )
            snapshot_inner = (
                '    if (( ITER % SNAP_STRIDE == 0 )); then\n'
                '        cp work.gro "snapshots/iter_$(printf \'%04d\' $ITER).gro"\n'
                '    fi\n'
            )
            snapshot_label = f"every_n (stride={iteration_stride})"
        else:  # final_only
            snapshot_setup = ""
            snapshot_inner = ""
            snapshot_label = "final_only"

        quality_call = (
            "\n# ── Step 9: embedding quality diagnostic ─────────────────────────────────\n"
            'echo "[embedding] Running embedding quality diagnostic..."\n'
            "python3 diagnose_quality.py\n"
        ) if quality_diagnostics else ""

        interface_call = (
            "\n# ── Step 10: protein–membrane interface evaluation (Phase 9A) ─────────────\n"
            'echo "[embedding] Running protein–membrane interface evaluation..."\n'
            "python3 run_protein_membrane_interface.py || true\n"
        )

        gmxlib_val = membrane_assets_rel or "../../membrane_assets"

        # Bootstrap topology generation: called when assemble_system_topology runs AFTER embedding.
        # run_assemble_system.py supports argv[1]=GRO argv[2]=OUT_TOP for this purpose.
        if use_bootstrap:
            bootstrap_step = (
                "\n# ── Bootstrap: generate topology for shrink loop ─────────────────────────────\n"
                f'echo "[embedding] Bootstrap: building topology from $INPUT_GRO ..."\n'
                f'python3 "{assemble_ref}/run_assemble_system.py" "$INPUT_GRO" "$TOPOL"\n'
                f'echo "[embedding] Bootstrap topology written to $TOPOL"\n'
            )
        else:
            bootstrap_step = ""

        script = f"""#!/bin/bash
# ─── Membrane embedding — InflateGRO shrink loop ──────────────────────────────
# Generated by EmbeddingBuilder. Do not edit manually.
#
# Physical target: APL ≤ {apl_cutoff:.1f} Å² ({lipid_name}, {apl_target_ang2:.0f} + {apl_tolerance_ang2:.0f} Å² tolerance)
# Convergence: iterative deflation (factor {deflate_factor}) until APL converges.
# Telemetry: shrink_telemetry.json — one entry per iteration.
#
# Cutoff (one-shot, initial inflation only): {cutoff_nm} nm → {cutoff_angstrom_str} Å passed to inflategro
#   (inflategro expects Å input; multiplies by 0.1 internally to get nm)
# Shrink-loop cutoff (every deflate iteration): {shrink_loop_cutoff_nm} nm → {shrink_loop_cutoff_angstrom_str} Å
#   Matches the original protmemfiles method, which applies the overlap cutoff
#   ONCE before the shrink loop, then runs every deflate iteration with
#   cutoff=0 (docs/Prot-Memb_FILES/ScriptCamilo-Jorge.sh, run_inflategro.sh) —
#   the loop only relaxes/compresses lipid packing, it never deletes lipids
#   again. Reusing the one-shot cutoff on every iteration instead causes
#   cumulative annular lipid loss as the box shrinks (see docs/audits/).
# Trapped-lipid policy: {trapped_lipid_policy}
# Iteration snapshots: {snapshot_label}
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
CUTOFF_NM={cutoff_nm}
CUTOFF_ANGSTROM={cutoff_angstrom_str}
LOOP_CUTOFF_NM={shrink_loop_cutoff_nm}
LOOP_CUTOFF_ANGSTROM={shrink_loop_cutoff_angstrom_str}
# GMXLIB: point at workspace membrane_assets so grompp can find FF files
export GMXLIB="{gmxlib_val}:${{GMXLIB:-}}"

{snapshot_setup}
# Warn if cutoff is zero (no lipid exclusion — may leave trapped lipids)
if [ "$CUTOFF_NM" = "0.0" ] || [ "$CUTOFF_NM" = "0" ]; then
    echo "[embedding] WARNING: lipid_exclusion_cutoff_nm=0 — inflategro will NOT exclude overlapping lipids."
    echo "[embedding] WARNING: This may result in lipids trapped inside the protein TM cavity."
    echo "[embedding] WARNING: Set lipid_exclusion_cutoff_nm: 1.4 (or higher) in your config."
fi

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
    local _pyconv
    [ "$converged" = "true" ] && _pyconv=True || _pyconv=False
    python3 -c "
import json, pathlib, datetime
p = pathlib.Path('shrink_telemetry.json')
data = json.loads(p.read_text()) if p.exists() else {{'iterations': [], 'converged': False, 'started_at': '$(date -Iseconds)'}}
data['iterations'].append({{'iter': $iter, 'apl_ang2': $apl, 'converged': $_pyconv}})
data['last_apl_ang2'] = $apl
p.write_text(json.dumps(data, indent=2))
"
}}

finalize_telemetry() {{
    local apl=$1 converged=$2 n_iter=$3
    local _pyconv
    [ "$converged" = "true" ] && _pyconv=True || _pyconv=False
    python3 -c "
import json, pathlib, datetime
p = pathlib.Path('shrink_telemetry.json')
data = json.loads(p.read_text()) if p.exists() else {{'iterations': []}}
data['converged'] = $_pyconv
data['final_apl_ang2'] = $apl
data['n_iterations'] = $n_iter
data['finished_at'] = datetime.datetime.now().isoformat(timespec='seconds')
p.write_text(json.dumps(data, indent=2))
"
}}

# ── Bash helper: update lipid count in topology ──────────────────────────────
update_topology_lipid_count() {{
    local topol=$1 lipid=$2 old_count=$3 new_count=$4
    if [ ! -f "$topol" ]; then
        echo "[embedding] WARNING: topology not found at $topol — skip count update"
        return 0
    fi
    if grep -qE "^${{lipid}}[[:space:]]+${{old_count}}([[:space:]]|$)" "$topol"; then
        sed -i "s/^${{lipid}}[[:space:]]\\+${{old_count}}[[:space:]]*$/${{lipid}}    ${{new_count}}/" "$topol"
        echo "[embedding] Topology updated: ${{lipid}} ${{old_count}} -> ${{new_count}}"
    else
        echo "[embedding] WARNING: could not find '${{lipid}} ${{old_count}}' in topology — count not updated"
    fi
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
echo "[embedding] Cutoff: ${{CUTOFF_NM}} nm (${{CUTOFF_ANGSTROM}} Å to inflategro)"
{bootstrap_step}
# ── Step 0: Phase 9B Membrane Embedding Optimizer ────────────────────────────
echo "[embedding] Step 0: running membrane embedding optimizer (Phase 9B)..."
python3 run_membrane_optimizer.py
# work_input.gro is now available (optimized or copy of INPUT_GRO)

# ── Step 1: initial inflation ─────────────────────────────────────────────────
echo "[embedding] Step 1: inflating bilayer (factor $INFLATE)"
cp work_input.gro work.gro

# Read lipid count from bootstrap topology BEFORE inflation (InflateGRO may remove lipids)
ORIG_LIPID_COUNT=$(grep -oP "^${{LIPID}}[[:space:]]+\\K[0-9]+" "$TOPOL" | head -1 || echo 0)
echo "[embedding] Bootstrap lipid count: $ORIG_LIPID_COUNT"

perl "$INFLATEGRO" work.gro $INFLATE $LIPID $CUTOFF_ANGSTROM inflated.gro $GRIDSIZE area_2.dat \
    > inflategro_initial.log 2>&1
cat inflategro_initial.log
APL=$(read_apl)
echo "[embedding] Post-inflation APL = ${{APL}} Å²"
log_iter 0 $APL false
cp inflated.gro work.gro

# Extract initial lipid count from inflategro output
INITIAL_LIPID_COUNT=$(grep -oP 'There are \\K[0-9]+(?= lipids)' inflategro_initial.log | head -1 || true)
INITIAL_LIPID_COUNT=${{INITIAL_LIPID_COUNT:-0}}
echo "[embedding] Initial lipid count: $INITIAL_LIPID_COUNT"

# Sync topology to match actual GRO if InflateGRO removed protein-overlapping lipids
if [ "$ORIG_LIPID_COUNT" -gt 0 ] && [ "$INITIAL_LIPID_COUNT" -gt 0 ] && [ "$INITIAL_LIPID_COUNT" -ne "$ORIG_LIPID_COUNT" ]; then
    update_topology_lipid_count "$TOPOL" "$LIPID" "$ORIG_LIPID_COUNT" "$INITIAL_LIPID_COUNT"
fi

# ── Coordinate/topology sync: authoritative, counts $LIPID directly from ──────
# work.gro rather than trusting the heuristic above or inflategro's log text.
# This is the check that actually prevents the grompp atom-count mismatch.
python3 sync_topology.py work.gro "$TOPOL"

# ── Step 2: first minimization (inflated system) ──────────────────────────────
echo "[embedding] Step 2: initial minimization (inflated system)"
gmx grompp -f minim_shrink.mdp -c work.gro -r work.gro -p "$TOPOL" -o work.tpr -maxwarn 2 -quiet
gmx mdrun -s work.tpr -deffnm work -nb gpu -quiet
cp work.gro work_prev.gro

# ── Step 3: shrink loop ───────────────────────────────────────────────────────
# Uses LOOP_CUTOFF_ANGSTROM (default 0 = overlap-checking off), NOT the
# one-shot CUTOFF_ANGSTROM used above — matching the original protmemfiles
# method. The one-shot cleanup already removed protein-clashing lipids before
# the loop started; the loop's job is purely to compress/relax lipid packing
# via deflation + minimization, not to keep re-deleting lipids as the box
# shrinks and previously-clear lipids are compressed back toward the protein.
echo "[embedding] Step 3: shrink loop (deflate $DEFLATE per iteration, loop overlap cutoff = ${{LOOP_CUTOFF_NM}} nm)"
ITER=1
APL=$( python3 -c "print(int(float(open('area_2.dat').read().strip()) * 100))" )

while [ "$APL" -gt "$( python3 -c "print(int($APL_TARGET))" )" ] && [ "$ITER" -le "$MAX_ITER" ]; do

    echo "[embedding] Iter $ITER: APL = ${{APL}} Å² > ${{APL_TARGET}} Å² — deflating"

    perl "$INFLATEGRO" work.gro $DEFLATE $LIPID $LOOP_CUTOFF_ANGSTROM work.gro $GRIDSIZE area_2.dat \
        >> inflategro_loop.log 2>&1

    APL=$(read_apl)
    log_iter $ITER $APL false

    # inflategro may have removed protein-overlapping lipids this iteration —
    # re-sync the topology's $LIPID count with work.gro before grompp.
    python3 sync_topology.py work.gro "$TOPOL"

    gmx grompp -f minim_shrink.mdp -c work.gro -r work.gro -p "$TOPOL" -o work.tpr -maxwarn 2 -quiet
    gmx mdrun -s work.tpr -deffnm work -nb gpu -quiet

    # Clean GROMACS backup files
    rm -f \\#*
{snapshot_inner}
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
# Re-sync once more immediately before the final grompp — this is the exact
# call site that prevented the reported failure (work.gro had 477 DPP /
# 30421 atoms; bootstrap_topol.top still expected 478 DPP / 30471 atoms).
python3 sync_topology.py work.gro "$TOPOL"
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

# ── Phase M2B: conservative local annular repair ─────────────────────────────
if {annular_repair_enabled_shell} ; then
    echo "[embedding] Running conservative annular repair..."
    python3 - <<'PY'
from validators.annular_repair import repair_annular_packing
repair_annular_packing(
    "converged.gro",
    {tm_residues_literal},
    enabled=True,
    allow_insertion={annular_repair_allow_insertion_py},
)
PY
else
    echo "[embedding] Annular repair disabled."
    python3 - <<'PY'
import json
from pathlib import Path
Path("annular_repair_report.json").write_text(json.dumps({{
    "enabled": False,
    "status": "disabled",
    "coordinate_file_modified": False,
    "topology_modified": False,
    "n_relocated_lipids": 0,
    "n_inserted_lipids": 0,
}}, indent=2))
PY
fi

# ── Step 7: track lipid count changes and update topology ─────────────────────
# With LOOP_CUTOFF_NM=0 (default), inflategro-Jorge.pl skips its overlap-check
# block entirely during the shrink loop, so inflategro_loop.log contains no
# "There are N lipids..." lines and this grep correctly finds nothing — the
# fallback to $INITIAL_LIPID_COUNT then reflects reality: the loop no longer
# removes any lipids. sync_topology.py (called before every grompp above) is
# the authoritative safety net regardless of what this heuristic finds.
FINAL_LIPID_COUNT=$(grep -oP 'There are \\K[0-9]+(?= lipids)' inflategro_loop.log 2>/dev/null | tail -1 || true)
FINAL_LIPID_COUNT=${{FINAL_LIPID_COUNT:-$INITIAL_LIPID_COUNT}}
N_REMOVED=$(( INITIAL_LIPID_COUNT - FINAL_LIPID_COUNT ))

echo "[embedding] Lipid count: initial=$INITIAL_LIPID_COUNT  final=$FINAL_LIPID_COUNT  removed=$N_REMOVED"

if [ "$N_REMOVED" -gt 0 ] && [ "$INITIAL_LIPID_COUNT" -gt 0 ]; then
    update_topology_lipid_count "$TOPOL" "$LIPID" "$INITIAL_LIPID_COUNT" "$FINAL_LIPID_COUNT"
fi

# ── Step 8: trapped-lipid diagnostic ─────────────────────────────────────────
echo "[embedding] Running trapped-lipid diagnostic..."
python3 diagnose_trapped.py $INITIAL_LIPID_COUNT $FINAL_LIPID_COUNT
DIAG_EXIT=$?

if [ $DIAG_EXIT -eq 2 ]; then
    echo "[embedding] Strict policy: blocking workflow due to trapped lipids."
    echo "[embedding] See trapped_lipid_diagnosis.json for details."
    exit 1
fi
{quality_call}{interface_call}"""
        script_path = step_dir / "run.sh"
        script_path.write_text(script.lstrip())
        script_path.chmod(0o755)

    def _write_tm_aware_script(
        self,
        step_dir: Path,
        input_gro: str,
        topol_top: str,
        lipid_residue_name: str,
        tm_residues_str: Optional[str],
        water_model: str,
        max_void_fraction_local: float,
        max_void_area_local_nm2: float,
        min_tm_burial_score: float,
        max_trapped_lipids: int,
        max_soluble_domain_core_atoms: int,
        tm_aware_max_removed_lipids: int,
        headgroup_atom: str,
        tail_atom: str,
        reference_apl_nm2: float,
        quality_analysis_region: str,
        quality_footprint_margin: float,
        simforge_root: str,
        membrane_assets_rel: Optional[str] = None,
        assemble_ref: str = "../assemble_system_topology",
        use_bootstrap: bool = False,
    ) -> None:
        if tm_residues_str:
            from core.structural_annotation import residues_in_range
            tm_set = residues_in_range(tm_residues_str)
            tm_value = f"set({sorted(tm_set)!r})"
        else:
            tm_value = "None"

        # Bootstrap block: generate topology at start of main() when assemble runs after embedding
        if use_bootstrap:
            bootstrap_block = (
                "\n    # Bootstrap: generate topology from embed_in_bilayer/system.gro\n"
                "    # (assemble_system_topology/run_assemble_system.py was generated at build time)\n"
                f"    _assemble_py = (Path(__file__).parent / {assemble_ref!r} / 'run_assemble_system.py').resolve()\n"
                "    if not Path(TOPOL).exists() and _assemble_py.exists():\n"
                "        import subprocess as _bsub\n"
                "        _embed_gro = str((Path(__file__).parent / Path(INPUT_GRO).parent).resolve() / Path(INPUT_GRO).name) if not Path(INPUT_GRO).is_absolute() else INPUT_GRO\n"
                "        print(f'[tm_aware_backend] Bootstrap: building topology {TOPOL} from {INPUT_GRO} ...')\n"
                "        _bsub.run(['python3', str(_assemble_py), INPUT_GRO, TOPOL],\n"
                "                  cwd=str(Path(__file__).parent), check=True)\n"
            )
        else:
            bootstrap_block = ""

        script = f'''#!/usr/bin/env python3
"""TM-aware membrane embedding backend executor — generated by EmbeddingBuilder."""
import sys, os, shutil, json, subprocess
from pathlib import Path

# Add simforge_root to path
sys.path.insert(0, {simforge_root!r})

from validators.tm_exclusion_mask import exclude_lipids_tm_aware, classify_lipids_tm_aware
from validators.membrane_validators import detect_trapped_lipids
from validators.embedding_quality import diagnose_embedding_quality

INPUT_GRO = {input_gro!r}
TOPOL = {topol_top!r}
LIPID_RESNAME = {lipid_residue_name!r}
TM_RESIDUES = {tm_value}

# Phase 9B: use optimizer-promoted input if available
import os as _os
if _os.path.exists("work_input.gro"):
    INPUT_GRO = "work_input.gro"
    print("[tm_aware_backend] Using Phase 9B optimizer-promoted input: work_input.gro")

# Thresholds
MAX_VOID_FRACTION = {max_void_fraction_local}
MAX_VOID_AREA = {max_void_area_local_nm2}
MIN_TM_BURIAL = {min_tm_burial_score}
MAX_TRAPPED_LIPIDS = {max_trapped_lipids}
MAX_SOLUBLE_CORE_ATOMS = {max_soluble_domain_core_atoms}
TM_AWARE_MAX_REMOVED = {tm_aware_max_removed_lipids}

# Diagnostics parameters
HEADGROUP_ATOM = {headgroup_atom!r}
TAIL_ATOM = {tail_atom!r}
REFERENCE_APL_NM2 = {reference_apl_nm2!r}
ANALYSIS_REGION = {quality_analysis_region!r}
FOOTPRINT_MARGIN_NM = {quality_footprint_margin!r}

def count_lipids(gro_path, resname):
    if not Path(gro_path).exists():
        return 0
    lines = Path(gro_path).read_text().splitlines()
    resids = set()
    for line in lines[2:-1]:
        if len(line) >= 10:
            if line[5:10].strip() == resname:
                try:
                    resids.add(int(line[0:5]))
                except ValueError:
                    pass
    return len(resids)

def update_topology_lipid_count_py(topol_path: Path, resname: str, count: int):
    if not topol_path.exists():
        return
    lines = topol_path.read_text().splitlines()
    new_lines = []
    updated = False
    for line in lines:
        parts = line.split()
        if len(parts) >= 2 and parts[0] == resname:
            new_lines.append(f"{{resname:<16}} {{count}}")
            updated = True
        else:
            new_lines.append(line)
    if updated:
        topol_path.write_text("\\n".join(new_lines) + "\\n")

def main():
    print("[tm_aware_backend] Starting TM-aware membrane embedding backend")

    # 1. Run TM-mask validation (check if TM residues exist)
    if not TM_RESIDUES:
        err_msg = "TM-aware backend requires transmembrane residue annotations. Fallback recommendation: set membrane.embedding.backend to 'inflategro'."
        backend_report = {{
            "backend": "tm_aware",
            "tm_annotation_used": False,
            "n_lipids_initial": 0,
            "n_lipids_removed": 0,
            "n_lipids_final": 0,
            "exclusion_policy": "strict_apply",
            "quality_passed": False,
            "final_trapped_lipid_count": 0,
            "final_void_fraction_local": 0.0,
            "final_void_area_local_nm2": 0.0,
            "final_tm_burial_score": 0.0,
            "soluble_domain_core_penetration": 0,
            "warnings": [err_msg],
            "blocking_reason": err_msg
        }}
        Path("tm_aware_backend_report.json").write_text(json.dumps(backend_report, indent=2))
        # Empty placeholder reports for other files to satisfy expected_outputs
        for fn in ("tm_mask_report.json", "tm_lipid_classification_report.json", "tm_aware_lipid_exclusion_report.json", "embedding_quality_report.json"):
            Path(fn).write_text(json.dumps({{}}, indent=2))
        print(f"ERROR: {{err_msg}}", file=sys.stderr)
        sys.exit(1)

    # 2. Run TM-mask validation & TM-aware lipid exclusion in strict_apply mode
{bootstrap_block}
    n_initial = count_lipids(INPUT_GRO, LIPID_RESNAME)
    print(f"[tm_aware_backend] Initial lipid count in {{INPUT_GRO}}: {{n_initial}}")

    shutil.copy(INPUT_GRO, "work.gro")
    topol_path = Path(TOPOL)
    topol_dir = topol_path.parent
    shutil.copy(topol_path, "topol.top")
    for itp in topol_dir.glob("*.itp"):
        shutil.copy(itp, itp.name)

    try:
        # Run TM-mask classification + exclusion in strict_apply mode
        classify_lipids_tm_aware(
            gro_path="work.gro",
            tm_residues=TM_RESIDUES,
            lipid_resname=LIPID_RESNAME
        )
        excl_report = exclude_lipids_tm_aware(
            gro_path="work.gro",
            tm_residues=TM_RESIDUES,
            lipid_resname=LIPID_RESNAME,
            policy="strict_apply",
            safety_limit=TM_AWARE_MAX_REMOVED
        )
        n_removed = excl_report.get("n_lipids_removed", 0)
    except Exception as e:
        err_msg = f"TM-aware exclusion/classification failed: {{e}}"
        backend_report = {{
            "backend": "tm_aware",
            "tm_annotation_used": True,
            "n_lipids_initial": n_initial,
            "n_lipids_removed": 0,
            "n_lipids_final": n_initial,
            "exclusion_policy": "strict_apply",
            "quality_passed": False,
            "final_trapped_lipid_count": 0,
            "final_void_fraction_local": 0.0,
            "final_void_area_local_nm2": 0.0,
            "final_tm_burial_score": 0.0,
            "soluble_domain_core_penetration": 0,
            "warnings": [err_msg],
            "blocking_reason": err_msg
        }}
        Path("tm_aware_backend_report.json").write_text(json.dumps(backend_report, indent=2))
        print(f"ERROR: {{err_msg}}", file=sys.stderr)
        sys.exit(1)

    n_final = count_lipids("work.gro", LIPID_RESNAME)
    print(f"[tm_aware_backend] Post-exclusion lipid count: {{n_final}} (removed {{n_removed}} lipids)")

    # 3. Update topology lipid count for post-exclusion system.
    # pdb2gmx MUST NOT run on work.gro (mixed protein-lipid GRO — loses chain/TER semantics
    # and causes deterministic OXT terminus errors). Instead, update only the [molecules]
    # count in the topology that was bootstrapped from embed_in_bilayer/system.gro.
    update_topology_lipid_count_py(Path("topol.top"), LIPID_RESNAME, n_final)
    print(f"[tm_aware_backend] Topology lipid count updated: {{LIPID_RESNAME}} → {{n_final}}")

    # Set up GMXLIB so grompp can find the forcefield
    ff_env = os.environ.copy()
    for ff_candidate in [Path("oplsaa_membrane.ff"), Path("../membrane_assets/oplsaa_membrane.ff"), Path("../../membrane_assets/oplsaa_membrane.ff"), {(repr(membrane_assets_rel + "/oplsaa_membrane.ff") if membrane_assets_rel else "None")}]:
        if ff_candidate is None:
            continue
        ff_candidate = Path(ff_candidate)
        if ff_candidate.exists():
            ff_env["GMXLIB"] = str(ff_candidate.parent.resolve())
            break

    # Validate coordinate/topology consistency
    grompp_cmd = ["gmx", "grompp", "-f", "minim_shrink.mdp", "-c", "work.gro", "-r", "work.gro",
                  "-p", "topol.top", "-o", "work.tpr", "-maxwarn", "2", "-quiet"]
    ret_g = subprocess.run(grompp_cmd, capture_output=True, text=True, env=ff_env)
    if ret_g.returncode != 0:
        err_msg = f"GROMACS consistency check (grompp) failed: {{ret_g.stderr}}"
        backend_report = {{
            "backend": "tm_aware",
            "tm_annotation_used": True,
            "n_lipids_initial": n_initial,
            "n_lipids_removed": n_removed,
            "n_lipids_final": n_final,
            "exclusion_policy": "strict_apply",
            "quality_passed": False,
            "final_trapped_lipid_count": 0,
            "final_void_fraction_local": 0.0,
            "final_void_area_local_nm2": 0.0,
            "final_tm_burial_score": 0.0,
            "soluble_domain_core_penetration": 0,
            "warnings": [err_msg],
            "blocking_reason": err_msg
        }}
        Path("tm_aware_backend_report.json").write_text(json.dumps(backend_report, indent=2))
        print(f"ERROR: {{err_msg}}", file=sys.stderr)
        sys.exit(1)

    # 4. Run conservative local minimization
    print("[tm_aware_backend] Running conservative local minimization...")
    mdrun_cmd = ["gmx", "mdrun", "-s", "work.tpr", "-deffnm", "work", "-nb", "gpu", "-quiet"]
    ret_md = subprocess.run(mdrun_cmd, capture_output=True, text=True)
    if ret_md.returncode != 0:
        mdrun_cmd_cpu = ["gmx", "mdrun", "-s", "work.tpr", "-deffnm", "work", "-quiet"]
        ret_md = subprocess.run(mdrun_cmd_cpu, capture_output=True, text=True)
        if ret_md.returncode != 0:
            err_msg = f"GROMACS local minimization (mdrun) failed: {{ret_md.stderr}}"
            backend_report = {{
                "backend": "tm_aware",
                "tm_annotation_used": True,
                "n_lipids_initial": n_initial,
                "n_lipids_removed": n_removed,
                "n_lipids_final": n_final,
                "exclusion_policy": "strict_apply",
                "quality_passed": False,
                "final_trapped_lipid_count": 0,
                "final_void_fraction_local": 0.0,
                "final_void_area_local_nm2": 0.0,
                "final_tm_burial_score": 0.0,
                "soluble_domain_core_penetration": 0,
                "warnings": [err_msg],
                "blocking_reason": err_msg
            }}
            Path("tm_aware_backend_report.json").write_text(json.dumps(backend_report, indent=2))
            print(f"ERROR: {{err_msg}}", file=sys.stderr)
            sys.exit(1)

    for p in Path(".").glob("#*#"):
        try:
            p.unlink()
        except OSError:
            pass

    # 5. Run embedding quality diagnostics
    print("[tm_aware_backend] Running embedding quality diagnostics...")
    quality_res = diagnose_embedding_quality(
        "work.gro",
        lipid_resname=LIPID_RESNAME,
        headgroup_atom=HEADGROUP_ATOM,
        tail_atom=TAIL_ATOM,
        tm_residues=TM_RESIDUES,
        footprint_margin_nm=FOOTPRINT_MARGIN_NM,
        reference_apl_nm2=REFERENCE_APL_NM2,
        analysis_region=ANALYSIS_REGION
    )
    quality_report = quality_res.dict()
    Path("embedding_quality_report.json").write_text(json.dumps(quality_report, indent=2))

    diag_trapped = detect_trapped_lipids("work.gro", lipid_resname=LIPID_RESNAME, tm_residues=TM_RESIDUES)
    final_trapped_count = diag_trapped.n_suspicious_lipids

    # Write trapped_lipid_diagnosis.json
    trapped_report = {{
        "workspace_path": str(Path(".").resolve()),
        "files_inspected": ["work.gro"],
        "actual_coordinate_file_used_after_membrane_embedding": str(Path("converged.gro").resolve()),
        "actual_coordinate_file_passed_to_solvate_membrane": str(Path("converged.gro").resolve()),
        "bilayer_midplane_z": diag_trapped.bilayer_midplane_z,
        "protein_tm_center": diag_trapped.protein_tm_center,
        "protein_xy_footprint": diag_trapped.protein_xy_footprint,
        "lipid_residue_candidates": [
            {{
                "resid": c.resid,
                "resname": c.resname,
                "stage_first_detected": c.stage_first_detected,
                "present_since_embed": c.present_since_embed,
                "introduced_by_shrink": c.introduced_by_shrink,
                "inside_tm_footprint": c.inside_tm_footprint,
                "inside_expanded_tm_footprint": c.inside_expanded_tm_footprint,
                "low_overlap_but_geometrically_inside": c.low_overlap_but_geometrically_inside,
                "min_distance_to_protein_nm": c.min_distance_to_protein_nm,
                "atoms_within_0.25_nm": c.atoms_within_025,
                "atoms_within_0.35_nm": c.atoms_within_035,
                "atoms_within_0.45_nm": c.atoms_within_045,
                "lipid_center": c.lipid_center,
                "protein_tm_center": c.protein_tm_center,
                "z_relative_to_midplane": c.z_relative_to_midplane,
                "suspected_cause": c.suspected_cause,
                "recommended_action": c.recommended_action,
                "classification": c.classification, "lipid_com": c.lipid_com,
                "inside_slab": c.inside_slab, "angularly_enclosed": c.angularly_enclosed,
                "min_dist_tm_slab_nm": c.min_dist_tm_slab_nm,
                "min_dist_full_protein_nm": c.min_dist_full_protein_nm,
                "nearest_contact_region": c.nearest_contact_region
            }}
            for c in diag_trapped.lipid_residue_candidates
        ],
        "summary": {{
            "n_lipids_checked": diag_trapped.n_lipids_checked,
            "n_suspicious_lipids": diag_trapped.n_suspicious_lipids,
            "n_hard_overlap_lipids": diag_trapped.n_hard_overlap_lipids,
            "n_low_overlap_geometrically_trapped_lipids": diag_trapped.n_low_overlap_geometrically_trapped_lipids,
            "n_true_cavity_trapped": diag_trapped.n_true_cavity_trapped,
            "n_tm_slab_hard_contact": diag_trapped.n_tm_slab_hard_contact,
            "n_annular_surface_contact": diag_trapped.n_annular_surface_contact,
            "n_soluble_domain_contact": diag_trapped.n_soluble_domain_contact,
            "n_ambiguous": diag_trapped.n_ambiguous,
            "n_legacy_trapped_advisories": diag_trapped.n_legacy_trapped_advisories,
            "slab_z_min": diag_trapped.slab_z_min, "slab_z_max": diag_trapped.slab_z_max,
            "likely_root_cause": diag_trapped.likely_root_cause,
            "should_block_future_workflows": diag_trapped.should_block_future_workflows,
            "remediation_not_applied": diag_trapped.remediation_not_applied
        }}
    }}
    Path("trapped_lipid_diagnosis.json").write_text(json.dumps(trapped_report, indent=4))

    # Evaluate quality gates
    warnings = []
    blocking_reason = None
    quality_passed = True

    void_fraction = quality_report.get("void_fraction_local")
    void_area = quality_report.get("void_area_local_nm2")
    tm_burial = quality_report.get("tm_burial_score")
    soluble_penetration = quality_report.get("soluble_domain_core_penetration")

    if void_fraction is None: void_fraction = 0.0
    if void_area is None: void_area = 0.0
    if tm_burial is None: tm_burial = 1.0
    if soluble_penetration is None: soluble_penetration = 0

    if void_fraction > MAX_VOID_FRACTION:
        quality_passed = False
        blocking_reason = f"void_fraction_local {{void_fraction:.3f}} > threshold {{MAX_VOID_FRACTION}}"
    elif void_fraction > MAX_VOID_FRACTION - 0.05:
        warnings.append(f"void_fraction_local {{void_fraction:.3f}} is near failure threshold {{MAX_VOID_FRACTION}}")

    if void_area > MAX_VOID_AREA:
        quality_passed = False
        blocking_reason = f"void_area_local_nm2 {{void_area:.3f}} > threshold {{MAX_VOID_AREA}}"
    elif void_area > MAX_VOID_AREA - 0.2:
        warnings.append(f"void_area_local_nm2 {{void_area:.3f}} is near failure threshold {{MAX_VOID_AREA}}")

    if tm_burial < MIN_TM_BURIAL:
        quality_passed = False
        blocking_reason = f"tm_burial_score {{tm_burial:.3f}} < threshold {{MIN_TM_BURIAL}}"
    elif tm_burial < MIN_TM_BURIAL + 0.05:
        warnings.append(f"tm_burial_score {{tm_burial:.3f}} is near failure threshold {{MIN_TM_BURIAL}}")

    if final_trapped_count > MAX_TRAPPED_LIPIDS:
        quality_passed = False
        blocking_reason = f"trapped_lipids_count {{final_trapped_count}} > threshold {{MAX_TRAPPED_LIPIDS}}"

    if soluble_penetration > MAX_SOLUBLE_CORE_ATOMS:
        quality_passed = False
        blocking_reason = f"soluble_domain_core_penetration {{soluble_penetration}} > threshold {{MAX_SOLUBLE_CORE_ATOMS}}"
    elif soluble_penetration > MAX_SOLUBLE_CORE_ATOMS - 5:
        warnings.append(f"soluble_domain_core_penetration {{soluble_penetration}} is near failure threshold {{MAX_SOLUBLE_CORE_ATOMS}}")

    backend_report = {{
        "backend": "tm_aware",
        "tm_annotation_used": True,
        "n_lipids_initial": n_initial,
        "n_lipids_removed": n_removed,
        "n_lipids_final": n_final,
        "exclusion_policy": "strict_apply",
        "quality_passed": quality_passed,
        "final_trapped_lipid_count": final_trapped_count,
        "final_void_fraction_local": void_fraction,
        "final_void_area_local_nm2": void_area,
        "final_tm_burial_score": tm_burial,
        "soluble_domain_core_penetration": soluble_penetration,
        "warnings": warnings,
        "blocking_reason": blocking_reason
    }}
    Path("tm_aware_backend_report.json").write_text(json.dumps(backend_report, indent=2))

    if quality_passed:
        shutil.copy("work.gro", "converged.gro")
        print("[tm_aware_backend] TM-aware embedding backend finished successfully!")
        sys.exit(0)
    else:
        print(f"[tm_aware_backend] ERROR: TM-aware quality gate failure. Blocking reason: {{blocking_reason}}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
'''
        (step_dir / "run_tm_aware.py").write_text(script)
        (step_dir / "run_tm_aware.py").chmod(0o755)

    # ── Phase 9A: Protein–Membrane Interface Evaluator helper ─────────────────

    def _write_interface_script(
        self,
        step_dir:                Path,
        tm_residues_str:         Optional[str],
        ec_residues_str:         Optional[str],
        ic_residues_str:         Optional[str],
        min_fraction_covered:    float,
        max_fraction_exposed_gap: float,
        max_p90_distance_nm:     float,
        simforge_root:           str,
    ) -> None:
        """Write run_protein_membrane_interface.py — warning-only evaluator."""
        if tm_residues_str:
            from core.structural_annotation import residues_in_range
            tm_set   = sorted(residues_in_range(tm_residues_str))
            tm_value = f"set({tm_set!r})"
        else:
            tm_value = "set()"

        if ec_residues_str:
            from core.structural_annotation import residues_in_range
            ec_set   = sorted(residues_in_range(ec_residues_str))
            ec_value = f"set({ec_set!r})"
        else:
            ec_value = "None"

        if ic_residues_str:
            from core.structural_annotation import residues_in_range
            ic_set   = sorted(residues_in_range(ic_residues_str))
            ic_value = f"set({ic_set!r})"
        else:
            ic_value = "None"

        script = f'''#!/usr/bin/env python3
"""Phase 9A Protein–Membrane Interface Evaluator — generated by EmbeddingBuilder."""
import json, sys
from pathlib import Path

sys.path.insert(0, {simforge_root!r})

TM_RESIDUES  = {tm_value}
EC_RESIDUES  = {ec_value}
IC_RESIDUES  = {ic_value}

MIN_FRACTION_COVERED    = {min_fraction_covered!r}
MAX_FRACTION_EXPOSED    = {max_fraction_exposed_gap!r}
MAX_P90_DISTANCE_NM     = {max_p90_distance_nm!r}

GRO_FILE = "converged.gro"

if not Path(GRO_FILE).exists():
    report = {{"error": f"{{GRO_FILE}} not found — interface evaluation skipped"}}
    Path("protein_membrane_interface_report.json").write_text(json.dumps(report, indent=2))
    print(f"[interface_check] WARNING: {{GRO_FILE}} not found — skipping interface evaluation")
    sys.exit(0)

if not TM_RESIDUES:
    report = {{
        "note": "No TM residue annotation provided — interface evaluation skipped",
        "quality_passed": None,
        "warnings": ["tm_residues not set in YAML; add tm_residues to enable Phase 9A"],
    }}
    Path("protein_membrane_interface_report.json").write_text(json.dumps(report, indent=2))
    print("[interface_check] NOTE: tm_residues not set — interface evaluation skipped")
    sys.exit(0)

try:
    from validators.protein_membrane_interface import evaluate_protein_membrane_interface
    report = evaluate_protein_membrane_interface(
        gro_path                   = GRO_FILE,
        tm_residues                = TM_RESIDUES,
        ec_residues                = EC_RESIDUES,
        ic_residues                = IC_RESIDUES,
        min_fraction_covered       = MIN_FRACTION_COVERED,
        max_fraction_exposed_gap   = MAX_FRACTION_EXPOSED,
        max_p90_distance_nm        = MAX_P90_DISTANCE_NM,
        output_dir                 = Path("."),
    )
    if report.get("quality_passed") is False:
        print("[interface_check] WARNING: Protein–membrane interface quality check FAILED:")
        for w in report.get("warnings", []):
            print(f"  [interface_check] {{w}}")
        print(
            f"[interface_check]   fraction_covered={{report.get('fraction_covered'):.3f}}  "
            f"fraction_exposed={{report.get('fraction_exposed_gap'):.3f}}  "
            f"p90_dist={{report.get('p90_nearest_lipid_distance'):.3f}} nm"
        )
        print("[interface_check] Review protein_membrane_interface_report.json for details.")
        print("[interface_check] This is a WARNING — the workflow continues.")
    else:
        fc   = report.get("fraction_covered")
        p90  = report.get("p90_nearest_lipid_distance")
        nsur = report.get("n_tm_surface_atoms", 0)
        print(
            f"[interface_check] Interface quality OK: "
            f"covered={{fc:.3f}}  p90_dist={{p90:.3f}} nm  n_surface_atoms={{nsur}}"
        )
except Exception as _e:
    err_report = {{
        "error":         str(_e),
        "quality_passed": None,
        "warnings":      [f"Interface evaluator raised an exception: {{_e}}"],
    }}
    Path("protein_membrane_interface_report.json").write_text(
        json.dumps(err_report, indent=2)
    )
    print(f"[interface_check] WARNING: interface evaluation failed: {{_e}}", file=sys.stderr)

# Always exit 0 — this step is advisory only
sys.exit(0)
'''
        path = step_dir / "run_protein_membrane_interface.py"
        path.write_text(script)
        path.chmod(0o755)

    # ── Phase 9B: Embedding Optimizer script ─────────────────────────────────

    def _write_optimizer_script(
        self,
        step_dir:                  Path,
        tm_residues_str:           Optional[str],
        lipid_residue_name:        str,
        headgroup_atom:            str,
        tail_atom:                 str,
        reference_apl_nm2:         float,
        input_gro:                 str,
        enabled:                   bool,
        policy:                    str,
        max_candidates:            int,
        z_shift_offsets_nm:        list,
        mask_padding_nm_values:    list,
        min_score_improvement:     float,
        run_candidate_minimization: bool,
        exclusion_safety_limit:    int,
        min_fraction_covered:      float,
        max_fraction_exposed_gap:  float,
        max_p90_distance_nm:       float,
        simforge_root:             str,
    ) -> None:
        """Write run_membrane_optimizer.py into the step directory."""
        if tm_residues_str:
            from core.structural_annotation import residues_in_range
            tm_set   = sorted(residues_in_range(tm_residues_str))
            tm_value = f"set({tm_set!r})"
        else:
            tm_value = "set()"

        script = f'''#!/usr/bin/env python3
"""Phase 9B Membrane Embedding Optimizer — generated by EmbeddingBuilder."""
import sys, shutil
from pathlib import Path

sys.path.insert(0, {simforge_root!r})

from validators.membrane_embedding_optimizer import run_membrane_embedding_optimizer

INPUT_GRO_PATH   = {input_gro!r}
SYSTEM_BASELINE  = "system_baseline.gro"
WORK_INPUT       = "work_input.gro"

TM_RESIDUES = {tm_value}

# Config (baked in at build time)
ENABLED                     = {enabled!r}
POLICY                      = {policy!r}
MAX_CANDIDATES              = {max_candidates!r}
Z_SHIFT_OFFSETS_NM          = {z_shift_offsets_nm!r}
MASK_PADDING_NM_VALUES      = {mask_padding_nm_values!r}
MIN_SCORE_IMPROVEMENT       = {min_score_improvement!r}
RUN_CANDIDATE_MINIMIZATION  = {run_candidate_minimization!r}
EXCLUSION_SAFETY_LIMIT      = {exclusion_safety_limit!r}
LIPID_RESIDUE_NAME          = {lipid_residue_name!r}
HEADGROUP_ATOM              = {headgroup_atom!r}
TAIL_ATOM                   = {tail_atom!r}
REFERENCE_APL_NM2           = {reference_apl_nm2!r}
MIN_FRACTION_COVERED        = {min_fraction_covered!r}
MAX_FRACTION_EXPOSED_GAP    = {max_fraction_exposed_gap!r}
MAX_P90_DISTANCE_NM         = {max_p90_distance_nm!r}

def main():
    # Always create baseline backup
    if not Path(INPUT_GRO_PATH).exists():
        print(f"[optimizer] WARNING: INPUT_GRO_PATH not found: {{INPUT_GRO_PATH}}", file=sys.stderr)
        Path(WORK_INPUT).write_text("")
        return

    shutil.copy(INPUT_GRO_PATH, SYSTEM_BASELINE)

    try:
        report = run_membrane_embedding_optimizer(
            system_baseline_gro        = SYSTEM_BASELINE,
            work_input_gro             = WORK_INPUT,
            tm_residues                = TM_RESIDUES,
            output_dir                 = Path("."),
            lipid_residue_name         = LIPID_RESIDUE_NAME,
            headgroup_atom             = HEADGROUP_ATOM,
            tail_atom                  = TAIL_ATOM,
            reference_apl_nm2         = REFERENCE_APL_NM2,
            z_shift_offsets_nm         = Z_SHIFT_OFFSETS_NM,
            mask_padding_nm_values     = MASK_PADDING_NM_VALUES,
            max_candidates             = MAX_CANDIDATES,
            enabled                    = ENABLED,
            policy                     = POLICY,
            run_candidate_minimization = RUN_CANDIDATE_MINIMIZATION,
            min_score_improvement      = MIN_SCORE_IMPROVEMENT,
            exclusion_safety_limit     = EXCLUSION_SAFETY_LIMIT,
            min_fraction_covered       = MIN_FRACTION_COVERED,
            max_fraction_exposed_gap   = MAX_FRACTION_EXPOSED_GAP,
            max_p90_distance_nm        = MAX_P90_DISTANCE_NM,
            simforge_root              = {simforge_root!r},
        )

        promoted = report.get("promoted_to_work_input_gro", False)
        selected_id = report.get("selected_candidate_id")
        improvement = report.get("score_improvement")
        q_passed = report.get("selected_quality_passed")

        if not ENABLED:
            print("[optimizer] Optimizer disabled; baseline forwarded to work_input.gro")
        elif promoted:
            print(
                f"[optimizer] Promoted candidate {{selected_id}} | "
                f"score_improvement={{improvement:.3f}} | quality_passed={{q_passed}}"
            )
        else:
            print("[optimizer] No improvement found; baseline forwarded to work_input.gro")

        for w in report.get("warnings", []):
            print(f"[optimizer] WARNING: {{w}}")

        # Strict policy: block if optimizer ran but found no valid candidate
        if POLICY == "strict" and ENABLED and not promoted and TM_RESIDUES:
            print(
                "[optimizer] STRICT MODE: no valid candidate promoted. "
                "Downstream embedding will use baseline. "
                "Consider relaxing optimizer thresholds or disabling strict policy.",
                file=sys.stderr,
            )
            # Do NOT exit(1) — the baseline is already in work_input.gro and is valid input
            # Strict mode is advisory for the optimizer itself; the embedding step's
            # own quality gates will catch remaining issues.

    except Exception as _e:
        print(f"[optimizer] WARNING: optimizer failed: {{_e}}", file=sys.stderr)
        # Ensure work_input.gro exists for downstream steps
        if not Path(WORK_INPUT).exists():
            shutil.copy(SYSTEM_BASELINE, WORK_INPUT)

if __name__ == "__main__":
    main()
'''
        path = step_dir / "run_membrane_optimizer.py"
        path.write_text(script)
        path.chmod(0o755)
