# core/variant_compiler.py
"""
Multi-system / variant workflow compiler.

Handles YAML files with a top-level `variants` key:

    variants:
      - id: WT
        file: WT.pdb
      - id: P419T
        file: P419T.pdb

    shared_workflow:
      forcefield: amber99sb
      water_model: tip3p
      production_ns: 100

Each variant gets its own workspace under:
    simforge_runs/<project>/<variant_id>/

A comparative workspace is also created:
    simforge_runs/<project>/comparative/

The variant compiler:
  1. Validates the shared_workflow config
  2. Generates a concrete SystemState per variant (patching the protein file)
  3. Compiles each variant independently
  4. Records a VariantManifest for comparative analysis

Architecture note: variants are compile-time expansion, NOT runtime.
The DAG stays single-system; comparison happens in a separate analysis step.
"""

from __future__ import annotations

import yaml
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from core.compiler_models import CompilationResult


# ═══════════════════════════════════════════════════════════════════════════════
# Variant manifest model
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class VariantSpec:
    variant_id:  str
    file:        str
    label:       Optional[str]  = None
    extra_patch: dict           = field(default_factory=dict)


@dataclass
class VariantManifest:
    project_name:    str
    variants:        list[VariantSpec]
    shared_workflow: dict
    workspaces:      dict[str, str]   = field(default_factory=dict)  # id → path
    errors:          dict[str, str]   = field(default_factory=dict)  # id → error


# ═══════════════════════════════════════════════════════════════════════════════
# YAML detection
# ═══════════════════════════════════════════════════════════════════════════════

def is_variant_yaml(yaml_path: Path) -> bool:
    """Return True if the YAML has a top-level `variants` key."""
    try:
        raw = yaml.safe_load(yaml_path.read_text())
        return isinstance(raw, dict) and "variants" in raw
    except Exception:
        return False


def parse_variant_yaml(yaml_path: Path) -> VariantManifest:
    """
    Parse a variants YAML into a VariantManifest.

    Expected structure:
        project:
          name: slc16a11_variants

        variants:
          - id: WT
            file: structures/WT.pdb
          - id: P419T
            file: structures/P419T.pdb
            label: "Pro419Thr mutant"

        shared_workflow:
          forcefield: amber99sb
          water_model: tip3p
          temperature_K: 310.0
          duration_ns: 100
          analysis:
            - rmsd
            - rmsf
            - radius_of_gyration
    """
    raw      = yaml.safe_load(yaml_path.read_text())
    base_dir = yaml_path.parent

    project_name = raw.get("project", {}).get("name", yaml_path.stem)

    raw_variants = raw.get("variants", [])
    if not raw_variants:
        raise ValueError("variants YAML must have at least one entry in `variants`")

    specs: list[VariantSpec] = []
    for entry in raw_variants:
        if "id" not in entry or "file" not in entry:
            raise ValueError(f"Each variant must have `id` and `file`. Got: {entry}")
        file_path = Path(entry["file"])
        if not file_path.is_absolute():
            file_path = base_dir / file_path
        specs.append(VariantSpec(
            variant_id  = entry["id"],
            file        = str(file_path),
            label       = entry.get("label"),
            extra_patch = {k: v for k, v in entry.items() if k not in ("id", "file", "label")},
        ))

    shared = raw.get("shared_workflow", {})
    if not shared:
        raise ValueError("variants YAML must have a `shared_workflow` section")

    return VariantManifest(
        project_name    = project_name,
        variants        = specs,
        shared_workflow = shared,
    )


# ═══════════════════════════════════════════════════════════════════════════════
# Variant YAML generator (expands each variant into a single-system YAML)
# ═══════════════════════════════════════════════════════════════════════════════

def _build_single_yaml(
    spec:        VariantSpec,
    shared:      dict,
    project_name: str,
) -> dict:
    """
    Build a standard single-system YAML dict for one variant.
    Merges shared_workflow fields into the standard SimForge YAML schema.
    """
    sw = shared  # alias for brevity

    # Environment block
    env: dict = {}
    if "water_model" in sw:
        env.setdefault("solvent", {})["water_model"] = sw["water_model"]
    if "temperature_K" in sw:
        env["temperature_K"] = sw["temperature_K"]
    if "duration_ns" in sw:
        env["duration_ns"] = sw["duration_ns"]
    env.setdefault("ions", {})["concentration"] = sw.get("ion_concentration", 0.154)

    # Forcefields
    ff: dict = {}
    if "forcefield" in sw:
        ff["protein"] = sw["forcefield"]
    if "ligand_forcefield" in sw:
        ff["ligands"] = sw["ligand_forcefield"]

    # Analysis
    analyses_raw = sw.get("analysis", ["rmsd", "rmsf", "energy"])
    analyses: list[dict] = []
    for a in analyses_raw:
        if isinstance(a, str):
            analyses.append({"type": a})
        elif isinstance(a, dict):
            analyses.append(a)

    # Simulation objectives
    objectives = sw.get("simulation_objectives", ["stability"])

    doc = {
        "project": {"name": f"{project_name}_{spec.variant_id}"},
        "components": [
            {
                "id":   "protein_1",
                "role": "protein",
                "file": spec.file,
            }
        ],
        "environment":           env,
        "forcefields":           ff,
        "simulation_objectives": objectives,
        "analysis":              analyses,
    }

    if "hardware" in sw:
        doc["hardware"] = sw["hardware"]

    return doc


# ═══════════════════════════════════════════════════════════════════════════════
# Compiler
# ═══════════════════════════════════════════════════════════════════════════════

def compile_variants(
    manifest:   VariantManifest,
    output_dir: str = "simforge_runs",
    no_build:   bool = False,
) -> VariantManifest:
    """
    Compile all variants.

    Each variant is compiled + materialized into:
        output_dir/<project_name>/<variant_id>/

    manifest.workspaces and manifest.errors are populated in place.
    Returns the updated manifest.
    """
    import tempfile
    import yaml as yaml_lib
    from core.compiler import SimulationCompiler
    from builders.workspace_builder import WorkspaceBuilder

    root_dir = Path(output_dir) / manifest.project_name

    for spec in manifest.variants:
        try:
            single_yaml = _build_single_yaml(
                spec, manifest.shared_workflow, manifest.project_name
            )

            # Write temp YAML for parsing
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".yaml", delete=False, prefix=f"simforge_{spec.variant_id}_"
            ) as tmp:
                yaml_lib.dump(single_yaml, tmp, default_flow_style=False)
                tmp_path = Path(tmp.name)

            try:
                from core.parser import parse_yaml
                state = parse_yaml(tmp_path)
            finally:
                tmp_path.unlink(missing_ok=True)

            result = SimulationCompiler().compile_from_state(state)

            if not no_build:
                workspace = WorkspaceBuilder().build(
                    result,
                    output_dir=str(root_dir),
                )
                manifest.workspaces[spec.variant_id] = str(workspace)

        except Exception as exc:
            manifest.errors[spec.variant_id] = str(exc)

    return manifest


# ═══════════════════════════════════════════════════════════════════════════════
# Comparative workspace generator
# ═══════════════════════════════════════════════════════════════════════════════

def build_comparative_workspace(
    manifest:   VariantManifest,
    output_dir: str = "simforge_runs",
) -> Path:
    """
    Create a comparative/ directory pointing to all variant workspaces.
    Writes a README and a comparative_plan.json for future analysis.
    """
    import json

    root_dir    = Path(output_dir) / manifest.project_name
    comp_dir    = root_dir / "comparative_analysis"
    comp_dir.mkdir(parents=True, exist_ok=True)

    # comparative_plan.json
    plan = {
        "project":    manifest.project_name,
        "variants":   [
            {
                "id":        s.variant_id,
                "label":     s.label or s.variant_id,
                "file":      s.file,
                "workspace": manifest.workspaces.get(s.variant_id, ""),
            }
            for s in manifest.variants
        ],
        "analyses": [
            "comparative_rmsd",
            "comparative_rmsf",
            "delta_radius_of_gyration",
            "clustering_comparison",
            "pca_comparison",
            "secondary_structure_comparison",
        ],
        "status": "pending",
    }

    (comp_dir / "comparative_plan.json").write_text(
        json.dumps(plan, indent=2)
    )

    # README
    readme_lines: list[str] = [
        f"# Comparative Analysis — {manifest.project_name}",
        "",
        "This directory contains the comparative analysis workspace for all variants.",
        "",
        "## Variants",
        "",
    ]
    for s in manifest.variants:
        ws = manifest.workspaces.get(s.variant_id, "not compiled")
        label = s.label or s.variant_id
        readme_lines.append(f"- **{s.variant_id}** ({label}): `{ws}`")

    if manifest.errors:
        readme_lines += [
            "",
            "## Compilation Errors",
            "",
        ]
        for vid, err in manifest.errors.items():
            readme_lines.append(f"- **{vid}**: {err}")

    readme_lines += [
        "",
        "## Planned Analyses",
        "",
        "- Comparative RMSD",
        "- Comparative RMSF",
        "- ΔRg (radius of gyration difference)",
        "- Clustering comparison",
        "- PCA comparison",
        "- Secondary structure comparison",
        "",
        "Run each variant simulation first, then execute comparative analysis.",
        "",
        "## Directory Structure",
        "",
        "```",
        f"{manifest.project_name}/",
    ]
    for s in manifest.variants:
        readme_lines.append(f"  {s.variant_id}/")
    readme_lines += [
        "  comparative_analysis/",
        "    comparative_plan.json",
        "    README.md",
        "    (reports generated post-run)",
        "```",
    ]

    (comp_dir / "README.md").write_text("\n".join(readme_lines))

    return comp_dir
