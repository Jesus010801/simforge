# builders/workspace_builder.py

from __future__ import annotations

from pathlib import Path
import json
import shutil
import uuid
from datetime import datetime

from core.compiler_models import (
    CompilationResult,
)
from builders.step_builders.minimization_builder import (
    MinimizationBuilder,
)
from builders.builder_registry import (
    STEP_BUILDERS,
)
from core.workspace_fingerprint import (
    compute_build_signature,
    compute_builder_signature,
    compute_template_hash,
    SIMFORGE_VERSION,
)

# ═══════════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════════

def _stamp_step(step_dir: Path, generated_at: str) -> None:
    """Add provenance fields (template_hash, generated_at) to step metadata.json."""
    meta_path = step_dir / "metadata.json"
    if not meta_path.exists():
        return
    try:
        meta = json.loads(meta_path.read_text())
        scripts = list(step_dir.glob("*.sh")) + list(step_dir.glob("*.py"))
        meta["template_hash"] = compute_template_hash(scripts)
        meta["generated_at"]  = generated_at
        meta_path.write_text(json.dumps(meta, indent=4))
    except Exception:
        pass


# ═══════════════════════════════════════════════════════════════════════════════
# Workspace Builder
# ═══════════════════════════════════════════════════════════════════════════════

class WorkspaceBuilder:
    """
    Construye un workspace físico reproducible
    desde un CompilationResult.
    """

    def build(
        self,
        result: CompilationResult,
        output_dir: str = "simforge_runs",
        workspace_name: str | None = None,
        yaml_source: str = "",
        workspace_path: "str | Path | None" = None,
    ) -> Path:
        """
        Args:
            result:          CompilationResult from the compiler.
            output_dir:      Root directory for all workspaces.
            workspace_name:  Override the workspace directory name.
                             Defaults to inferred_system_type.
            yaml_source:     Absolute path to the original YAML config.
            workspace_path:  Explicit run directory (overrides output_dir +
                             workspace_name).  Used by the CLI to pass a
                             pre-created timestamped run directory so that
                             each compile lands in its own immutable location.
        """
        if workspace_path is not None:
            root = Path(workspace_path)
        else:
            system_name = workspace_name or (
                result.state.inferred_system_type
                or "simforge_system"
            )
            root = Path(output_dir) / system_name

        # ────────────────────────────────────────────────────────────────────
        # Core directories
        # ────────────────────────────────────────────────────────────────────

        workflow_dir = root / "workflow"

        reports_dir = root / "reports"

        metadata_dir = root / "metadata"

        steps_dir = root / "steps"

        workflow_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        reports_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        metadata_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        steps_dir.mkdir(parents=True, exist_ok=True)

        # ────────────────────────────────────────────────────────────────────
        # inputs/ — stage all external source files into the workspace
        # ────────────────────────────────────────────────────────────────────

        inputs_dir = root / "inputs"
        inputs_dir.mkdir(exist_ok=True)
        self._stage_inputs(result, inputs_dir)

        # ────────────────────────────────────────────────────────────────────
        # Step folders — pass 1: create dirs + build step_dir_map
        # ────────────────────────────────────────────────────────────────────

        step_dir_map: dict[str, Path] = {
            "__workspace_root__": root.resolve(),
        }

        for i, step in enumerate(result.execution_order, start=1):
            step_dir = steps_dir / f"{i:02d}_{step.step_id}"
            step_dir.mkdir(exist_ok=True)
            step_dir_map[step.step_id] = step_dir.resolve()

        # ────────────────────────────────────────────────────────────────────
        # Step materialization — pass 2: builders receive full map
        # ────────────────────────────────────────────────────────────────────

        generated_at = datetime.now().isoformat(timespec="seconds")

        for i, step in enumerate(result.execution_order, start=1):
            step_dir = steps_dir / f"{i:02d}_{step.step_id}"
            builder = STEP_BUILDERS.get(step.stage.value)
            if builder is not None:
                builder.build(step, step_dir, step_dir_map)
                _stamp_step(step_dir, generated_at)

        # ────────────────────────────────────────────────────────────────────
        # Mermaid workflow
        # ────────────────────────────────────────────────────────────────────

        mermaid_path = (
            workflow_dir
            / "workflow.mmd"
        )

        mermaid_path.write_text(
            result.mermaid_graph
        )

        # ────────────────────────────────────────────────────────────────────
        # User workflow
        # ────────────────────────────────────────────────────────────────────

        workflow_txt = (
            workflow_dir
            / "workflow.txt"
        )

        lines = []

        for i, step in enumerate(
            result.user_view,
            start=1,
        ):

            lines.append(
                f"{i:02d}. {step}"
            )

        workflow_txt.write_text(
            "\n".join(lines)
        )

        # ────────────────────────────────────────────────────────────────────
        # Summary metadata
        # ────────────────────────────────────────────────────────────────────

        summary_json = (
            metadata_dir
            / "summary.json"
        )

        summary_data = {
            "system_type": (
                result.state.inferred_system_type
            ),

            "workflow_steps": (
                len(result.plan.steps)
            ),

            "blocking_issues": (
                len(result.plan.blocking_issues)
            ),

            "special_protocols": (
                result.plan.special_protocols
            ),
        }

        summary_json.write_text(
            json.dumps(
                summary_data,
                indent=4,
            )
        )

        # ────────────────────────────────────────────────────────────────────
        # Execution manifest — fuente de verdad para el executor
        # ────────────────────────────────────────────────────────────────────

        manifest_entries = []
        for i, step in enumerate(result.execution_order, start=1):
            manifest_entries.append({
                "step_id":         step.step_id,
                "dir_name":        f"{i:02d}_{step.step_id}",
                "stage":           step.stage.value,
                "step_type":       step.step_type.value,
                "automation_level": step.effective_automation_level().value,
                "blocking":        step.blocking,
                "depends_on":      step.depends_on,
            })

        compiled_at       = datetime.now().isoformat(timespec="seconds")
        builder_sig       = compute_builder_signature()
        build_sig         = compute_build_signature(yaml_source) if yaml_source else builder_sig

        manifest_data = {
            "compile_id":         str(uuid.uuid4()),
            "compiled_at":        compiled_at,
            "simforge_version":   SIMFORGE_VERSION,
            "yaml_source":        str(Path(yaml_source).resolve()) if yaml_source else "",
            "builder_signature":  builder_sig,
            "build_signature":    build_sig,
            "system_type":        result.state.inferred_system_type,
            "n_steps":            len(manifest_entries),
            "steps":              manifest_entries,
        }

        manifest_path = metadata_dir / "execution_manifest.json"
        manifest_path.write_text(json.dumps(manifest_data, indent=4))

        # ────────────────────────────────────────────────────────────────────
        # Compile report (professional scientific summary)
        # ────────────────────────────────────────────────────────────────────
        try:
            from core.report_generator import generate_compile_report
            generate_compile_report(result, root)
        except Exception:
            pass   # report failure must never block workspace creation

        return root

    # ────────────────────────────────────────────────────────────────────────
    # Input staging
    # ────────────────────────────────────────────────────────────────────────

    def _stage_inputs(
        self,
        result: "CompilationResult",
        inputs_dir: Path,
    ) -> None:
        """
        Copy all external source files declared in component definitions
        into workspace/inputs/.

        Raises FileNotFoundError at compile time (not runtime) if any
        declared source file is missing — workspaces must be self-contained.

        Files are stored as inputs/{component_id}{original_extension}.
        """
        for comp in result.state.components:
            if not comp.file:
                continue

            src = Path(comp.file)

            if not src.exists():
                raise FileNotFoundError(
                    f"Source file for component '{comp.id}' not found: {src}\n"
                    f"  Check the 'file:' field in your YAML config.\n"
                    f"  Path was resolved relative to the YAML directory."
                )

            ext = src.suffix or ".pdb"
            dst = inputs_dir / f"{comp.id}{ext}"

            if src.resolve() != dst.resolve():
                shutil.copy2(src, dst)