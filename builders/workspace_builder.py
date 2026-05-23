# builders/workspace_builder.py

from __future__ import annotations

from pathlib import Path
import json
import shutil

from core.compiler_models import (
    CompilationResult,
)
from builders.step_builders.minimization_builder import (
    MinimizationBuilder,
)

from builders.builder_registry import (
    STEP_BUILDERS,
)

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
    ) -> Path:

        system_name = (
            result.state.inferred_system_type
            or "simforge_system"
        )

        root = (
            Path(output_dir)
            / system_name
        )

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

        if steps_dir.exists():
            shutil.rmtree(steps_dir)
        steps_dir.mkdir(parents=True)

        # ────────────────────────────────────────────────────────────────────
        # Step folders — pass 1: create dirs + build step_dir_map
        # ────────────────────────────────────────────────────────────────────

        step_dir_map: dict[str, Path] = {}

        for i, step in enumerate(result.execution_order, start=1):
            step_dir = steps_dir / f"{i:02d}_{step.step_id}"
            step_dir.mkdir(exist_ok=True)
            step_dir_map[step.step_id] = step_dir.resolve()

        # ────────────────────────────────────────────────────────────────────
        # Step materialization — pass 2: builders receive full map
        # ────────────────────────────────────────────────────────────────────

        for i, step in enumerate(result.execution_order, start=1):
            step_dir = steps_dir / f"{i:02d}_{step.step_id}"
            builder = STEP_BUILDERS.get(step.stage.value)
            if builder is not None:
                builder.build(step, step_dir, step_dir_map)

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
                "step_id":   step.step_id,
                "dir_name":  f"{i:02d}_{step.step_id}",
                "stage":     step.stage.value,
                "step_type": step.step_type.value,
                "blocking":  step.blocking,
                "depends_on": step.depends_on,
            })

        manifest_data = {
            "compiled_at": __import__("datetime").datetime.now().isoformat(timespec="seconds"),
            "system_type": result.state.inferred_system_type,
            "n_steps":     len(manifest_entries),
            "steps":       manifest_entries,
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