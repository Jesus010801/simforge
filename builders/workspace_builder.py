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
        # Step folders
        # ────────────────────────────────────────────────────────────────────

        for i, step in enumerate(
            result.execution_order,
            start=1,
        ):

            step_dir = (
                steps_dir
                / f"{i:02d}_{step.step_id}"
            )

            step_dir.mkdir(
                exist_ok=True
            )
            # ────────────────────────────────────────────────────────────────
            # Step materialization
            # ────────────────────────────────────────────────────────────────
            builder = STEP_BUILDERS.get(
                step.stage.value
            )

            if builder is not None:

                builder.build(
                    step,
                    step_dir,
                )

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

        return root