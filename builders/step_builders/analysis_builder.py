# builders/step_builders/analysis_builder.py

from __future__ import annotations

from pathlib import Path
import json

from core.execution_models import (
    SimulationStep,
)


# ═══════════════════════════════════════════════════════════════════════════════
# Analysis Builder
# ═══════════════════════════════════════════════════════════════════════════════

class AnalysisBuilder:
    """
    Genera scripts y configuración
    para análisis post-simulación.
    """

    def build(
        self,
        step: SimulationStep,
        step_dir: Path,
    ) -> None:

        analysis_type = (
            step.step_id
            .replace("analysis_", "")
        )

        # ────────────────────────────────────────────────────────────────────
        # Script
        # ────────────────────────────────────────────────────────────────────

        script = f"""
#!/bin/bash

echo "Running analysis: {analysis_type}"

mkdir -p outputs

# Placeholder analysis command

echo "Analysis complete"
"""

        script_path = (
            step_dir / "run_analysis.sh"
        )

        script_path.write_text(
            script.strip()
        )

        # ────────────────────────────────────────────────────────────────────
        # Config
        # ────────────────────────────────────────────────────────────────────

        config = {
            "analysis_type": (
                analysis_type
            ),

            "input_trajectory": (
                "md.xtc"
            ),

            "input_topology": (
                "md.tpr"
            ),

            "expected_outputs": [
                "plots/",
                "tables/",
                "outputs/",
            ],
        }

        config_path = (
            step_dir
            / "analysis_config.json"
        )

        config_path.write_text(
            json.dumps(
                config,
                indent=4,
            )
        )

        # ────────────────────────────────────────────────────────────────────
        # Output folders
        # ────────────────────────────────────────────────────────────────────

        (step_dir / "plots").mkdir(
            exist_ok=True
        )

        (step_dir / "tables").mkdir(
            exist_ok=True
        )

        (step_dir / "outputs").mkdir(
            exist_ok=True
        )