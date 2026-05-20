# pipelines/inhibition_pipeline.py

from __future__ import annotations

from pipelines.base_pipeline import (
    BasePipeline,
)

from core.models import SystemState

from core.execution_models import (
    SimulationPlan,
)

from core.decision_engine import (
    build_simulation_plan,
)


# ═══════════════════════════════════════════════════════════════════════════════
# Inhibition Study Pipeline
# ═══════════════════════════════════════════════════════════════════════════════

class InhibitionPipeline(BasePipeline):
    """
    Pipeline especializado para estudios de inhibición.

    FUTURO:
        - apo MD
        - pocket detection
        - docking
        - pose validation
        - complex MD
        - binding analysis
    """

    pipeline_type = "inhibition-study"

    def build_plan(
        self,
        state: SystemState,
    ) -> SimulationPlan:

        # ────────────────────────────────────────────────────────────────────
        # TEMPORALMENTE:
        # reutiliza workflow MD estándar
        # ────────────────────────────────────────────────────────────────────

        return build_simulation_plan(
            state
        )