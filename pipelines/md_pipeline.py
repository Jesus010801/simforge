# pipelines/md_pipeline.py

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
# Molecular Dynamics Pipeline
# ═══════════════════════════════════════════════════════════════════════════════

class MDPipeline(BasePipeline):

    pipeline_type = "molecular-dynamics"

    def build_plan(
        self,
        state: SystemState,
    ) -> SimulationPlan:

        return build_simulation_plan(
            state
        )