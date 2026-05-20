# pipelines/base_pipeline.py

from __future__ import annotations

from abc import ABC, abstractmethod

from core.models import SystemState

from core.execution_models import (
    SimulationPlan,
)


# ═══════════════════════════════════════════════════════════════════════════════
# Base Pipeline
# ═══════════════════════════════════════════════════════════════════════════════

class BasePipeline(ABC):
    """
    Pipeline científico abstracto.

    Cada pipeline define:
        - cómo interpretar el sistema
        - qué workflow construir
        - qué análisis ejecutar
    """

    pipeline_type: str = "base"

    @abstractmethod
    def build_plan(
        self,
        state: SystemState,
    ) -> SimulationPlan:
        """
        Construye un SimulationPlan.
        """

        raise NotImplementedError