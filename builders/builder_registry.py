# builders/builder_registry.py

from __future__ import annotations

from builders.step_builders.minimization_builder import (
    MinimizationBuilder,
)

from builders.step_builders.equilibration_builder import (
    EquilibrationBuilder,
)

from builders.step_builders.production_builder import (
    ProductionBuilder,
)

from builders.step_builders.analysis_builder import (
    AnalysisBuilder,
)




# ═══════════════════════════════════════════════════════════════════════════════
# Builder Registry
# ═══════════════════════════════════════════════════════════════════════════════

STEP_BUILDERS = {
    "minimization": MinimizationBuilder(),

    "equilibration": EquilibrationBuilder(),

    "production": ProductionBuilder(),

    "analysis": AnalysisBuilder(),
}
