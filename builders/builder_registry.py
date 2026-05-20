# builders/builder_registry.py

from __future__ import annotations

from builders.step_builders.preparation_builder import (
    PreparationBuilder,
)

from builders.step_builders.parametrization_builder import (
    ParametrizationBuilder,
)

from builders.step_builders.assembly_builder import (
    AssemblyBuilder,
)

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

from builders.step_builders.validation_builder import (
    ValidationBuilder,
)

from builders.step_builders.enhanced_sampling_builder import (
    EnhancedSamplingBuilder,
)


# ═══════════════════════════════════════════════════════════════════════════════
# Builder Registry
# ═══════════════════════════════════════════════════════════════════════════════

STEP_BUILDERS = {
    "preparation":       PreparationBuilder(),
    "parametrization":   ParametrizationBuilder(),
    "validation":        ValidationBuilder(),
    "assembly":          AssemblyBuilder(),
    "minimization":      MinimizationBuilder(),
    "equilibration":     EquilibrationBuilder(),
    "enhanced_sampling": EnhancedSamplingBuilder(),
    "production":        ProductionBuilder(),
    "analysis":          AnalysisBuilder(),
}