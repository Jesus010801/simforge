"""Shared scientific intermediates and the dependency DAG (Phase 13).

Observables declare ``IntermediateRequest``s via ``ObservableSpec.dependencies``;
intermediates may declare their own.  The planner builds one deduplicated,
acyclic, deterministically ordered graph; the executor reuses or computes each
node once through the persistent, identity-addressed store and hands the
results to consumers via ``AnalysisContext.intermediates``.
"""
from analysis.campaign.intermediates.base import (  # noqa: F401
    INTERMEDIATE_SCHEMA, IntermediateContext, IntermediateRequest, IntermediateResult,
    IntermediateScope, IntermediateSpec, IntermediateStatus,
)
from analysis.campaign.intermediates.executor import (  # noqa: F401
    consumer_dependencies, resolve_plan,
)
from analysis.campaign.intermediates.planner import (  # noqa: F401
    DependencyCycleError, DependencyPlan, plan_dependencies,
)
from analysis.campaign.intermediates.store import (  # noqa: F401
    IntermediateConflict, IntermediateStore,
)
