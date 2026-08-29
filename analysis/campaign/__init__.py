"""SimForge comparative MD **study campaign** analysis layer.

A researcher points SimForge at a directory tree of finished MD runs; SimForge
resolves files -> systems -> conditions -> replicates -> molecular components ->
trajectory validity -> observable-specific preprocessing -> explicitly selected
analyses -> reproducible, provenance-rich results.

Public entry points:

    from analysis.campaign.orchestration.study_analyzer import run_inspect, run_analyze
    from analysis.campaign.manifest import build_manifest, load_manifest, write_manifest
    from analysis.campaign.observables import registry

Nothing is analysed by default — the caller must name analyses explicitly.
"""

__all__ = ["models"]
