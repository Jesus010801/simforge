# core/workflow_hints.py
"""
WorkflowHints — structured hints from semantic layer to workflow policy.

Populated by semantic_inference.apply_preset() and _inject_membrane_objectives().
Consumed by decision_engine._build_workflow_policy().

This is the bridge that closes the gap between "what the user wants scientifically"
and "what parameters the simulation should use". It eliminates the knowledge
duplication that previously existed between SIMULATION_PRESETS (text hints in
global_reasoning.notes) and individual pipeline classes (MembraneWorkflowOPLSAA).

Pipeline flow:
    YAML simulation_profile
        → semantic_inference.apply_preset()
        → state.workflow_hints (this model)
        → decision_engine._build_workflow_policy()
        → WorkflowPolicy (concrete values)
        → _build_*_params()
        → step.params
        → Builder → MDP / scripts
"""
from __future__ import annotations
from pydantic import BaseModel


class WorkflowHints(BaseModel):
    """
    Structured scientific constraints from the semantic layer.

    Each field is a boolean signal that the decision engine translates into
    concrete simulation parameters. Fields map 1:1 to SIMULATION_PRESETS hints.
    """

    # ── Membrane system ───────────────────────────────────────────────────────
    membrane_required:      bool = False  # system embeds a bilayer
    semiisotropic_coupling: bool = False  # NPT/production use semiisotropic P-coupling
    conservative_timestep:  bool = False  # dt=0.001 ps (OPLS-AA lipids, 0.002 is unstable)
    lipid_aware_restraints: bool = False  # POSRES on protein only, lipids unrestrained
    membrane_equilibration: bool = False  # extended NPT equilibration for bilayer relaxation

    # ── Sampling / timing ─────────────────────────────────────────────────────
    extended_equilibration: bool = False  # longer NVT/NPT (IDRs, disordered loops)
    extended_production:    bool = False  # longer production run (allosteric, slow motions)
    long_production:        bool = False  # alias for extended_production (idr_sampling preset)
    enhanced_sampling:      bool = False  # enables REST2 / metadynamics

    def any_membrane_hint(self) -> bool:
        return (
            self.membrane_required
            or self.semiisotropic_coupling
            or self.conservative_timestep
            or self.membrane_equilibration
        )
