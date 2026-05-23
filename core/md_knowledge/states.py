"""core/md_knowledge/states.py — Simulation state ontology."""
from __future__ import annotations
from enum import Enum


class SimulationState(str, Enum):
    STABLE_EQUILIBRATED     = "stable_equilibrated"
    METASTABLE              = "metastable"
    DRIFTING                = "drifting"
    UNSTABLE                = "unstable"
    PARTIALLY_CONVERGED     = "partially_converged"
    LIGAND_DISSOCIATED      = "ligand_dissociated"
    INSUFFICIENT_SAMPLING   = "insufficient_sampling"
    NONPHYSICAL_BEHAVIOR    = "nonphysical_behavior"
    CONFORMATIONAL_TRANSITION = "conformational_transition"


STATE_DESCRIPTIONS: dict[SimulationState, str] = {
    SimulationState.STABLE_EQUILIBRATED: (
        "System has reached thermodynamic equilibrium. Key observables plateau with low "
        "variance. Structure is stable and dynamics are physically meaningful."
    ),
    SimulationState.METASTABLE: (
        "System occupies a local free-energy minimum but may not have sampled the global one. "
        "Observables are locally stable but history-dependent. "
        "Consider extended sampling or replica exchange."
    ),
    SimulationState.DRIFTING: (
        "One or more observables show monotonic drift without plateau. "
        "Either equilibration is incomplete or a slow conformational change is occurring. "
        "Extend simulation before drawing thermodynamic conclusions."
    ),
    SimulationState.UNSTABLE: (
        "Large fluctuations, energy spikes, or physically implausible values detected. "
        "Likely force-field / topology issue or insufficient minimization/equilibration."
    ),
    SimulationState.PARTIALLY_CONVERGED: (
        "Some observables plateau while others drift. Mixed convergence; "
        "statistical averages are unreliable for the drifting quantities."
    ),
    SimulationState.LIGAND_DISSOCIATED: (
        "Ligand has left the binding pocket during simulation. "
        "May indicate a weakly-bound complex, incorrect protonation state, "
        "or insufficient restraints during equilibration."
    ),
    SimulationState.INSUFFICIENT_SAMPLING: (
        "Simulation is too short to assess convergence reliably. "
        "Standard deviations are large relative to system complexity."
    ),
    SimulationState.NONPHYSICAL_BEHAVIOR: (
        "Values outside physically plausible ranges detected "
        "(T >> 400 K, negative pressure extremes, exploding RMSD). "
        "Simulation should be discarded or restarted."
    ),
    SimulationState.CONFORMATIONAL_TRANSITION: (
        "Observable shows a clear, sustained jump to a new plateau. "
        "May represent a biologically relevant conformational change "
        "or force-field artefact — inspect trajectory visually."
    ),
}


# Severity ranking: higher = more problematic
STATE_SEVERITY: dict[SimulationState, int] = {
    SimulationState.NONPHYSICAL_BEHAVIOR:     5,
    SimulationState.UNSTABLE:                 4,
    SimulationState.LIGAND_DISSOCIATED:        3,
    SimulationState.DRIFTING:                  2,
    SimulationState.CONFORMATIONAL_TRANSITION: 2,
    SimulationState.PARTIALLY_CONVERGED:       1,
    SimulationState.INSUFFICIENT_SAMPLING:     1,
    SimulationState.METASTABLE:                1,
    SimulationState.DRIFTING:                  2,
    SimulationState.STABLE_EQUILIBRATED:       0,
}
