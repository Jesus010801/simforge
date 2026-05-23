"""
core/md_knowledge — Scientific knowledge base for MD simulation interpretation.

Provides context-aware, evidence-based interpretation of GROMACS trajectories.
Philosophy: soft ranges + temporal patterns + multi-observable evidence, NOT rigid thresholds.
"""
from core.md_knowledge.states import SimulationState, STATE_DESCRIPTIONS
from core.md_knowledge.patterns import TemporalPattern, PatternResult, detect_temporal_pattern
from core.md_knowledge.contexts import SystemContext, SYSTEM_CONTEXTS
from core.md_knowledge.heuristics import OBSERVABLE_HEURISTICS, get_heuristic
from core.md_knowledge.evidence import Evidence, EvidenceBundle, accumulate_evidence
from core.md_knowledge.interpreter import InterpretationResult, interpret_simulation

__all__ = [
    "SimulationState", "STATE_DESCRIPTIONS",
    "TemporalPattern", "PatternResult", "detect_temporal_pattern",
    "SystemContext", "SYSTEM_CONTEXTS",
    "OBSERVABLE_HEURISTICS", "get_heuristic",
    "Evidence", "EvidenceBundle", "accumulate_evidence",
    "InterpretationResult", "interpret_simulation",
]
