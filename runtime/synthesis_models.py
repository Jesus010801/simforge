"""Data models for the scientific synthesis layer."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class SignalProfile:
    """Normalized observable signals for one system."""
    system: str
    # normalized: 0 = lowest across study, 1 = highest (sigmoid-based)
    normalized:    dict[str, float] = field(default_factory=dict)
    # direction_score: 0 = unfavorable, 1 = favorable (accounts for polarity)
    direction_score: dict[str, float] = field(default_factory=dict)
    # consistency: 0 = high inter-replica variance, 1 = fully consistent
    consistency:   dict[str, float] = field(default_factory=dict)


@dataclass
class RuleMatch:
    """Result of evaluating one interaction rule for one system."""
    state:                  str
    description:            str
    score:                  float        # raw rule score 0-1
    confidence:             float        # score × consensus_multiplier
    n_conditions_met:       int
    n_conditions_available: int
    supporting:             list[str]    # observables supporting this state
    opposing:               list[str]    # observables contradicting this state


@dataclass
class SystemSynthesis:
    """Complete scientific synthesis for one system."""
    system:            str
    primary_state:     str
    primary_confidence: float
    binding_score:     float        # 0-1 composite binding quality
    stability_score:   float        # 0-1 structural stability
    composite_score:   float        # overall 0-1
    rule_matches:      list[RuleMatch]
    evidence:          list[str]    # key findings
    explanation:       str          # one-paragraph narrative


@dataclass
class ConsensusResult:
    """Replica consensus for one system × observable."""
    system:      str
    observable:  str
    n_total:     int
    n_agreeing:  int
    label:       str    # "full" | "strong" | "moderate" | "weak" | "conflicting"
    multiplier:  float  # confidence multiplier 0.5-1.0


@dataclass
class TemporalEvent:
    """A detected time-resolved event in a simulation trajectory."""
    system:     str
    replica:    str
    observable: str
    event_type: str     # "late_destabilization" | "abrupt_transition" | "contact_loss" | "ligand_drift"
    time_ns:    float   # approximate onset time
    description: str


@dataclass
class SynthesisResult:
    """Complete scientific synthesis for a comparative MD study."""
    systems:   dict[str, SystemSynthesis] = field(default_factory=dict)
    consensus: list[ConsensusResult]      = field(default_factory=list)
    events:    list[TemporalEvent]        = field(default_factory=list)
    ranking:   list[tuple[str, float]]   = field(default_factory=list)  # (system, score)
    narrative: str = ""
