"""Study-layer data models for comparative multi-system MD analysis."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class ObservableSeries:
    """Time-series for one observable measured in one replica."""
    observable:        str
    replica:           str
    system:            str
    xvg_path:          Path
    time_ns:           list[float]
    values:            list[float]
    mean:              float = 0.0
    std:               float = 0.0
    median:            float = 0.0
    drift:             float = 0.0          # units/ns, linear fit over full trajectory
    plateau_std:       float = 0.0          # std of last 20% window
    convergence_score: float = 0.0          # 0 = not converged, 1 = well converged


@dataclass
class Replica:
    """One replica within a system."""
    label:           str
    system:          str
    observables:     dict[str, ObservableSeries] = field(default_factory=dict)
    is_outlier:      bool = False
    outlier_reasons: list[str] = field(default_factory=list)


@dataclass
class AggregateMetrics:
    """Cross-replica aggregate statistics for one observable within one system."""
    observable:       str
    system:           str
    n_replicas:       int
    mean:             float     # mean of per-replica means
    std:              float     # inter-replica std (variability across replicas)
    median:           float
    min_val:          float
    max_val:          float
    mean_convergence: float     # mean convergence score across replicas
    units:            str = ""


@dataclass
class SystemGroup:
    """One experimental system (e.g., 'AA') with one or more replicas."""
    name:      str
    replicas:  dict[str, Replica]       = field(default_factory=dict)
    aggregate: dict[str, AggregateMetrics] = field(default_factory=dict)

    @property
    def n_replicas(self) -> int:
        return len(self.replicas)

    @property
    def observables(self) -> list[str]:
        obs: set[str] = set()
        for r in self.replicas.values():
            obs.update(r.observables.keys())
        return sorted(obs)


@dataclass
class ComparativeFinding:
    message:    str
    level:      str           = "info"      # info | highlight | warning
    system:     Optional[str] = None
    observable: Optional[str] = None
    replica:    Optional[str] = None


@dataclass
class ComparativeSummary:
    findings:         list[ComparativeFinding]          = field(default_factory=list)
    outlier_replicas: list[tuple[str, str, str]]        = field(default_factory=list)  # (system, replica, reason)
    system_ranking:   dict[str, float]                  = field(default_factory=dict)  # system → stability score


@dataclass
class Study:
    """Top-level study object representing a full comparative MD dataset."""
    path:                 Path
    systems:              dict[str, SystemGroup]  = field(default_factory=dict)
    observables_detected: list[str]               = field(default_factory=list)
    observable_display:   dict[str, str]          = field(default_factory=dict)
    observable_units:     dict[str, str]          = field(default_factory=dict)
    summary:              Optional[ComparativeSummary] = None
    n_xvg_discovered:     int = 0
    n_xvg_parsed:         int = 0
    n_xvg_ungrouped:      int = 0
    parse_errors:         list[str] = field(default_factory=list)
