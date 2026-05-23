from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

from runtime.artifacts import ArtifactRef, checksum, make_ref


# ─── Base ────────────────────────────────────────────────────────────────────

@dataclass
class BaseArtifact:
    path:    Path
    step_id: str

    @property
    def exists(self) -> bool:
        return self.path.exists()

    def as_ref(self) -> ArtifactRef:
        raise NotImplementedError


# ─── Topology ────────────────────────────────────────────────────────────────

@dataclass
class TopologyState(BaseArtifact):
    """
    Tracks a GROMACS topology file (topol.top) through its mutation chain.

    The topology is created by pdb2gmx and then modified in-place by:
      solvate_membrane → clean_water → add_ions

    Each modification is recorded so callers can reason about who last
    changed the file.
    """
    forcefield:  str              = ""
    n_molecules: dict[str, int]   = field(default_factory=dict)
    mutations:   list[str]        = field(default_factory=list)  # step_ids in order

    def record_mutation(self, step_id: str) -> None:
        self.mutations.append(step_id)

    def as_ref(self) -> ArtifactRef:
        return make_ref(self.path, "topology", self.step_id)


# ─── Coordinates ─────────────────────────────────────────────────────────────

@dataclass
class CoordinateArtifact(BaseArtifact):
    """GRO/PDB file produced by a simulation or preparation step."""
    n_atoms:    int    = 0
    box_nm:     tuple[float, float, float] | None = None
    source_step: str  = ""  # which step produced this GRO

    def as_ref(self) -> ArtifactRef:
        return make_ref(self.path, "coordinates", self.step_id)


# ─── Trajectory ──────────────────────────────────────────────────────────────

@dataclass
class TrajectoryArtifact(BaseArtifact):
    """XTC or TRR trajectory file."""
    format:      Literal["xtc", "trr"] = "xtc"
    n_frames:    int   = 0
    duration_ns: float = 0.0
    dt_ps:       float = 0.0

    def as_ref(self) -> ArtifactRef:
        return make_ref(self.path, "trajectory", self.step_id)

    def size_gb(self) -> float:
        try:
            return self.path.stat().st_size / 1e9
        except FileNotFoundError:
            return 0.0


# ─── Analysis ─────────────────────────────────────────────────────────────────

@dataclass
class AnalysisArtifact(BaseArtifact):
    """XVG or other analysis output."""
    analysis_type: str = ""   # "rmsd", "hbond", "distance", etc.
    units:         str = ""

    def as_ref(self) -> ArtifactRef:
        return make_ref(self.path, "analysis", self.step_id)


# ─── Checkpoint ──────────────────────────────────────────────────────────────

@dataclass
class CheckpointArtifact(BaseArtifact):
    """GROMACS .cpt checkpoint file."""
    simulation_time_ps: float = 0.0

    def as_ref(self) -> ArtifactRef:
        return make_ref(self.path, "checkpoint", self.step_id)
