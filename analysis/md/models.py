"""Data models for MD run discovery and analysis planning."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class DiscoveredFile:
    path: Path
    role: str
    size_bytes: int = 0

    def to_dict(self) -> dict:
        return {"path": str(self.path), "role": self.role, "size_bytes": self.size_bytes}


@dataclass
class AnalysisWarning:
    code: str
    message: str
    severity: str = "warn"  # warn | info

    def to_dict(self) -> dict:
        return {"code": self.code, "message": self.message, "severity": self.severity}


@dataclass
class DiscoveredMDRun:
    run_dir: Path
    tpr: Optional[Path] = None
    trajectory: Optional[Path] = None        # primary: .xtc preferred over .trr
    trajectory_candidates: list[Path] = field(default_factory=list)
    structure: Optional[Path] = None         # .gro or .pdb
    edr: Optional[Path] = None
    topology: Optional[Path] = None          # .top
    index: Optional[Path] = None             # .ndx
    log: Optional[Path] = None               # primary log
    xvg_files: list[Path] = field(default_factory=list)
    all_files: list[DiscoveredFile] = field(default_factory=list)
    warnings: list[AnalysisWarning] = field(default_factory=list)
    user_overrides: dict[str, str] = field(default_factory=dict)  # role -> original user path string

    def to_dict(self) -> dict:
        return {
            "run_dir": str(self.run_dir),
            "tpr": str(self.tpr) if self.tpr else None,
            "trajectory": str(self.trajectory) if self.trajectory else None,
            "trajectory_candidates": [str(p) for p in self.trajectory_candidates],
            "structure": str(self.structure) if self.structure else None,
            "edr": str(self.edr) if self.edr else None,
            "topology": str(self.topology) if self.topology else None,
            "index": str(self.index) if self.index else None,
            "log": str(self.log) if self.log else None,
            "xvg_files": [str(p) for p in self.xvg_files],
            "all_files": [f.to_dict() for f in self.all_files],
            "warnings": [w.to_dict() for w in self.warnings],
            "user_overrides": self.user_overrides,
        }
