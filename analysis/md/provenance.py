"""Provenance record for an MD analysis run."""
from __future__ import annotations

import datetime
import sys
from dataclasses import dataclass, field
from pathlib import Path

from analysis.md.models import DiscoveredMDRun


@dataclass
class AnalysisArtifact:
    path: Path
    type: str

    def to_dict(self) -> dict:
        return {"path": str(self.path), "type": self.type}


@dataclass
class Provenance:
    timestamp: str
    simforge_version: str
    run_dir: str
    output_dir: str
    dry_run: bool
    argv: list[str]
    warnings: list[dict] = field(default_factory=list)
    discovered_files: dict[str, object] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "timestamp": self.timestamp,
            "simforge_version": self.simforge_version,
            "run_dir": self.run_dir,
            "output_dir": self.output_dir,
            "dry_run": self.dry_run,
            "argv": self.argv,
            "warnings": self.warnings,
            "discovered_files": self.discovered_files,
        }


def build_provenance(run: DiscoveredMDRun, output_dir: Path, dry_run: bool) -> Provenance:
    try:
        from simforge import __version__ as _v
        version = _v
    except Exception:
        version = "unknown"

    discovered: dict[str, object] = {
        "tpr":                  str(run.tpr)        if run.tpr        else None,
        "trajectory":           str(run.trajectory) if run.trajectory else None,
        "structure":            str(run.structure)  if run.structure  else None,
        "edr":                  str(run.edr)        if run.edr        else None,
        "topology":             str(run.topology)   if run.topology   else None,
        "index":                str(run.index)      if run.index      else None,
        "log":                  str(run.log)        if run.log        else None,
        "xvg_files":            [str(p) for p in run.xvg_files],
        "trajectory_candidates":[str(p) for p in run.trajectory_candidates],
    }

    return Provenance(
        timestamp=datetime.datetime.utcnow().isoformat() + "Z",
        simforge_version=version,
        run_dir=str(run.run_dir),
        output_dir=str(output_dir.resolve()),
        dry_run=dry_run,
        argv=sys.argv[:],
        warnings=[w.to_dict() for w in run.warnings],
        discovered_files=discovered,
    )
