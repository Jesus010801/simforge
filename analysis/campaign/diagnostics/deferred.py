"""Reserved detectors — registered so reports state they were NOT run.

``molecule_split`` needs per-molecule connectivity (bonds / molecule blocks
from the .tpr) plus sampled per-atom coordinates; a centre-of-mass heuristic
cannot establish splitting and is deliberately not shipped.
"""
from __future__ import annotations

from analysis.campaign.diagnostics.registry import DetectorSpec, register

_RESERVED = [
    ("molecule_split", 3,
     "molecules split across periodic boundaries",
     "needs .tpr molecule/bond connectivity and sampled per-atom coordinates; "
     "not derivable robustly from centre series — deferred"),
    ("membrane_split_z", 4, "bilayer split across the periodic z boundary",
     "membrane-specific (Level 4) — deferred"),
    ("leaflet_discontinuity", 4, "leaflet assignment discontinuities",
     "requires leaflet assignment (future intermediate) — deferred"),
    ("membrane_normal_rotation", 4, "membrane normal rotating away from z",
     "membrane-specific (Level 4) — deferred"),
]

for _id, _level, _desc, _why in _RESERVED:
    register(DetectorSpec(id=_id, version="0", level=_level, description=_desc, deferred=_why))
