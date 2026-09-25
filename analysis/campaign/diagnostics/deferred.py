"""Reserved detectors — registered so reports state they were NOT run.

Membrane / leaflet detectors are domain-specific (Level 4) and deferred.
"""
from __future__ import annotations

from analysis.campaign.diagnostics.registry import DetectorSpec, register

# ``molecule_split`` (reserved until Phase 13.5) is now implemented as
# ``molecule_periodic_image_change`` (diagnostics/molecules.py).
_RESERVED = [
    ("membrane_split_z", 4, "bilayer split across the periodic z boundary",
     "membrane-specific (Level 4) — deferred"),
    ("leaflet_discontinuity", 4, "leaflet assignment discontinuities",
     "requires leaflet assignment (future intermediate) — deferred"),
    ("membrane_normal_rotation", 4, "membrane normal rotating away from z",
     "membrane-specific (Level 4) — deferred"),
]

for _id, _level, _desc, _why in _RESERVED:
    register(DetectorSpec(id=_id, version="0", level=_level, description=_desc, deferred=_why))
