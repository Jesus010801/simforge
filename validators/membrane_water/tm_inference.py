"""
validators/membrane_water/tm_inference.py

Best-effort transmembrane (TM) residue inference from membrane geometry when
the caller does not supply an explicit `tm_residues` set. This is new
capability: nothing else in the repo infers TM residues geometrically today
(everything else requires an explicit user-supplied list — see
core/structural_annotation.py). It is deliberately conservative and reports
a confidence score rather than a hard yes/no, so a low-confidence inference
can force safe cleanup behavior downstream instead of silently degrading to
"delete everything in the membrane Z-range."

Heuristic (documented, not claimed to be exact): a protein residue is a TM
candidate if its heavy-atom centroid Z falls inside the *local* membrane
core band (from leaflets.MembraneModel, not a global scalar). Candidates are
grouped into contiguous residue-number segments (small gaps tolerated).
Confidence combines segment-length plausibility (real TM helices/strands are
typically 15-35 residues), segment count, and overall candidate fraction of
the protein (a plausible TM annotation is a minority of the protein, not all
of it — flagging 100% of residues is itself a red flag, not a confident
answer).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from validators.membrane_water.leaflets import MembraneModel
from validators.membrane_water.species import SpeciesData

_MAX_SEGMENT_GAP = 2
_PLAUSIBLE_SEGMENT_LEN = (12, 40)
_MARGIN_NM = 0.2


@dataclass
class TMAnnotation:
    source: str                       # "explicit" | "inferred"
    residue_count: int
    segments: list[tuple[int, int]]
    confidence: float
    residues: "set[int]" = field(default_factory=set)


def explicit_annotation(tm_residues: "set[int]") -> TMAnnotation:
    residues = set(tm_residues)
    segments = _segments_from_residues(sorted(residues))
    return TMAnnotation(
        source="explicit",
        residue_count=len(residues),
        segments=segments,
        confidence=1.0,
        residues=residues,
    )


def infer_tm_annotation(species: SpeciesData, membrane_model: MembraneModel) -> TMAnnotation:
    per_residue: dict[int, list[float]] = {}
    per_residue_xy: dict[int, list[tuple[float, float]]] = {}
    for i in species.protein_idx.tolist():
        atom = species.gro.atoms[i]
        per_residue.setdefault(atom.residue_number, []).append(atom.z)
        per_residue_xy.setdefault(atom.residue_number, []).append((atom.x, atom.y))

    if not per_residue:
        return TMAnnotation(source="inferred", residue_count=0, segments=[], confidence=0.0, residues=set())

    resids = sorted(per_residue.keys())
    centroid_xy = np.array([
        (float(np.mean([p[0] for p in per_residue_xy[r]])), float(np.mean([p[1] for p in per_residue_xy[r]])))
        for r in resids
    ])
    centroid_z = np.array([float(np.mean(per_residue[r])) for r in resids])

    lower, upper = membrane_model.core_bounds_batch(centroid_xy)
    in_band = (centroid_z >= lower - _MARGIN_NM) & (centroid_z <= upper + _MARGIN_NM)

    candidates = sorted(int(r) for r, ok in zip(resids, in_band) if ok)
    segments = _segments_from_residues(candidates)

    total_protein_residues = len(resids)
    candidate_fraction = len(candidates) / total_protein_residues if total_protein_residues else 0.0

    confidence = _score_confidence(segments, candidate_fraction)

    return TMAnnotation(
        source="inferred",
        residue_count=len(candidates),
        segments=segments,
        confidence=round(confidence, 4),
        residues=set(candidates),
    )


def _segments_from_residues(residues: list[int]) -> list[tuple[int, int]]:
    if not residues:
        return []
    segments: list[tuple[int, int]] = []
    start = prev = residues[0]
    for r in residues[1:]:
        if r - prev <= 1 + _MAX_SEGMENT_GAP:
            prev = r
            continue
        segments.append((start, prev))
        start = prev = r
    segments.append((start, prev))
    return segments


def _score_confidence(segments: list[tuple[int, int]], candidate_fraction: float) -> float:
    if not segments:
        return 0.0

    lo, hi = _PLAUSIBLE_SEGMENT_LEN
    length_scores = []
    for start, end in segments:
        length = end - start + 1
        if lo <= length <= hi:
            length_scores.append(1.0)
        elif length < lo:
            length_scores.append(max(0.0, length / lo))
        else:
            length_scores.append(max(0.0, hi / length))
    length_score = float(np.mean(length_scores))

    segment_count_score = min(1.0, len(segments) / 4.0)  # saturates around 4+ TM segments

    if 0.05 <= candidate_fraction <= 0.65:
        fraction_score = 1.0
    elif candidate_fraction < 0.05:
        fraction_score = candidate_fraction / 0.05
    else:
        fraction_score = max(0.0, 1.0 - (candidate_fraction - 0.65) / 0.35)

    return 0.5 * length_score + 0.25 * segment_count_score + 0.25 * fraction_score
