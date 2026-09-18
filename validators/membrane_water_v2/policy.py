"""Swappable cleanup policy, independent of atlas generation."""
from .models import CleanupDecision
REMOVE = {"membrane_defect_candidate", "steric_clash"}


def decide(mapping, policy="conservative"):
    if mapping.hard_clash:
        return CleanupDecision(mapping.molecule_id, "REMOVE", "confirmed hard excluded-volume clash", mapping.region_id, policy)
    if policy == "conservative":
        remove = mapping.classification == "membrane_defect_candidate" and not mapping.uncertainty
    elif policy == "aggressive":
        remove = mapping.classification in REMOVE and not mapping.uncertainty
    else:
        raise ValueError(f"Unknown water cleanup policy: {policy}")
    return CleanupDecision(mapping.molecule_id, "REMOVE" if remove else "PRESERVE",
        "supported membrane defect" if remove else ("uncertain or permitted solvent region" if mapping.uncertainty else "policy preserves this region"),
        mapping.region_id, policy)

