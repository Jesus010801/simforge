"""
Kabsch alignment and RMSD utilities for comparing atom coordinate sets.

No RDKit dependency. Works on plain numpy arrays.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def compute_centroid(coords: NDArray[np.float64]) -> NDArray[np.float64]:
    """Return the centroid of an (N, 3) coordinate array."""
    return coords.mean(axis=0)


def kabsch_rotation(
    mobile: NDArray[np.float64],
    reference: NDArray[np.float64],
) -> NDArray[np.float64]:
    """
    Compute the Kabsch optimal rotation matrix that minimises RMSD.

    Both arrays must be pre-centred (centroid at origin).
    Returns a (3, 3) rotation matrix R such that mobile @ R ≈ reference.
    """
    H = mobile.T @ reference
    U, _, Vt = np.linalg.svd(H)
    # Correct for reflection: ensure det(R) == +1
    d = np.linalg.det(U @ Vt)
    sign_matrix = np.diag([1.0, 1.0, d])
    # Right-multiplication convention: aligned = mobile @ R
    R = U @ sign_matrix @ Vt
    return R


def align(
    mobile: NDArray[np.float64],
    reference: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64], float]:
    """
    Superimpose *mobile* onto *reference* using the Kabsch algorithm.

    Returns:
        aligned  — transformed mobile coordinates (N, 3)
        R        — rotation matrix applied (3, 3)
        rmsd     — RMSD after alignment
    """
    if mobile.shape != reference.shape:
        raise ValueError(
            f"Shape mismatch: mobile {mobile.shape} vs reference {reference.shape}"
        )

    c_mob = compute_centroid(mobile)
    c_ref = compute_centroid(reference)

    mob_centered = mobile - c_mob
    ref_centered = reference - c_ref

    R = kabsch_rotation(mob_centered, ref_centered)
    aligned = mob_centered @ R + c_ref
    rmsd = compute_rmsd(aligned, reference)
    return aligned, R, rmsd


def compute_rmsd(
    coords_a: NDArray[np.float64],
    coords_b: NDArray[np.float64],
) -> float:
    """Return RMSD between two (N, 3) coordinate arrays (no alignment)."""
    if coords_a.shape != coords_b.shape:
        raise ValueError(
            f"Shape mismatch: {coords_a.shape} vs {coords_b.shape}"
        )
    diff = coords_a - coords_b
    return float(np.sqrt((diff ** 2).sum() / len(coords_a)))


def elements_match(elements_a: list[str], elements_b: list[str]) -> bool:
    """
    Return True if both ordered element lists are identical (case-insensitive).

    Useful for verifying that two coordinate sets correspond to the same atoms
    before computing RMSD.
    """
    if len(elements_a) != len(elements_b):
        return False
    return all(a.upper() == b.upper() for a, b in zip(elements_a, elements_b))
