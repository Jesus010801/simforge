"""Unit tests for utils/geometry_alignment.py."""

import numpy as np
import pytest

from utils.geometry_alignment import (
    align,
    compute_centroid,
    compute_rmsd,
    elements_match,
    kabsch_rotation,
)


def _random_coords(n: int, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.standard_normal((n, 3))


# ── centroid ─────────────────────────────────────────────────────────────────

def test_centroid_origin():
    coords = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    c = compute_centroid(coords)
    np.testing.assert_allclose(c, [0.0, 0.0, 0.0], atol=1e-12)


def test_centroid_shifted():
    coords = np.ones((4, 3)) * 3.0
    np.testing.assert_allclose(compute_centroid(coords), [3.0, 3.0, 3.0])


# ── RMSD ─────────────────────────────────────────────────────────────────────

def test_rmsd_identical():
    coords = _random_coords(10)
    assert compute_rmsd(coords, coords) == pytest.approx(0.0, abs=1e-12)


def test_rmsd_known_value():
    a = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    b = np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    # every atom displaced by 1 Å along x → RMSD = 1.0
    assert compute_rmsd(a, b) == pytest.approx(1.0, abs=1e-9)


def test_rmsd_shape_mismatch():
    with pytest.raises(ValueError, match="Shape mismatch"):
        compute_rmsd(np.zeros((3, 3)), np.zeros((4, 3)))


# ── Kabsch rotation ───────────────────────────────────────────────────────────

def test_kabsch_identity():
    """Kabsch of identical centred coords should return identity matrix."""
    coords = _random_coords(8)
    centred = coords - compute_centroid(coords)
    R = kabsch_rotation(centred, centred)
    np.testing.assert_allclose(R, np.eye(3), atol=1e-10)


def test_kabsch_known_rotation():
    """45° rotation around Z recovered for non-coplanar 3D points."""
    rng = np.random.default_rng(42)
    mobile = rng.standard_normal((12, 3))
    mobile -= compute_centroid(mobile)

    angle = np.pi / 4
    R_true = np.array([
        [np.cos(angle), -np.sin(angle), 0.0],
        [np.sin(angle),  np.cos(angle), 0.0],
        [0.0,            0.0,           1.0],
    ])
    # right-multiply: reference = mobile @ R_true
    reference = mobile @ R_true
    R_est = kabsch_rotation(mobile, reference)
    np.testing.assert_allclose(R_est, R_true, atol=1e-9)


def test_kabsch_det_positive():
    """Rotation matrix must have det +1 (no reflection)."""
    mobile = _random_coords(12, seed=7)
    ref = _random_coords(12, seed=99)
    mob_c = mobile - compute_centroid(mobile)
    ref_c = ref - compute_centroid(ref)
    R = kabsch_rotation(mob_c, ref_c)
    assert np.linalg.det(R) == pytest.approx(1.0, abs=1e-9)


# ── align (full pipeline) ─────────────────────────────────────────────────────

def test_align_identical_coords():
    coords = _random_coords(20, seed=3)
    aligned, R, rmsd = align(coords, coords)
    assert rmsd == pytest.approx(0.0, abs=1e-9)


def test_align_pure_translation():
    coords = _random_coords(15, seed=5)
    translated = coords + np.array([3.0, -2.0, 1.0])
    aligned, R, rmsd = align(coords, translated)
    assert rmsd == pytest.approx(0.0, abs=1e-9)


def test_align_pure_rotation():
    """After Kabsch alignment of a rotated copy, RMSD must be ~0."""
    coords = _random_coords(10, seed=11)
    angle = np.pi / 4
    Rz = np.array([
        [np.cos(angle), -np.sin(angle), 0.0],
        [np.sin(angle),  np.cos(angle), 0.0],
        [0.0,            0.0,           1.0],
    ])
    rotated = coords @ Rz.T
    aligned, R, rmsd = align(coords, rotated)
    assert rmsd == pytest.approx(0.0, abs=1e-9)


def test_align_shape_mismatch():
    with pytest.raises(ValueError, match="Shape mismatch"):
        align(np.zeros((5, 3)), np.zeros((6, 3)))


def test_align_rmsd_positive_for_different():
    a = _random_coords(10, seed=1)
    b = _random_coords(10, seed=2)
    _, _, rmsd = align(a, b)
    assert rmsd > 0.0


# ── elements_match ────────────────────────────────────────────────────────────

def test_elements_match_same():
    assert elements_match(["C", "N", "O"], ["C", "N", "O"])


def test_elements_match_case_insensitive():
    assert elements_match(["c", "n"], ["C", "N"])


def test_elements_match_different():
    assert not elements_match(["C", "N"], ["C", "O"])


def test_elements_match_different_length():
    assert not elements_match(["C", "N", "O"], ["C", "N"])


def test_elements_match_empty():
    assert elements_match([], [])
