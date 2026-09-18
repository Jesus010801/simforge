"""One lattice convention: column cell vectors, Cartesian nm, periodic XYZ.

Rotated orthogonal cells are supported. Skew cells are represented faithfully
but rejected by the numerical backend rather than silently treated as cuboids.
"""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class PeriodicCell:
    matrix: np.ndarray

    def __post_init__(self):
        h = np.asarray(self.matrix, dtype=float)
        if h.shape != (3, 3) or not np.isfinite(h).all() or np.linalg.det(h) <= 0:
            raise ValueError("Cell must be a finite, right-handed, nonsingular 3x3 matrix")
        object.__setattr__(self, "matrix", h.copy())

    @classmethod
    def from_gro(cls, values):
        if len(values) == 3:
            return cls(np.diag(values))
        if len(values) == 9:
            a, b, c, ay, az, bx, bz, cx, cy = values
            return cls(np.array([[a, bx, cx], [ay, b, cy], [az, bz, c]]))
        raise ValueError("GRO box requires 3 or 9 values")

    @property
    def lengths(self):
        return np.linalg.norm(self.matrix, axis=0)

    @property
    def basis(self):
        return self.matrix / self.lengths

    def require_orthogonal(self):
        if not np.allclose(self.basis.T @ self.basis, np.eye(3), atol=1e-7):
            raise ValueError("Skew cell backend unavailable: preserve water; do not use slab cleanup")

    def fractional(self, xyz):
        return np.asarray(xyz) @ np.linalg.inv(self.matrix).T

    def cartesian(self, fractional):
        return np.asarray(fractional) @ self.matrix.T

    def wrap(self, xyz):
        return self.cartesian(self.fractional(xyz) % 1.0)

    def minimum_image(self, delta):
        self.require_orthogonal()
        f = self.fractional(delta)
        return self.cartesian(f - np.floor(f + 0.5))

    def local(self, xyz):
        self.require_orthogonal()
        return (self.fractional(xyz) % 1.0) * self.lengths


@dataclass(frozen=True)
class Lattice:
    cell: PeriodicCell
    shape: tuple[int, int, int]
    phase: tuple[float, float, float] = (0.5, 0.5, 0.5)

    @classmethod
    def create(cls, cell, spacing, phase=(0.5, 0.5, 0.5)):
        cell.require_orthogonal()
        if not np.isfinite(spacing) or spacing <= 0:
            raise ValueError("spacing must be positive")
        return cls(cell, tuple(np.maximum(2, np.ceil(cell.lengths / spacing).astype(int))), tuple(phase))

    @property
    def steps(self):
        return self.cell.lengths / np.array(self.shape)

    @property
    def error_bound(self):
        return float(np.linalg.norm(self.steps) / 2)

    def points(self, flat_indices):
        ijk = np.array(np.unravel_index(flat_indices, self.shape)).T
        return self.cell.cartesian(((ijk + self.phase) / self.shape) % 1)

    def indices(self, xyz):
        q = np.floor(self.cell.fractional(xyz) * self.shape - self.phase + 0.5).astype(int)
        return q % self.shape

    def flat_indices(self, xyz):
        return np.ravel_multi_index(self.indices(xyz).T, self.shape)
