import numpy as np
import pytest
from validators.membrane_water_v2.cell import PeriodicCell, Lattice
from validators.membrane_water_v2.identity import load_system, OTHER


def write_gro(path, atoms, box=(6., 6., 8.)):
    lines = ["atlas test", str(len(atoms))]
    for i, (r, rn, an, x, y, z) in enumerate(atoms):
        lines.append(f"{r % 100000:5d}{rn:<5}{an:>5}{(i+1)%100000:5d}{x:8.3f}{y:8.3f}{z:8.3f}")
    lines.append(" ".join(str(x) for x in box))
    path.write_text("\n".join(lines) + "\n")
    return path


def water(r, xyz, model="SOL", sites=3):
    x, y, z = xyz
    atoms = [(r, model, "OW", x, y, z), (r, model, "HW1", x+.01, y, z), (r, model, "HW2", x, y+.01, z)]
    return atoms + [(r, model, f"MW{i}", x, y, z) for i in range(sites-3)]


def test_identity_and_species(tmp_path):
    atoms = [(1, "EUG", "C1", 1, 1, 1)] + water(1, (2,2,2)) + water(1, (3,3,3), "TIP4", 4)
    s = load_system(write_gro(tmp_path / "a.gro", atoms))
    assert [w.uid for w in s.waters] == [0, 1]
    assert [len(w.atom_indices) for w in s.waters] == [3, 4]
    assert s.geometry.species.tolist() == [OTHER]
    assert all(w.valid for w in s.waters)


def test_cell_period_and_rotation():
    h = np.diag([6.13, 7.27, 8.19])
    c = PeriodicCell(h)
    g = Lattice.create(c, .2)
    np.testing.assert_allclose(g.steps * g.shape, c.lengths)
    p = np.array([[.01, 2., 3.]])
    np.testing.assert_array_equal(g.indices(p), g.indices(p + h[:, 0]))
    theta = .71
    q = np.array([[np.cos(theta), 0, np.sin(theta)], [0,1,0], [-np.sin(theta),0,np.cos(theta)]])
    cr = PeriodicCell(q @ h)
    np.testing.assert_allclose(cr.minimum_image((p + h[:,0]) @ q.T), c.minimum_image(p) @ q.T)


def test_skew_is_explicitly_unsupported():
    c = PeriodicCell.from_gro([6, 6, 8, 0, 0, 1, 0, 0, 0])
    with pytest.raises(ValueError, match="Skew"):
        Lattice.create(c, .2)
