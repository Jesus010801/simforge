"""Unit tests for utils/gro_parser.py."""

from pathlib import Path

import pytest

from utils.gro_parser import GroAtom, GroFile, parse_gro, write_gro


# ── fixtures ──────────────────────────────────────────────────────────────────

def _gro_atom(resnum, resname, atomname, atomnum, x, y, z, vx=None, vy=None, vz=None):
    """Build a correctly-formatted .gro atom line (fixed-width)."""
    line = f"{resnum:5d}{resname:<5s}{atomname:>5s}{atomnum:5d}{x:8.3f}{y:8.3f}{z:8.3f}"
    if vx is not None:
        line += f"{vx:8.4f}{vy:8.4f}{vz:8.4f}"
    return line


MINIMAL_GRO = "\n".join([
    "Simple test system",
    "    3",
    _gro_atom(1, "LIG", "C1", 1, 1.0, 2.0, 3.0),
    _gro_atom(1, "LIG", "N2", 2, 4.0, 5.0, 6.0),
    _gro_atom(1, "LIG", "O3", 3, 7.0, 8.0, 9.0),
    "10.00000   10.00000   10.00000",
]) + "\n"

GRO_WITH_VELOCITIES = "\n".join([
    "Velocities test",
    "    2",
    _gro_atom(1, "SOL", "OW",  1, 0.100, 0.200, 0.300, 0.1000, 0.2000, 0.3000),
    _gro_atom(1, "SOL", "HW1", 2, 0.110, 0.210, 0.310, 0.1100, 0.2100, 0.3100),
    "5.00000   5.00000   5.00000",
]) + "\n"


@pytest.fixture
def minimal_gro_file(tmp_path: Path) -> Path:
    p = tmp_path / "minimal.gro"
    p.write_text(MINIMAL_GRO)
    return p


@pytest.fixture
def velocities_gro_file(tmp_path: Path) -> Path:
    p = tmp_path / "velocities.gro"
    p.write_text(GRO_WITH_VELOCITIES)
    return p


# ── parse ─────────────────────────────────────────────────────────────────────

def test_parse_title(minimal_gro_file):
    gro = parse_gro(minimal_gro_file)
    assert gro.title == "Simple test system"


def test_parse_atom_count(minimal_gro_file):
    gro = parse_gro(minimal_gro_file)
    assert gro.atom_count == 3


def test_parse_atom_names(minimal_gro_file):
    gro = parse_gro(minimal_gro_file)
    names = [a.atom_name for a in gro.atoms]
    assert names == ["C1", "N2", "O3"]


def test_parse_residue_name(minimal_gro_file):
    gro = parse_gro(minimal_gro_file)
    assert all(a.residue_name == "LIG" for a in gro.atoms)


def test_parse_coordinates(minimal_gro_file):
    gro = parse_gro(minimal_gro_file)
    atom = gro.atoms[0]
    assert atom.x == pytest.approx(1.0)
    assert atom.y == pytest.approx(2.0)
    assert atom.z == pytest.approx(3.0)


def test_parse_box(minimal_gro_file):
    gro = parse_gro(minimal_gro_file)
    assert gro.box == pytest.approx([10.0, 10.0, 10.0])


def test_parse_no_velocities(minimal_gro_file):
    gro = parse_gro(minimal_gro_file)
    assert gro.atoms[0].vx is None


def test_parse_velocities(velocities_gro_file):
    gro = parse_gro(velocities_gro_file)
    atom = gro.atoms[0]
    assert atom.vx == pytest.approx(0.1)
    assert atom.vy == pytest.approx(0.2)
    assert atom.vz == pytest.approx(0.3)


def test_parse_atom_numbers(minimal_gro_file):
    gro = parse_gro(minimal_gro_file)
    assert [a.atom_number for a in gro.atoms] == [1, 2, 3]


# ── error handling ────────────────────────────────────────────────────────────

def test_parse_too_few_lines(tmp_path):
    p = tmp_path / "bad.gro"
    p.write_text("title\n")
    with pytest.raises(ValueError, match="Too few lines"):
        parse_gro(p)


def test_parse_bad_atom_count(tmp_path):
    p = tmp_path / "bad.gro"
    # Three lines so we pass the "too few lines" check, but line 2 is not a number
    p.write_text("title\nnotanumber\n1.0 1.0 1.0\n")
    with pytest.raises(ValueError, match="atom count"):
        parse_gro(p)


def test_parse_count_mismatch(tmp_path):
    p = tmp_path / "bad.gro"
    content = "title\n    5\n    1LIG      C1    1   1.000   2.000   3.000\n10.0 10.0 10.0\n"
    p.write_text(content)
    with pytest.raises(ValueError, match="Declared 5 atoms"):
        parse_gro(p)


# ── write / round-trip ────────────────────────────────────────────────────────

def test_write_round_trip(minimal_gro_file, tmp_path):
    gro = parse_gro(minimal_gro_file)
    out = tmp_path / "out.gro"
    write_gro(gro, out)
    gro2 = parse_gro(out)

    assert gro2.title == gro.title
    assert gro2.atom_count == gro.atom_count
    assert gro2.box == pytest.approx(gro.box)
    for a, b in zip(gro.atoms, gro2.atoms):
        assert a.atom_name == b.atom_name
        assert a.x == pytest.approx(b.x, abs=1e-3)
        assert a.y == pytest.approx(b.y, abs=1e-3)
        assert a.z == pytest.approx(b.z, abs=1e-3)


def test_write_velocities_round_trip(velocities_gro_file, tmp_path):
    gro = parse_gro(velocities_gro_file)
    out = tmp_path / "out_vel.gro"
    write_gro(gro, out)
    gro2 = parse_gro(out)
    assert gro2.atoms[0].vx == pytest.approx(gro.atoms[0].vx, abs=1e-4)


def test_write_preserves_atom_count(minimal_gro_file, tmp_path):
    gro = parse_gro(minimal_gro_file)
    out = tmp_path / "out.gro"
    write_gro(gro, out)
    lines = out.read_text().splitlines()
    assert int(lines[1].strip()) == 3


def test_write_preserves_box(minimal_gro_file, tmp_path):
    gro = parse_gro(minimal_gro_file)
    out = tmp_path / "out.gro"
    write_gro(gro, out)
    gro2 = parse_gro(out)
    assert gro2.box == pytest.approx([10.0, 10.0, 10.0], abs=1e-4)


# ── GroFile helpers ───────────────────────────────────────────────────────────

def test_grofile_atom_count_property():
    gro = GroFile(title="t", atoms=[
        GroAtom(1, "RES", "CA", 1, 0.0, 0.0, 0.0),
        GroAtom(1, "RES", "CB", 2, 1.0, 0.0, 0.0),
    ], box=[10.0, 10.0, 10.0])
    assert gro.atom_count == 2
