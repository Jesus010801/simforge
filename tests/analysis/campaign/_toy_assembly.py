"""Toy periodic fixtures for Phase 13.5: a covalent dimer (two topological
molecules) + a 2-atom ligand, in a rectangular or triclinic cell.  Frames:
A (both subunits in one image), B (subunit 2 translated by an exact lattice
vector), A.  A and B are the same physical configuration."""
import subprocess, sys
from pathlib import Path
import numpy as np

TOP = """[ defaults ]
1 2 yes 0.5 0.8333
[ atomtypes ]
CX 12.011 0.0 A 0.34 0.36
[ moleculetype ]
SUB 3
[ atoms ]
1 CX 1 SUB C1 1 0 12.011
2 CX 1 SUB C2 1 0 12.011
3 CX 1 SUB C3 1 0 12.011
4 CX 1 SUB C4 1 0 12.011
[ bonds ]
1 2 1 0.15 1000
2 3 1 0.15 1000
3 4 1 0.15 1000
[ moleculetype ]
LIG 3
[ atoms ]
1 CX 1 LIG L1 1 0 12.011
2 CX 1 LIG L2 1 0 12.011
[ bonds ]
1 2 1 0.15 1000
[ system ]
toy
[ molecules ]
SUB 2
LIG 1
"""
MDP = "integrator=md\nnsteps=0\ncutoff-scheme=Verlet\nrvdw=1.0\nrcoulomb=1.0\ncoulombtype=cut-off\n"


def gro(frames, box9, path):
    names = [("SUB", "C1"), ("SUB", "C2"), ("SUB", "C3"), ("SUB", "C4")] * 2 + [("LIG", "L1"), ("LIG", "L2")]
    res = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3]
    out = []
    for t, x in frames:
        out.append(f"toy t= {t:.5f} step= {int(t)}")
        out.append(f"{len(x):5d}")
        for i, ((rn, an), r, p) in enumerate(zip(names, res, x), 1):
            out.append(f"{r:5d}{rn:<5s}{an:>5s}{i:5d}{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}")
        out.append(" ".join(f"{v:.5f}" for v in box9))
    Path(path).write_text("\n".join(out) + "\n")


def build(d, triclinic):
    d = Path(d); d.mkdir(parents=True, exist_ok=True)
    if triclinic:
        a, b, c = np.array([5.0, 0, 0]), np.array([2.0, 5.0, 0]), np.array([1.0, 1.5, 5.0])
        box9 = [5, 5, 5, 0, 0, 2.0, 0, 1.0, 1.5]
    else:
        a, b, c = np.array([5.0, 0, 0]), np.array([0, 5.0, 0]), np.array([0, 0, 5.0])
        box9 = [5, 5, 5, 0, 0, 0, 0, 0, 0]
    s1 = np.array([[2.0, 2.5, 2.5], [2.15, 2.5, 2.5], [2.3, 2.5, 2.5], [2.45, 2.5, 2.5]])
    s2 = s1 + np.array([0, 0.45, 0.0])                      # adjacent: a compact dimer
    lig = np.array([[2.2, 3.25, 2.6], [2.35, 3.25, 2.6]])
    A = np.vstack([s1, s2, lig])
    B = np.vstack([s1, s2 + b, lig])                         # subunit 2 moved by lattice vector b
    gro([(0.0, A)], box9, d / "md.gro")
    gro([(0.0, A), (10.0, B), (20.0, A)], box9, d / "frames.gro")
    (d / "topol.top").write_text(TOP); (d / "md.mdp").write_text(MDP)
    run = lambda *a, **k: subprocess.run(list(a), cwd=d, capture_output=True, text=True, **k)
    r = run("gmx", "grompp", "-f", "md.mdp", "-c", "md.gro", "-p", "topol.top", "-o", "md.tpr", "-maxwarn", "10")
    assert r.returncode == 0, r.stderr[-800:]
    r = run("gmx", "trjconv", "-f", "frames.gro", "-s", "md.tpr", "-o", "md.xtc", input="0\n")
    assert r.returncode == 0, r.stderr[-800:]
    ndx = ("[ Receptor ]\n1 2 3 4 5 6 7 8\n[ Ligand ]\n9 10\n[ Complex ]\n1 2 3 4 5 6 7 8 9 10\n"
           "[ System ]\n1 2 3 4 5 6 7 8 9 10\n[ Receptor_Backbone ]\n1 2 3 4 5 6 7 8\n"
           "[ Ann_span ]\n3 4 5 6\n[ Ann_one ]\n1 2 3\n")
    (d / "index.ndx").write_text(ndx)
    return d

