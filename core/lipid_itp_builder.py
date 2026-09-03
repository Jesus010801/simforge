"""
core/lipid_itp_builder.py
─────────────────────────
Generate standalone GROMACS ITP files for lipid molecules from RTP residue entries.

The RTP format defines residue templates (atoms, bonds, impropers).
This module generates complete ITP files including angles, proper dihedrals,
and 1-4 pairs derived from the bond graph — without calling pdb2gmx.

Usage:
    from core.lipid_itp_builder import build_lipid_itp_from_rtp
    itp_text = build_lipid_itp_from_rtp(ff_dir, residue_name="DPP")
    Path("dpp.itp").write_text(itp_text)
"""

from __future__ import annotations

from pathlib import Path
from typing import NamedTuple


# ─── Mapping of known lipid residue names → static ITP filenames ──────────────
# The ITP is shipped inside oplsaa_membrane.ff/ alongside forcefield.itp.
LIPID_ITP_REGISTRY: dict[str, str] = {
    "DPP":  "dpp.itp",   # DPPC united-atom in OPLS-AA membrane ff
    "DPPC": "dpp.itp",   # common alias
}


class _RtpAtom(NamedTuple):
    name: str
    atype: str
    charge: float
    cgnr: int


class _RtpResidue(NamedTuple):
    name: str
    atoms: list
    bonds: list        # list of (name1, name2)
    impropers: list    # list of (name1, name2, name3, name4, dih_type)


def _parse_rtp_residue(rtp_path: Path, residue_name: str) -> _RtpResidue:
    """Parse a named residue entry from a GROMACS .rtp file."""
    text = rtp_path.read_text()
    lines = text.splitlines()

    res_start: int | None = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("[") and "]" in stripped:
            inner = stripped[stripped.index("[") + 1 : stripped.index("]")].strip()
            if inner.upper() == residue_name.upper():
                res_start = i
                break

    if res_start is None:
        raise ValueError(f"Residue {residue_name!r} not found in {rtp_path}")

    _SUBSECTIONS = {
        "atoms", "bonds", "angles", "dihedrals", "impropers",
        "exclusions", "pairs", "cmap",
    }

    atoms: list[_RtpAtom] = []
    bonds: list[tuple[str, str]] = []
    impropers: list[tuple[str, str, str, str, str]] = []
    current: str | None = None

    for line in lines[res_start + 1 :]:
        stripped = line.strip()
        if not stripped or stripped.startswith(";"):
            continue
        if stripped.startswith("[") and "]" in stripped:
            inner = stripped[stripped.index("[") + 1 : stripped.index("]")].strip().lower()
            if inner not in _SUBSECTIONS:
                break   # hit next residue header
            current = inner
            continue
        parts = stripped.split(";")[0].split()
        if not parts:
            continue
        if current == "atoms" and len(parts) >= 4:
            atoms.append(_RtpAtom(
                name=parts[0], atype=parts[1],
                charge=float(parts[2]), cgnr=int(parts[3]),
            ))
        elif current == "bonds" and len(parts) >= 2:
            bonds.append((parts[0], parts[1]))
        elif current == "impropers" and len(parts) >= 5:
            impropers.append((parts[0], parts[1], parts[2], parts[3], parts[4]))

    return _RtpResidue(name=residue_name, atoms=atoms, bonds=bonds, impropers=impropers)


def _parse_masses(ff_dir: Path) -> dict[str, float]:
    """Parse atomtype → mass from atomtypes.atp."""
    masses: dict[str, float] = {}
    atp = ff_dir / "atomtypes.atp"
    if not atp.exists():
        return masses
    for line in atp.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith(";"):
            continue
        parts = stripped.split(";")[0].split()
        if len(parts) >= 2:
            try:
                masses[parts[0]] = float(parts[1])
            except ValueError:
                pass
    return masses


def _build_adj(n: int, bonds: list[tuple[str, str]], names: list[str]) -> dict[int, list[int]]:
    """Build sorted adjacency list (1-indexed) from atom name pairs."""
    idx = {name: i + 1 for i, name in enumerate(names)}
    adj: dict[int, list[int]] = {i + 1: [] for i in range(n)}
    for a, b in bonds:
        if a not in idx or b not in idx:
            continue
        i, j = idx[a], idx[b]
        adj[i].append(j)
        adj[j].append(i)
    for k in adj:
        adj[k] = sorted(set(adj[k]))
    return adj


def _angles(adj: dict[int, list[int]]) -> list[tuple[int, int, int]]:
    """Generate all angles i-j-k (center at j) from the bond graph."""
    result: list[tuple[int, int, int]] = []
    for j in sorted(adj):
        nbrs = adj[j]
        for a_idx in range(len(nbrs)):
            for b_idx in range(a_idx + 1, len(nbrs)):
                result.append((nbrs[a_idx], j, nbrs[b_idx]))
    return sorted(result)


def _proper_dihedrals(adj: dict[int, list[int]]) -> list[tuple[int, int, int, int]]:
    """Generate all proper dihedrals i-j-k-l from the bond graph.

    Each unique dihedral is emitted once (no duplicates from forward/reverse traversal).
    """
    seen: set[tuple[int, int, int, int]] = set()
    result: list[tuple[int, int, int, int]] = []
    for j in sorted(adj):
        for k in sorted(adj[j]):
            if k <= j:
                continue  # process each bond j-k once (j < k)
            for i in sorted(adj[j]):
                if i == k:
                    continue
                for l in sorted(adj[k]):
                    if l == j or l == i:
                        continue
                    fwd = (i, j, k, l)
                    rev = (l, k, j, i)
                    key = min(fwd, rev)
                    if key not in seen:
                        seen.add(key)
                        result.append(fwd)
    return sorted(result)


def _pairs_14(dihedrals: list[tuple[int, int, int, int]]) -> list[tuple[int, int]]:
    """Generate sorted 1-4 pairs from proper dihedral endpoints."""
    pairs: set[tuple[int, int]] = set()
    for i, _, _, l in dihedrals:
        pairs.add((min(i, l), max(i, l)))
    return sorted(pairs)


def build_lipid_itp_from_rtp(
    ff_dir: Path,
    residue_name: str = "DPP",
    nrexcl: int = 3,
    proper_dih_func: int = 3,
    improper_dih_func: int = 1,
) -> str:
    """
    Build a complete GROMACS ITP for a single-residue lipid from its RTP template.

    Parameters
    ----------
    ff_dir          : forcefield directory containing aminoacids.rtp and atomtypes.atp
    residue_name    : residue name as it appears in the RTP file
    nrexcl          : nonbonded exclusion count (3 for OPLS-AA)
    proper_dih_func : GROMACS dihedral function type for proper dihedrals
                      (3 = Ryckaert-Bellemans, used by OPLS-AA)
    improper_dih_func : GROMACS dihedral function type for improper dihedrals
                        (1 = periodic improper)

    Returns
    -------
    str : complete ITP file text, ready to write to disk or include in a topology
    """
    rtp = _parse_rtp_residue(ff_dir / "aminoacids.rtp", residue_name)
    masses = _parse_masses(ff_dir)

    atom_names = [a.name for a in rtp.atoms]
    name_to_idx = {name: i + 1 for i, name in enumerate(atom_names)}
    adj = _build_adj(len(atom_names), rtp.bonds, atom_names)

    bond_list = sorted({
        (min(name_to_idx[a], name_to_idx[b]), max(name_to_idx[a], name_to_idx[b]))
        for a, b in rtp.bonds
        if a in name_to_idx and b in name_to_idx
    })
    angle_list   = _angles(adj)
    dihedral_list = _proper_dihedrals(adj)
    pair_list    = _pairs_14(dihedral_list)
    improper_list = [
        (name_to_idx[a1], name_to_idx[a2], name_to_idx[a3], name_to_idx[a4], dt)
        for a1, a2, a3, a4, dt in rtp.impropers
        if all(n in name_to_idx for n in (a1, a2, a3, a4))
    ]

    L = [
        f"; {residue_name} lipid moleculetype — OPLS-AA membrane forcefield",
        f"; Generated from aminoacids.rtp by core/lipid_itp_builder.py",
        f"; {len(rtp.atoms)} atoms  {len(bond_list)} bonds  {len(angle_list)} angles"
        f"  {len(dihedral_list)} proper dihedrals  {len(pair_list)} pairs"
        f"  {len(improper_list)} impropers",
        f"; Forcefield atom types and bonded parameters are looked up via forcefield.itp",
        "",
        "[ moleculetype ]",
        f"; Name    nrexcl",
        f"{residue_name:<8}      {nrexcl}",
        "",
        "[ atoms ]",
        ";   nr       type  resnr residue  atom   cgnr     charge       mass",
    ]
    for nr, atom in enumerate(rtp.atoms, start=1):
        mass = masses.get(atom.atype, 0.0)
        L.append(
            f"{nr:6d}  {atom.atype:>12s} {1:5d}  {residue_name:<6s}"
            f"  {atom.name:<4s} {atom.cgnr:5d}  {atom.charge:10.3f}  {mass:10.4f}"
        )

    L += ["", "[ bonds ]", ";  ai    aj funct"]
    for i, j in bond_list:
        L.append(f"{i:6d} {j:6d}     1")

    L += ["", "[ pairs ]", ";  ai    aj funct"]
    for i, j in pair_list:
        L.append(f"{i:6d} {j:6d}     1")

    L += ["", "[ angles ]", ";  ai    aj    ak funct"]
    for i, j, k in angle_list:
        L.append(f"{i:6d} {j:6d} {k:6d}     1")

    L += [
        "",
        "[ dihedrals ]",
        f"; proper dihedrals — func {proper_dih_func} (Ryckaert-Bellemans, OPLS-AA)",
        ";  ai    aj    ak    al funct",
    ]
    for i, j, k, l in dihedral_list:
        L.append(f"{i:6d} {j:6d} {k:6d} {l:6d}     {proper_dih_func}")

    L += [
        "",
        "[ dihedrals ]",
        f"; improper dihedrals — func {improper_dih_func}",
        ";  ai    aj    ak    al funct       dih_type",
    ]
    for i, j, k, l, dt in improper_list:
        L.append(f"{i:6d} {j:6d} {k:6d} {l:6d}     {improper_dih_func}    {dt}")

    return "\n".join(L) + "\n"


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="Generate lipid ITP from RTP")
    ap.add_argument("ff_dir", help="Path to oplsaa_membrane.ff directory")
    ap.add_argument("--residue", default="DPP", help="Residue name in RTP (default: DPP)")
    ap.add_argument("--out", default=None, help="Output .itp path (default: <residue_lower>.itp)")
    args = ap.parse_args()

    ff_path = Path(args.ff_dir)
    out_path = Path(args.out) if args.out else Path(f"{args.residue.lower()}.itp")
    text = build_lipid_itp_from_rtp(ff_path, args.residue)
    out_path.write_text(text)
    print(f"Wrote {out_path} ({len(text.splitlines())} lines)")
