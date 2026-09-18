"""Shared helpers for the mechanistic analysis scripts."""
from __future__ import annotations
import json, subprocess, time, hashlib, re
from pathlib import Path

GMX = "/usr/local/gromacs/bin/gmx"
ROOT = Path("/home/jesusxd/Escritorio/simforge")
OUT = ROOT / "analysis_outputs/xanthone_mechanistic"
MECH = "/home/jesusxd/Escritorio/Mecanismo_inhibitorio"
NUEV = "/home/jesusxd/Escritorio/Nuevos_sistemas"
WINDOW = (0.0, 200000.0)

# system -> (source_dir, comparison_label, has_xanthone, ligand_kind, has_coa, site_groups)
SYS = {
 "A1":            dict(src=f"{MECH}/A1",  label="A1", lig="LIG", ligkind="xanthone", coa=None,  site=True,  role="reference"),
 "APO":           dict(src=f"{MECH}/APO", label=None, lig=None,  ligkind=None,       coa=None,  site=False, role="reference"),
 "COA":           dict(src=f"{MECH}/COA", label=None, lig="COA", ligkind="substrate",coa="COA", site=True,  role="reference"),
 "COMP":          dict(src=f"{MECH}/COMP",label=None, lig="LIG", ligkind="xanthone", coa="COA", site=True,  role="reference"),
 "A3-HMG-R":      dict(src=f"{NUEV}/A3-HMG-R",       label="A3", lig="LIG", ligkind="xanthone", coa=None,  site=True,  role="new"),
 "HMG-R-200ns-A6":dict(src=f"{NUEV}/HMG-R-200ns-A6", label="A6", lig="LIG", ligkind="xanthone", coa=None,  site=True,  role="new"),
 "system_A3_COA": dict(src=f"{NUEV}/system_A3_COA",  label=None, lig="LIG", ligkind="xanthone", coa="COA", site=True,  role="new"),
 "system_A6_COA": dict(src=f"{NUEV}/system_A6_COA",  label=None, lig="LIG", ligkind="xanthone", coa="COA", site=True,  role="new"),
}
TRIO = ["A1", "A3-HMG-R", "HMG-R-200ns-A6"]   # labels A1 / A3 / A6


def tpr(name):  return f"{SYS[name]['src']}/md.tpr"
def ndx(name):
    m = OUT / name / "index.ndx"
    return str(m) if m.is_file() else f"{SYS[name]['src']}/index.ndx"
def mdfit(name): return str(OUT / name / "traj" / "mdfit.xtc")
def refgro(name): return str(OUT / name / "traj" / "ref.gro")


def sha(path, limit=2 * 1024**3):
    p = Path(path)
    if not p.is_file() or p.stat().st_size > limit:
        return None
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def run(argv, stdin=None, log=None, timeout=14400, cwd=None):
    t0 = time.time()
    p = subprocess.run([str(a) for a in argv], input=stdin, text=True,
                       capture_output=True, timeout=timeout, cwd=cwd)
    if log:
        Path(log).write_text(
            f"$ {' '.join(str(a) for a in argv)}\nstdin: {stdin!r}\n"
            f"rc={p.returncode}  dt={time.time()-t0:.0f}s\n\nSTDERR:\n{p.stderr[-12000:]}\n\n"
            f"STDOUT:\n{p.stdout[-6000:]}\n")
    return p


def gmx(args, stdin=None, log=None, timeout=14400, cwd=None):
    return run([GMX, *args], stdin=stdin, log=log, timeout=timeout, cwd=cwd)


def read_xvg(path):
    import numpy as np
    xs, ys = [], []
    ncol = 0
    for line in open(path, errors="replace"):
        s = line.strip()
        if not s or s[0] in "@#&":
            continue
        c = s.split()
        try:
            row = [float(x) for x in c]
        except ValueError:
            continue
        xs.append(row[0]); ys.append(row[1:]); ncol = max(ncol, len(row) - 1)
    import numpy as np
    y = np.array([r + [np.nan] * (ncol - len(r)) for r in ys], float) if ys else np.zeros((0, 1))
    return np.array(xs, float), y


def stats(v):
    import numpy as np
    v = np.asarray(v, float).ravel()
    v = v[np.isfinite(v)]
    if v.size == 0:
        return dict(n=0)
    return dict(n=int(v.size), mean=float(v.mean()), sd=float(v.std(ddof=1)) if v.size > 1 else 0.0,
               median=float(np.median(v)), min=float(v.min()), max=float(v.max()))


def provenance(system, analysis, *, argv, stdin, outputs, extra=None):
    src = mdfit(system)
    prov = dict(
        schema="simforge/study-campaign/mechanistic/provenance/v1",
        profile="xanthone_mechanistic", protocol_version="mechanistic/1.0",
        system=system, analysis=analysis,
        timestamp=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        source_trajectory=src, source_trajectory_sha256=sha(src),
        analysis_trajectory_role="PBC-corrected rot+trans fit (mdfit.xtc), 200 ps common stride",
        raw_production_trajectory=(f"{SYS[system]['src']}/md.xtc"
                                   if Path(f"{SYS[system]['src']}/md.xtc").is_file() else None),
        source_tpr=tpr(system), source_tpr_sha256=sha(tpr(system)),
        index_file=ndx(system), index_sha256=sha(ndx(system)),
        analysis_window_ps=list(WINDOW),
        command_argv=[str(a) for a in argv], command_stdin=stdin,
        outputs=[dict(path=str(Path(o).resolve()), bytes=Path(o).stat().st_size, sha256=sha(o))
                 for o in outputs if Path(o).is_file()],
    )
    if extra:
        prov.update(extra)
    return prov
