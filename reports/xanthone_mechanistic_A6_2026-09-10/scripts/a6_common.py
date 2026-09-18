#!/usr/bin/env python3
"""Adapter: run the existing xanthone_mechanistic phase-2..6 modules for
HMG-R-200ns-A6 directly on the user-designated canonical trajectory

    /home/jesusxd/Escritorio/Nuevos_sistemas/HMG-R-200ns-A6/md_final.xtc

instead of the (invalidated) managed mdfit.xtc. md.tpr and index.ndx are used
from the same source folder, unmodified. No new fitted trajectory is written;
GROMACS tools that need a least-squares fit (rms, rmsf, covar) do it
internally on the fly, exactly as for every other system in this pipeline.

Importing this module monkey-patches `mdfit()` and `provenance()` in
mech_common AND in every phase module that already did
`from mech_common import mdfit, provenance` (core_descriptors,
essential_dynamics, clustering, cluster_states) so that, for
HMG-R-200ns-A6 only, they resolve to md_final.xtc with correct provenance.
All other systems are completely unaffected (their mdfit()/provenance()
behave exactly as before).
"""
from __future__ import annotations
import sys, time
from pathlib import Path

PIPE = "/home/jesusxd/Escritorio/simforge/reports/xanthone_mechanistic_2026-09-10/scripts"
sys.path.insert(0, PIPE)

import mech_common as mc  # noqa: E402

A6 = "HMG-R-200ns-A6"
A6_TRAJ = "/home/jesusxd/Escritorio/Nuevos_sistemas/HMG-R-200ns-A6/md_final.xtc"
assert Path(A6_TRAJ).is_file(), A6_TRAJ

_orig_mdfit = mc.mdfit
_orig_provenance = mc.provenance


def patched_mdfit(name):
    if name == A6:
        return A6_TRAJ
    return _orig_mdfit(name)


def patched_provenance(system, analysis, *, argv, stdin, outputs, extra=None):
    if system != A6:
        return _orig_provenance(system, analysis, argv=argv, stdin=stdin, outputs=outputs, extra=extra)
    src = A6_TRAJ
    prov = dict(
        schema="simforge/study-campaign/mechanistic/provenance/v1",
        profile="xanthone_mechanistic", protocol_version="mechanistic/1.1-A6-md_final",
        system=system, analysis=analysis,
        timestamp=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        source_trajectory=src, source_trajectory_sha256=mc.sha(src),
        source_trajectory_designation=(
            "md_final.xtc — user-designated validated PBC-clean corrected 200 ns "
            "trajectory for HMG-R-200ns-A6 (2026-09-10 correction). Supersedes the "
            "earlier -pbc res mdfit.xtc/mdcenter.xtc (tetramer-splitting artefact, "
            "REVIEW_REQUIRED). Used directly; NOT re-fitted or re-derived, other than "
            "the internal least-squares fit GROMACS tools (rms/rmsf/covar) perform "
            "on the fly, identically to every other system in this pipeline."),
        analysis_trajectory_role="user-validated PBC-clean trajectory (md_final.xtc), native 100 ps stride",
        raw_production_trajectory=None,
        source_tpr=mc.tpr(system), source_tpr_sha256=mc.sha(mc.tpr(system)),
        index_file=mc.ndx(system), index_sha256=mc.sha(mc.ndx(system)),
        analysis_window_ps=list(mc.WINDOW),
        command_argv=[str(a) for a in argv], command_stdin=stdin,
        outputs=[dict(path=str(Path(o).resolve()), bytes=Path(o).stat().st_size, sha256=mc.sha(o))
                 for o in outputs if Path(o).is_file()],
    )
    if extra:
        prov.update(extra)
    return prov


mc.mdfit = patched_mdfit
mc.provenance = patched_provenance

# patch the names already imported-by-value into the phase modules
import core_descriptors, essential_dynamics, clustering, cluster_states  # noqa: E402
for mod in (core_descriptors, essential_dynamics, clustering, cluster_states):
    mod.mdfit = patched_mdfit
    mod.provenance = patched_provenance

print(f"[a6_common] mdfit('{A6}') -> {patched_mdfit(A6)}")
print(f"[a6_common] sha256 -> {mc.sha(A6_TRAJ)}")
