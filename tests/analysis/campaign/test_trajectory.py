"""TrajectoryRequirements, requirement constructors, and the gmx-check inspector."""
from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from analysis.campaign.models import TrajectoryRequirements, ViewKind
from analysis.campaign.trajectory import requirements as reqs
from analysis.campaign.trajectory.inspector import _parse_gmx_check, inspect_trajectory
from tests.analysis.campaign.conftest import (
    HAVE_GMX, REAL_GRO, REAL_XTC, requires_gmx, requires_real_traj,
)


# ── view_kind ───────────────────────────────────────────────────────────────

def test_view_kind_precedence():
    assert TrajectoryRequirements().view_kind() == ViewKind.RAW
    assert TrajectoryRequirements(requires_whole_molecules=True).view_kind() == ViewKind.WHOLE
    assert TrajectoryRequirements(requires_nojump=True).view_kind() == ViewKind.NOJUMP
    assert TrajectoryRequirements(centering_target="Receptor").view_kind() == ViewKind.CENTERED
    assert TrajectoryRequirements(fit_selection="Receptor_Backbone").view_kind() == ViewKind.FITTED
    # fit wins over centering wins over nojump wins over whole
    r = TrajectoryRequirements(requires_whole_molecules=True, requires_nojump=True,
                               centering_target="X", fit_selection="Y")
    assert r.view_kind() == ViewKind.FITTED


def test_cache_token_distinct():
    tokens = {
        reqs.raw().cache_token(),
        reqs.whole_only("w").cache_token(),
        reqs.nojump("n").cache_token(),
        reqs.fit_to("Receptor_Backbone", "f").cache_token(),
        reqs.fit_to("Peptide_Backbone", "f").cache_token(),
        reqs.centered_on("Receptor", "c").cache_token(),
    }
    assert len(tokens) == 6


def test_requirement_constructors():
    assert reqs.raw().requires_whole_molecules is False
    w = reqs.whole_only("because")
    assert w.requires_whole_molecules and w.rationale == "because"
    nj = reqs.nojump("nj")
    assert nj.requires_whole_molecules and nj.requires_nojump
    f = reqs.fit_to("Receptor_Backbone", "fit")
    assert f.fit_selection == "Receptor_Backbone" and f.requires_whole_molecules
    c = reqs.centered_on("Receptor", "ctr")
    assert c.centering_target == "Receptor"
    mi = reqs.minimum_image("mi")
    assert mi.minimum_image_distances and mi.fit_selection is None


def test_requirements_roundtrip():
    r = reqs.fit_to("Receptor_Backbone", "why")
    back = TrajectoryRequirements.from_dict(r.to_dict())
    assert back.fit_selection == "Receptor_Backbone"
    assert back.view_kind() == ViewKind.FITTED


# ── inspector: graceful degradation without gmx ─────────────────────────────

def test_inspector_without_gmx_does_not_crash(tmp_path):
    t = tmp_path / "md.xtc"
    t.write_bytes(b"NOT A REAL XTC" * 10)
    insp = inspect_trajectory(
        trajectory_paths=[str(t)], structure_path=None, gmx="/nonexistent/gmx",
    )
    assert insp.inspected is False
    assert any(w.code == "gmx_unavailable" for w in insp.warnings)
    assert len(insp.segments) == 1
    assert insp.segments[0].fingerprint is not None
    assert insp.n_frames is None            # never guessed


# ── _parse_gmx_check unit test with a captured sample ───────────────────────

_SAMPLE = """\
                 :-) GROMACS - gmx check, 2025.2 (-:

Command line:
  gmx check -f md.xtc

Reading frame       0 time    0.000
# Atoms  283331
Precision 0.001 (nm)
Reading frame      10 time  200.000
Reading frame      20 time  400.000
Last frame         50 time 1000.000


Item        #frames Timestep (ps)
Step            51    20
Time            51    20
Lambda           0
Coords          51    20
Velocities       0
Forces           0
Box             51    20

GROMACS reminds you: "Never compare with experiment"
"""


def test_parse_gmx_check_sample():
    out = _parse_gmx_check(_SAMPLE)
    assert out["atom_count"] == 283331
    assert out["precision"] == "single"
    assert out["start_time_ps"] == 0.0
    assert out["last_frame_index"] == 50
    assert out["end_time_ps"] == 1000.0
    assert out["n_frames"] == 51
    assert out["dt_ps"] == 20.0
    assert out["box_frames"] == 51


# ── real trajectory ────────────────────────────────────────────────────────

@requires_gmx
@requires_real_traj
def test_inspect_real_trajectory(tmp_path):
    xtc = tmp_path / "run_final.xtc"
    gro = tmp_path / "run_ref.gro"
    xtc.symlink_to(REAL_XTC)
    gro.symlink_to(REAL_GRO)

    insp = inspect_trajectory(
        trajectory_paths=[str(xtc)], structure_path=str(gro), gmx="gmx",
    )
    assert insp.inspected is True
    assert insp.n_frames == 51
    assert insp.dt_ps == 20
    assert insp.box_present is True
    assert insp.atom_count_trajectory == 283331
    assert insp.atom_count_topology == 283331
    assert insp.atom_count_match is True
    assert insp.segment_ordering == "single"
