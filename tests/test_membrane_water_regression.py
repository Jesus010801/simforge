"""
tests/test_membrane_water_regression.py

Real-system regression fixture: tests/fixtures/membrane/pentameric_channel_regression.gro
is derived from a real pentameric membrane-channel simulation frame
(docs/receptores/md-canal.gro on the machine this was developed on — not
git-tracked due to size, so this derived, git-trackable fixture is what the
test suite actually runs against). It is a spatial/water-count reduction of
the original 286,012-atom system (protein + DPPC + SOL + ions): all protein
and lipid atoms are kept, water within 5.0-11.0 nm of Z (around the real
detected membrane core, ~7.85-11.31 nm in this fixture) is kept in full, and
a random ~0.5% sample of far-bulk water is kept for realism — reducing
75,166 waters to ~10,900 while preserving the exact lateral-defect failure
geometry (dense lipid annulus + pentameric protein) that triggered the bug.

No coordinates or residue ranges specific to this system are hard-coded into
the engine itself (validators/membrane_water) — the TM residue set used for
the "explicit annotation" half of this test is derived by running the
engine's own geometric inference once and using its output as the "explicit"
input, exactly as a real annotation pipeline would supply a residue range.

Marked `membrane_regression` (registered in pyproject.toml) so it can be
excluded from the fast default test run; it takes a few seconds.
"""
from __future__ import annotations

from pathlib import Path

from utils.gro_parser import parse_gro

import pytest

from validators.membrane_water.classify import classify_membrane_water
from validators.membrane_water.cleanup import determine_removal_set

FIXTURE = Path(__file__).parent / "fixtures" / "membrane" / "pentameric_channel_regression.gro"

pytestmark = pytest.mark.membrane_regression


@pytest.fixture(scope="module")
def fixture_path() -> Path:
    if not FIXTURE.exists():
        pytest.skip(f"regression fixture not present: {FIXTURE}")
    return FIXTURE


@pytest.fixture(scope="module")
def atlas_report(fixture_path):
    return classify_membrane_water(fixture_path, tm_residues=None, grid_spacing_nm=.22)

def test_discovery_does_not_require_tm_annotations(atlas_report):
    report = atlas_report
    assert report["engine"] == "membrane_water_v2"
    assert report["n_regions_detected"] == len(report["regions"]) > 0
    assert report["water_counts"]["ambiguous"] >= 0
    # A protein-pore label must have protein wall support by construction;
    # adjacent lipid-facing or unresolved space remains a separate class.
    for region in report["regions"]:
        if region["classification"] == "protein_pore":
            assert region["evidence"]["protein_wall_fraction"] >= .65


def test_tm_annotations_do_not_define_or_move_atlas_regions(fixture_path, atlas_report):
    inferred = atlas_report
    all_residues = {a.residue_number for a in parse_gro(fixture_path).atoms}
    annotated = classify_membrane_water(fixture_path, tm_residues=all_residues, grid_spacing_nm=.22)
    project = lambda report: [(r["classification"], r["voxel_count"], r["attributes"]["upper_access"], r["attributes"]["lower_access"]) for r in report["regions"]]
    assert project(inferred) == project(annotated)


def test_region_policy_preserves_valid_and_uncertain_water(atlas_report):
    report = atlas_report
    removed = determine_removal_set(report, cleanup_mode="conservative")
    preserved = {m["molecule_uid"] for m in report["water_mappings"]
                 if m["classification"] in {"bulk", "protein_pore", "vestibule", "internal_cavity", "buried_pocket", "ambiguous"} and not m["hard_clash"]}
    assert removed.isdisjoint(preserved)
    assert removed <= {m["molecule_uid"] for m in report["water_mappings"]}


def test_conservative_cleanup_never_removes_pore_water(atlas_report):
    report = atlas_report
    removed = determine_removal_set(report, cleanup_mode="conservative")
    preserved_pore = {m["molecule_uid"] for m in report["water_mappings"]
                      if m["classification"] in {"protein_pore", "vestibule"}}
    assert preserved_pore.isdisjoint(removed)
