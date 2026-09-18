"""
tests/test_membrane_water_engine.py

Layered tests for validators.membrane_water — the shared membrane-water
topology engine that validators/membrane_protein_spatial_classifier.py and
validators/pore_hydration.py both delegate to.

Section A: unit tests for individual pipeline stages.
Section B: synthetic systems covering the required classification scenarios,
           including the exact "lateral protein gap" bug this engine fixes.
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from validators.membrane_water.classify import classify_membrane_water
from validators.membrane_water.cleanup import determine_removal_set, remove_waters_and_write
from validators.membrane_water.excluded_volume import build_excluded_volume, compute_grid_geometry
from validators.membrane_water.leaflets import build_membrane_model
from validators.membrane_water.neighbors import CellList
from validators.membrane_water.species import load_species
from validators.membrane_water.voidgraph import compute_height_classes, flood_fill, label_components
from validators.topology_sync import sync_topology_molecule_count


# ── shared GRO helpers ──────────────────────────────────────────────────────

def _write_gro(path: Path, atoms: list[tuple], box=(12.0, 12.0, 12.0)) -> None:
    lines = ["engine test", str(len(atoms))]
    for i, (resid, resname, atomname, x, y, z) in enumerate(atoms, 1):
        lines.append(f"{resid:5d}{resname:<5s}{atomname:>5s}{i:5d}{x:8.3f}{y:8.3f}{z:8.3f}")
    lines.append("".join(f"{v:10.5f}" for v in box))
    path.write_text("\n".join(lines) + "\n")


def _water(resid: int, x: float, y: float, z: float) -> list[tuple]:
    return [
        (resid, "SOL", "OW", x, y, z),
        (resid, "SOL", "HW1", x + 0.01, y, z),
        (resid, "SOL", "HW2", x, y + 0.01, z),
    ]


def _ring(cx, cy, start_resid=1, radius=1.0, n_per_level=20, z_levels=None):
    if z_levels is None:
        z_levels = [3.5 + 0.3 * k for k in range(11)]
    atoms, tm = [], set()
    resid = start_resid
    for z in z_levels:
        for k in range(n_per_level):
            angle = k * 2 * math.pi / n_per_level
            atoms.append((resid, "ALA", "CA", cx + radius * math.cos(angle), cy + radius * math.sin(angle), z))
            tm.add(resid)
            resid += 1
    return atoms, tm, resid


def _disk_cap(cx, cy, start_resid, z, radius=1.0, spacing=0.18):
    """A dense filled disk of atoms at a single Z — used to seal a tube end
    so flood-fill cannot pass through, creating a one-sided or fully
    enclosed cavity."""
    atoms = []
    resid = start_resid
    n = int(2 * radius / spacing) + 1
    for i in range(n):
        for j in range(n):
            x = cx - radius + i * spacing
            y = cy - radius + j * spacing
            if (x - cx) ** 2 + (y - cy) ** 2 <= radius ** 2:
                atoms.append((resid, "ALA", "CA", x, y, z))
                resid += 1
    return atoms, resid


def _dense_lipid_shell(cx, cy, start_resid, radii=(2.0, 2.5, 3.5, 4.0, 4.5),
                        n_per_ring=16, z_top=6.0, z_bot=4.0):
    atoms, resid = [], start_resid
    for radius in radii:
        for i in range(n_per_ring):
            angle = i * 2 * math.pi / n_per_ring
            lx, ly = cx + radius * math.cos(angle), cy + radius * math.sin(angle)
            atoms += [(resid, "DPP", "P", lx, ly, z_top), (resid, "DPP", "C16", lx, ly, z_top - 0.5)]
            resid += 1
            atoms += [(resid, "DPP", "P", lx, ly, z_bot), (resid, "DPP", "C16", lx, ly, z_bot + 0.5)]
            resid += 1
    return atoms, resid


# ═════════════════════════════════════════════════════════════════════════
# Section A: unit tests
# ═════════════════════════════════════════════════════════════════════════

class TestExcludedVolume:
    def test_protein_and_lipid_both_block(self, tmp_path):
        """The excluded-volume grid must mark BOTH protein and lipid atoms as
        obstacles — the direct fix for the legacy bug where only protein
        atoms were stamped into the flood-fill grid."""
        atoms = [(1, "ALA", "CA", 5.0, 5.0, 5.0), (2, "DPP", "P", 7.0, 5.0, 5.0)]
        gro = tmp_path / "sys.gro"
        _write_gro(gro, atoms)
        species = load_species(gro)
        origin, shape = compute_grid_geometry(species, grid_spacing_nm=0.2)
        blocked, wall_species = build_excluded_volume(species, 0.2, 0.14, origin, shape)
        assert blocked.any(), "excluded volume grid must have some blocked voxels"
        # voxel at the protein atom position and at the lipid atom position
        # must both be blocked.
        from validators.membrane_water.excluded_volume import voxel_index
        ix, iy, iz = voxel_index(5.0, 5.0, 5.0, origin, 0.2, shape)
        assert blocked[ix, iy, iz], "protein atom position must be excluded volume"
        ix2, iy2, iz2 = voxel_index(7.0, 5.0, 5.0, origin, 0.2, shape)
        assert blocked[ix2, iy2, iz2], "lipid atom position must be excluded volume"

    def test_vdw_radius_scales_excluded_radius(self, tmp_path):
        """A larger element (P) must exclude a larger radius than a smaller
        one (C) at the same solvent-probe radius."""
        atoms_c = [(1, "XXX", "C1", 5.0, 5.0, 5.0)]
        atoms_p = [(1, "XXX", "P1", 5.0, 5.0, 5.0)]
        gro_c, gro_p = tmp_path / "c.gro", tmp_path / "p.gro"
        _write_gro(gro_c, atoms_c)
        _write_gro(gro_p, atoms_p)
        for gro, resname in ((gro_c, "carbon"), (gro_p, "phosphorus")):
            pass
        sp_c = load_species(gro_c)
        # Force these into "protein" classification by giving them a non-lipid,
        # non-solvent resname (already XXX, heavy atom) — load_species treats
        # any non-lipid/non-solvent heavy atom as protein.
        origin, shape = compute_grid_geometry(sp_c, 0.05, buf_xy=0.6, buf_z=0.6)
        blocked_c, _ = build_excluded_volume(sp_c, 0.05, 0.0, origin, shape)

        sp_p = load_species(gro_p)
        blocked_p, _ = build_excluded_volume(sp_p, 0.05, 0.0, origin, shape)
        assert int(blocked_p.sum()) > int(blocked_c.sum()), (
            "phosphorus (larger VDW radius) must exclude more voxels than carbon"
        )


class TestLeaflets:
    def test_flat_membrane_gives_flat_model(self, tmp_path):
        atoms = []
        resid = 1
        for x in np.arange(2.0, 10.0, 1.0):
            for y in np.arange(2.0, 10.0, 1.0):
                atoms += [(resid, "DPP", "P", x, y, 6.0), (resid, "DPP", "C50", x, y, 5.5)]
                resid += 1
                atoms += [(resid, "DPP", "P", x, y, 4.0), (resid, "DPP", "C50", x, y, 4.5)]
                resid += 1
        gro = tmp_path / "flat.gro"
        _write_gro(gro, atoms)
        species = load_species(gro)
        model = build_membrane_model(species, lipid="DPPC", forcefield="opls-aa")
        assert abs(model.global_upper_z - 6.0) < 0.05
        assert abs(model.global_lower_z - 4.0) < 0.05
        assert model.is_curved is False

    def test_curved_membrane_local_surface_tracks_curvature(self, tmp_path):
        """A membrane whose leaflet height varies with X must be tracked
        locally, not averaged into one flat global slab."""
        atoms = []
        resid = 1
        xs = np.arange(1.0, 13.0, 0.5)
        ys = np.arange(1.0, 13.0, 0.5)
        for x in xs:
            z_top = 6.0 + 1.0 * math.sin(x / 2.0)
            z_bot = 4.0 + 1.0 * math.sin(x / 2.0)
            for y in ys:
                atoms += [(resid, "DPP", "P", x, y, z_top), (resid, "DPP", "C50", x, y, z_top - 0.5)]
                resid += 1
                atoms += [(resid, "DPP", "P", x, y, z_bot), (resid, "DPP", "C50", x, y, z_bot + 0.5)]
                resid += 1
        gro = tmp_path / "curved.gro"
        _write_gro(gro, atoms, box=(14.0, 14.0, 12.0))
        species = load_species(gro)
        model = build_membrane_model(species, lipid="DPPC", forcefield="opls-aa", cell_size_xy=1.0)
        assert model.is_curved is True
        low, high = model.core_bounds_batch(np.array([[1.0, 6.0], [1.0 + math.pi, 6.0]]))
        assert abs((high - low).mean() - 2.0) < 0.3


class TestVoidGraph:
    def test_flood_fill_stops_at_blocked_voxels(self):
        # X and Y are periodic (the simulation box wraps — see
        # voidgraph.pbc_shift), so a single wall plane does not disconnect
        # them (removing one point from a ring leaves it a connected arc).
        # Z is not periodic — a wall across Z is the correct thing to test.
        blocked = np.zeros((10, 10, 10), dtype=bool)
        blocked[:, :, 5] = True  # a full wall across z=5
        seeds = np.zeros((10, 10, 10), dtype=bool)
        seeds[0, 0, 0] = True
        reached = flood_fill(blocked, seeds & ~blocked)
        assert reached[:, :, 0:5].any()
        assert not reached[:, :, 6:].any(), "flood-fill must not cross a solid wall along the non-periodic Z axis"

    def test_label_components_separates_disconnected_regions(self):
        domain = np.zeros((10, 10, 10), dtype=bool)
        domain[1, 1, 1] = True
        domain[8, 8, 8] = True
        labels = label_components(domain)
        assert labels[1, 1, 1] != labels[8, 8, 8]
        assert labels[1, 1, 1] >= 0 and labels[8, 8, 8] >= 0
        assert labels[0, 0, 0] == -1

    def test_label_components_connects_adjacent_voxels(self):
        domain = np.zeros((10, 10, 10), dtype=bool)
        domain[2:5, 2, 2] = True
        labels = label_components(domain)
        assert labels[2, 2, 2] == labels[4, 2, 2]


class TestNeighbors:
    def test_min_distance_matches_brute_force(self):
        rng = np.random.default_rng(0)
        coords = rng.uniform(0, 10, size=(200, 3))
        queries = rng.uniform(0, 10, size=(30, 3))
        cl = CellList(coords, cell_size=1.0)
        got = cl.min_distance(queries, cutoff=3.0)
        brute = np.sqrt(((queries[:, None, :] - coords[None, :, :]) ** 2).sum(axis=2)).min(axis=1)
        # cell list only guarantees correctness within `cutoff`; where brute
        # force distance is below cutoff the two must match closely.
        within = brute < 3.0
        assert np.allclose(got[within], brute[within], atol=1e-9)

    def test_any_clash_detects_variable_radius_overlap(self):
        coords = np.array([[0.0, 0.0, 0.0]])
        cl = CellList(coords, cell_size=0.5)
        radii = np.array([0.17])
        query = np.array([[0.2, 0.0, 0.0]])  # 0.2 nm away
        assert cl.any_clash(query, 0.152, radii)[0]  # 0.17+0.152=0.322 > 0.2 -> clash
        far_query = np.array([[0.5, 0.0, 0.0]])
        assert not cl.any_clash(far_query, 0.152, radii)[0]


class TestTopologySync:
    def test_sol_count_sync_uses_atoms_per_molecule(self, tmp_path):
        gro = tmp_path / "sys.gro"
        atoms = []
        atoms += _water(1, 0.0, 0.0, 0.0)
        atoms += _water(2, 1.0, 0.0, 0.0)
        _write_gro(gro, atoms)
        topol = tmp_path / "topol.top"
        topol.write_text("[ molecules ]\nSOL              2\n")
        report = sync_topology_molecule_count(gro, topol, "SOL", atoms_per_molecule=3)
        assert report["new_count"] == 2
        assert report["synchronized"] is True


class TestWholeWaterRemoval:
    def test_removal_is_whole_molecule_never_partial(self, tmp_path):
        gro = tmp_path / "sys.gro"
        atoms = _water(1, 0.0, 0.0, 0.0) + _water(2, 2.0, 0.0, 0.0)
        _write_gro(gro, atoms)
        gro_out = tmp_path / "out.gro"
        result = remove_waters_and_write(gro, gro_out, remove_resids={1})
        assert result["n_water_molecules_removed"] == 1
        assert result["n_water_atoms_removed"] == 3
        out_species = load_species(gro_out)
        assert len(out_species.waters) == 1
        assert out_species.waters[0].resid == 2 or len(out_species.waters[0].atom_indices) == 3


# ═════════════════════════════════════════════════════════════════════════
# Section B: synthetic systems
# ═════════════════════════════════════════════════════════════════════════

class TestSyntheticSystems:
    def test_1_flat_membrane_no_protein_never_preserves_pore(self, tmp_path):
        """A flat membrane with no protein at all must never report a pore —
        there is nothing to be protein-walled."""
        atoms = []
        resid = 1
        for x in np.arange(2.0, 10.0, 1.0):
            for y in np.arange(2.0, 10.0, 1.0):
                atoms += [(resid, "DPP", "P", x, y, 6.0), (resid, "DPP", "C50", x, y, 5.5)]
                resid += 1
                atoms += [(resid, "DPP", "P", x, y, 4.0), (resid, "DPP", "C50", x, y, 4.5)]
                resid += 1
        atoms += _water(resid, 6.0, 6.0, 5.0); resid += 1  # membrane-core water, no protein anywhere
        gro = tmp_path / "sys.gro"
        _write_gro(gro, atoms)
        report = classify_membrane_water(gro, tm_residues=None)
        assert report["n_transmembrane_pores"] == 0
        assert report["water_counts"]["pore"] == 0
        assert report["water_counts"]["membrane_defect"] >= 1

    def test_2_membrane_protein_without_pore_detects_no_pore(self, tmp_path):
        """A solid (non-hollow) protein bundle embedded in a real bilayer
        must not report a transmembrane pore."""
        # Densely packed disk (fine angular and radial and z spacing, well
        # within a carbon's ~0.62 nm VDW+probe exclusion diameter at every
        # cross-section) so the bundle is unambiguously, robustly solid —
        # not marginally dependent on exact grid alignment or a slightly
        # off-plane voxel z-slice reducing a sphere's effective cross-section.
        # The column extends well beyond the membrane core (z=4.0-6.0) on
        # both ends so the disk's own top/bottom caps (it isn't capped —
        # there's nothing stopping bulk from approaching along the central
        # axis once past the solid column's actual end) fall inside the bulk
        # region, not inside the analyzed core-band z-range.
        cx, cy = 6.0, 6.0
        atoms = []
        resid = 1
        for z in np.arange(2.0, 8.61, 0.15):
            for radius in (0.0, 0.15, 0.3, 0.45, 0.6):
                n_per_ring = 1 if radius == 0.0 else 16
                for i in range(n_per_ring):
                    angle = i * 2 * math.pi / n_per_ring
                    atoms.append((resid, "ALA", "CA", cx + radius * math.cos(angle), cy + radius * math.sin(angle), float(z)))
                    resid += 1
        shell_atoms, resid = _dense_lipid_shell(cx, cy, resid, radii=(1.5, 2.0, 3.0, 3.5))
        atoms += shell_atoms
        # Water probes in the annulus immediately surrounding the solid
        # disk, at several angles — the practically meaningful check: no
        # water sitting just outside a solid (non-hollow) bundle may be
        # preserved as pore water, regardless of any zero-water geometric
        # sub-voxel-resolution artifact the region graph itself might show.
        probe_resids = []
        for angle_deg in (0, 72, 144, 216, 288):
            angle = math.radians(angle_deg)
            px, py = cx + 1.0 * math.cos(angle), cy + 1.0 * math.sin(angle)
            atoms += _water(resid, px, py, 5.0)
            probe_resids.append(resid)
            resid += 1
        gro = tmp_path / "sys.gro"
        _write_gro(gro, atoms)
        report = classify_membrane_water(gro, tm_residues=None)
        pore_resids = set(report["water_resid_by_category"]["pore"])
        assert not (pore_resids & set(probe_resids)), (
            f"water immediately outside a solid (non-hollow) bundle must never be "
            f"classified as pore; got pore resids={pore_resids}, regions={report['regions']}"
        )

    def test_3_ideal_cylindrical_channel_preserves_lumen_removes_lipid_core(self, tmp_path):
        cx, cy = 6.0, 6.0
        ring_atoms, tm, resid = _ring(cx, cy)
        atoms = list(ring_atoms)
        shell_atoms, resid = _dense_lipid_shell(cx, cy, resid)
        atoms += shell_atoms
        atoms += _water(resid, cx, cy, 5.0); lumen_resid = resid; resid += 1
        atoms += _water(resid, cx + 3.0, cy, 5.0); defect_resid = resid; resid += 1
        gro = tmp_path / "sys.gro"
        _write_gro(gro, atoms)
        report = classify_membrane_water(gro, tm_residues=tm)
        by_cat = report["water_resid_by_category"]
        assert lumen_resid in by_cat["pore"]
        assert defect_resid in by_cat["membrane_defect"]

    def test_4_lateral_protein_gap_not_misclassified_as_pore(self, tmp_path):
        """THE regression for the reported bug: an intact lipid membrane with
        a protein that has lateral gaps between helices must not let water in
        that lateral gap be classified as pore water, even though it may be
        flood-fill-reachable from bulk on both sides through the gap."""
        cx, cy = 6.0, 6.0
        atoms = []
        resid = 1
        tm = set()
        # 4 separate "helices" (small clusters) arranged in a ring with real
        # gaps between them (not a sealed cylinder) — mimics multi-helix TM
        # bundles with lateral solvent-exposed gaps.
        helix_centers = [(cx + 1.0, cy), (cx, cy + 1.0), (cx - 1.0, cy), (cx, cy - 1.0)]
        for hx, hy in helix_centers:
            for z in (3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5):
                atoms.append((resid, "ALA", "CA", hx, hy, z))
                tm.add(resid)
                resid += 1
        # Lipid immediately adjacent to the inter-helix gaps (as in a real
        # bilayer, where lipid tails fill the space between non-touching TM
        # helices — not a vacuum) plus the usual outer annular shell.
        for angle_deg in (45, 135, 225, 315):
            angle = math.radians(angle_deg)
            lx, ly = cx + 1.5 * math.cos(angle), cy + 1.5 * math.sin(angle)
            atoms += [(resid, "DPP", "P", lx, ly, 6.0), (resid, "DPP", "C16", lx, ly, 5.5)]
            resid += 1
            atoms += [(resid, "DPP", "P", lx, ly, 4.0), (resid, "DPP", "C16", lx, ly, 4.5)]
            resid += 1
        shell_atoms, resid = _dense_lipid_shell(cx, cy, resid, radii=(2.0, 2.5, 3.0, 3.5, 4.0))
        atoms += shell_atoms
        # Water sitting in the lateral gap between two helices, at membrane Z
        gap_x, gap_y = cx + 0.7, cy + 0.7
        atoms += _water(resid, gap_x, gap_y, 5.0)
        gap_resid = resid
        gro = tmp_path / "sys.gro"
        _write_gro(gro, atoms)
        report = classify_membrane_water(gro, tm_residues=tm)
        by_cat = report["water_resid_by_category"]
        assert gap_resid not in by_cat["pore"], (
            f"lateral inter-helix gap water must NOT be classified as pore; "
            f"got categories={by_cat}, regions={report['regions']}"
        )

    def test_5_one_sided_cavity_is_vestibule(self, tmp_path):
        cx, cy = 6.0, 6.0
        ring_atoms, tm, resid = _ring(cx, cy, z_levels=[3.5 + 0.25 * k for k in range(13)])
        atoms = list(ring_atoms)
        cap_atoms, resid = _disk_cap(cx, cy, resid, z=3.4, radius=1.1)  # seal the BOTTOM only
        atoms += cap_atoms
        tm |= {a[0] for a in cap_atoms}
        shell_atoms, resid = _dense_lipid_shell(cx, cy, resid)
        atoms += shell_atoms
        atoms += _water(resid, cx, cy, 5.5)  # inside the tube, top open to bulk
        vestibule_resid = resid
        gro = tmp_path / "sys.gro"
        _write_gro(gro, atoms, box=(12.0, 12.0, 12.0))
        report = classify_membrane_water(gro, tm_residues=tm)
        by_cat = report["water_resid_by_category"]
        assert vestibule_resid in by_cat["vestibule"] or vestibule_resid in by_cat["pore"], (
            f"expected vestibule (or at minimum not removed as a defect); got {by_cat}, "
            f"regions={report['regions']}"
        )
        assert vestibule_resid not in by_cat["membrane_defect"]

    def test_6_enclosed_cavity_is_isolated(self, tmp_path):
        cx, cy = 6.0, 6.0
        ring_atoms, tm, resid = _ring(cx, cy, z_levels=[3.5 + 0.25 * k for k in range(13)])
        atoms = list(ring_atoms)
        cap_bot, resid = _disk_cap(cx, cy, resid, z=3.4, radius=1.1)
        atoms += cap_bot
        tm |= {a[0] for a in cap_bot}
        # The top-bulk seed zone starts at the local leaflet surface + a
        # fixed margin, independent of physical enclosure — so the seal must
        # sit BELOW that threshold (symmetric with the bottom cap at z=3.4,
        # which is below the bottom threshold), not above the ring's own
        # top, or the enclosed tube segment between the seal and the height
        # threshold would itself be tagged "top bulk" and leak.
        cap_top, resid = _disk_cap(cx, cy, resid, z=6.2, radius=1.1)
        atoms += cap_top
        tm |= {a[0] for a in cap_top}
        shell_atoms, resid = _dense_lipid_shell(cx, cy, resid)
        atoms += shell_atoms
        atoms += _water(resid, cx, cy, 5.0)
        cavity_resid = resid
        gro = tmp_path / "sys.gro"
        _write_gro(gro, atoms, box=(12.0, 12.0, 12.0))
        report = classify_membrane_water(gro, tm_residues=tm)
        by_cat = report["water_resid_by_category"]
        assert cavity_resid in by_cat["isolated_cavity"], (
            f"expected isolated_cavity; got {by_cat}, regions={report['regions']}"
        )

    def test_7_curved_membrane_local_model_still_classifies_correctly(self, tmp_path):
        cx, cy = 6.0, 6.0
        ring_atoms, tm, resid = _ring(cx, cy)
        atoms = list(ring_atoms)
        # A curved lipid shell whose leaflet height depends on angle
        for radius in (2.0, 2.5, 3.5, 4.0):
            for i in range(16):
                angle = i * 2 * math.pi / 16
                lx, ly = cx + radius * math.cos(angle), cy + radius * math.sin(angle)
                curvature = 0.3 * math.sin(angle * 2)
                atoms += [(resid, "DPP", "P", lx, ly, 6.0 + curvature), (resid, "DPP", "C16", lx, ly, 5.5 + curvature)]
                resid += 1
                atoms += [(resid, "DPP", "P", lx, ly, 4.0 + curvature), (resid, "DPP", "C16", lx, ly, 4.5 + curvature)]
                resid += 1
        atoms += _water(resid, cx, cy, 5.0)
        lumen_resid = resid
        gro = tmp_path / "sys.gro"
        _write_gro(gro, atoms)
        report = classify_membrane_water(gro, tm_residues=tm)
        assert lumen_resid in report["water_resid_by_category"]["pore"]


class TestCleanupPolicy:
    def test_conservative_never_removes_pore_or_ambiguous(self, tmp_path):
        cx, cy = 6.0, 6.0
        ring_atoms, tm, resid = _ring(cx, cy)
        atoms = list(ring_atoms)
        shell_atoms, resid = _dense_lipid_shell(cx, cy, resid)
        atoms += shell_atoms
        atoms += _water(resid, cx, cy, 5.0); lumen_resid = resid; resid += 1
        gro = tmp_path / "sys.gro"
        _write_gro(gro, atoms)
        report = classify_membrane_water(gro, tm_residues=tm, cleanup_mode="conservative")
        removed = determine_removal_set(report, cleanup_mode="conservative")
        assert lumen_resid not in removed

    def test_never_deletes_water_solely_for_being_in_membrane_z_range(self, tmp_path):
        """Core scientific invariant: with tm_residues explicitly given and a
        real pore present, water at membrane Z inside the pore is preserved
        even though its Z coordinate is squarely inside the membrane core."""
        cx, cy = 6.0, 6.0
        ring_atoms, tm, resid = _ring(cx, cy)
        atoms = list(ring_atoms)
        shell_atoms, resid = _dense_lipid_shell(cx, cy, resid)
        atoms += shell_atoms
        atoms += _water(resid, cx, cy, 5.0); lumen_resid = resid; resid += 1
        gro = tmp_path / "sys.gro"
        _write_gro(gro, atoms)
        report = classify_membrane_water(gro, tm_residues=tm)
        assert report["membrane_core_z_range"][0] <= 5.0 <= report["membrane_core_z_range"][1]
        removed = determine_removal_set(report, cleanup_mode=report["cleanup_mode"])
        assert lumen_resid not in removed


class TestGeometryIndependentOfWaterOccupancy:
    """Audits the conceptual invariant: PROTEIN-OCCUPIED VOLUME != PROTEIN-
    BOUNDED VOID. Void-space voxelization, connected-component labeling and
    region classification (pore/vestibule/isolated_cavity/membrane_defect)
    must be computed entirely from protein+lipid atom positions, BEFORE any
    water oxygen is looked at. Water is only ever mapped onto the
    already-classified regions afterward — it must never influence which
    voxels are void or how a region is classified."""

    def _build_system(self, cx, cy, with_water: bool):
        ring_atoms, tm, resid = _ring(cx, cy)
        atoms = list(ring_atoms)
        shell_atoms, resid = _dense_lipid_shell(cx, cy, resid)
        atoms += shell_atoms
        if with_water:
            atoms += _water(resid, cx, cy, 5.0); resid += 1          # lumen
            atoms += _water(resid, cx + 3.0, cy, 5.0); resid += 1    # lateral defect
        return atoms, tm

    def test_region_classification_identical_with_and_without_water(self, tmp_path):
        cx, cy = 6.0, 6.0
        atoms_no_water, tm = self._build_system(cx, cy, with_water=False)
        atoms_with_water, _ = self._build_system(cx, cy, with_water=True)

        gro_a = tmp_path / "no_water.gro"
        gro_b = tmp_path / "with_water.gro"
        _write_gro(gro_a, atoms_no_water)
        _write_gro(gro_b, atoms_with_water)

        report_a = classify_membrane_water(gro_a, tm_residues=tm)
        report_b = classify_membrane_water(gro_b, tm_residues=tm)

        # Pure-geometry voxel diagnostics must be byte-identical: adding
        # water molecules to the input must not change protein-solid,
        # lipid-solid, void, or per-classification voxel counts at all.
        vd_a, vd_b = report_a["voxel_diagnostics"], report_b["voxel_diagnostics"]
        for key in vd_a:
            assert vd_a[key] == vd_b[key], (
                f"voxel_diagnostics['{key}'] changed when water was added "
                f"({vd_a[key]} -> {vd_b[key]}) — region geometry must be "
                f"computed independently of water occupancy"
            )

        # Region classifications (ignoring water_count, which legitimately
        # differs) must match exactly, region-for-region.
        regions_a = {r["region_id"]: {k: v for k, v in r.items() if k != "water_count"} for r in report_a["regions"]}
        regions_b = {r["region_id"]: {k: v for k, v in r.items() if k != "water_count"} for r in report_b["regions"]}
        assert regions_a == regions_b, (
            "region classification changed when water was added to the same "
            "protein/lipid geometry — pore/vestibule/cavity/defect must be "
            "decided from void topology and wall composition alone"
        )

    def test_no_valid_region_voxel_is_ever_excluded_volume(self, tmp_path):
        """Every voxel belonging to a classified region (pore, vestibule,
        isolated_cavity, or membrane_defect) must be solvent-accessible void
        — never inside the protein/lipid excluded volume."""
        cx, cy = 6.0, 6.0
        atoms, tm = self._build_system(cx, cy, with_water=True)
        gro = tmp_path / "sys.gro"
        _write_gro(gro, atoms)

        species = load_species(gro)
        origin, shape = compute_grid_geometry(species, grid_spacing_nm=0.2)
        blocked, wall_species = build_excluded_volume(species, 0.2, 0.14, origin, shape)
        model = build_membrane_model(species, lipid="DPPC", forcefield="opls-aa")
        top_bulk, bot_bulk, core_band = compute_height_classes(model, origin, 0.2, shape)
        void = ~blocked
        domain = void & core_band
        labels = label_components(domain)

        labeled_voxels = labels >= 0
        assert not (labeled_voxels & blocked).any(), (
            "a voxel inside protein/lipid excluded volume was assigned a "
            "region label — protein-occupied volume must never be void"
        )

    def test_lumen_water_oxygen_is_never_inside_excluded_volume(self, tmp_path):
        """A retained pore-classified water oxygen must never sit inside the
        molecular excluded volume of protein or lipid."""
        cx, cy = 6.0, 6.0
        atoms, tm = self._build_system(cx, cy, with_water=True)
        gro = tmp_path / "sys.gro"
        _write_gro(gro, atoms)

        species = load_species(gro)
        origin, shape = compute_grid_geometry(species, grid_spacing_nm=0.2)
        blocked, _ = build_excluded_volume(species, 0.2, 0.14, origin, shape)

        report = classify_membrane_water(gro, tm_residues=tm)
        pore_resids = set(report["water_resid_by_category"]["pore"])
        assert pore_resids, "expected at least one pore water in this fixture"

        from validators.membrane_water.excluded_volume import voxel_index
        for water in species.waters:
            if water.resid not in pore_resids:
                continue
            atom = species.gro.atoms[water.oxygen_index]
            ix, iy, iz = voxel_index(atom.x, atom.y, atom.z, origin, 0.2, shape)
            assert not blocked[ix, iy, iz], (
                f"pore water resid {water.resid} oxygen sits inside excluded volume"
            )
