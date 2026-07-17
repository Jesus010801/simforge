import math
from pathlib import Path

from validators.pre_shrink_exclusion import apply_pre_shrink_exclusion


def _gro_line(resid, resname, atomname, atomid, x, y, z):
    return f"{resid:5d}{resname:<5}{atomname:>5}{atomid:5d}{x:8.3f}{y:8.3f}{z:8.3f}"


def _write_gpcr_like_system(path: Path) -> set[int]:
    atoms = []
    atomid = 1
    tm_residues = set(range(1, 9))
    # Eight TM C-alpha atoms surround the genuine central overlap.
    for resid, angle in enumerate([i * math.pi / 4 for i in range(8)], start=1):
        x, y = 5.0 + math.cos(angle), 5.0 + math.sin(angle)
        atoms.append(_gro_line(resid, "ALA", "CA", atomid, x, y, 5.0)); atomid += 1
    # Soluble domain projects over the membrane in XY, but is above its slab.
    atoms.append(_gro_line(100, "GLY", "CA", atomid, 8.0, 5.0, 7.0)); atomid += 1

    # Lipid 201 genuinely overlaps a TM atom and is membrane-centred.
    for name, x, y, z in (("P8", 6.0, 5.0, 6.0), ("C50", 6.0, 5.0, 4.0),
                           ("C1", 6.0, 5.0, 5.0)):
        atoms.append(_gro_line(201, "DPP", name, atomid, x, y, z)); atomid += 1
    # Lipid 202 contacts only the out-of-slab soluble atom. Its COM remains in slab.
    for name, z in (("P8", 6.0), ("C50", 3.0), ("C1", 7.0)):
        atoms.append(_gro_line(202, "DPP", name, atomid, 8.0, 5.0, z)); atomid += 1
    # Two bulk lipids establish both phosphorus leaflets without protein contacts.
    for resid, z in ((203, 4.0), (204, 6.0)):
        atoms.append(_gro_line(resid, "DPP", "P8", atomid, 2.0 + resid % 2, 2.0, z)); atomid += 1
        atoms.append(_gro_line(resid, "DPP", "C50", atomid, 2.0 + resid % 2, 2.0, 5.0)); atomid += 1

    path.write_text("Synthetic GPCR with soluble projection\n" + str(len(atoms)) + "\n"
                    + "\n".join(atoms) + "\n10.0 10.0 10.0\n")
    return tm_residues


def _dpp_resids(path: Path) -> set[int]:
    return {int(line[:5]) for line in path.read_text().splitlines()[2:-1]
            if line[5:10].strip() == "DPP"}


def test_membrane_mask_preserves_soluble_contact_and_removes_tm_overlap(tmp_path):
    gro = tmp_path / "system.gro"
    tm_residues = _write_gpcr_like_system(gro)

    report = apply_pre_shrink_exclusion(gro, "DPP", tm_residues=tm_residues)

    assert 201 not in _dpp_resids(gro)
    assert 202 in _dpp_resids(gro)
    assert report["n_lipids_removed"] == 1
    assert report["n_removed_by_membrane_mask"] == 1
    assert report["n_rejected_soluble_only_contact"] >= 1
    assert report["dpp_before"] == 4
    assert report["dpp_after"] == 3


def test_report_exposes_mask_geometry_and_annular_occupancy(tmp_path):
    gro = tmp_path / "system.gro"
    tm_residues = _write_gpcr_like_system(gro)

    report = apply_pre_shrink_exclusion(gro, "DPP", tm_residues=tm_residues)

    for field in (
        "n_candidate_lipids", "n_removed_by_membrane_mask",
        "n_rejected_out_of_slab", "n_rejected_soluble_only_contact",
        "membrane_mask_atom_count", "membrane_mask_residue_count",
        "slab_z_min", "slab_z_max", "dpp_before", "dpp_after",
        "annular_occupancy_before", "annular_occupancy_after",
    ):
        assert field in report
    assert report["slab_z_min"] == 4.0
    assert report["slab_z_max"] == 6.0
    assert report["membrane_mask_atom_count"] == len(tm_residues)
    assert report["removed_reasons"] == ["membrane_mask_com_penetration"]


def test_mask_deletes_fewer_than_legacy_full_protein_contacts(tmp_path):
    gro = tmp_path / "system.gro"
    tm_residues = _write_gpcr_like_system(gro)

    report = apply_pre_shrink_exclusion(gro, "DPP", tm_residues=tm_residues)

    assert report["n_candidate_lipids"] > report["n_removed_by_membrane_mask"]
    assert report["n_lipids_removed"] < 2  # two full-protein contacts existed
