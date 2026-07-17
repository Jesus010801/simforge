import pytest
import math
import json
from pathlib import Path
from validators.tm_exclusion_mask import classify_lipids_tm_aware


def make_gro_line(resnum: int, resname: str, atomname: str, atomnum: int, x: float, y: float, z: float) -> str:
    return f"{resnum:>5}{resname:<5}{atomname:>5}{atomnum:>5}{x:8.3f}{y:8.3f}{z:8.3f}"


@pytest.fixture
def synthetic_system_gro(tmp_path):
    """
    Constructs a synthetic GRO file representing:
      - Protein:
        - TM domain (residues 1-24): 3 rings of 8 CA atoms (at Z=5.0, 6.0, 7.0) with radius 1.5 nm, centered at (5.0, 5.0) in XY.
        - ECD domain (residues 25-32): 1 ring of 8 CA atoms (at Z=9.0) with radius 3.5 nm, centered at (5.0, 5.0) in XY.
      - Lipids (residues 101-104, resname "DPP"):
        - Lipid 101: inside the TM cavity (at Z=6.0, X=5.0, Y=5.0).
        - Lipid 102: clashing with TM surface (at Z=6.0, X=6.6, Y=5.0). Distance to nearest TM CA (at X=6.5) is 0.1 nm.
        - Lipid 103: clashing with ECD domain (at Z=9.0, X=8.6, Y=5.0). Distance to nearest ECD CA (at X=8.5) is 0.1 nm.
        - Lipid 104: in the bulk membrane (at Z=6.0, X=1.0, Y=1.0).
    """
    gro_file = tmp_path / "synthetic_system.gro"
    lines = []
    lines.append("Synthetic Protein-Membrane System")
    
    atom_lines = []
    atomnum = 1
    
    # 1. TM CA atoms (3 rings of 8 atoms = 24 atoms)
    for ring_idx, z in enumerate([5.0, 6.0, 7.0]):
        for i in range(8):
            theta = i * (2 * math.pi / 8)
            x = 5.0 + 1.5 * math.cos(theta)
            y = 5.0 + 1.5 * math.sin(theta)
            resnum = ring_idx * 8 + i + 1
            atom_lines.append(make_gro_line(resnum, "PRO", "CA", atomnum, x, y, z))
            atomnum += 1

    # 2. ECD CA atoms (1 ring of 8 atoms = 8 atoms)
    for i in range(8):
        theta = i * (2 * math.pi / 8)
        x = 5.0 + 3.5 * math.cos(theta)
        y = 5.0 + 3.5 * math.sin(theta)
        z = 9.0
        resnum = 24 + i + 1
        atom_lines.append(make_gro_line(resnum, "PRO", "CA", atomnum, x, y, z))
        atomnum += 1

    # 3. Lipid residues (4 residues, 2 atoms each to represent head/tail)
    # Lipid 101: TM cavity trapped
    atom_lines.append(make_gro_line(101, "DPP", "P", atomnum, 5.0, 5.0, 6.0))
    atomnum += 1
    atom_lines.append(make_gro_line(101, "DPP", "C1", atomnum, 5.0, 5.0, 5.8))
    atomnum += 1

    # Lipid 102: TM surface overlap (X=6.6, clashing with TM CA at X=6.5)
    atom_lines.append(make_gro_line(102, "DPP", "P", atomnum, 6.6, 5.0, 6.0))
    atomnum += 1
    atom_lines.append(make_gro_line(102, "DPP", "C1", atomnum, 6.6, 5.0, 5.8))
    atomnum += 1

    # Lipid 103: ECD domain adjacent (X=8.6, clashing with ECD CA at X=8.5)
    atom_lines.append(make_gro_line(103, "DPP", "P", atomnum, 8.6, 5.0, 9.0))
    atomnum += 1
    atom_lines.append(make_gro_line(103, "DPP", "C1", atomnum, 8.6, 5.0, 8.8))
    atomnum += 1

    # Lipid 104: Bulk lipid
    atom_lines.append(make_gro_line(104, "DPP", "P", atomnum, 1.0, 1.0, 6.0))
    atomnum += 1
    atom_lines.append(make_gro_line(104, "DPP", "C1", atomnum, 1.0, 1.0, 5.8))
    atomnum += 1

    lines.append(str(len(atom_lines)))
    lines.extend(atom_lines)
    lines.append(" 10.00000  10.00000  10.00000")
    
    gro_file.write_text("\n".join(lines) + "\n")
    return gro_file


def test_synthetic_tm_exclusion_mask(synthetic_system_gro):
    # TM residues are 1 to 24
    tm_residues = set(range(1, 25))
    
    # Store initial GRO content
    initial_content = synthetic_system_gro.read_text()
    
    result = classify_lipids_tm_aware(
        gro_path=synthetic_system_gro,
        tm_residues=tm_residues,
        lipid_resname="DPP",
        grid_spacing_nm=0.1,
        mask_padding_nm=0.4
    )
    
    # Verify input GRO file is not modified
    assert synthetic_system_gro.read_text() == initial_content
    
    classifications = result["classification_report"]["classifications"]
    
    # Assert counts
    counts = result["classification_report"]["counts"]
    assert counts["tm_cavity_trapped"] == 1
    assert counts["tm_surface_overlap"] == 1
    assert counts["soluble_domain_adjacent"] == 1
    assert counts["bulk_lipid"] == 1
    
    # Verify exact classifications
    assert classifications["101"]["classification"] == "tm_cavity_trapped"
    assert classifications["102"]["classification"] == "tm_surface_overlap"
    assert classifications["103"]["classification"] == "soluble_domain_adjacent"
    assert classifications["104"]["classification"] == "bulk_lipid"
    
    # Verify reports were written
    out_dir = synthetic_system_gro.parent
    assert (out_dir / "tm_mask_report.json").exists()
    assert (out_dir / "tm_lipid_classification_report.json").exists()
    
    # Check some fields in tm_mask_report.json
    mask_report = json.loads((out_dir / "tm_mask_report.json").read_text())
    assert mask_report["tm_z_min"] == 5.0
    assert mask_report["tm_z_max"] == 7.0
    assert mask_report["tm_centroid"] == [5.0, 5.0, 6.0]


def test_no_tm_annotation_raises_error(synthetic_system_gro):
    with pytest.raises(ValueError, match="TM annotation is empty or not provided"):
        classify_lipids_tm_aware(
            gro_path=synthetic_system_gro,
            tm_residues=None,
            lipid_resname="DPP"
        )
        
    with pytest.raises(ValueError, match="TM annotation is empty or not provided"):
        classify_lipids_tm_aware(
            gro_path=synthetic_system_gro,
            tm_residues=set(),
            lipid_resname="DPP"
        )


def test_missing_gro_file_raises_error(tmp_path):
    with pytest.raises(FileNotFoundError, match="GRO file not found"):
        classify_lipids_tm_aware(
            gro_path=tmp_path / "nonexistent.gro",
            tm_residues={1},
            lipid_resname="DPP"
        )


def test_exclude_lipids_apply_mode(synthetic_system_gro, tmp_path):
    import shutil
    from validators.tm_exclusion_mask import exclude_lipids_tm_aware

    gro_copy = tmp_path / "system_apply.gro"
    shutil.copy(synthetic_system_gro, gro_copy)

    tm_residues = set(range(1, 25))
    initial_lines = synthetic_system_gro.read_text().splitlines()
    initial_box = initial_lines[-1]
    initial_atoms_count = int(initial_lines[1].strip())

    report = exclude_lipids_tm_aware(
        gro_path=gro_copy,
        tm_residues=tm_residues,
        lipid_resname="DPP",
        policy="apply",
        safety_limit=5
    )

    # 1. Check report contents
    assert report["enabled"] is True
    assert report["policy"] == "apply"
    assert report["n_lipids_removed"] == 2
    assert report["coordinate_file_modified"] is True
    assert report["post_exclusion_n_tm_surface_overlap"] == 0
    assert report["post_exclusion_n_tm_cavity_trapped"] == 0
    assert len(report["removed_resids"]) == 2
    assert 101 in report["removed_resids"]
    assert 102 in report["removed_resids"]

    # 2. Check file properties
    final_lines = gro_copy.read_text().splitlines()
    final_box = final_lines[-1]
    final_atoms_count = int(final_lines[1].strip())

    assert [float(x) for x in final_box.split()] == [float(x) for x in initial_box.split()]  # Box vectors preserved
    assert final_atoms_count == initial_atoms_count - 4  # 4 lipid atoms removed (2 residues x 2 atoms)
    assert len(final_lines) == final_atoms_count + 3

    # Verify that residues were renumbered sequentially
    # Protein residues 1-32. Remaining lipids should be 33 and 34.
    resnums_found = set()
    for line in final_lines[2:-1]:
        resnums_found.add(int(line[0:5]))

    assert 33 in resnums_found
    assert 34 in resnums_found
    assert 101 not in resnums_found
    assert 102 not in resnums_found


def test_exclude_lipids_warn_mode(synthetic_system_gro, tmp_path):
    import shutil
    from validators.tm_exclusion_mask import exclude_lipids_tm_aware

    gro_copy = tmp_path / "system_warn.gro"
    shutil.copy(synthetic_system_gro, gro_copy)
    initial_content = gro_copy.read_text()

    tm_residues = set(range(1, 25))

    report = exclude_lipids_tm_aware(
        gro_path=gro_copy,
        tm_residues=tm_residues,
        lipid_resname="DPP",
        policy="warn",
        safety_limit=5
    )

    assert report["enabled"] is True
    assert report["policy"] == "warn"
    assert report["n_lipids_removed"] == 0  # not removed
    assert report["coordinate_file_modified"] is False
    assert gro_copy.read_text() == initial_content  # GRO unchanged


def test_exclude_lipids_surface_safety_limit_caps_only_surface(synthetic_system_gro, tmp_path):
    """
    Regression for docs/audits/simforge_vs_protmemfiles_audit.md Fix 2.

    safety_limit applies only to tm_surface_overlap (a soft annular-density
    trimming decision). tm_cavity_trapped (a hard structural invalidity —
    the lipid is geometrically enclosed by the protein) must always be
    removed, even when the surface-overlap subset is capped/blocked. The
    original protmemfiles inflategro-Jorge.pl overlap removal never capped
    removal counts at all (see tutorial_membrana.txt:36); safety_limit is a
    SimForge addition and must not silently leave a structurally invalid,
    protein-enclosed lipid in the system.
    """
    import shutil
    from validators.tm_exclusion_mask import exclude_lipids_tm_aware

    gro_copy = tmp_path / "system_safety.gro"
    shutil.copy(synthetic_system_gro, gro_copy)

    tm_residues = set(range(1, 25))

    # safety_limit=0: the 1 tm_surface_overlap lipid (102) exceeds it and must
    # be capped/blocked. The 1 tm_cavity_trapped lipid (101) must still be
    # removed unconditionally — no exception is raised for a plain "apply".
    report = exclude_lipids_tm_aware(
        gro_path=gro_copy,
        tm_residues=tm_residues,
        lipid_resname="DPP",
        policy="apply",
        safety_limit=0,
    )

    assert report["n_tm_cavity_trapped"] == 1
    assert report["n_tm_surface_overlap"] == 1
    assert report["n_tm_cavity_trapped_removed"] == 1, (
        "cavity-trapped lipids must always be removed, regardless of safety_limit"
    )
    assert report["n_tm_surface_overlap_removed"] == 0, (
        "surface-overlap removal must be capped by safety_limit"
    )
    assert report["surface_overlap_safety_status"] == "exceeded"
    assert report["safety_status"] == "exceeded"  # back-compat alias
    assert report["n_lipids_removed"] == 1
    assert report["coordinate_file_modified"] is True
    assert report["removed_resids"] == [101]
    assert 102 not in report["removed_resids"]

    # Verify the GRO file itself: exactly one lipid residue (the cavity-trapped
    # one) was removed — the capped surface-overlap lipid, plus the two
    # untouched lipids (soluble-adjacent 103, bulk 104), all survive. Residues
    # are renumbered sequentially after removal (gmx editconf -resnr 1), so
    # check by remaining DPP residue count rather than original resid values.
    final_lines = gro_copy.read_text().splitlines()
    initial_atom_count = int(synthetic_system_gro.read_text().splitlines()[1].strip())
    final_atom_count = int(final_lines[1].strip())
    assert final_atom_count == initial_atom_count - 2, "exactly 1 lipid residue (2 atoms) removed"

    dpp_resnums = {int(line[0:5]) for line in final_lines[2:-1] if line[5:10].strip() == "DPP"}
    assert len(dpp_resnums) == 3, "3 of the original 4 DPP lipid residues must remain"

    # Verify report was written to disk with the same fields.
    report_file = tmp_path / "tm_aware_lipid_exclusion_report.json"
    assert report_file.exists()
    disk_report = json.loads(report_file.read_text())
    assert disk_report["surface_overlap_safety_status"] == "exceeded"
    assert disk_report["n_tm_cavity_trapped_removed"] == 1
    assert disk_report["n_tm_surface_overlap_removed"] == 0


def test_exclude_lipids_strict_apply_still_blocks_on_capped_surface(synthetic_system_gro, tmp_path):
    """
    strict_apply must still block when a capped tm_surface_overlap subset
    remains after removal — the cap changes *what* is removed, not whether
    strict_apply validates the final state.
    """
    import shutil
    from validators.tm_exclusion_mask import exclude_lipids_tm_aware

    gro_copy = tmp_path / "system_safety_strict.gro"
    shutil.copy(synthetic_system_gro, gro_copy)
    tm_residues = set(range(1, 25))

    with pytest.raises(ValueError, match="TM-aware strict_apply validation failed"):
        exclude_lipids_tm_aware(
            gro_path=gro_copy,
            tm_residues=tm_residues,
            lipid_resname="DPP",
            policy="strict_apply",
            safety_limit=0,
        )

    report_file = tmp_path / "tm_aware_lipid_exclusion_report.json"
    assert report_file.exists()
    report = json.loads(report_file.read_text())
    # Cavity-trapped was still removed even though strict_apply then blocked
    # on the residual (capped) surface-overlap lipid.
    assert report["n_tm_cavity_trapped_removed"] == 1
    assert report["n_tm_surface_overlap_removed"] == 0
    assert report["post_exclusion_n_tm_surface_overlap"] >= 1


def test_exclude_lipids_strict_apply_blocks(synthetic_system_gro, tmp_path):
    import shutil
    from unittest.mock import patch
    from validators.tm_exclusion_mask import exclude_lipids_tm_aware, classify_lipids_tm_aware

    gro_copy = tmp_path / "system_strict.gro"
    shutil.copy(synthetic_system_gro, gro_copy)
    tm_residues = set(range(1, 25))

    # We want classify_lipids_tm_aware to return surface overlap count > 0 on the post-exclusion call
    original_classify = classify_lipids_tm_aware
    call_count = 0
    def mock_classify(*args, **kwargs):
        nonlocal call_count
        call_count += 1
        res = original_classify(*args, **kwargs)
        if call_count > 1:
            res["classification_report"]["counts"]["tm_surface_overlap"] = 1
        return res

    with patch("validators.tm_exclusion_mask.classify_lipids_tm_aware", side_effect=mock_classify):
        with pytest.raises(ValueError, match="TM-aware strict_apply validation failed"):
            exclude_lipids_tm_aware(
                gro_path=gro_copy,
                tm_residues=tm_residues,
                lipid_resname="DPP",
                policy="strict_apply",
                safety_limit=5
            )

