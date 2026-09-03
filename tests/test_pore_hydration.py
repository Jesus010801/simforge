from __future__ import annotations

import json
import math
from pathlib import Path

from validators.pore_hydration import (
    classify_pore_hydration,
    clean_water_channel_aware,
)


def _write_gro(path: Path, atoms: list[tuple], box=(12.0, 12.0, 10.0)) -> None:
    lines = ["pore hydration test", str(len(atoms))]
    for atom_id, (resid, resname, atomname, x, y, z) in enumerate(atoms, 1):
        lines.append(
            f"{resid:5d}{resname:<5s}{atomname:>5s}{atom_id:5d}"
            f"{x:8.3f}{y:8.3f}{z:8.3f}"
        )
    lines.append("".join(f"{value:10.5f}" for value in box))
    path.write_text("\n".join(lines) + "\n")


def _water(resid: int, x: float, y: float, z: float) -> list[tuple]:
    return [
        (resid, "SOL", "OW", x, y, z),
        (resid, "SOL", "HW1", x + 0.01, y, z),
        (resid, "SOL", "HW2", x, y + 0.01, z),
    ]


def _channel(path: Path, include_clash: bool = False) -> tuple[Path, set[int]]:
    atoms: list[tuple] = []
    tm_residues: set[int] = set()
    resid = 1
    for z in (3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5):
        for index in range(20):
            angle = index * 2 * math.pi / 20
            atoms.append((
                resid, "ALA", "CA",
                6.0 + math.cos(angle), 6.0 + math.sin(angle), z,
            ))
            tm_residues.add(resid)
            resid += 1
    for z_head, z_tail in ((6.0, 5.5), (4.0, 4.5)):
        for index in range(12):
            angle = index * 2 * math.pi / 12
            x, y = 6.0 + 3.0 * math.cos(angle), 6.0 + 3.0 * math.sin(angle)
            atoms.extend([
                (resid, "DPP", "P", x, y, z_head),
                (resid, "DPP", "C50", x, y, z_tail),
            ])
            resid += 1
    atoms += _water(resid, 6.0, 6.0, 5.0)  # lumen
    resid += 1
    atoms += _water(resid, 9.0, 6.0, 5.0)  # external lipid core
    resid += 1
    atoms += _water(resid, 6.0, 6.0, 8.0)  # bulk
    resid += 1
    if include_clash:
        atoms += _water(resid, 7.0, 6.0, 5.0)  # on protein ring
    _write_gro(path, atoms)
    return path, tm_residues


def test_channel_lumen_water_preserved_and_external_core_removed(tmp_path):
    gro, tm = _channel(tmp_path / "channel.gro")
    report = classify_pore_hydration(gro, tm_residues=tm, output_dir=tmp_path)
    assert report["channel_like_system_detected"] is True
    assert report["n_channel_lumen_water"] >= 1
    assert report["n_lipid_core_outside_protein"] >= 1
    assert report["n_bulk_water"] >= 1

    out = tmp_path / "clean.gro"
    clean = clean_water_channel_aware(
        gro, out, tm_residues=tm, output_dir=tmp_path
    )
    assert clean["n_water_molecules_removed"] == 1
    assert out.read_text().count(" OW") == 2


def test_clashing_water_is_flagged(tmp_path):
    gro, tm = _channel(tmp_path / "channel.gro", include_clash=True)
    report = classify_pore_hydration(gro, tm_residues=tm, output_dir=tmp_path)
    assert report["n_clash_water"] >= 1


def test_compact_gpcr_like_system_does_not_force_pore_hydration(tmp_path):
    atoms: list[tuple] = []
    tm = set()
    resid = 1
    for z in (4.0, 5.0, 6.0):
        for radius in (0.0, 0.35, 0.7):
            for index in range(12):
                angle = index * 2 * math.pi / 12
                atoms.append((
                    resid, "ALA", "CA",
                    6 + radius * math.cos(angle), 6 + radius * math.sin(angle), z,
                ))
                tm.add(resid)
                resid += 1
    for z_head, z_tail in ((6.0, 5.5), (4.0, 4.5)):
        atoms.extend([
            (resid, "DPP", "P", 9.0, 6.0, z_head),
            (resid, "DPP", "C50", 9.0, 6.0, z_tail),
        ])
        resid += 1
    _write_gro(tmp_path / "gpcr.gro", atoms)
    report = classify_pore_hydration(
        tmp_path / "gpcr.gro", tm_residues=tm, output_dir=tmp_path
    )
    assert report["channel_like_system_detected"] is False
    assert report["n_channel_lumen_water"] == 0


def test_report_always_written_when_disabled(tmp_path):
    gro, tm = _channel(tmp_path / "channel.gro")
    report = classify_pore_hydration(
        gro, tm_residues=tm, output_dir=tmp_path, enabled=False
    )
    saved = json.loads((tmp_path / "pore_hydration_report.json").read_text())
    assert report["enabled"] is False
    assert saved["coordinate_file_modified"] is False


def test_channel_cleanup_updates_topology_and_reports(tmp_path):
    gro, tm = _channel(tmp_path / "channel.gro")
    topol_in = tmp_path / "topol.in"
    topol_out = tmp_path / "topol.top"
    topol_in.write_text("[ molecules ]\nDPP 24\nSOL 3\n")
    clean_water_channel_aware(
        gro,
        tmp_path / "clean.gro",
        tm_residues=tm,
        output_dir=tmp_path,
        topol_in=topol_in,
        topol_out=topol_out,
    )
    assert "SOL              2" in topol_out.read_text()
    clean_report = json.loads((tmp_path / "clean_water_report.json").read_text())
    assert clean_report["pore_aware"] is True
    assert clean_report["n_pore_waters_preserved"] >= 1
