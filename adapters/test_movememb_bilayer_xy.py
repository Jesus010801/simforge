# adapters/test_movememb_bilayer_xy.py
"""
Phase 3 tests: bilayer XY authority in MoveMembAdapter.

The selected bilayer is authoritative for simulation-cell XY dimensions.
Output system.gro must use bilayer box X/Y; protein GRO provides Z.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from adapters.movememb_adapter import MoveMembAdapter, _parse_gro


# ── GRO construction helpers (mirrors test_movememb_tm_alignment.py) ──────────

def _gro_line(resnum: int, resname: str, atomname: str, serial: int,
              x: float, y: float, z: float) -> str:
    return f"{resnum:5d}{resname:<5}{atomname:>5}{serial:5d}{x:8.3f}{y:8.3f}{z:8.3f}"


def _write_gro(path: Path, title: str, atoms: list[tuple],
               box: tuple[float, float, float] = (10.0, 10.0, 20.0)) -> Path:
    lines = [title, str(len(atoms))]
    for i, (rn, rname, aname, x, y, z) in enumerate(atoms, start=1):
        lines.append(_gro_line(rn, rname, aname, i, x, y, z))
    lines.append(f"{box[0]:.5f}  {box[1]:.5f}  {box[2]:.5f}")
    path.write_text("\n".join(lines) + "\n")
    return path


def _protein_gro(tmp_path: Path, prot_box=(10.0, 10.0, 15.0)) -> Path:
    """
    Small protein with 3×3 nm XY footprint, 5 nm Z.
    Box is passed in explicitly to simulate what match_box_to_bilayer produces.
    """
    atoms = [
        (1, "ALA", "CA", 3.5, 3.5, 5.0),
        (2, "ALA", "CA", 6.5, 6.5, 5.0),
        (3, "ALA", "CA", 5.0, 5.0, 7.5),
    ]
    return _write_gro(tmp_path / "protein_boxed.gro", "boxed protein", atoms, box=prot_box)


def _bilayer_gro(tmp_path: Path, bil_box=(12.84, 12.89, 12.0),
                 midplane_z=6.0) -> Path:
    """Minimal bilayer with given box XY and midplane at midplane_z."""
    half = 3.0
    atoms = [
        (1, "DPP", "O33", 6.0, 6.0, midplane_z - half),
        (2, "DPP", "O33", 6.0, 6.0, midplane_z + half),
    ]
    return _write_gro(tmp_path / "bilayer.gro", "DPPC bilayer", atoms, box=bil_box)


# ── Box vector tests ──────────────────────────────────────────────────────────

class TestMoveMembAdapterBilayerXY:
    """Bilayer XY is authoritative; protein GRO provides Z."""

    def test_output_box_xy_equals_bilayer_box_xy(self, tmp_path):
        prot = _protein_gro(tmp_path, prot_box=(10.0, 10.0, 15.0))
        bil  = _bilayer_gro(tmp_path, bil_box=(12.84, 12.89, 12.0))
        out  = tmp_path / "system.gro"

        result = MoveMembAdapter().run(
            protein_gro=prot, bilayer_gro=bil, gro_out=out,
        )
        assert result.success, result.error_message

        _, _, box = _parse_gro(out)
        assert abs(box[0] - 12.84) < 1e-4, f"Expected box X=12.84, got {box[0]}"
        assert abs(box[1] - 12.89) < 1e-4, f"Expected box Y=12.89, got {box[1]}"

    def test_output_box_z_equals_protein_gro_z(self, tmp_path):
        prot = _protein_gro(tmp_path, prot_box=(10.0, 10.0, 15.0))
        bil  = _bilayer_gro(tmp_path, bil_box=(12.84, 12.89, 12.0))
        out  = tmp_path / "system.gro"

        MoveMembAdapter().run(protein_gro=prot, bilayer_gro=bil, gro_out=out)
        _, _, box = _parse_gro(out)
        assert abs(box[2] - 15.0) < 1e-4, f"Expected box Z=15.0 (from protein GRO), got {box[2]}"

    def test_protein_box_xy_is_not_used_for_output(self, tmp_path):
        """Box XY must be bilayer's, not protein GRO's."""
        prot = _protein_gro(tmp_path, prot_box=(7.5, 8.0, 15.0))
        bil  = _bilayer_gro(tmp_path, bil_box=(12.84, 12.89, 12.0))
        out  = tmp_path / "system.gro"

        MoveMembAdapter().run(protein_gro=prot, bilayer_gro=bil, gro_out=out)
        _, _, box = _parse_gro(out)
        assert abs(box[0] - 7.5) > 0.1, "Output box X must NOT come from protein GRO"
        assert abs(box[1] - 8.0) > 0.1, "Output box Y must NOT come from protein GRO"
        assert abs(box[0] - 12.84) < 1e-4
        assert abs(box[1] - 12.89) < 1e-4

    def test_metadata_xy_authority_is_selected_bilayer(self, tmp_path):
        prot = _protein_gro(tmp_path)
        bil  = _bilayer_gro(tmp_path)
        out  = tmp_path / "system.gro"

        m = MoveMembAdapter().run(
            protein_gro=prot, bilayer_gro=bil, gro_out=out,
        ).metadata
        assert m["xy_authority"] == "selected_bilayer"

    def test_metadata_bilayer_xy_fields_present(self, tmp_path):
        prot = _protein_gro(tmp_path)
        bil  = _bilayer_gro(tmp_path, bil_box=(12.84, 12.89, 12.0))
        out  = tmp_path / "system.gro"

        m = MoveMembAdapter().run(
            protein_gro=prot, bilayer_gro=bil, gro_out=out,
        ).metadata

        assert abs(m["bilayer_x"] - 12.84) < 1e-4
        assert abs(m["bilayer_y"] - 12.89) < 1e-4
        assert "protein_x" in m
        assert "protein_y" in m
        assert "required_margin_xy" in m
        assert "available_margin_x" in m
        assert "available_margin_y" in m
        assert "xy_fit_status" in m

    def test_protein_smaller_than_bilayer_passes(self, tmp_path):
        """Protein 3×3 nm in 12.84×12.89 nm bilayer → plenty of margin → pass."""
        prot = _protein_gro(tmp_path)
        bil  = _bilayer_gro(tmp_path, bil_box=(12.84, 12.89, 12.0))
        out  = tmp_path / "system.gro"

        result = MoveMembAdapter().run(
            protein_gro=prot, bilayer_gro=bil, gro_out=out, required_margin_nm=1.0,
        )
        assert result.success
        assert result.metadata["xy_fit_status"] == "pass"
        assert result.metadata["available_margin_x"] > 1.0
        assert result.metadata["available_margin_y"] > 1.0

    def test_protein_larger_than_bilayer_recorded_as_fail(self, tmp_path):
        """
        Protein footprint larger than bilayer → xy_fit_status=fail in metadata.
        The adapter records the failure but does NOT raise/exit — the gate
        (box_match_helper) is responsible for blocking the pipeline.
        """
        # Large protein: 10×10 nm, bilayer only 8×8 nm, required margin 1 nm
        atoms = [
            (1, "ALA", "CA", 0.5, 0.5, 5.0),
            (2, "ALA", "CA", 10.5, 10.5, 5.0),
        ]
        prot = _write_gro(
            tmp_path / "prot_large.gro", "large protein", atoms, box=(12.0, 12.0, 15.0)
        )
        bil = _bilayer_gro(tmp_path, bil_box=(8.0, 8.0, 8.0))
        out = tmp_path / "system.gro"

        result = MoveMembAdapter().run(
            protein_gro=prot, bilayer_gro=bil, gro_out=out, required_margin_nm=1.0,
        )
        assert result.success  # adapter still writes the file — gate blocks downstream
        assert result.metadata["xy_fit_status"] == "fail"
        assert result.metadata["available_margin_x"] < 1.0

    def test_metadata_required_margin_default_is_one(self, tmp_path):
        prot = _protein_gro(tmp_path)
        bil  = _bilayer_gro(tmp_path)
        out  = tmp_path / "system.gro"

        m = MoveMembAdapter().run(
            protein_gro=prot, bilayer_gro=bil, gro_out=out,
        ).metadata
        assert m["required_margin_xy"] == pytest.approx(1.0)

    def test_custom_required_margin_stored_in_metadata(self, tmp_path):
        prot = _protein_gro(tmp_path)
        bil  = _bilayer_gro(tmp_path)
        out  = tmp_path / "system.gro"

        m = MoveMembAdapter().run(
            protein_gro=prot, bilayer_gro=bil, gro_out=out, required_margin_nm=2.0,
        ).metadata
        assert m["required_margin_xy"] == pytest.approx(2.0)

    def test_different_bilayers_produce_different_box_xy(self, tmp_path):
        """Two bilayer patches of different sizes produce different output box XY."""
        prot = _protein_gro(tmp_path)
        out1 = tmp_path / "sys1.gro"
        out2 = tmp_path / "sys2.gro"
        bil1 = _write_gro(
            tmp_path / "bilayer1.gro", "large bilayer",
            [(1, "DPP", "O33", 6.0, 6.0, 3.0), (2, "DPP", "O33", 6.0, 6.0, 9.0)],
            box=(12.84, 12.89, 12.0),
        )
        bil2 = _write_gro(
            tmp_path / "bilayer2.gro", "small bilayer",
            [(1, "DPP", "O33", 4.0, 4.0, 3.0), (2, "DPP", "O33", 4.0, 4.0, 6.0)],
            box=(8.0, 8.5, 8.0),
        )

        MoveMembAdapter().run(protein_gro=prot, bilayer_gro=bil1, gro_out=out1)
        MoveMembAdapter().run(protein_gro=prot, bilayer_gro=bil2, gro_out=out2)

        _, _, box1 = _parse_gro(out1)
        _, _, box2 = _parse_gro(out2)
        assert abs(box1[0] - box2[0]) > 1.0, "Different bilayers must produce different box X"
        assert abs(box1[1] - box2[1]) > 1.0, "Different bilayers must produce different box Y"

    def test_z_field_present_and_not_from_bilayer(self, tmp_path):
        """Z must come from protein GRO, not from bilayer GRO box Z."""
        prot = _protein_gro(tmp_path, prot_box=(10.0, 10.0, 18.5))
        bil  = _bilayer_gro(tmp_path, bil_box=(12.84, 12.89, 9.0))  # Z=9, not 18.5
        out  = tmp_path / "system.gro"

        MoveMembAdapter().run(protein_gro=prot, bilayer_gro=bil, gro_out=out)
        _, _, box = _parse_gro(out)
        # box Z must be protein's 18.5, not bilayer's 9.0
        assert abs(box[2] - 18.5) < 1e-4
        assert abs(box[2] - 9.0)  > 1.0
