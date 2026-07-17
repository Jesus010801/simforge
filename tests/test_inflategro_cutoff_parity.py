"""
tests/test_inflategro_cutoff_parity.py

Regression tests for docs/audits/simforge_vs_protmemfiles_audit.md Fix 1.

These tests invoke the real, vendored inflategro-Jorge.pl script directly —
the same file `builders/workspace_builder.py` stages into every generated
workspace's membrane_assets/ (see `_PROT_MEMB_FILES` in that module) and the
same file the generated embedding step's shrink loop calls as `$INFLATEGRO`.
Using the original script itself as the oracle (rather than reimplementing
its overlap-removal logic in Python) is the most direct possible proof that
SimForge's default cutoff now reproduces the original protmemfiles recipe.
"""
from __future__ import annotations

import math
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

_HAS_PERL = bool(shutil.which("perl"))
_INFLATEGRO_SCRIPT = Path("docs/Prot-Memb_FILES/inflategro-Jorge.pl")


def _make_gro_line(resnum: int, resname: str, atomname: str, atomnum: int,
                    x: float, y: float, z: float) -> str:
    return f"{resnum:>5}{resname:<5}{atomname:>5}{atomnum:>5}{x:8.3f}{y:8.3f}{z:8.3f}"


def _write_synthetic_system(path: Path) -> dict[str, float]:
    """
    Build a synthetic protein-ring + lipid system where lipids sit at known
    radial gaps (in nm) from the protein Cα ring surface, spanning both
    leaflets:

      gap 0.05 nm  -> inside both the old (0.14 nm) and new (1.4 nm) cutoff
      gap 0.30 nm  -> BETWEEN the old and new cutoff (the regression target)
      gap 0.80 nm  -> BETWEEN the old and new cutoff (the regression target)
      gap 2.00 nm  -> outside both cutoffs (bulk lipid, must always survive)

    Returns a dict mapping resid -> gap_nm for reference in assertions.
    """
    atoms: list[str] = []
    n = 1

    # Protein: ring of 8 CA atoms, radius 1.0 nm, centered at (5, 5),
    # replicated at both leaflet Z heights so both leaflets clash-check.
    for i in range(8):
        theta = i * (2 * math.pi / 8)
        x = 5.0 + 1.0 * math.cos(theta)
        y = 5.0 + 1.0 * math.sin(theta)
        for z, resbase in [(3.5, 1), (6.5, 9)]:
            atoms.append(_make_gro_line(resbase + i, "PRO", "CA", n, x, y, z))
            n += 1

    gaps = {
        101: (0.05, 6.5),
        102: (0.30, 6.5),
        103: (0.80, 6.5),
        104: (2.00, 6.5),
        105: (0.05, 3.5),
        106: (0.30, 3.5),
        107: (0.80, 3.5),
        108: (2.00, 3.5),
    }
    for resid, (gap, z) in gaps.items():
        r = 1.0 + gap
        x = 5.0 + r
        y = 5.0
        atoms.append(_make_gro_line(resid, "DPP", "P", n, x, y, z))
        n += 1
        atoms.append(_make_gro_line(resid, "DPP", "C1", n, x, y, z - 0.2 if z > 5 else z + 0.2))
        n += 1

    lines = ["synthetic clash-removal system", str(len(atoms))]
    lines.extend(atoms)
    lines.append(f"{20.0:10.5f}{20.0:10.5f}{20.0:10.5f}")
    path.write_text("\n".join(lines) + "\n")

    return {resid: gap for resid, (gap, _z) in gaps.items()}


def _remaining_dpp_resids(gro_path: Path) -> set[int]:
    lines = gro_path.read_text().splitlines()
    n = int(lines[1].strip())
    resids = set()
    for line in lines[2: 2 + n]:
        if line[5:10].strip() == "DPP":
            resids.add(int(line[0:5]))
    return resids


@pytest.mark.skipif(not _HAS_PERL, reason="perl not found on PATH")
@pytest.mark.skipif(not _INFLATEGRO_SCRIPT.exists(), reason="docs/Prot-Memb_FILES/inflategro-Jorge.pl not found")
class TestInflategroCutoffParity:

    def test_corrected_cutoff_removes_lipids_between_old_and_new_default(self, tmp_path):
        """
        At the corrected default (cutoff_nm=1.4 -> Perl arg "14"), every lipid
        within 1.4 nm of the protein ring — including those between the old
        buggy default (0.14 nm) and the corrected one — must be removed.
        Only the 2.0 nm bulk lipids survive.
        """
        gro_in = tmp_path / "system.gro"
        gaps = _write_synthetic_system(gro_in)

        script = _INFLATEGRO_SCRIPT.resolve()
        out_gro = tmp_path / "system_inflated.gro"
        ret = subprocess.run(
            ["perl", str(script), "system.gro", "1", "DPP", "14",
             "system_inflated.gro", "5", "area.dat"],
            cwd=tmp_path, capture_output=True, text=True, timeout=60,
        )
        assert ret.returncode == 0, f"inflategro-Jorge.pl failed:\n{ret.stdout}\n{ret.stderr}"

        remaining = _remaining_dpp_resids(out_gro)

        for resid, gap in gaps.items():
            if gap < 1.4:
                assert resid not in remaining, (
                    f"lipid {resid} at gap={gap} nm (< 1.4 nm cutoff) must be removed "
                    f"by the corrected default; remaining DPP resids: {sorted(remaining)}"
                )
            else:
                assert resid in remaining, (
                    f"lipid {resid} at gap={gap} nm (>= 1.4 nm cutoff, bulk) must survive; "
                    f"remaining DPP resids: {sorted(remaining)}"
                )

    def test_old_default_leaves_mid_range_lipids_trapped(self, tmp_path):
        """
        Demonstrates the bug the fix corrects: at the previous default
        (cutoff_nm=0.14 -> Perl arg "1.4"), lipids at 0.30 nm and 0.80 nm —
        squarely between the old and new cutoff — are NOT removed and remain
        clashing with the protein.
        """
        gro_in = tmp_path / "system.gro"
        gaps = _write_synthetic_system(gro_in)

        script = _INFLATEGRO_SCRIPT.resolve()
        out_gro = tmp_path / "system_inflated.gro"
        ret = subprocess.run(
            ["perl", str(script), "system.gro", "1", "DPP", "1.4",
             "system_inflated.gro", "5", "area.dat"],
            cwd=tmp_path, capture_output=True, text=True, timeout=60,
        )
        assert ret.returncode == 0, f"inflategro-Jorge.pl failed:\n{ret.stdout}\n{ret.stderr}"

        remaining = _remaining_dpp_resids(out_gro)

        # gap=0.05 nm: removed even under the old (weaker) cutoff.
        assert 101 not in remaining
        assert 105 not in remaining
        # gap=0.30 nm and 0.80 nm: the regression — old cutoff leaves these
        # trapped against the protein.
        for resid in (102, 103, 106, 107):
            assert resid in remaining, (
                f"lipid {resid} unexpectedly removed at the old 0.14 nm cutoff — "
                "fixture assumption violated"
            )


def _write_deflate_probe_system(path: Path, gap_nm: float = 2.0) -> None:
    """
    Build a synthetic protein-ring + single-lipid-per-leaflet system where
    each lipid starts safely outside the one-shot cutoff (gap_nm from the
    ring surface), positioned so that repeated radial deflation (scale 0.95)
    compresses it back toward the protein over successive iterations —
    reproducing the shrink loop's physical behavior.
    """
    atoms: list[str] = []
    n = 1
    for i in range(8):
        theta = i * (2 * math.pi / 8)
        x = 5.0 + 1.0 * math.cos(theta)
        y = 5.0 + 1.0 * math.sin(theta)
        for z, resbase in [(3.5, 1), (6.5, 9)]:
            atoms.append(_make_gro_line(resbase + i, "PRO", "CA", n, x, y, z))
            n += 1

    r = 1.0 + gap_nm
    for resid, z in [(101, 6.5), (102, 3.5)]:
        x = 5.0 + r
        y = 5.0
        atoms.append(_make_gro_line(resid, "DPP", "P", n, x, y, z))
        n += 1
        atoms.append(_make_gro_line(resid, "DPP", "C1", n, x, y, z - 0.2 if z > 5 else z + 0.2))
        n += 1

    lines = ["deflate-loop probe system", str(len(atoms))]
    lines.extend(atoms)
    lines.append(f"{20.0:10.5f}{20.0:10.5f}{20.0:10.5f}")
    path.write_text("\n".join(lines) + "\n")


@pytest.mark.skipif(not _HAS_PERL, reason="perl not found on PATH")
@pytest.mark.skipif(not _INFLATEGRO_SCRIPT.exists(), reason="docs/Prot-Memb_FILES/inflategro-Jorge.pl not found")
class TestShrinkLoopCutoffParity:
    """
    Regression for the GLP-1R membrane-void bug: reusing the one-shot cutoff
    on every shrink-loop deflate iteration cumulatively deletes lipids that
    deflation itself compresses back within range, even though the protein
    never moved. See docs/audits/glp1r_membrane_void_shrink_loop_audit.md.

    Uses the real vendored inflategro-Jorge.pl as ground truth, simulating
    the shrink loop's repeated deflate calls directly.
    """

    def test_reusing_one_shot_cutoff_in_loop_deletes_the_lipid(self, tmp_path):
        """
        Reproduces the bug: cutoff=14 (1.4 nm) reapplied on every deflate
        iteration removes a lipid that started safely outside the cutoff.

        Once both lipids are gone, inflategro-Jorge.pl itself divides by zero
        computing per-leaflet lipid counts (there is nothing left to divide
        by) — that crash is itself confirmation that lipid depletion reached
        zero, since the script writes its output file only after that point,
        so loop.gro is left holding the last successfully-written (already
        fully depleted) state.
        """
        gro = tmp_path / "loop.gro"
        _write_deflate_probe_system(gro, gap_nm=2.0)

        script = _INFLATEGRO_SCRIPT.resolve()
        depleted_early = False
        for _ in range(10):
            ret = subprocess.run(
                ["perl", str(script), "loop.gro", "0.95", "DPP", "14",
                 "loop.gro", "5", "area.dat"],
                cwd=tmp_path, capture_output=True, text=True, timeout=60,
            )
            if ret.returncode != 0:
                assert "division by zero" in ret.stderr.lower(), (
                    f"unexpected inflategro-Jorge.pl failure:\n{ret.stdout}\n{ret.stderr}"
                )
                depleted_early = True
                break

        remaining = _remaining_dpp_resids(gro)
        assert not remaining, (
            "fixture assumption: reusing the one-shot cutoff every deflate "
            f"iteration should delete both lipids within 10 iterations; "
            f"remaining: {sorted(remaining)} (depleted_early={depleted_early})"
        )

    def test_zero_cutoff_in_loop_preserves_the_lipid(self, tmp_path):
        """Proves the fix: cutoff=0 during every deflate iteration (matching
        the original protmemfiles method) leaves the lipid untouched through
        the same 10 iterations, even though the box is compressing."""
        gro = tmp_path / "loop.gro"
        _write_deflate_probe_system(gro, gap_nm=2.0)

        script = _INFLATEGRO_SCRIPT.resolve()
        for _ in range(10):
            ret = subprocess.run(
                ["perl", str(script), "loop.gro", "0.95", "DPP", "0",
                 "loop.gro", "5", "area.dat"],
                cwd=tmp_path, capture_output=True, text=True, timeout=60,
            )
            assert ret.returncode == 0, f"inflategro-Jorge.pl failed:\n{ret.stdout}\n{ret.stderr}"

        remaining = _remaining_dpp_resids(gro)
        assert remaining == {101, 102}, (
            "with cutoff=0 (the fix), both lipids must survive every deflate "
            f"iteration untouched; remaining: {sorted(remaining)}"
        )
