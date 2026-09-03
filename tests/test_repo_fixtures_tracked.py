"""Guard against 'removed a file but left a reference' regressions.

A file can exist on the maintainer's disk (untracked) and make the suite pass
locally while breaking a fresh clone / CI. These tests fail if a config or a
pipeline asset that the test suite depends on is not actually tracked by git.

Skipped entirely when not run from a git working tree (e.g. from an sdist).
"""
from __future__ import annotations

import subprocess
from functools import lru_cache
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parent.parent


@lru_cache(maxsize=1)
def _tracked_files() -> frozenset[str] | None:
    if not (_ROOT / ".git").exists():
        return None
    try:
        out = subprocess.run(
            ["git", "-C", str(_ROOT), "ls-files"],
            capture_output=True, text=True, timeout=30,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if out.returncode != 0:
        return None
    return frozenset(out.stdout.split())


@pytest.fixture(scope="session")
def tracked() -> frozenset[str]:
    t = _tracked_files()
    if t is None:
        pytest.skip("not a git working tree — cannot check tracked files")
    return t


def _rel(p: Path) -> str:
    return p.resolve().relative_to(_ROOT).as_posix()


# ── config `file:` references ───────────────────────────────────────────────

def _config_yamls() -> list[Path]:
    # configs/systems/*.yaml are documentation templates (protein_only.* is a
    # placeholder), not loaded by any test — not checked here.
    return sorted((_ROOT / "configs").glob("*.yaml"))


@pytest.mark.parametrize("cfg", _config_yamls(), ids=lambda p: p.name)
def test_config_referenced_inputs_are_tracked(cfg, tracked):
    if _rel(cfg) not in tracked:
        pytest.skip(f"{cfg.name} is not tracked (local-only config)")

    import yaml

    doc = yaml.safe_load(cfg.read_text()) or {}
    for comp in (doc.get("components") or []):
        ref = comp.get("file")
        if not ref:
            continue
        cand = (cfg.parent / ref).resolve()
        if not cand.exists():
            cand = (_ROOT / ref).resolve()  # some configs use repo-rel paths
        assert cand.exists(), f"{cfg.name} references missing file: {ref}"
        assert _rel(cand) in tracked, (
            f"{cfg.name} references {ref} which exists on disk but is NOT "
            f"tracked by git — a fresh clone / CI would fail to build it"
        )


# ── membrane pipeline assets ───────────────────────────────────────────────

_MEMBRANE_ASSETS = [
    "docs/Prot-Memb_FILES/dppc512_whole.gro",
    "docs/Prot-Memb_FILES/dppc128.gro",
    "docs/Prot-Memb_FILES/inflategro-Jorge.pl",
    "docs/Prot-Memb_FILES/oplsaa_membrane.ff/forcefield.itp",
    "docs/Prot-Memb_FILES/oplsaa_membrane.ff/ffnonbonded.itp",
    "docs/Prot-Memb_FILES/oplsaa_membrane.ff/ffbonded.itp",
    "docs/Prot-Memb_FILES/oplsaa_membrane.ff/dpp.itp",
]


@pytest.mark.parametrize("asset", _MEMBRANE_ASSETS)
def test_membrane_pipeline_asset_is_tracked(asset, tracked):
    assert (_ROOT / asset).exists(), f"membrane asset missing on disk: {asset}"
    assert asset in tracked, (
        f"membrane asset {asset} is on disk but NOT tracked by git — the "
        f"membrane test suite would error on a fresh clone"
    )


def test_oplsaa_membrane_ff_is_complete(tracked):
    """A GROMACS .ff directory needs its .itp payload, not just .rtp/.atp."""
    ff = _ROOT / "docs/Prot-Memb_FILES/oplsaa_membrane.ff"
    itps = {p.name for p in ff.glob("*.itp")}
    required = {"forcefield.itp", "ffnonbonded.itp", "ffbonded.itp"}
    missing = required - itps
    assert not missing, f"oplsaa_membrane.ff is missing {missing}"
    for name in required:
        assert _rel(ff / name) in tracked, f"{name} not tracked"
