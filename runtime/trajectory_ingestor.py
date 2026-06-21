"""Phase A — Scientific Interpretation Layer: trajectory ingestor.

Auto-discovers MD output files from any directory structure (SimForge workspaces,
raw GROMACS output, or mixed layouts) and labels XVG files by content type.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from runtime.xvg_parser import parse_xvg, XVGData


# ─────────────────────────────────────────────────────────────────────────────
# TrajectoryManifest
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class TrajectoryManifest:
    """Auto-discovered MD output files from any directory."""
    root:               Path
    xvg_files:          dict[str, Path]   # label → path
    xtc_files:          list[Path]
    edr_files:          list[Path]
    log_files:          list[Path]
    gro_files:          list[Path]
    tpr_files:          list[Path]
    simforge_workspace: bool              # True if looks like a SimForge workspace


# ─────────────────────────────────────────────────────────────────────────────
# XVG labeling heuristics
# ─────────────────────────────────────────────────────────────────────────────

_LABEL_RULES: list[tuple[list[str], str]] = [
    # (substrings to search in stem/title, canonical label)
    # Order matters: first match wins; longer/more-specific patterns listed first.
    (["rmsd"],                                                        "rmsd"),
    (["rmsf"],                                                        "rmsf"),
    (["potential energy", "epot", "potential_energy"],                "potential_energy"),
    (["energy"],                                                      "potential_energy"),
    (["temperature"],                                                 "temperature"),
    (["pressure"],                                                    "pressure"),
    (["sasa", "solvent accessible", "solvent_accessible", "solvacc"], "sasa"),
    (["gyrate", "gyration", "radius of gyration", "radius_of_gyration"], "gyration"),
    (["density"],                                                     "density"),
    (["hbnum", "hydrogen bond", "hydrogen_bond", "hbond", "hbonds"], "hydrogen_bonds"),
    (["distance"],                                                    "distance"),
]

# Short prefixes checked against the file stem only — kept separate because
# single letters / two-letter strings are unsafe as general substrings
# (e.g. "rg" appears inside "energy"; "hb" inside "inhibit").
# These only fire when no _LABEL_RULES rule matched first.
_PREFIX_RULES: list[tuple[str, str]] = [
    ("hb",  "hydrogen_bonds"),  # hb_.xvg, hb-1.xvg, hbnum.xvg, …
    ("rg",  "gyration"),        # rg.xvg, rg_prot.xvg, rg-.xvg, …
]


def _label_xvg(path: Path, title: str) -> str:
    """Determine a semantic label for an XVG file from its stem and title.

    Matching order:
      1. Substring rules (_LABEL_RULES) — checked in declaration order.
      2. Stem-prefix rules (_PREFIX_RULES) — only for short tokens that would
         produce false positives as general substrings.
      3. Fallback: the raw file stem.
    """
    name_lower  = path.stem.lower()
    title_lower = title.lower()

    for keywords, label in _LABEL_RULES:
        for kw in keywords:
            if kw in name_lower or kw in title_lower:
                return label

    for prefix, label in _PREFIX_RULES:
        if name_lower.startswith(prefix):
            return label

    # Fallback: use the file stem as-is
    return path.stem


def _peek_title(path: Path) -> str:
    """Read only the title directive from an XVG file (fast, no full parse)."""
    try:
        with path.open(errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if not line.startswith("@"):
                    break
                if "title" in line.lower() and '"' in line:
                    parts = line.split('"')
                    if len(parts) >= 2:
                        return parts[1]
    except OSError:
        pass
    return ""


# ─────────────────────────────────────────────────────────────────────────────
# SimForge workspace detection
# ─────────────────────────────────────────────────────────────────────────────

def _is_simforge_workspace(path: Path) -> bool:
    """Heuristic: SimForge workspace has steps/ and metadata/ subdirectories."""
    return (path / "steps").exists() and (path / "metadata").exists()


# ─────────────────────────────────────────────────────────────────────────────
# discover_trajectory — public API
# ─────────────────────────────────────────────────────────────────────────────

def discover_trajectory(path: Path) -> TrajectoryManifest:
    """Walk a directory and classify all MD output files.

    Works with:
      - SimForge workspaces (steps/ subdirs)
      - Raw GROMACS output directories (flat or nested)
      - Any mix of the above

    For XVG files with duplicate labels, later-found files overwrite earlier
    ones; the last file discovered is used (typically the most complete run).
    To keep all, callers can inspect manifest.xvg_files directly.
    """
    if not path.exists() or not path.is_dir():
        return TrajectoryManifest(
            root               = path,
            xvg_files          = {},
            xtc_files          = [],
            edr_files          = [],
            log_files          = [],
            gro_files          = [],
            tpr_files          = [],
            simforge_workspace = False,
        )

    is_sf = _is_simforge_workspace(path)

    xvg_files: dict[str, Path] = {}
    xtc_files: list[Path]      = []
    edr_files: list[Path]      = []
    log_files: list[Path]      = []
    gro_files: list[Path]      = []
    tpr_files: list[Path]      = []

    # Walk the entire tree (sorted for determinism)
    for p in sorted(path.rglob("*")):
        if not p.is_file():
            continue

        suffix = p.suffix.lower()

        if suffix == ".xvg":
            title = _peek_title(p)
            label = _label_xvg(p, title)
            # If label already exists, make a unique key
            if label in xvg_files:
                # Keep the one with more data (larger file)
                existing = xvg_files[label]
                if p.stat().st_size > existing.stat().st_size:
                    xvg_files[label] = p
            else:
                xvg_files[label] = p

        elif suffix == ".xtc":
            xtc_files.append(p)
        elif suffix == ".edr":
            edr_files.append(p)
        elif suffix == ".log":
            log_files.append(p)
        elif suffix == ".gro":
            gro_files.append(p)
        elif suffix == ".tpr":
            tpr_files.append(p)

    return TrajectoryManifest(
        root               = path,
        xvg_files          = xvg_files,
        xtc_files          = xtc_files,
        edr_files          = edr_files,
        log_files          = log_files,
        gro_files          = gro_files,
        tpr_files          = tpr_files,
        simforge_workspace = is_sf,
    )


# ─────────────────────────────────────────────────────────────────────────────
# load_xvg_files — public API
# ─────────────────────────────────────────────────────────────────────────────

def load_xvg_files(manifest: TrajectoryManifest) -> dict[str, XVGData]:
    """Parse all discovered XVG files. Skip files that fail to parse."""
    result: dict[str, XVGData] = {}
    for label, path in manifest.xvg_files.items():
        try:
            data = parse_xvg(path)
            result[label] = data
        except Exception:
            pass   # skip unreadable files — degrade gracefully
    return result
