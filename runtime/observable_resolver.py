"""ObservableResolver — naming-convention-agnostic XVG classification.

Accepts heterogeneous filename fragments (rmsd_protein, protein_rmsd,
rmsdProtein, rmsd-ligand, contacts_lig_active, mindist, …) and maps them
to canonical observable names, display labels, units, and semantic groups.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Callable


@dataclass
class ResolvedObservable:
    canonical: str    # e.g. "protein_rmsd"
    display:   str    # e.g. "Protein RMSD"
    units:     str    # e.g. "nm"
    group:     str    # structural | energetic | interaction | other


# ─── Normalization ────────────────────────────────────────────────────────────

_CAMEL_RE = re.compile(r'([a-z])([A-Z])')


def _normalize(s: str) -> str:
    """Lowercase, split camelCase, collapse separators to underscore."""
    s = _CAMEL_RE.sub(r'\1_\2', s)
    return re.sub(r'[-\s]+', '_', s).lower().strip('_')


# ─── Rule table ──────────────────────────────────────────────────────────────
# Each entry: (canonical, display, units, group, match_predicate)
# Rules are evaluated in order; first match wins.

_LIG_KEYWORDS = ("lig", "drug", "substrate", "small_mol", "inhibitor", "ligand")

_RULES: list[tuple[str, str, str, str, Callable[[str], bool]]] = [
    # Ligand RMSD — must come before protein_rmsd so "rmsd_lig" doesn't fall through
    ("ligand_rmsd", "Ligand RMSD", "nm", "structural",
     lambda s: "rmsd" in s and any(k in s for k in _LIG_KEYWORDS)),

    # Protein RMSD — anything with rmsd that isn't ligand-specific
    ("protein_rmsd", "Protein RMSD", "nm", "structural",
     lambda s: "rmsd" in s),

    # RMSF
    ("rmsf", "RMSF", "nm", "structural",
     lambda s: "rmsf" in s),

    # Minimum distance — before generic distance
    ("mindist", "Min Distance", "nm", "interaction",
     lambda s: "mindist" in s or "min_dist" in s or "minimum_dist" in s or "min_distance" in s),

    # Catalytic distance — before generic distance
    ("catalytic_distance", "Catalytic Distance", "nm", "interaction",
     lambda s: "catalytic" in s or "dist_cat" in s or "cat_dist" in s or "active_dist" in s),

    # Contacts
    ("contacts", "Contacts", "", "interaction",
     lambda s: "contact" in s),

    # Radius of gyration
    ("radius_of_gyration", "Radius of Gyration", "nm", "structural",
     lambda s: "gyration" in s or bool(re.search(r'\brg\b', s)) or "gyr" in s),

    # Hydrogen bonds
    ("hydrogen_bonds", "Hydrogen Bonds", "", "interaction",
     lambda s: "hbond" in s or "h_bond" in s or "hydrogen" in s or "hb_" in s),

    # Kinetic energy — before generic energy
    ("kinetic_energy", "Kinetic Energy", "kJ/mol", "energetic",
     lambda s: "kinetic" in s or "ekin" in s),

    # Total energy — before generic energy
    ("total_energy", "Total Energy", "kJ/mol", "energetic",
     lambda s: "total_energy" in s or "etot" in s),

    # Potential energy
    ("potential_energy", "Potential Energy", "kJ/mol", "energetic",
     lambda s: "potential" in s or "epot" in s or "energy" in s),

    # Temperature
    ("temperature", "Temperature", "K", "energetic",
     lambda s: "temp" in s),

    # Pressure
    ("pressure", "Pressure", "bar", "energetic",
     lambda s: "pressure" in s or "press" in s),

    # Generic distance fallback
    ("distance", "Distance", "nm", "interaction",
     lambda s: "dist" in s),
]


class ObservableResolver:
    """Resolve XVG hints to canonical observable metadata."""

    def resolve(
        self,
        hint:      str,
        xvg_title: str = "",
        xvg_ylabel: str = "",
    ) -> ResolvedObservable:
        """
        Combine filename hint, XVG title, and y-axis label, then return the
        best-matching ResolvedObservable. Falls back to a generic entry if no
        rule matches.
        """
        combined = " ".join([
            _normalize(hint),
            _normalize(xvg_title),
            _normalize(xvg_ylabel),
        ])

        for canonical, display, units, group, predicate in _RULES:
            try:
                if predicate(combined):
                    return ResolvedObservable(canonical=canonical, display=display,
                                              units=units, group=group)
            except Exception:
                continue

        norm = _normalize(hint) or "unknown"
        return ResolvedObservable(canonical=norm, display=hint or "Unknown",
                                  units="", group="other")

    def resolve_from_path(
        self,
        xvg_path:  "str | Path",
        xvg_title: str = "",
        xvg_ylabel: str = "",
    ) -> ResolvedObservable:
        return self.resolve(Path(xvg_path).stem, xvg_title, xvg_ylabel)
