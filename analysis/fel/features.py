"""Build a validated two-feature time-aligned table from two XVGSeries."""
from __future__ import annotations

import math
from typing import Optional

from analysis.fel.models import FeatureTable, FELWarning, XVGSeries

_NS_FACTORS: dict[str, float] = {"ps": 1e-3, "ns": 1.0, "us": 1e3}


def _to_ns(raw: list[float], unit: Optional[str]) -> list[float]:
    factor = _NS_FACTORS.get(unit or "ps", 1e-3)
    return [t * factor for t in raw]


def build_feature_table(
    sx: XVGSeries,
    sy: XVGSeries,
    *,
    override_unit_x: Optional[str] = None,
    override_unit_y: Optional[str] = None,
    x_label: Optional[str] = None,
    y_label: Optional[str] = None,
    time_tolerance: float = 1e-6,
) -> FeatureTable:
    """Build and validate a FeatureTable from two aligned XVGSeries.

    Both series must have the same number of frames. Their time axes (after
    conversion to nanoseconds) must agree within `time_tolerance` ns at every
    frame. If they differ only by a ps↔ns unit mismatch the conversion is
    applied automatically and recorded in the provenance note.

    Raises
    ------
    ValueError
        For length mismatch or irreconcilable time-axis misalignment.
    """
    ux = override_unit_x or sx.x_unit
    uy = override_unit_y or sy.x_unit

    tx_ns = _to_ns(sx.time_raw, ux)
    ty_ns = _to_ns(sy.time_raw, uy)

    warnings: list[FELWarning] = []
    conversion_note = _build_conversion_note(ux, uy)

    # ── Length check ──────────────────────────────────────────────────────────
    if len(tx_ns) != len(ty_ns):
        raise ValueError(
            f"Feature series have different lengths: "
            f"X has {len(tx_ns)} frames, Y has {len(ty_ns)} frames. "
            "Both XVG files must cover the same trajectory."
        )

    # ── Time alignment check ──────────────────────────────────────────────────
    max_diff = max(abs(a - b) for a, b in zip(tx_ns, ty_ns))
    if max_diff > time_tolerance:
        # Try detecting a factor-of-1000 ps/ns mismatch
        corrected = _try_unit_correction(sx, sy, ux, uy, time_tolerance)
        if corrected is not None:
            tx_ns, ty_ns, ux, uy, conversion_note = corrected
            warnings.append(FELWarning(
                "time_unit_corrected",
                f"Time axes differed by ~1000×; applied ps↔ns conversion. "
                f"X unit: {ux}, Y unit: {uy}.",
                severity="info",
            ))
        else:
            raise ValueError(
                f"Time axes of X and Y series do not align within tolerance "
                f"{time_tolerance} ns (max diff = {max_diff:.4g} ns). "
                "Use --time-unit-x / --time-unit-y to specify units explicitly, "
                "or check that both XVG files cover the same trajectory."
            )

    # ── Warn about skipped lines ──────────────────────────────────────────────
    for label, s in (("X", sx), ("Y", sy)):
        if s.n_skipped > 0:
            pct = 100.0 * s.n_skipped / max(1, s.n_frames() + s.n_skipped)
            warnings.append(FELWarning(
                f"skipped_lines_{label.lower()}",
                f"{s.n_skipped} malformed line(s) skipped in {s.path.name} ({pct:.1f}%)",
                severity="warn" if pct > 5 else "info",
            ))

    conv_x = _NS_FACTORS.get(ux or "ps", 1e-3)
    conv_y = _NS_FACTORS.get(uy or "ps", 1e-3)

    return FeatureTable(
        time_ns              = tx_ns,
        x                    = list(sx.data_values),
        y                    = list(sy.data_values),
        x_label              = x_label or sx.column_name,
        y_label              = y_label or sy.column_name,
        x_unit_original      = sx.x_unit,
        y_unit_original      = sy.x_unit,
        time_conversion_note = conversion_note,
        n_frames             = len(tx_ns),
        warnings             = warnings,
        x_unit               = sx.data_unit,
        y_unit               = sy.data_unit,
        time_unit_normalized = "ns",
        time_conversion_x_to_ns = conv_x,
        time_conversion_y_to_ns = conv_y,
    )


def _build_conversion_note(ux: Optional[str], uy: Optional[str]) -> str:
    notes: list[str] = []
    for label, u in (("X", ux), ("Y", uy)):
        detected = u or "ps (assumed — unit not declared in header)"
        notes.append(f"{label}: {detected}")
    return "Time units — " + ", ".join(notes) + " → all converted to ns."


def _try_unit_correction(
    sx: XVGSeries,
    sy: XVGSeries,
    ux: Optional[str],
    uy: Optional[str],
    tol: float,
) -> Optional[tuple]:
    """Try swapping ps↔ns on one series to reconcile time axes.

    Returns (tx_ns, ty_ns, new_ux, new_uy, note) or None.
    """
    _SWAPS = [("ps", "ns"), ("ns", "ps")]
    for axis, (from_u, to_u) in ((0, ("ps", "ns")), (0, ("ns", "ps")),
                                  (1, ("ps", "ns")), (1, ("ns", "ps"))):
        if axis == 0 and ux == from_u:
            tx_ns2 = _to_ns(sx.time_raw, to_u)
            ty_ns2 = _to_ns(sy.time_raw, uy)
        elif axis == 1 and uy == from_u:
            tx_ns2 = _to_ns(sx.time_raw, ux)
            ty_ns2 = _to_ns(sy.time_raw, to_u)
        else:
            continue
        if len(tx_ns2) == len(ty_ns2):
            diff = max(abs(a - b) for a, b in zip(tx_ns2, ty_ns2))
            if diff <= tol:
                new_ux = to_u if axis == 0 else ux
                new_uy = to_u if axis == 1 else uy
                note = (
                    f"Time unit auto-corrected: axis {'X' if axis == 0 else 'Y'} "
                    f"reinterpreted as {to_u} instead of {from_u}."
                )
                return tx_ns2, ty_ns2, new_ux, new_uy, note
    return None
