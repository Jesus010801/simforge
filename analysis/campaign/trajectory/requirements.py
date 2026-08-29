"""Helpers around :class:`TrajectoryRequirements`.

The typed model lives in ``analysis.campaign.models``; this module only holds
convenience constructors so observables can express intent declaratively.
"""
from __future__ import annotations

from analysis.campaign.models import TrajectoryRequirements


def raw() -> TrajectoryRequirements:
    return TrajectoryRequirements(rationale="observable operates on raw coordinates")


def whole_only(reason: str) -> TrajectoryRequirements:
    return TrajectoryRequirements(requires_whole_molecules=True, rationale=reason)


def nojump(reason: str) -> TrajectoryRequirements:
    return TrajectoryRequirements(
        requires_whole_molecules=True, requires_nojump=True, rationale=reason,
    )


def fit_to(fit_selection: str, reason: str) -> TrajectoryRequirements:
    """Whole molecules + least-squares rot/trans fit to ``fit_selection``."""
    return TrajectoryRequirements(
        requires_whole_molecules=True,
        fit_selection=fit_selection,
        rationale=reason,
    )


def centered_on(target: str, reason: str) -> TrajectoryRequirements:
    return TrajectoryRequirements(
        requires_whole_molecules=True, centering_target=target, rationale=reason,
    )


def minimum_image(reason: str) -> TrajectoryRequirements:
    """Observable needs whole molecules but does its own minimum-image work;
    no global fitting (would corrupt distance/contact analysis)."""
    return TrajectoryRequirements(
        requires_whole_molecules=True, minimum_image_distances=True, rationale=reason,
    )
