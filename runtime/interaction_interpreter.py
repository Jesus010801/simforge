"""Multi-observable interaction state interpreter.

Evaluates explicit weighted rules against normalized observable signals to
classify molecular interaction states (stable_binding, weak_binding,
ligand_destabilization, etc.) with confidence scores and evidence.

Design principles:
- No black-box scoring — every conclusion traces back to named observables.
- Sigmoid normalization keeps middle values at 0.5, avoids extreme artifacts.
- Rules require a minimum fraction of their conditions to be present in data.
- If no rule fires above threshold → uncertain_behavior.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from runtime.study_models import Study

from runtime.synthesis_models import SignalProfile, RuleMatch, SystemSynthesis


# ─── Observable polarity ─────────────────────────────────────────────────────
# Lower is favorable (stability/proximity): direction_score = 1 - normalized
_LOWER_IS_BETTER = frozenset({
    "protein_rmsd", "ligand_rmsd", "rmsf", "mindist",
    "catalytic_distance", "distance", "radius_of_gyration",
})
# Higher is favorable (interaction strength): direction_score = normalized
_HIGHER_IS_BETTER = frozenset({"contacts", "hydrogen_bonds"})


def _direction_score(obs: str, normalized_value: float) -> float:
    if obs in _LOWER_IS_BETTER:
        return 1.0 - normalized_value
    if obs in _HIGHER_IS_BETTER:
        return normalized_value
    return 0.5   # neutral


# ─── Sigmoid normalization ────────────────────────────────────────────────────

def _sigmoid(x: float) -> float:
    return 1.0 / (1.0 + math.exp(-x))


def _normalize_sigmoid(value: float, study_mean: float, study_std: float) -> float:
    """Sigmoid normalization centered on study mean.

    Returns ≈0.5 at study_mean, approaches 1.0 two std above, 0.0 below.
    """
    if study_std < 1e-10:
        return 0.5
    return _sigmoid((value - study_mean) / study_std)


# ─── Rule definitions ─────────────────────────────────────────────────────────

@dataclass
class RuleCondition:
    observable: str
    expected:   str    # "good" | "bad"  (good = favorable value expected for this state)
    weight:     float


@dataclass
class InteractionRule:
    state:            str
    description:      str
    conditions:       list[RuleCondition]
    min_data_fraction: float = 0.50   # min fraction of total weight that must be available


_RULES: list[InteractionRule] = [
    # ── Binding states ────────────────────────────────────────────────────────
    InteractionRule(
        state="stable_binding",
        description="Strong, persistent molecular interaction",
        conditions=[
            RuleCondition("contacts",           "good", 0.40),
            RuleCondition("ligand_rmsd",        "good", 0.25),
            RuleCondition("protein_rmsd",       "good", 0.20),
            RuleCondition("catalytic_distance", "good", 0.15),
        ],
        min_data_fraction=0.45,
    ),
    InteractionRule(
        state="interaction_persistent",
        description="Consistent interaction profile across observables",
        conditions=[
            RuleCondition("contacts",           "good", 0.45),
            RuleCondition("catalytic_distance", "good", 0.30),
            RuleCondition("mindist",            "good", 0.25),
        ],
        min_data_fraction=0.45,
    ),
    InteractionRule(
        state="weak_binding",
        description="Interaction maintained but with elevated ligand mobility",
        conditions=[
            RuleCondition("ligand_rmsd",  "bad",  0.35),
            RuleCondition("contacts",     "bad",  0.30),
            RuleCondition("protein_rmsd", "good", 0.20),   # protein still stable
            RuleCondition("mindist",      "bad",  0.15),
        ],
        min_data_fraction=0.55,  # requires ligand + contact evidence
    ),
    InteractionRule(
        state="transient_binding",
        description="Interaction present but not consistently maintained",
        conditions=[
            RuleCondition("contacts",           "bad", 0.45),
            RuleCondition("mindist",            "bad", 0.30),
            RuleCondition("catalytic_distance", "bad", 0.25),
        ],
        min_data_fraction=0.50,
    ),
    InteractionRule(
        state="ligand_destabilization",
        description="Ligand losing interaction stability",
        conditions=[
            RuleCondition("ligand_rmsd",        "bad", 0.35),
            RuleCondition("contacts",           "bad", 0.35),
            RuleCondition("mindist",            "bad", 0.15),
            RuleCondition("catalytic_distance", "bad", 0.15),
        ],
        min_data_fraction=0.50,
    ),
    InteractionRule(
        state="possible_dissociation",
        description="Evidence consistent with ligand separation from binding site",
        conditions=[
            RuleCondition("ligand_rmsd", "bad", 0.40),
            RuleCondition("mindist",     "bad", 0.35),
            RuleCondition("contacts",    "bad", 0.25),
        ],
        min_data_fraction=0.60,
    ),
    # ── Structural states ─────────────────────────────────────────────────────
    InteractionRule(
        state="structurally_stable",
        description="Well-equilibrated, converged protein conformation",
        conditions=[
            RuleCondition("protein_rmsd",       "good", 0.50),
            RuleCondition("rmsf",               "good", 0.30),
            RuleCondition("radius_of_gyration", "good", 0.20),
        ],
        min_data_fraction=0.45,
    ),
    InteractionRule(
        state="flexible_but_stable",
        description="Elevated flexibility without indication of unfolding",
        conditions=[
            RuleCondition("rmsf",               "bad",  0.45),
            RuleCondition("protein_rmsd",       "good", 0.35),
            RuleCondition("radius_of_gyration", "good", 0.20),
        ],
        min_data_fraction=0.55,
    ),
    InteractionRule(
        state="conformational_rearrangement",
        description="System undergoing significant structural reorganization",
        conditions=[
            RuleCondition("protein_rmsd",       "bad", 0.50),
            RuleCondition("rmsf",               "bad", 0.30),
            RuleCondition("radius_of_gyration", "bad", 0.20),
        ],
        min_data_fraction=0.60,  # needs multiple structural indicators
    ),
]

# Composite score weights
_BINDING_WEIGHTS:   dict[str, float] = {
    "contacts": 0.40, "ligand_rmsd": 0.30, "mindist": 0.20, "catalytic_distance": 0.10,
}
_STABILITY_WEIGHTS: dict[str, float] = {
    "protein_rmsd": 0.55, "rmsf": 0.30, "radius_of_gyration": 0.15,
}

# Minimum score for a rule to be considered "active"
_SCORE_THRESHOLD = 0.55


# ─── Signal normalization ─────────────────────────────────────────────────────

def build_signal_profiles(study: "Study") -> dict[str, SignalProfile]:
    """Compute per-system normalized signals and direction scores."""
    from runtime.study_models import Study  # avoid circular import at module level

    profiles: dict[str, SignalProfile] = {}

    # Study-wide stats per observable
    for obs in study.observables_detected:
        sys_means = [
            sg.aggregate[obs].mean
            for sg in study.systems.values()
            if obs in sg.aggregate
        ]
        if not sys_means:
            continue

        st_mean = sum(sys_means) / len(sys_means)
        st_std  = (
            (sum((v - st_mean) ** 2 for v in sys_means) / max(len(sys_means) - 1, 1)) ** 0.5
        )

        # Study-wide max inter-replica std (for consistency normalization)
        all_stds = [
            sg.aggregate[obs].std
            for sg in study.systems.values()
            if obs in sg.aggregate
        ]
        max_std = max(all_stds) if all_stds else 1.0

        for sys_name, sg in study.systems.items():
            if obs not in sg.aggregate:
                continue

            agg = sg.aggregate[obs]
            if sys_name not in profiles:
                profiles[sys_name] = SignalProfile(system=sys_name)

            p = profiles[sys_name]
            norm = _normalize_sigmoid(agg.mean, st_mean, st_std)
            p.normalized[obs]     = norm
            p.direction_score[obs] = _direction_score(obs, norm)
            # Consistency: 0 = most variable, 1 = most consistent
            p.consistency[obs] = 1.0 - (agg.std / (max_std + 1e-10))

    return profiles


# ─── Rule evaluation ─────────────────────────────────────────────────────────

def _condition_score(cond: RuleCondition, profile: SignalProfile) -> Optional[float]:
    """Return [0,1] score for one condition, or None if observable not available."""
    if cond.observable not in profile.direction_score:
        return None
    ds = profile.direction_score[cond.observable]
    return ds if cond.expected == "good" else (1.0 - ds)


def evaluate_rules(
    profile:           SignalProfile,
    study:             "Study",
    consensus_map:     dict[str, float] | None = None,
) -> list[RuleMatch]:
    """Evaluate all interaction rules for one system, return sorted list."""
    matches: list[RuleMatch] = []

    for rule in _RULES:
        total_weight = sum(c.weight for c in rule.conditions)
        if total_weight < 1e-10:
            continue

        avail_weight = 0.0
        raw_score    = 0.0
        supporting:  list[str] = []
        opposing:    list[str] = []
        n_available  = 0

        for cond in rule.conditions:
            cs = _condition_score(cond, profile)
            if cs is None:
                continue
            n_available += 1
            avail_weight += cond.weight
            raw_score    += cs * cond.weight

            obs_label = cond.observable.replace("_", " ")
            dir_label = (
                "favorable" if cond.expected == "good" else "unfavorable"
            )
            if cs >= 0.60:
                supporting.append(f"{dir_label} {obs_label}")
            elif cs <= 0.40:
                opposing.append(f"{obs_label} not {dir_label}")

        if avail_weight < rule.min_data_fraction * total_weight:
            continue

        score = raw_score / avail_weight if avail_weight > 0 else 0.0

        # Apply inter-replica consistency boost/penalty
        obs_in_rule = [c.observable for c in rule.conditions]
        consistencies = [
            profile.consistency[o]
            for o in obs_in_rule
            if o in profile.consistency
        ]
        consistency_factor = (
            sum(consistencies) / len(consistencies) if consistencies else 0.8
        )
        # Blend: low consistency slightly reduces confidence
        confidence_mult = 0.7 + 0.3 * consistency_factor

        # Optionally use cross-observable consensus
        if consensus_map:
            key = profile.system
            cval = consensus_map.get(key, 1.0)
            confidence_mult *= cval

        matches.append(RuleMatch(
            state                  = rule.state,
            description            = rule.description,
            score                  = score,
            confidence             = min(1.0, score * confidence_mult),
            n_conditions_met       = n_available,
            n_conditions_available = n_available,
            supporting             = supporting,
            opposing               = opposing,
        ))

    matches.sort(key=lambda m: -m.score)
    return matches


# ─── Composite scores ─────────────────────────────────────────────────────────

def _weighted_score(
    direction_scores: dict[str, float],
    weights:          dict[str, float],
) -> Optional[float]:
    """Return weighted mean of available direction scores, or None."""
    total_w = avail_w = 0.0
    score = 0.0
    for obs, w in weights.items():
        total_w += w
        if obs in direction_scores:
            score  += direction_scores[obs] * w
            avail_w += w
    if avail_w < 0.3 * total_w:  # need at least 30% weight coverage
        return None
    return score / avail_w


# ─── Main interpretation ──────────────────────────────────────────────────────

def interpret_system(
    profile:        SignalProfile,
    study:          "Study",
    consensus_map:  dict[str, float] | None = None,
) -> SystemSynthesis:
    """Build a SystemSynthesis from a SignalProfile."""
    matches = evaluate_rules(profile, study, consensus_map)

    # Composite scores
    binding_score  = _weighted_score(profile.direction_score, _BINDING_WEIGHTS)
    stability_score = _weighted_score(profile.direction_score, _STABILITY_WEIGHTS)

    # Mean convergence from aggregate
    sg = study.systems.get(profile.system)
    conv_scores = (
        [agg.mean_convergence for agg in sg.aggregate.values()]
        if sg else []
    )
    mean_conv = sum(conv_scores) / len(conv_scores) if conv_scores else 0.5

    # Composite: binding × 0.5 + stability × 0.3 + convergence × 0.2
    # Redistribute weights if components missing
    components: list[tuple[float, float]] = []
    if binding_score  is not None: components.append((binding_score,  0.50))
    if stability_score is not None: components.append((stability_score, 0.30))
    components.append((mean_conv, 0.20))

    total_w = sum(w for _, w in components)
    composite = sum(v * w for v, w in components) / total_w if total_w else mean_conv

    # Primary state
    active = [m for m in matches if m.score >= _SCORE_THRESHOLD]
    if active:
        top = active[0]
        # Check if top two are too close (ambiguous)
        if len(active) >= 2 and (active[0].score - active[1].score) < 0.08:
            # Tiebreaker: among tied rules, prefer the one with the most
            # conditions satisfied (more evidence = more specific classification).
            # This avoids floating-point ordering artifacts.
            best = max(active, key=lambda m: m.n_conditions_available)
            second_best_n = sorted(
                active, key=lambda m: -m.n_conditions_available
            )[1].n_conditions_available
            if best.n_conditions_available > second_best_n:
                primary_state = best.state
                primary_conf  = best.confidence * 0.85  # slight penalty for near-tie
            else:
                primary_state = "uncertain_behavior"
                primary_conf  = max(active[0].confidence, 0.3)
        else:
            primary_state = top.state
            primary_conf  = top.confidence
    else:
        primary_state = "uncertain_behavior"
        primary_conf  = 0.30

    # Build evidence list
    evidence: list[str] = []
    for obs, ds in sorted(profile.direction_score.items(), key=lambda x: -abs(x[1] - 0.5)):
        agg = sg.aggregate.get(obs) if sg else None
        if agg is None:
            continue
        units = study.observable_units.get(obs, "")
        disp  = study.observable_display.get(obs, obs)
        val   = f"{agg.mean:.3g}{' ' + units if units else ''}"
        if ds >= 0.65:
            evidence.append(f"favorable {disp} ({val})")
        elif ds <= 0.35:
            evidence.append(f"unfavorable {disp} ({val})")

    # One-line explanation
    state_desc = next(
        (r.description for r in _RULES if r.state == primary_state),
        primary_state.replace("_", " "),
    )
    ev_summary = "; ".join(evidence[:3]) if evidence else "limited observable data"
    explanation = (
        f"{profile.system} classified as {primary_state.replace('_', ' ')} "
        f"(confidence {primary_conf:.2f}): {state_desc.lower()}. "
        f"Evidence: {ev_summary}."
    )

    return SystemSynthesis(
        system            = profile.system,
        primary_state     = primary_state,
        primary_confidence = primary_conf,
        binding_score     = binding_score  if binding_score  is not None else 0.5,
        stability_score   = stability_score if stability_score is not None else 0.5,
        composite_score   = composite,
        rule_matches      = matches,
        evidence          = evidence,
        explanation       = explanation,
    )


def interpret_all(
    study:         "Study",
    consensus_map: dict[str, float] | None = None,
) -> dict[str, SystemSynthesis]:
    """Interpret all systems in a study."""
    profiles = build_signal_profiles(study)
    return {
        sys_name: interpret_system(profiles[sys_name], study, consensus_map)
        for sys_name in study.systems
        if sys_name in profiles
    }
