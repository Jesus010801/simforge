"""Tests for core/md_knowledge scientific knowledge base."""
from __future__ import annotations

import math
import pytest

from core.md_knowledge.states import SimulationState, STATE_DESCRIPTIONS, STATE_SEVERITY
from core.md_knowledge.patterns import (
    TemporalPattern, PatternResult, detect_temporal_pattern, MIN_POINTS,
)
from core.md_knowledge.contexts import SystemContext, SYSTEM_CONTEXTS, SoftRange
from core.md_knowledge.heuristics import OBSERVABLE_HEURISTICS, get_heuristic
from core.md_knowledge.evidence import Evidence, EvidenceBundle, accumulate_evidence
from core.md_knowledge.interpreter import ObservableResult, InterpretationResult, interpret_simulation


# ── Fixtures ──────────────────────────────────────────────────────────────────

def _plateau(n=50, mean=0.15, std=0.01):
    """Generate a plateau-like signal with small Gaussian noise."""
    import random
    rng = random.Random(42)
    return [mean + rng.gauss(0, std) for _ in range(n)]


def _drift(n=50, start=0.10, end=0.40):
    """Generate a linearly drifting signal."""
    return [start + (end - start) * i / (n - 1) for i in range(n)]


def _jump(n=50, pre=0.15, post=0.35, jump_at=0.5):
    """Generate a signal that jumps at jump_at fraction of n."""
    import random
    rng = random.Random(7)
    idx = int(n * jump_at)
    return [pre + rng.gauss(0, 0.005) if i < idx else post + rng.gauss(0, 0.005)
            for i in range(n)]


def _times(n=50, total_ns=10.0):
    return [total_ns * i / (n - 1) for i in range(n)]


# ═══════════════════════════════════════════════════════════════════════════════
# SimulationState
# ═══════════════════════════════════════════════════════════════════════════════

class TestSimulationState:
    def test_all_states_have_descriptions(self):
        for s in SimulationState:
            assert s in STATE_DESCRIPTIONS, f"Missing description for {s}"
            assert len(STATE_DESCRIPTIONS[s]) > 20

    def test_all_states_have_severity(self):
        for s in SimulationState:
            assert s in STATE_SEVERITY, f"Missing severity for {s}"
            assert isinstance(STATE_SEVERITY[s], int)

    def test_stable_severity_zero(self):
        assert STATE_SEVERITY[SimulationState.STABLE_EQUILIBRATED] == 0

    def test_nonphysical_highest_severity(self):
        sev = STATE_SEVERITY[SimulationState.NONPHYSICAL_BEHAVIOR]
        assert all(
            STATE_SEVERITY[s] <= sev
            for s in SimulationState
        )

    def test_string_values(self):
        assert SimulationState.STABLE_EQUILIBRATED.value == "stable_equilibrated"
        assert SimulationState.LIGAND_DISSOCIATED.value == "ligand_dissociated"


# ═══════════════════════════════════════════════════════════════════════════════
# TemporalPattern / detect_temporal_pattern
# ═══════════════════════════════════════════════════════════════════════════════

class TestDetectTemporalPattern:
    def test_insufficient_when_too_few_points(self):
        result = detect_temporal_pattern([0.1, 0.2])
        assert result.pattern == TemporalPattern.INSUFFICIENT
        assert result.confidence == 1.0

    def test_plateau_clean_signal(self):
        vals = _plateau(n=100, mean=0.15, std=0.003)
        ts   = _times(100, 20.0)
        result = detect_temporal_pattern(vals, ts, plateau_slope_threshold=0.005)
        assert result.pattern == TemporalPattern.PLATEAU
        assert result.confidence > 0.5

    def test_drift_linear_signal(self):
        vals = _drift(n=100, start=0.10, end=0.40)
        ts   = _times(100, 20.0)
        result = detect_temporal_pattern(vals, ts, drift_slope_threshold=0.002)
        assert result.pattern == TemporalPattern.DRIFT
        assert result.slope > 0

    def test_jump_plateau_detected(self):
        vals = _jump(n=100, pre=0.15, post=0.35)
        ts   = _times(100, 20.0)
        result = detect_temporal_pattern(vals, ts)
        assert result.pattern == TemporalPattern.JUMP_PLATEAU
        assert result.confidence > 0.5

    def test_slope_sign(self):
        # Negative drift
        vals = _drift(n=80, start=0.40, end=0.10)
        ts   = _times(80, 20.0)
        result = detect_temporal_pattern(vals, ts, drift_slope_threshold=0.002)
        assert result.slope < 0

    def test_n_points_recorded(self):
        vals = _plateau(n=60)
        result = detect_temporal_pattern(vals)
        assert result.n_points == 60

    def test_mean_last_approximates_plateau_mean(self):
        vals = _plateau(n=100, mean=0.20, std=0.002)
        result = detect_temporal_pattern(vals)
        assert abs(result.mean_last - 0.20) < 0.02

    def test_plateau_std_is_low_for_clean_signal(self):
        vals = _plateau(n=100, mean=0.15, std=0.003)
        result = detect_temporal_pattern(vals)
        # std of last 20% should be close to 0.003
        assert result.plateau_std < 0.015

    def test_without_times_uses_indices(self):
        vals = _plateau(n=50)
        result = detect_temporal_pattern(vals)   # no times_ns
        assert result.n_points == 50


# ═══════════════════════════════════════════════════════════════════════════════
# SystemContext / SoftRange
# ═══════════════════════════════════════════════════════════════════════════════

class TestSoftRange:
    def _make(self):
        return SoftRange(
            excellent_low=0.0, excellent_high=0.20,
            acceptable_low=0.0, acceptable_high=0.35,
            warning_low=0.0,   warning_high=0.50,
            unit="nm",
        )

    def test_classify_excellent(self):
        sr = self._make()
        assert sr.classify(0.10) == "excellent"

    def test_classify_acceptable(self):
        sr = self._make()
        assert sr.classify(0.28) == "acceptable"

    def test_classify_warning(self):
        sr = self._make()
        assert sr.classify(0.42) == "warning"

    def test_classify_critical(self):
        sr = self._make()
        assert sr.classify(0.80) == "critical"

    def test_score_centre_is_one(self):
        sr = self._make()
        score = sr.score(0.10)  # centre of excellent range
        assert score >= 0.8

    def test_score_outside_is_low(self):
        sr = self._make()
        score = sr.score(1.0)
        assert score < 0.3


class TestSystemContexts:
    def test_all_contexts_defined(self):
        for ctx in SystemContext:
            assert ctx in SYSTEM_CONTEXTS, f"Missing profile for {ctx}"

    def test_globular_rmsd_tight(self):
        profile = SYSTEM_CONTEXTS[SystemContext.GLOBULAR_PROTEIN]
        assert profile.rmsd_nm.excellent_high < 0.25

    def test_idr_rmsd_loose(self):
        profile = SYSTEM_CONTEXTS[SystemContext.IDR]
        assert profile.rmsd_nm.excellent_high > 0.5

    def test_membrane_rmsd_intermediate(self):
        profile = SYSTEM_CONTEXTS[SystemContext.MEMBRANE_PROTEIN]
        glob = SYSTEM_CONTEXTS[SystemContext.GLOBULAR_PROTEIN]
        # membrane protein allows higher RMSD than globular
        assert profile.rmsd_nm.excellent_high > glob.rmsd_nm.excellent_high

    def test_temperature_range_all_contexts(self):
        for ctx, profile in SYSTEM_CONTEXTS.items():
            assert profile.temperature_k.excellent_low >= 280
            assert profile.temperature_k.excellent_high <= 320


# ═══════════════════════════════════════════════════════════════════════════════
# OBSERVABLE_HEURISTICS / get_heuristic
# ═══════════════════════════════════════════════════════════════════════════════

class TestObservableHeuristics:
    EXPECTED_OBSERVABLES = {
        "rmsd", "rg", "temperature", "pressure", "potential_energy",
        "rmsf", "sasa", "hbonds", "ligand_pocket_distance", "secondary_structure",
    }

    def test_all_observables_registered(self):
        assert self.EXPECTED_OBSERVABLES <= set(OBSERVABLE_HEURISTICS.keys())

    def test_each_observable_has_all_contexts(self):
        for obs in self.EXPECTED_OBSERVABLES:
            for ctx in SystemContext:
                h = get_heuristic(obs, ctx)
                assert h is not None, f"Missing heuristic for {obs} + {ctx}"

    def test_heuristic_has_rules(self):
        h = get_heuristic("rmsd", SystemContext.GLOBULAR_PROTEIN)
        assert h is not None
        assert len(h.rules) >= 3

    def test_heuristic_has_citations_for_rmsd(self):
        h = get_heuristic("rmsd", SystemContext.GLOBULAR_PROTEIN)
        assert len(h.citations) >= 1

    def test_fallback_to_unknown(self):
        h = get_heuristic("rmsd", SystemContext.UNKNOWN)
        assert h is not None

    def test_ligand_heuristic_has_dissociation_rule(self):
        h = get_heuristic("ligand_pocket_distance", SystemContext.PROTEIN_LIGAND_COMPLEX)
        msgs = [r.interpretation for r in h.rules]
        assert any("dissociat" in m.lower() for m in msgs)


# ═══════════════════════════════════════════════════════════════════════════════
# EvidenceBundle / accumulate_evidence
# ═══════════════════════════════════════════════════════════════════════════════

class TestEvidenceBundle:
    def _plateau_pair(self, obs_name="rmsd", mean=0.15):
        vals = _plateau(n=80, mean=mean, std=0.005)
        pr = detect_temporal_pattern(vals, _times(80, 20.0), plateau_slope_threshold=0.005)
        return obs_name, (pr, pr.mean_last)

    def test_bundle_adds_evidence(self):
        obs_name, pair = self._plateau_pair()
        bundle = accumulate_evidence({obs_name: pair})
        assert len(bundle.items) == 1

    def test_plateau_votes_stable(self):
        obs_name, pair = self._plateau_pair()
        bundle = accumulate_evidence({obs_name: pair})
        ev = bundle.items[0]
        assert ev.state_vote == SimulationState.STABLE_EQUILIBRATED

    def test_drift_votes_drifting(self):
        vals = _drift(n=80, start=0.10, end=0.40)
        pr = detect_temporal_pattern(vals, _times(80, 20.0), drift_slope_threshold=0.002)
        bundle = accumulate_evidence({"rmsd": (pr, pr.mean_last)})
        ev = bundle.items[0]
        assert ev.state_vote == SimulationState.DRIFTING

    def test_ligand_dissociated_when_far(self):
        vals = _plateau(n=80, mean=1.0, std=0.01)
        pr = detect_temporal_pattern(vals, _times(80, 20.0))
        bundle = accumulate_evidence({"ligand_pocket_distance": (pr, 1.0)})
        ev = bundle.items[0]
        assert ev.state_vote == SimulationState.LIGAND_DISSOCIATED
        assert ev.severity == "error"

    def test_high_temperature_flagged(self):
        vals = _plateau(n=80, mean=500.0, std=1.0)
        pr = detect_temporal_pattern(vals, _times(80, 20.0))
        bundle = accumulate_evidence({"temperature": (pr, 500.0)})
        ev = bundle.items[0]
        assert ev.state_vote == SimulationState.NONPHYSICAL_BEHAVIOR
        assert ev.severity == "error"

    def test_warnings_property(self):
        obs_name, pair = self._plateau_pair()
        bundle = accumulate_evidence({obs_name: pair})
        # A plateau signal has no warnings
        assert isinstance(bundle.warnings, list)

    def test_dominant_state_most_votes(self):
        _, p1 = self._plateau_pair("rmsd")
        _, p2 = self._plateau_pair("rg")
        bundle = accumulate_evidence({"rmsd": p1, "rg": p2})
        assert bundle.dominant_state() == SimulationState.STABLE_EQUILIBRATED

    def test_compute_confidence_positive(self):
        _, p1 = self._plateau_pair()
        bundle = accumulate_evidence({"rmsd": p1})
        assert 0.0 < bundle.compute_confidence() <= 1.0

    def test_error_severity_lowers_confidence(self):
        # Ligand dissociated (error) vs plateau rmsd (info)
        _, rmsd_pair = self._plateau_pair("rmsd")
        vals_lig = _plateau(n=80, mean=1.2, std=0.01)
        pr_lig = detect_temporal_pattern(vals_lig, _times(80, 20.0))
        bundle_no_error   = accumulate_evidence({"rmsd": rmsd_pair})
        bundle_with_error = accumulate_evidence({"rmsd": rmsd_pair, "ligand_pocket_distance": (pr_lig, 1.2)})
        assert bundle_with_error.overall_confidence <= bundle_no_error.overall_confidence


# ═══════════════════════════════════════════════════════════════════════════════
# interpret_simulation
# ═══════════════════════════════════════════════════════════════════════════════

class TestInterpretSimulation:
    def _make_obs(self, name, vals, total_ns=20.0):
        return ObservableResult(
            name=name,
            values=vals,
            times_ns=_times(len(vals), total_ns),
        )

    def test_empty_observables_returns_insufficient(self):
        result = interpret_simulation([], SystemContext.GLOBULAR_PROTEIN)
        assert result.state == SimulationState.INSUFFICIENT_SAMPLING
        assert result.confidence == 0.0

    def test_stable_rmsd_gives_stable_state(self):
        vals = _plateau(n=100, mean=0.15, std=0.004)
        obs  = self._make_obs("rmsd", vals)
        result = interpret_simulation([obs], SystemContext.GLOBULAR_PROTEIN)
        assert result.state == SimulationState.STABLE_EQUILIBRATED
        assert result.confidence > 0.4

    def test_drifting_rmsd_gives_drifting_state(self):
        vals = _drift(n=100, start=0.10, end=0.50)
        obs  = self._make_obs("rmsd", vals)
        result = interpret_simulation([obs], SystemContext.GLOBULAR_PROTEIN)
        assert result.state == SimulationState.DRIFTING

    def test_jump_rmsd_gives_conformational_transition(self):
        vals = _jump(n=100, pre=0.15, post=0.40)
        obs  = self._make_obs("rmsd", vals)
        result = interpret_simulation([obs], SystemContext.GLOBULAR_PROTEIN)
        assert result.state == SimulationState.CONFORMATIONAL_TRANSITION

    def test_ligand_dissociation_overrides(self):
        rmsd_obs = self._make_obs("rmsd", _plateau(n=80, mean=0.15, std=0.004))
        lig_obs  = self._make_obs("ligand_pocket_distance", _plateau(n=80, mean=1.2, std=0.02))
        result = interpret_simulation([rmsd_obs, lig_obs], SystemContext.PROTEIN_LIGAND_COMPLEX)
        assert result.state == SimulationState.LIGAND_DISSOCIATED

    def test_high_temperature_gives_nonphysical(self):
        temp_obs = self._make_obs("temperature", _plateau(n=80, mean=500.0, std=1.0))
        result = interpret_simulation([temp_obs], SystemContext.GLOBULAR_PROTEIN)
        assert result.state == SimulationState.NONPHYSICAL_BEHAVIOR

    def test_idr_accepts_high_rmsd(self):
        # RMSD=1.5 nm for an IDR — should not be catastrophic
        vals = _plateau(n=100, mean=1.5, std=0.05)
        obs  = self._make_obs("rmsd", vals)
        result = interpret_simulation([obs], SystemContext.IDR)
        # State may be stable or partially converged but NOT problematic
        assert result.state not in (
            SimulationState.NONPHYSICAL_BEHAVIOR,
            SimulationState.UNSTABLE,
        )

    def test_globular_high_rmsd_is_problematic(self):
        # RMSD=1.5 nm for a globular protein — should get a lower confidence
        vals = _plateau(n=100, mean=1.5, std=0.05)
        obs  = self._make_obs("rmsd", vals)
        result_idr  = interpret_simulation([obs], SystemContext.IDR)
        result_glob = interpret_simulation([obs], SystemContext.GLOBULAR_PROTEIN)
        # Globular confidence should be lower for the same high RMSD
        assert result_glob.confidence <= result_idr.confidence

    def test_summary_non_empty(self):
        vals = _plateau(n=80, mean=0.15, std=0.004)
        obs  = self._make_obs("rmsd", vals)
        result = interpret_simulation([obs], SystemContext.GLOBULAR_PROTEIN)
        assert len(result.summary) > 20

    def test_quality_tier_excellent_for_well_converged(self):
        vals = _plateau(n=200, mean=0.12, std=0.003)
        ts   = _times(200, 50.0)
        obs  = ObservableResult("rmsd", vals, ts)
        result = interpret_simulation([obs], SystemContext.GLOBULAR_PROTEIN)
        assert result.quality_tier in ("EXCELLENT", "GOOD", "ACCEPTABLE")

    def test_quality_tier_failed_for_nonphysical(self):
        obs = self._make_obs("temperature", _plateau(n=80, mean=600.0, std=2.0))
        result = interpret_simulation([obs], SystemContext.GLOBULAR_PROTEIN)
        assert result.quality_tier == "FAILED"

    def test_state_description_available(self):
        obs = self._make_obs("rmsd", _plateau(n=80, mean=0.15, std=0.004))
        result = interpret_simulation([obs], SystemContext.GLOBULAR_PROTEIN)
        assert len(result.state_description) > 10

    def test_context_stored(self):
        obs = self._make_obs("rmsd", _plateau(n=50, mean=0.15, std=0.005))
        result = interpret_simulation([obs], SystemContext.MEMBRANE_PROTEIN)
        assert result.context == SystemContext.MEMBRANE_PROTEIN

    def test_pattern_results_populated(self):
        obs = self._make_obs("rmsd", _plateau(n=80, mean=0.15, std=0.004))
        result = interpret_simulation([obs], SystemContext.GLOBULAR_PROTEIN)
        assert "rmsd" in result.pattern_results


# ═══════════════════════════════════════════════════════════════════════════════
# quality_classifier integration with context
# ═══════════════════════════════════════════════════════════════════════════════

class TestQualityClassifierWithContext:
    """Verify classify_run respects the context= argument."""

    def _xvg(self, values, times_ps):
        from pathlib import Path
        from runtime.xvg_parser import XVGData, XVGSeries
        series = XVGSeries(name="test", values=values)
        return XVGData(title="test", xlabel="t", ylabel="v",
                       time_ps=times_ps, series=[series], source=Path("/tmp/test.xvg"))

    def _plateau_xvg(self, mean=0.15, total_ps=20000, n=200):
        import random
        rng = random.Random(99)
        times = [total_ps * i / (n - 1) for i in range(n)]
        vals  = [mean + rng.gauss(0, 0.003) for _ in range(n)]
        return self._xvg(vals, times)

    def test_no_context_uses_legacy_classifier(self):
        from runtime.quality_classifier import classify_run
        rmsd = self._plateau_xvg(mean=0.15, total_ps=10000)
        report = classify_run(rmsd_data=rmsd)
        # Should not crash; quality is one of the known enums
        from runtime.quality_classifier import RunQuality
        assert report.quality in RunQuality

    def test_with_context_returns_report(self):
        from runtime.quality_classifier import classify_run, RunQuality
        rmsd = self._plateau_xvg(mean=0.15, total_ps=20000)
        report = classify_run(rmsd_data=rmsd, context="globular_protein")
        assert report.quality in RunQuality
        assert 0.0 <= report.confidence <= 1.0

    def test_context_in_metrics(self):
        from runtime.quality_classifier import classify_run
        rmsd = self._plateau_xvg(mean=0.15, total_ps=20000)
        report = classify_run(rmsd_data=rmsd, context="globular_protein")
        assert report.metrics.get("context") == "globular_protein"

    def test_unknown_context_falls_back(self):
        from runtime.quality_classifier import classify_run
        rmsd = self._plateau_xvg(mean=0.15, total_ps=20000)
        # Should not raise
        report = classify_run(rmsd_data=rmsd, context="nonexistent_context_xyz")
        assert report is not None
        assert report.metrics.get("context") == "unknown"

    def test_high_rmsd_bad_for_globular(self):
        from runtime.quality_classifier import classify_run, RunQuality
        rmsd = self._plateau_xvg(mean=1.5, total_ps=20000)
        report = classify_run(rmsd_data=rmsd, context="globular_protein")
        # High RMSD should result in lower confidence than low RMSD
        report_good = classify_run(rmsd_data=self._plateau_xvg(mean=0.12, total_ps=20000),
                                   context="globular_protein")
        assert report.confidence <= report_good.confidence
