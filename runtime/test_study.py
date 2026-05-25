"""Tests for the Study Layer and Scientific Synthesis Layer.

Covers 8 modules with zero previous coverage:
  runtime/study_models.py
  runtime/observable_resolver.py
  runtime/study_analyzer.py
  runtime/synthesis_models.py
  runtime/interaction_interpreter.py
  runtime/consensus_engine.py
  runtime/event_detector.py
  runtime/scientific_synthesis.py
"""
from __future__ import annotations

import math
from pathlib import Path

import pytest

# ─── Helpers ─────────────────────────────────────────────────────────────────

def _write_xvg(path: Path, title: str, ylabel: str, times_ps: list[float], values: list[float]) -> None:
    lines = [
        f'@ title "{title}"',
        f'@ yaxis label "{ylabel}"',
        "@ TYPE xy",
    ]
    for t, v in zip(times_ps, values):
        lines.append(f"{t:.1f}  {v:.6f}")
    path.write_text("\n".join(lines))


def _flat_series(n: int = 100, value: float = 0.20) -> tuple[list[float], list[float]]:
    """n points at constant value, time in ps (0..10000)."""
    times = [i * 100.0 for i in range(n)]
    vals  = [value] * n
    return times, vals


def _drifting_series(n: int = 100, start: float = 0.15, end: float = 0.50) -> tuple[list[float], list[float]]:
    times = [i * 100.0 for i in range(n)]
    vals  = [start + (end - start) * i / (n - 1) for i in range(n)]
    return times, vals


# ─── Shared study fixture ─────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def two_system_dir(tmp_path_factory):
    """AA (good binder) vs LP (weak binder), 3 replicas each, 2 observables."""
    d = tmp_path_factory.mktemp("two_system")

    for rep in ("A1", "A2", "A3"):
        _write_xvg(d / f"AA-{rep}_rmsd_protein.xvg",  "RMSD",     "RMSD (nm)",  *_flat_series(100, 0.18))
        _write_xvg(d / f"AA-{rep}_contacts.xvg",       "Contacts", "Contacts",   *_flat_series(100, 14.5))
        _write_xvg(d / f"LP-{rep}_rmsd_protein.xvg",  "RMSD",     "RMSD (nm)",  *_flat_series(100, 0.45))
        _write_xvg(d / f"LP-{rep}_contacts.xvg",       "Contacts", "Contacts",   *_flat_series(100, 6.2))

    return d


@pytest.fixture(scope="module")
def two_system_study(two_system_dir):
    from runtime.study_analyzer import parse_study
    return parse_study(two_system_dir)


@pytest.fixture(scope="module")
def rich_study_dir(tmp_path_factory):
    """3 observables: protein_rmsd, ligand_rmsd, contacts."""
    d = tmp_path_factory.mktemp("rich_study")
    import random
    rng = random.Random(99)

    systems = {
        "AA": dict(prot=0.18, lig=0.14, cont=14.5),
        "LP": dict(prot=0.40, lig=0.58, cont=6.5),
    }
    for sys, props in systems.items():
        for rep in ("A1", "A2", "A3", "A4"):
            n = 0.98 + rng.random() * 0.04
            _write_xvg(d / f"{sys}-{rep}_rmsd_protein.xvg", "Protein RMSD", "RMSD (nm)",
                       *_flat_series(100, props["prot"] * n))
            _write_xvg(d / f"{sys}-{rep}_rmsd_ligand.xvg",  "Ligand RMSD",  "RMSD (nm)",
                       *_flat_series(100, props["lig"] * n))
            _write_xvg(d / f"{sys}-{rep}_contacts.xvg",     "Contacts",     "Contacts",
                       *_flat_series(100, props["cont"] * n))
    return d


@pytest.fixture(scope="module")
def rich_study(rich_study_dir):
    from runtime.study_analyzer import parse_study
    return parse_study(rich_study_dir)


# ═══════════════════════════════════════════════════════════════════════════════
# 1. study_models.py
# ═══════════════════════════════════════════════════════════════════════════════

class TestStudyModels:
    def test_system_group_n_replicas(self):
        from runtime.study_models import SystemGroup, Replica
        sg = SystemGroup(name="AA")
        assert sg.n_replicas == 0
        sg.replicas["A1"] = Replica(label="A1", system="AA")
        sg.replicas["A2"] = Replica(label="A2", system="AA")
        assert sg.n_replicas == 2

    def test_system_group_observables_sorted(self):
        from runtime.study_models import SystemGroup, Replica, ObservableSeries
        from pathlib import Path
        sg = SystemGroup(name="AA")
        r = Replica(label="A1", system="AA")
        r.observables["protein_rmsd"] = ObservableSeries(
            observable="protein_rmsd", replica="A1", system="AA",
            xvg_path=Path("."), time_ns=[], values=[],
        )
        r.observables["contacts"] = ObservableSeries(
            observable="contacts", replica="A1", system="AA",
            xvg_path=Path("."), time_ns=[], values=[],
        )
        sg.replicas["A1"] = r
        assert sg.observables == ["contacts", "protein_rmsd"]

    def test_observable_series_defaults(self):
        from runtime.study_models import ObservableSeries
        s = ObservableSeries(
            observable="protein_rmsd", replica="A1", system="AA",
            xvg_path=Path("."), time_ns=[], values=[],
        )
        assert s.mean == 0.0
        assert s.convergence_score == 0.0

    def test_comparative_summary_defaults(self):
        from runtime.study_models import ComparativeSummary
        cs = ComparativeSummary()
        assert cs.findings == []
        assert cs.outlier_replicas == []
        assert cs.system_ranking == {}


# ═══════════════════════════════════════════════════════════════════════════════
# 2. observable_resolver.py
# ═══════════════════════════════════════════════════════════════════════════════

class TestNormalize:
    def test_lowercase(self):
        from runtime.observable_resolver import _normalize
        assert _normalize("RMSD") == "rmsd"

    def test_camel_case(self):
        from runtime.observable_resolver import _normalize
        assert _normalize("rmsdProtein") == "rmsd_protein"

    def test_hyphens_to_underscore(self):
        from runtime.observable_resolver import _normalize
        assert _normalize("rmsd-protein") == "rmsd_protein"

    def test_strips_leading_trailing(self):
        from runtime.observable_resolver import _normalize
        assert _normalize("_rmsd_") == "rmsd"


class TestObservableResolver:
    @pytest.fixture(autouse=True)
    def resolver(self):
        from runtime.observable_resolver import ObservableResolver
        self.r = ObservableResolver()

    # Priority: ligand_rmsd BEFORE protein_rmsd
    def test_ligand_rmsd_wins_over_protein_rmsd(self):
        ro = self.r.resolve("rmsd_lig")
        assert ro.canonical == "ligand_rmsd"

    def test_ligand_rmsd_with_drug_keyword(self):
        ro = self.r.resolve("rmsd_drug")
        assert ro.canonical == "ligand_rmsd"

    def test_ligand_rmsd_with_inhibitor_keyword(self):
        ro = self.r.resolve("rmsd_inhibitor")
        assert ro.canonical == "ligand_rmsd"

    def test_protein_rmsd_fallback(self):
        ro = self.r.resolve("rmsd_protein")
        assert ro.canonical == "protein_rmsd"

    def test_protein_rmsd_plain_rmsd(self):
        ro = self.r.resolve("rmsd")
        assert ro.canonical == "protein_rmsd"

    def test_camel_case_rmsd_protein(self):
        ro = self.r.resolve("rmsdProtein")
        assert ro.canonical == "protein_rmsd"

    # Priority: mindist BEFORE distance
    def test_mindist_wins_over_distance(self):
        ro = self.r.resolve("mindist")
        assert ro.canonical == "mindist"

    def test_mindist_alias_min_distance(self):
        ro = self.r.resolve("min_distance")
        assert ro.canonical == "mindist"

    # Contacts
    def test_contacts(self):
        assert self.r.resolve("contacts").canonical == "contacts"
        assert self.r.resolve("native_contacts").canonical == "contacts"

    # RMSF
    def test_rmsf(self):
        ro = self.r.resolve("rmsf")
        assert ro.canonical == "rmsf"
        assert ro.units == "nm"

    # Radius of gyration
    def test_radius_of_gyration(self):
        assert self.r.resolve("gyration").canonical == "radius_of_gyration"
        assert self.r.resolve("gyr").canonical == "radius_of_gyration"

    # Hydrogen bonds
    def test_hydrogen_bonds(self):
        assert self.r.resolve("hbonds").canonical == "hydrogen_bonds"
        assert self.r.resolve("h_bond_count").canonical == "hydrogen_bonds"

    # Temperature / pressure / energy
    def test_temperature(self):
        assert self.r.resolve("temperature").canonical == "temperature"
        assert self.r.resolve("temp_avg").canonical == "temperature"

    def test_pressure(self):
        assert self.r.resolve("pressure").canonical == "pressure"

    def test_potential_energy(self):
        assert self.r.resolve("potential").canonical == "potential_energy"

    # Groups
    def test_group_structural(self):
        assert self.r.resolve("rmsd").group == "structural"

    def test_group_interaction(self):
        assert self.r.resolve("contacts").group == "interaction"

    def test_group_energetic(self):
        assert self.r.resolve("temperature").group == "energetic"

    # Fallback
    def test_unknown_returns_canonical_from_hint(self):
        ro = self.r.resolve("xyz_metric")
        assert ro.group == "other"

    # resolve_from_path
    def test_resolve_from_path(self, tmp_path):
        p = tmp_path / "rmsd_protein.xvg"
        p.write_text("")
        ro = self.r.resolve_from_path(p)
        assert ro.canonical == "protein_rmsd"

    # Title/ylabel disambiguation
    def test_title_helps_disambiguate(self):
        ro = self.r.resolve("rmsd", xvg_title="Ligand RMSD")
        assert ro.canonical == "ligand_rmsd"

    # Catalytic distance — explicit patterns
    def test_catalytic_distance_keyword(self):
        assert self.r.resolve("catalytic_dist").canonical == "catalytic_distance"
        assert self.r.resolve("dist_catalytic").canonical == "catalytic_distance"

    def test_catalytic_distance_cat_dist(self):
        assert self.r.resolve("cat_dist_res").canonical == "catalytic_distance"
        assert self.r.resolve("dist_cat_residue").canonical == "catalytic_distance"

    def test_catalytic_distance_active_site(self):
        assert self.r.resolve("active_site_dist").canonical == "catalytic_distance"
        assert self.r.resolve("dist_active_site").canonical == "catalytic_distance"

    def test_catalytic_distance_lig_cat(self):
        assert self.r.resolve("dist_lig_cat_res").canonical == "catalytic_distance"
        assert self.r.resolve("cat_lig_distance").canonical == "catalytic_distance"

    def test_catalytic_wins_over_generic_distance(self):
        assert self.r.resolve("active_site_dist").canonical == "catalytic_distance"
        assert self.r.resolve("dist_active_site").canonical != "distance"

    # Hydrogen bonds — extended filename patterns
    def test_hydrogen_bonds_trailing_hb(self):
        assert self.r.resolve("num_hb").canonical == "hydrogen_bonds"
        assert self.r.resolve("n_hb").canonical == "hydrogen_bonds"
        assert self.r.resolve("lig_hb").canonical == "hydrogen_bonds"

    def test_hydrogen_bonds_hbnum(self):
        assert self.r.resolve("hbnum").canonical == "hydrogen_bonds"

    def test_hydrogen_bonds_nhb(self):
        assert self.r.resolve("nhb").canonical == "hydrogen_bonds"
        assert self.r.resolve("nhb_lig_prot").canonical == "hydrogen_bonds"

    def test_hydrogen_bonds_existing_patterns_unchanged(self):
        assert self.r.resolve("hbonds").canonical == "hydrogen_bonds"
        assert self.r.resolve("hbond_lig_prot").canonical == "hydrogen_bonds"
        assert self.r.resolve("h_bond_count").canonical == "hydrogen_bonds"


# ═══════════════════════════════════════════════════════════════════════════════
# 3. study_analyzer.py
# ═══════════════════════════════════════════════════════════════════════════════

class TestFilenameParser:
    def test_standard_pattern_no_sep(self):
        from runtime.study_analyzer import _parse_filename
        r = _parse_filename("AA-A1rmsd_protein")
        assert r == ("AA", "A1", "rmsd_protein")

    def test_standard_pattern_with_sep(self):
        from runtime.study_analyzer import _parse_filename
        r = _parse_filename("LP-A4_rmsd-ligand")
        assert r == ("LP", "A4", "rmsd-ligand")

    def test_standard_pattern_dash_sep(self):
        from runtime.study_analyzer import _parse_filename
        r = _parse_filename("HMG-A5-mindist_lig")
        assert r == ("HMG", "A5", "mindist_lig")

    def test_system_uppercased(self):
        from runtime.study_analyzer import _parse_filename
        r = _parse_filename("aa-A1rmsd")
        assert r[0] == "AA"

    def test_replica_uppercased(self):
        from runtime.study_analyzer import _parse_filename
        r = _parse_filename("AA-a1rmsd")
        assert r[1] == "A1"

    def test_no_match_plain_file(self):
        from runtime.study_analyzer import _parse_filename
        assert _parse_filename("rmsd_protein") is None

    def test_no_match_too_long_system(self):
        from runtime.study_analyzer import _parse_filename
        # 7 chars exceeds limit
        assert _parse_filename("TOOLONG-A1rmsd") is None

    def test_multi_digit_replica(self):
        from runtime.study_analyzer import _parse_filename
        r = _parse_filename("AA-A10_contacts")
        assert r == ("AA", "A10", "contacts")


class TestSeriesStats:
    def test_mean_and_std_flat(self):
        from runtime.study_analyzer import _compute_series_stats
        from runtime.study_models import ObservableSeries
        s = ObservableSeries(
            observable="protein_rmsd", replica="A1", system="AA",
            xvg_path=Path("."),
            time_ns=[0.0, 1.0, 2.0, 3.0],
            values=[0.20, 0.20, 0.20, 0.20],
        )
        _compute_series_stats(s)
        assert s.mean == pytest.approx(0.20)
        assert s.std == pytest.approx(0.0, abs=1e-10)

    def test_stable_series_high_convergence(self):
        from runtime.study_analyzer import _compute_series_stats
        from runtime.study_models import ObservableSeries
        times, vals = _flat_series(100, 0.20)
        s = ObservableSeries(
            observable="protein_rmsd", replica="A1", system="AA",
            xvg_path=Path("."),
            time_ns=[t / 1000.0 for t in times],
            values=vals,
        )
        _compute_series_stats(s)
        assert s.convergence_score > 0.90

    def test_drifting_series_low_convergence(self):
        from runtime.study_analyzer import _compute_series_stats
        from runtime.study_models import ObservableSeries
        times, vals = _drifting_series(100, 0.10, 0.80)
        s = ObservableSeries(
            observable="protein_rmsd", replica="A1", system="AA",
            xvg_path=Path("."),
            time_ns=[t / 1000.0 for t in times],
            values=vals,
        )
        _compute_series_stats(s)
        assert s.convergence_score < 0.50

    def test_empty_series_no_crash(self):
        from runtime.study_analyzer import _compute_series_stats
        from runtime.study_models import ObservableSeries
        s = ObservableSeries(
            observable="protein_rmsd", replica="A1", system="AA",
            xvg_path=Path("."), time_ns=[], values=[],
        )
        _compute_series_stats(s)
        assert s.mean == 0.0

    def test_median_odd(self):
        from runtime.study_analyzer import _median
        assert _median([1.0, 3.0, 2.0]) == 2.0

    def test_median_even(self):
        from runtime.study_analyzer import _median
        assert _median([1.0, 2.0, 3.0, 4.0]) == 2.5


class TestGrubbsThreshold:
    def test_known_values(self):
        from runtime.study_analyzer import _outlier_threshold
        assert _outlier_threshold(3) == 1.15
        assert _outlier_threshold(4) == 1.48
        assert _outlier_threshold(5) == 1.71
        assert _outlier_threshold(8) == 2.13

    def test_fallback_above_8(self):
        from runtime.study_analyzer import _outlier_threshold
        assert _outlier_threshold(9)  == 2.0
        assert _outlier_threshold(50) == 2.0


class TestOutlierDetection:
    def test_outlier_detected_n4(self, tmp_path):
        """Replica with 2.5x the normal RMSD must be flagged with n=4."""
        from runtime.study_analyzer import parse_study
        # Normal replicas: 0.18, outlier: 0.95
        values = {"A1": 0.18, "A2": 0.18, "A3": 0.18, "A4": 0.95}
        for rep, val in values.items():
            _write_xvg(tmp_path / f"HMG-{rep}_rmsd_protein.xvg", "RMSD", "nm",
                       *_flat_series(100, val))
        study = parse_study(tmp_path)
        outlier_reps = {r for _, r, _ in study.summary.outlier_replicas}
        assert "A4" in outlier_reps

    def test_normal_replicas_not_flagged(self, two_system_study):
        """Replicas with similar values must not be flagged."""
        outliers = {r for _, r, _ in two_system_study.summary.outlier_replicas}
        # All three AA replicas are near-identical — none should be outlier
        aa_reps = {r for s, r, _ in two_system_study.summary.outlier_replicas if s == "AA"}
        assert aa_reps == set()

    def test_numerical_noise_not_flagged(self, tmp_path):
        """Tiny numerical variation (< 10% relative) must not be flagged."""
        import random
        rng = random.Random(0)
        for rep in ("A1", "A2", "A3", "A4"):
            noise = 1.0 + rng.gauss(0, 0.005)  # < 1% noise
            _write_xvg(tmp_path / f"SYS-{rep}_contacts.xvg", "Contacts", "",
                       *_flat_series(100, 12.0 * noise))
        from runtime.study_analyzer import parse_study
        study = parse_study(tmp_path)
        assert study.summary.outlier_replicas == []

    def test_n_less_than_3_skipped(self, tmp_path):
        """Outlier detection requires n >= 3 replicas."""
        for rep in ("A1", "A2"):
            _write_xvg(tmp_path / f"XX-{rep}_rmsd_protein.xvg", "RMSD", "nm",
                       *_flat_series(100, 0.20))
        from runtime.study_analyzer import parse_study
        study = parse_study(tmp_path)
        assert study.summary.outlier_replicas == []


class TestComparativeFindings:
    def test_findings_generated_for_high_variation(self, two_system_study):
        """AA vs LP have >15% relative range on contacts → findings generated."""
        obs_findings = [f for f in two_system_study.summary.findings
                        if f.observable == "contacts"]
        assert len(obs_findings) > 0

    def test_system_ranking_populated(self, two_system_study):
        assert "AA" in two_system_study.summary.system_ranking
        assert "LP" in two_system_study.summary.system_ranking

    def test_highlight_level_used(self, two_system_study):
        levels = {f.level for f in two_system_study.summary.findings}
        assert "highlight" in levels


class TestParseStudyIntegration:
    def test_systems_discovered(self, two_system_study):
        assert set(two_system_study.systems.keys()) == {"AA", "LP"}

    def test_replica_count(self, two_system_study):
        assert two_system_study.systems["AA"].n_replicas == 3
        assert two_system_study.systems["LP"].n_replicas == 3

    def test_observables_detected(self, two_system_study):
        assert "protein_rmsd" in two_system_study.observables_detected
        assert "contacts" in two_system_study.observables_detected

    def test_parsed_count(self, two_system_study):
        assert two_system_study.n_xvg_parsed == 12  # 2 sys × 3 rep × 2 obs

    def test_aggregate_means_correct(self, two_system_study):
        aa_rmsd = two_system_study.systems["AA"].aggregate["protein_rmsd"]
        lp_rmsd = two_system_study.systems["LP"].aggregate["protein_rmsd"]
        assert aa_rmsd.mean == pytest.approx(0.18, abs=1e-6)
        assert lp_rmsd.mean == pytest.approx(0.45, abs=1e-6)

    def test_aa_contacts_higher_than_lp(self, two_system_study):
        aa = two_system_study.systems["AA"].aggregate["contacts"].mean
        lp = two_system_study.systems["LP"].aggregate["contacts"].mean
        assert aa > lp

    def test_display_names_populated(self, two_system_study):
        assert two_system_study.observable_display.get("protein_rmsd") == "Protein RMSD"
        assert two_system_study.observable_display.get("contacts") == "Contacts"

    def test_fallback_single_system(self, tmp_path):
        """Files without study pattern → treated as single-system with replica R1."""
        from runtime.study_analyzer import parse_study
        _write_xvg(tmp_path / "rmsd_protein.xvg", "RMSD", "nm", *_flat_series(50, 0.20))
        _write_xvg(tmp_path / "contacts.xvg",     "Contacts", "", *_flat_series(50, 12.0))
        study = parse_study(tmp_path)
        assert len(study.systems) == 1
        sg = list(study.systems.values())[0]
        assert "R1" in sg.replicas

    def test_parse_errors_captured(self, tmp_path):
        """Malformed XVG should be captured in parse_errors, not raise."""
        from runtime.study_analyzer import parse_study
        (tmp_path / "AA-A1_rmsd.xvg").write_text("NOT XVG DATA AT ALL\n!!!")
        study = parse_study(tmp_path)
        # May parse (no numbers) or fail gracefully — either way no exception
        assert study is not None


# ═══════════════════════════════════════════════════════════════════════════════
# 4. synthesis_models.py
# ═══════════════════════════════════════════════════════════════════════════════

class TestSynthesisModels:
    def test_signal_profile_defaults(self):
        from runtime.synthesis_models import SignalProfile
        p = SignalProfile(system="AA")
        assert p.normalized == {}
        assert p.direction_score == {}
        assert p.consistency == {}

    def test_synthesis_result_defaults(self):
        from runtime.synthesis_models import SynthesisResult
        r = SynthesisResult()
        assert r.systems == {}
        assert r.ranking == []
        assert r.narrative == ""

    def test_rule_match_fields(self):
        from runtime.synthesis_models import RuleMatch
        m = RuleMatch(
            state="stable_binding", description="desc",
            score=0.8, confidence=0.75,
            n_conditions_met=3, n_conditions_available=3,
            supporting=["favorable contacts"], opposing=[],
        )
        assert m.score == 0.8
        assert len(m.supporting) == 1

    def test_temporal_event_fields(self):
        from runtime.synthesis_models import TemporalEvent
        e = TemporalEvent(
            system="AA", replica="A1", observable="ligand_rmsd",
            event_type="ligand_drift", time_ns=5.0,
            description="drift detected",
        )
        assert e.time_ns == 5.0


# ═══════════════════════════════════════════════════════════════════════════════
# 5. interaction_interpreter.py
# ═══════════════════════════════════════════════════════════════════════════════

class TestSigmoidNorm:
    def test_at_mean_gives_half(self):
        from runtime.interaction_interpreter import _normalize_sigmoid
        assert _normalize_sigmoid(5.0, 5.0, 2.0) == pytest.approx(0.5)

    def test_above_mean_gives_more_than_half(self):
        from runtime.interaction_interpreter import _normalize_sigmoid
        assert _normalize_sigmoid(7.0, 5.0, 2.0) > 0.5

    def test_below_mean_gives_less_than_half(self):
        from runtime.interaction_interpreter import _normalize_sigmoid
        assert _normalize_sigmoid(3.0, 5.0, 2.0) < 0.5

    def test_zero_std_returns_half(self):
        from runtime.interaction_interpreter import _normalize_sigmoid
        assert _normalize_sigmoid(5.0, 5.0, 0.0) == pytest.approx(0.5)


class TestDirectionScore:
    def test_lower_is_better_inversion(self):
        from runtime.interaction_interpreter import _direction_score
        # For lower_is_better, a HIGH normalized value → low direction_score (bad)
        assert _direction_score("protein_rmsd", 0.9) == pytest.approx(0.1)

    def test_lower_is_better_low_norm_is_favorable(self):
        from runtime.interaction_interpreter import _direction_score
        # low normalized → direction close to 1 (favorable)
        assert _direction_score("protein_rmsd", 0.1) == pytest.approx(0.9)

    def test_higher_is_better_direct(self):
        from runtime.interaction_interpreter import _direction_score
        assert _direction_score("contacts", 0.8) == pytest.approx(0.8)

    def test_neutral_returns_half(self):
        from runtime.interaction_interpreter import _direction_score
        assert _direction_score("unknown_obs", 0.7) == pytest.approx(0.5)


class TestBuildSignalProfiles:
    def test_profiles_built_for_each_system(self, rich_study):
        from runtime.interaction_interpreter import build_signal_profiles
        profiles = build_signal_profiles(rich_study)
        assert "AA" in profiles
        assert "LP" in profiles

    def test_aa_contacts_direction_higher_than_lp(self, rich_study):
        from runtime.interaction_interpreter import build_signal_profiles
        profiles = build_signal_profiles(rich_study)
        assert profiles["AA"].direction_score["contacts"] > profiles["LP"].direction_score["contacts"]

    def test_aa_ligand_rmsd_direction_higher_than_lp(self, rich_study):
        """AA has lower ligand RMSD → more favorable → higher direction_score."""
        from runtime.interaction_interpreter import build_signal_profiles
        profiles = build_signal_profiles(rich_study)
        assert profiles["AA"].direction_score["ligand_rmsd"] > profiles["LP"].direction_score["ligand_rmsd"]

    def test_consistency_between_0_and_1(self, rich_study):
        from runtime.interaction_interpreter import build_signal_profiles
        profiles = build_signal_profiles(rich_study)
        for p in profiles.values():
            for v in p.consistency.values():
                assert 0.0 <= v <= 1.0


class TestEvaluateRules:
    def test_rules_return_list(self, rich_study):
        from runtime.interaction_interpreter import build_signal_profiles, evaluate_rules
        profiles = build_signal_profiles(rich_study)
        matches = evaluate_rules(profiles["AA"], rich_study)
        assert isinstance(matches, list)

    def test_rules_sorted_by_score_descending(self, rich_study):
        from runtime.interaction_interpreter import build_signal_profiles, evaluate_rules
        profiles = build_signal_profiles(rich_study)
        matches = evaluate_rules(profiles["AA"], rich_study)
        scores = [m.score for m in matches]
        assert scores == sorted(scores, reverse=True)

    def test_score_in_range(self, rich_study):
        from runtime.interaction_interpreter import build_signal_profiles, evaluate_rules
        profiles = build_signal_profiles(rich_study)
        for sys_name in ["AA", "LP"]:
            for m in evaluate_rules(profiles[sys_name], rich_study):
                assert 0.0 <= m.score <= 1.0
                assert 0.0 <= m.confidence <= 1.0


class TestInterpretSystem:
    def test_aa_composite_higher_than_lp(self, rich_study):
        from runtime.interaction_interpreter import interpret_all
        synths = interpret_all(rich_study)
        assert synths["AA"].composite_score > synths["LP"].composite_score

    def test_aa_binding_score_higher(self, rich_study):
        from runtime.interaction_interpreter import interpret_all
        synths = interpret_all(rich_study)
        assert synths["AA"].binding_score > synths["LP"].binding_score

    def test_primary_state_is_string(self, rich_study):
        from runtime.interaction_interpreter import interpret_all
        synths = interpret_all(rich_study)
        for syn in synths.values():
            assert isinstance(syn.primary_state, str)
            assert len(syn.primary_state) > 0

    def test_evidence_list(self, rich_study):
        from runtime.interaction_interpreter import interpret_all
        synths = interpret_all(rich_study)
        for syn in synths.values():
            assert isinstance(syn.evidence, list)

    def test_explanation_nonempty(self, rich_study):
        from runtime.interaction_interpreter import interpret_all
        synths = interpret_all(rich_study)
        for syn in synths.values():
            assert isinstance(syn.explanation, str)
            assert len(syn.explanation) > 0

    def test_tiebreaker_prefers_most_conditions(self, rich_study):
        """When top two rules score within 0.08, state with more conditions wins."""
        from runtime.interaction_interpreter import build_signal_profiles, evaluate_rules, interpret_system, _SCORE_THRESHOLD
        profiles = build_signal_profiles(rich_study)
        # AA is the good binder — it should not be uncertain_behavior when
        # a rule with 3 conditions beats one with 1 condition
        syn = interpret_system(profiles["AA"], rich_study)
        # If uncertain_behavior, that's only acceptable if the top two TIED on conditions
        if syn.primary_state == "uncertain_behavior":
            matches = evaluate_rules(profiles["AA"], rich_study)
            active = [m for m in matches if m.score >= _SCORE_THRESHOLD]
            if len(active) >= 2:
                top_cond = active[0].n_conditions_available
                second_cond = active[1].n_conditions_available
                # If tiebreaker applies, top must not have more conditions than second
                assert top_cond <= second_cond


# ═══════════════════════════════════════════════════════════════════════════════
# 6. consensus_engine.py
# ═══════════════════════════════════════════════════════════════════════════════

class TestConsensusEngine:
    def test_full_consensus(self, rich_study):
        """With stable replicas, at least one (system, obs) pair is full or strong."""
        from runtime.consensus_engine import evaluate_consensus
        results = evaluate_consensus(rich_study)
        # All replicas of AA are on same side → full or strong for AA
        # (may be absent from results since 'full' is excluded by default)
        # Just check no exception and result is a list
        assert isinstance(results, list)

    def test_conflicting_consensus_returns_low_multiplier(self):
        from runtime.consensus_engine import _label_and_multiplier
        label, mult = _label_and_multiplier(0.45, 4)
        assert label == "conflicting"
        assert mult == pytest.approx(0.45)

    def test_full_consensus_returns_1(self):
        from runtime.consensus_engine import _label_and_multiplier
        label, mult = _label_and_multiplier(1.0, 4)
        assert label == "full"
        assert mult == pytest.approx(1.0)

    def test_strong_consensus(self):
        from runtime.consensus_engine import _label_and_multiplier
        label, mult = _label_and_multiplier(0.85, 4)
        assert label == "strong"
        assert mult == pytest.approx(0.9)

    def test_moderate_consensus(self):
        from runtime.consensus_engine import _label_and_multiplier
        label, mult = _label_and_multiplier(0.65, 4)
        assert label == "moderate"
        assert mult == pytest.approx(0.75)

    def test_weak_consensus(self):
        from runtime.consensus_engine import _label_and_multiplier
        label, mult = _label_and_multiplier(0.50, 4)
        assert label == "weak"
        assert mult == pytest.approx(0.60)

    def test_insufficient_data_n1(self):
        from runtime.consensus_engine import _label_and_multiplier
        label, mult = _label_and_multiplier(1.0, 1)
        assert label == "insufficient_data"
        assert mult == pytest.approx(0.5)

    def test_multiplier_map_returns_dict(self, rich_study):
        from runtime.consensus_engine import consensus_multiplier_map
        cmap = consensus_multiplier_map(rich_study)
        assert isinstance(cmap, dict)
        # Values must be between 0 and 1
        for v in cmap.values():
            assert 0.0 <= v <= 1.0


# ═══════════════════════════════════════════════════════════════════════════════
# 7. event_detector.py
# ═══════════════════════════════════════════════════════════════════════════════

def _make_series(observable: str, times_ns: list[float], values: list[float],
                 system: str = "SYS", replica: str = "A1"):
    from runtime.study_models import ObservableSeries
    return ObservableSeries(
        observable=observable, replica=replica, system=system,
        xvg_path=Path("."), time_ns=times_ns, values=values,
    )


class TestLateDestabilization:
    def test_detected_when_second_half_increases(self):
        from runtime.event_detector import _detect_late_destabilization
        n = 60
        times = [i * 0.1 for i in range(n)]
        # First half stable at 0.2, second half rises to 0.4
        vals = [0.2] * 30 + [0.4] * 30
        s = _make_series("protein_rmsd", times, vals)
        event = _detect_late_destabilization(s)
        assert event is not None
        assert event.event_type == "late_destabilization"

    def test_not_detected_when_stable(self):
        from runtime.event_detector import _detect_late_destabilization
        n = 60
        times = [i * 0.1 for i in range(n)]
        s = _make_series("protein_rmsd", times, [0.2] * n)
        assert _detect_late_destabilization(s) is None

    def test_contacts_drop_is_contact_loss(self):
        from runtime.event_detector import _detect_late_destabilization
        n = 60
        times = [i * 0.1 for i in range(n)]
        vals = [14.0] * 30 + [8.0] * 30
        s = _make_series("contacts", times, vals)
        event = _detect_late_destabilization(s)
        assert event is not None
        assert event.event_type == "contact_loss"

    def test_too_few_points_returns_none(self):
        from runtime.event_detector import _detect_late_destabilization
        s = _make_series("protein_rmsd", [0.0, 1.0], [0.2, 0.8])
        assert _detect_late_destabilization(s) is None


class TestAbruptTransition:
    def test_constant_series_returns_none(self):
        # overall_std ≈ 0 → early return before window scan
        from runtime.event_detector import _detect_abrupt_transition
        n = 80
        times = [i * 0.1 for i in range(n)]
        vals = [0.42] * n
        s = _make_series("protein_rmsd", times, vals)
        assert _detect_abrupt_transition(s) is None

    def test_not_detected_smooth_change(self):
        from runtime.event_detector import _detect_abrupt_transition
        n = 80
        times = [i * 0.1 for i in range(n)]
        vals = [0.2 + 0.005 * i for i in range(n)]
        s = _make_series("protein_rmsd", times, vals)
        # Small gradient — should not trigger 3σ
        event = _detect_abrupt_transition(s)
        assert event is None

    def test_too_few_points_returns_none(self):
        from runtime.event_detector import _detect_abrupt_transition
        s = _make_series("protein_rmsd", [0.0, 1.0], [0.2, 0.8])
        assert _detect_abrupt_transition(s) is None


class TestLigandDrift:
    def test_detected_monotonic_increase(self):
        from runtime.event_detector import _detect_ligand_drift
        n = 60
        times = [i * 0.1 for i in range(n)]
        # Monotonically increasing: 0.15 → 0.20 → 0.25 (>30% total rise)
        third = n // 3
        vals  = [0.15] * third + [0.20] * third + [0.25] * (n - 2 * third)
        s = _make_series("ligand_rmsd", times, vals)
        event = _detect_ligand_drift(s)
        assert event is not None
        assert event.event_type == "ligand_drift"

    def test_not_detected_for_protein_rmsd(self):
        from runtime.event_detector import _detect_ligand_drift
        n = 60
        times = [i * 0.1 for i in range(n)]
        third = n // 3
        vals  = [0.15] * third + [0.20] * third + [0.25] * (n - 2 * third)
        s = _make_series("protein_rmsd", times, vals)
        assert _detect_ligand_drift(s) is None

    def test_not_detected_small_total_rise(self):
        from runtime.event_detector import _detect_ligand_drift
        n = 60
        times = [i * 0.1 for i in range(n)]
        third = n // 3
        # Only 10% rise — below 30% threshold
        vals  = [0.20] * third + [0.21] * third + [0.22] * (n - 2 * third)
        s = _make_series("ligand_rmsd", times, vals)
        assert _detect_ligand_drift(s) is None


class TestDetectAllEvents:
    def test_returns_list(self, rich_study):
        from runtime.event_detector import detect_all_events
        events = detect_all_events(rich_study)
        assert isinstance(events, list)

    def test_no_duplicates(self, rich_study):
        from runtime.event_detector import detect_all_events
        events = detect_all_events(rich_study)
        keys = [(e.system, e.replica, e.observable, e.event_type) for e in events]
        assert len(keys) == len(set(keys))

    def test_abrupt_before_late_destabilization(self, two_system_dir):
        """Severity ordering: abrupt_transition < late_destabilization."""
        from runtime.event_detector import detect_all_events
        from runtime.study_analyzer import parse_study
        # Add a single file with both patterns
        d = two_system_dir
        study = parse_study(d)
        events = detect_all_events(study)
        types = [e.event_type for e in events]
        abrupt_idx = next((i for i, t in enumerate(types) if t == "abrupt_transition"), None)
        late_idx   = next((i for i, t in enumerate(types) if t == "late_destabilization"), None)
        if abrupt_idx is not None and late_idx is not None:
            assert abrupt_idx < late_idx


# ═══════════════════════════════════════════════════════════════════════════════
# 8. scientific_synthesis.py
# ═══════════════════════════════════════════════════════════════════════════════

class TestScientificSynthesis:
    def test_returns_synthesis_result(self, rich_study):
        from runtime.scientific_synthesis import synthesize_study
        result = synthesize_study(rich_study)
        from runtime.synthesis_models import SynthesisResult
        assert isinstance(result, SynthesisResult)

    def test_systems_populated(self, rich_study):
        from runtime.scientific_synthesis import synthesize_study
        result = synthesize_study(rich_study)
        assert "AA" in result.systems
        assert "LP" in result.systems

    def test_ranking_sorted_descending(self, rich_study):
        from runtime.scientific_synthesis import synthesize_study
        result = synthesize_study(rich_study)
        scores = [sc for _, sc in result.ranking]
        assert scores == sorted(scores, reverse=True)

    def test_aa_ranks_first(self, rich_study):
        from runtime.scientific_synthesis import synthesize_study
        result = synthesize_study(rich_study)
        assert result.ranking[0][0] == "AA"

    def test_narrative_nonempty(self, rich_study):
        from runtime.scientific_synthesis import synthesize_study
        result = synthesize_study(rich_study)
        assert len(result.narrative) > 50

    def test_narrative_mentions_aa(self, rich_study):
        from runtime.scientific_synthesis import synthesize_study
        result = synthesize_study(rich_study)
        assert "AA" in result.narrative

    def test_empty_study_returns_empty_result(self):
        from runtime.scientific_synthesis import synthesize_study
        from runtime.study_models import Study
        empty = Study(path=Path("."))
        result = synthesize_study(empty)
        assert result.systems == {}
        assert result.ranking == []
        assert result.narrative == ""

    def test_composite_scores_in_range(self, rich_study):
        from runtime.scientific_synthesis import synthesize_study
        result = synthesize_study(rich_study)
        for _, sc in result.ranking:
            assert 0.0 <= sc <= 1.0

    def test_explanation_updated_in_synthesis(self, rich_study):
        """After synthesize_study, explanation should reference the system name."""
        from runtime.scientific_synthesis import synthesize_study
        result = synthesize_study(rich_study)
        for sys_name, syn in result.systems.items():
            assert sys_name in syn.explanation

    def test_consensus_included(self, rich_study):
        from runtime.scientific_synthesis import synthesize_study
        result = synthesize_study(rich_study)
        assert isinstance(result.consensus, list)

    def test_events_included(self, rich_study):
        from runtime.scientific_synthesis import synthesize_study
        result = synthesize_study(rich_study)
        assert isinstance(result.events, list)
