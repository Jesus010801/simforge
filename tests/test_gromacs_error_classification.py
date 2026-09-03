# tests/test_gromacs_error_classification.py
"""
Phase 11 — GROMACS error classification tests.

Verifies that:
1. TOPOLOGY_PREPROCESSING_ERROR exists in ErrorCategory
2. OXT/terminus signal matches the expected patterns
3. The planner returns a non-retryable plan for OXT errors
4. The plan does not appear in the retryable categories
5. DiagnosisResult.category is TOPOLOGY_PREPROCESSING_ERROR for OXT logs
6. The plan sets is_applicable=False and max_retries=0
"""
from __future__ import annotations

import datetime
from pathlib import Path

import pytest


# ═══════════════════════════════════════════════════════════════════════════════
# 1. ErrorCategory has TOPOLOGY_PREPROCESSING_ERROR
# ═══════════════════════════════════════════════════════════════════════════════

class TestErrorCategoryEnum:

    def test_topology_preprocessing_error_exists(self):
        from executors.remediation_models import ErrorCategory
        assert hasattr(ErrorCategory, "TOPOLOGY_PREPROCESSING_ERROR"), (
            "ErrorCategory must have TOPOLOGY_PREPROCESSING_ERROR member"
        )

    def test_topology_preprocessing_error_value(self):
        from executors.remediation_models import ErrorCategory
        assert ErrorCategory.TOPOLOGY_PREPROCESSING_ERROR.value == "topology_preprocessing_error"

    def test_topology_preprocessing_error_is_str_enum(self):
        from executors.remediation_models import ErrorCategory
        assert isinstance(ErrorCategory.TOPOLOGY_PREPROCESSING_ERROR, str)


# ═══════════════════════════════════════════════════════════════════════════════
# 2. OXT signal matches expected log lines
# ═══════════════════════════════════════════════════════════════════════════════

_OXT_LOG = (
    "Fatal error:\n"
    "Atom OXT in residue LYS 439 was not found in rtp entry LYSH with 22 atoms\n"
    "while sorting atoms\n"
)

_OXT_LOG_2 = (
    "Fatal error:\n"
    "was not found in rtp entry CTHR with 18 atoms\n"
)

_LINCS_LOG = (
    "LINCS warning: step 1000\n"
    "relative constraint deviation after LINCS: 0.012\n"
)

_NORMAL_LOG = "Step 1000: Epot = -123456.789 kJ/mol\n"


class TestOXTSignalDetection:

    def _detect(self, text: str):
        from executors.signal_detector import _SIGNAL_PATTERNS
        for signal in _SIGNAL_PATTERNS:
            hit, line = signal.match(text)
            if hit:
                return signal
        return None

    def test_oxt_log_matches_topology_error(self):
        signal = self._detect(_OXT_LOG)
        from executors.remediation_models import ErrorCategory
        assert signal is not None, "OXT log must match a signal"
        assert signal.category == ErrorCategory.TOPOLOGY_PREPROCESSING_ERROR

    def test_rtp_entry_log_matches_topology_error(self):
        signal = self._detect(_OXT_LOG_2)
        from executors.remediation_models import ErrorCategory
        assert signal is not None
        assert signal.category == ErrorCategory.TOPOLOGY_PREPROCESSING_ERROR

    def test_lincs_log_does_not_match_topology_error(self):
        signal = self._detect(_LINCS_LOG)
        from executors.remediation_models import ErrorCategory
        assert signal is not None
        assert signal.category != ErrorCategory.TOPOLOGY_PREPROCESSING_ERROR

    def test_normal_log_does_not_match_topology_error(self):
        signal = self._detect(_NORMAL_LOG)
        if signal:
            from executors.remediation_models import ErrorCategory
            assert signal.category != ErrorCategory.TOPOLOGY_PREPROCESSING_ERROR

    def test_oxt_signal_has_fatal_severity(self):
        from executors.signal_detector import _SIGNAL_PATTERNS
        from executors.remediation_models import ErrorCategory, ErrorSeverity
        oxt_signals = [
            s for s in _SIGNAL_PATTERNS
            if s.category == ErrorCategory.TOPOLOGY_PREPROCESSING_ERROR
        ]
        assert oxt_signals, "At least one TOPOLOGY_PREPROCESSING_ERROR signal must exist"
        assert all(s.severity == ErrorSeverity.FATAL for s in oxt_signals), (
            "TOPOLOGY_PREPROCESSING_ERROR signals must have FATAL severity"
        )

    def test_oxt_signal_high_confidence(self):
        from executors.signal_detector import _SIGNAL_PATTERNS
        from executors.remediation_models import ErrorCategory
        oxt_signals = [
            s for s in _SIGNAL_PATTERNS
            if s.category == ErrorCategory.TOPOLOGY_PREPROCESSING_ERROR
        ]
        assert all(s.confidence >= 0.90 for s in oxt_signals), (
            "TOPOLOGY_PREPROCESSING_ERROR signals must have confidence >= 0.90"
        )

    def test_oxt_signal_comes_before_nonzero_exit_in_list(self):
        """OXT signal must be checked before the generic NONZERO_EXIT fallback."""
        from executors.signal_detector import _SIGNAL_PATTERNS
        from executors.remediation_models import ErrorCategory
        categories = [s.category for s in _SIGNAL_PATTERNS]
        topology_idx = next(
            (i for i, c in enumerate(categories)
             if c == ErrorCategory.TOPOLOGY_PREPROCESSING_ERROR),
            None,
        )
        nonzero_idx = next(
            (i for i, c in enumerate(categories)
             if c == ErrorCategory.NONZERO_EXIT),
            None,
        )
        assert topology_idx is not None, "TOPOLOGY_PREPROCESSING_ERROR signal not in list"
        assert nonzero_idx is not None, "NONZERO_EXIT signal not in list"
        assert topology_idx < nonzero_idx, (
            "TOPOLOGY_PREPROCESSING_ERROR signal must appear before NONZERO_EXIT fallback"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# 3. Planner returns non-retryable plan for OXT error
# ═══════════════════════════════════════════════════════════════════════════════

class TestTopologyPreprocessingPlanner:

    def _make_diag(self, text: str, step_dir: Path):
        from executors.remediation_models import (
            DiagnosisResult, ErrorCategory, ErrorSeverity,
        )
        from executors.signal_detector import _SIGNAL_PATTERNS
        for signal in _SIGNAL_PATTERNS:
            hit, line = signal.match(text)
            if hit:
                return DiagnosisResult(
                    step_id       = "generate_topology",
                    step_dir      = str(step_dir),
                    category      = signal.category,
                    severity      = signal.severity,
                    confidence    = signal.confidence,
                    primary_signal = line,
                    explanation   = signal.explanation,
                )
        raise ValueError("No signal matched")

    def test_oxt_plan_is_not_applicable(self, tmp_path):
        from executors.signal_detector import AdaptiveReasoner
        diag = self._make_diag(_OXT_LOG, tmp_path)
        plan = AdaptiveReasoner().plan_remediation(diag, tmp_path)
        assert plan.is_applicable is False, (
            "Topology preprocessing error plan must not be applicable (non-retryable)"
        )

    def test_oxt_plan_max_retries_zero(self, tmp_path):
        from executors.signal_detector import AdaptiveReasoner
        diag = self._make_diag(_OXT_LOG, tmp_path)
        plan = AdaptiveReasoner().plan_remediation(diag, tmp_path)
        assert plan.max_retries == 0, (
            "Topology preprocessing error must have max_retries=0 (deterministic failure)"
        )

    def test_oxt_plan_requires_human(self, tmp_path):
        from executors.signal_detector import AdaptiveReasoner
        diag = self._make_diag(_OXT_LOG, tmp_path)
        plan = AdaptiveReasoner().plan_remediation(diag, tmp_path)
        assert plan.requires_human is True, (
            "Topology preprocessing error requires human intervention"
        )

    def test_oxt_plan_action_is_log_only(self, tmp_path):
        from executors.signal_detector import AdaptiveReasoner
        from executors.remediation_models import ActionType
        diag = self._make_diag(_OXT_LOG, tmp_path)
        plan = AdaptiveReasoner().plan_remediation(diag, tmp_path)
        assert len(plan.actions) >= 1
        assert plan.actions[0].action_type == ActionType.LOG_ONLY, (
            "Non-retryable error must use LOG_ONLY action (no file modifications)"
        )

    def test_oxt_plan_strategy_mentions_phase11(self, tmp_path):
        from executors.signal_detector import AdaptiveReasoner
        diag = self._make_diag(_OXT_LOG, tmp_path)
        plan = AdaptiveReasoner().plan_remediation(diag, tmp_path)
        combined = (plan.strategy + plan.expected_outcome).lower()
        assert "pdb" in combined or "topology" in combined, (
            "Plan strategy must mention the topology fix (use PDB input)"
        )


# ═══════════════════════════════════════════════════════════════════════════════
# 4. AdaptiveReasoner.diagnose classifies OXT correctly from a step record
# ═══════════════════════════════════════════════════════════════════════════════

class TestAdaptiveReasonerDiagnosis:

    def _make_record(self, stderr: str):
        from executors.execution_state import StepExecutionRecord
        record = StepExecutionRecord(
            step_id    = "generate_topology",
            step_dir   = "/fake/step",
            stdout     = "",
            stderr     = stderr,
            exit_code  = 1,
            started_at = datetime.datetime.now(),
            ended_at   = datetime.datetime.now(),
        )
        return record

    def test_diagnose_oxt_as_topology_preprocessing_error(self, tmp_path):
        from executors.signal_detector import AdaptiveReasoner
        from executors.remediation_models import ErrorCategory
        record = self._make_record(_OXT_LOG)
        diag = AdaptiveReasoner().diagnose(record, tmp_path)
        assert diag.category == ErrorCategory.TOPOLOGY_PREPROCESSING_ERROR, (
            f"Expected TOPOLOGY_PREPROCESSING_ERROR, got {diag.category}"
        )

    def test_diagnose_oxt_explanation_mentions_pdb(self, tmp_path):
        from executors.signal_detector import AdaptiveReasoner
        record = self._make_record(_OXT_LOG)
        diag = AdaptiveReasoner().diagnose(record, tmp_path)
        assert "PDB" in diag.explanation or "pdb2gmx" in diag.explanation.lower(), (
            "OXT diagnosis explanation must mention the PDB fix"
        )

    def test_diagnose_oxt_not_retryable_via_plan(self, tmp_path):
        from executors.signal_detector import AdaptiveReasoner
        record = self._make_record(_OXT_LOG)
        diag = AdaptiveReasoner().diagnose(record, tmp_path)
        plan = AdaptiveReasoner().plan_remediation(diag, tmp_path)
        assert plan.is_applicable is False
        assert plan.max_retries == 0

    def test_diagnose_lincs_is_still_recoverable(self, tmp_path):
        from executors.signal_detector import AdaptiveReasoner
        from executors.remediation_models import ErrorCategory
        record = self._make_record(_LINCS_LOG)
        diag = AdaptiveReasoner().diagnose(record, tmp_path)
        assert diag.category in (
            ErrorCategory.LINCS_WARNING,
            ErrorCategory.LINCS_FATAL,
        ), f"LINCS log must classify as LINCS error, got {diag.category}"

    def test_diagnose_rtp_entry_log_topology_error(self, tmp_path):
        from executors.signal_detector import AdaptiveReasoner
        from executors.remediation_models import ErrorCategory
        record = self._make_record(_OXT_LOG_2)
        diag = AdaptiveReasoner().diagnose(record, tmp_path)
        assert diag.category == ErrorCategory.TOPOLOGY_PREPROCESSING_ERROR
