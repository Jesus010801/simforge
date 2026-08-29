"""Multi-level study validation, run *before* any observable executes.

* **study-level**   – is the study coherent? (systems found, topologies present)
* **system-level**  – does this trajectory/topology pair look usable?
* **component-level** – are the semantic components identified?
* **observable-level** – for each requested analysis × system, are its
  ``required_components`` RESOLVED and its trajectory usable?

Statuses are explicit (:class:`ValidationState` + ``review_required``).  A
system with an unresolved peptide does not block ``rmsd-receptor`` but does
block ``rmsd-peptide-*`` — reported per system, never a global crash.
"""
from __future__ import annotations

from analysis.campaign.models import (
    CampaignWarning, ClassificationState, ComponentType, ObservableSystemStatus,
    Severity, StudyManifest, SystemRecord, ValidationReport, ValidationState,
)
from analysis.campaign.observables import registry

_COMPONENT_LABEL = {
    ComponentType.RECEPTOR: "receptor",
    ComponentType.PEPTIDE: "peptide",
    ComponentType.COMPLEX: "complex",
    ComponentType.LIGAND: "ligand",
    ComponentType.PROTEIN_PARTNER: "protein partner",
}


# Warnings that do NOT by themselves reduce a system's usability.
# (``topology_stage_unverified`` is deliberately NOT here — an unconfirmed
# production topology is a real caveat and should surface as a WARNING.)
_BENIGN_WARNING_CODES = {
    "no_replicate_marker", "component_inference", "assignment_invariant_repaired",
    "unbalanced_replicates", "dry_run",
}


def _system_usable(rec: SystemRecord) -> tuple[str, str]:
    """Central warning-propagation policy → (ValidationState, reason)."""
    if not rec.production_trajectory_paths:
        return ValidationState.INVALID, (
            "no production trajectory identified (only equilibration/preparation runs)")
    if not rec.topology_path:
        return ValidationState.INVALID, "no topology (.tpr/.top) associated"

    insp = rec.trajectory_inspection
    if insp and insp.inspected:
        if insp.atom_count_match is False:
            return ValidationState.INVALID, "trajectory/topology atom-count mismatch"

    # escalate from the system's own warnings
    worst = ValidationState.VALID
    reason = ""
    for w in rec.warnings:
        if w.code in _BENIGN_WARNING_CODES:
            continue
        if w.severity == Severity.ERROR:
            return ValidationState.INVALID, w.message
        if w.severity == Severity.REVIEW and worst != ValidationState.INVALID:
            worst, reason = ValidationState.WARNING, w.message
        elif w.severity == Severity.WARN and worst == ValidationState.VALID:
            worst, reason = ValidationState.WARNING, w.message

    if insp and insp.inspected:
        if insp.box_present is False:
            return ValidationState.WARNING, "some production frames have no box vectors"
        if insp.segment_overlap:
            return ValidationState.WARNING, "production trajectory segments overlap in time"
    if not rec.structure_path and worst == ValidationState.VALID:
        return ValidationState.WARNING, "no reference structure; limited semantic analysis"
    return worst, reason


def _component_resolved(rec: SystemRecord, ctype: str) -> bool:
    c = rec.component(ctype)
    return bool(c and c.classification_state == ClassificationState.RESOLVED)


def validate_study(
    manifest: StudyManifest, requested_analyses: list[str],
) -> ValidationReport:
    registry.ensure_loaded()
    report = ValidationReport()
    lines: list[str] = []

    # ── study level ────────────────────────────────────────────────────────
    if not manifest.systems:
        report.study_state = ValidationState.INVALID
        report.warnings.append(CampaignWarning(
            "no_systems", "no systems discovered", Severity.ERROR))
        report.lines = ["Study: INVALID — no systems discovered"]
        return report

    stats = manifest.stats()
    lines.append("Study:")
    lines.append(f"  systems:      {stats['n_systems']}")
    lines.append(f"  conditions:   {stats['n_conditions']}")
    lines.append(f"  partners:     {stats['n_partners']}  {stats['partners']}")
    lines.append(f"  water models: {stats['n_water_models']}  {stats['water_models']}")
    lines.append(f"  replicates:   {stats['replicate_counts']}"
                 f"  ({'balanced' if stats['balanced'] else 'UNBALANCED'})")

    # ── system level ──────────────────────────────────────────────────────
    n_valid = n_warn = n_invalid = 0
    recep_ok = pep_ok = memb_ok = 0
    for rec in manifest.systems:
        state, reason = _system_usable(rec)
        report.system_states[rec.system_id] = state
        # propagate every non-benign system warning into the study report
        for w in rec.warnings:
            if w.code in _BENIGN_WARNING_CODES:
                continue
            if w.severity in (Severity.WARN, Severity.REVIEW, Severity.ERROR):
                report.warnings.append(CampaignWarning(
                    w.code, f"{rec.system_id}: {w.message}", w.severity, scope=rec.system_id))
        if state == ValidationState.VALID:
            n_valid += 1
        elif state == ValidationState.WARNING:
            n_warn += 1
            report.warnings.append(CampaignWarning(
                "system_warning", f"{rec.system_id}: {reason}", Severity.WARN,
                scope=rec.system_id))
        else:
            n_invalid += 1
            report.warnings.append(CampaignWarning(
                "system_invalid", f"{rec.system_id}: {reason}", Severity.ERROR,
                scope=rec.system_id))
        if _component_resolved(rec, ComponentType.RECEPTOR):
            recep_ok += 1
        if _component_resolved(rec, ComponentType.PEPTIDE):
            pep_ok += 1
        if rec.membrane_present:
            memb_ok += 1

    n = len(manifest.systems)
    lines += [
        "",
        "Files:",
        f"  systems with topology:  {sum(1 for s in manifest.systems if s.topology_path)}/{n}",
        f"  systems with structure: {sum(1 for s in manifest.systems if s.structure_path)}/{n}",
        "",
        "Structural classification:",
        f"  receptor resolved: {recep_ok}/{n}",
        f"  peptide resolved:  {pep_ok}/{n}",
        f"  membrane present:  {memb_ok}/{n}",
        "",
        "System usability:",
        f"  valid: {n_valid}   warning: {n_warn}   invalid: {n_invalid}",
    ]

    # ── trajectory validation summary ────────────────────────────────────
    inspected = [s for s in manifest.systems
                 if s.trajectory_inspection and s.trajectory_inspection.inspected]
    if inspected:
        box_ok = sum(1 for s in inspected if s.trajectory_inspection.box_present)
        atom_ok = sum(1 for s in inspected if s.trajectory_inspection.atom_count_match)
        lines += [
            "",
            "Trajectory validation:",
            f"  inspected:            {len(inspected)}/{n}",
            f"  box present:          {box_ok}/{len(inspected)}",
            f"  atom-count matches:   {atom_ok}/{len(inspected)}",
        ]

    # ── observable level ────────────────────────────────────────────────
    if requested_analyses:
        lines += ["", "Requested analyses:"]
    for analysis_id in requested_analyses:
        if not registry.is_registered(analysis_id):
            report.warnings.append(CampaignWarning(
                "unknown_analysis",
                f"'{analysis_id}' is not a registered analysis. Valid: {registry.ids()}",
                Severity.ERROR))
            lines.append(f"  {analysis_id}: ERROR — unknown analysis id")
            continue
        spec = registry.get(analysis_id)
        n_ok = n_skip = n_review = 0
        for rec in manifest.systems:
            sys_state = report.system_states[rec.system_id]
            if sys_state == ValidationState.INVALID:
                report.observable_statuses.append(ObservableSystemStatus(
                    rec.system_id, analysis_id, ValidationState.INVALID,
                    "system is not usable"))
                n_skip += 1
                continue
            missing = [
                _COMPONENT_LABEL.get(ct, ct) for ct in spec.required_components
                if not _component_resolved(rec, ct)
            ]
            if missing:
                # is it review_required (ambiguous) or genuinely absent (skip)?
                any_ambiguous = any(
                    (rec.component(ct) and rec.component(ct).classification_state
                     in (ClassificationState.AMBIGUOUS, ClassificationState.REVIEW_REQUIRED))
                    for ct in spec.required_components
                )
                st = ClassificationState.REVIEW_REQUIRED if any_ambiguous else ValidationState.INVALID
                reason = (
                    f"required component(s) not resolved: {missing}"
                    + (" (ambiguous — assign in study_manifest.yaml)" if any_ambiguous
                       else " (not present in this system)")
                )
                remediation = (
                    "Edit study_manifest.yaml, set the ambiguity 'resolution' field "
                    "(e.g. 'receptor=A;peptide=B'), and re-run with --manifest."
                    if any_ambiguous else
                    "This system does not contain the component; the analysis is skipped."
                )
                report.observable_statuses.append(ObservableSystemStatus(
                    rec.system_id, analysis_id, st, reason, remediation))
                if any_ambiguous:
                    n_review += 1
                else:
                    n_skip += 1
                continue
            report.observable_statuses.append(ObservableSystemStatus(
                rec.system_id, analysis_id, ValidationState.VALID))
            n_ok += 1
        lines.append(
            f"  {analysis_id:<28} ready:{n_ok}  skipped:{n_skip}  review_required:{n_review}"
        )
        report.counts[analysis_id] = {"ready": n_ok, "skipped": n_skip,
                                      "review_required": n_review}

    # ── overall study state ─────────────────────────────────────────────
    if n_invalid == n:
        report.study_state = ValidationState.INVALID
    elif n_invalid or n_warn:
        report.study_state = ValidationState.WARNING
    else:
        report.study_state = ValidationState.VALID

    report.counts["systems"] = {"valid": n_valid, "warning": n_warn, "invalid": n_invalid}
    report.lines = lines
    return report
