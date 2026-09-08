"""Duration evidence for legacy plans; no tools or source writes."""
import math
import re
from pathlib import Path

from analysis.campaign.trajectory.stage_classifier import _find_mdp, _parse_mdp, mdp_timing


def duration_evidence(rec, directory):
    directory = Path(directory)
    sources, issues = [], []
    intended = {}
    for trajectory in sorted(rec.production_trajectory_paths):
        mdp = _find_mdp(Path(trajectory))
        if mdp:
            _, ps = mdp_timing(_parse_mdp(mdp))
            if ps is not None and math.isfinite(ps) and ps > 0:
                intended[str(mdp)] = ps / 1000
            else:
                issues.append(f'invalid or unbounded production MDP duration: {mdp}')
    for path, ns in sorted(intended.items()):
        sources.append(f'{path}: nsteps * dt = {ns:g} ns (intended, not measured)')
    if len(set(intended.values())) > 1:
        issues.append('conflicting production MDP durations')
    ti = rec.trajectory_inspection
    tpr_ns = ti.tpr_duration_ps / 1000 if ti and ti.tpr_duration_ps is not None else None
    if tpr_ns is not None and math.isfinite(tpr_ns) and tpr_ns > 0:
        sources.append(f'{ti.topology_path}: gmx dump nsteps * dt = {tpr_ns:g} ns (TPR intended, not measured)')
        if any(abs(tpr_ns - ns) > max(.01, ns * .01) for ns in intended.values()):
            issues.append('TPR intended duration conflicts with production MDP duration')
    else:
        tpr_ns = None
    measured = ti.total_duration_ps if ti else None
    failed = ti and any(w.code in ('gmx_check_failed', 'gmx_unavailable', 'gmx_check_incomplete') for w in ti.warnings)
    if measured is not None and math.isfinite(measured) and measured > 0 and not failed:
        value, confidence = measured / 1000, .95
        sources.insert(0, f'gmx check: production trajectory span = {value:g} ns; ' + ', '.join(sorted(rec.production_trajectory_paths)))
        if len(rec.production_trajectory_paths) == 1 and any(abs(value - ns) > max(.01, ns * .01) for ns in intended.values()):
            issues.append('measured trajectory duration conflicts with intended production MDP duration')
        if tpr_ns is not None and len(rec.production_trajectory_paths) == 1 and abs(value - tpr_ns) > max(.01, tpr_ns * .01):
            issues.append('measured trajectory duration conflicts with intended TPR duration')
    elif tpr_ns is not None:
        value, confidence = tpr_ns, .85
    elif intended:
        value, confidence = max(intended.values()), .65
    else:
        # Compatibility for records supplied without paired MDP files. Measured
        # artifacts are explicitly excluded: their end time is not a duration.
        values = [a.end_time_ps / 1000 for a in rec.trajectory_artifacts
                  if a.path in rec.production_trajectory_paths and a.n_frames is None
                  and a.end_time_ps is not None and math.isfinite(a.end_time_ps) and a.end_time_ps > 0]
        if values:
            value, confidence = max(values), .65
            sources.append('production artifact: intended duration from stage classifier (not measured)')
        else:
            value, confidence = None, 0
    hint = re.search(r'(\d+(?:\.\d+)?)ns', directory.name, re.I)
    if hint:
        sources.append(f'folder duration hint: {directory.name}')
        if value is None:
            value, confidence = float(hint[1]), .35
        elif confidence > .5 and abs(float(hint[1]) - value) > 5:
            issues.append('folder duration conflicts with production metadata')
    if not rec.production_trajectory_paths:
        issues.append('missing production trajectory')
    elif any(not Path(p).is_file() for p in rec.production_trajectory_paths):
        issues.append('production trajectory path does not exist')
    if failed:
        issues.append('production trajectory inspection failed; measured duration unavailable')
    return {'value': value, 'confidence': confidence, 'provenance': sources}, issues
