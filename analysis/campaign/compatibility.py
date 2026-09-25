"""Scientific compatibility of existing results (Trajectory Review Phase 8).

Decides whether an already stored result may satisfy a new request::

    requested observable -> definition evidence -> candidate result
        -> per-dimension checks -> COMPATIBLE | INCOMPATIBLE |
           INSUFFICIENT_EVIDENCE | EXPLICIT_IMPORT_ONLY

What equality meant before this phase: a native ``definition_token`` hashes the
whole definition evidence — observable, selection atom-set / annotation
evidence, scientific parameters, coordinate semantics, .tpr / reference digest,
backend + version — *and* ``trajectory_view.cache_key`` (which encodes the
source trajectory digests and preprocessing).  It therefore conflates "same
definition" with "same input".  Here the two are separated additively:

* ``definition_identity`` — token of the evidence without ``trajectory_view``;
* input identity — the view reference plus the analysed file's content digest.

Nothing is inferred from filenames, titles or directory names.  Stored output
files are re-fingerprinted before any reuse.  Imported series are never
auto-reused.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

from analysis.campaign.models import AnalysisResult

RESOLVER_VERSION = "simforge/compatibility/v1"


class CompatibilityStatus:
    COMPATIBLE = "compatible"
    INCOMPATIBLE = "incompatible"
    INSUFFICIENT_EVIDENCE = "insufficient_evidence"
    EXPLICIT_IMPORT_ONLY = "explicit_import_only"


class CheckOutcome:
    MATCH = "match"
    MISMATCH = "mismatch"
    UNKNOWN = "unknown"          # evidence missing -> never a basis for reuse


class ProvenanceTier:
    NATIVE = "tier1_native"                 # full SimForge provenance
    DECLARED = "tier2_declared"             # external / header evidence only
    NAME_ONLY = "tier3_name_only"           # a file that merely looks right


@dataclass
class CompatibilityCheck:
    dimension: str
    outcome: str
    detail: str = ""
    evidence: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {"dimension": self.dimension, "outcome": self.outcome,
                "detail": self.detail, "evidence": self.evidence}


@dataclass
class Candidate:
    """One stored result that might satisfy a request."""
    source: str                               # provenance path / file path
    tier: str
    result: Optional[AnalysisResult] = None
    system_id: Optional[str] = None
    in_current_location: bool = False

    def to_dict(self) -> dict:
        r = self.result
        return {"source": self.source, "tier": self.tier, "system_id": self.system_id,
                "analysis_id": r.analysis_id if r else None,
                "definition_token": r.definition_token if r else None,
                "in_current_location": self.in_current_location}


@dataclass
class CompatibilityResult:
    status: str
    candidate: Candidate
    requested: dict = field(default_factory=dict)
    checks: list[CompatibilityCheck] = field(default_factory=list)
    reason: str = ""
    resolver: str = RESOLVER_VERSION

    def to_dict(self) -> dict:
        return {"resolver": self.resolver, "status": self.status, "reason": self.reason,
                "candidate": self.candidate.to_dict(), "requested": self.requested,
                "checks": [c.to_dict() for c in self.checks]}

    def summary(self) -> str:
        return "\n".join(f"{c.dimension:18s} {c.outcome.upper()}"
                         + (f"  ({c.detail})" if c.detail and c.outcome != CheckOutcome.MATCH else "")
                         for c in self.checks) + f"\n-> {self.status.upper()}: {self.reason}"


# ═══════════════════════════════════════════════════════════════════════════════
# Identities
# ═══════════════════════════════════════════════════════════════════════════════

INPUT_KEYS = ("trajectory_view",)


def definition_identity(evidence: Optional[dict]) -> Optional[str]:
    """Scientific-definition identity: the evidence without input identity."""
    if not evidence:
        return None
    from analysis.campaign.results import definition_token
    return definition_token({k: v for k, v in evidence.items() if k not in INPUT_KEYS})


def _flatten(d: Any, prefix: str = "") -> dict[str, Any]:
    if isinstance(d, dict):
        out: dict[str, Any] = {}
        for k, v in d.items():
            out.update(_flatten(v, f"{prefix}.{k}" if prefix else str(k)))
        return out
    return {prefix: d}


def evidence_differences(a: Optional[dict], b: Optional[dict], *, limit: int = 8) -> list[str]:
    """Human-readable paths where two definition evidences differ."""
    fa, fb = _flatten({k: v for k, v in (a or {}).items() if k not in INPUT_KEYS}), \
        _flatten({k: v for k, v in (b or {}).items() if k not in INPUT_KEYS})
    diffs = [k for k in sorted(set(fa) | set(fb)) if fa.get(k) != fb.get(k)]
    return diffs[:limit]


# ═══════════════════════════════════════════════════════════════════════════════
# Output-file evidence
# ═══════════════════════════════════════════════════════════════════════════════

_NUM = re.compile(r"^[-+]?\d*\.?(\d*)(?:[eE]([-+]?\d+))?$")


def output_resolution(path: str | Path, column: int) -> Optional[dict]:
    """Printed numerical resolution of one column (from the file itself)."""
    decimals: set[int] = set()
    sig: set[int] = set()
    try:
        for line in Path(path).read_text(errors="replace").splitlines():
            s = line.strip()
            if not s or s[0] in "#@":
                continue
            tok = s.split()[column] if len(s.split()) > column else None
            m = _NUM.match(tok or "")
            if not m:
                continue
            if m.group(2) is not None:
                mant = tok.lower().split("e")[0].lstrip("+-").replace(".", "").lstrip("0")
                sig.add(len(mant) or 1)
            else:
                decimals.add(len(m.group(1)))
    except OSError:
        return None
    if sig and not decimals:
        return {"format": "exponent", "significant_digits": min(sig)}
    if decimals and not sig:
        n = min(decimals)
        return {"format": "fixed", "decimals": n, "resolution": 10.0 ** (-n)}
    return None


def storage_integrity(result: AnalysisResult) -> tuple[bool, str]:
    """Every referenced output file exists and still has its recorded fingerprint."""
    from analysis.campaign.fingerprint import fingerprint_file
    if not result.arrays:
        return False, "result has no ResultArray"
    for arr in result.arrays:
        st = arr.storage
        if st is None:
            return False, f"array {arr.name!r} has no storage"
        p = Path(st.path)
        if not p.is_file():
            return False, f"output file missing: {p}"
        if st.format not in ("xvg", "csv", "npy", "npz") or p.suffix.lstrip(".").lower() != st.format:
            return False, f"storage format {st.format!r} does not match {p.name}"
        if st.fingerprint is None:
            return False, f"no recorded fingerprint for {p.name}"
        now = fingerprint_file(p, strong=st.fingerprint.mode == "strong")
        if now.digest != st.fingerprint.digest:
            return False, f"output file changed after the analysis: {p.name}"
    return True, "all output files present with recorded fingerprints"


# ═══════════════════════════════════════════════════════════════════════════════
# Discovery (locations SimForge wrote; no filesystem crawl)
# ═══════════════════════════════════════════════════════════════════════════════

def discover_candidates(output_root: str | Path, system_id: str, analysis_id: str,
                        *, current_dir: Optional[Path] = None) -> list[Candidate]:
    """Stored results for (system, analysis) under a campaign output root.

    Native candidates come from ``provenance.json``; an output-looking file
    without provenance in the same directory is reported as a name-only
    candidate (never reusable).
    """
    root = Path(output_root)
    obs = root / "systems" / system_id / "observables" / analysis_id
    out: list[Candidate] = []
    if not obs.is_dir():
        return out
    prov = obs / "provenance.json"
    current = current_dir.resolve() if current_dir else None
    if prov.is_file():
        try:
            d = json.loads(prov.read_text())
            res = AnalysisResult.from_dict(d["result"])
            out.append(Candidate(source=str(prov.resolve()), tier=ProvenanceTier.NATIVE,
                                 result=res, system_id=d.get("system_id"),
                                 in_current_location=current == obs.resolve()))
        except (OSError, ValueError, KeyError):
            out.append(Candidate(source=str(prov.resolve()), tier=ProvenanceTier.NAME_ONLY))
    else:
        for f in sorted(obs.glob("*.xvg")):
            out.append(Candidate(source=str(f.resolve()), tier=ProvenanceTier.NAME_ONLY))
    return out


# ═══════════════════════════════════════════════════════════════════════════════
# Evaluation
# ═══════════════════════════════════════════════════════════════════════════════

def evaluate(candidate: Candidate, *, requested_evidence: Optional[dict],
             requested_view_ref: Optional[str], requested_schema: list,
             requested_timeline: Optional[list[float]] = None,
             view_invariant: bool = False) -> CompatibilityResult:
    """Check one candidate against the request, dimension by dimension."""
    req_id = definition_identity(requested_evidence)
    requested = {"definition_identity": req_id, "view_ref": requested_view_ref,
                 "timeline_frames": len(requested_timeline) if requested_timeline else None}
    res = CompatibilityResult(status=CompatibilityStatus.INSUFFICIENT_EVIDENCE,
                              candidate=candidate, requested=requested)
    add = res.checks.append

    if candidate.tier == ProvenanceTier.NAME_ONLY or candidate.result is None:
        add(CompatibilityCheck("provenance", CheckOutcome.UNKNOWN,
                               "only a file name / unreadable record — no scientific evidence"))
        res.reason = "no provenance: a matching file name never establishes equivalence"
        return res
    r = candidate.result
    arrays = r.arrays or []
    if any(a.attrs.get("externally_supplied") for a in arrays) or not r.definition_evidence:
        add(CompatibilityCheck("provenance", CheckOutcome.UNKNOWN,
                               "externally supplied / no native definition evidence"))
        res.status = CompatibilityStatus.EXPLICIT_IMPORT_ONLY
        res.reason = ("externally supplied result: usable only by explicit user choice; "
                      "scientific equivalence is not verified")
        return res
    if r.status not in ("success",):
        add(CompatibilityCheck("provenance", CheckOutcome.MISMATCH,
                               f"candidate status is {r.status!r}, not a completed computation"))
        res.status, res.reason = CompatibilityStatus.INCOMPATIBLE, "candidate is not a completed result"
        return res
    add(CompatibilityCheck("provenance", CheckOutcome.MATCH, ProvenanceTier.NATIVE))

    # A. scientific definition (input-independent)
    cand_id = definition_identity(r.definition_evidence)
    if req_id is None:
        add(CompatibilityCheck("definition", CheckOutcome.UNKNOWN,
                               "the request has no definition evidence (not applicable now?)"))
    elif cand_id == req_id:
        add(CompatibilityCheck("definition", CheckOutcome.MATCH, "",
                               {"definition_identity": cand_id}))
    else:
        add(CompatibilityCheck("definition", CheckOutcome.MISMATCH,
                               "differs at " + ", ".join(evidence_differences(
                                   requested_evidence, r.definition_evidence)),
                               {"requested": req_id, "candidate": cand_id}))
    # backend / version (part of the definition, reported on its own)
    rb, cb = (requested_evidence or {}).get("backend"), r.definition_evidence.get("backend")
    add(CompatibilityCheck("backend", CheckOutcome.MATCH if rb == cb else CheckOutcome.MISMATCH,
                           "" if rb == cb else f"requested {rb} vs candidate {cb}"))

    # C. coordinate view
    cand_view = (r.definition_evidence.get("trajectory_view") or {}).get("cache_key")
    if cand_view and cand_view == requested_view_ref:
        add(CompatibilityCheck("view", CheckOutcome.MATCH, "", {"view_ref": cand_view}))
    else:
        add(CompatibilityCheck(
            "view", CheckOutcome.MISMATCH,
            f"candidate view {cand_view} vs requested {requested_view_ref}"
            + (" (observable declares view invariance — not used without equal sources)"
               if view_invariant else "")))

    # B. source data (the analysed file's content)
    digests = {a.axis("time").alignment.trajectory_digest for a in arrays
               if a.axis("time") and a.axis("time").alignment}
    add(CompatibilityCheck(
        "source", CheckOutcome.MATCH if cand_view == requested_view_ref and digests
        else CheckOutcome.UNKNOWN if not digests else CheckOutcome.MISMATCH,
        "view identity encodes the source trajectory fingerprints"
        if cand_view == requested_view_ref else "view differs",
        {"analysed_file_digests": sorted(d for d in digests if d)}))

    # D. timeline
    add(_timeline_check(r, requested_timeline))

    # E. output schema / integrity / resolution
    want = {(s.quantity, s.unit, s.axis_signature()) for s in requested_schema}
    have = {(a.quantity, a.unit, a.axis_signature()) for a in arrays}
    add(CompatibilityCheck("output_schema", CheckOutcome.MATCH if want == have else CheckOutcome.MISMATCH,
                           "" if want == have else f"requested {sorted(want)} vs {sorted(have)}"))
    ok, why = storage_integrity(r)
    add(CompatibilityCheck("output_integrity", CheckOutcome.MATCH if ok else CheckOutcome.MISMATCH, why))
    reso = {a.name: output_resolution(a.storage.path, a.storage.column or 0)
            for a in arrays if a.storage} if ok else {}
    add(CompatibilityCheck(
        "output_resolution",
        CheckOutcome.MATCH if ok and all(reso.values()) else CheckOutcome.UNKNOWN,
        "same backend, version and command options -> same printed resolution"
        if ok and all(reso.values()) else "printed resolution could not be established",
        {"measured": reso}))

    outcomes = {c.dimension: c.outcome for c in res.checks}
    if any(o == CheckOutcome.MISMATCH for o in outcomes.values()):
        res.status = CompatibilityStatus.INCOMPATIBLE
        bad = [c for c in res.checks if c.outcome == CheckOutcome.MISMATCH]
        res.reason = "; ".join(f"{c.dimension}: {c.detail}" for c in bad)
    elif any(o == CheckOutcome.UNKNOWN for o in outcomes.values()):
        res.status = CompatibilityStatus.INSUFFICIENT_EVIDENCE
        res.reason = "; ".join(f"{c.dimension}: {c.detail}" for c in res.checks
                               if c.outcome == CheckOutcome.UNKNOWN)
    else:
        res.status = CompatibilityStatus.COMPATIBLE
        res.reason = "same definition, source, view, timeline and output semantics"
    return res


def same_time(a: float, b: float) -> bool:
    """Two recorded timestamps denote the same frame time (printed-precision
    tolerance: 1e-3 ps or float32 relative precision, whichever is larger)."""
    return abs(a - b) <= max(1e-3, 6e-6 * abs(b))


def _timeline_check(r: AnalysisResult, requested: Optional[list[float]]) -> CompatibilityCheck:
    """Exact sample match only (no slicing, no interpolation, order preserved)."""
    from analysis.campaign.results import read_column
    for arr in r.arrays:
        t = arr.axis("time")
        if t is None or t.values_ref is None:
            return CompatibilityCheck("timeline", CheckOutcome.UNKNOWN, "no stored time axis")
        if t.alignment is None or t.alignment.mode != "one_row_per_frame":
            return CompatibilityCheck("timeline", CheckOutcome.UNKNOWN,
                                      "candidate rows were not aligned one-per-frame")
        if requested is None:
            return CompatibilityCheck("timeline", CheckOutcome.UNKNOWN,
                                      "requested timeline not available without a scan")
        try:
            times = read_column(t.values_ref)
        except (OSError, ValueError) as exc:
            return CompatibilityCheck("timeline", CheckOutcome.UNKNOWN, f"time axis unreadable: {exc}")
        if len(times) != len(requested):
            return CompatibilityCheck("timeline", CheckOutcome.MISMATCH,
                                      f"{len(times)} samples vs {len(requested)} requested frames")
        bad = next((i for i, (a, b) in enumerate(zip(times, requested))
                    if not same_time(a, b)), None)
        if bad is not None:
            return CompatibilityCheck("timeline", CheckOutcome.MISMATCH,
                                      f"sample {bad}: {times[bad]} ps vs requested {requested[bad]} ps")
    return CompatibilityCheck("timeline", CheckOutcome.MATCH,
                              f"{len(requested)} samples match the requested frames in order",
                              {"n": len(requested), "first_ps": requested[0] if requested else None,
                               "last_ps": requested[-1] if requested else None})


def resolve_compatible_result(candidates: list[Candidate], **request
                              ) -> tuple[Optional[CompatibilityResult], list[CompatibilityResult]]:
    """(chosen compatible result or None, every evaluation).

    Several fully compatible candidates are chosen deterministically — the
    one in the current output location first, then by source path — never by
    a scientific "better" judgement.
    """
    reports = [evaluate(c, **request) for c in candidates]
    ok = [r for r in reports if r.status == CompatibilityStatus.COMPATIBLE]
    ok.sort(key=lambda r: (not r.candidate.in_current_location, r.candidate.source))
    return (ok[0] if ok else None), reports
