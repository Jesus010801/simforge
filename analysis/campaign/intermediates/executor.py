"""Resolve a :class:`DependencyPlan`: reuse or compute each node once.

Per node, in topological order:

1. a planning error, or a dependency that is not usable → BLOCKED (only the
   nodes and consumers that need it; everything else proceeds);
2. applicability against the current system state (RESOLVED components,
   ACTIVE annotations) → else BLOCKED;
3. the node's own coordinate view, planned by the Phase 5 policy for the
   node's purpose (via ``view_resolver``) — never inherited from a consumer;
4. definition evidence + dependency identities + input identity → identity;
5. the store: a valid entry → CACHED; otherwise compute in a build dir,
   publish atomically → COMPUTED (dry run: PLANNED, nothing written).
"""
from __future__ import annotations

import shutil
import traceback
from pathlib import Path
from typing import Callable, Optional

from analysis.campaign.intermediates import registry
from analysis.campaign.intermediates.base import (
    INTERMEDIATE_SCHEMA, IntermediateContext, IntermediateResult, IntermediateScope,
    IntermediateStatus,
)
from analysis.campaign.intermediates.planner import DependencyPlan
from analysis.campaign.intermediates.store import IntermediateConflict, IntermediateStore

#: (requirements, purpose) → TrajectoryView, planned by the policy (orchestration supplies it)
ViewResolver = Callable[..., object]


def identities(spec, evidence: dict, dep_results: list[IntermediateResult], ctx) -> dict:
    """definition / input / artifact identity (sha256 tokens)."""
    from analysis.campaign.results import definition_token
    definition = {"schema": INTERMEDIATE_SCHEMA, "intermediate": spec.id,
                  "version": spec.version, "evidence": evidence,
                  "dependencies": [[d.intermediate_id, d.version, d.artifact_identity]
                                   for d in dep_results]}
    if spec.scope == IntermediateScope.VIEW:
        v = ctx.trajectory_view
        inputs = {"scope": spec.scope, "view_ref": v.cache_key,
                  "sources": [f.digest for f in v.source_fingerprints]}
    else:
        fp = ctx.system.source_fingerprints.get("structure")
        inputs = {"scope": spec.scope, "structure": fp.digest if fp else None}
    d, i = definition_token(definition), definition_token(inputs)
    return {"definition": definition, "inputs": inputs, "definition_identity": d,
            "input_identity": i,
            "artifact_identity": definition_token({"definition": d, "inputs": i})}


def resolve_plan(plan: DependencyPlan, *, system, semantic_index, topology_path: str,
                 structure_path: Optional[str], store: IntermediateStore,
                 view_resolver: Optional[ViewResolver] = None, gmx: str = "gmx",
                 dry_run: bool = False, time_index_cache_dir=None
                 ) -> dict[str, IntermediateResult]:
    results: dict[str, IntermediateResult] = {}
    for key in plan.order:
        node = plan.nodes[key]
        req = node.request

        def blocked(reason, status=IntermediateStatus.BLOCKED, spec=None):
            results[key] = IntermediateResult(
                intermediate_id=req.id, version=spec.version if spec else (req.version or 0),
                request=req.to_dict(), status=status, reason=reason,
                dependencies=[results[d].reference() for d in node.depends_on if d in results])
        if node.error:
            blocked(node.error)
            continue
        spec = registry.get(req.id)
        bad = [results[d] for d in node.depends_on if not results[d].usable]
        if bad:
            blocked("; ".join(f"dependency {b.request.get('key')} {b.status}: {b.reason}"
                              for b in bad), spec=spec)
            continue
        appl = spec.applicability(system, req.params, semantic_index)
        if not appl.get("applicable"):
            blocked("; ".join(appl.get("reasons") or ["not applicable"]), spec=spec)
            continue
        view = None
        reqs = spec.trajectory_requirements(req.params)
        if reqs is not None:
            if view_resolver is None:
                blocked("no coordinate view resolver supplied", spec=spec)
                continue
            view = view_resolver(reqs, spec.purpose)
            if not view.safe or not view.cache_key:
                blocked("coordinate view not admitted: " + "; ".join(
                    w.message for w in view.warnings) or "unsafe view", spec=spec)
                continue
        deps = [results[d] for d in node.depends_on]
        ctx = IntermediateContext(system=system, semantic_index=semantic_index,
                                  trajectory_view=view, topology_path=topology_path,
                                  structure_path=structure_path, parameters=req.params,
                                  work_dir=Path("."), gmx=gmx, dry_run=dry_run,
                                  time_index_cache_dir=time_index_cache_dir,
                                  dependencies={d: results[d] for d in node.depends_on})
        try:
            evidence = spec.definition_evidence(ctx)
        except Exception as exc:  # noqa: BLE001 — never crash the run
            blocked(f"definition evidence failed: {exc}", IntermediateStatus.FAILED, spec)
            continue
        if evidence is None:
            blocked("definition evidence unavailable", spec=spec)
            continue
        ids = identities(spec, evidence, deps, ctx)
        refs = [d.reference() for d in deps]
        base = dict(intermediate_id=spec.id, version=spec.version, request=req.to_dict(),
                    definition_identity=ids["definition_identity"],
                    input_identity=ids["input_identity"],
                    artifact_identity=ids["artifact_identity"], definition_evidence=evidence,
                    view_ref=view.cache_key if view is not None else None, dependencies=refs)
        expect = [d.artifact_identity for d in deps]
        cached, why = store.load(ids["artifact_identity"], expect_dependencies=expect)
        if cached is not None:
            results[key] = cached
            continue
        if dry_run:
            results[key] = IntermediateResult(status=IntermediateStatus.PLANNED,
                                              reason=f"would compute ({why})", **base)
            continue
        with store.lock(ids["artifact_identity"]):
            cached, why = store.load(ids["artifact_identity"], expect_dependencies=expect)
            if cached is not None:                       # another builder published it
                results[key] = cached
                continue
            build = store.new_build_dir(ids["artifact_identity"])
            ctx.work_dir = build
            try:
                arrays, artifacts, prov = spec.compute(ctx)
                res = IntermediateResult(status=IntermediateStatus.COMPUTED, arrays=arrays,
                                         artifacts=artifacts, provenance=prov, **base)
                res.provenance["cache_miss_reason"] = why
                results[key] = store.publish(build, res)
            except IntermediateConflict:
                raise
            except Exception as exc:  # noqa: BLE001
                shutil.rmtree(build, ignore_errors=True)
                results[key] = IntermediateResult(
                    status=IntermediateStatus.FAILED,
                    reason=f"{type(exc).__name__}: {exc}",
                    provenance={"traceback": traceback.format_exc()[-2000:]}, **base)
    return results


def consumer_dependencies(plan: DependencyPlan, results: dict, consumer: str
                          ) -> tuple[dict, Optional[str], Optional[str]]:
    """``(resolved direct deps by key, blocking status, reason)`` for one consumer."""
    direct = {k: results[k] for k in plan.consumers.get(consumer, [])}
    bad = [r for r in direct.values() if not r.usable and r.status != IntermediateStatus.PLANNED]
    if not bad:
        return direct, None, None
    status = (IntermediateStatus.FAILED if any(b.status == IntermediateStatus.FAILED for b in bad)
              else IntermediateStatus.BLOCKED)
    return direct, status, "; ".join(f"dependency {b.request.get('key')} {b.status}: {b.reason}"
                                     for b in bad)
