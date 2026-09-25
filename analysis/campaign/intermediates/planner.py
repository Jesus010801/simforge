"""Dependency DAG planning — pure: no data access, nothing executed.

consumers (observables, in request order) → their ``IntermediateRequest``s →
recursively each intermediate's own requests.  Identical requests (same id +
parameters) are one node however many consumers ask.  Unknown intermediates,
invalid parameters and broken version/schema contracts become *node errors*
(their consumers are blocked; unrelated consumers are not); a cycle makes
the whole graph invalid (:class:`DependencyCycleError`, with the path).

Order is deterministic: depth-first post-order following consumers and their
declared requests in order — every node after all of its dependencies.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from analysis.campaign.intermediates import registry
from analysis.campaign.intermediates.base import IntermediateRequest


class DependencyCycleError(ValueError):
    def __init__(self, path: list[str]):
        self.path = path
        super().__init__("dependency cycle: " + " -> ".join(path))


@dataclass
class DependencyNode:
    request: IntermediateRequest
    depends_on: list[str] = field(default_factory=list)       # node keys
    consumers: list[str] = field(default_factory=list)        # who asked (provenance only)
    error: str = ""                                            # planning error → blocked

    @property
    def key(self) -> str:
        return self.request.key


@dataclass
class DependencyPlan:
    nodes: dict[str, DependencyNode] = field(default_factory=dict)
    order: list[str] = field(default_factory=list)            # topological
    consumers: dict[str, list[str]] = field(default_factory=dict)   # consumer → direct keys

    def closure(self, key: str) -> list[str]:
        """Every node ``key`` transitively needs, itself included, in order."""
        need, stack = set(), [key]
        while stack:
            k = stack.pop()
            if k not in need:
                need.add(k)
                stack.extend(self.nodes[k].depends_on)
        return [k for k in self.order if k in need]

    def explain(self, statuses: Optional[dict[str, str]] = None) -> str:
        """Human-readable tree: which consumer needs what, and sharing."""
        statuses = statuses or {}
        lines, seen = [], {}
        for consumer, keys in self.consumers.items():
            lines.append(consumer)
            for k in keys:
                self._explain(k, 1, lines, seen, consumer, statuses)
        return "\n".join(lines)

    def _explain(self, key, depth, lines, seen, consumer, statuses):
        node = self.nodes[key]
        state = statuses.get(key, "planned")
        note = f"  (shared with {seen[key]})" if key in seen and seen[key] != consumer else ""
        lines.append("  " * depth + f"└── {key}  [{node.error or state}]{note}")
        seen.setdefault(key, consumer)
        for d in node.depends_on:
            self._explain(d, depth + 1, lines, seen, consumer, statuses)

    def to_dict(self) -> dict:
        return {"order": self.order, "consumers": self.consumers,
                "nodes": {k: {"request": n.request.to_dict(), "depends_on": n.depends_on,
                              "consumers": n.consumers, "error": n.error}
                          for k, n in self.nodes.items()}}


def _contract_error(req: IntermediateRequest) -> str:
    if not registry.is_registered(req.id):
        return f"unknown intermediate {req.id!r}"
    spec = registry.get(req.id)
    if req.version is not None and req.version != spec.version:
        return (f"contract mismatch: consumer requires {req.id} v{req.version}, "
                f"registered v{spec.version}")
    problems = spec.validate(req.params)
    return "; ".join(problems)


def plan_dependencies(consumers: list[tuple[str, list[IntermediateRequest]]]) -> DependencyPlan:
    plan = DependencyPlan()
    done: set[str] = set()

    def visit(req: IntermediateRequest, stack: list[str], consumer: str) -> str:
        key = req.key
        if key in stack:
            raise DependencyCycleError(stack[stack.index(key):] + [key])
        node = plan.nodes.get(key)
        if node is None:
            node = plan.nodes[key] = DependencyNode(req)
            node.error = _contract_error(req)
        elif req.version is not None and node.request.version not in (None, req.version):
            node.error = node.error or (f"consumers pin different versions of {req.id} "
                                        f"({node.request.version} vs {req.version})")
        if consumer not in node.consumers:
            node.consumers.append(consumer)
        if key in done:
            return key
        if not node.error:
            for dep in registry.get(req.id).dependencies(req.params):
                dk = visit(dep, stack + [key], consumer)
                if dk not in node.depends_on:
                    node.depends_on.append(dk)
        done.add(key)
        plan.order.append(key)
        return key

    for consumer, requests in consumers:
        keys = []
        for req in requests:
            k = visit(req, [], consumer)
            if k not in keys:
                keys.append(k)
        plan.consumers[consumer] = keys
    return plan
