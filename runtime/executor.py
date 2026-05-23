from __future__ import annotations

import asyncio
import time
from dataclasses import dataclass, field
from pathlib import Path

from rich.console import Console

from executors.base_executor import BaseExecutor
from executors.execution_state import StepExecutionRecord, StepStatus

from runtime.artifacts import ArtifactRegistry, SemanticRole
from runtime.cache import StepCacheManager
from runtime.events import BoundPublisher, EventBus, EventSeverity, EventType
from runtime.execution_backend import LocalSubprocessBackend
from runtime.journal import JournalWriter
from runtime.metrics import SystemMetricsCollector
from runtime.stream import AsyncProcessRunner

MAX_REMEDIATION_DEPTH = 3

_console = Console(highlight=False)


def _json_dumps(obj: object) -> str:
    """JSON serializer that handles Path objects."""
    import json
    def _default(o):
        if isinstance(o, Path):
            return str(o)
        raise TypeError(f"Object of type {type(o)} is not JSON serializable")
    return json.dumps(obj, indent=2, default=_default)

# file-extension → semantic role
_ROLE_MAP: dict[str, SemanticRole] = {
    ".top":  "topology",
    ".itp":  "topology",
    ".gro":  "coordinates",
    ".xtc":  "trajectory",
    ".trr":  "trajectory",
    ".cpt":  "checkpoint",
    ".xvg":  "analysis",
    ".edr":  "analysis",
    ".mdp":  "parameter",
    ".log":  "log",
}

_SCRIPT_CANDIDATES: list[str] = [
    "run.sh",
    "run_md.sh",
    "run_nvt.sh",
    "run_npt.sh",
    "run_analysis.sh",
    "commands.sh",
]


def _semantic_role(path: Path) -> SemanticRole:
    return _ROLE_MAP.get(path.suffix.lower(), "other")


# ── Resume planning model ─────────────────────────────────────────────────────

@dataclass
class ResumePlan:
    """Classification of all steps before execution starts."""
    reusable:      list[str]       # step_ids with valid cache + intact artifacts
    resume_point:  str | None      # first step that must actually run
    pending:       list[str]       # steps after the resume point
    invalidated:   list[str] = field(default_factory=list)  # stale cache entries cleared
    recoverable:   list[str] = field(default_factory=list)  # steps with .cpt but no final .gro


# ═══════════════════════════════════════════════════════════════════════════════

class RuntimeExecutor(BaseExecutor):
    """
    Observable, event-driven executor.

    Extends BaseExecutor with:
      - Resume planning phase (_pre_run_hook): classifies all steps as
        REUSABLE (cached + artifacts intact) or PENDING before execution starts.
        Reusable steps are silently skipped via _should_skip(); the user sees
        a clear resume summary instead of per-step "cache hit" noise.
      - asyncio subprocess streaming (AsyncProcessRunner)
      - structured ExecutionEvent emission via EventBus
      - append-only JSONL journal
      - content-addressed step caching
      - file provenance via ArtifactRegistry
      - live system metrics (CPU/RAM/GPU)
    """

    def __init__(
        self,
        workspace_path: str | Path,
        dry_run:        bool     = True,
        bus:            EventBus | None = None,
        metrics_interval_s: int = 15,
    ) -> None:
        super().__init__(workspace_path, dry_run)
        wp = Path(workspace_path)

        self._bus      = bus or EventBus()
        self._journal  = JournalWriter(wp)
        self._journal.register(self._bus)

        self._artifacts = ArtifactRegistry(wp)
        self._cache     = StepCacheManager(wp)

        self._workspace_id       = wp.name
        self._metrics_interval_s = metrics_interval_s
        self._cached_steps: set[str] = set()      # populated by _pre_run_hook
        self._resumable_steps: dict[str, str] = {}  # step_id → resume script name
        self._backend = LocalSubprocessBackend()

    # ── EventBus access (for external subscribers) ────────────────────────────

    @property
    def bus(self) -> EventBus:
        return self._bus

    # ── Resume planning (BaseExecutor hooks) ──────────────────────────────────

    def _pre_run_hook(self) -> None:
        """
        Runs once after state is initialized, before the step loop.

        Classifies every step as reusable (cached + outputs present) or pending.
        Populates _cached_steps so _should_skip() can fast-path them.
        Prints a resume summary when at least one step can be reused.
        """
        if self.dry_run or self.state is None:
            return
        plan = self._plan_resume()
        if plan.reusable:
            self._print_resume_plan(plan)

    def _should_skip(self, record: StepExecutionRecord) -> bool:
        """True when the step's cache is valid and all artifacts are on disk."""
        if record.step_id not in self._cached_steps:
            return False
        # Emit a structured event so subscribers (journal, tests) can observe it.
        pub = self._bus.make_publisher(self._workspace_id, record.step_id)
        pub.emit(
            EventType.CACHE_HIT,
            message=f"Step '{record.step_id}' — reusing cached artifacts",
        )
        return True

    def _plan_resume(self) -> ResumePlan:
        """
        Inspect every step's fingerprint and expected outputs.

        Side effects:
          - Adds valid steps to self._cached_steps.
          - Invalidates stale cache entries (fingerprint match but outputs gone).
          - Populates self._resumable_steps for steps with checkpoints.
        """
        reusable:     list[str] = []
        pending:      list[str] = []
        invalidated:  list[str] = []
        recoverable:  list[str] = []
        resume_point: str | None = None

        for record in self.state.steps:
            step_dir         = Path(record.step_dir)
            meta             = self._read_metadata(step_dir)
            expected_outputs = meta.get("expected_outputs", [])
            input_paths      = [
                (step_dir / rel).resolve()
                for rel in meta.get("required_inputs", [])
            ]
            fp = self._cache.fingerprint(meta.get("params", {}), input_paths)

            if self._cache.is_cached(record.step_id, fp):
                missing = [
                    name for name in expected_outputs
                    if not (step_dir / name).exists()
                ]
                if not missing:
                    reusable.append(record.step_id)
                    self._cached_steps.add(record.step_id)
                    continue
                else:
                    # Fingerprint matches but artifacts are gone — stale entry.
                    self._cache.invalidate(record.step_id)
                    invalidated.append(record.step_id)
                    self._log(
                        f"[CACHE INVALIDATED] {record.step_id} — "
                        f"missing outputs: {missing}"
                    )

            # ── Checkpoint detection ──────────────────────────────────────────
            # Check for interrupted MD steps: .cpt present, final .gro missing.
            # Priority: npt > nvt > md (prefer the most-advanced checkpoint).
            cpt_info = self._detect_checkpoint(step_dir)
            if cpt_info is not None:
                _cpt_path, resume_script = cpt_info
                recoverable.append(record.step_id)
                self._resumable_steps[record.step_id] = resume_script
                record.status = StepStatus.RECOVERABLE
                self._log(
                    f"[CHECKPOINT] {record.step_id} — "
                    f"will resume with {resume_script}"
                )

            if resume_point is None:
                resume_point = record.step_id
            else:
                pending.append(record.step_id)

        return ResumePlan(
            reusable=reusable,
            resume_point=resume_point,
            pending=pending,
            invalidated=invalidated,
            recoverable=recoverable,
        )

    @staticmethod
    def _detect_checkpoint(step_dir: Path) -> tuple[Path, str] | None:
        """
        Return (cpt_path, resume_script_name) if the step has a checkpoint but
        no corresponding final .gro (i.e., the run was interrupted mid-flight).

        Checks in priority order: npt → nvt → md, so that a partially-done
        equilibration resumes from the furthest-advanced phase.
        Returns None if no recoverable checkpoint is found.
        """
        for deffnm, resume_script in [
            ("npt", "run_npt_resume.sh"),
            ("nvt", "run_nvt_resume.sh"),
            ("md",  "run_md_resume.sh"),
        ]:
            cpt = step_dir / f"{deffnm}.cpt"
            gro = step_dir / f"{deffnm}.gro"
            if cpt.exists() and not gro.exists():
                return cpt, resume_script
        return None

    def _print_resume_plan(self, plan: ResumePlan) -> None:
        _console.print()
        _console.print("  [dim]Analyzing existing workflow state...[/dim]")
        _console.print()

        _console.print("  [bold]Reusable artifacts:[/bold]")
        for step_id in plan.reusable:
            _console.print(f"    [green]✓[/green]  [dim]{step_id}[/dim]")

        if plan.recoverable:
            _console.print()
            _console.print("  [bold blue]Recoverable (checkpoint detected):[/bold blue]")
            for step_id in plan.recoverable:
                script = self._resumable_steps.get(step_id, "?")
                _console.print(
                    f"    [blue]⟳[/blue]  {step_id}  "
                    f"[dim](→ {script})[/dim]"
                )

        if plan.invalidated:
            _console.print()
            _console.print("  [bold yellow]Invalidated (outputs missing):[/bold yellow]")
            for step_id in plan.invalidated:
                _console.print(f"    [yellow]![/yellow]  {step_id}")

        if plan.resume_point:
            _console.print()
            _console.print("  [bold]Resume point:[/bold]")
            _console.print(f"    [cyan]↻[/cyan]  [bold]{plan.resume_point}[/bold]")

            if plan.pending:
                _console.print()
                _console.print("  [bold]Pending downstream:[/bold]")
                for step_id in plan.pending:
                    _console.print(f"    [dim]○  {step_id}[/dim]")

            _console.print()
            _console.print(
                f"  [dim]Resuming execution from[/dim] "
                f"[bold]{plan.resume_point}[/bold]...\n"
            )
        else:
            _console.print()
            _console.print(
                "  [bold green]✓[/bold green]  "
                "All steps have valid cached artifacts.\n"
            )

    # ── BaseExecutor hooks ────────────────────────────────────────────────────

    def _run_step(self, record: StepExecutionRecord) -> None:
        """Bridge: sync BaseExecutor → async implementation."""
        asyncio.run(self._run_step_async(record))

    def _post_run_hook(self) -> None:
        """Generate scientific summary after the step loop completes."""
        if self.dry_run:
            return
        try:
            from runtime.scientific_summary import generate_summary
            summary = generate_summary(self.workspace_path)
            meta_dir = self.workspace_path / "metadata"
            meta_dir.mkdir(parents=True, exist_ok=True)
            (meta_dir / "scientific_summary.json").write_text(
                _json_dumps(summary.as_dict())
            )
            (meta_dir / "scientific_summary.md").write_text(summary.as_markdown())
            self._log(
                f"[SUMMARY] Scientific summary written "
                f"(converged={summary.converged}, analyses={len(summary.analyses)})"
            )
        except Exception as exc:
            self._log(f"[SUMMARY] Could not generate scientific summary: {exc}")

    # ── Async implementation ──────────────────────────────────────────────────

    async def _run_step_async(self, record: StepExecutionRecord) -> None:
        step_dir = Path(record.step_dir)
        meta     = self._read_metadata(step_dir)
        script   = self._find_script(step_dir, record.step_id)
        pub      = self._bus.make_publisher(self._workspace_id, record.step_id)

        # ── Manual / external steps → skip immediately ────────────────────────
        if script is None:
            step_type = meta.get("step_type", "automatic")
            if step_type in ("manual", "external", "validation"):
                record.status        = StepStatus.SKIPPED
                record.error_message = f"Step tipo '{step_type}' — requiere acción manual"
                pub.emit(EventType.STEP_SKIPPED, message=record.error_message, step_type=step_type)
            else:
                record.status        = StepStatus.SKIPPED
                record.error_message = "No se encontró script ejecutable"
                pub.emit(EventType.STEP_SKIPPED, message=record.error_message)
            return

        expected_outputs: list[str] = meta.get("expected_outputs", [])

        # ── Dry run ───────────────────────────────────────────────────────────
        if self.dry_run:
            self._dry_run_step(record, step_dir, script, expected_outputs, meta, pub)
            return

        # ── Cache check (safety net for steps not caught by _should_skip) ────
        # _pre_run_hook populates _cached_steps; _should_skip intercepts most
        # cache hits before _run_step is called. This block handles edge cases
        # (e.g., cache became valid between planning and execution).
        input_paths = [
            (step_dir / rel).resolve()
            for rel in meta.get("required_inputs", [])
        ]
        fp = self._cache.fingerprint(meta.get("params", {}), input_paths)

        if self._cache.is_cached(record.step_id, fp):
            missing_artifacts = [
                name for name in expected_outputs
                if not (step_dir / name).exists()
            ]
            if missing_artifacts:
                self._cache.invalidate(record.step_id)
                self._log(
                    f"[CACHE INVALIDATED] {record.step_id} — "
                    f"missing outputs: {missing_artifacts}"
                )
                pub.emit(
                    EventType.CACHE_MISS,
                    fingerprint=fp,
                    message=f"Cache invalidated: missing outputs {missing_artifacts}",
                    missing_artifacts=missing_artifacts,
                )
                _console.print(
                    f"  [yellow]⚠[/yellow]  {record.step_id} "
                    f"[dim](cache invalidated — outputs missing, re-running)[/dim]"
                )
            else:
                record.status = StepStatus.DONE
                pub.emit(
                    EventType.CACHE_HIT,
                    message=f"Step '{record.step_id}' unchanged — skipping",
                    fingerprint=fp,
                )
                _console.print(
                    f"  [dim]⚡ {record.step_id} (cache hit — skipped)[/dim]"
                )
                return

        pub.emit(EventType.CACHE_MISS, fingerprint=fp)

        # ── Real execution ────────────────────────────────────────────────────
        pub.emit(
            EventType.STEP_STARTED,
            message=f"Executing {record.step_id}",
            script=str(script.relative_to(step_dir)),
            expected_outputs=expected_outputs,
            dry_run=False,
        )

        runner  = AsyncProcessRunner(pub)
        metrics = SystemMetricsCollector(
            pub,
            workspace_dir=str(self.workspace_path),
            interval_s=self._metrics_interval_s,
        )

        metrics_task = asyncio.create_task(metrics.run())
        try:
            result = await runner.run(
                cmd=["bash", str(script)],
                cwd=str(step_dir),
            )
        except Exception as exc:
            record.status        = StepStatus.FAILED
            record.error_message = str(exc)
            pub.emit(
                EventType.STEP_FAILED,
                message=str(exc),
                severity=EventSeverity.ERROR,
            )
            return
        finally:
            metrics_task.cancel()
            try:
                await metrics_task
            except asyncio.CancelledError:
                pass

        if result.ns_per_day:
            metrics.record_performance(result.ns_per_day, result.hours_per_ns or 0.0)

        record.stdout    = "\n".join(result.stdout_lines[-500:])
        record.stderr    = "\n".join(result.stderr_lines[-200:])
        record.exit_code = result.returncode

        if result.returncode != 0:
            record.status        = StepStatus.FAILED
            record.error_message = (
                f"Exit code {result.returncode}. "
                f"Last stderr: {record.stderr[-400:]}"
            )
            pub.emit(
                EventType.STEP_FAILED,
                message=record.error_message,
                severity=EventSeverity.ERROR,
                returncode=result.returncode,
                wall_time_s=result.wall_time_s,
            )
            # Adaptive remediation loop (Sprint 1c)
            if record.n_remediations() < MAX_REMEDIATION_DEPTH and not self.dry_run:
                remediated = await self._attempt_remediation(record, step_dir, pub)
                if remediated:
                    return  # retry succeeded — record is now DONE
            return

        # ── Verify expected outputs ───────────────────────────────────────────
        found, missing = self._check_expected_outputs(step_dir, expected_outputs)
        record.outputs_found   = found
        record.outputs_missing = missing

        if missing:
            record.status        = StepStatus.FAILED
            record.error_message = f"Outputs esperados no encontrados: {missing}"
            pub.emit(
                EventType.STEP_FAILED,
                message=record.error_message,
                severity=EventSeverity.ERROR,
                missing_outputs=missing,
            )
            # Adaptive remediation loop (Sprint 1c)
            if record.n_remediations() < MAX_REMEDIATION_DEPTH and not self.dry_run:
                remediated = await self._attempt_remediation(record, step_dir, pub)
                if remediated:
                    return  # retry succeeded — record is now DONE
            return

        # ── Register produced artifacts ───────────────────────────────────────
        for name in found:
            p = step_dir / name
            if p.exists():
                self._artifacts.register(p, _semantic_role(p), record.step_id)
                pub.emit(
                    EventType.ARTIFACT_CREATED,
                    message=name,
                    path=str(p),
                    semantic_role=_semantic_role(p),
                )

        # ── Cache the fingerprint ─────────────────────────────────────────────
        self._cache.record(record.step_id, fp)

        record.status = StepStatus.DONE
        pub.emit(
            EventType.STEP_COMPLETED,
            message=f"Step '{record.step_id}' completed",
            wall_time_s=result.wall_time_s,
            ns_per_day=result.ns_per_day,
            outputs_found=found,
        )

    # ── Dry-run helper ────────────────────────────────────────────────────────

    def _dry_run_step(
        self,
        record:    StepExecutionRecord,
        step_dir:  Path,
        script:    Path,
        expected:  list[str],
        meta:      dict,
        pub:       BoundPublisher,
    ) -> None:
        pub.emit(
            EventType.STEP_STARTED,
            message=f"[dry-run] {record.step_id}",
            script=script.name,
            dry_run=True,
        )
        try:
            lines   = script.read_text().splitlines()
            preview = [l for l in lines if l.strip() and not l.startswith("#")][:4]
            for line in preview:
                pub.emit(EventType.STDOUT, message=line.strip())
        except Exception:
            pass

        record.status        = StepStatus.DONE
        record.stdout        = "[dry-run: no output]"
        record.exit_code     = 0
        record.outputs_found = expected

        pub.emit(EventType.STEP_COMPLETED, message=f"[dry-run] {record.step_id}", dry_run=True)

    # ── Adaptive remediation (Sprint 1c) ─────────────────────────────────────

    async def _attempt_remediation(
        self,
        record:   "StepExecutionRecord",
        step_dir: Path,
        pub:      BoundPublisher,
    ) -> bool:
        """
        Try to remediate a failed step using RemediationExecutor + adaptive reasoning.

        depth is tracked via record.n_remediations() so recursive re-runs cannot
        exceed MAX_REMEDIATION_DEPTH.

        Returns True if the step eventually ends in DONE, False otherwise.
        """
        depth = record.n_remediations()

        # depth 3+ → escalation, no more retries
        if depth >= MAX_REMEDIATION_DEPTH:
            pub.emit(
                EventType.STEP_FAILED,
                message=(
                    f"Max remediation depth ({MAX_REMEDIATION_DEPTH}) reached for "
                    f"{record.step_id} — escalating to FAILED with rich diagnosis."
                ),
                severity=EventSeverity.ERROR,
                remediation_depth=depth,
            )
            return False

        try:
            from executors.remediation_executor import RemediationExecutor

            rem = RemediationExecutor(
                workspace_path     = self.workspace_path,
                dry_run            = False,
                max_global_retries = 1,
                backend            = self._backend,
            )
            rem._bus = self._bus  # share the bus so events flow through

            # diagnose + plan
            plan = rem.diagnose_and_plan(record)
            if not plan.is_applicable:
                # fatal — no automatic fix possible
                return False

            # apply
            rem_record = rem.apply(plan, record)

            pub.emit(
                EventType.REMEDIATION_APPLIED,
                message=(
                    f"Remediation applied to {record.step_id} "
                    f"(depth {depth + 1}/{MAX_REMEDIATION_DEPTH}): "
                    f"{plan.strategy}"
                ),
                step_id=record.step_id,
                depth=depth + 1,
                files_modified=rem_record.files_modified,
            )

            # Emit a STEP_STARTED event so observers know we're retrying
            pub.emit(
                EventType.STEP_STARTED,
                message=(
                    f"Retrying {record.step_id} after remediation "
                    f"(attempt {depth + 2})"
                ),
                retry=True,
                remediation_depth=depth + 1,
            )

            # Reset record for re-run
            from datetime import datetime as _dt
            record.status        = StepStatus.RUNNING
            record.started_at    = _dt.now()
            record.stdout        = ""
            record.stderr        = ""
            record.exit_code     = None
            record.error_message = None
            record.outputs_found   = []
            record.outputs_missing = []

            # Recursive re-run — depth is now record.n_remediations()
            await self._run_step_async(record)

            return record.status == StepStatus.DONE

        except Exception as exc:
            pub.emit(
                EventType.WARNING,
                message=f"Remediation attempt failed: {exc}",
                severity=EventSeverity.WARNING,
            )
            return False

    # ── Script discovery ──────────────────────────────────────────────────────

    def _find_script(self, step_dir: Path, step_id: str | None = None) -> Path | None:
        # Prefer resume script for recoverable steps.
        if step_id and step_id in self._resumable_steps:
            resume_name = self._resumable_steps[step_id]
            resume_path = step_dir / resume_name
            if resume_path.exists():
                return resume_path
        for name in _SCRIPT_CANDIDATES:
            p = step_dir / name
            if p.exists():
                return p
        return None
