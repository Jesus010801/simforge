"""Tests for the runtime/ package — no GROMACS required."""
from __future__ import annotations

import asyncio
import json
import time
from pathlib import Path

import pytest

from runtime.events import BoundPublisher, EventBus, EventSeverity, EventType, ExecutionEvent
from runtime.executor import RuntimeExecutor
from runtime.journal import JournalWriter
from runtime.artifacts import ArtifactRegistry, checksum, make_ref
from runtime.cache import StepCacheManager
from runtime.stream import AsyncProcessRunner
from executors.execution_state import StepStatus


# ══════════════════════════════════════════════════════════
# Events
# ══════════════════════════════════════════════════════════

class TestExecutionEvent:
    def test_round_trip(self):
        ev = ExecutionEvent(
            event_type=EventType.STEP_STARTED,
            workspace_id="ws1",
            step_id="energy_minimization",
            message="hello",
            data={"key": 42},
        )
        assert ExecutionEvent.from_dict(ev.to_dict()).message == "hello"
        assert ExecutionEvent.from_dict(ev.to_dict()).data["key"] == 42

    def test_default_severity_is_info(self):
        ev = ExecutionEvent(EventType.STDOUT, "ws", "step")
        assert ev.severity == EventSeverity.INFO


class TestEventBus:
    def test_global_handler_receives_all_events(self):
        bus = EventBus()
        received = []
        bus.subscribe(received.append)
        bus.publish(ExecutionEvent(EventType.STDOUT, "ws", "s1"))
        bus.publish(ExecutionEvent(EventType.HEARTBEAT, "ws", "s1"))
        assert len(received) == 2

    def test_typed_handler_receives_only_matching(self):
        bus = EventBus()
        received = []
        bus.subscribe(received.append, event_type=EventType.STEP_COMPLETED)
        bus.publish(ExecutionEvent(EventType.STDOUT, "ws", "s1"))
        bus.publish(ExecutionEvent(EventType.STEP_COMPLETED, "ws", "s1"))
        assert len(received) == 1
        assert received[0].event_type == EventType.STEP_COMPLETED

    def test_unsubscribe_global(self):
        bus = EventBus()
        received = []
        handler = received.append
        bus.subscribe(handler)
        bus.unsubscribe(handler)
        bus.publish(ExecutionEvent(EventType.STDOUT, "ws", "s1"))
        assert received == []

    def test_handler_exception_does_not_crash_bus(self):
        bus = EventBus()
        def bad(_): raise RuntimeError("boom")
        bus.subscribe(bad)
        bus.publish(ExecutionEvent(EventType.STDOUT, "ws", "s1"))  # must not raise

    def test_bound_publisher_emit(self):
        bus = EventBus()
        received = []
        bus.subscribe(received.append)
        pub = BoundPublisher(bus, "myws", "step_one")
        pub.emit(EventType.HEARTBEAT, message="alive")
        assert received[0].workspace_id == "myws"
        assert received[0].step_id == "step_one"
        assert received[0].message == "alive"

    def test_bound_publisher_rebind(self):
        bus = EventBus()
        received = []
        bus.subscribe(received.append)
        pub = BoundPublisher(bus, "ws", "step_a")
        pub2 = pub.rebind("step_b")
        pub2.emit(EventType.STDOUT)
        assert received[0].step_id == "step_b"
        assert received[0].workspace_id == "ws"


# ══════════════════════════════════════════════════════════
# Journal
# ══════════════════════════════════════════════════════════

class TestJournalWriter:
    def test_writes_jsonl(self, tmp_path):
        bus = EventBus()
        j   = JournalWriter(tmp_path)
        j.register(bus)
        bus.publish(ExecutionEvent(EventType.STEP_STARTED, "ws", "s1", message="go"))
        bus.publish(ExecutionEvent(EventType.STEP_COMPLETED, "ws", "s1"))
        lines = j.path().read_text().splitlines()
        assert len(lines) == 2
        assert json.loads(lines[0])["event_type"] == "STEP_STARTED"

    def test_heartbeat_not_written(self, tmp_path):
        bus = EventBus()
        j   = JournalWriter(tmp_path)
        j.register(bus)
        bus.publish(ExecutionEvent(EventType.HEARTBEAT, "ws", "s1"))
        assert not j.path().exists() or j.path().read_text().strip() == ""

    def test_stdout_not_written(self, tmp_path):
        bus = EventBus()
        j   = JournalWriter(tmp_path)
        j.register(bus)
        bus.publish(ExecutionEvent(EventType.STDOUT, "ws", "s1", message="line"))
        assert not j.path().exists() or j.path().read_text().strip() == ""

    def test_read_all_returns_events(self, tmp_path):
        bus = EventBus()
        j   = JournalWriter(tmp_path)
        j.register(bus)
        bus.publish(ExecutionEvent(EventType.STEP_STARTED, "ws", "s1"))
        bus.publish(ExecutionEvent(EventType.WARNING, "ws", "s1", message="warn"))
        events = j.read_all()
        assert len(events) == 2
        assert events[1].message == "warn"

    def test_read_all_empty_if_no_file(self, tmp_path):
        j = JournalWriter(tmp_path)
        assert j.read_all() == []


# ══════════════════════════════════════════════════════════
# Artifacts
# ══════════════════════════════════════════════════════════

class TestChecksum:
    def test_same_content_same_hash(self, tmp_path):
        f = tmp_path / "a.txt"
        f.write_text("hello")
        assert checksum(f) == checksum(f)

    def test_different_content_different_hash(self, tmp_path):
        a = tmp_path / "a.txt"; a.write_text("hello")
        b = tmp_path / "b.txt"; b.write_text("world")
        assert checksum(a) != checksum(b)

    def test_missing_file_returns_empty_string(self, tmp_path):
        assert checksum(tmp_path / "nope.txt") == ""


class TestArtifactRegistry:
    def test_register_and_lineage(self, tmp_path):
        f = tmp_path / "topol.top"; f.write_text("[ system ]")
        reg = ArtifactRegistry(tmp_path)
        ref = reg.register(f, "topology", "generate_topology")
        assert ref.semantic_role == "topology"
        assert ref.step_id == "generate_topology"
        lin = reg.lineage(f)
        assert lin is not None
        assert lin.ref.checksum == ref.checksum

    def test_record_modification(self, tmp_path):
        f = tmp_path / "topol.top"; f.write_text("v1")
        reg = ArtifactRegistry(tmp_path)
        reg.register(f, "topology", "step_a")
        f.write_text("v2")
        reg.record_modification(f, "step_b", "added ions")
        lin = reg.lineage(f)
        assert len(lin.modifications) == 1
        assert lin.modifications[0].step_id == "step_b"
        assert lin.modifications[0].description == "added ions"

    def test_persisted_and_reloaded(self, tmp_path):
        f = tmp_path / "coord.gro"; f.write_text("GROMACS gro file")
        reg = ArtifactRegistry(tmp_path)
        reg.register(f, "coordinates", "solvate")
        reg2 = ArtifactRegistry(tmp_path)
        lin = reg2.lineage(f)
        assert lin is not None
        assert lin.ref.semantic_role == "coordinates"


# ══════════════════════════════════════════════════════════
# Cache
# ══════════════════════════════════════════════════════════

class TestStepCacheManager:
    def test_miss_then_hit(self, tmp_path):
        c  = StepCacheManager(tmp_path)
        fp = c.fingerprint({"nsteps": 50000}, [])
        assert not c.is_cached("energy_minimization", fp)
        c.record("energy_minimization", fp)
        assert c.is_cached("energy_minimization", fp)

    def test_different_params_different_fingerprint(self, tmp_path):
        c   = StepCacheManager(tmp_path)
        fp1 = c.fingerprint({"nsteps": 50000}, [])
        fp2 = c.fingerprint({"nsteps": 100000}, [])
        assert fp1 != fp2

    def test_different_input_file_different_fingerprint(self, tmp_path):
        f1 = tmp_path / "a.gro"; f1.write_text("v1")
        f2 = tmp_path / "a.gro"  # same path, different content
        c  = StepCacheManager(tmp_path)
        fp1 = c.fingerprint({}, [f1])
        f2.write_text("v2")
        fp2 = c.fingerprint({}, [f2])
        assert fp1 != fp2

    def test_invalidate(self, tmp_path):
        c  = StepCacheManager(tmp_path)
        fp = c.fingerprint({}, [])
        c.record("step_x", fp)
        c.invalidate("step_x")
        assert not c.is_cached("step_x", fp)

    def test_persisted_across_instances(self, tmp_path):
        c1 = StepCacheManager(tmp_path)
        fp = c1.fingerprint({"nsteps": 1000}, [])
        c1.record("nvt", fp)
        c2 = StepCacheManager(tmp_path)
        assert c2.is_cached("nvt", fp)


# ══════════════════════════════════════════════════════════
# AsyncProcessRunner
# ══════════════════════════════════════════════════════════

class TestAsyncProcessRunner:
    def _pub(self) -> tuple[EventBus, BoundPublisher, list]:
        bus = EventBus()
        received = []
        bus.subscribe(received.append)
        pub = BoundPublisher(bus, "ws", "step")
        return bus, pub, received

    def test_simple_echo(self):
        _, pub, received = self._pub()
        runner = AsyncProcessRunner(pub)
        result = asyncio.run(runner.run(["echo", "hello world"]))
        assert result.returncode == 0
        stdout_events = [e for e in received if e.event_type == EventType.STDOUT]
        assert any("hello world" in e.message for e in stdout_events)

    def test_exit_code_nonzero(self):
        _, pub, _ = self._pub()
        runner = AsyncProcessRunner(pub)
        result = asyncio.run(runner.run(["bash", "-c", "exit 42"]))
        assert result.returncode == 42

    def test_stderr_emitted(self):
        _, pub, received = self._pub()
        runner = AsyncProcessRunner(pub)
        asyncio.run(runner.run(["bash", "-c", "echo err >&2"]))
        stderr_events = [e for e in received if e.event_type == EventType.STDERR]
        assert any("err" in e.message for e in stderr_events)

    def test_wall_time_measured(self):
        _, pub, _ = self._pub()
        runner = AsyncProcessRunner(pub)
        result = asyncio.run(runner.run(["bash", "-c", "sleep 0.1"]))
        assert result.wall_time_s >= 0.1

    def test_performance_line_parsed(self):
        _, pub, received = self._pub()
        runner = AsyncProcessRunner(pub)
        fake_output = 'echo "Performance:  42.5 ns/day  0.56 hours/ns"'
        result = asyncio.run(runner.run(["bash", "-c", fake_output]))
        assert result.ns_per_day == pytest.approx(42.5)
        perf_events = [e for e in received if e.event_type == EventType.PERFORMANCE]
        assert len(perf_events) == 1
        assert perf_events[0].data["ns_per_day"] == pytest.approx(42.5)


# ══════════════════════════════════════════════════════════
# Streaming robustness
# ══════════════════════════════════════════════════════════

class TestStreamRobustness:
    """
    AsyncProcessRunner must survive large output, oversized lines, and
    mixed streams without losing the process returncode.

    Motivation: GROMACS mdrun production runs emit log blocks that can
    exceed asyncio's default 64 KB StreamReader limit, causing
    LimitOverrunError and exit_code=None in older code.
    """

    def _pub(self) -> tuple[EventBus, BoundPublisher, list]:
        bus = EventBus()
        received = []
        bus.subscribe(received.append)
        pub = BoundPublisher(bus, "ws", "step")
        return bus, pub, received

    # ── Oversized lines ───────────────────────────────────────────────────────

    def test_line_128kb_no_crash(self):
        """A stdout line of 128 KB (2× old limit) must not raise LimitOverrunError."""
        _, pub, _ = self._pub()
        runner = AsyncProcessRunner(pub)
        result = asyncio.run(runner.run(
            ["python3", "-c", "print('x' * 131072)"]
        ))
        assert result.returncode == 0, (
            f"returncode={result.returncode} — process must complete despite large line"
        )

    def test_line_1mb_no_crash(self):
        """A 1 MB line — stress test for the 10 MB buffer limit."""
        _, pub, _ = self._pub()
        runner = AsyncProcessRunner(pub)
        result = asyncio.run(runner.run(
            ["python3", "-c", "print('y' * 1_048_576)"]
        ))
        assert result.returncode == 0

    def test_multiple_oversized_lines(self):
        """Several lines > 64 KB each — all must be handled, returncode preserved."""
        _, pub, _ = self._pub()
        runner = AsyncProcessRunner(pub)
        script = (
            "import sys\n"
            "for _ in range(5):\n"
            "    print('z' * 100_000)\n"
            "print('done')\n"
        )
        result = asyncio.run(runner.run(["python3", "-c", script]))
        assert result.returncode == 0

    # ── Returncode isolation ──────────────────────────────────────────────────

    def test_nonzero_exit_captured_after_large_output(self):
        """Large stdout must not prevent capturing a non-zero exit code."""
        _, pub, _ = self._pub()
        runner = AsyncProcessRunner(pub)
        script = "for i in $(seq 200); do echo line $i; done; exit 7"
        result = asyncio.run(runner.run(["bash", "-c", script]))
        assert result.returncode == 7, (
            f"Expected exit 7, got {result.returncode}"
        )

    def test_returncode_captured_after_large_stderr(self):
        """Large stderr output must not cause returncode to be lost."""
        _, pub, _ = self._pub()
        runner = AsyncProcessRunner(pub)
        script = "python3 -c \"import sys; [sys.stderr.write('e'*200+'\\n') for _ in range(200)]\""
        result = asyncio.run(runner.run(["bash", "-c", script]))
        assert result.returncode == 0

    # ── High-volume output ────────────────────────────────────────────────────

    def test_10000_lines_captured(self):
        """10 000 lines of output — all captured, no memory crash."""
        _, pub, _ = self._pub()
        runner = AsyncProcessRunner(pub)
        script = "import sys\nfor i in range(10_000): print(f'line {i}')"
        result = asyncio.run(runner.run(["python3", "-c", script]))
        assert result.returncode == 0
        assert len(result.stdout_lines) == 10_000

    def test_mixed_stdout_stderr_high_volume(self):
        """Interleaved stdout/stderr at high volume — no deadlock or crash."""
        _, pub, _ = self._pub()
        runner = AsyncProcessRunner(pub)
        script = (
            "import sys\n"
            "for i in range(1000):\n"
            "    print(f'out {i}')\n"
            "    sys.stderr.write(f'err {i}\\n')\n"
        )
        result = asyncio.run(runner.run(["python3", "-c", script]))
        assert result.returncode == 0
        assert len(result.stdout_lines) == 1000

    # ── Heartbeat lifecycle ───────────────────────────────────────────────────

    def test_heartbeat_emitted_for_slow_process(self):
        """Heartbeat events fire while the process is running."""
        _, pub, received = self._pub()
        runner = AsyncProcessRunner(pub, heartbeat_s=1)  # short interval for test
        result = asyncio.run(runner.run(["bash", "-c", "sleep 2.5"]))
        assert result.returncode == 0
        heartbeats = [e for e in received if e.event_type == EventType.HEARTBEAT]
        assert len(heartbeats) >= 2, (
            f"Expected ≥2 heartbeats for 2.5s run, got {len(heartbeats)}"
        )

    def test_heartbeat_does_not_prevent_completion(self):
        """A short process completes even if heartbeat interval > process duration."""
        _, pub, _ = self._pub()
        runner = AsyncProcessRunner(pub, heartbeat_s=60)  # very long interval
        result = asyncio.run(runner.run(["echo", "quick"]))
        assert result.returncode == 0

    # ── stdout_lines buffer ───────────────────────────────────────────────────

    def test_stdout_lines_accumulate_correctly(self):
        """Lines are split on newline and accumulated in stdout_lines."""
        _, pub, _ = self._pub()
        runner = AsyncProcessRunner(pub)
        result = asyncio.run(runner.run(
            ["bash", "-c", "echo alpha; echo beta; echo gamma"]
        ))
        assert result.returncode == 0
        assert "alpha" in result.stdout_lines
        assert "beta" in result.stdout_lines
        assert "gamma" in result.stdout_lines


# ══════════════════════════════════════════════════════════
# Cache integrity — RuntimeExecutor integration tests
# ══════════════════════════════════════════════════════════

def _make_fake_workspace(
    tmp_path: Path,
    step_id: str = "test_step",
    expected_outputs: list[str] | None = None,
    params: dict | None = None,
) -> Path:
    """
    Build a minimal SimForge workspace with one step whose script
    creates `expected_outputs` via `touch`. No GROMACS required.
    """
    if expected_outputs is None:
        expected_outputs = ["result.gro"]
    if params is None:
        params = {"nsteps": 1000}

    step_dir = tmp_path / "steps" / f"01_{step_id}"
    step_dir.mkdir(parents=True)
    meta_dir = tmp_path / "metadata"
    meta_dir.mkdir()

    # Script creates expected output files
    touches = "\n".join(f'touch "{o}"' for o in expected_outputs)
    (step_dir / "run.sh").write_text(f"#!/bin/bash\n{touches}\n")

    (step_dir / "metadata.json").write_text(json.dumps({
        "step_type":        "automatic",
        "blocking":         False,
        "required_inputs":  [],
        "expected_outputs": expected_outputs,
        "params":           params,
        "stage":            "minimization",
        "engine":           "test",
    }))

    (meta_dir / "execution_manifest.json").write_text(json.dumps({
        "system_type": "test",
        "steps": [
            {
                "step_id":   step_id,
                "dir_name":  f"01_{step_id}",
                "depends_on": [],
            }
        ],
    }))

    return tmp_path


class TestCacheIntegrity:
    """
    Materialized artifact validation on cache hits.

    A cache hit MUST be rejected (and the step re-run) whenever any
    expected_output is missing on disk, even if the fingerprint matches.
    """

    def _run(self, workspace: Path) -> tuple[list[ExecutionEvent], object]:
        events: list[ExecutionEvent] = []
        executor = RuntimeExecutor(workspace, dry_run=False)
        executor.bus.subscribe(events.append)
        state = executor.run()
        return events, state

    # ── Happy path: outputs present ───────────────────────────────────────────

    def test_first_run_populates_cache(self, tmp_path):
        ws = _make_fake_workspace(tmp_path)
        events, state = self._run(ws)
        assert state.steps[0].status == StepStatus.DONE
        # Second run should hit cache
        events2, state2 = self._run(ws)
        hits = [e for e in events2 if e.event_type == EventType.CACHE_HIT]
        assert len(hits) == 1, "Expected one cache hit on second run with intact outputs"

    def test_cache_hit_skips_step_when_outputs_exist(self, tmp_path):
        ws = _make_fake_workspace(tmp_path, expected_outputs=["em.gro", "em.edr", "em.log"])
        self._run(ws)  # first run — populate cache and create outputs
        events, state = self._run(ws)
        hits = [e for e in events if e.event_type == EventType.CACHE_HIT]
        assert len(hits) == 1
        assert state.steps[0].status == StepStatus.DONE

    # ── Stale cache: outputs deleted after successful run ─────────────────────

    def test_deleted_output_invalidates_cache(self, tmp_path):
        ws = _make_fake_workspace(tmp_path, expected_outputs=["result.gro"])
        self._run(ws)
        step_dir = ws / "steps" / "01_test_step"
        (step_dir / "result.gro").unlink()

        events, state = self._run(ws)

        # Must NOT be a cache hit
        hits = [e for e in events if e.event_type == EventType.CACHE_HIT]
        assert hits == [], "Cache hit must be rejected when output is missing"

        # Step must have re-run successfully and recreated the output
        assert state.steps[0].status == StepStatus.DONE
        assert (step_dir / "result.gro").exists(), "Output must be recreated after re-run"

    def test_partial_artifact_loss_invalidates_cache(self, tmp_path):
        ws = _make_fake_workspace(
            tmp_path,
            expected_outputs=["em.gro", "em.edr", "em.log"],
        )
        self._run(ws)
        step_dir = ws / "steps" / "01_test_step"
        (step_dir / "em.edr").unlink()  # delete only one of three outputs

        events, state = self._run(ws)

        hits = [e for e in events if e.event_type == EventType.CACHE_HIT]
        assert hits == [], "Partial artifact loss must also invalidate the cache"
        assert state.steps[0].status == StepStatus.DONE
        # All three outputs must exist after re-run
        for name in ("em.gro", "em.edr", "em.log"):
            assert (step_dir / name).exists()

    def test_cache_invalidated_entry_is_removed_from_disk(self, tmp_path):
        ws = _make_fake_workspace(tmp_path, expected_outputs=["result.gro"])
        self._run(ws)

        # Corrupt the artifact store
        (ws / "steps" / "01_test_step" / "result.gro").unlink()

        # Run once to trigger invalidation
        self._run(ws)

        # The cache file should no longer contain a stale entry for this step
        # (it gets replaced with a valid entry after the successful re-run)
        cache_path = ws / "metadata" / "step_cache.json"
        assert cache_path.exists()
        cache = json.loads(cache_path.read_text())
        # After re-run the entry is valid again (not absent — it was re-recorded)
        assert "test_step" in cache

    def test_stale_cache_emits_cache_miss_not_cache_hit(self, tmp_path):
        ws = _make_fake_workspace(tmp_path, expected_outputs=["result.gro"])
        self._run(ws)
        (ws / "steps" / "01_test_step" / "result.gro").unlink()

        events, _ = self._run(ws)

        miss_events = [e for e in events if e.event_type == EventType.CACHE_MISS]
        hit_events  = [e for e in events if e.event_type == EventType.CACHE_HIT]
        assert len(miss_events) >= 1, "CACHE_MISS must be emitted on invalidation"
        assert len(hit_events) == 0,  "No CACHE_HIT may be emitted when outputs are missing"

    # ── No expected_outputs declared ──────────────────────────────────────────

    def test_empty_expected_outputs_still_hits_cache(self, tmp_path):
        ws = _make_fake_workspace(tmp_path, expected_outputs=[])
        self._run(ws)  # first run — no outputs to check
        events, state = self._run(ws)
        hits = [e for e in events if e.event_type == EventType.CACHE_HIT]
        assert len(hits) == 1, "Empty expected_outputs → all 0 outputs exist → cache valid"
