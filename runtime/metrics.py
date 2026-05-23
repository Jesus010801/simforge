from __future__ import annotations

import asyncio
import json
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from typing import Optional

try:
    import psutil
    _HAS_PSUTIL = True
except ImportError:
    _HAS_PSUTIL = False

from runtime.events import BoundPublisher, EventType


@dataclass
class GROMACSPerformance:
    ns_per_day:   float
    hours_per_ns: float
    sampled_at:   float = field(default_factory=time.time)


@dataclass
class SystemSnapshot:
    timestamp:    float
    cpu_percent:  float
    ram_used_gb:  float
    ram_total_gb: float
    disk_used_gb: float | None  = None
    gpu_util:     int   | None  = None   # 0-100
    vram_used_mb: int   | None  = None
    vram_total_mb: int  | None  = None
    performance:  GROMACSPerformance | None = None

    def to_dict(self) -> dict:
        d = {
            "timestamp":    self.timestamp,
            "cpu_percent":  self.cpu_percent,
            "ram_used_gb":  self.ram_used_gb,
            "ram_total_gb": self.ram_total_gb,
        }
        if self.disk_used_gb is not None:
            d["disk_used_gb"] = self.disk_used_gb
        if self.gpu_util is not None:
            d["gpu_util"]      = self.gpu_util
            d["vram_used_mb"]  = self.vram_used_mb
            d["vram_total_mb"] = self.vram_total_mb
        return d


def _query_nvidia_smi() -> tuple[int, int, int] | None:
    """Returns (gpu_util%, vram_used_mb, vram_total_mb) or None."""
    if not shutil.which("nvidia-smi"):
        return None
    try:
        out = subprocess.check_output(
            ["nvidia-smi",
             "--query-gpu=utilization.gpu,memory.used,memory.total",
             "--format=csv,noheader,nounits"],
            timeout=5,
            text=True,
        ).strip().split("\n")[0]
        parts = [p.strip() for p in out.split(",")]
        return int(parts[0]), int(parts[1]), int(parts[2])
    except Exception:
        return None


class SystemMetricsCollector:
    """
    Collects CPU/RAM/GPU metrics in a background asyncio task and emits
    METRICS_SNAPSHOT events at `interval_s` seconds.

    Usage:
        collector = SystemMetricsCollector(publisher, workspace_dir, interval_s=15)
        task = asyncio.create_task(collector.run())
        ...
        task.cancel()
    """

    def __init__(
        self,
        publisher:    BoundPublisher,
        workspace_dir: "Optional[str]" = None,
        interval_s:   int = 15,
    ) -> None:
        self._pub         = publisher
        self._workspace   = workspace_dir
        self._interval    = interval_s
        self._snapshots:  list[SystemSnapshot] = []
        self._last_perf:  GROMACSPerformance | None = None

    def record_performance(self, ns_per_day: float, hours_per_ns: float) -> None:
        self._last_perf = GROMACSPerformance(ns_per_day=ns_per_day, hours_per_ns=hours_per_ns)

    async def run(self) -> None:
        while True:
            snap = self._collect()
            self._snapshots.append(snap)
            self._pub.emit(
                EventType.METRICS_SNAPSHOT,
                message=f"CPU {snap.cpu_percent:.1f}% RAM {snap.ram_used_gb:.1f}/{snap.ram_total_gb:.1f} GB",
                **snap.to_dict(),
            )
            await asyncio.sleep(self._interval)

    def _collect(self) -> SystemSnapshot:
        cpu   = 0.0
        ru_gb = 0.0
        rt_gb = 0.0
        disk  = None
        if _HAS_PSUTIL:
            cpu   = psutil.cpu_percent(interval=1)
            vm    = psutil.virtual_memory()
            ru_gb = vm.used  / 1e9
            rt_gb = vm.total / 1e9
            if self._workspace:
                try:
                    du = psutil.disk_usage(self._workspace)
                    disk = du.used / 1e9
                except Exception:
                    pass

        gpu_info = _query_nvidia_smi()
        snap = SystemSnapshot(
            timestamp    = time.time(),
            cpu_percent  = cpu,
            ram_used_gb  = ru_gb,
            ram_total_gb = rt_gb,
            disk_used_gb = disk,
            gpu_util     = gpu_info[0] if gpu_info else None,
            vram_used_mb = gpu_info[1] if gpu_info else None,
            vram_total_mb= gpu_info[2] if gpu_info else None,
            performance  = self._last_perf,
        )
        return snap

    def snapshots(self) -> list[SystemSnapshot]:
        return list(self._snapshots)
