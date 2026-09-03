"""
executors/membrane_integration_runner.py
=========================================
Compile → build workspace → execute each step against real GROMACS.
Stop at the first real failure, classify it, produce membrane_build_report.yaml.

Usage (library):
    from executors.membrane_integration_runner import MembraneIntegrationRunner
    runner = MembraneIntegrationRunner("configs/membrane_orient_test.yaml")
    report = runner.run()
    # report is a dict; also written to <workspace>/membrane_build_report.yaml

Usage (CLI):
    python -m executors.membrane_integration_runner configs/membrane_orient_test.yaml
    python -m executors.membrane_integration_runner configs/membrane_orient_test.yaml \\
        --workspace /tmp/my_run --stop-after assemble_system_topology
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional

import yaml


# ── Failure classification ──────────────────────────────────────────────────
_FAILURE_CATEGORY: dict[str, str] = {
    "orient_protein":           "file_path",
    "match_box_to_bilayer":     "geometry",
    "embed_in_bilayer":         "geometry",
    "generate_protein_topology": "topology",
    "assemble_system_topology": "topology",
    "membrane_embedding":       "force_field",
    "solvate_membrane":         "gromacs_command",
    "clean_water":              "gromacs_command",
    "add_ions":                 "gromacs_command",
    "energy_minimization":      "physical_instability",
    "equilibration":            "physical_instability",
    "production_md":            "physical_instability",
}

# Steps that are considered "assembly" (no MD — fast for validation)
_FAST_STEPS = {
    "orient_protein",
    "match_box_to_bilayer",
    "embed_in_bilayer",
    "generate_protein_topology",
    "assemble_system_topology",
    "solvate_membrane",
    "clean_water",
    "add_ions",
}


class MembraneIntegrationRunner:
    """
    Execute a protein-membrane SimForge workflow step-by-step against real GROMACS.

    Parameters
    ----------
    yaml_path:      Path to the YAML config file.
    workspace_dir:  Where to build the workspace.  A temp dir is used if None.
    stop_after:     Stop execution after this step_id (inclusive).  None = run all.
    timeout_fast:   Per-step timeout (s) for fast steps (editconf, pdb2gmx, etc.).
    timeout_slow:   Per-step timeout (s) for long-running steps (shrink loop, MD).
    """

    def __init__(
        self,
        yaml_path:     str | Path,
        workspace_dir: str | Path | None = None,
        stop_after:    str | None = None,
        timeout_fast:  int = 300,
        timeout_slow:  int = 7200,
    ) -> None:
        self.yaml_path    = Path(yaml_path).resolve()
        self.workspace_dir = Path(workspace_dir).resolve() if workspace_dir else None
        self.stop_after   = stop_after
        self.timeout_fast = timeout_fast
        self.timeout_slow = timeout_slow

    # ── Public entry point ──────────────────────────────────────────────────

    def run(self) -> dict:
        """
        Compile + build workspace, execute steps, return report dict.
        Also writes <workspace>/membrane_build_report.yaml.
        """
        print(f"[runner] Compiling {self.yaml_path}")
        workspace = self._compile_and_build()
        print(f"[runner] Workspace: {workspace}")

        report = self._make_empty_report(workspace)
        steps_dir = workspace / "steps"

        step_dirs = sorted(steps_dir.iterdir()) if steps_dir.exists() else []
        for step_dir in step_dirs:
            step_id = self._step_id_from_dir(step_dir)
            if not step_id:
                continue

            print(f"\n[runner] ── {step_dir.name} ──")
            status, detail = self._run_step(step_dir, step_id)
            self._update_report(report, step_id, step_dir, status, detail)

            if status == "fail":
                cat = _FAILURE_CATEGORY.get(step_id, "unknown")
                report["first_failure_stage"]    = step_id
                report["first_failure_category"] = cat
                report["first_failure_detail"]   = detail
                print(f"[runner] FAILED at {step_id} (category: {cat})")
                print(f"[runner] Detail: {detail[:400]}")
                break

            if self.stop_after and step_id == self.stop_after:
                print(f"[runner] Stopping after {step_id} as requested")
                break

        report["final_system_ready"] = self._is_system_ready(report)
        self._enrich_report_from_artifacts(report, workspace)

        report_path = workspace / "membrane_build_report.yaml"
        report_path.write_text(yaml.dump(report, allow_unicode=True, sort_keys=False))
        print(f"\n[runner] Report: {report_path}")
        return report

    # ── Compile + build ─────────────────────────────────────────────────────

    def _compile_and_build(self) -> Path:
        from core.compiler import SimulationCompiler
        from builders.workspace_builder import WorkspaceBuilder

        # Resolve paths relative to the YAML file
        import os
        orig_cwd = os.getcwd()
        yaml_dir = self.yaml_path.parent
        os.chdir(yaml_dir)
        try:
            result = SimulationCompiler().compile(str(self.yaml_path))
        finally:
            os.chdir(orig_cwd)

        if self.workspace_dir is None:
            import tempfile
            tmpdir = tempfile.mkdtemp(prefix="simforge_membrane_")
            self.workspace_dir = Path(tmpdir)

        ws = WorkspaceBuilder().build(
            result,
            workspace_path=self.workspace_dir,
            yaml_source=str(self.yaml_path),
        )
        return ws

    # ── Step execution ──────────────────────────────────────────────────────

    def _run_step(self, step_dir: Path, step_id: str) -> tuple[str, str]:
        """
        Returns ("pass"|"fail"|"skip", detail_string).
        """
        meta_path = step_dir / "metadata.json"
        if not meta_path.exists():
            return "skip", "no metadata.json"

        meta = json.loads(meta_path.read_text())
        automation_level = meta.get("automation_level", "")
        step_type        = meta.get("step_type", "")

        if automation_level in ("manual", "guided") or step_type in ("manual",):
            print(f"  [skip] {step_id} is {automation_level or step_type} — skipping")
            return "skip", f"manual step: {automation_level or step_type}"

        script = self._find_script(step_dir)
        if script is None:
            return "skip", "no executable script found"

        timeout = self.timeout_fast if step_id in _FAST_STEPS else self.timeout_slow
        t0 = time.monotonic()

        try:
            proc = subprocess.run(
                ["bash", str(script)],
                cwd=str(step_dir),
                capture_output=True,
                text=True,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired:
            return "fail", f"timeout after {timeout}s"
        except Exception as exc:
            return "fail", str(exc)

        elapsed = time.monotonic() - t0
        combined = (proc.stdout or "") + "\n" + (proc.stderr or "")

        if proc.returncode != 0:
            last_err = (proc.stderr or "")[-800:]
            print(f"  [fail] exit={proc.returncode} in {elapsed:.1f}s")
            return "fail", f"exit={proc.returncode}\n{last_err}"

        expected = meta.get("expected_outputs", [])
        missing  = [f for f in expected if not (step_dir / f).exists()]
        if missing:
            print(f"  [fail] missing outputs: {missing}")
            return "fail", f"missing outputs: {missing}"

        print(f"  [pass] {elapsed:.1f}s  outputs={expected}")
        return "pass", combined[:300]

    def _find_script(self, step_dir: Path) -> Path | None:
        for name in ("run.sh", "run_md.sh", "run_nvt.sh"):
            s = step_dir / name
            if s.exists():
                return s
        return None

    # ── Report construction ─────────────────────────────────────────────────

    def _make_empty_report(self, workspace: Path) -> dict:
        return {
            "workspace":              str(workspace),
            "yaml_source":            str(self.yaml_path),
            "selected_bilayer":       None,
            "bilayer_x":              None,
            "bilayer_y":              None,
            "bilayer_z":              None,
            "protein_x":              None,
            "protein_y":              None,
            "protein_z":              None,
            "xy_fit_status":          None,
            "topology_status":        "not_run",
            "solvation_status":       "not_run",
            "ionization_status":      "not_run",
            "em_status":              "not_run",
            "nvt_status":             "not_run",
            "npt_status":             "not_run",
            "first_failure_stage":    None,
            "first_failure_category": None,
            "first_failure_detail":   None,
            "remediation_applied":    [],
            "final_system_ready":     False,
        }

    def _update_report(
        self,
        report:  dict,
        step_id: str,
        step_dir: Path,
        status:  str,
        detail:  str,
    ) -> None:
        _map = {
            "generate_protein_topology": "topology_status",
            "assemble_system_topology":  "topology_status",
            "solvate_membrane":          "solvation_status",
            "add_ions":                  "ionization_status",
            "energy_minimization":       "em_status",
            "equilibration":             "nvt_status",
            "production_md":             "npt_status",
        }
        key = _map.get(step_id)
        if key:
            report[key] = "pass" if status == "pass" else ("skipped" if status == "skip" else "fail")

    def _enrich_report_from_artifacts(self, report: dict, workspace: Path) -> None:
        """Read JSON artifacts produced by steps to fill geometry fields."""
        steps_dir = workspace / "steps"
        if not steps_dir.exists():
            return

        # box_match_report.json → bilayer geometry + protein geometry
        for step_dir in sorted(steps_dir.iterdir()):
            if "match_box" not in step_dir.name:
                continue
            rpt = step_dir / "box_match_report.json"
            if rpt.exists():
                try:
                    data = json.loads(rpt.read_text())
                    rb   = data.get("recommended_box", {})
                    bxy  = data.get("bilayer_xy", {})
                    pg   = data.get("protein_geometry", {})
                    meta_path = step_dir / "metadata.json"
                    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
                    report["selected_bilayer"] = meta.get("bilayer_file", "dppc512_whole.gro")
                    report["bilayer_x"]        = bxy.get("bilayer_x_nm") or rb.get("box_x_nm")
                    report["bilayer_y"]        = bxy.get("bilayer_y_nm") or rb.get("box_y_nm")
                    report["bilayer_z"]        = rb.get("box_z_nm")
                    report["protein_x"]        = pg.get("x_extent_nm")
                    report["protein_y"]        = pg.get("y_extent_nm")
                    report["protein_z"]        = pg.get("z_extent_nm")
                    report["xy_fit_status"]    = bxy.get("xy_fit_status")
                except Exception:
                    pass
            break

        # orientation_report.json → orient status already captured via step pass/fail
        # overlap_report.json → embed geometry (informational only)

    def _is_system_ready(self, report: dict) -> bool:
        """True when at least through add_ions without failures."""
        return (
            report["first_failure_stage"] is None
            and report.get("topology_status") == "pass"
            and report.get("solvation_status") in ("pass", "not_run")
            and report.get("ionization_status") in ("pass", "not_run")
        )

    # ── Utilities ───────────────────────────────────────────────────────────

    @staticmethod
    def _step_id_from_dir(step_dir: Path) -> str | None:
        """Extract step_id from e.g. '01_orient_protein' → 'orient_protein'."""
        name = step_dir.name
        parts = name.split("_", 1)
        return parts[1] if len(parts) == 2 and parts[0].isdigit() else None


# ── CLI ─────────────────────────────────────────────────────────────────────

def main() -> None:
    import argparse
    p = argparse.ArgumentParser(description="Run a SimForge membrane workflow end-to-end")
    p.add_argument("yaml",          help="Path to YAML config file")
    p.add_argument("--workspace",   default=None, help="Workspace directory (default: temp)")
    p.add_argument("--stop-after",  default=None, help="Stop after this step_id")
    p.add_argument("--timeout-fast", type=int, default=300)
    p.add_argument("--timeout-slow", type=int, default=7200)
    args = p.parse_args()

    runner = MembraneIntegrationRunner(
        yaml_path=args.yaml,
        workspace_dir=args.workspace,
        stop_after=args.stop_after,
        timeout_fast=args.timeout_fast,
        timeout_slow=args.timeout_slow,
    )
    report = runner.run()

    print("\n── membrane_build_report ───────────────────────────────────────")
    for k, v in report.items():
        if k not in ("workspace", "yaml_source", "first_failure_detail"):
            print(f"  {k:<30} {v}")
    if report.get("first_failure_detail"):
        print(f"\nFirst failure detail:\n{report['first_failure_detail'][:600]}")

    sys.exit(0 if report["final_system_ready"] else 1)


if __name__ == "__main__":
    main()
