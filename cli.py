#!/usr/bin/env python3
# cli.py — SimForge command-line interface

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.table import Table
from rich import box
from rich.panel import Panel
from rich.text import Text
from rich.tree import Tree

app     = Console()
cli     = typer.Typer(
    name="simforge",
    help="Molecular simulation workflow compiler.",
    add_completion=True,
)

# ── Ligand sub-app ────────────────────────────────────────────────────────────
_ligand_app = typer.Typer(
    name="ligand",
    help="Ligand preparation and parameterization utilities.",
    no_args_is_help=True,
)
cli.add_typer(_ligand_app)


# ═══════════════════════════════════════════════════════════════════════════════
# Compile helpers
# ═══════════════════════════════════════════════════════════════════════════════

_LIGAND_ROLES = {"ligand", "competitive_ligand", "substrate", "cofactor", "inhibitor"}

_STAGE_COLORS = {
    "preparation":       "cyan",
    "parametrization":   "blue",
    "assembly":          "blue",
    "minimization":      "yellow",
    "equilibration":     "yellow",
    "production":        "green",
    "analysis":          "magenta",
    "enhanced_sampling": "red",
    "membrane_embedding":"red",
    "validation":        "dim",
    "review":            "dim",
}


def _stage_ok(n: int, total: int, label: str, elapsed: float, note: str = "") -> None:
    note_str = f"  [dim]{note}[/dim]" if note else ""
    app.print(
        f"  [green]✓[/green] [dim][{n}/{total}][/dim] "
        f"{label:<38} [dim]{elapsed:.2f}s[/dim]{note_str}"
    )


def _collect_scientific_warnings(state, plan) -> list[tuple[str, str]]:
    """Returns list of (message, color) tuples."""
    warns: list[tuple[str, str]] = []

    # Membrane-associated without membrane enabled
    for comp in state.components:
        ctx = getattr(comp, "biological_context", None) or []
        if "membrane_associated" in ctx:
            mem = getattr(state.environment, "membrane", None)
            if not (mem and getattr(mem, "enabled", False)):
                warns.append((
                    f"membrane_associated protein ({comp.id}) without membrane in config",
                    "yellow",
                ))

    # High-flexibility ligands
    for comp in state.components:
        if comp.role in _LIGAND_ROLES:
            desc = getattr(comp, "descriptors", None)
            if desc and desc.flexibility_class in ("flexible", "very_flexible"):
                warns.append((
                    f"High-flexibility ligand {comp.id} "
                    f"({desc.n_rotatable_bonds} rotatable bonds) — consider longer equilibration",
                    "yellow",
                ))
            if desc and abs(desc.net_charge) > 2:
                warns.append((
                    f"Ligand {comp.id} has net charge {desc.net_charge:+d} — verify ion neutralization",
                    "yellow",
                ))

    # Blocking issues from plan
    for issue in plan.blocking_issues:
        color = "red" if issue.severity.value in ("critical", "error") else "yellow"
        warns.append((issue.message, color))

    # Enhanced sampling
    if plan.workflow_policy.enhanced_sampling:
        warns.append((
            f"Enhanced sampling: {plan.workflow_policy.sampling_method} — "
            f"verify replica count and temperature range",
            "cyan",
        ))

    # Long run
    if plan.workflow_policy.production_time_ns >= 100:
        warns.append((
            f"Long production ({plan.workflow_policy.production_time_ns} ns) — "
            f"verify disk space and GPU availability",
            "yellow",
        ))

    return warns


def _build_dag_tree(steps: list) -> Tree:
    """Build a Rich Tree from execution order using primary-parent algorithm."""
    step_ids   = [s.step_id for s in steps]
    id_to_step = {s.step_id: s for s in steps}

    # Each node's primary parent = its latest dependency in execution order
    primary_parent: dict[str, str | None] = {}
    for step in steps:
        deps = [d for d in (step.depends_on or []) if d in id_to_step]
        primary_parent[step.step_id] = (
            None if not deps
            else max(deps, key=lambda d: step_ids.index(d))
        )

    children: dict[str, list[str]] = {sid: [] for sid in step_ids}
    roots: list[str] = []
    for sid, parent in primary_parent.items():
        if parent is None:
            roots.append(sid)
        else:
            children[parent].append(sid)

    def _render(node, sid: str) -> None:
        step  = id_to_step.get(sid)
        color = _STAGE_COLORS.get(step.stage.value if step else "", "white")
        stage = f"  [dim]{step.stage.value}[/dim]" if step else ""
        child_node = node.add(f"[{color}]{sid}[/{color}]{stage}")
        for cid in children.get(sid, []):
            _render(child_node, cid)

    tree = Tree("[bold dim]Workflow DAG[/bold dim]")
    for rid in roots:
        _render(tree, rid)
    return tree


def _show_system_summary(result) -> None:
    state  = result.state
    plan   = result.plan
    policy = plan.workflow_policy

    comp_lines = []
    for c in state.components:
        desc = getattr(c, "descriptors", None)
        extra = ""
        if desc:
            if desc.n_heavy_atoms:
                extra += f"  {desc.n_heavy_atoms} heavy atoms"
            if desc.flexibility_class != "unknown":
                extra += f"  [{desc.flexibility_class}]"
            if desc.net_charge:
                extra += f"  charge {desc.net_charge:+d}"
        comp_lines.append(f"  [dim]•[/dim] {c.id}  [dim]{c.role}[/dim]{extra}")

    dur   = policy.production_time_ns
    gpu_h = max(0.1, dur * 0.2)
    cpu_h = dur * 4.0
    physio = " [dim](37 °C — physiological)[/dim]" if abs(policy.temperature_K - 309.65) < 0.1 else ""

    lines = [
        f"  System type     [bold]{state.inferred_system_type}[/bold]",
        f"  Temperature     {policy.temperature_K} K{physio}",
        f"  Pressure        {policy.pressure_bar} bar",
        f"  Duration        {dur} ns  ({int(dur / policy.timestep_ps * 1000):,} steps)",
        f"  Workflow        {len(plan.steps)} steps",
        "",
        "  [bold]Components[/bold]",
    ] + comp_lines + [
        "",
        f"  [bold]Runtime estimate[/bold] [dim](rough — depends on hardware & system size)[/dim]",
        f"  GPU ≈ {gpu_h:.1f} h     CPU ≈ {cpu_h:.0f} h",
    ]

    if policy.enhanced_sampling:
        lines.append(f"\n  [red]Enhanced sampling:[/red] {policy.sampling_method}")

    app.print(Panel("\n".join(lines), title="System Summary", border_style="cyan", padding=(0, 2)))


def _show_scientific_warnings(warnings: list[tuple[str, str]]) -> None:
    if not warnings:
        return
    table = Table(box=box.SIMPLE, show_header=False, padding=(0, 1))
    table.add_column("icon", width=2, no_wrap=True)
    table.add_column("message")
    for msg, color in warnings:
        icon = "⚠" if color == "yellow" else ("✗" if color == "red" else "◆")
        table.add_row(f"[{color}]{icon}[/{color}]", f"[{color}]{msg}[/{color}]")
    app.print(Panel(table, title="Scientific Warnings", border_style="yellow", padding=(0, 1)))


def _run_geometry_advisory(state, plan) -> "list | None":
    """
    Run geometry analysis on all protein/structure components in the state.
    Returns a list of (component_id, GeometryReport) tuples, or None if nothing
    to analyze.
    """
    from core.geometry_advisor import GeometryAdvisor

    protein_roles = {"protein", "receptor", "enzyme", "antibody", "antigen"}
    components = [
        c for c in state.components
        if c.role in protein_roles and getattr(c, "file", None)
    ]
    if not components:
        return None

    # Extract box_distance from solvate step params (or use standard default)
    padding_nm = 1.2
    box_type   = "triclinic"
    for step in getattr(plan, "steps", []):
        if "solvate" in step.step_id:
            padding_nm = step.params.get("box_distance", padding_nm)
            box_type   = step.params.get("box_type", box_type)
            break

    advisor = GeometryAdvisor()
    results = []
    for comp in components:
        pdb = Path(comp.file)
        if not pdb.exists():
            continue
        report = advisor.analyze(pdb, padding_nm=padding_nm, box_type=box_type)
        results.append((comp.id, report))
    return results or None


def _show_geometry_advisories(results: "list") -> None:
    """Display geometry advisory reports from _run_geometry_advisory()."""
    from core.geometry_advisor import Advisory

    _LEVEL_COLOR = {"INFO": "cyan", "WARNING": "yellow", "SUGGEST": "green"}
    _LEVEL_ICON  = {"INFO": "◆",    "WARNING": "⚠",      "SUGGEST": "→"}

    table = Table(box=box.SIMPLE, show_header=False, padding=(0, 1))
    table.add_column("icon",    width=2,  no_wrap=True)
    table.add_column("message", min_width=50)
    table.add_column("detail",  style="dim")

    has_warnings = False
    for comp_id, report in results:
        for adv in report.advisories:
            color  = _LEVEL_COLOR.get(adv.level, "white")
            icon   = _LEVEL_ICON.get(adv.level, "·")
            table.add_row(
                f"[{color}]{icon}[/{color}]",
                f"[{color}]{adv.message}[/{color}]",
                adv.detail.replace("\n", "  "),
            )
            if adv.level in ("WARNING", "SUGGEST"):
                has_warnings = True

    border = "yellow" if has_warnings else "cyan"
    app.print(Panel(
        table,
        title=f"Geometry Advisory — {', '.join(cid for cid, _ in results)}",
        border_style=border,
        padding=(0, 1),
    ))


def _show_parse_error(exc: Exception) -> None:
    """Display a parse error with semantic suggestions for objective errors."""
    msg = str(exc)
    body = f"[red]Parse error:[/red] {msg}\n\nCheck YAML syntax and component file paths."

    # Detect objective-related errors and add suggestions
    if "objective" in msg.lower() or "objetivo" in msg.lower():
        from core.semantic_objectives import suggest_objectives, SIMULATION_PRESETS
        import re
        quoted = re.findall(r"'([^']+)'", msg)
        if quoted:
            unknown = quoted[0]
            suggestions = suggest_objectives(unknown)
            presets = list(SIMULATION_PRESETS.keys())[:4]
            body += (
                f"\n\n[yellow]Unknown objective:[/yellow] '{unknown}'\n"
                + (f"  Closest matches: {suggestions}\n" if suggestions else "")
                + f"  Or use a preset:  simulation_profile: "
                + " | ".join(presets)
            )

    app.print(Panel(body, title="❌ Compile Failed", border_style="red"))


def _show_semantic_normalization_notes(state) -> None:
    """Show normalization notes, membrane inference, and workflow selection hints."""

    # Surface objective normalization notes from global_reasoning
    semantic_notes = [
        n for n in (state.global_reasoning.notes or [])
        if n.startswith("[semantic]")
    ]
    inference_notes = [
        n for n in (state.global_reasoning.notes or [])
        if n.startswith("[inference]")
    ]
    # Surface unknown objective warnings
    unknown_warns = [
        w for w in (state.warnings or [])
        if "Unknown simulation objective" in w.message
    ]

    is_membrane = state.inferred_system_type in (
        "protein-membrane", "protein-membrane-ligand"
    )
    hints = state.workflow_hints

    has_anything = semantic_notes or unknown_warns or inference_notes or is_membrane
    if not has_anything:
        return

    from rich.table import Table
    table = Table(box=box.SIMPLE, show_header=False, padding=(0, 1))
    table.add_column("icon", width=2, no_wrap=True)
    table.add_column("message")

    has_warn = False

    # ── Membrane workflow trazabilidad ────────────────────────────────────────
    if is_membrane:
        mem = state.environment.membrane
        lipid = (getattr(mem, "type", None) or "DPPC").upper()
        table.add_row("[green]✓[/green]", f"[green]Membrane system detected → [bold]MembraneWorkflowOPLSAA[/bold] selected[/green]")
        table.add_row("[green]✓[/green]", f"[green]Lipid environment: [bold]{lipid}[/bold][/green]")
        if hints.semiisotropic_coupling:
            table.add_row("[green]✓[/green]", "[green]Semiisotropic pressure coupling enabled[/green]")
        if hints.conservative_timestep:
            table.add_row("[green]✓[/green]", "[green]Conservative timestep (dt=0.001 ps) enabled[/green]")
        if hints.membrane_equilibration:
            table.add_row("[green]✓[/green]", "[green]Membrane equilibration protocol enabled[/green]")

    # ── Inference fallback warnings ───────────────────────────────────────────
    for note in inference_notes:
        clean = note.removeprefix("[inference] ")
        table.add_row("[yellow]⚠[/yellow]", f"[yellow]{clean}[/yellow]")
        has_warn = True

    # ── Objective normalization notes ─────────────────────────────────────────
    for note in semantic_notes:
        clean = note.removeprefix("[semantic] ")
        table.add_row("[cyan]◆[/cyan]", f"[cyan]{clean}[/cyan]")

    for w in unknown_warns:
        from core.semantic_objectives import suggest_objectives, SIMULATION_PRESETS
        import re
        quoted = re.findall(r"'([^']+)'", w.message)
        unknown = quoted[0] if quoted else "?"
        suggestions = suggest_objectives(unknown)
        presets = [k for k in SIMULATION_PRESETS if unknown.replace("_", "") in k.replace("_", "")]
        msg = f"[yellow]Unknown objective:[/yellow] '{unknown}'"
        if suggestions:
            msg += f"\n      Closest matches: {suggestions}"
        if presets:
            msg += f"\n      Or use preset:   simulation_profile: {presets[0]}"
        table.add_row("[yellow]⚠[/yellow]", msg)
        has_warn = True

    border = "yellow" if has_warn else ("green" if is_membrane else "cyan")
    title = "Semantic Normalization"
    app.print(Panel(table, title=title, border_style=border, padding=(0, 1)))


def _inspect_variants(yaml_path: Path) -> None:
    from core.variant_compiler import parse_variant_yaml

    app.print(f"\n[bold cyan]SimForge Inspect[/bold cyan]  [magenta]variants[/magenta]  "
              f"[dim]{yaml_path}[/dim]\n")

    try:
        manifest = parse_variant_yaml(yaml_path)
    except Exception as e:
        app.print(f"[red]Parse failed:[/red] {e}")
        raise typer.Exit(1)

    app.print(f"  Project:   [bold]{manifest.project_name}[/bold]")
    app.print(f"  Variants:  {len(manifest.variants)}\n")

    for spec in manifest.variants:
        exists = "[green]✓[/green]" if Path(spec.file).exists() else "[red]✗[/red]"
        label  = spec.label or spec.variant_id
        app.print(f"  {exists}  [bold]{spec.variant_id}[/bold]  [dim]({label})[/dim]")
        app.print(f"       [dim]{spec.file}[/dim]")

    sw = manifest.shared_workflow
    app.print("\n  [bold]Shared workflow:[/bold]")
    for k, v in sw.items():
        if k != "analysis":
            app.print(f"    {k}: {v}")
    if "analysis" in sw:
        app.print(f"    analysis: {sw['analysis']}")
    app.print()


def _write_compile_report(result, workspace: Path, yaml_path: Path, warnings: list) -> None:
    """Write workspace/metadata/compile_report.md."""
    from datetime import datetime

    state  = result.state
    plan   = result.plan
    policy = plan.workflow_policy

    comps_md = "\n".join(
        f"- `{c.id}` ({c.role})" + (f" — `{c.file}`" if c.file else "")
        for c in state.components
    )

    steps_md = "| # | Step | Stage | Engine |\n|---|------|-------|--------|\n"
    for i, step in enumerate(result.execution_order, 1):
        steps_md += f"| {i} | `{step.step_id}` | {step.stage.value} | {step.engine} |\n"

    warns_md = "\n".join(f"- {msg}" for msg, _ in warnings) or "None"

    env   = state.environment
    water = getattr(getattr(env, "solvent", None), "water_model", "?")
    ions  = getattr(getattr(env, "ions", None), "concentration", "?")

    report = f"""# SimForge Compile Report

**Generated:** {datetime.now().isoformat(timespec="seconds")}
**Config:** `{yaml_path}`
**SimForge version:** 0.1.0

---

## System

| Field | Value |
|-------|-------|
| Type | {state.inferred_system_type} |
| Temperature | {policy.temperature_K} K |
| Pressure | {policy.pressure_bar} bar |
| Duration | {policy.production_time_ns} ns |
| Water model | {water} |
| Ion concentration | {ions} M |
| Hardware | {getattr(policy, "hardware", "auto")} |
| Enhanced sampling | {"yes — " + policy.sampling_method if policy.enhanced_sampling else "no"} |

## Components

{comps_md}

## Workflow ({len(result.execution_order)} steps)

{steps_md}
## Scientific Warnings

{warns_md}

## Reproducibility

- Config snapshot: `{yaml_path}`
- Execution manifest: `metadata/execution_manifest.json`
- Planning session: `metadata/planning_session.json`
"""

    (workspace / "metadata" / "compile_report.md").write_text(report)


# ═══════════════════════════════════════════════════════════════════════════════
# compile
# ═══════════════════════════════════════════════════════════════════════════════

@cli.command()
def compile(
    yaml_path:  Path = typer.Argument(..., help="YAML config file to compile."),
    output_dir: str  = typer.Option("simforge_runs", "--output-dir", "-o", help="Root output directory."),
    no_build:   bool = typer.Option(False, "--no-build", help="Only compile — skip workspace materialization."),
    no_plan:    bool = typer.Option(False, "--no-plan",  help="Skip scientific planning dialogue."),
):
    """Compile a YAML config into an executable workspace."""
    from core.compiler import SimulationCompiler
    from core.parser import parse_yaml
    from builders.workspace_builder import WorkspaceBuilder

    if not yaml_path.exists():
        app.print(Panel(
            f"[red]File not found:[/red] {yaml_path}\n\n"
            f"Suggested fix: check the path or run [dim]simforge init[/dim]",
            title="❌ Compile Failed", border_style="red",
        ))
        raise typer.Exit(1)

    # ── Variants dispatch ─────────────────────────────────────────────────────
    try:
        from core.variant_compiler import is_variant_yaml
        if is_variant_yaml(yaml_path):
            _compile_variants(yaml_path, output_dir, no_build)
            return
    except Exception:
        pass  # not a variants YAML — continue normally

    t_start = time.perf_counter()
    app.print(Panel(
        f"[bold cyan]SimForge Compiler[/bold cyan]  [dim]{yaml_path}[/dim]",
        border_style="cyan", padding=(0, 2),
    ))

    # ── [1/6] Parse YAML ──────────────────────────────────────────────────────
    t = time.perf_counter()
    try:
        with app.status("  [dim][1/6][/dim] Parsing YAML..."):
            state = parse_yaml(str(yaml_path))
    except Exception as e:
        _show_parse_error(e)
        raise typer.Exit(1)
    _stage_ok(1, 6, "Parsing YAML", time.perf_counter() - t,
              note=f"{len(state.components)} component(s)  system: {state.inferred_system_type}")
    _show_semantic_normalization_notes(state)

    # ── [2/6] Scientific planning ──────────────────────────────────────────────
    t = time.perf_counter()
    session_record = _run_planning_dialogue(state, yaml_path, skip=no_plan)
    if session_record is None:
        raise typer.Exit(1)
    _stage_ok(2, 6, "Scientific planning", time.perf_counter() - t,
              note="skipped" if (no_plan or session_record.skipped) else
                   f"{len(session_record.answers or [])} answer(s) applied")

    # ── [3/6] Build execution plan ────────────────────────────────────────────
    t = time.perf_counter()
    try:
        with app.status("  [dim][3/6][/dim] Building execution plan..."):
            result = SimulationCompiler().compile_from_state(state)
    except Exception as e:
        app.print(Panel(
            f"[red]Compilation error:[/red] {e}",
            title="❌ Compile Failed", border_style="red",
        ))
        raise typer.Exit(1)
    _stage_ok(3, 6, "Building execution plan", time.perf_counter() - t,
              note=f"{len(result.plan.steps)} steps planned")

    # ── [4/6] Workflow DAG + warnings ─────────────────────────────────────────
    t = time.perf_counter()
    warnings = _collect_scientific_warnings(result.state, result.plan)
    _stage_ok(4, 6, "Computing workflow DAG", time.perf_counter() - t,
              note=f"{len(result.execution_order)} nodes  "
                   f"{len(warnings)} warning(s)  "
                   f"{len(result.plan.blocking_issues)} blocking issue(s)")

    # ── Geometry advisory (between DAG and workspace build) ───────────────────
    geo_results = _run_geometry_advisory(result.state, result.plan)

    # ── [5/6] Materialize workspace ───────────────────────────────────────────
    workspace = None
    if not no_build:
        from core.project_manager import ProjectManager
        from datetime import datetime as _dt

        project_name = (
            result.state.inferred_system_type
            or Path(yaml_path).stem
            or "simforge_system"
        )
        project_dir = ProjectManager.project_dir(output_dir, project_name)
        run_ts      = _dt.now().strftime("%Y-%m-%d_%H-%M-%S")
        run_dir     = ProjectManager.create_run_dir(project_dir, run_ts)

        t = time.perf_counter()
        try:
            with app.status("  [dim][5/6][/dim] Materializing workspace..."):
                workspace = WorkspaceBuilder().build(
                    result,
                    output_dir=output_dir,
                    yaml_source=str(Path(yaml_path).resolve()),
                    workspace_path=run_dir,
                )
        except Exception as e:
            app.print(Panel(
                f"[red]Workspace error:[/red] {e}",
                title="❌ Compile Failed", border_style="red",
            ))
            raise typer.Exit(1)
        _stage_ok(5, 6, "Materializing workspace", time.perf_counter() - t,
                  note=str(workspace))

        # ── [6/6] Write reports ───────────────────────────────────────────────
        t = time.perf_counter()
        with app.status("  [dim][6/6][/dim] Writing reports..."):
            _save_planning_session(session_record, workspace / "metadata" / "planning_session.json")
            _write_compile_report(result, workspace, yaml_path, warnings)
            # Write run provenance into the run directory itself
            run_info = {
                "run_id":       run_ts,
                "run_path":     str(workspace.resolve()),
                "compiled_at":  run_ts,
                "yaml_source":  str(Path(yaml_path).resolve()),
                "system_type":  result.state.inferred_system_type,
                "n_steps":      len(result.execution_order),
            }
            (workspace / "metadata" / "run_info.json").write_text(
                __import__("json").dumps(run_info, indent=4)
            )
            ProjectManager.update_project_registry(project_dir, run_info)
        _stage_ok(6, 6, "Writing reports", time.perf_counter() - t,
                  note="compile_report.md  execution_manifest.json")
    else:
        app.print("  [dim]--no-build: steps 5–6 skipped[/dim]")

    # ── System summary ────────────────────────────────────────────────────────
    app.print()
    _show_system_summary(result)

    # ── Scientific warnings ───────────────────────────────────────────────────
    _show_scientific_warnings(warnings)

    # ── Geometry advisory ─────────────────────────────────────────────────────
    if geo_results:
        _show_geometry_advisories(geo_results)

    # ── DAG preview ───────────────────────────────────────────────────────────
    app.print(Panel(
        _build_dag_tree(result.execution_order),
        title="Workflow DAG", border_style="dim", padding=(0, 2),
    ))

    # ── Blocking issues ───────────────────────────────────────────────────────
    if result.plan.blocking_issues:
        app.print("[red]Blocking issues:[/red]")
        for issue in result.plan.blocking_issues:
            app.print(f"  [red]✗[/red] [{issue.severity.value}] {issue.source}: {issue.message}")
        app.print()

    # ── Final panel ───────────────────────────────────────────────────────────
    total = time.perf_counter() - t_start
    if workspace:
        prior_runs = ProjectManager.get_run_history(project_dir)
        run_history = (
            f"  Run #{len(prior_runs)} of {len(prior_runs)} in project\n"
            if len(prior_runs) > 1
            else ""
        )
        app.print(Panel(
            f"[bold green]✓[/bold green] Project  → [dim]{project_dir}[/dim]\n"
            f"  Run      → [cyan]{workspace}[/cyan]\n"
            f"{run_history}"
            f"  Report   → [dim]{workspace}/metadata/compile_report.md[/dim]\n"
            f"  Manifest → [dim]{workspace}/metadata/execution_manifest.json[/dim]\n\n"
            f"Compiled in [bold]{total:.2f}s[/bold]\n\n"
            f"Next:  [dim]simforge run {workspace}[/dim]",
            title="✓ Done", border_style="green", padding=(0, 2),
        ))
    else:
        app.print(f"\n[dim]--no-build: workspace not materialized  ({total:.2f}s)[/dim]\n")


# ═══════════════════════════════════════════════════════════════════════════════
# variants
# ═══════════════════════════════════════════════════════════════════════════════

def _compile_variants(yaml_path: Path, output_dir: str, no_build: bool) -> None:
    """Handle a variants YAML: compile each variant + build comparative workspace."""
    from core.variant_compiler import (
        parse_variant_yaml, compile_variants, build_comparative_workspace,
    )

    app.print(Panel(
        f"[bold cyan]SimForge Variants[/bold cyan]  [dim]{yaml_path}[/dim]",
        border_style="magenta", padding=(0, 2),
    ))

    try:
        manifest = parse_variant_yaml(yaml_path)
    except Exception as e:
        app.print(f"[red]Variants parse failed:[/red] {e}")
        raise typer.Exit(1)

    app.print(f"  Project:  [bold]{manifest.project_name}[/bold]")
    app.print(f"  Variants: {len(manifest.variants)}")
    for spec in manifest.variants:
        label   = spec.label or spec.variant_id
        exists  = "[green]✓[/green]" if Path(spec.file).exists() else "[red]✗[/red]"
        app.print(f"    {exists}  [cyan]{spec.variant_id}[/cyan]  ({label})  [dim]{spec.file}[/dim]")
    app.print()

    with app.status("  Compiling variants..."):
        manifest = compile_variants(manifest, output_dir=output_dir, no_build=no_build)

    app.print()
    for spec in manifest.variants:
        if spec.variant_id in manifest.errors:
            app.print(f"  [red]✗[/red]  {spec.variant_id}  — {manifest.errors[spec.variant_id]}")
        elif spec.variant_id in manifest.workspaces:
            app.print(
                f"  [green]✓[/green]  [bold]{spec.variant_id}[/bold]  →  "
                f"[cyan]{manifest.workspaces[spec.variant_id]}[/cyan]"
            )

    if not no_build and manifest.workspaces:
        comp_dir = build_comparative_workspace(manifest, output_dir=output_dir)
        app.print(Panel(
            f"  Comparative workspace:  [cyan]{comp_dir}[/cyan]\n\n"
            f"  Run each variant with:  [dim]simforge run <workspace>[/dim]\n"
            f"  Then compare results:   [dim]see {comp_dir}/README.md[/dim]",
            title="✓ Variants compiled",
            border_style="magenta",
        ))


# ═══════════════════════════════════════════════════════════════════════════════
# run
# ═══════════════════════════════════════════════════════════════════════════════

def _execute_workspace(
    workspace:     Path,
    dry_run:       bool,
    executor_type: str,
    no_confirm:    bool = False,
) -> None:
    """Shared execution logic for `run` and `dry-run` commands."""
    if not workspace.exists():
        app.print(f"[red]Error:[/red] Workspace not found: {workspace}")
        raise typer.Exit(1)

    manifest_file = workspace / "metadata" / "execution_manifest.json"
    if not manifest_file.exists():
        app.print(f"[red]Error:[/red] No execution manifest found. Run [cyan]simforge compile[/cyan] first.")
        raise typer.Exit(1)

    if dry_run:
        app.print(f"\n[bold cyan]SimForge[/bold cyan]  [yellow]dry-run[/yellow]  "
                  f"[green]{workspace}[/green]\n")
    else:
        app.print(f"\n[bold cyan]SimForge[/bold cyan]  [bold]run[/bold]  "
                  f"[green]{workspace}[/green]\n")

    # Pre-run summary for real execution
    if not dry_run:
        _print_prerun_summary(workspace, manifest_file)
        if not no_confirm:
            confirmed = typer.confirm("\nProceed with real execution?")
            if not confirmed:
                app.print("[dim]Aborted.[/dim]")
                raise typer.Exit(0)

    if executor_type == "gromacs":
        from executors.gromacs_executor import GROMACSExecutor
        executor = GROMACSExecutor(workspace, dry_run=dry_run)
    else:
        from runtime.executor import RuntimeExecutor
        executor = RuntimeExecutor(workspace, dry_run=dry_run)

    app.print()

    try:
        state = executor.run()
    except Exception as e:
        app.print(f"[red]Execution error:[/red] {e}")
        raise typer.Exit(1)

    done    = sum(1 for s in state.steps if s.status.value == "done")
    failed  = sum(1 for s in state.steps if s.status.value == "failed")
    skipped = sum(1 for s in state.steps if s.status.value in ("skipped", "blocked"))
    total   = len(state.steps)

    app.print()
    if failed == 0:
        app.print(Panel(
            f"  [green]✓[/green]  {done}/{total} steps completed  "
            f"[dim](skipped: {skipped})[/dim]",
            border_style="green",
            title="Run complete",
        ))
    else:
        failed_lines = "\n".join(
            f"  [red]·[/red]  {s.step_id}  exit={s.exit_code}"
            for s in state.steps if s.status.value == "failed"
        )
        app.print(Panel(
            f"  [red]✗[/red]  {failed} step(s) failed  |  {done} done  |  {skipped} skipped\n"
            + failed_lines,
            border_style="red",
            title="Run failed",
        ))

    if not dry_run and done > 0:
        _print_output_paths(workspace)

    if failed:
        raise typer.Exit(1)


def _print_prerun_summary(workspace: Path, manifest_file: Path) -> None:
    """Show estimated cost/time before a real run."""
    manifest    = json.loads(manifest_file.read_text())
    n_steps     = manifest.get("n_steps", "?")
    system_type = manifest.get("system_type", "unknown")

    runtime_est = None
    report_path = workspace / "metadata" / "compile_report.md"
    if report_path.exists():
        for line in report_path.read_text().splitlines():
            if line.startswith(">") and "estimate" in line.lower():
                runtime_est = line.lstrip("> ").strip()
                break

    lines = [
        f"  Workspace     [cyan]{workspace}[/cyan]",
        f"  System type   {system_type}",
        f"  Steps         {n_steps}",
    ]
    if runtime_est:
        lines.append(f"  Runtime est.  {runtime_est}")
    lines += [
        "",
        "  [dim]This will execute real GROMACS commands and may run for hours.[/dim]",
    ]
    app.print(Panel("\n".join(lines), title="Pre-run summary", border_style="yellow"))


def _print_output_paths(workspace: Path) -> None:
    """Show where outputs landed after a run."""
    app.print("[bold]Results:[/bold]")
    app.print(f"  Workspace:  [cyan]{workspace}[/cyan]")

    steps_dir = workspace / "steps"
    if steps_dir.exists():
        for step_dir in sorted(steps_dir.iterdir(), reverse=True):
            if not step_dir.is_dir():
                continue
            xtc = step_dir / "md.xtc"
            gro = step_dir / "md.gro"
            if xtc.exists():
                app.print(f"  Trajectory: [cyan]{xtc}[/cyan]")
            if gro.exists():
                app.print(f"  Structure:  [cyan]{gro}[/cyan]")
            if xtc.exists() or gro.exists():
                break

    analysis_dir = workspace / "analysis"
    if analysis_dir.exists():
        app.print(f"  Analysis:   [cyan]{analysis_dir}[/cyan]")

    report = workspace / "metadata" / "compile_report.md"
    if report.exists():
        app.print(f"  Report:     [cyan]{report}[/cyan]")
    app.print()


@cli.command()
def run(
    workspace:     Path = typer.Argument(..., help="Workspace directory to execute."),
    executor_type: str  = typer.Option("shell", "--executor", "-e", help="Executor type: shell | gromacs."),
    no_confirm:    bool = typer.Option(False, "--no-confirm", help="Skip confirmation prompt (for automation)."),
):
    """Execute a compiled workspace with real GROMACS commands."""
    _execute_workspace(workspace, dry_run=False, executor_type=executor_type, no_confirm=no_confirm)


@cli.command(name="dry-run")
def dry_run_cmd(
    workspace:     Path = typer.Argument(..., help="Workspace directory to simulate."),
    executor_type: str  = typer.Option("shell", "--executor", "-e", help="Executor type: shell | gromacs."),
):
    """Simulate execution without running real GROMACS commands."""
    _execute_workspace(workspace, dry_run=True, executor_type=executor_type)


# ═══════════════════════════════════════════════════════════════════════════════
# validate
# ═══════════════════════════════════════════════════════════════════════════════

@cli.command()
def validate(
    yaml_path: Path = typer.Argument(..., help="YAML config to validate."),
):
    """Parse and validate a config without generating a workspace."""
    from core.parser import parse_yaml
    from core.compiler import SimulationCompiler

    if not yaml_path.exists():
        app.print(f"[red]Error:[/red] File not found: {yaml_path}")
        raise typer.Exit(1)

    app.print(f"\n[bold cyan]SimForge Validate[/bold cyan]  [dim]{yaml_path}[/dim]\n")

    try:
        with app.status("  Parsing..."):
            state = parse_yaml(str(yaml_path))
        app.print(f"  [green]✓[/green] Parse        {len(state.components)} component(s)  system: {state.inferred_system_type}")
    except Exception as e:
        app.print(Panel(
            f"[red]Parse failed:[/red] {e}",
            title="❌ Validation Failed", border_style="red",
        ))
        raise typer.Exit(1)

    # Show per-component validation issues immediately (file not found, dir errors)
    comp_issues = [
        c for c in state.components
        if c.validation and (not c.validation.is_valid or c.validation.validation_error)
    ]
    for comp in comp_issues:
        app.print(Panel(
            f"[red]✗[/red]  [bold]{comp.id}[/bold]\n"
            f"  {comp.validation.validation_error}",
            title=f"Component issue: {comp.id}",
            border_style="red",
        ))

    try:
        with app.status("  Compiling plan..."):
            result = SimulationCompiler().compile_from_state(state)
        app.print(f"  [green]✓[/green] Plan         {len(result.plan.steps)} steps  {len(result.execution_order)} in DAG")
    except Exception as e:
        app.print(Panel(
            f"[red]Compilation failed:[/red] {e}",
            title="❌ Validation Failed", border_style="red",
        ))
        raise typer.Exit(1)

    warnings = _collect_scientific_warnings(result.state, result.plan)
    _show_scientific_warnings(warnings)

    if result.plan.blocking_issues:
        app.print(f"\n  [red]✗[/red] {len(result.plan.blocking_issues)} blocking issue(s):")
        for issue in result.plan.blocking_issues:
            app.print(f"    [red]•[/red] [{issue.severity.value}] {issue.message}")
        raise typer.Exit(1)

    app.print(
        f"\n  [bold green]✓ Valid[/bold green]  "
        f"{len(warnings)} warning(s)  0 blocking issues\n"
    )


# ═══════════════════════════════════════════════════════════════════════════════
# inspect
# ═══════════════════════════════════════════════════════════════════════════════

@cli.command()
def inspect(
    yaml_path: Path = typer.Argument(..., help="YAML config to inspect."),
):
    """Show parsed IR, descriptors, and DAG without writing any files."""
    from core.parser import parse_yaml
    from core.compiler import SimulationCompiler

    if not yaml_path.exists():
        app.print(f"[red]Error:[/red] File not found: {yaml_path}")
        raise typer.Exit(1)

    # Variants dispatch
    try:
        from core.variant_compiler import is_variant_yaml, parse_variant_yaml
        if is_variant_yaml(yaml_path):
            _inspect_variants(yaml_path)
            return
    except Exception:
        pass

    app.print(f"\n[bold cyan]SimForge Inspect[/bold cyan]  [dim]{yaml_path}[/dim]\n")

    try:
        with app.status("  Parsing + compiling..."):
            state  = parse_yaml(str(yaml_path))
            result = SimulationCompiler().compile_from_state(state)
    except Exception as e:
        app.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(1)

    # System summary
    _show_system_summary(result)

    # Component descriptors table
    desc_table = Table(box=box.SIMPLE, show_header=True, header_style="bold dim")
    desc_table.add_column("Component")
    desc_table.add_column("Role")
    desc_table.add_column("Heavy atoms", justify="right")
    desc_table.add_column("Flexibility")
    desc_table.add_column("Charge", justify="right")
    desc_table.add_column("H-bond D/A", justify="right")
    for comp in state.components:
        d = getattr(comp, "descriptors", None)
        desc_table.add_row(
            comp.id,
            comp.role,
            str(d.n_heavy_atoms) if d and d.n_heavy_atoms else "—",
            d.flexibility_class if d else "—",
            f"{d.net_charge:+d}" if d else "—",
            f"{d.hb_donors}/{d.hb_acceptors}" if d else "—",
        )
    app.print(Panel(desc_table, title="Component Descriptors", border_style="dim"))

    # DAG preview
    app.print(Panel(
        _build_dag_tree(result.execution_order),
        title="Workflow DAG", border_style="dim", padding=(0, 2),
    ))

    # Key step params (IR snapshot)
    params_table = Table(box=box.SIMPLE, show_header=True, header_style="bold dim")
    params_table.add_column("Step", style="white")
    params_table.add_column("Key IR params", overflow="fold", style="dim")
    _KEY_PARAMS = {"emtol", "nsteps", "temperature", "box_type", "box_distance", "hardware", "concentration"}
    for step in result.execution_order:
        shown = {k: v for k, v in (step.params or {}).items() if k in _KEY_PARAMS}
        if shown:
            params_table.add_row(step.step_id, "  ".join(f"{k}={v}" for k, v in shown.items()))
    app.print(Panel(params_table, title="Step Parameters (IR snapshot)", border_style="dim"))

    # Warnings
    _show_scientific_warnings(_collect_scientific_warnings(state, result.plan))
    app.print()


# ═══════════════════════════════════════════════════════════════════════════════
# clean / recompile
# ═══════════════════════════════════════════════════════════════════════════════

@cli.command()
def clean(
    workspace: Path = typer.Argument(..., help="Workspace directory to clean."),
    yes:       bool = typer.Option(False, "--yes", "-y", help="Skip confirmation."),
):
    """Delete stale scripts and execution state from a workspace."""
    import shutil

    if not workspace.exists():
        app.print(f"[red]Error:[/red] Workspace not found: {workspace}")
        raise typer.Exit(1)

    steps_dir  = workspace / "steps"
    state_file = workspace / "execution_state.json"

    targets = [p for p in [steps_dir, state_file] if p.exists()]
    if not targets:
        app.print(f"  [dim]Nothing to clean in {workspace}[/dim]")
        return

    app.print(f"\n[bold cyan]SimForge Clean[/bold cyan]  [dim]{workspace}[/dim]\n")
    for t in targets:
        app.print(f"  [dim]Will delete:[/dim] {t.name}")

    if not yes:
        confirm = typer.confirm("\n  Proceed?")
        if not confirm:
            app.print("  Aborted.")
            return

    for t in targets:
        if t.is_dir():
            shutil.rmtree(t)
        else:
            t.unlink()
        app.print(f"  [green]✓[/green] Deleted {t.name}")

    app.print(f"\n  Workspace cleaned. Run [dim]simforge recompile <yaml>[/dim] to regenerate.\n")


@cli.command()
def recompile(
    yaml_path:  Path = typer.Argument(..., help="Original YAML config."),
    output_dir: str  = typer.Option("simforge_runs", "--output-dir", help="Root output directory."),
    yes:        bool = typer.Option(False, "--yes", "-y", help="Skip clean confirmation."),
):
    """
    Force-regenerate all workspace scripts from the current builders.

    Equivalent to: simforge clean <workspace> && simforge compile <yaml>
    """
    import shutil
    from core.parser import parse_yaml
    from core.compiler import SimulationCompiler
    from builders.workspace_builder import WorkspaceBuilder

    if not yaml_path.exists():
        app.print(f"[red]Error:[/red] YAML not found: {yaml_path}")
        raise typer.Exit(1)

    # Infer workspace path by compiling first (just to get the name)
    try:
        with app.status("  Parsing config..."):
            state = parse_yaml(str(yaml_path))
    except Exception as e:
        app.print(f"[red]Parse error:[/red] {e}")
        raise typer.Exit(1)

    system_name   = state.inferred_system_type or "simforge_system"
    workspace_dir = Path(output_dir) / system_name

    app.print(f"\n[bold cyan]SimForge Recompile[/bold cyan]  [dim]{yaml_path}[/dim]\n")
    app.print(f"  Workspace: [cyan]{workspace_dir}[/cyan]")

    # Clean existing steps + state
    steps_dir  = workspace_dir / "steps"
    state_file = workspace_dir / "execution_state.json"
    to_delete  = [p for p in [steps_dir, state_file] if p.exists()]

    if to_delete and not yes:
        for p in to_delete:
            app.print(f"  [dim]Will delete:[/dim] {p.name}")
        if not typer.confirm("\n  Proceed with clean?"):
            app.print("  Aborted.")
            return

    for p in to_delete:
        if p.is_dir():
            shutil.rmtree(p)
        else:
            p.unlink()

    # Full recompile
    try:
        with app.status("  Compiling..."):
            result = SimulationCompiler().compile(str(yaml_path))
        with app.status("  Building workspace..."):
            workspace = WorkspaceBuilder().build(
                result,
                output_dir=output_dir,
                yaml_source=str(yaml_path.resolve()),
            )
    except Exception as e:
        app.print(f"[red]Recompile failed:[/red] {e}")
        raise typer.Exit(1)

    from core.workspace_fingerprint import compute_builder_signature
    app.print(f"\n  [green]✓[/green] Recompiled  {len(result.execution_order)} steps")
    app.print(f"  builder_signature: [dim]{compute_builder_signature()}[/dim]")
    app.print(f"  Workspace: [cyan]{workspace}[/cyan]")
    app.print(f"\nRun:  [dim]simforge run {workspace}[/dim]\n")


# ═══════════════════════════════════════════════════════════════════════════════
# status
# ═══════════════════════════════════════════════════════════════════════════════

_STATUS_STYLE = {
    "done":    ("green",  "✓"),
    "failed":  ("red",    "✗"),
    "running": ("yellow", "▶"),
    "skipped": ("dim",    "–"),
    "blocked": ("dim",    "⊘"),
    "pending": ("dim",    "·"),
}

@cli.command()
def status(
    workspace: Path = typer.Argument(Path("."), help="Workspace directory to inspect (default: current directory)."),
    verbose: bool   = typer.Option(False, "--verbose", "-v", help="Show timing and error details."),
):
    """Show execution status of a workspace."""

    state_file = workspace / "execution_state.json"
    manifest_file = workspace / "metadata" / "execution_manifest.json"

    if not state_file.exists():
        if manifest_file.exists():
            manifest = json.loads(manifest_file.read_text())
            app.print(f"\n[cyan]{workspace}[/cyan]  [dim]compiled, not yet executed[/dim]")
            app.print(f"  {manifest['n_steps']} steps  system: {manifest.get('system_type', '?')}")
            # Freshness check
            from core.workspace_fingerprint import check_workspace_freshness
            is_fresh, msg = check_workspace_freshness(manifest_file)
            if not is_fresh and "Legacy" not in msg:
                app.print(f"\n  [yellow]⚠[/yellow]  [bold]Stale workspace[/bold] — builders have changed since compilation.")
                app.print(f"  [dim]Run: simforge recompile {manifest.get('yaml_source', '<yaml>')}[/dim]")
            app.print(f"\nRun:  [dim]simforge run {workspace}[/dim]\n")
        else:
            app.print(f"[red]Error:[/red] No state or manifest found in {workspace}")
            raise typer.Exit(1)
        return

    raw   = json.loads(state_file.read_text())
    steps = raw.get("steps", [])

    # Enrich steps with stage info from manifest (not stored in execution_state)
    if manifest_file.exists():
        manifest_steps = {e["step_id"]: e for e in json.loads(manifest_file.read_text()).get("steps", [])}
        for s in steps:
            if s["step_id"] in manifest_steps:
                s.setdefault("stage", manifest_steps[s["step_id"]].get("stage", ""))

    done    = sum(1 for s in steps if s["status"] == "done")
    failed  = sum(1 for s in steps if s["status"] == "failed")
    total   = len(steps)
    dry_run = raw.get("dry_run", False)

    header = f"[bold cyan]{workspace}[/bold cyan]"
    if dry_run:
        header += "  [dim](dry-run)[/dim]"
    if raw.get("is_complete"):
        header += "  [green]complete[/green]"
    elif raw.get("was_interrupted"):
        header += "  [yellow]interrupted[/yellow]"

    app.print(f"\n{header}\n")

    table = Table(box=box.SIMPLE, show_header=True, header_style="bold")
    table.add_column("",        width=2)
    table.add_column("Step",    style="white")
    table.add_column("Stage",   style="cyan", no_wrap=True)
    table.add_column("Status",  no_wrap=True)
    if verbose:
        table.add_column("Time", style="dim", justify="right")
        table.add_column("Note", style="dim")

    for s in steps:
        status_val = s["status"]
        color, icon = _STATUS_STYLE.get(status_val, ("dim", "?"))
        status_text = Text(f"{icon} {status_val}", style=color)

        if verbose:
            elapsed = f"{s.get('elapsed_s', 0):.1f}s" if s.get("elapsed_s") else ""
            note    = s.get("error_message", "")[:60] if status_val == "failed" else ""
            table.add_row("", s["step_id"], s.get("stage", ""), status_text, elapsed, note)
        else:
            table.add_row("", s["step_id"], s.get("stage", ""), status_text)

    app.print(table)

    # Summary line
    pct = int(done / total * 100) if total else 0
    bar_filled = int(pct / 5)
    bar = "█" * bar_filled + "░" * (20 - bar_filled)
    app.print(f"  [{bar}] {done}/{total} steps  ({pct}%)")
    if failed:
        app.print(f"\n  [red]{failed} step(s) failed.[/red]  Run [dim]simforge status {workspace} --verbose[/dim] for details.\n")
    else:
        app.print()

    # Orientation gate report (shown when orient_protein has been executed)
    try:
        from runtime.orientation_gate import read_orientation_report
        orient = read_orientation_report(workspace)
        if orient is not None:
            conf    = orient.get("confidence", 0.0)
            passed  = orient.get("passed", False)
            warns   = orient.get("warnings", [])
            errors  = orient.get("errors",   [])
            status_color = "green" if passed else "red"
            status_label = "yes" if passed else "no"
            lines = [
                f"  confidence: {conf:.2f}   "
                f"passed: [{status_color}]{status_label}[/{status_color}]",
            ]
            for w in warns:
                lines.append(f"  [yellow]⚠[/yellow]  {w}")
            for e in errors:
                lines.append(f"  [red]✗[/red]  {e}")
            app.print(Panel(
                "\n".join(lines),
                title="Orientation",
                border_style="yellow" if (warns or errors) else "green",
                padding=(0, 2),
            ))
    except Exception:
        pass

    # Box-match report (shown when match_box_to_bilayer has been executed)
    try:
        from runtime.box_match_gate import read_box_match_report
        bm = read_box_match_report(workspace)
        if bm is not None:
            rb      = bm.get("recommended_box", {})
            est     = bm.get("estimates", {})
            passed  = bm.get("passed", False)
            conf    = bm.get("confidence", 0.0)
            warns   = bm.get("warnings", [])
            errors  = bm.get("errors",   [])
            sc      = "green" if passed else "red"
            sl      = "yes" if passed else "no"
            lines   = [
                f"  confidence: {conf:.2f}   passed: [{sc}]{sl}[/{sc}]"
                f"   lipid: {bm.get('lipid_type', '?')}",
            ]
            if rb:
                lines.append(
                    f"  box: {rb.get('box_x_nm', '?'):.3f} × "
                    f"{rb.get('box_y_nm', '?'):.3f} × "
                    f"{rb.get('box_z_nm', '?'):.3f} nm"
                )
            if est.get("n_lipids_estimate"):
                coverage_pct = int(est.get("protein_xy_coverage", 0) * 100)
                lines.append(
                    f"  N lipids estimate: {est['n_lipids_estimate']}   "
                    f"protein XY coverage: {coverage_pct}%"
                )
            for w in warns:
                lines.append(f"  [yellow]⚠[/yellow]  {w}")
            for e in errors:
                lines.append(f"  [red]✗[/red]  {e}")
            app.print(Panel(
                "\n".join(lines),
                title="Box–Bilayer Match",
                border_style="yellow" if (warns or errors) else "green",
                padding=(0, 2),
            ))
    except Exception:
        pass


# ═══════════════════════════════════════════════════════════════════════════════
# init
# ═══════════════════════════════════════════════════════════════════════════════

def _ask(prompt: str, options: list[tuple[str, str]], default: int = 1) -> str:
    """Numbered menu prompt. Returns the key of the selected option."""
    app.print(f"\n  [bold]{prompt}[/bold]")
    for i, (key, label) in enumerate(options, 1):
        marker = "[green]●[/green]" if i == default else " "
        app.print(f"    [{i}] {marker} {label}")
    while True:
        raw = input(f"\n  Select [1-{len(options)}] (Enter = {default}): ").strip()
        if raw == "":
            return options[default - 1][0]
        if raw.isdigit() and 1 <= int(raw) <= len(options):
            return options[int(raw) - 1][0]
        app.print("  [red]Opción inválida.[/red]")


def _ask_multi(prompt: str, options: list[tuple[str, str]], defaults: list[int]) -> list[str]:
    """Multi-select numbered menu. Returns list of selected keys."""
    app.print(f"\n  [bold]{prompt}[/bold]")
    app.print("  [dim](separados por coma, ej: 1,2,3)[/dim]")
    for i, (key, label) in enumerate(options, 1):
        marker = "[green]●[/green]" if i in defaults else " "
        app.print(f"    [{i}] {marker} {label}")
    default_str = ",".join(str(d) for d in defaults)
    while True:
        raw = input(f"\n  Select (Enter = {default_str}): ").strip()
        if raw == "":
            return [options[d - 1][0] for d in defaults]
        parts = [p.strip() for p in raw.split(",")]
        if all(p.isdigit() and 1 <= int(p) <= len(options) for p in parts):
            return [options[int(p) - 1][0] for p in parts]
        app.print("  [red]Entrada inválida.[/red]")


def _ask_float(prompt: str, default: float) -> float:
    """Free-text float prompt with default."""
    raw = input(f"\n  {prompt} (Enter = {default}): ").strip()
    if raw == "":
        return default
    try:
        return float(raw)
    except ValueError:
        app.print("  [yellow]Valor inválido, usando default.[/yellow]")
        return default


def _ask_str(prompt: str, default: str = "") -> str:
    """Free-text string prompt."""
    suffix = f" (Enter = {default})" if default else ""
    raw = input(f"\n  {prompt}{suffix}: ").strip()
    return raw if raw else default


@cli.command()
def init(
    output: Optional[Path] = typer.Option(None, "--output", "-o", help="Output YAML path (default: configs/<name>.yaml)."),
):
    """Interactive wizard — genera un YAML de configuración para SimForge."""

    app.print(Panel(
        "[bold cyan]SimForge Init[/bold cyan]  [dim]generador de configuración[/dim]",
        border_style="cyan",
        padding=(0, 2),
    ))

    # ── 1. Nombre del proyecto ────────────────────────────────────────────────
    name = _ask_str("Nombre del proyecto")
    if not name:
        app.print("[red]El nombre es obligatorio.[/red]")
        raise typer.Exit(1)
    name = name.replace(" ", "_")

    # ── 2. Tipo de sistema ────────────────────────────────────────────────────
    system_type = _ask("Tipo de sistema", [
        ("protein_in_water",  "Proteína / péptido en agua"),
        ("protein_ligand",    "Proteína + ligando (inhibición, binding)"),
        ("protein_membrane",  "Proteína + membrana"),
    ])

    # ── 2b. Configuración de membrana ─────────────────────────────────────────
    lipid_type: str | None = None
    if system_type == "protein_membrane":
        lipid_type = _ask("Tipo de lípido para la bicapa", [
            ("DPPC", "DPPC  (dipalmitoylphosphatidylcholine — más común, OPLS-AA)"),
            ("POPC", "POPC  (1-palmitoyl-2-oleoyl-sn-glycero-3-phosphocholine)"),
            ("POPE", "POPE  (phosphatidylethanolamine)"),
            ("DMPC", "DMPC  (dimyristoylphosphatidylcholine — bicapas finas)"),
        ])

    # ── 2c. Structural annotation (protein_membrane only) ─────────────────────
    _structural_ann = None
    if system_type == "protein_membrane":
        want_ann = _ask("¿Añadir anotación EC/IC/TM ahora?", [
            ("yes", "Sí — definir topología de membrana (habilita orient_protein AUTOMATED)"),
            ("no",  "No — añadir después con 'simforge annotate-structure <config>'"),
        ], default=2)
        if want_ann == "yes":
            from core.structural_annotation import (
                StructuralAnnotation,
                MembraneTopologyAnnotation,
                OrientationAnnotation,
                OrientationEvidence,
            )
            app.print(Panel(
                "[bold cyan]Structural Biology Annotation[/bold cyan]  "
                "[dim]EC · IC · TM topology[/dim]",
                border_style="cyan",
                padding=(0, 2),
            ))

            app.print("\n[bold]── Regiones extracelulares ──────────────────────────────────────────[/bold]")
            _ec = _ask_range_list("Residuos en la cara extracelular", [])

            app.print("\n[bold]── Regiones intracelulares ──────────────────────────────────────────[/bold]")
            if _ec:
                app.print(
                    "  [dim]Necesario para orientación automática: sin IC no se puede calcular el eje EC→IC.[/dim]"
                )
            _ic = _ask_range_list("Residuos en la cara intracelular", [])

            if _ec and not _ic:
                app.print(
                    "\n  [yellow]⚠[/yellow]  [yellow]orient_protein quedará GUIDED:[/yellow] "
                    "hay regiones EC pero faltan regiones IC.\n"
                    "     Define al menos una región IC para habilitar orientación automática."
                )

            app.print("\n[bold]── Segmentos transmembrana ──────────────────────────────────────────[/bold]")
            _tm = _ask_tm_segments([])

            app.print("\n[bold]── Orientación geométrica ───────────────────────────────────────────[/bold]")
            app.print(
                "  [dim]¿En qué dirección del eje Z queda la cara extracelular "
                "después de editconf -princ?[/dim]"
            )
            _ec_side = _ask("Cara extracelular", [
                ("+z", "+Z  (convención GROMACS estándar)"),
                ("-z", "-Z  (proteína invertida respecto a la convención)"),
            ])
            _ic_side = "-z" if _ec_side == "+z" else "+z"

            app.print("\n[bold]── Fuente y confianza ───────────────────────────────────────────────[/bold]")
            _source = _ask("Fuente de la anotación", [
                ("user_annotation", "Anotación manual (tú lo sabes)"),
                ("opm_database",    "OPM database (Orientations of Proteins in Membranes)"),
                ("pdbtm",           "PDBTM (PDB Transmembrane database)"),
                ("uniprot",         "UniProt (feature annotations)"),
                ("predicted",       "Predicción computacional (TMH, DeepTMHMM, etc.)"),
            ])
            _conf = max(0.0, min(1.0, _ask_float("Confianza [0.0–1.0]", 0.9)))
            _ref  = _ask_str("Referencia (OPM accession, PubMed ID, URL — opcional)")

            _topology = MembraneTopologyAnnotation(
                extracellular_regions=_ec,
                intracellular_regions=_ic,
                transmembrane_segments=_tm,
            )
            _orientation = OrientationAnnotation(
                extracellular_side=_ec_side,
                intracellular_side=_ic_side,
                source=_source,
                confidence=_conf,
                reference=_ref or None,
            ) if (_ec or _ic) else None
            _evidence_note = (
                "Anotación creada con simforge init"
                + (f"; referencia: {_ref}" if _ref else "")
            )
            _structural_ann = StructuralAnnotation(
                membrane_topology=_topology if (_ec or _ic or _tm) else None,
                orientation=_orientation,
                evidence=[OrientationEvidence(
                    source=_source, confidence=_conf,
                    reference=_ref or None, notes=_evidence_note,
                )] if (_ec or _ic) else [],
            )

            # Show validation and resulting automation level
            _ann_warns = _structural_ann.validation_warnings()
            if _ann_warns:
                app.print("\n[yellow]Advertencias:[/yellow]")
                for w in _ann_warns:
                    app.print(f"  [yellow]⚠[/yellow]  {w}")

            _overlap_warns = _topology.overlap_warnings() if (_ec or _ic or _tm) else []
            if _overlap_warns:
                app.print("\n[red]Solapamiento de regiones:[/red]")
                for w in _overlap_warns:
                    app.print(f"  [red]✗[/red]  {w}")

            app.print()
            if _structural_ann.is_complete_for_orient() or _structural_ann.is_partial_for_orient():
                app.print(
                    "  [green]● orient_protein → AUTOMATED[/green]  "
                    "[dim](EC + IC + orientación definidas)[/dim]"
                )
            else:
                app.print(
                    "  [red]● orient_protein → GUIDED[/red]  "
                    "[dim](topología incompleta — añade IC/EC para AUTOMATED)[/dim]"
                )

    # ── 3. Componentes ────────────────────────────────────────────────────────
    app.print("\n  [bold]Archivos PDB[/bold]")
    protein_file = _ask_str("  Archivo PDB de la proteína/péptido")
    if not protein_file:
        app.print("[red]El PDB de proteína es obligatorio.[/red]")
        raise typer.Exit(1)

    extra_components: list[dict] = []

    if system_type == "protein_ligand":
        n_extra = int(_ask("¿Cuántos ligandos/sustratos?", [
            ("1", "1"), ("2", "2"), ("3", "3"),
        ]))
        for i in range(1, n_extra + 1):
            role = _ask(f"Rol del componente {i}", [
                ("competitive_ligand", "Ligando competitivo"),
                ("substrate",          "Sustrato"),
                ("ligand",             "Ligando (genérico)"),
                ("cofactor",           "Cofactor"),
            ])
            pdb = _ask_str(f"  Archivo PDB del {role}_{i}")
            extra_components.append({"role": role, "file": pdb, "index": i})

    # ── 4. Forcefield ─────────────────────────────────────────────────────────
    ff_protein = _ask("Forcefield (proteína)", [
        ("opls-aa",   "OPLS-AA  (péptidos, proteínas pequeñas)"),
        ("charmm36",  "CHARMM36 (proteínas, membrana)"),
        ("amber99sb", "AMBER99SB-ILDN"),
    ])

    ff_ligands = None
    if system_type == "protein_ligand":
        ff_ligands = _ask("Forcefield (ligandos)", [
            ("gaff2",  "GAFF2   (AMBER ecosystem)"),
            ("cgenff", "CGenFF  (CHARMM ecosystem)"),
            ("opls",   "OPLS-AA (ligandos simples)"),
        ])

    # ── 5. Modelo de agua ─────────────────────────────────────────────────────
    water_model = _ask("Modelo de agua", [
        ("spce",   "SPC/E   (recomendado para OPLS-AA)"),
        ("tip3p",  "TIP3P   (recomendado para CHARMM/AMBER)"),
        ("tip4p",  "TIP4P-Ew"),
    ])

    # ── 6. Temperatura ────────────────────────────────────────────────────────
    temp_choice = _ask("Temperatura", [
        ("300",    "300 K  (estándar)"),
        ("309.65", "309.65 K  (37 °C — fisiológico humano)"),
        ("custom", "Otra..."),
    ])
    if temp_choice == "custom":
        temperature = _ask_float("Temperatura (K)", 300.0)
    else:
        temperature = float(temp_choice)

    # ── 7. Duración ───────────────────────────────────────────────────────────
    dur_choice = _ask("Duración de la simulación de producción", [
        ("1",    " 1 ns  (exploración rápida)"),
        ("10",   "10 ns  (estabilidad)"),
        ("50",   "50 ns  (binding, inhibición)"),
        ("100",  "100 ns"),
        ("custom", "Otra..."),
    ])
    if dur_choice == "custom":
        duration = _ask_float("Duración (ns)", 10.0)
    else:
        duration = float(dur_choice)

    # ── 8. Análisis ───────────────────────────────────────────────────────────
    all_analyses = [
        ("rmsd",             "RMSD  (estabilidad global)"),
        ("rmsf",             "RMSF  (flexibilidad por residuo)"),
        ("energy",           "Energía  (Ep, Ek, T, P)"),
        ("hydrogen_bonds",   "Puentes de hidrógeno"),
        ("radius_of_gyration", "Radio de giro"),
        ("distance_analysis", "Distancias entre grupos"),
    ]
    default_analyses = [1, 2, 3]
    if system_type == "protein_ligand":
        default_analyses = [1, 2, 3, 4, 6]

    selected_analyses = _ask_multi("Análisis a incluir", all_analyses, default_analyses)

    # ── 9. Hardware ───────────────────────────────────────────────────────────
    hardware = _ask("Hardware disponible para mdrun", [
        ("auto", "Auto-detectar en runtime (recomendado)"),
        ("gpu",  "GPU  (NVIDIA, siempre usar flags GPU)"),
        ("cpu",  "CPU  (sin GPU)"),
    ])

    # ── 10. Objetivos de simulación ───────────────────────────────────────────
    objectives_map = {
        "protein_in_water": ["stability"],
        "protein_ligand":   ["competitive_binding", "active_site_stability"],
        "protein_membrane": ["stability", "membrane_protein_dynamics"],
    }
    objectives = objectives_map[system_type]

    # ── Construir YAML ────────────────────────────────────────────────────────
    components_yaml = f"  - id: protein_1\n    role: protein\n    file: {protein_file}\n"
    for ec in extra_components:
        role_id = f"{ec['role']}_{ec['index']}"
        components_yaml += (
            f"\n  - id: {role_id}\n"
            f"    role: {ec['role']}\n"
            f"    file: {ec['file']}\n"
        )

    forcefields_yaml = f"  protein: {ff_protein}\n"
    if ff_ligands:
        forcefields_yaml += f"  ligands: {ff_ligands}\n"

    objectives_yaml = "\n".join(f"  - {o}" for o in objectives)

    analyses_yaml = ""
    for a in selected_analyses:
        analyses_yaml += f"\n  - type: {a}"
        if a == "distance_analysis" and system_type == "protein_ligand" and extra_components:
            first = extra_components[0]
            second = extra_components[1] if len(extra_components) > 1 else extra_components[0]
            g1 = f"{first['role']}_{first['index']}"
            g2 = f"{second['role']}_{second['index']}"
            analyses_yaml += f"\n    selection:\n      group1: {g1}\n      group2: {g2}"

    hardware_line = f"\nhardware: {hardware}\n" if hardware != "auto" else ""

    membrane_section = ""
    if system_type == "protein_membrane" and lipid_type:
        membrane_section = f"""  membrane:
    enabled: true
    type: {lipid_type}
"""

    yaml_content = f"""project:
  name: {name}

components:
{components_yaml}
environment:
{membrane_section}  solvent:
    water_model: {water_model}
  ions:
    concentration: 0.154
  temperature_K: {temperature}
  duration_ns: {duration}

forcefields:
{forcefields_yaml}
simulation_objectives:
{objectives_yaml}

analysis:{analyses_yaml}
{hardware_line}"""

    # ── Escribir archivo ──────────────────────────────────────────────────────
    if output is None:
        configs_dir = Path("configs")
        configs_dir.mkdir(exist_ok=True)
        output = configs_dir / f"{name}.yaml"

    output.write_text(yaml_content)

    # Inject structural_annotation if collected during the wizard
    if _structural_ann is not None:
        from ruamel.yaml import YAML as RuamelYAML
        _ryaml = RuamelYAML()
        _ryaml.preserve_quotes   = True
        _ryaml.default_flow_style = False
        _ryaml.best_sequence_indent = 2
        _ryaml.best_map_flow_style  = False
        _ryaml.width = 4096
        with open(output) as _f:
            _doc = _ryaml.load(_f)
        _doc["structural_annotation"] = _annotation_to_yaml_dict(_structural_ann)
        with open(output, "w") as _f:
            _ryaml.dump(_doc, _f)
        app.print(f"  [green]✓[/green] structural_annotation incluida")
        if _structural_ann.membrane_topology:
            _mt = _structural_ann.membrane_topology
            app.print(
                f"    EC: {_mt.extracellular_regions}  "
                f"IC: {_mt.intracellular_regions}  "
                f"TM: {_mt.tm_segment_count()} segmentos"
            )

    app.print(f"\n[bold green]✓[/bold green] Config generado → [cyan]{output}[/cyan]")
    app.print(f"\nSiguiente:  [dim]simforge compile {output}[/dim]\n")


# ═══════════════════════════════════════════════════════════════════════════════
# Scientific planning dialogue
# ═══════════════════════════════════════════════════════════════════════════════

_KIND_LABELS = {
    "missing_parameter":   ("?",  "dim"),
    "policy_selection":    ("◆",  "cyan"),
    "risk_confirmation":   ("⚠",  "yellow"),
    "scientific_tradeoff": ("⟁",  "magenta"),
    "protocol_selection":  ("◎",  "blue"),
}

_KIND_TITLES = {
    "missing_parameter":   "Missing Parameter",
    "policy_selection":    "Policy Selection",
    "risk_confirmation":   "Risk Confirmation",
    "scientific_tradeoff": "Scientific Tradeoff",
    "protocol_selection":  "Protocol Selection",
}


def _run_planning_dialogue(
    state,
    yaml_path: Path,
    skip: bool = False,
) -> Optional[object]:
    """
    Detecta preguntas de planning, muestra el diálogo y aplica patches al state.

    Retorna PlanningSessionRecord, o None si el usuario eligió abortar.
    """
    from core.scientific_planner import detect_questions, apply_patches
    from core.planning_models import PlanningSessionRecord, PlanningAnswer
    from datetime import datetime

    questions = detect_questions(state)

    record = PlanningSessionRecord(
        created_at  = datetime.now().isoformat(timespec="seconds"),
        config_path = str(yaml_path),
        skipped     = skip or len(questions) == 0,
    )

    if not questions or skip:
        if questions and skip:
            app.print(f"[dim]--no-plan: skipping {len(questions)} planning question(s)[/dim]\n")
        return record

    app.print(f"[bold]Scientific Planning[/bold]  [dim]{len(questions)} question(s) detected[/dim]\n")

    answers: list[PlanningAnswer] = []
    total = len(questions)

    for idx, q in enumerate(questions, 1):
        icon, color = _KIND_LABELS.get(q.kind.value, ("?", "dim"))
        title = _KIND_TITLES.get(q.kind.value, q.kind.value)

        # ── Header ───────────────────────────────────────────────────────────
        app.rule(f"[{color}]{icon} [{idx}/{total}] {title}[/{color}]")
        app.print()

        # ── Context ───────────────────────────────────────────────────────────
        for line in q.context.splitlines():
            app.print(f"  {line}")
        app.print()

        # ── Question ──────────────────────────────────────────────────────────
        app.print(f"  [bold]{q.question}[/bold]")
        app.print()

        # ── Options ───────────────────────────────────────────────────────────
        for i, opt in enumerate(q.options, 1):
            marker = "[bold green]●[/bold green]" if opt.is_default else " "
            app.print(f"  [{i}] {marker} [white]{opt.label}[/white]")
            app.print(f"       [dim]{opt.description}[/dim]")
        app.print()

        # ── Prompt ────────────────────────────────────────────────────────────
        default_idx = next(
            (i for i, o in enumerate(q.options, 1) if o.is_default), 1
        )

        while True:
            raw = input(f"  Select [1-{len(q.options)}] (Enter = {default_idx}): ").strip()
            if raw == "":
                choice = default_idx
                break
            if raw.isdigit() and 1 <= int(raw) <= len(q.options):
                choice = int(raw)
                break
            app.print(f"  [red]Invalid choice.[/red] Enter a number between 1 and {len(q.options)}.")

        selected = q.options[choice - 1]
        app.print(f"\n  [green]→[/green] {selected.label}\n")

        # ── Abort check ───────────────────────────────────────────────────────
        if selected.state_patch.get("abort"):
            app.print("[yellow]Compile aborted by planning decision.[/yellow]")
            app.print(f"  Reason: {selected.description}\n")
            return None

        answer = PlanningAnswer(
            question_id    = q.id,
            selected_key   = selected.key,
            selected_label = selected.label,
            state_patch    = {k: v for k, v in selected.state_patch.items() if k != "abort"},
        )
        answers.append(answer)

    # ── Apply patches ─────────────────────────────────────────────────────────
    all_patches: list[dict] = [a.state_patch for a in answers if a.state_patch]
    apply_patches(state, answers)

    # ── Show summary ──────────────────────────────────────────────────────────
    if all_patches:
        app.rule("[dim]Planning complete[/dim]")
        app.print()
        for a in answers:
            if a.state_patch:
                app.print(f"  [green]✓[/green] {a.selected_label}")
        app.print()

    # ── Build record ──────────────────────────────────────────────────────────
    record.questions_shown = [q.model_dump(mode="json") for q in questions]
    record.answers         = answers
    record.patches_applied = all_patches

    return record


def _save_planning_session(record, path: Path) -> None:
    """Persiste la sesión de planning en metadata/ del workspace."""
    if record is None:
        return
    try:
        path.write_text(
            json.dumps(record.model_dump(mode="json"), indent=4)
        )
    except Exception:
        pass   # no bloquear compile por fallo de persistencia


# ═══════════════════════════════════════════════════════════════════════════════
# summary
# ═══════════════════════════════════════════════════════════════════════════════

@cli.command()
def summary(
    workspace: Path = typer.Argument(Path("."), help="Workspace directory to summarize (default: current directory)."),
    json_out:  bool = typer.Option(False, "--json", help="Output raw JSON instead of rich text."),
):
    """Show (or generate) the scientific summary for a workspace."""
    from runtime.scientific_summary import generate_summary

    if not workspace.exists():
        app.print(f"[red]Error:[/red] Workspace not found: {workspace}")
        raise typer.Exit(1)

    # Try reading cached summary first; regenerate if absent
    cached_json = workspace / "metadata" / "scientific_summary.json"
    if cached_json.exists():
        try:
            import json as _json
            data = _json.loads(cached_json.read_text())
            if json_out:
                import sys
                _json.dump(data, sys.stdout, indent=2)
                app.print()
                return
        except Exception:
            pass  # fall through to regenerate

    # Regenerate
    sm = generate_summary(workspace)

    if json_out:
        import json as _json
        import sys
        _json.dump(sm.as_dict(), sys.stdout, indent=2)
        app.print()
        return

    # Rich display
    converged_str = "[green]Yes[/green]" if sm.converged else "[red]No[/red]"
    runtime_str   = f"{sm.runtime_ns:.3f} ns" if sm.runtime_ns else "unknown"

    app.print(Panel(
        f"  Workspace  [cyan]{workspace}[/cyan]\n"
        f"  Converged  {converged_str}\n"
        f"  Runtime    {runtime_str}\n"
        f"  Analyses   {len(sm.analyses)}",
        title="Scientific Summary",
        border_style="cyan",
    ))

    # Orientation section (membrane workflows only)
    try:
        from runtime.orientation_gate import read_orientation_report
        orient = read_orientation_report(workspace)
        if orient is not None:
            conf   = orient.get("confidence", 0.0)
            passed = orient.get("passed", False)
            warns  = orient.get("warnings", [])
            errors = orient.get("errors",   [])
            geom   = orient.get("geometry", {})
            status_color = "green" if passed else "red"
            status_label = "yes" if passed else "no"
            orient_lines = [
                f"  confidence: {conf:.2f}   "
                f"passed: [{status_color}]{status_label}[/{status_color}]",
            ]
            if geom.get("tilt_angle_deg") is not None:
                orient_lines.append(f"  tilt: {geom['tilt_angle_deg']:.1f}°   "
                    f"EC COM z: {geom.get('ec_com_z_nm', 'N/A')} nm   "
                    f"IC COM z: {geom.get('ic_com_z_nm', 'N/A')} nm")
            for w in warns:
                orient_lines.append(f"  [yellow]⚠[/yellow]  {w}")
            for e in errors:
                orient_lines.append(f"  [red]✗[/red]  {e}")
            app.print(Panel(
                "\n".join(orient_lines),
                title="Orientation",
                border_style="yellow" if (warns or errors) else "green",
                padding=(0, 2),
            ))
    except Exception:
        pass

    # Box-match section (membrane workflows only)
    try:
        from runtime.box_match_gate import read_box_match_report
        bm = read_box_match_report(workspace)
        if bm is not None:
            rb     = bm.get("recommended_box", {})
            pg     = bm.get("protein_geometry", {})
            est    = bm.get("estimates", {})
            passed = bm.get("passed", False)
            conf   = bm.get("confidence", 0.0)
            warns  = bm.get("warnings", [])
            errors = bm.get("errors",   [])
            sc     = "green" if passed else "red"
            sl     = "yes" if passed else "no"
            bm_lines = [
                f"  confidence: {conf:.2f}   passed: [{sc}]{sl}[/{sc}]"
                f"   lipid: {bm.get('lipid_type', '?')}",
            ]
            if rb:
                bm_lines.append(
                    f"  recommended box: {rb.get('box_x_nm', '?'):.3f} × "
                    f"{rb.get('box_y_nm', '?'):.3f} × "
                    f"{rb.get('box_z_nm', '?'):.3f} nm"
                )
            if pg:
                bm_lines.append(
                    f"  protein footprint: {pg.get('x_extent_nm', '?'):.2f} × "
                    f"{pg.get('y_extent_nm', '?'):.2f} nm   "
                    f"height: {pg.get('z_extent_nm', '?'):.2f} nm"
                )
            if est.get("n_lipids_estimate"):
                cov = int(est.get("protein_xy_coverage", 0) * 100)
                bm_lines.append(
                    f"  N lipids estimate: {est['n_lipids_estimate']}   "
                    f"XY coverage: {cov}%   "
                    f"solvent: {est.get('solvent_volume_nm3', '?'):.0f} nm³"
                )
            for w in warns:
                bm_lines.append(f"  [yellow]⚠[/yellow]  {w}")
            for e in errors:
                bm_lines.append(f"  [red]✗[/red]  {e}")
            app.print(Panel(
                "\n".join(bm_lines),
                title="Box–Bilayer Match",
                border_style="yellow" if (warns or errors) else "green",
                padding=(0, 2),
            ))
    except Exception:
        pass

    if sm.rmsd_verdict:
        app.print(Panel(sm.rmsd_verdict,   title="RMSD Convergence",  border_style="dim"))
    if sm.energy_verdict:
        app.print(Panel(sm.energy_verdict, title="Energy Stability",   border_style="dim"))

    if sm.warnings:
        for w in sm.warnings:
            app.print(f"  [yellow]⚠[/yellow]  {w}")


# ═══════════════════════════════════════════════════════════════════════════════
# analyze
# ═══════════════════════════════════════════════════════════════════════════════

_QUALITY_STYLE = {
    "converged":            ("[green]",  "✓", "green"),
    "partially_converged":  ("[yellow]", "~", "yellow"),
    "not_converged":        ("[red]",    "✗", "red"),
    "problematic":          ("[red]",    "⚠", "red"),
    "insufficient_data":    ("[dim]",    "?", "dim"),
}


@cli.command()
def analyze(
    path:    str           = typer.Argument(".", help="Path to simulation directory or XVG files (default: current directory)"),
    output:  Optional[str] = typer.Option(None,        "--output",  "-o", help="Output file (default: stdout)"),
    format:  str           = typer.Option("markdown",  "--format",  "-f", help="Output format: markdown|json"),
    context: Optional[str] = typer.Option(None,        "--context", "-c",
        help="System context for context-aware interpretation: globular_protein | membrane_protein | "
             "idr | protein_ligand_complex | multimeric_complex | enzyme | peptide | "
             "membrane_system | flexible_domain_protein"),
):
    """Analyze an existing MD simulation and classify its scientific quality."""
    from runtime.scientific_summary import analyze_trajectory
    import json as _json

    sim_path = Path(path)
    if not sim_path.exists():
        app.print(f"[red]Error:[/red] Path not found: {sim_path}")
        raise typer.Exit(1)

    ctx_label = f" [{context}]" if context else ""
    with app.status(f"  Analyzing [cyan]{sim_path}[/cyan]{ctx_label}..."):
        try:
            from runtime.trajectory_ingestor import discover_trajectory
            manifest = discover_trajectory(sim_path)
            _summary, report = analyze_trajectory(sim_path, context=context)
        except Exception as exc:
            app.print(f"[red]Analysis error:[/red] {exc}")
            raise typer.Exit(1)

    if context:
        app.print(f"[dim]  Context: {context}[/dim]")

    # ── JSON output ──────────────────────────────────────────────────────────
    if format == "json":
        out = report.as_dict()
        out["discovered_files"] = {lbl: str(p.name) for lbl, p in manifest.xvg_files.items()}
        out_text = _json.dumps(out, indent=2)
        if output:
            Path(output).write_text(out_text)
            app.print(f"[dim]Written to {output}[/dim]")
        else:
            import sys
            print(out_text, file=sys.stdout)
        return

    # ── Markdown / Rich output ────────────────────────────────────────────────
    q_val = report.quality.value
    color_tag, icon, border = _QUALITY_STYLE.get(q_val, ("[white]", "?", "white"))

    badge = f"{color_tag}{icon} {q_val.replace('_', ' ').upper()}[/{color_tag[1:]}"
    confidence_bar_filled = max(0, int(report.confidence * 20))
    confidence_bar = "█" * confidence_bar_filled + "░" * (20 - confidence_bar_filled)
    conf_pct = f"{report.confidence:.0%}"

    header_lines = [
        f"  Path         [cyan]{sim_path}[/cyan]",
        f"  Quality      {badge}",
        f"  Confidence   [{confidence_bar}] {conf_pct}",
    ]
    if report.metrics:
        if "rmsd_total_ns" in report.metrics:
            header_lines.append(f"  RMSD data    {report.metrics['rmsd_total_ns']} ns")
        if "energy_total_ns" in report.metrics:
            header_lines.append(f"  Energy data  {report.metrics['energy_total_ns']} ns")

    # ── Discovered XVG files ─────────────────────────────────────────────────
    if manifest.xvg_files:
        from rich.table import Table as _Table
        from runtime.trajectory_ingestor import load_xvg_files as _load_xvg
        _xvg_map = _load_xvg(manifest)
        ftable = _Table(box=box.SIMPLE, show_header=False, padding=(0, 1))
        ftable.add_column("label", style="cyan", no_wrap=True)
        ftable.add_column("arrow", width=2, style="dim")
        ftable.add_column("file",  style="dim")
        ftable.add_column("unit",  style="dim")
        for lbl, fpath in sorted(manifest.xvg_files.items()):
            _d = _xvg_map.get(lbl)
            if _d and _d.x_unit:
                _unit_tag = f"[time: {_d.x_unit}]"
            elif _d and _d.time_ps:
                _unit_tag = "[time: ps (assumed)]"
            else:
                _unit_tag = ""
            ftable.add_row(lbl, "←", fpath.name, _unit_tag)
        app.print(Panel(ftable, title="Discovered XVG files", border_style="dim", padding=(0, 1)))
    else:
        app.print(Panel(
            "  [dim]No XVG files found in this directory.[/dim]",
            title="Discovered XVG files", border_style="dim", padding=(0, 1),
        ))

    app.print(Panel(
        "\n".join(header_lines),
        title="MD Simulation Quality",
        border_style=border,
        padding=(0, 2),
    ))

    # Evidence
    if report.evidence:
        ev_text = "\n".join(f"  [dim]•[/dim] {e}" for e in report.evidence)
        app.print(Panel(ev_text, title="Evidence", border_style="dim", padding=(0, 1)))

    # Warnings
    if report.warnings:
        w_text = "\n".join(f"  [yellow]⚠[/yellow] {w}" for w in report.warnings)
        app.print(Panel(w_text, title="Warnings", border_style="yellow", padding=(0, 1)))

    # Recommendations
    if report.recommendations:
        r_text = "\n".join(f"  [cyan]→[/cyan] {r}" for r in report.recommendations)
        app.print(Panel(r_text, title="Recommendations", border_style="cyan", padding=(0, 1)))

    # Key metrics table
    if report.metrics:
        from rich.table import Table as _Table
        mt = _Table(box=box.SIMPLE, show_header=True, header_style="bold dim")
        mt.add_column("Metric", style="white")
        mt.add_column("Value", justify="right")
        for k, v in report.metrics.items():
            mt.add_row(k.replace("_", " "), str(v))
        app.print(Panel(mt, title="Key Metrics", border_style="dim", padding=(0, 1)))

    # Optionally write to file (markdown)
    if output:
        Path(output).write_text(report.as_markdown())
        app.print(f"\n[dim]Report written to {output}[/dim]")


# ═══════════════════════════════════════════════════════════════════════════════
# study
# ═══════════════════════════════════════════════════════════════════════════════

_FIND_LEVEL_COLOR = {"highlight": "cyan",   "info": "dim",  "warning": "yellow"}
_FIND_LEVEL_ICON  = {"highlight": "→",      "info": "·",    "warning": "⚠"}

_STATE_STYLE: dict[str, tuple[str, str]] = {
    "stable_binding":               ("green",    "●"),
    "interaction_persistent":       ("green",    "●"),
    "structurally_stable":          ("green",    "●"),
    "weak_binding":                 ("yellow",   "◐"),
    "flexible_but_stable":          ("yellow",   "◐"),
    "transient_binding":            ("yellow",   "◑"),
    "ligand_destabilization":       ("red",      "○"),
    "conformational_rearrangement": ("red",      "○"),
    "possible_dissociation":        ("red",      "✗"),
    "uncertain_behavior":           ("dim",      "?"),
}


def _generate_markdown_report(study_obj, summary_obj, synthesis, path: Path) -> str:
    from datetime import date as _date

    lines: list[str] = []

    # ── Header ────────────────────────────────────────────────────────────────
    lines += [
        "# SimForge Study Report",
        f"**Path:** `{path.resolve()}`  ",
        f"**Date:** {_date.today().isoformat()}  ",
        f"**Files:** {study_obj.n_xvg_parsed} XVG parsed / {study_obj.n_xvg_discovered} discovered",
        "",
    ]

    # ── Systems overview ──────────────────────────────────────────────────────
    lines += ["---", "", "## Systems Overview", ""]
    lines += ["| System | Replicas | Observables |", "|--------|----------|-------------|"]
    for sys_name, sg in sorted(study_obj.systems.items()):
        obs_str = ", ".join(study_obj.observable_display.get(o, o) for o in sg.observables)
        lines.append(f"| {sys_name} | {sg.n_replicas} | {obs_str} |")
    lines.append("")

    # ── Observable statistics ──────────────────────────────────────────────────
    if study_obj.observables_detected:
        sys_names = sorted(study_obj.systems.keys())
        lines += ["## Observable Statistics", ""]
        header_cols = " | ".join(f"{s} mean ± std" for s in sys_names)
        sep_cols    = "|".join("---" for _ in sys_names)
        lines += [f"| Observable | {header_cols} |", f"|------------|{sep_cols}|"]
        for obs in study_obj.observables_detected:
            disp  = study_obj.observable_display.get(obs, obs)
            units = study_obj.observable_units.get(obs, "")
            row   = f"| {disp}"
            for sys_name in sys_names:
                sg = study_obj.systems.get(sys_name)
                if sg and obs in sg.aggregate:
                    agg = sg.aggregate[obs]
                    u   = f" {units}" if units else ""
                    row += f" | {agg.mean:.3g} ± {agg.std:.3g}{u}"
                else:
                    row += " | —"
            lines.append(row + " |")
        lines.append("")

    # ── Scientific synthesis ───────────────────────────────────────────────────
    if synthesis:
        lines += ["## Scientific Synthesis", ""]

        # Interaction state classification
        if synthesis.systems:
            lines += ["### Interaction State Classification", ""]
            lines += [
                "| System | State | Conf. | Binding | Stability | Evidence |",
                "|--------|-------|-------|---------|-----------|----------|",
            ]
            for sys_name, syn in sorted(synthesis.systems.items()):
                _, icon = _STATE_STYLE.get(syn.primary_state, ("dim", "?"))
                state_label = f"{icon} {syn.primary_state.replace('_', ' ')}"
                ev_str = "; ".join(syn.evidence[:2]) if syn.evidence else "—"
                lines.append(
                    f"| {sys_name} | {state_label} | {syn.primary_confidence:.2f}"
                    f" | {syn.binding_score:.2f} | {syn.stability_score:.2f} | {ev_str} |"
                )
            lines.append("")

        # Composite ranking
        if synthesis.ranking and len(synthesis.ranking) > 1:
            lines += ["### Composite Ranking", ""]
            lines += ["| Rank | System | Composite | Binding | Stability |",
                      "|------|--------|-----------|---------|-----------|"]
            for i, (sys_name, score) in enumerate(synthesis.ranking, 1):
                syn = synthesis.systems.get(sys_name)
                b = f"{syn.binding_score:.3f}"   if syn else "—"
                s = f"{syn.stability_score:.3f}" if syn else "—"
                lines.append(f"| {i} | {sys_name} | {score:.3f} | {b} | {s} |")
            lines.append("")

        # Temporal events
        if synthesis.events:
            lines += ["### Temporal Events", ""]
            lines += [
                "| System | Replica | Observable | Event | Time (ns) |",
                "|--------|---------|------------|-------|-----------|",
            ]
            for evt in synthesis.events:
                obs_d = study_obj.observable_display.get(evt.observable, evt.observable)
                etype = evt.event_type.replace("_", " ")
                lines.append(
                    f"| {evt.system} | {evt.replica} | {obs_d}"
                    f" | {etype} | {evt.time_ns:.1f} |"
                )
            lines += ["", "> **Descriptions**", ""]
            for evt in synthesis.events:
                lines.append(f"- {evt.description}")
            lines.append("")

        # Narrative
        if synthesis.narrative:
            lines += ["### Scientific Narrative", "", synthesis.narrative, ""]

    # ── Comparative findings ───────────────────────────────────────────────────
    if summary_obj and summary_obj.findings:
        sorted_findings = sorted(
            summary_obj.findings,
            key=lambda f: (0 if f.level != "info" else 1, f.level),
        )
        lines += ["## Comparative Findings", ""]
        _ICON = {"highlight": "→", "info": "·", "warning": "⚠"}
        for f in sorted_findings[:12]:
            icon = _ICON.get(f.level, "·")
            lines.append(f"{icon} {f.message}")
        lines.append("")

    # ── Outlier replicas ───────────────────────────────────────────────────────
    if summary_obj and summary_obj.outlier_replicas:
        lines += ["## Outlier Replicas", ""]
        lines += ["| System | Replica | Reason |", "|--------|---------|--------|"]
        for sys_n, rep_n, reason in summary_obj.outlier_replicas:
            lines.append(f"| {sys_n} | {rep_n} | {reason} |")
        lines.append("")

    lines += ["---", "*Generated by [SimForge](https://github.com/simforge)*", ""]
    return "\n".join(lines)


def _show_interaction_states(synthesis, study_obj) -> None:
    table = Table(box=box.SIMPLE, show_header=True, header_style="bold dim")
    table.add_column("System",    style="bold cyan")
    table.add_column("State",     min_width=26)
    table.add_column("Conf",      justify="right", style="dim")
    table.add_column("Binding",   justify="right", style="dim")
    table.add_column("Stability", justify="right", style="dim")
    table.add_column("Evidence",  overflow="fold")

    for sys_name, syn in sorted(synthesis.systems.items()):
        color, icon = _STATE_STYLE.get(syn.primary_state, ("dim", "?"))
        state_label = syn.primary_state.replace("_", " ")
        ev_str      = "; ".join(syn.evidence[:2]) if syn.evidence else "—"
        table.add_row(
            sys_name,
            f"[{color}]{icon} {state_label}[/{color}]",
            f"{syn.primary_confidence:.0%}",
            f"{syn.binding_score:.2f}",
            f"{syn.stability_score:.2f}",
            f"[dim]{ev_str}[/dim]",
        )

    app.print(Panel(
        table,
        title="Interaction State Classification",
        border_style="magenta", padding=(0, 1),
    ))


def _show_scientific_ranking(synthesis) -> None:
    lines: list[str] = []
    for i, (sys_name, score) in enumerate(synthesis.ranking, 1):
        syn = synthesis.systems.get(sys_name)
        color, _ = _STATE_STYLE.get(syn.primary_state if syn else "", ("dim", "?"))
        bar = "█" * int(score * 20) + "░" * (20 - int(score * 20))
        state_tag = (
            f"  [{color}]{syn.primary_state.replace('_', ' ')}[/{color}]" if syn else ""
        )
        lines.append(f"  [{i}] [cyan]{sys_name:<8}[/cyan] [{bar}] {score:.2f}{state_tag}")

    app.print(Panel(
        "\n".join(lines),
        title="Scientific composite ranking  (binding × 0.5 + stability × 0.3 + convergence × 0.2)",
        border_style="magenta", padding=(0, 1),
    ))


def _show_temporal_events(synthesis) -> None:
    _EVT_COLOR = {
        "abrupt_transition":    "red",
        "late_destabilization": "yellow",
        "contact_loss":         "yellow",
        "ligand_drift":         "magenta",
    }
    table = Table(box=box.SIMPLE, show_header=False, padding=(0, 1))
    table.add_column("icon", width=2)
    table.add_column("message")

    for evt in synthesis.events[:8]:
        color = _EVT_COLOR.get(evt.event_type, "dim")
        table.add_row(
            f"[{color}]▶[/{color}]",
            f"[{color}][bold]{evt.system}[/bold] {evt.replica}[/{color}]  "
            f"[dim]{evt.description}[/dim]",
        )

    border = "red" if any(e.event_type == "abrupt_transition" for e in synthesis.events) else "yellow"
    app.print(Panel(
        table,
        title=f"Time-resolved events  ({len(synthesis.events)})",
        border_style=border, padding=(0, 1),
    ))


@cli.command()
def study(
    path:   Path          = typer.Argument(Path("."), help="Directory containing XVG files (default: current directory)."),
    output: Optional[str] = typer.Option(None, "--output", "-o", help="Write JSON summary to file."),
    report: Optional[str] = typer.Option(None, "--report", "-r", help="Export Markdown report to file."),
):
    """Analyze a multi-system comparative MD study: auto-discovers systems, replicas,
    and observables, computes aggregate statistics, detects outliers, and produces
    a comparative summary.

    Expected filename convention:
      SYSTEM-REPLICAobservable.xvg    e.g.  AA-A1rmsd_protein.xvg
      SYSTEM-REPLICA_observable.xvg   e.g.  LP-A4_rmsd-ligand.xvg
    """
    from runtime.study_analyzer import parse_study

    if not path.exists():
        app.print(f"[red]Error:[/red] Path not found: {path}")
        raise typer.Exit(1)

    app.print(Panel(
        f"[bold cyan]SimForge Study[/bold cyan]  [dim]{path.resolve()}[/dim]",
        border_style="cyan", padding=(0, 2),
    ))

    with app.status("  Scanning and analyzing XVG files..."):
        study_obj = parse_study(path)

    # ── Discovery summary ─────────────────────────────────────────────────────
    app.print(
        f"\n  Discovered  [bold]{study_obj.n_xvg_discovered}[/bold] XVG files  "
        f"|  [green]{study_obj.n_xvg_parsed}[/green] parsed  "
        f"|  [dim]{study_obj.n_xvg_ungrouped} ungrouped[/dim]\n"
    )

    if study_obj.parse_errors:
        for err in study_obj.parse_errors[:5]:
            app.print(f"  [yellow]⚠[/yellow] {err}")
        if len(study_obj.parse_errors) > 5:
            app.print(f"  [dim]... and {len(study_obj.parse_errors) - 5} more errors[/dim]")
        app.print()

    if not study_obj.systems:
        app.print(Panel(
            "No multi-system structure detected in this directory.\n\n"
            "Expected naming pattern:\n"
            "  [dim]SYSTEM-REPLICAobservable.xvg[/dim]   e.g.  [dim]AA-A1rmsd_protein.xvg[/dim]\n"
            "  [dim]SYSTEM-REPLICA_observable.xvg[/dim]  e.g.  [dim]LP-A4_contacts.xvg[/dim]\n\n"
            "For single-simulation analysis use:  [dim]simforge analyze[/dim]",
            title="No study structure detected",
            border_style="yellow",
        ))
        raise typer.Exit(0)

    # ── Systems ───────────────────────────────────────────────────────────────
    sys_table = Table(box=box.SIMPLE, show_header=False, padding=(0, 1))
    sys_table.add_column("bullet", width=2)
    sys_table.add_column("system", style="bold cyan")
    sys_table.add_column("replicas", style="dim")
    sys_table.add_column("observables", style="dim")

    for sys_name, sg in sorted(study_obj.systems.items()):
        obs_str = "  ".join(
            study_obj.observable_display.get(o, o) for o in sg.observables
        )
        sys_table.add_row(
            "[cyan]●[/cyan]", sys_name,
            f"{sg.n_replicas} replica{'s' if sg.n_replicas != 1 else ''}",
            obs_str,
        )

    app.print(Panel(
        sys_table,
        title=f"Systems detected  ({len(study_obj.systems)})",
        border_style="cyan", padding=(0, 1),
    ))

    # ── Observables ───────────────────────────────────────────────────────────
    if study_obj.observables_detected:
        obs_table = Table(box=box.SIMPLE, show_header=False, padding=(0, 1))
        obs_table.add_column("bullet", width=2)
        obs_table.add_column("display", style="white")
        obs_table.add_column("units",   style="dim")
        obs_table.add_column("group",   style="dim")

        from runtime.observable_resolver import ObservableResolver as _OR
        _res = _OR()
        for obs_name in study_obj.observables_detected:
            ro = _res.resolve(obs_name)
            obs_table.add_row(
                "[blue]●[/blue]",
                study_obj.observable_display.get(obs_name, obs_name),
                study_obj.observable_units.get(obs_name, ""),
                ro.group,
            )

        app.print(Panel(
            obs_table,
            title=f"Observables detected  ({len(study_obj.observables_detected)})",
            border_style="blue", padding=(0, 1),
        ))

    # ── Aggregate statistics table ────────────────────────────────────────────
    obs_cols = study_obj.observables_detected[:6]
    if obs_cols and len(study_obj.systems) > 0:
        agg_table = Table(box=box.SIMPLE, show_header=True, header_style="bold dim")
        agg_table.add_column("System", style="bold cyan")
        for obs_name in obs_cols:
            display = study_obj.observable_display.get(obs_name, obs_name)
            units   = study_obj.observable_units.get(obs_name, "")
            hdr     = f"{display}" + (f"\nmean ± std ({units})" if units else "\nmean ± std")
            agg_table.add_column(hdr, justify="right")

        for sys_name, sg in sorted(study_obj.systems.items()):
            row: list[str] = [sys_name]
            for obs_name in obs_cols:
                if obs_name in sg.aggregate:
                    ag = sg.aggregate[obs_name]
                    row.append(f"{ag.mean:.3g} ± {ag.std:.3g}")
                else:
                    row.append("[dim]—[/dim]")
            agg_table.add_row(*row)

        app.print(Panel(
            agg_table,
            title="Aggregate statistics  (mean ± inter-replica std)",
            border_style="dim", padding=(0, 1),
        ))

    summary_obj = study_obj.summary

    # ── Scientific synthesis ──────────────────────────────────────────────────
    synthesis = None
    if study_obj.systems and study_obj.observables_detected:
        try:
            from runtime.scientific_synthesis import synthesize_study as _synth
            with app.status("  Running scientific synthesis..."):
                synthesis = _synth(study_obj)
        except Exception as _e:
            app.print(f"  [dim]Scientific synthesis skipped: {_e}[/dim]")

    # ── Outliers ──────────────────────────────────────────────────────────────
    if summary_obj and summary_obj.outlier_replicas:
        out_table = Table(box=box.SIMPLE, show_header=False, padding=(0, 1))
        out_table.add_column("icon", width=2)
        out_table.add_column("message")
        for sys_n, rep_n, reason in summary_obj.outlier_replicas:
            out_table.add_row(
                "[yellow]⚠[/yellow]",
                f"[bold]{sys_n}[/bold] replica [bold]{rep_n}[/bold]  [dim]{reason}[/dim]",
            )
        app.print(Panel(
            out_table,
            title=f"Potential outliers  ({len(summary_obj.outlier_replicas)})",
            border_style="yellow", padding=(0, 1),
        ))

    # ── Comparative findings ──────────────────────────────────────────────────
    if summary_obj and summary_obj.findings:
        find_table = Table(box=box.SIMPLE, show_header=False, padding=(0, 1))
        find_table.add_column("icon", width=2)
        find_table.add_column("message")

        sorted_findings = sorted(
            summary_obj.findings,
            key=lambda f: (0 if f.level != "info" else 1, f.level),
        )

        for finding in sorted_findings[:12]:
            color = _FIND_LEVEL_COLOR.get(finding.level, "dim")
            icon  = _FIND_LEVEL_ICON.get(finding.level, "·")
            find_table.add_row(
                f"[{color}]{icon}[/{color}]",
                f"[{color}]{finding.message}[/{color}]",
            )

        app.print(Panel(
            find_table,
            title="Comparative analysis",
            border_style="cyan", padding=(0, 1),
        ))

    # ── Interaction state classification ──────────────────────────────────────
    if synthesis and synthesis.systems:
        _show_interaction_states(synthesis, study_obj)

    # ── Time-resolved events ──────────────────────────────────────────────────
    if synthesis and synthesis.events:
        _show_temporal_events(synthesis)

    # ── Scientific ranking ────────────────────────────────────────────────────
    if synthesis and synthesis.ranking and len(synthesis.ranking) > 1:
        _show_scientific_ranking(synthesis)
    elif summary_obj and len(summary_obj.system_ranking) > 1:
        ranking = sorted(summary_obj.system_ranking.items(), key=lambda x: -x[1])
        rank_lines: list[str] = []
        for i, (sn, score) in enumerate(ranking, 1):
            bar_filled = int(score * 20)
            bar = "█" * bar_filled + "░" * (20 - bar_filled)
            rank_lines.append(f"  [{i}] [cyan]{sn:<8}[/cyan] [{bar}] {score:.2f}")
        app.print(Panel(
            "\n".join(rank_lines),
            title="System stability ranking  (convergence-based)",
            border_style="dim", padding=(0, 1),
        ))

    # ── Scientific narrative ──────────────────────────────────────────────────
    if synthesis and synthesis.narrative:
        app.print(Panel(
            synthesis.narrative,
            title="Scientific Narrative",
            border_style="magenta", padding=(1, 2),
        ))

    # ── JSON output ───────────────────────────────────────────────────────────
    if output:
        import json as _json
        out_data: dict = {
            "path":             str(path.resolve()),
            "n_xvg_discovered": study_obj.n_xvg_discovered,
            "n_xvg_parsed":     study_obj.n_xvg_parsed,
            "systems": {
                sn: {
                    "n_replicas":  sg.n_replicas,
                    "observables": sg.observables,
                    "aggregate": {
                        obs: {"mean": ag.mean, "std": ag.std, "n_replicas": ag.n_replicas}
                        for obs, ag in sg.aggregate.items()
                    },
                }
                for sn, sg in study_obj.systems.items()
            },
            "outliers": [
                {"system": s, "replica": r, "reason": reason}
                for s, r, reason in (summary_obj.outlier_replicas if summary_obj else [])
            ],
            **({"synthesis": {
                "ranking": [
                    {"system": sn, "composite_score": sc}
                    for sn, sc in synthesis.ranking
                ],
                "systems": {
                    sn: {
                        "primary_state":      syn.primary_state,
                        "primary_confidence": syn.primary_confidence,
                        "binding_score":      syn.binding_score,
                        "stability_score":    syn.stability_score,
                        "composite_score":    syn.composite_score,
                        "evidence":           syn.evidence,
                    }
                    for sn, syn in synthesis.systems.items()
                },
                "events": [
                    {
                        "system": e.system, "replica": e.replica,
                        "observable": e.observable, "event_type": e.event_type,
                        "time_ns": e.time_ns, "description": e.description,
                    }
                    for e in synthesis.events
                ],
                "narrative": synthesis.narrative,
            }} if synthesis else {}),
        }
        Path(output).write_text(_json.dumps(out_data, indent=2))
        app.print(f"\n[dim]JSON written to {output}[/dim]")

    # ── Markdown report ───────────────────────────────────────────────────────
    if report:
        md = _generate_markdown_report(study_obj, summary_obj, synthesis, path)
        Path(report).write_text(md, encoding="utf-8")
        app.print(f"\n[dim]Markdown report written to {report}[/dim]")


# ═══════════════════════════════════════════════════════════════════════════════
# annotate-structure
# ═══════════════════════════════════════════════════════════════════════════════


def _ask_range_list(label: str, existing: list[str]) -> list[str]:
    """Prompt the user to enter a list of residue range strings (e.g. '1-50').

    Shows current values, lets user add or replace entries.
    Empty input on first prompt keeps existing list.
    """
    app.print(f"\n  [bold]{label}[/bold]")
    app.print("  [dim]Formato: '1-50', '1-20,45-60', o '5,10,15'. Línea vacía para terminar.[/dim]")
    if existing:
        app.print(f"  [dim]Actuales: {existing}[/dim]")
        raw = input("  ¿Reemplazar? (s/n, Enter=n): ").strip().lower()
        if raw != "s":
            return existing

    from core.structural_annotation import parse_residue_range
    entries: list[str] = []
    while True:
        raw = input(f"  {'Rango ' + str(len(entries)+1)}: ").strip()
        if not raw:
            break
        try:
            parse_residue_range(raw)
            entries.append(raw)
        except ValueError as e:
            app.print(f"  [red]{e}[/red]")
    return entries


def _ask_tm_segments(existing: list) -> list:
    """Prompt for transmembrane segments (string or {residues, label, helix_type})."""
    from core.structural_annotation import parse_residue_range, TransmembraneSegment

    app.print("\n  [bold]Segmentos transmembrana[/bold] [dim](opcional)[/dim]")
    app.print("  [dim]Formato: '51-75' o entrada detallada con label/tipo.[/dim]")

    if existing:
        names = [
            (s.residues if isinstance(s, TransmembraneSegment) else s)
            for s in existing
        ]
        app.print(f"  [dim]Actuales: {names}[/dim]")
        raw = input("  ¿Reemplazar? (s/n, Enter=n): ").strip().lower()
        if raw != "s":
            return existing

    mode = _ask("Modo de entrada", [
        ("simple",   "Simple — solo rangos de residuos"),
        ("detailed", "Detallado — rangos + label + tipo de hélice"),
    ])

    segments: list = []
    i = 1
    while True:
        raw_r = input(f"  TM{i} residuos (Enter para terminar): ").strip()
        if not raw_r:
            break
        try:
            parse_residue_range(raw_r)
        except ValueError as e:
            app.print(f"  [red]{e}[/red]")
            continue

        if mode == "detailed":
            label     = input(f"  TM{i} label (ej: TM1, Enter=TM{i}): ").strip() or f"TM{i}"
            helix_raw = _ask("Tipo de hélice", [
                ("alpha", "α-hélice (más común)"),
                ("310",   "Hélice 3₁₀"),
                ("pi",    "Hélice π"),
                ("none",  "No especificar"),
            ])
            segments.append(TransmembraneSegment(
                residues=raw_r,
                label=label,
                helix_type=helix_raw if helix_raw != "none" else None,
            ))
        else:
            segments.append(raw_r)
        i += 1
    return segments


def _annotation_to_yaml_dict(ann) -> dict:
    """Convert a StructuralAnnotation to a plain dict suitable for YAML serialization."""
    from core.structural_annotation import TransmembraneSegment

    def _tm_entry(s):
        if isinstance(s, TransmembraneSegment):
            d: dict = {"residues": s.residues}
            if s.label:
                d["label"] = s.label
            if s.helix_type:
                d["helix_type"] = s.helix_type
            return d
        return s

    result: dict = {}

    if ann.membrane_topology:
        mt = ann.membrane_topology
        topo: dict = {}
        if mt.extracellular_regions:
            topo["extracellular_regions"] = mt.extracellular_regions
        if mt.intracellular_regions:
            topo["intracellular_regions"] = mt.intracellular_regions
        if mt.transmembrane_segments:
            topo["transmembrane_segments"] = [_tm_entry(s) for s in mt.transmembrane_segments]
        result["membrane_topology"] = topo

    if ann.orientation:
        o = ann.orientation
        result["orientation"] = {
            "extracellular_side": o.extracellular_side,
            "intracellular_side": o.intracellular_side,
            "source":             o.source,
            "confidence":         o.confidence,
        }
        if o.reference:
            result["orientation"]["reference"] = o.reference

    if ann.evidence:
        evlist = []
        for ev in ann.evidence:
            d = {"source": ev.source, "confidence": ev.confidence}
            if ev.reference:
                d["reference"] = ev.reference
            if ev.notes:
                d["notes"] = ev.notes
            evlist.append(d)
        result["evidence"] = evlist

    return result


@cli.command("annotate-structure")
def annotate_structure(
    yaml_path: Path = typer.Argument(..., help="Config YAML a anotar."),
    output: Optional[Path] = typer.Option(
        None, "--output", "-o",
        help="Ruta de salida (default: sobreescribe el YAML de entrada).",
    ),
):
    """Wizard interactivo — define la Structural Biology Annotation Layer de un YAML.

    Permite especificar topología de membrana, dominios biológicos, orientación
    geométrica y evidencia. El resultado se escribe de vuelta al YAML como
    la clave top-level `structural_annotation:`.
    """
    from core.structural_annotation import (
        StructuralAnnotation,
        MembraneTopologyAnnotation,
        OrientationAnnotation,
        OrientationEvidence,
    )
    from ruamel.yaml import YAML as RuamelYAML
    import yaml as pyyaml

    if not yaml_path.exists():
        app.print(f"[red]Archivo no encontrado: {yaml_path}[/red]")
        raise typer.Exit(1)

    app.print(Panel(
        "[bold cyan]SimForge · Annotate Structure[/bold cyan]  "
        "[dim]Structural Biology Annotation Layer[/dim]",
        border_style="cyan",
    ))

    # ── Cargar YAML y mostrar estado actual ───────────────────────────────────
    with open(yaml_path) as f:
        raw = pyyaml.safe_load(f)

    project_name = (raw.get("project") or {}).get("name", yaml_path.stem)
    app.print(f"\n  Proyecto: [bold]{project_name}[/bold]")

    existing_ann: Optional[StructuralAnnotation] = None
    sa_raw = raw.get("structural_annotation")
    if sa_raw:
        try:
            existing_ann = StructuralAnnotation(**sa_raw)
            app.print("  [dim]structural_annotation existente detectada.[/dim]")
        except Exception:
            app.print("  [yellow]structural_annotation existente tiene errores — se reemplazará.[/yellow]")

    existing_mt = existing_ann.membrane_topology if existing_ann else None
    existing_orient = existing_ann.orientation if existing_ann else None
    existing_evidence = existing_ann.evidence if existing_ann else []

    # ── 1. Topología extracellular ────────────────────────────────────────────
    app.print("\n[bold]── 1 · Regiones extracelulares ─────────────────────────────────────[/bold]")
    ec = _ask_range_list("Residuos en la cara extracelular", existing_mt.extracellular_regions if existing_mt else [])

    # ── 2. Topología intracelular ─────────────────────────────────────────────
    app.print("\n[bold]── 2 · Regiones intracelulares ──────────────────────────────────────[/bold]")
    if ec:
        app.print(
            "  [dim]Necesario para orientación automática: sin IC no se puede calcular el eje EC→IC.[/dim]"
        )
    ic = _ask_range_list("Residuos en la cara intracelular", existing_mt.intracellular_regions if existing_mt else [])

    if ec and not ic:
        app.print(
            "\n  [yellow]⚠[/yellow]  [yellow]orient_protein quedará GUIDED:[/yellow] "
            "hay regiones EC pero faltan regiones IC.\n"
            "     Sin al menos una región intracelular no se puede calcular el eje EC→IC.\n"
            "     Puedes continuar el wizard y añadir IC, o volver a ejecutarlo después."
        )

    # ── 3. Segmentos transmembrana ────────────────────────────────────────────
    app.print("\n[bold]── 3 · Segmentos transmembrana ──────────────────────────────────────[/bold]")
    tm = _ask_tm_segments(existing_mt.transmembrane_segments if existing_mt else [])

    # ── 4. Orientación geométrica ─────────────────────────────────────────────
    app.print("\n[bold]── 4 · Orientación geométrica ───────────────────────────────────────[/bold]")
    app.print("  [dim]¿En qué dirección del eje Z queda la cara extracelular después de editconf -princ?[/dim]")
    ec_side = _ask("Cara extracelular", [
        ("+z", "+Z  (convención GROMACS estándar)"),
        ("-z", "-Z  (proteína invertida respecto a la convención)"),
    ], default=1 if (not existing_orient or existing_orient.extracellular_side == "+z") else 2)

    ic_side = "-z" if ec_side == "+z" else "+z"

    # ── 5. Fuente y confianza ─────────────────────────────────────────────────
    app.print("\n[bold]── 5 · Fuente y confianza ───────────────────────────────────────────[/bold]")
    source = _ask("Fuente de la anotación", [
        ("user_annotation", "Anotación manual (tú lo sabes)"),
        ("opm_database",    "OPM database (Orientations of Proteins in Membranes)"),
        ("pdbtm",           "PDBTM (PDB Transmembrane database)"),
        ("uniprot",         "UniProt (feature annotations)"),
        ("predicted",       "Predicción computacional (TMH, DeepTMHMM, etc.)"),
    ])

    default_conf = existing_orient.confidence if existing_orient else 0.9
    confidence = _ask_float("Confianza [0.0–1.0]", default_conf)
    confidence = max(0.0, min(1.0, confidence))

    ref = _ask_str("Referencia (OPM accession, PubMed ID, URL — opcional)")

    # ── Construir StructuralAnnotation ────────────────────────────────────────
    topology = MembraneTopologyAnnotation(
        extracellular_regions=ec,
        intracellular_regions=ic,
        transmembrane_segments=tm,
    )
    orientation = OrientationAnnotation(
        extracellular_side=ec_side,
        intracellular_side=ic_side,
        source=source,
        confidence=confidence,
        reference=ref or None,
    ) if (ec or ic) else None

    evidence_note = (
        "Anotación creada con simforge annotate-structure"
        + (f"; referencia: {ref}" if ref else "")
    )
    evidence = [OrientationEvidence(
        source=source,
        confidence=confidence,
        reference=ref or None,
        notes=evidence_note,
    )] if (ec or ic) else existing_evidence

    ann = StructuralAnnotation(
        membrane_topology=topology if (ec or ic or tm) else None,
        orientation=orientation,
        evidence=evidence,
    )

    # ── Mostrar validación ────────────────────────────────────────────────────
    warns = ann.validation_warnings()
    if warns:
        app.print("\n[yellow]Advertencias de consistencia:[/yellow]")
        for w in warns:
            app.print(f"  [yellow]⚠[/yellow]  {w}")

    overlap_warns = topology.overlap_warnings() if (ec or ic or tm) else []
    if overlap_warns:
        app.print("\n[red]Solapamiento de regiones:[/red]")
        for w in overlap_warns:
            app.print(f"  [red]✗[/red]  {w}")

    # ── Mostrar AutomationLevel resultante ────────────────────────────────────
    app.print("\n")
    if ann.is_complete_for_orient():
        app.print("  [green]● orient_protein → AUTOMATED[/green]  "
                  "[dim](EC + IC + orientación definidas)[/dim]")
    elif ann.is_partial_for_orient():
        app.print("  [green]● orient_protein → AUTOMATED[/green]  "
                  "[dim](EC + IC definidas; eje asumido como estándar +Z)[/dim]")
    else:
        if ec and not ic:
            guided_reason = "falta intracellular_regions — no se puede calcular el eje EC→IC"
            guided_hint   = "Define al menos una región IC (ej. '200-250') para habilitar orientación automática."
        elif ic and not ec:
            guided_reason = "falta extracellular_regions — no se puede calcular el eje EC→IC"
            guided_hint   = "Define al menos una región EC (ej. '1-50') para habilitar orientación automática."
        else:
            guided_reason = "sin topología de membrana definida"
            guided_hint   = "Define regiones EC e IC con 'simforge annotate-structure <config.yaml>'."
        app.print(f"  [red]● orient_protein → GUIDED[/red]  [dim]({guided_reason})[/dim]")
        app.print(f"  [dim]  → {guided_hint}[/dim]")

    # ── Confirmar escritura ───────────────────────────────────────────────────
    app.print()
    confirm = input("  ¿Escribir structural_annotation al YAML? (s/n, Enter=s): ").strip().lower()
    if confirm == "n":
        app.print("[dim]Cancelado.[/dim]")
        raise typer.Exit(0)

    # ── Escribir con ruamel.yaml (preserva formato y comentarios) ─────────────
    ryaml = RuamelYAML()
    ryaml.preserve_quotes = True
    ryaml.default_flow_style = False
    ryaml.best_sequence_indent = 2
    ryaml.best_map_flow_style  = False
    ryaml.width = 4096  # avoid line-wrapping in long lists

    with open(yaml_path) as f:
        doc = ryaml.load(f)

    doc["structural_annotation"] = _annotation_to_yaml_dict(ann)

    out_path = output or yaml_path
    with open(out_path, "w") as f:
        ryaml.dump(doc, f)

    app.print(f"\n[green]✓[/green] structural_annotation escrita en [bold]{out_path}[/bold]")
    if ann.membrane_topology:
        mt = ann.membrane_topology
        app.print(f"  EC: {mt.extracellular_regions}  IC: {mt.intracellular_regions}  "
                  f"TM: {mt.tm_segment_count()} segmentos")


# ═══════════════════════════════════════════════════════════════════════════════
# Ligand commands
# ═══════════════════════════════════════════════════════════════════════════════

_RDKIT_MISSING_MSG = (
    "LigParGen export requires RDKit. "
    "Activate the rdkit_env environment or install RDKit.\n"
    "  conda activate rdkit_env\n"
    "  # or: conda install -c conda-forge rdkit"
)

_RMSD_WARN_THRESHOLD = 0.05  # Å — mirrors export._RMSD_WARN_THRESHOLD


@_ligand_app.command("export-ligpargen")
def export_ligpargen_cmd(
    input_file: Path = typer.Argument(
        ...,
        help="Input ligand file (.sdf, .mol, .pdb).",
        exists=False,  # validated manually for a cleaner error message
    ),
    output_dir: Path = typer.Option(
        Path("ligpargen_export"),
        "--output-dir", "-o",
        help="Directory to write output files (default: ./ligpargen_export).",
    ),
    mol_name: str = typer.Option(
        "LIG",
        "--mol-name",
        help="4-char GROMACS-compatible molecule name (default: LIG).",
    ),
    smiles: bool = typer.Option(
        False,
        "--smiles",
        help=(
            "Export canonical SMILES (.smi + metadata + charge advisory). "
            "Also available alongside --legacy output."
        ),
    ),
    legacy: bool = typer.Option(
        False,
        "--legacy/--no-legacy",
        help=(
            "Write the experimentally validated legacy PDB (ATOM records + CONECT). "
            "Also writes companion SMILES, metadata JSON, and charge advisory."
        ),
    ),
) -> None:
    """Export a ligand for LigParGen parameterization.

    --legacy produces a PDB that has been experimentally accepted by the
    LigParGen web server, plus companion files in the output directory:

    \b
        LIG_ligpargen_legacy.pdb   ← upload to LigParGen
        LIG_ligpargen.smi          ← canonical SMILES companion
        LIG_meta.json              ← metadata including formal charge
        LIG_charge.txt             ← charge advisory (ALWAYS check this)

    --smiles exports canonical SMILES and is also accepted by LigParGen.

    IMPORTANT — formal charge:
    The reported formal charge must match the charge you select in LigParGen.
    A charge mismatch (e.g. submitting charge=0 for a +1 molecule) is the
    most common cause of failed or incorrect parameterization.

    Examples:

        simforge ligand export-ligpargen E20.pdb --legacy   # validated PDB
        simforge ligand export-ligpargen E20.pdb --smiles   # SMILES mode
    """
    if not input_file.exists():
        app.print(f"[red]Error:[/red] Input file not found: {input_file}")
        raise typer.Exit(1)

    from ligand import export as _lex

    if smiles:
        fmt_label = "SMILES (.smi)"
    elif legacy:
        fmt_label = "legacy PDB (ATOM records) — experimentally validated"
    else:
        fmt_label = "modern PDB (HETATM records)"

    app.print(Panel(
        f"[bold cyan]LigParGen Export[/bold cyan]  "
        f"[dim]{input_file.name}[/dim]  →  [dim]{fmt_label}[/dim]",
        border_style="cyan", padding=(0, 2),
    ))

    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        if smiles:
            result = _lex.export_for_ligpargen_smiles(
                path=input_file,
                output_dir=output_dir,
                mol_name=mol_name,
            )
        elif legacy:
            result = _lex.export_for_ligpargen_legacy(
                path=input_file,
                output_dir=output_dir,
                mol_name=mol_name,
            )
        else:
            result = _lex.export_for_ligpargen(
                path=input_file,
                output_dir=output_dir,
                mol_name=mol_name,
            )
    except ImportError:
        app.print(f"[red]Error:[/red] {_RDKIT_MISSING_MSG}")
        raise typer.Exit(1)

    if not result.success:
        app.print(f"\n  [red]✗  Export failed:[/red] {result.error}")
        raise typer.Exit(1)

    # ── Success output ────────────────────────────────────────────────────────
    app.print(f"\n  [green]✓[/green]  [bold]Exported:[/bold]  {result.exported_path}")
    app.print(f"     Molecule name:  [bold]{result.molecule_name}[/bold]")

    # ── Charge display (shown for all modes when available) ───────────────────
    if result.formal_charge is not None:
        from ligand.export import _charge_label as _clabel
        clabel = _clabel(result.formal_charge)
        app.print(f"     Formal charge:  [bold]{clabel}[/bold]")
        app.print(f"     LigParGen charge: [bold]{clabel}[/bold]")

    if smiles:
        app.print(f"     SMILES:         [bold]{result.smiles}[/bold]")
        app.print(f"     Heavy atoms:    {result.atom_count}")
        meta_path = output_dir / f"{result.molecule_name}_meta.json"
        charge_path = output_dir / f"{result.molecule_name}_charge.txt"
        if meta_path.exists():
            app.print(f"     Metadata:       {meta_path}")
        if charge_path.exists():
            app.print(f"     Charge file:    {charge_path}")
    else:
        app.print(f"     Atom count:     {result.atom_count}  (including explicit H)")
        app.print(f"     Format:         {fmt_label}")

        if result.heavy_atom_rmsd is not None:
            rmsd_str = f"{result.heavy_atom_rmsd:.4f} Å"
            if result.heavy_atom_rmsd > _RMSD_WARN_THRESHOLD:
                app.print(
                    f"\n  [yellow]⚠[/yellow]  Heavy-atom coordinate RMSD: [yellow]{rmsd_str}[/yellow]"
                    f"  (> {_RMSD_WARN_THRESHOLD} Å — coordinates may have shifted "
                    f"during hydrogen addition; verify in PyMOL or Avogadro)"
                )
            else:
                app.print(f"     Heavy-atom RMSD: {rmsd_str}  [dim](within tolerance)[/dim]")

        if legacy:
            name = result.molecule_name
            app.print(f"\n  [dim]Companion files written:[/dim]")
            app.print(f"     {output_dir / f'{name}_ligpargen.smi'}")
            app.print(f"     {output_dir / f'{name}_meta.json'}")
            app.print(f"     {output_dir / f'{name}_charge.txt'}")

    # ── Charge warning ────────────────────────────────────────────────────────
    if result.formal_charge is not None and result.formal_charge != 0:
        from ligand.export import _charge_label as _clabel
        clabel = _clabel(result.formal_charge)
        app.print(Panel(
            f"[bold]WARNING[/bold]\n"
            f"This molecule carries a formal charge of [bold]{clabel}[/bold].\n"
            f"In LigParGen select charge [bold]{clabel}[/bold] instead of 0.",
            border_style="yellow", padding=(0, 2),
        ))

    app.print(f"\n  [dim]Next step:[/dim] submit [bold]{result.exported_path.name}[/bold] "
              f"to https://zarbi.chem.yale.edu/ligpargen/\n")


# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    cli()
