"""
CLI implementations for the SimForge platform-consolidation commands:

    simforge build            SYSTEM.yaml   -> assembled GROMACS system
    simforge doctor                          -> environment report
    simforge validate-system  SYSTEM_DIR     -> structured system validation
    simforge inspect-run      RUN_DIR        -> reproducibility manifest summary

These are thin presentation wrappers around ``core.system_spec``,
``core.system_build``, ``core.doctor`` and ``core.system_validate``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Optional

import typer
from rich import box
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

app = Console()

_STATUS_COLOR = {
    "PASS": "green", "AVAILABLE": "green",
    "WARN": "yellow", "SKIP": "dim", "WRONG_ENVIRONMENT": "yellow",
    "FAIL": "red", "MISSING": "red",
}
_STATUS_ICON = {
    "PASS": "✓", "AVAILABLE": "✓",
    "WARN": "⚠", "SKIP": "–", "WRONG_ENVIRONMENT": "⚠",
    "FAIL": "✗", "MISSING": "✗",
}


def _fmt(status: str) -> str:
    c = _STATUS_COLOR.get(status, "white")
    return f"[{c}]{_STATUS_ICON.get(status, '·')} {status}[/{c}]"


# ═══════════════════════════════════════════════════════════════════════════════
# simforge doctor
# ═══════════════════════════════════════════════════════════════════════════════

def doctor_fn(
    forcefield: str = typer.Option("oplsaa", "--forcefield", "--ff", help="Base force field to verify."),
    json_out: bool = typer.Option(False, "--json", help="Emit a machine-readable report."),
) -> None:
    """Inspect the runtime environment and scientific software dependencies."""
    from core.doctor import run_doctor

    report = run_doctor(forcefield=forcefield)

    if json_out:
        print(json.dumps(report.to_dict(), indent=2))
        raise typer.Exit(0 if report.ready else 1)

    table = Table(box=box.SIMPLE, show_header=True, header_style="bold dim")
    table.add_column("Component")
    table.add_column("Level", style="dim")
    table.add_column("Version")
    table.add_column("Status")
    table.add_column("Detail", style="dim", overflow="fold")
    for c in report.checks:
        table.add_row(
            c.name, c.level, c.version or "—", _fmt(c.status), c.detail or "",
        )

    app.print(Panel(table, title="SIMFORGE ENVIRONMENT", border_style="cyan", padding=(0, 1)))
    color = "green" if report.ready else "red"
    app.print(f"\n  Overall status: [{color}][bold]{report.overall}[/bold][/{color}]\n")
    raise typer.Exit(0 if report.ready else 1)


# ═══════════════════════════════════════════════════════════════════════════════
# simforge validate-system
# ═══════════════════════════════════════════════════════════════════════════════

def validate_system_fn(
    system_dir: Path = typer.Argument(..., help="Built GROMACS system directory (contains topol.top)."),
    forcefield: str = typer.Option("oplsaa", "--forcefield", "--ff"),
    coordinate: Optional[Path] = typer.Option(None, "--coordinate", "-c", help="Explicit coordinate .gro (else auto-detected)."),
    no_grompp: bool = typer.Option(False, "--no-grompp", help="Skip the gmx grompp dry-run."),
    json_report: Optional[Path] = typer.Option(None, "--json-report", help="Also write validation_report.json here."),
) -> None:
    """Validate a built system without modifying it. Non-zero exit on blocking failure."""
    from core.system_validate import validate_system_dir

    report = validate_system_dir(
        system_dir, forcefield=forcefield, coordinate=coordinate,
        run_grompp=(False if no_grompp else None),
    )

    table = Table(box=box.SIMPLE, show_header=True, header_style="bold dim")
    table.add_column("Check")
    table.add_column("Severity", style="dim")
    table.add_column("Status")
    table.add_column("Message", overflow="fold")
    for c in report.checks:
        table.add_row(c.check_id, c.severity, _fmt(c.status), c.message)

    app.print(Panel(table, title=f"SYSTEM VALIDATION — {report.system_dir}",
                    border_style="cyan", padding=(0, 1)))

    color = {"READY": "green", "READY (with warnings)": "yellow"}.get(report.readiness, "red")
    app.print(f"\n  Scientific readiness: [{color}][bold]{report.readiness}[/bold][/{color}]\n")

    out = json_report or (Path(report.system_dir) / "validation_report.json")
    try:
        out.write_text(json.dumps(report.to_dict(), indent=2, default=str))
        app.print(f"  [dim]report → {out}[/dim]\n")
    except Exception:
        pass

    raise typer.Exit(0 if report.ok else 1)


# ═══════════════════════════════════════════════════════════════════════════════
# simforge build
# ═══════════════════════════════════════════════════════════════════════════════

def build_fn(
    spec_path: Path = typer.Argument(..., help="Unified system specification YAML."),
    out_dir: Optional[Path] = typer.Option(None, "--out", "-o", help="Output directory (default: ./<name>_system)."),
    no_grompp: bool = typer.Option(False, "--no-grompp", help="Skip gmx grompp validation during assembly + check."),
    skip_preflight: bool = typer.Option(False, "--skip-preflight", help="Skip the deterministic atomtype preflight."),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Show full atomtype rename maps."),
) -> None:
    """Build a GROMACS-ready system from a single declarative YAML specification."""
    from core.system_spec import SystemSpecError, load_system_spec
    from core.system_build import build_system

    if not spec_path.exists():
        app.print(f"[red]Error:[/red] spec file not found: {spec_path}")
        raise typer.Exit(1)

    try:
        spec = load_system_spec(spec_path)
    except SystemSpecError as exc:
        app.print(Panel(
            f"[red]Invalid system specification:[/red]\n\n{exc}",
            title="❌ Build Failed", border_style="red",
        ))
        raise typer.Exit(1)

    resolved_out = out_dir or Path.cwd() / f"{spec.name or spec_path.stem}_system"

    # ── header ────────────────────────────────────────────────────────────
    comp_desc = " + ".join(f"[cyan]{c.id}[/cyan]" for c in spec.components)
    app.print(Panel(
        f"[bold cyan]SimForge Build[/bold cyan]  "
        f"[dim]{spec.protein_structure_path.name}[/dim] ({spec.protein_mode}) + {comp_desc}"
        f"  →  [dim]{resolved_out}[/dim]",
        border_style="cyan", padding=(0, 2),
    ))
    for note in spec.discovery_notes:
        app.print(f"  [green]↪[/green] {note}")
    for w in spec.warnings:
        app.print(f"  [yellow]⚠[/yellow] {w}")

    command = ["simforge", "build", str(spec_path)]
    if out_dir:
        command += ["--out", str(out_dir)]

    with app.status("  Assembling system..."):
        outcome = build_system(
            spec, resolved_out, run_grompp=not no_grompp,
            command=command, skip_preflight=skip_preflight,
        )

    if outcome.preflight_errors:
        app.print(Panel(
            "[red]Atomtype preflight failed — assembly not attempted:[/red]\n\n"
            + "\n\n".join(f"  • {e}" for e in outcome.preflight_errors),
            title="❌ Build Failed (preflight)", border_style="red",
        ))
        raise typer.Exit(1)

    if outcome.errors and not outcome.assembly_result:
        app.print(Panel(
            "[red]Assembly failed:[/red]\n\n" + "\n".join(f"  • {e}" for e in outcome.errors),
            title="❌ Build Failed", border_style="red",
        ))
        raise typer.Exit(1)

    _render_build_summary(outcome, verbose=verbose)

    if outcome.errors:
        app.print(Panel(
            "[red]Validation errors:[/red]\n" + "\n".join(f"  • {e}" for e in outcome.errors),
            border_style="red", title="Build finished with errors",
        ))
        raise typer.Exit(1)

    app.print(Panel(
        f"[green]✓[/green] System   → [cyan]{outcome.out_dir}[/cyan]\n"
        f"  Topology → [dim]{outcome.out_dir}/topol.top[/dim]\n"
        f"  Coords   → [dim]{outcome.out_dir}/complex.gro[/dim]\n"
        f"  Provenance → [dim]{outcome.out_dir}/provenance.json[/dim]\n"
        f"  Validation → [dim]{outcome.out_dir}/validation_report.json[/dim]\n\n"
        f"Next:  [dim]gmx solvate -cp {outcome.out_dir}/complex.gro -cs spc216.gro "
        f"-o {outcome.out_dir}/solvated.gro -p {outcome.out_dir}/topol.top[/dim]",
        title="✓ Build complete", border_style="green", padding=(0, 2),
    ))


def _render_build_summary(outcome, *, verbose: bool) -> None:
    r = outcome.assembly_result
    spec = outcome.spec

    lines = [
        "  [bold]Protein[/bold]",
        f"    input mode: {spec.protein_mode}",
    ]
    if spec.protein_mode == "preparameterized":
        lines.append(f"    chains: {len(spec.protein_topology_paths)}")
    lines.append(f"    atoms: {r.protein_atom_count}")
    for mol in r.protein_molecules:
        lines.append(f"      [dim]{mol['name']} × {mol['count']}[/dim]")

    lines.append("")
    lines.append("  [bold]Components[/bold]")
    total_renames = 0
    for comp in r.components:
        role = f"  [dim]{comp['role']}[/dim]" if comp.get("role") else ""
        lines.append(
            f"    {comp['id']:<8} {comp['molecule_name']:<8} "
            f"{comp['atom_count']:>4} atoms{role}"
        )
        renames = comp.get("atomtype_renames") or {}
        total_renames += len(renames)
        if verbose and renames:
            for old, new in renames.items():
                lines.append(f"        [dim]{old} → {new}[/dim]")

    lines.append("")
    if total_renames:
        lines.append(f"  [bold]Atomtype conflicts[/bold]")
        lines.append(
            f"    {total_renames} resolved automatically"
            + ("" if verbose else "  [dim](use --verbose for the rename map)[/dim]")
        )
        lines.append("")
    lines.append(f"  [bold]Total atoms[/bold]  {r.total_atom_count}")

    lines.append("")
    v = outcome.validation
    vc = {"READY": "green", "READY (with warnings)": "yellow"}.get(v.get("readiness"), "red")
    lines.append(f"  [bold]Validation[/bold]  [{vc}]{v.get('readiness', '?')}[/{vc}]")
    for chk in v.get("checks", []):
        if chk["status"] in ("FAIL", "WARN"):
            lines.append(f"    {_fmt(chk['status'])}  {chk['check_id']}: {chk['message']}")

    app.print(Panel("\n".join(lines), title="System", border_style="cyan", padding=(0, 2)))


# ═══════════════════════════════════════════════════════════════════════════════
# simforge campaign-build
# ═══════════════════════════════════════════════════════════════════════════════

def campaign_build_fn(
    campaign_path: Path = typer.Argument(..., help="Declarative campaign YAML (shared protein + per-system components)."),
    out_dir: Optional[Path] = typer.Option(None, "--out", "-o", help="Root output directory (default: ./<campaign>_campaign)."),
    no_grompp: bool = typer.Option(False, "--no-grompp", help="Skip gmx grompp validation for every system."),
) -> None:
    """Build several systems that share a protein, from one campaign specification."""
    from core.system_spec import SystemSpecError
    from core.system_campaign import run_campaign

    if not campaign_path.exists():
        app.print(f"[red]Error:[/red] campaign file not found: {campaign_path}")
        raise typer.Exit(1)

    try:
        from core.system_campaign import expand_campaign
        name, specs, _ = expand_campaign(campaign_path)
    except SystemSpecError as exc:
        app.print(Panel(f"[red]Invalid campaign:[/red]\n\n{exc}",
                        title="❌ Campaign Failed", border_style="red"))
        raise typer.Exit(1)

    resolved_out = out_dir or Path.cwd() / f"{name}_campaign"
    app.print(Panel(
        f"[bold cyan]SimForge Campaign[/bold cyan]  [dim]{name}[/dim]  "
        f"{len(specs)} system(s)  →  [dim]{resolved_out}[/dim]",
        border_style="cyan", padding=(0, 2),
    ))

    with app.status("  Building systems..."):
        result = run_campaign(
            campaign_path, resolved_out, run_grompp=not no_grompp,
            command=["simforge", "campaign-build", str(campaign_path)],
        )

    table = Table(box=box.SIMPLE, show_header=True, header_style="bold dim")
    table.add_column("System")
    table.add_column("Atoms", justify="right")
    table.add_column("Validation")
    table.add_column("Status")
    for sys_name, o in result.outcomes.items():
        v = o.validation.get("readiness", "—")
        table.add_row(
            sys_name,
            str(getattr(o.assembly_result, "total_atom_count", "—")),
            v,
            _fmt("PASS" if o.success else "FAIL"),
        )
    for sys_name, err in result.errors.items():
        table.add_row(sys_name, "—", "—", f"[red]✗ {err[:40]}[/red]")
    app.print(table)

    app.print(f"\n  [bold]{result.succeeded}/{result.total}[/bold] systems built  "
              f"[dim]summary → {resolved_out}/campaign_summary.yaml[/dim]\n")
    raise typer.Exit(0 if result.success else 1)


# ═══════════════════════════════════════════════════════════════════════════════
# simforge inspect-run
# ═══════════════════════════════════════════════════════════════════════════════

def inspect_run_fn(
    run_dir: Path = typer.Argument(..., help="A directory containing provenance.json (a built system)."),
    json_out: bool = typer.Option(False, "--json", help="Print the raw provenance record."),
) -> None:
    """Summarize a built system's provenance / reproducibility manifest."""
    prov_path = run_dir / "provenance.json"
    if not prov_path.exists():
        app.print(f"[red]Error:[/red] no provenance.json in {run_dir}")
        raise typer.Exit(1)

    prov = json.loads(prov_path.read_text())
    if json_out:
        print(json.dumps(prov, indent=2))
        return

    env = prov.get("environment", {})
    git = env.get("git", {})
    app.print(Panel(
        f"[bold]Generated[/bold]   {prov.get('generated_utc')}\n"
        f"[bold]Command[/bold]     {' '.join(prov.get('command') or []) or '—'}\n"
        f"[bold]Spec[/bold]        {prov.get('spec_path') or '—'}\n"
        f"[bold]SimForge[/bold]    {env.get('simforge_version')}  "
        f"(git {str(git.get('commit'))[:10]}{' +dirty' if git.get('dirty') else ''})\n"
        f"[bold]Python[/bold]      {env.get('python_version')}  [dim]{env.get('python_executable')}[/dim]\n"
        f"[bold]GROMACS[/bold]     {env.get('gromacs_version') or 'not found'}\n"
        f"[bold]RDKit[/bold]       {env.get('rdkit_version') or 'not found'}\n"
        f"[bold]Force field[/bold] {prov.get('forcefield')}   [bold]Water[/bold] {prov.get('water_model')}\n"
        f"[bold]Protein[/bold]     {prov.get('protein_mode')}  "
        f"(pdb2gmx {'used' if prov.get('pdb2gmx_used') else 'skipped'})",
        title=f"Reproducibility Manifest — {run_dir}", border_style="cyan", padding=(0, 2),
    ))

    if prov.get("discovery_notes"):
        app.print("  [bold]Auto-discovery[/bold]")
        for n in prov["discovery_notes"]:
            app.print(f"    [green]↪[/green] {n}")

    it = Table(box=box.SIMPLE, show_header=True, header_style="bold dim", title="Inputs")
    it.add_column("Role")
    it.add_column("File")
    it.add_column("SHA256", style="dim")
    for label, meta in (prov.get("inputs") or {}).items():
        it.add_row(label, Path(meta.get("path", "")).name, str(meta.get("sha256"))[:16])
    app.print(it)

    ot = Table(box=box.SIMPLE, show_header=True, header_style="bold dim", title="Outputs")
    ot.add_column("Role")
    ot.add_column("File")
    ot.add_column("SHA256", style="dim")
    for label, meta in (prov.get("outputs") or {}).items():
        ot.add_row(label, Path(meta.get("path", "")).name, str(meta.get("sha256"))[:16])
    app.print(ot)

    renames = prov.get("atomtype_renames") or {}
    if renames:
        n = sum(len(v) for v in renames.values())
        app.print(f"\n  [bold]Atomtype renames[/bold]: {n} across {len(renames)} component(s)")

    v = prov.get("validation") or {}
    vc = {"READY": "green", "READY (with warnings)": "yellow"}.get(v.get("readiness"), "red")
    app.print(f"  [bold]Validation[/bold]: [{vc}]{v.get('readiness', 'unknown')}[/{vc}]\n")
