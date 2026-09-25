"""``simforge trajectory review`` / ``review-observables`` — thin entry points.

Headless: preparation writes a ReviewDataset and exits.  No viewer, browser
or server is started (those are later layers).
"""
from __future__ import annotations

import json
import tempfile
from pathlib import Path
from typing import Annotated, Optional

import typer
from rich.console import Console

_console = Console()

_SHOW_HELP = ("Observable instance(s), repeatable / comma-separated: an id "
              "(rmsd-receptor) or id(key=value,...) e.g. "
              "'rg(selection=component:receptor)', "
              "'com-distance(selection_a=component:ligand,selection_b=annotation:site)'.")


def _requests(show, request_file):
    from analysis.review.request import RequestError, load_request_file, parse_show
    try:
        reqs = parse_show(show or [])
        if request_file:
            reqs += [r for r in load_request_file(request_file) if r not in reqs]
    except (RequestError, OSError) as exc:
        _console.print(f"[red]Error:[/red] {exc}")
        raise typer.Exit(2)
    return reqs


def review_fn(
    target: Annotated[Path, typer.Argument(help="RUN_DIR (campaign discovery) or a TRAJECTORY file.")],
    topology: Annotated[Optional[Path], typer.Option("--topology", "-s",
        help="TRAJECTORY mode: run input / topology (.tpr recommended).")] = None,
    structure: Annotated[Optional[Path], typer.Option("--structure",
        help="TRAJECTORY mode: reference structure (.gro/.pdb).")] = None,
    index: Annotated[Optional[Path], typer.Option("--index", help="TRAJECTORY mode: user .ndx.")] = None,
    manifest: Annotated[Optional[str], typer.Option("--manifest",
        help="Use an existing (hand-corrected) study_manifest.yaml.")] = None,
    system: Annotated[Optional[str], typer.Option("--system",
        help="System id when the run directory holds several.")] = None,
    show: Annotated[Optional[list[str]], typer.Option("--show", help=_SHOW_HELP)] = None,
    request_file: Annotated[Optional[Path], typer.Option("--request",
        help="YAML/JSON file with an 'observables' list (observable + parameters).")] = None,
    display: Annotated[Optional[str], typer.Option("--display",
        help="Display view: raw (default; source coordinates) | whole | center:<Group> | "
             "fit:<Group>. Judged by the preprocessing policy for purpose 'display'. "
             "An explicit choice replaces a profile's preference.")] = None,
    profile: Annotated[Optional[str], typer.Option("--profile",
        help="Review profile: a built-in id (see `trajectory review-profiles`) or a YAML "
             "file. Its requests are merged with --show/--request.")] = None,
    hide: Annotated[Optional[list[str]], typer.Option("--hide",
        help="Remove a profile request (instance id or observable id; repeatable).")] = None,
    display_intent: Annotated[bool, typer.Option("--display-intent",
        help="Explicitly choose an interpretation-changing display transformation "
             "(recorded as user_flag intent; applies to the display view only).")] = False,
    output: Annotated[Optional[str], typer.Option("--output", "-o",
        help="Analysis root (default: <RUN_DIR>/simforge_analysis, or ./simforge_analysis "
             "in TRAJECTORY mode). Sessions go to <root>/review/.")] = None,
    reuse: Annotated[bool, typer.Option("--reuse/--no-reuse",
        help="Reuse scientifically compatible existing results (Phase 8).")] = True,
    reuse_from: Annotated[Optional[list[Path]], typer.Option("--reuse-from",
        help="Additional campaign output root(s) to search for reusable results.")] = None,
    cache_dir: Annotated[Optional[Path], typer.Option("--cache-dir",
        help="Diagnostics probe / time-index cache (default: <root>/review/cache).")] = None,
    prepare_only: Annotated[bool, typer.Option("--prepare-only",
        help="Prepare the session and exit (the only mode in this version).")] = True,
    dry_run: Annotated[bool, typer.Option("--dry-run",
        help="Resolve system, annotations, display policy, applicability and expected "
             "reuse; run no GROMACS observable and build no derived trajectory.")] = False,
    strict: Annotated[bool, typer.Option("--strict",
        help="Exit non-zero if any requested observable or the display view is not available.")] = False,
    force: Annotated[bool, typer.Option("--force",
        help="Re-prepare even if an identical valid session exists (results may still be reused).")] = False,
    gmx: Annotated[str, typer.Option("--gmx", help="GROMACS binary.")] = "gmx",
    as_json: Annotated[bool, typer.Option("--json", help="Emit the ReviewDataset as JSON.")] = False,
    serve: Annotated[bool, typer.Option("--serve",
        help="After preparing, serve the session on a local dashboard (see `trajectory serve`).")] = False,
    no_browser: Annotated[bool, typer.Option("--no-browser",
        help="With --serve: do not open a browser; just print the URL.")] = False,
    viewer: Annotated[str, typer.Option("--viewer",
        help="With --serve: viewer (null | fake | vmd).")] = "null",
    vmd: Annotated[Optional[str], typer.Option("--vmd",
        help="VMD executable (default: `vmd` on PATH).")] = None,
    vmd_headless: Annotated[bool, typer.Option("--vmd-headless",
        help="Run VMD in text mode (no graphics; testing / development).")] = False,
) -> None:
    """Prepare a reproducible trajectory review session (headless by default)."""
    from analysis.review.prepare import ReviewError, prepare_review
    from analysis.review.request import DisplayRequest, RequestError, ReviewRequest
    if serve and (dry_run or as_json):
        _console.print("[red]Error:[/red] --serve needs a real preparation (no --dry-run / --json)")
        raise typer.Exit(2)
    if not prepare_only and not serve:
        _console.print("[red]Error:[/red] nothing to do: use --prepare-only (default) or --serve")
        raise typer.Exit(2)
    reqs = _requests(show, request_file)
    from analysis.review.profiles import ProfileError, get_profile
    try:
        prof = get_profile(profile) if profile else None
        if display is not None:
            disp = DisplayRequest.parse(display, intent=display_intent)
        elif display_intent and prof is not None and prof.display:
            # the user explicitly authorises the profile's display preference
            disp = DisplayRequest.parse(prof.display, intent=True)
        elif display_intent:
            disp = DisplayRequest.parse("raw", intent=True)
        else:
            disp = None
        req = ReviewRequest(reqs, disp)
        prep = prepare_review(target, req, topology=topology, structure=structure, index=index,
                              manifest=manifest, system=system, output=output, reuse=reuse,
                              reuse_from=tuple(reuse_from or ()), cache_dir=cache_dir,
                              dry_run=dry_run, force=force, gmx=gmx, profile=prof,
                              hide=tuple(hide or ()))
    except (ReviewError, RequestError, ProfileError) as exc:
        _console.print(f"[red]Error:[/red] {exc}")
        raise typer.Exit(2)
    failures = prep.strict_failures()
    if as_json:
        print(json.dumps({"preparation": {
            "dataset_path": str(prep.dataset_path) if prep.dataset_path else None,
            "reused_session": prep.reused_session, "dry_run": prep.dry_run,
            "replaced": prep.replaced, "strict_failures": failures},
            "review_dataset": prep.dataset.to_dict()}, indent=2))
    else:
        render_preparation(prep)
    if strict and failures:
        raise typer.Exit(1)
    if serve:
        _serve(prep.dataset_path, no_browser=no_browser, viewer=viewer, vmd=vmd,
               vmd_headless=vmd_headless)
    raise typer.Exit(0)


def _serve(dataset_path, *, no_browser: bool, viewer: str, host: str = "127.0.0.1",
           port: int = 0, vmd: Optional[str] = None, vmd_headless: bool = False) -> None:
    from analysis.review.runtime.server import ReviewServeError, ReviewServer, serve_forever
    from analysis.review.runtime.vmd import VMDError
    opts = {"vmd": vmd, "headless": vmd_headless} if viewer == "vmd" else {}
    if viewer == "vmd":
        from analysis.review.runtime.vmd import find_vmd
        try:
            find_vmd(vmd)                        # actionable error before anything starts
        except VMDError as exc:
            _console.print(f"[red]Error:[/red] {exc}")
            raise typer.Exit(1)
    try:
        server = ReviewServer(dataset_path, host=host, port=port, viewer=viewer,
                              viewer_options=opts)
    except (ReviewServeError, VMDError, ValueError) as exc:
        _console.print(f"[red]Refusing to serve:[/red] {exc}")
        raise typer.Exit(1)
    serve_forever(server, open_browser=not no_browser,
                  echo=lambda m: print(m, flush=True))      # the URL must reach pipes now


def serve_fn(
    dataset: Annotated[Path, typer.Argument(help="review_dataset.json or its session directory.")],
    no_browser: Annotated[bool, typer.Option("--no-browser",
        help="Do not open a browser; print the URL only (headless use).")] = False,
    host: Annotated[str, typer.Option("--host",
        help="Bind address (default loopback only).")] = "127.0.0.1",
    port: Annotated[int, typer.Option("--port", help="Port (0 = any free port).")] = 0,
    viewer: Annotated[str, typer.Option("--viewer",
        help="Viewer: null (default) | fake | vmd (loads the session's display trajectory).")] = "null",
    vmd: Annotated[Optional[str], typer.Option("--vmd",
        help="VMD executable (default: `vmd` on PATH; nothing else is searched).")] = None,
    vmd_headless: Annotated[bool, typer.Option("--vmd-headless",
        help="Run VMD in text mode (no graphics; testing / development).")] = False,
) -> None:
    """Serve an already prepared, valid review session on a local dashboard.

    Presentation only: nothing is discovered, analysed or re-prepared; a stale
    or invalid session is refused."""
    if host not in ("127.0.0.1", "localhost", "::1"):
        _console.print(f"[yellow]warning:[/yellow] binding to {host} exposes the session "
                       f"beyond this machine (token still required)")
    _serve(dataset, no_browser=no_browser, viewer=viewer, host=host, port=port, vmd=vmd,
           vmd_headless=vmd_headless)


_STATE_STYLE = {"available": "green", "planned": "cyan", "not_applicable": "dim",
                "review_required": "yellow", "ambiguous": "yellow", "unsupported": "red", "refused": "red",
                "failed": "red"}


def render_preparation(prep) -> None:
    ds = prep.dataset
    s = ds.summary()
    tl = ds.timeline
    head = "Review session (dry run — nothing written)" if prep.dry_run else "Review session"
    _console.print(f"\n[bold cyan]{head}[/bold cyan]  [dim]{ds.session_id}[/dim]")
    if prep.reused_session:
        _console.print("[green]identical valid session already prepared — reused as is[/green]")
    for why in prep.replaced:
        _console.print(f"[yellow]previous session re-prepared:[/yellow] {why}")
    comps = ", ".join(f"{c['type']}" for c in ds.system.get("components", [])
                      if c["state"] == "resolved")
    _console.print(f"System      {ds.system.get('system_id')}  [dim]({comps})[/dim]")
    src = (ds.sources.get("trajectory") or {}).get("path")
    _console.print(f"Trajectory  {src}")
    if tl.get("n_frames") is not None:
        dup = f", {tl['n_duplicate_groups']} duplicate-time groups" if tl.get("n_duplicate_groups") else ""
        _console.print(f"Frames      {tl['n_frames']}  ({tl['start_time_ps']} → "
                       f"{tl['end_time_ps']} ps, timeline {tl['state']}{dup})")
    else:
        _console.print(f"Frames      [yellow]{tl.get('state')}[/yellow] {tl.get('reason', '')}")
    _render_profile((ds.provenance or {}).get("profile"))
    dv = ds.display_view
    colour = _STATE_STYLE.get(dv.get("status"), "white")
    _console.print(f"Display     {(dv.get('request') or {}).get('mode')} → "
                   f"[{colour}]{dv.get('status')}[/{colour}] ({dv.get('kind')}, frames "
                   f"{(tl.get('display') or {}).get('relation', 'n/a')})")
    if dv.get("reason"):
        _console.print(f"            [dim]{dv['reason']}[/dim]")
    dg = ds.diagnostics
    counts = dg.get("counts") or {}
    _console.print(f"Diagnostics {dg.get('execution')}/{dg.get('status')}"
                   + (f"  warnings={counts.get('warnings', 0)} review={counts.get('review_required', 0)} "
                      f"errors={counts.get('errors', 0)}" if counts else f"  {dg.get('reason', '')}"))
    _console.print("Observables")
    for o in ds.observables:
        colour = _STATE_STYLE.get(o.state, "white")
        how = o.availability if o.state == "available" else o.state
        if o.state == "planned" and o.reuse:
            how = {"would_reuse": "planned: would reuse",
                   "would_recompute": "planned: would compute"}.get(o.reuse.get("decision"), how)
        elif o.state == "planned":
            how = "planned: would compute"
        sync = ""
        if o.results:
            sync = "  sync=" + ",".join("yes" if r.sync.get("syncable") else "no" for r in o.results)
        _console.print(f"  [{colour}]{how:<24}[/{colour}] {o.instance_id}{sync}")
        if o.state not in ("available", "planned") and o.reason:
            _console.print(f"      [dim]{o.reason[:300]}[/dim]")
    _console.print(f"Summary     computed={s['observables']['computed']} "
                   f"reused={s['observables']['cached']} blocked={s['observables']['blocked']} "
                   f"not_applicable={s['observables']['not_applicable']} "
                   f"failed={s['observables']['failed']}")
    if prep.dataset_path:
        _console.print(f"ReviewDataset {prep.dataset_path}")


def _render_profile(p) -> None:
    if not p:
        return
    colour = {"applicable": "green", "partially_applicable": "yellow"}.get(p["status"], "red")
    _console.print(f"Profile     {p['id']} v{p['version']} ({p['source']}) → "
                   f"[{colour}]{p['status']}[/{colour}]  [dim]{p['definition_identity'][:12]}[/dim]")
    for r in p["requirements"]:
        mark = "✓" if r["ok"] else "✗"
        _console.print(f"            {mark} {r['requirement']} {r['component']}"
                       + (f" — {r['reason']}" if r["reason"] else ""))
    _console.print(f"            included: {len(p['requests'])} request(s) from the profile")
    for o in p["optional"]:
        if o["outcome"] == "included":
            _console.print(f"            optional included: {o['instance_id']}")
        else:
            _console.print(f"            optional {o['outcome']}: {o['request']} — {o['reason']}")
    ov = p["overrides"]
    if ov["added"] or ov["hidden"]:
        _console.print(f"            user added {ov['added']} · hidden {ov['hidden']}")
    _console.print(f"            display from: {ov['display_source']}")


def review_profiles_fn(
    target: Annotated[Optional[Path], typer.Argument(
        help="Optional RUN_DIR / TRAJECTORY: also check each profile's requirements "
             "against that system (factual, unranked).")] = None,
    topology: Annotated[Optional[Path], typer.Option("--topology", "-s")] = None,
    structure: Annotated[Optional[Path], typer.Option("--structure")] = None,
    system: Annotated[Optional[str], typer.Option("--system")] = None,
    output: Annotated[Optional[str], typer.Option("--output", "-o")] = None,
    profile: Annotated[Optional[list[str]], typer.Option("--profile",
        help="Also list these profile files / ids.")] = None,
    gmx: Annotated[str, typer.Option("--gmx")] = "gmx",
    as_json: Annotated[bool, typer.Option("--json")] = False,
) -> None:
    """List review profiles; with a system, show whether each one's requirements hold."""
    from analysis.review.profiles import ProfileError, get_profile, list_profiles, resolve_profile
    try:
        profiles = list_profiles() + [get_profile(x) for x in (profile or [])]
    except ProfileError as exc:
        _console.print(f"[red]Error:[/red] {exc}")
        raise typer.Exit(2)
    rec = None
    if target is not None:
        from analysis.campaign.manifest import build_manifest
        from analysis.review.prepare import ReviewError, _select_system, resolve_inputs
        try:
            root, _, _ = resolve_inputs(target, topology=topology, structure=structure,
                                        output=output)
            rec = _select_system(build_manifest(root, inspect_trajectories=False, gmx=gmx), system)
        except ReviewError as exc:
            _console.print(f"[red]Error:[/red] {exc}")
            raise typer.Exit(2)
    rows = []
    for p in profiles:
        row = p.to_dict()
        if rec is not None:
            res = resolve_profile(p, rec)
            row["system"] = {"system_id": rec.system_id, "requirements_status": res.status,
                             "requirements": res.requirements, "optional": res.optional}
        rows.append(row)
    if as_json:
        print(json.dumps(rows, indent=2))
        raise typer.Exit(0)
    for r in rows:
        _console.print(f"\n[bold cyan]{r['id']}[/bold cyan] v{r['version']} "
                       f"[dim]({r['source']}; {r['definition_identity'][:12]})[/dim]")
        _console.print(f"  {' '.join(r['description'].split())}")
        req = r["requirements"]
        _console.print(f"  requires: {req['components'] or '—'}"
                       + (f"  absent: {req['absent_components']}" if req["absent_components"] else ""))
        _console.print(f"  display: {r['display'] or 'raw (default)'}"
                       + (f"  composition: {r['composition']}" if len(r["composition"]) > 1 else ""))
        for o in r["observables"]:
            _console.print(f"    • {o['observable']} {o['parameters'] or ''}")
        for o in r["optional_observables"]:
            _console.print(f"    ◦ if {o['when']}: {o['observable']} {o['parameters']}")
        if "system" in r:
            st = r["system"]["requirements_status"]
            text = {"applicable": "requirements satisfied",
                    "not_applicable": "requirements not satisfied",
                    "review_required": "requirements need review"}.get(st, st)
            why = "; ".join(x["reason"] for x in r["system"]["requirements"] if x["reason"])
            _console.print(f"  on {r['system']['system_id']}: {text}" + (f" — {why}" if why else ""))


# ═══════════════════════════════════════════════════════════════════════════════
# review-observables: read-only capability listing
# ═══════════════════════════════════════════════════════════════════════════════

def review_observables_fn(
    target: Annotated[Path, typer.Argument(help="RUN_DIR or TRAJECTORY file.")],
    topology: Annotated[Optional[Path], typer.Option("--topology", "-s")] = None,
    structure: Annotated[Optional[Path], typer.Option("--structure")] = None,
    index: Annotated[Optional[Path], typer.Option("--index")] = None,
    manifest: Annotated[Optional[str], typer.Option("--manifest")] = None,
    system: Annotated[Optional[str], typer.Option("--system")] = None,
    show: Annotated[Optional[list[str]], typer.Option("--show", help=_SHOW_HELP)] = None,
    output: Annotated[Optional[str], typer.Option("--output", "-o")] = None,
    gmx: Annotated[str, typer.Option("--gmx")] = "gmx",
    as_json: Annotated[bool, typer.Option("--json")] = False,
) -> None:
    """List observables, their outputs and applicability to this system — nothing is computed."""
    from analysis.review.prepare import ReviewError
    try:
        listing = list_capabilities(target, requests=_requests(show, None), topology=topology,
                                    structure=structure, index=index, manifest=manifest,
                                    system=system, output=output, gmx=gmx)
    except ReviewError as exc:
        _console.print(f"[red]Error:[/red] {exc}")
        raise typer.Exit(2)
    if as_json:
        print(json.dumps(listing, indent=2))
        raise typer.Exit(0)
    _console.print(f"\n[bold cyan]Observables for[/bold cyan] {listing['system_id']}\n")
    _console.print("[bold]Selectable references[/bold]")
    for r in listing["selections"]:
        colour = "green" if r["usable"] else "yellow"
        _console.print(f"  [{colour}]{r['ref']:<40}[/{colour}] {r['state']}"
                       + (f"  [dim]{r['reason']}[/dim]" if r.get("reason") else ""))
    _console.print("\n[bold]Observables[/bold]")
    for o in listing["observables"]:
        colour = {"yes": "green", "no": "yellow", "needs parameters": "cyan"}[o["applicable"]]
        _console.print(f"  [cyan]{o['id']:<28}[/cyan] {o['quantity']} [{o['unit']}] "
                       f"axes={o['axis_signature']}  applicable: [{colour}]{o['applicable']}[/{colour}]")
        if o["parameters"]:
            _console.print(f"      [dim]parameters: {o['parameters']}[/dim]")
        for why in o["reasons"]:
            _console.print(f"      [dim]{why}[/dim]")
    for i in listing["instances"]:
        colour = "green" if i["applicable"] else "yellow"
        _console.print(f"  [{colour}]{i['instance_id']}[/{colour}] applicable={i['applicable']}"
                       + "".join(f"\n      [dim]{w}[/dim]" for w in i["reasons"]))


def list_capabilities(target, *, requests=(), topology=None, structure=None, index=None,
                      manifest=None, system=None, output=None, gmx="gmx") -> dict:
    """Pure capability listing: discovery + annotations + a temporary semantic
    index; no diagnostics, no views, no observable execution."""
    import shutil
    from analysis.campaign.manifest import build_manifest, load_manifest
    from analysis.campaign.models import AnnotationState, ClassificationState
    from analysis.campaign.observables import registry
    from analysis.campaign.observables.selections import COMPONENT_GROUPS
    from analysis.campaign.orchestration.study_analyzer import prepare_system
    from analysis.review.prepare import _select_system, resolve_inputs
    from analysis.review.request import request_problems
    registry.ensure_loaded()
    root, _, _ = resolve_inputs(target, topology=topology, structure=structure, index=index,
                                output=output)
    mf = load_manifest(manifest) if manifest else build_manifest(
        root, inspect_trajectories=False, gmx=gmx)
    rec = _select_system(mf, system)
    tmp = Path(tempfile.mkdtemp(prefix="simforge-review-list-"))
    try:
        prepare_system(rec, tmp, gmx=gmx, diagnostics="off")
        selections = [{"ref": "component:system", "state": "resolved", "usable": True}]
        for c in rec.components:
            if c.component_type in COMPONENT_GROUPS:
                ok = c.classification_state == ClassificationState.RESOLVED
                selections.append({"ref": f"component:{c.component_type}",
                                   "state": c.classification_state, "usable": ok})
        for a in rec.annotations:
            selections.append({"ref": f"annotation:{a.annotation_id}", "state": a.state,
                               "usable": a.state == AnnotationState.ACTIVE and bool(a.group_name),
                               "reason": "; ".join(a.reasons)})
        observables = []
        for oid in registry.ids():
            spec = registry.get(oid)
            (arr, *_) = spec.output_schema({}) or [None]
            schema = spec.parameters_schema()
            if schema:
                appl, reasons = "needs parameters", []
            else:
                a = spec.applicability(rec, {}, rec.semantic_index)
                appl, reasons = ("yes" if a["applicable"] else "no"), a["reasons"]
            observables.append({
                "id": oid, "display_name": spec.display_name, "purpose": spec.purpose,
                "quantity": arr.quantity if arr else None, "unit": arr.unit if arr else None,
                "axis_signature": list(arr.axis_signature()) if arr else [],
                "parameters": schema, "applicable": appl, "reasons": reasons,
                "required_components": list(spec.required_components),
                "required_annotations": list(spec.required_annotations)})
        instances = []
        for r in requests:
            probs = request_problems(r)
            if probs:
                instances.append({"instance_id": r.instance_id, "applicable": False,
                                  "reasons": probs})
                continue
            a = registry.get(r.observable).applicability(rec, {r.observable: r.params},
                                                         rec.semantic_index)
            instances.append({"instance_id": r.instance_id, "applicable": a["applicable"],
                              "reasons": a["reasons"]})
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return {"system_id": rec.system_id, "selections": selections, "observables": observables,
            "instances": instances}
