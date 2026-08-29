# SimForge — Study Campaign Analysis

`simforge study inspect | analyze | analyses` is a **trajectory-first** layer for
turning a directory tree of finished MD runs into validated, per-system
observables with full provenance.

It discovers systems from the **trajectories and topologies** it finds on disk —
not from a filename convention. File names are treated as the weakest possible
evidence and never define scientific identity.

```
simforge study inspect   <dir>                     discover + validate, run NOTHING
simforge study analyses                            list registered analyses
simforge study analyze   <dir> --analysis <id> …   run ONLY the analyses you name
```

Related commands: `simforge build` / `doctor` / `validate-system` / `inspect-run`
(see `docs/system_build.md`).

---

## 1. Overview

```
Study            a directory tree of MD artifacts
└── Condition    one combination of every condition dimension (receptor, partner,
    │            water model, plus any extra axis found in the path)
    └── Replicate   one independent simulation of that condition (rep01, rep02, …)
        └── System  the concrete files: trajectory(s) + topology + structure + index
```

A **System** is one usable trajectory/topology pair. **Conditions** group
replicate systems that share every dimension. Discovery assigns each system a
deterministic, filesystem-safe **canonical id** built only from resolved
identity, e.g. `glp1r__semaglutide__tip3p__rep02` (unknown fields become the
literal token `unknown`; collisions get a stable content-hash suffix). Trajectory
and topology filenames never appear in it.

### Trajectory-first vs. the legacy XVG mode

`simforge study <xvgdir>` (no sub-command) is the **legacy comparative XVG
mode** — it still exists and is unchanged. It scans for pre-computed
`SYSTEM-REPLICAobservable.xvg` files (e.g. `AA-A1rmsd_protein.xvg`) and builds a
comparative report from them. It does not read trajectories.

The `inspect` / `analyze` / `analyses` sub-commands are the new layer: they read
trajectories, infer molecular components, build safe derived trajectories and run
GROMACS observables themselves.

---

## 2. `simforge study inspect <dir>`

Discovery → manifest → component inference → trajectory inspection → validation.
**Runs no scientific analysis.**

```
simforge study inspect ./glp1r_study
simforge study inspect ./glp1r_study -o /scratch/glp1r_analysis
simforge study inspect ./glp1r_study --no-trajectory-inspection
simforge study inspect ./glp1r_study --json
```

| Option | Effect |
|---|---|
| `-o, --output DIR` | Output directory (default `<study>/simforge_analysis`). |
| `--manifest FILE` | Load an existing `study_manifest.yaml` instead of rediscovering. |
| `--no-trajectory-inspection` | Skip `gmx check` (faster; no frame/box/atom-count metadata). |
| `--strong-fingerprint` | Full-content SHA-256 of every source file (slow for multi-GB trajectories; default is a fast size+mtime+edge hash). |
| `--gmx BIN` | GROMACS binary (default `gmx`). |
| `--json` | Machine-readable output. |

What it does:

* **Discovers** candidate systems (one per directory containing `.xtc`/`.trr`).
* **Associates** each with the nearest topology (`.tpr` preferred, then `.top`),
  reference structure (`.gro`/`.pdb`) and optional user `.ndx`, recording the
  evidence for every choice.
* **Infers molecular components** (receptor / peptide / membrane / water / ions …).
* **Inspects trajectories** with `gmx check`: frame count, timestep, box
  presence, precision, trajectory-vs-topology atom-count match, segment
  ordering/overlap.
* **Validates** at study / system / component level and writes a report.

Exit code is non-zero only when the study as a whole is `INVALID` (no systems, or
every system unusable).

---

## 3. `simforge study analyses`

Lists every registered analysis: id, category, display name, description and the
components it requires.

```
simforge study analyses
simforge study analyses --json
```

Phase 1 registers exactly the four RMSD observables (§8).

---

## 4. `simforge study analyze <dir> --analysis <id> …`

Runs the **same** discovery + validation pipeline as `inspect`, then executes
**only the analyses you explicitly select**, per system.

```
# explicit, repeatable --analysis (recommended for scripts)
simforge study analyze ./glp1r_study \
    --analysis rmsd-receptor --analysis rmsd-peptide-receptor-frame

# see the exact gmx commands without running them
simforge study analyze ./glp1r_study --analysis rmsd-complex --dry-run

# re-run after hand-correcting a classification
simforge study analyze ./glp1r_study --manifest ./glp1r_study/simforge_analysis/study_manifest.yaml \
    --analysis rmsd-peptide-intrinsic
```

| Option | Effect |
|---|---|
| `--analysis ID` | Analysis to run. **Repeatable.** Nothing runs unless you pass this (or pick interactively). |
| `-o, --output DIR` | Output directory (default `<study>/simforge_analysis`). |
| `--manifest FILE` | Use an existing (optionally hand-corrected) manifest. |
| `--non-interactive` | Never prompt; require `--analysis`. |
| `--dry-run` | Discover, validate and construct the GROMACS analysis commands, but do not execute them (status `planned`). Trajectory preprocessing (`gmx trjconv`) and the observable itself (`gmx rms`) are **not** run. The structure-only semantic-index build (`gmx select`) still runs so selection resolution can be validated. |
| `--force` | Ignore cached trajectory views / results. |
| `--gmx BIN` | GROMACS binary. |
| `--json` | Machine-readable output. |

**Selection rules:**

* **Nothing runs by default** — there is no implicit "run all".
* TTY and no `--analysis` → an **interactive menu**; enter numbers (`1,3`),
  `all`, or blank to cancel.
* Non-interactive (no TTY, `--non-interactive`, or `--json`) and no `--analysis`
  → prints guidance and exits `2`, running nothing.
* An unknown analysis id is an error (`simforge study analyses` lists valid ids).

Per system, each requested analysis ends in one status:

| Status | Meaning |
|---|---|
| `success` | Ran; `rmsd.xvg` produced. |
| `planned` | `--dry-run`: command constructed, not executed. |
| `cached` | Reused a prior result. |
| `skipped` | Requirement genuinely absent for this system (no such component, unusable system, or no safe trajectory view). |
| `review_required` | A component the analysis needs is **ambiguous** — resolve it in the manifest and re-run. |
| `failed` | GROMACS/tool error during execution. |
| `error` | Invalid request (unknown id, bad parameters). |

Analyses that need an unresolved distinction are **never run silently** — they
come back `skipped` or `review_required`, per system.

---

## 5. Discovery & evidence

**File names never define scientific identity.** Directory structure is recorded
as *evidence*, not treated as truth. Association uses progressively weaker
evidence:

1. **SimForge run metadata** — `metadata.json` / `provenance.json` /
   `simforge_run.json` in the sim dir or an ancestor (strongest).
2. **Same-directory** topology / structure / index files.
3. **Ancestor-directory** files (`ancestor_directory(+N)`), nearest first.
4. Stem preferences as a last tie-breaker (`md` / `prod` / `production` for
   `.tpr`, `topol` / `system` for `.top`, etc.).

Each decision is stored in `discovery_evidence` / `association_evidence`. When
two same-directory `.tpr` files exist, one is picked by the stem rule and an
`ambiguous_topology` warning is emitted.

### Trajectory-stage classification

Every trajectory in a system directory is classified into an MD-workflow
**stage** and only `production` trajectories are analysed:

| stage | how it is recognised |
|---|---|
| `production` | stem `md` / `prod` / `production` / `run`; paired `.mdp` with a long run (≥ 2 ns) or `part`/`seg` continuation markers; the coarsest output spacing / longest duration in the directory; or *by elimination* when it is the only trajectory and nothing points to equilibration |
| `equilibration_nvt` | stem `nvt`; paired `.mdp` with `gen_vel = yes` and no pressure coupling |
| `equilibration_npt` | stem `npt`; paired `.mdp` with `gen_vel = yes` + pressure coupling, or a short run with pressure coupling next to a much longer sibling |
| `equilibration_other` | stem `heat` / `anneal` / `equil` / `relax` / `posre` |
| `minimization` | stem `em` / `min`; `.mdp` `integrator = steep`/`cg` |
| `preparation` | stem `ions` / `solv` / `genion` |
| `derived` | stem contains `nojump` / `nopbc` / `whole` / `center` / `fit` / `aligned` / `wrapped` — a PBC/fit derivative of another trajectory |
| `unknown` | no evidence either way |

Evidence combined: filename stem (never decisive on its own), the paired `.mdp`
(`mdout.mdp` / `<stem>.mdp` — `gen_vel`, `pcoupl`, `nsteps × dt`, output
spacing), and directory-relative frame spacing / duration. A `.mdp` that says
`gen_vel = yes` overrides a stem of `md`.

`SystemRecord` records `production_trajectory_paths` (analysed) and, separately,
`trajectory_artifacts` (**every** trajectory with its `stage`,
`stage_confidence` and `stage_evidence`) and `workflow_trajectory_paths`
(equilibration / minimisation). Production frame count and duration are derived
**only** from the production trajectory(s) — equilibration frames are never
summed in.

* **Production continuation segments** — multiple `production` trajectories in
  one directory are treated as **one ordered collection only** when *all* carry
  `part`/`seg`/`chunk` markers. Otherwise an `ambiguous_production_trajectories`
  (`review_required`) warning is raised and they are **not** concatenated
  (same-directory membership alone is not evidence of one run).
* **Stage-matched topology / structure** — the production trajectory is paired
  with a stem-matching `.tpr` / `.gro` (`md.xtc ↔ md.tpr ↔ md.gro`). Only
  multiple *production-stage* topologies with no stem match trigger
  `ambiguous_production_topology` (`review_required`); a clean stem match is
  never flagged.
* **Equilibration-only directory** → no production trajectory →
  `no_production_trajectory` warning, system `invalid` for observables.
* **Replicate inference** — a path segment matching `rep`/`replica`/`r`/
  `run`/`sim`/`seed` + number → `repNN`. A **bare integer** directory
  (`300`, `310` — a temperature/salt axis) is only treated as the replicate
  when there is *no* explicit `rep` marker elsewhere in the path; otherwise it
  becomes a condition dimension. With no marker at all, `rep01` is assumed
  (`no_replicate_marker`, INFO — a documented harmless assumption).
* **Water-model inference** — a path segment matching a known token (`tip3p`,
  `spce`, `tip4p`, `opc`, …) sets `water_model` as a condition dimension.

Remaining path segments become the `partner` hint and extra `dimension_N` axes
(directory semantics — not structurally confirmed).

---

## 6. The StudyManifest

Written to the output dir as **both** `study_manifest.yaml` and
`study_manifest.json` (identical content). It records every discovery decision:
`study_root`, `generated_utc`, `simforge_version`, `schema_version`,
`discovery_settings`, `condition_dimensions`, `conditions[]`, `systems[]`
(files, fingerprints, components, semantic index, trajectory inspection,
evidence, warnings), `unassigned_files`, `ambiguities[]`, `warnings[]`, plus a
`stats` block.

It is a **reproducible contract**: pass it back with `--manifest` and the run
skips rescanning entirely — discovery, association, component inference and
trajectory metadata are all read from the file.

### Hand-correcting an ambiguous classification

When component inference cannot separate, say, receptor from peptide, the
manifest lists it under `ambiguities[]`:

```yaml
ambiguities:
  - system_id: glp1r__unknown__tip3p__rep01
    kind: component_classification
    message: "protein chains (unresolved): chains ['A', 'B'] with residue counts [463, 31]
              — no membrane anchor and no dominant signal to separate receptor from partner"
    options: ["chain A", "chain B"]
    resolution: null          # <- you fill this in
```

Set `resolution` to a `;`-separated list of `component=chain` assignments, then
re-run with `--manifest`:

```yaml
    resolution: "receptor=A;peptide=B"
```

```
simforge study analyze ./glp1r_study \
    --manifest ./glp1r_study/simforge_analysis/study_manifest.yaml \
    --analysis rmsd-peptide-receptor-frame
```

On load, the named components are set to those chains, their
`classification_state` becomes `resolved`, and the system is marked
`user_overridden: true`.

---

## 7. Molecular components & ambiguity

Each system carries a list of `components`, every one with **evidence**
(weighted `ComponentEvidence` entries), a `confidence` score and a
`classification_state`:

| State | Meaning |
|---|---|
| `resolved` | Inferred with sufficient confidence. |
| `ambiguous` | Multiple interpretations are compatible. |
| `unknown` | No evidence at all. |
| `review_required` | Ambiguous **and** it blocks a requested analysis. |

Inference combines whatever is cheaply available: chain list and residue
sequences from the reference structure, the `.top` `[ molecules ]` section, and
the membrane-embedding signal (does a chain's backbone span the bilayer?).

> **No unconditional size rule.** "Largest chain = receptor" is at most one weak
> tie-breaker (weight 0.10) and never decides on its own. A clear
> transmembrane-embedded chain plus a short non-TM chain resolves cleanly; a
> soluble multi-chain complex with only a size difference stays `ambiguous`;
> genuinely inseparable chains become `review_required`.

Water, ions, membrane and cofactors are resolved from known residue-name sets.
Sequence identity / database lookups are **out of scope** for Phase 1 (no
network dependency).

An analysis whose `required_components` are not all `resolved` for a system is
**not run**: `skipped` if the component is simply absent, `review_required` if
it is ambiguous (with remediation text pointing at the manifest).

---

## 8. Analyses implemented in Phase 1 — the four RMSD observables

All four operate on **backbone** atoms, emit RMSD in **nm** vs. time in **ps**,
and record their fit selection, measurement selection and reference frame in the
result and provenance. Selections are explicit **named** semantic-index groups
(`Receptor_Backbone`, `Peptide_Backbone`, `Complex_Backbone`) — never numeric
GROMACS group ids.

**Reference coordinates** for all four = the production **`.tpr`** (the
start-of-production frame `grompp` wrote from the equilibrated system — the
standard "RMSD relative to the production start"), falling back to a
`.gro`/`.pdb` only when no `.tpr` is available. The chosen reference filename is
recorded in `reference_frame` in every result and provenance file.

| id | category | requires | fit selection | measured | reference |
|---|---|---|---|---|---|
| `rmsd-receptor` | structural | receptor | receptor backbone | receptor backbone | production `.tpr` |
| `rmsd-complex` | structural | complex | complex backbone | complex backbone | production `.tpr` |
| `rmsd-peptide-intrinsic` | structural | peptide | peptide backbone | peptide backbone | production `.tpr` |
| `rmsd-peptide-receptor-frame` | interaction | receptor + peptide | receptor backbone (trajectory pre-fitted) | peptide backbone, **no further fit** (`gmx rms -nofit`) | production `.tpr` (trajectory pre-fitted to its receptor backbone) |

The first three run `gmx rms` on a `whole` trajectory view and let `gmx rms`
perform the least-squares fit to the measurement group. The last runs on a
**receptor-fitted** view and measures with `-nofit`.

### `rmsd-peptide-intrinsic` vs. `rmsd-peptide-receptor-frame` — the key distinction

* **`rmsd-peptide-intrinsic`** — the peptide is fitted **to itself** (its own
  backbone, against the reference `.tpr` coordinates). All rigid-body motion of
  the peptide is removed. This is the peptide's **own internal conformational
  drift** — is it folding, unfolding, fraying? — with receptor motion
  completely irrelevant.

* **`rmsd-peptide-receptor-frame`** — the trajectory is fitted to the
  **receptor** backbone, then the peptide backbone RMSD is measured with **no
  further fitting**. The peptide is *not* re-aligned to itself, so this captures
  **how the peptide moves relative to the receptor**: translation, rocking and
  drift out of the binding site all show up. This is a **binding-pose stability**
  metric. A peptide that stays rigidly bound but in a shifted pose has low
  `intrinsic` RMSD and high `receptor-frame` RMSD.

---

## 9. Trajectory views & PBC policy

Each observable declares **typed** `TrajectoryRequirements` (whole molecules?
nojump? centering target? fit selection?). The orchestration layer resolves them
into the **minimal safe derived trajectory**, built with `gmx trjconv` under the
output dir — applying only the steps that observable asked for. There is **no
universal `-pbc mol` recipe.**

| View | Built by | Meaning |
|---|---|---|
| `raw` | — | single source trajectory, untouched |
| `whole` | `-pbc whole` | molecules made whole across PBC |
| `nojump` | `-pbc nojump` | unwrapped, diffusion-safe |
| `centered` | `-pbc mol -center` | a named semantic group centred in the box |
| `fitted` | `-fit rot+trans` | least-squares fit to a named group (reference = frame 0) |

* **Originals are never modified.** Derived trajectories are cached by a key
  from source fingerprints + the requirement token; observables share a view
  only when it is scientifically identical.
* Every `gmx trjconv` step is fed a **named** index group (never a numeric id)
  and recorded as a `PreprocessingOperation` (command, stdin, reason, return
  code, output tails).
* **If a safe view cannot be built** — a required semantic group is missing, a
  step fails, or the trajectory has multiple un-concatenated segments — the view
  is `safe = false` and the analysis **fails / is skipped**. It never falls back
  to raw coordinates.

---

## 10. Output layout

```
simforge_analysis/                     (or -o <dir>)
├── study_manifest.yaml
├── study_manifest.json
├── validation_report.json
├── validation_report.txt
├── analysis_results.json              (analyze only)
└── systems/
    └── <canonical_id>/
        ├── indices/
        │   └── semantic_index.ndx      generated; the user's .ndx is never overwritten
        ├── trajectories/
        │   └── <cache_key>/<view>.xtc  derived trajectory views (whole/nojump/…)
        └── observables/
            └── <analysis_id>/
                ├── rmsd.xvg
                ├── gmx_rms.log         exact argv + stdin + stdout/stderr
                └── provenance.json
```

---

## 11. Provenance

Every observable writes `provenance.json` (schema
`simforge/study-campaign/observable-provenance/v1`) recording:

* `timestamp`, `simforge_version`, the exact `argv`, and the collected
  `environment` (Python / platform / conda env, git commit + branch + dirty
  flag, GROMACS / RDKit / OpenBabel versions).
* `system_id` / `condition_id` / `replicate_id`, `analysis_id`, `parameters`.
* Semantic selections: `fit_selection`, `measure_selection`, `reference_frame`.
* `source_files` + `source_fingerprints` (path, size, mtime, digest, mode) for
  every trajectory / topology / structure / index.
* `semantic_index` (named groups, atom/residue counts, selection logic).
* `trajectory_view` — the full chain of preprocessing operations with commands.
* `component_evidence` — each component's type, state, confidence, chains and
  weighted evidence.
* `trajectory_inspection` — frames, timestep, box, atom-count match, segments.
* `result` — status, output files, data summary (mean / min / max / final /
  plateau-mean RMSD, frame count, GROMACS version).

The manifest plus these files make each result reproducible without rescanning.

---

## 12. Non-interactive / automation usage

```
# discovery + validation, JSON, no gmx check
simforge study inspect ./study --json --no-trajectory-inspection > inspect.json

# strict batch run — fails loudly rather than prompting
simforge study analyze ./study --non-interactive --json \
    --analysis rmsd-receptor --analysis rmsd-complex > results.json
```

* `inspect` exits non-zero when the study is `INVALID`; `analyze` exits non-zero
  when a requested run produced no successful/planned/cached result.
* `--non-interactive` with no `--analysis` exits `2` and runs nothing.
* `analysis_results.json` and the manifest/validation files are the stable
  machine artifacts.

---

## 13. Current limitations (Phase 1)

* **Only RMSD observables** are implemented — the four in §8. No RMSF, gyration,
  contacts, distances, SASA, hydrogen bonds, etc. yet.
* **No replicate aggregation.** Results are per system. Averaging across the
  replicates of a condition is Phase 2.
* **Production continuation segments are not auto-concatenated.** Multiple
  `production` trajectories (even with `part` markers) are recorded and ordered
  but not merged — the observable's `whole` view build fails safe and the
  analysis is `skipped`. Provide a single trajectory or a validated
  concatenation.
* **Trajectory-stage classification is heuristic.** It combines filename stem,
  the paired `.mdp` and directory-relative timing. It is deliberately
  conservative (a stage it cannot confirm as production-by-elimination gets a
  `production_stage_unverified` note). A production run with an mdp missing and
  an equilibration-suggesting name would need a manifest correction.
* **Component inference has no sequence / database step** — chains, residue
  classes, the `.top` molecule table and membrane embedding only. A single
  polymer chain is labelled `primary protein` with `study_role =
  "receptor (assumed)"` — it is used for receptor-based analyses but is **not**
  a verified biological receptor identity.
* **The canonical `receptor` slot is filled from the study-root directory name**
  — a **study label**, not a verified molecular identity (recorded as such in
  `discovery_evidence`). If the root is generically named (`data`, `runs`,
  `production`, …) it is left `unknown`.
* **`gmx` is required** for trajectory inspection, semantic-index generation and
  every observable. Without it, inspection is skipped and analyses fail /
  `review_required`.

---

## 14. No comparison, statistics or figures yet

This layer discovers, validates and computes **per-system observables** with
provenance. Cross-condition comparison, statistical testing and figures are
Phase 2 / Phase 3 and are not produced. For a comparative report today, the
legacy `simforge study <xvgdir>` XVG mode still applies.
