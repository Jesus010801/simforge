# xanthone_short controlled execution — final report

- **Date:** 2026-09-08
- **Command:** `simforge study run /home/jesusxd/Escritorio/Nuevos_sistemas --profile xanthone_short --missing-only`
- **Protocol / profile:** `xanthone_short` / `xanthone_short/1.0`
- **GROMACS:** 2025.2 (`/usr/local/gromacs/bin/gmx`)
- **Analysis window:** `-b 0 -e 25000` ps (0–25 ns) on every analysis
- **Contact cutoff:** `-d 0.6` nm
- **Catalytic distance:** `COM(LIG)` vs `COM(Catalytic_*)` via `gmx distance -select`
- **Output root:** `/home/jesusxd/Escritorio/simforge/analysis_outputs/xanthone_short/` (outside the scientific source tree)

## 1. Systems executed

| System | Eligible | Basis | Production (measured) | Analysis window |
|---|---|---|---|---|
| AA-A6 | yes | compound A6 (0.90, structural), target AA (0.85, index headers) | 25 000 ps | 0–25 000 ps |
| AG-A6 | yes | compound A6 (0.90), target AG (0.85) | 25 000 ps | 0–25 000 ps |
| LP-A6 | yes | compound A6 (0.90), target LP (0.85) | 25 000 ps | 0–25 000 ps |
| HMG-R-25ns-A6 | yes — **override** | compound A6 (0.90), target HMG-R (0.85); folder claims ≤50 ns scope | **200 000 ps** | 0–25 000 ps |

HMG-R-25ns-A6 ran under the explicit-`xanthone_short` override: identity and all four
semantic groups are high-confidence, the folder name asserts a ≤50 ns scope, and the
only discrepancy — `folder duration conflicts with production metadata` (MDP = 200 ns) —
is recorded as `waived`. The 200 ns production trajectory was **analysed over 0–25 ns and
never truncated or rewritten**.

## 2. Systems blocked (REVIEW_REQUIRED) — with missing evidence

| System | Why blocked |
|---|---|
| HMG-R-200ns-A6 | compound identity is a folder hint only (confidence 0.45 — no `topol.top`/`A6.itp`); no production trajectory found (only derived `mdcenter.xtc` / `mdfit.xtc`); classified `mechanistic` |
| A3-HMG-R | compound is **A3**, not A6; classified `mechanistic` (intended 200 ns); folder does not claim a short scope |
| system_A3_COA | cofactor system (`COA` at index group 14) — requires a dedicated mechanistic mapping; biological target not resolvable from generic `Protein` label |
| system_A6_COA | cofactor system (`COA` at index group 14) — same as above |

`HMG-R-200ns-A6` was **not** used as a short-study production source: its provenance does
not establish that, so it stays `REVIEW_REQUIRED` as instructed.

## 3. Analyses executed / skipped / blocked

| Status | Count | Detail |
|---|---|---|
| EXECUTED | 20 | 5 analyses × 4 systems (protein RMSD, ligand RMSD, ligand–active-site min-dist, ligand–active-site contacts, ligand–catalytic-site COM distance) |
| REVIEW_REQUIRED | 20 | 5 analyses × 4 blocked systems |
| FAILED | 0 | — |
| SKIP_EXISTING | 0 (first run) → 20 (rerun) | rerun re-detected the managed results and did **not** recompute |

No optional analyses (RMSF, SASA, H-bonds, MM-PBSA, interaction energy, PCA, FEL,
clustering, DCCM) were run. Existing validated MM-PBSA summaries in the source
directories were left untouched.

## 4. Resolved semantic groups (dynamic, from each `index.ndx`)

| Role | AA-A6 | AG-A6 | LP-A6 | HMG-R-25ns-A6 |
|---|---|---|---|---|
| protein | `Protein` (1) | `Protein` (1) | `Protein` (1) | `Protein` (1) |
| ligand | `LIG` (13) | `LIG` (13) | `LIG` (13) | `LIG` (13) |
| active_site | `ActiveSite_AA` (21) | `ActiveSite_AG` (21) | `ActiveSite_LIP` (21) | `ActiveSite_HMG` (21) |
| catalytic_site | `Catalytic_AA` (22) | `Catalytic_AG` (22) | `Catalytic_LP` (22) | `Catalytic_HMG` (22) |

Groups are passed to GROMACS **by name**; the numeric ids above are recorded in
provenance only and are never fed to the tools.

## 5. GROMACS commands (AA-A6 shown; others identical modulo path / group name)

```
gmx rms      -s <dir>/md.tpr -f <dir>/md.xtc -n <dir>/index.ndx -b 0 -e 25000 \
             -o .../protein_rmsd/rmsd_protein.xvg              # stdin: "Protein\nProtein\n"

gmx rms      -s <dir>/md.tpr -f <dir>/md.xtc -n <dir>/index.ndx -b 0 -e 25000 \
             -o .../ligand_rmsd/rmsd_ligand.xvg                # stdin: "Protein\nLIG\n"

gmx mindist  -s <dir>/md.tpr -f <dir>/md.xtc -n <dir>/index.ndx -b 0 -e 25000 -d 0.6 \
             -od .../active_site/mindist_lig_active.xvg \
             -on .../active_site/contacts_lig_active.xvg       # stdin: "LIG\nActiveSite_AA\n"

gmx distance -s <dir>/md.tpr -f <dir>/md.xtc -n <dir>/index.ndx -b 0 -e 25000 \
             -select 'com of group "LIG" plus com of group "Catalytic_AA"' \
             -oall .../catalytic_com_distance/dist_lig_catalytic.xvg
```

`<dir>` = `/home/jesusxd/Escritorio/Nuevos_sistemas/AA-A6`. Full argv, stdin, and the
`gmx.log` transcript are stored per analysis.

## 6. Time coverage & point counts (post-run validation, spec §10)

Every one of the 20 generated XVGs:

- readable XVG structure, `@TYPE xy`, correct axis labels
- time range **0.0 → 25 000.0 ps**
- **2501 points** (10 ps sampling — matches the historical A1–A5 corpus exactly)
- all values finite

Header `@` directives are **byte-identical** to the historical thesis families
(`AA-A5rmsd_protein.xvg`, `AA-A2contacts_lig_active.xvg`, `AA-A3dist_lig_catalytic.xvg`):
same `@ subtitle "Protein after lsq fit to Protein"`, same
`@ title "Number of Contacts < 0.6 nm"`, same `@ s0 legend "LIG-ActiveSite_AA"`.
No numerical similarity to A1–A5 is claimed or required — only methodological
compatibility, which holds.

## 7. Output checksums (SHA-256)

| System | file | bytes | sha256 |
|---|---|---|---|
| AA-A6 | rmsd_protein.xvg | 67325 | `cc1231288fe202dee3cf1ca0d5424e361098003a08e32535e5aa8cc54fa7a112` |
| AA-A6 | rmsd_ligand.xvg | 67319 | `49597986592ea499a05cf6fb298d9cb4612d91402b4b986050011d6f29f12573` |
| AA-A6 | mindist_lig_active.xvg | 68591 | `6d49c782b06aea8de7ecda4e17da74b1a55a2602564d3040aa6fda661c782430` |
| AA-A6 | contacts_lig_active.xvg | 58591 | `599a73df0a694e75a6808a425bd5f0d0d8153f107f0bf9c4dcfe7381f5d2aee4` |
| AA-A6 | dist_lig_catalytic.xvg | 53356 | `922ebb2c0e5adec32f3489faee16c0cbe503714b908d6f318a4a2f9f2a2e8872` |
| AG-A6 | rmsd_protein.xvg | 67339 | `e0f76cc1f27187ee7c59b80055927d48797449f0589a8f0b79492b941f436886` |
| AG-A6 | rmsd_ligand.xvg | 67331 | `03a8e663d067459375b2d84a1afdb46fef64edd7c0c5696b964e377899d08b53` |
| AG-A6 | mindist_lig_active.xvg | 68608 | `c6581ea17662e7e7cacc6f1f7d7931261331836caf450b3e80e9f77c17e1ffb2` |
| AG-A6 | contacts_lig_active.xvg | 58608 | `9422c4d09cc898ced7cb8e2bb11a631cf3c767c38f3c8b8bdb76ad1eae1da1e3` |
| AG-A6 | dist_lig_catalytic.xvg | 53380 | `5995fcae36bf00c5d7b09fceabbdc4c580f643a39bbdae3803f078148d1690b7` |
| LP-A6 | rmsd_protein.xvg | 67336 | `005e83c04ad1663b551adf9c5e9fab59e1a153401b28a498762e35d219584249` |
| LP-A6 | rmsd_ligand.xvg | 67346 | `d4468b5659236944a91364f4bbf21d876ec7cd8e7c75c3825868c6c5f25ed5c6` |
| LP-A6 | mindist_lig_active.xvg | 68619 | `f8698a52fb3ef5d7e44f685c4e341595feed23fd6bbcb9ced5f44c7f5b39b06a` |
| LP-A6 | contacts_lig_active.xvg | 58619 | `f8495c69b7055665922fc56133378b9d48ad5075a8882361f44364707ec067de` |
| LP-A6 | dist_lig_catalytic.xvg | 53387 | `485bf0afbfb8fe4e9b33fd8f0e369fa7716ade8932bdde0ae141185c7c210c97` |
| HMG-R-25ns-A6 | rmsd_protein.xvg | 67362 | `c750e11641c727725c7f15a1c1114d29cdbb2d3ff32ddf6868a30d68e0af0efe` |
| HMG-R-25ns-A6 | rmsd_ligand.xvg | 67351 | `bd68b142fec1d1ff370ebc27cf0759afd4ef4c9827d041220bc38c381307deaf` |
| HMG-R-25ns-A6 | mindist_lig_active.xvg | 68632 | `2971af170dffe87fd4b045162556809320881a3c96b791e8b81ceb87b6c7a1e2` |
| HMG-R-25ns-A6 | contacts_lig_active.xvg | 58632 | `d7737fc135d8cc4f8b0006c87e7e4caed828bd5e0c8c627d14a5db19882f9cea` |
| HMG-R-25ns-A6 | dist_lig_catalytic.xvg | 53419 | `293a66c34525fdc155712844c1202a3b0047cb0c49a13b32282fcb11ea9caa90` |

(`active_site_mindist` and `active_site_contacts` share one `gmx mindist` invocation, so
their provenance records both files.)

## 8. Provenance

Every generated result has a `provenance.<analysis>.json` recording: source trajectory,
source TPR, index file (all with SHA-256 / size / mtime — trajectory SHA-256 skipped
above 2 GiB, size+mtime retained), resolved semantic group names, resolved numeric
group ids, analysis window, contact cutoff / catalytic semantics, full command argv and
stdin, GROMACS return code and version, UTC timestamp, output SHA-256, protocol version
and profile, measured trajectory end, and intended-production-duration evidence.

## 9. Source-integrity verification (spec §11)

| Check | Before | After | Verdict |
|---|---|---|---|
| file count under `Nuevos_sistemas` | 782 | 782 | unchanged |
| files added | — | 0 | ✅ |
| files removed | — | 0 | ✅ |
| files modified (size or mtime) | — | 0 | ✅ |
| `index.ndx` SHA-256 (all 8 systems) | — | all identical | ✅ |

`source_integrity.json` (full before/after inventory) is written at the output root.
Nothing under `/home/jesusxd/Escritorio/Nuevos_sistemas` was created, modified, or removed.
A rerun produced the same verdict.

## 10. Rerun determinism

Second invocation with the same arguments: **20 SKIP_EXISTING + 20 REVIEW_REQUIRED, 0
EXECUTED, 0 FAILED**. No validated result was recomputed or overwritten (that requires
`--force`); source integrity unchanged again.

## 11. Test results

- `tests/analysis/campaign/` — **315 passed**
- `tests/analysis/campaign/test_xanthone_short.py` — **23 passed** (full spec §8 matrix:
  trajectory longer than window, exact 25 ns, missing active-site group, missing catalytic
  group, LIG at nonstandard id, existing-output skip, failed GROMACS command, output-path
  collision, rerun determinism, source files unchanged, provenance completeness, XVG
  post-validation, equilibration-TPR guard, unresolved-group guard)
- adversarial review fixed 3 latent defects: equilibration `*.tpr` fallback removed
  (`md.tpr` now mandatory), `build_command` refuses `None` group ids, blank source
  directory now marks a system ineligible

## 12. New / changed code

- new: `analysis/campaign/xanthone_short.py`, `tests/analysis/campaign/test_xanthone_short.py`
- changed: `analysis/campaign/cli.py` (+`study run`), `cli.py` (subcommand wiring),
  `docs/study_campaign.md`, `.gitignore` (`/analysis_outputs/`)
- the `study inspect|selections|plan|analyze|analyses` commands remain read-only;
  `study run` is the only command that invokes GROMACS against scientific data.

No mechanistic-profile execution was started.
