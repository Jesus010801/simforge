# MM-PBSA provenance / configuration equivalence audit (A1–A6)

> **SUPERSEDED 2026-09-09 by `mmpbsa_canonical_source_recovery.md` +
> `mmpbsa_rerun_decision.md`.** This audit only inspected `/Escritorio/MMPBSA/`
> and concluded the A1–A5 dataset column was a "cross-batch mix that does not
> reproduce from on-disk folders → recommend a 24-system rerun." That was
> **incomplete discovery.** The canonical source was subsequently found:
> `/home/jesusxd/Escritorio/heatmap.py` → `Tesis-maestría/Heat.pdf` → thesis
> Fig. `heatmap_mmpbsa` → `mmpbsa_dataset.csv` (all 20 A1–A5 values match). Its
> methodology (thesis Methods: last 5 ns, 200 ps intervals) is **equivalent to
> A6**. The `/Escritorio/MMPBSA/` folders that disagree are a **discarded
> preliminary `-dt 100` batch**. **Decision: retain A1–A5, no 24-system rerun.**
> The per-term configuration comparison below remains valid; the "not
> reproducible / rerun 24" conclusion does not.

**Date:** 2026-09-09
**Scope:** all A1–A6 × 4-target MM-PBSA results present on disk.
**Mandate:** decide comparability from **scientific configuration and provenance**,
not from numerical magnitude and not from whether the run was launched manually
or by SimForge. Report any real parameter mismatch. **No MM-PBSA was rerun and
no ML dataset was modified** in this pass.

## Sources on disk

| batch | location | systems |
|---|---|---|
| A1–A5 (per target) | `/home/jesusxd/Escritorio/MMPBSA/<TARGET>-<A#>/` | AA-{A1,A3,A4,A5}, AG-{A3,A4,A5}, HMGR-{A1,A2,A3,A4,A5}, LP-{A1,A2,A3,A4,A5} |
| A6 (25 ns) | `/home/jesusxd/Escritorio/Nuevos_sistemas/{AA-A6,AG-A6,HMG-R-25ns-A6,LP-A6}/` | the 4 A6 short systems |
| A3/A6 (200 ns, cofactor) | `Nuevos_sistemas/{A3-HMG-R,HMG-R-200ns-A6,system_A3_COA/*-MMPBSA,system_A6_COA/*-MMPBSA}` | mechanistic |

Not on disk (dataset rows exist): **AA-A2, AG-A1, AG-A2** — these came from a batch
no longer present in the folders inspected.

## Configuration comparison

### `mmpbsa.mdp` — two distinct files

| sha256 (first 16) | systems | notes |
|---|---|---|
| `954785fb180552bc…` | AA-A1/A3/A4/A5, **AA-A6**, AG-A3/A4/A5, **AG-A6**, HMGR-A1, **HMG-R-25ns-A6**, LP-A1/A2/A3/A4/A5, **LP-A6** (17) | the **majority** file — the one every A6 short run used |
| `1ee1293d34cf0cb1…` | HMGR-A2, HMGR-A3, HMGR-A4, HMGR-A5 (4) | a separate HMG-CoA-reductase batch |

The two `mmpbsa.mdp` files carry the same PB/SA physics that matters for
end-state ΔG:

| parameter | value (both files) |
|---|---|
| solute dielectric `pdie` | 2 |
| solvent dielectric `sdie` | 78.3 |
| reference dielectric `vdie` | 1 |
| PB solver | `lpbe` (linearised) |
| ion charge/radius/conc | ±1 / 0.95 & 1.81 Å / 0.150 M |
| grid spacing | 0.5 Å fine |
| `srfm` / `chgm` | `smol` / `spl4` |
| apolar model | SASA-only: `gamma` 0.0226778 kJ/mol/Å², `sasaconst` 3.84928 kJ/mol; SAV `press` 0.234304; WCA off |
| temperature | 309.65 K |

*(The 4-system difference between the two mdp files was not isolated to a
physics term in this audit; HMGR-A2..A5 also differ in residue count — see below —
so that batch is treated as separate regardless.)*

### g_mmpbsa binary

All inspected runs (A1–A5 and A6) used the **same build**:
`g_mmpbsa run, 2025.0-dev-20250210-6949615-unknown` (from `energy_MM.xvg` headers).
Manual vs automated launch is **not** a distinguishing factor — every run was a
manual `g_mmpbsa` invocation; SimForge did not run any of them, and that has no
bearing on comparability.

### Trajectory + interval + sampling

| item | A1–A5 (25 ns) | A6 (25 ns) | verdict |
|---|---|---|---|
| input trajectory | `mdfit.xtc` (PBC-corrected rot+trans fit) | `mdfit.xtc` | **same** |
| preprocessing chain | `md.xtc → -pbc res -ur compact -center → mdcenter.xtc → -fit rot+trans → mdfit.xtc` | identical | **same** |
| `-unit1 / -unit2` | `Protein` / `LIG` (by name) | `Protein` / `LIG` | **same** |
| analysis interval | `-b 20000 -e 25000` (last 5 ns) | `-b 20000 -e 25000` (last 5 ns) | **same** |
| stride | `-dt 100` → **51 snapshots** | `-dt 200` → **26 snapshots** | **DIFFERENT** |

**The one real configuration difference: snapshot stride (100 ps vs 200 ps),
i.e. 51 vs 26 frames over the identical 20–25 ns window.** Same interval, same
trajectory, same PB/SA physics — coarser sampling only. Effect on the end-state
ΔG *mean* is expected to be small (both sample the same 5 ns); it inflates the
*standard error* of the A6 means relative to A1–A5.

### Per-residue decomposition availability

| system(s) | `residues_energy_summary.csv` rows |
|---|---|
| AA-* (incl. AA-A6) | 496 |
| AG-* (incl. AG-A6) | 945 |
| LP-* (incl. LP-A6) | 450 |
| HMGR-A1, **HMG-R-25ns-A6** | **1615** |
| HMGR-A2..A5 | **815** |

A6 per-residue decomposition matches its A1 counterpart in every target family
(HMGR-A1 ↔ HMG-R-25ns-A6 both 1615) — so a catalytic-residue A1↔A6 comparison is
structurally ready. The HMGR-A2..A5 815-row batch is **not** aligned with these
and is out of scope for an A1/A6 comparison.

## The dataset-level mismatch (the real problem — reported, not fixed)

`mmpbsa_dataset_v2_A1-A6.csv` values vs the on-disk `energy_summary.csv` Totals:

| dataset row | dataset ΔG (kJ/mol) | on-disk folder | on-disk Total | agree? |
|---|---|---|---|---|
| A1 / alpha_amylase | −165.2 | AA-A1 | −165.2 | ✅ |
| A1 / hmgcoa_reductase | −190.9 | HMGR-A1 | −190.9 | ✅ |
| A1 / pancreatic_lipase | −188.9 | LP-A1 | **−168.9** | ❌ (Δ 20) |
| A1 / alpha_glucosidase | −174.9 | *(no AG-A1 folder)* | — | n/a |
| A3 / alpha_amylase | −157.3 | AA-A3 | −157.3 | ✅ |
| A3 / alpha_glucosidase | −115.7 | AG-A3 | −115.7 | ✅ |
| A3 / pancreatic_lipase | −60.9 | LP-A3 | −60.9 | ✅ |
| A5 / alpha_amylase | −162.6 | AA-A5 | **−242.6** | ❌ (Δ 80) |
| A5 / alpha_glucosidase | −181.7 | AG-A5 | **−231.7** | ❌ (Δ 50) |
| A5 / hmgcoa_reductase | −140.7 | HMGR-A5 | **−198.2** | ❌ (Δ 57) |
| A5 / pancreatic_lipase | −179.9 | LP-A5 | **−219.9** | ❌ (Δ 40) |
| A6 / alpha_amylase | −161.0 | AA-A6 | −161.0 | ✅ |
| A6 / alpha_glucosidase | −144.2 | AG-A6 | −144.2 | ✅ |
| A6 / hmgcoa_reductase | −266.2 | HMG-R-25ns-A6 | −266.2 | ✅ |
| A6 / pancreatic_lipase | −304.9 | LP-A6 | −304.9 | ✅ |

The **A6 dataset values ARE the on-disk A6 folder values.** The **A5 row and the
A1/lipase row are NOT** — they are 20–80 kJ/mol less negative than any on-disk
folder, with systematically smaller-magnitude vdW/electrostatic terms. Per the
earlier audit (`ML/.../targeted_A6_completion_2026-09-09/qc_and_unresolved_items.md`
Q4), those dataset numbers come from an **earlier (July-2025) MM-PBSA batch with
different ligand partial charges**, which is not on disk.

## Conclusions

1. **Methodological configuration is equivalent** between the A6 runs and the
   on-disk A1–A5 runs that share `mmpbsa.mdp` `954785fb1805…`: same fitted
   trajectory, same 20–25 ns interval, same PB/SA physics, same dielectrics,
   same temperature, same g_mmpbsa build, same `Protein`/`LIG` selection.
   **Manual execution is not a comparability criterion** and does not make A6
   "a different methodology".

2. **One real parameter difference exists: snapshot stride** — A6 used `-dt 200`
   (26 frames), A1–A5 used `-dt 100` (51 frames) over the same window. This is a
   sampling-density difference, not a physics or interval difference. It should
   be equalised in any future single-batch rerun; it is not on its own grounds
   to call the results incomparable.

3. **HMGR-A2..A5** used a different `mmpbsa.mdp` and a different residue count
   (815 vs 1615) — that is a genuinely separate batch and must not be pooled
   with HMGR-A1 / HMG-R-25ns-A6.

4. **The blocking problem is at the dataset level, not the method level:** the
   `mmpbsa_dataset_v2` A1–A5 column is a **cross-batch mix**. Several rows
   (all of A5, A1/lipase) do not reproduce from any on-disk folder and predate
   the current charge protocol. Until a **single-batch rerun of all 24
   systems** (A1–A6 × 4 targets — one topology/charge protocol, one window, one
   stride, one `mmpbsa.mdp`) is done and `mmpbsa_dataset` is regenerated from
   scratch, `mmpbsa_dg_kjmol` should be treated as **ordinal within a
   single batch only**.

**No MM-PBSA recompute performed. No `mmpbsa_dataset_v2` / `master_v2` cell
changed.** The "different batch = incompatible" framing is *partly* wrong (A6 vs
the `954785fb1805` A1–A5 subset is method-equivalent) and *partly* right (the
A5 / A1-lipase dataset values are from an unreproduced batch). Recommended
action — a 24-system single-batch rerun — is **reported, not executed**.
