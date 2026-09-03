# SimForge talk demo

A deterministic ~1-minute demo. No MD is run.

```bash
./demo/run_demo.sh
```

Story: **environment → declarative spec → guarded construction → GROMACS
validation → provenance.**

1. `simforge doctor` — the environment SimForge sits on top of.
2. `demo/system.yaml` — a whole protein + 2-component system in ~20 lines.
3. `simforge build` — assembles the topology + coordinates. The two components
   (`MTH`, `EOL`) each ship an `[ atomtypes ]` block that declares `opls_800`
   for a **different element**. SimForge namespaces the collision
   (`MTH_opls_800` / `EOL_opls_800`) instead of silently picking one.
4. `simforge validate-system` — 8 checks, ending in `gmx grompp -maxwarn 0`.
5. `simforge inspect-run` — the provenance manifest; the atomtype rename map is
   recorded there.

This is a **software demonstration** (6-residue peptide, methane + methanol),
not a scientific benchmark.

> The "Reproducibility advisory" panel in `expected_output.txt` appears because
> the capture was taken from a modified working tree. From a clean clone it
> does not appear — that panel *is* the feature working.

## Fallback

`./demo/run_demo.sh --capture` regenerates `demo/expected_output.txt`. If the
live demo fails on stage, show that file (or a screen recording of a prior
run). Do **not** run `mdrun`, the membrane builder, or the full test suite
live.
