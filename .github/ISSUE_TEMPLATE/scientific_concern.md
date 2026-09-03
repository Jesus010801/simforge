---
name: Scientific concern
about: You believe SimForge produced a scientifically wrong topology, system, or analysis
title: "[scientific] "
labels: scientific-concern
assignees: ""
---

## Summary

<!-- Briefly: what did SimForge produce that you believe is scientifically
     incorrect? e.g. wrong atom/charge counts, wrong protonation, misplaced
     ions/water, broken coordinate merge, implausible free-energy landscape. -->

## What you expected vs. what you got

**Expected:**

**Got:**

<!-- Be specific and quantitative where possible (atom counts, total charge,
     box dimensions, RMSD, energies, ...). -->

## Attachments (please attach all that apply)

- [ ] The spec YAML you ran
- [ ] `provenance.json` from the run
- [ ] `validation_report.json` (from `simforge validate-system`)
- [ ] Full `gmx grompp` output (with `-maxwarn 0`)
- [ ] `gmx mdrun` / analysis output, if relevant
- [ ] Any input structures/topologies needed to reproduce (or a link)

## Environment

- SimForge version / commit:
- GROMACS version:
- Force field / water model:
- OS / Python version:

## Additional context

<!-- Literature references, prior known-good results, or anything else that
     supports the concern. -->
