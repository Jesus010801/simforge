# Membrane-water atlas implementation checkpoints

The region-atlas design in the task is authoritative. Existing uncommitted
repository changes are retained; checkpoints are documented phases, not commits
that would accidentally include those changes.

1. Identity/PBC: independent molecular-instance IDs; water excluded from geometry;
   explicit amino-acid recognition; full cell matrix and exact periodic lattice.
   Rotated orthogonal cells supported, skew numerical backend explicitly rejected.
2. Exclusion/clearance: implemented with heavy-atom exclusion and nominal water-probe clearance on the shared cell lattice; hard clashes are a separate cleanup signal. Focused clearance tests pass.
3. Region discovery: implemented as periodic full-void components followed by independently labeled membrane-core components; water coordinates are not inputs to field construction or discovery. Topology fixtures A-P pass.
4. Graph/hierarchy: implemented as region nodes with explicit upper/lower reservoir portals and pathway summaries. Elementary components remain explicit; shared-portal/pathway subdivision within one connected volume is still limited and recorded as a known limitation.
5. Classification: implemented as per-region attributes, evidence, and uncertainty; species identity is explicit and unknown residues are OTHER.
6. Policy/cleanup: implemented as swappable conservative policy and atomic molecule decisions, with retained coordinate records copied unchanged and topology SOL synchronization by molecule identity.
7. Production integration: public spatial-classifier and pore-hydration entrypoints delegate to v2; assembly routing uses v2 without TM annotations; historical slab path is explicit-only; runtime core-water counts are advisory. Focused adapter tests pass.

Real-system audit artifacts are present under reports/membrane_water_v2_real_validation/. Retained atom records and complete water molecules validate; the copied topology passes grompp with one allowed warning for net charge 25.0002. Atlas topology is NOT converged: 0.20 nm yields 304 regions and 0.22 nm yields 230 (24.3% change). The static XY review shows a plausible pore but many isolated ambiguous/lipid-facing points; scientific visual approval is pending. Do not release cleanup as production-ready until region aggregation/classification stabilizes and the structures are reviewed in VMD.

The A-P suite currently exercises synthetic free-space topology masks rather than coordinates generated from molecular obstacle sets. It passes the requested component, portal, PBC-seam, multiple-pore, and vestibule relationships, but does not yet validate the entire exclusion-to-atlas pipeline for every named fixture. Bifurcated shared-junction decomposition is also incomplete: each path is not yet represented as a distinct child of one shared junction/vestibule. The numerical PBC backend supports orthogonal GRO cells, including rotated frames; skew cells fail safe and remain unsupported.

The legacy backbone/region helper files remain importable because existing unit tests exercise some low-level utilities, but public production routes no longer import them for classification; deleting or quarantining those private helpers is still migration cleanup.

Full repository tests: 2716 passed, 4 unrelated ligand hydrogenation failures, 6 skipped, 9 xfailed, 1 xpassed. The latest focused membrane-water/gate/regression suite passed 103 tests with 1 xpass.
