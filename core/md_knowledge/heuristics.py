"""core/md_knowledge/heuristics.py — Observable heuristics for 10 MD quantities × context."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from core.md_knowledge.contexts import SystemContext, SoftRange, SYSTEM_CONTEXTS


@dataclass
class InterpretationRule:
    """A single conditional interpretation rule for an observable."""
    condition:      str           # human-readable condition description
    interpretation: str
    severity:       str           # "info" | "warning" | "error"
    recommendation: str = ""


@dataclass
class ObservableHeuristic:
    observable:     str
    context:        SystemContext
    soft_range:     SoftRange
    rules:          list[InterpretationRule]
    citations:      list[str] = field(default_factory=list)
    notes:          list[str] = field(default_factory=list)


# ── Builder helpers ───────────────────────────────────────────────────────────

def _rmsd_rules(rmsd_range: SoftRange) -> list[InterpretationRule]:
    hi_exc = rmsd_range.excellent_high
    hi_acc = rmsd_range.acceptable_high
    hi_warn = rmsd_range.warning_high or 1.0
    return [
        InterpretationRule(
            condition=f"RMSD plateau < {hi_exc} nm",
            interpretation="Excellent structural convergence. System is well equilibrated.",
            severity="info",
        ),
        InterpretationRule(
            condition=f"RMSD plateau {hi_exc}–{hi_acc} nm",
            interpretation="Acceptable convergence; minor local flexibility or slow relaxation.",
            severity="info",
            recommendation="Visual inspection recommended; consider extending production.",
        ),
        InterpretationRule(
            condition=f"RMSD plateau {hi_acc}–{hi_warn} nm",
            interpretation="Elevated RMSD. May indicate partial unfolding, large domain motion, or poor equilibration.",
            severity="warning",
            recommendation="Inspect trajectory for secondary structure loss. Check whether RMSD drift is still ongoing.",
        ),
        InterpretationRule(
            condition=f"RMSD > {hi_warn} nm",
            interpretation="Very high RMSD. Likely unfolding, nonphysical behavior, or wrong reference structure.",
            severity="error",
            recommendation="Inspect trajectory immediately. Check minimization, solvation, and force-field parameters.",
        ),
        InterpretationRule(
            condition="RMSD shows monotonic drift (no plateau)",
            interpretation="System has not reached equilibrium. Thermodynamic averages are unreliable.",
            severity="warning",
            recommendation="Extend equilibration or simulation. Check thermostat/barostat parameters.",
        ),
        InterpretationRule(
            condition="RMSD shows jump then plateau",
            interpretation="Conformational transition detected. May be biologically relevant.",
            severity="info",
            recommendation="Verify with visual inspection; report transition time and magnitude.",
        ),
    ]


def _energy_rules() -> list[InterpretationRule]:
    return [
        InterpretationRule(
            condition="Potential energy stable with low variance",
            interpretation="System is thermodynamically stable.",
            severity="info",
        ),
        InterpretationRule(
            condition="Potential energy drifting",
            interpretation="System has not reached thermal equilibrium.",
            severity="warning",
            recommendation="Extend NVT/NPT equilibration.",
        ),
        InterpretationRule(
            condition="Large energy spikes",
            interpretation="Numerical instability — likely timestep too large or bad contacts.",
            severity="error",
            recommendation="Reduce timestep. Check for clashes in initial structure.",
        ),
    ]


def _temperature_rules() -> list[InterpretationRule]:
    return [
        InterpretationRule(
            condition="T within ±5 K of target",
            interpretation="Thermostat maintaining target temperature correctly.",
            severity="info",
        ),
        InterpretationRule(
            condition="T drifting or > ±20 K from target",
            interpretation="Thermostat not converged. Heat bath coupling may be too loose.",
            severity="warning",
            recommendation="Check tau_t parameter. Verify thermostat type is appropriate.",
        ),
        InterpretationRule(
            condition="T > 400 K",
            interpretation="Physically implausible temperature — system may be exploding.",
            severity="error",
            recommendation="Abort or restart with better minimization and reduced timestep.",
        ),
    ]


def _pressure_rules() -> list[InterpretationRule]:
    return [
        InterpretationRule(
            condition="|P - P_target| < 100 bar (time-average)",
            interpretation="Barostat maintaining target pressure within acceptable fluctuations.",
            severity="info",
            recommendation="Note: instantaneous pressure fluctuations of ±500 bar are normal in MD.",
        ),
        InterpretationRule(
            condition="Pressure drifting over simulation",
            interpretation="Density has not converged; NPT equilibration may be incomplete.",
            severity="warning",
            recommendation="Check barostat compressibility parameter. Extend NPT equilibration.",
        ),
    ]


def _rg_rules(rg_range: SoftRange) -> list[InterpretationRule]:
    return [
        InterpretationRule(
            condition=f"Rg within [{rg_range.excellent_low}, {rg_range.excellent_high}] nm (plateau)",
            interpretation="Compact, stable fold maintained throughout simulation.",
            severity="info",
        ),
        InterpretationRule(
            condition="Rg increasing monotonically",
            interpretation="Protein is unfolding or expanding. Check stability.",
            severity="warning",
            recommendation="Inspect trajectory for unfolding events. Verify solvation and pH.",
        ),
        InterpretationRule(
            condition="Rg decreasing monotonically",
            interpretation="Protein is collapsing. May indicate over-compaction or aggregation.",
            severity="warning",
            recommendation="Check periodic boundary conditions and protein-image distances.",
        ),
    ]


# ── OBSERVABLE_HEURISTICS dict ────────────────────────────────────────────────

OBSERVABLE_HEURISTICS: dict[str, dict[SystemContext, ObservableHeuristic]] = {}


def _register(observable: str, context: SystemContext, heuristic: ObservableHeuristic) -> None:
    OBSERVABLE_HEURISTICS.setdefault(observable, {})[context] = heuristic


def _build_all() -> None:
    for ctx, profile in SYSTEM_CONTEXTS.items():
        # ── RMSD ──────────────────────────────────────────────────────────────
        _register("rmsd", ctx, ObservableHeuristic(
            observable="rmsd",
            context=ctx,
            soft_range=profile.rmsd_nm,
            rules=_rmsd_rules(profile.rmsd_nm),
            citations=[
                "Knapp & Nilsson, J Chem Theory Comput, 2008",
                "Grossfield & Zuckerman, Annu Rev Biophys, 2009",
            ],
            notes=profile.notes,
        ))

        # ── Radius of gyration ────────────────────────────────────────────────
        _register("rg", ctx, ObservableHeuristic(
            observable="rg",
            context=ctx,
            soft_range=profile.rg_nm,
            rules=_rg_rules(profile.rg_nm),
            citations=["Lobanov et al., Mol Biol, 2008"],
        ))

        # ── Temperature ───────────────────────────────────────────────────────
        _register("temperature", ctx, ObservableHeuristic(
            observable="temperature",
            context=ctx,
            soft_range=profile.temperature_k,
            rules=_temperature_rules(),
        ))

        # ── Pressure ──────────────────────────────────────────────────────────
        _register("pressure", ctx, ObservableHeuristic(
            observable="pressure",
            context=ctx,
            soft_range=profile.pressure_bar,
            rules=_pressure_rules(),
            notes=["Instantaneous pressure fluctuations ±500 bar are normal in MD."],
        ))

        # ── Potential energy ──────────────────────────────────────────────────
        _register("potential_energy", ctx, ObservableHeuristic(
            observable="potential_energy",
            context=ctx,
            soft_range=profile.energy_stability_kj,
            rules=_energy_rules(),
        ))

        # ── RMSF ─────────────────────────────────────────────────────────────
        rmsf_range = SoftRange(
            excellent_low=0.0,
            excellent_high=0.15 if ctx not in (SystemContext.IDR, SystemContext.PEPTIDE) else 0.5,
            acceptable_low=0.0,
            acceptable_high=0.30 if ctx not in (SystemContext.IDR, SystemContext.PEPTIDE) else 1.0,
            warning_low=0.0,
            warning_high=0.50 if ctx not in (SystemContext.IDR, SystemContext.PEPTIDE) else 2.0,
            unit="nm",
            note="Per-residue flexibility; high RMSF in loops is normal",
        )
        _register("rmsf", ctx, ObservableHeuristic(
            observable="rmsf",
            context=ctx,
            soft_range=rmsf_range,
            rules=[
                InterpretationRule(
                    condition="RMSF < 0.1 nm for structured regions",
                    interpretation="Low backbone flexibility — well-folded core.",
                    severity="info",
                ),
                InterpretationRule(
                    condition="RMSF > 0.3 nm in structured regions",
                    interpretation="Elevated flexibility. May indicate local unfolding or hinge motion.",
                    severity="warning",
                    recommendation="Check RMSF vs B-factor correlation with crystal structure.",
                ),
                InterpretationRule(
                    condition="RMSF very high in loop regions",
                    interpretation="Expected for flexible loops and termini.",
                    severity="info",
                ),
            ],
            citations=["Skjaerven et al., BMC Bioinformatics, 2009"],
        ))

        # ── SASA ─────────────────────────────────────────────────────────────
        sasa_range = SoftRange(
            excellent_low=0.0, excellent_high=200.0,   # nm² — highly variable, use stability
            acceptable_low=0.0, acceptable_high=500.0,
            warning_low=0.0,   warning_high=1000.0,
            unit="nm²",
            note="SASA absolute value is system-size dependent; stability matters more",
        )
        _register("sasa", ctx, ObservableHeuristic(
            observable="sasa",
            context=ctx,
            soft_range=sasa_range,
            rules=[
                InterpretationRule(
                    condition="SASA stable at plateau",
                    interpretation="Protein surface exposure is stable — fold is maintained.",
                    severity="info",
                ),
                InterpretationRule(
                    condition="SASA increasing monotonically",
                    interpretation="Protein is unfolding or expanding — hydrophobic core exposed.",
                    severity="warning",
                    recommendation="Correlate with RMSD drift and secondary structure content.",
                ),
            ],
            citations=["Fraternali & Cavallo, Nucleic Acids Res, 2002"],
        ))

        # ── H-bonds ───────────────────────────────────────────────────────────
        hbond_range = SoftRange(
            excellent_low=0.0, excellent_high=500.0,   # count — system-size dependent
            acceptable_low=0.0, acceptable_high=1000.0,
            warning_low=0.0,   warning_high=2000.0,
            unit="count",
            note="Use stability (plateau) rather than absolute count for convergence",
        )
        _register("hbonds", ctx, ObservableHeuristic(
            observable="hbonds",
            context=ctx,
            soft_range=hbond_range,
            rules=[
                InterpretationRule(
                    condition="H-bond count stable",
                    interpretation="Hydrogen-bond network is stable — secondary/tertiary structure maintained.",
                    severity="info",
                ),
                InterpretationRule(
                    condition="H-bond count decreasing",
                    interpretation="Loss of internal H-bonds — potential unfolding or denaturation.",
                    severity="warning",
                    recommendation="Correlate with RMSD and secondary structure content.",
                ),
                InterpretationRule(
                    condition="H-bond count < 20% of expected for system size",
                    interpretation="Very few intramolecular H-bonds. Check topology or structure preparation.",
                    severity="error",
                ),
            ],
        ))

        # ── Ligand-pocket distance ────────────────────────────────────────────
        ligand_range = SoftRange(
            excellent_low=0.0, excellent_high=0.30,
            acceptable_low=0.0, acceptable_high=0.50,
            warning_low=0.0,   warning_high=0.80,
            unit="nm",
            note="Distance between ligand centroid and binding pocket centroid",
        )
        _register("ligand_pocket_distance", ctx, ObservableHeuristic(
            observable="ligand_pocket_distance",
            context=ctx,
            soft_range=ligand_range,
            rules=[
                InterpretationRule(
                    condition="Ligand-pocket distance < 0.3 nm (stable)",
                    interpretation="Ligand remains bound in the binding pocket.",
                    severity="info",
                ),
                InterpretationRule(
                    condition="Ligand-pocket distance 0.3–0.8 nm",
                    interpretation="Ligand is moving at the pocket periphery or partially unbound.",
                    severity="warning",
                    recommendation="Inspect trajectory. Evaluate binding free energy if relevant.",
                ),
                InterpretationRule(
                    condition="Ligand-pocket distance > 0.8 nm",
                    interpretation="Ligand has dissociated from the binding pocket.",
                    severity="error",
                    recommendation=(
                        "Check protonation state, force-field parameters, and initial pose. "
                        "Consider restraints during equilibration."
                    ),
                ),
                InterpretationRule(
                    condition="Ligand-pocket distance increasing monotonically",
                    interpretation="Ligand is actively dissociating during the simulation.",
                    severity="error",
                ),
            ],
            citations=["Buch et al., PNAS, 2011"],
        ))

        # ── Secondary structure retention ─────────────────────────────────────
        ss_range = SoftRange(
            excellent_low=0.90, excellent_high=1.0,
            acceptable_low=0.75, acceptable_high=1.0,
            warning_low=0.50,   warning_high=1.0,
            unit="fraction",
            note="Fraction of native secondary structure retained (helix + beta-strand)",
        )
        _register("secondary_structure", ctx, ObservableHeuristic(
            observable="secondary_structure",
            context=ctx,
            soft_range=ss_range,
            rules=[
                InterpretationRule(
                    condition="SS retention > 90%",
                    interpretation="Excellent: secondary structure fully preserved.",
                    severity="info",
                ),
                InterpretationRule(
                    condition="SS retention 75–90%",
                    interpretation="Minor secondary structure fluctuations. Termini and loops may be fraying.",
                    severity="info",
                    recommendation="Check which regions lose structure (termini fraying is normal).",
                ),
                InterpretationRule(
                    condition="SS retention 50–75%",
                    interpretation="Significant secondary structure loss — potential partial unfolding.",
                    severity="warning",
                    recommendation="Correlate with RMSD increase. Inspect trajectory.",
                ),
                InterpretationRule(
                    condition="SS retention < 50%",
                    interpretation="Severe unfolding or force-field incompatibility.",
                    severity="error",
                    recommendation=(
                        "Verify force field, pH/protonation, disulfide bonds, and "
                        "simulation parameters."
                    ),
                ),
            ],
            citations=["Best et al., J Chem Theory Comput, 2012"],
        ))


_build_all()


def get_heuristic(observable: str, context: SystemContext) -> ObservableHeuristic | None:
    """Return the heuristic for an observable in a given context, falling back to UNKNOWN."""
    by_ctx = OBSERVABLE_HEURISTICS.get(observable, {})
    return by_ctx.get(context) or by_ctx.get(SystemContext.UNKNOWN)
