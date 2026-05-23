"""core/md_knowledge/contexts.py — System context definitions with soft observable ranges."""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class SystemContext(str, Enum):
    GLOBULAR_PROTEIN         = "globular_protein"
    MEMBRANE_PROTEIN         = "membrane_protein"
    IDR                      = "idr"                    # intrinsically disordered
    PROTEIN_LIGAND_COMPLEX   = "protein_ligand_complex"
    MULTIMERIC_COMPLEX       = "multimeric_complex"
    ENZYME                   = "enzyme"
    PEPTIDE                  = "peptide"
    MEMBRANE_SYSTEM          = "membrane_system"        # pure lipid bilayer
    FLEXIBLE_DOMAIN_PROTEIN  = "flexible_domain_protein"
    UNKNOWN                  = "unknown"


@dataclass
class SoftRange:
    """A soft [low, high] range with descriptive labels for each zone."""
    excellent_low:  float
    excellent_high: float
    acceptable_low: float
    acceptable_high: float
    warning_low:    float | None  # None = open lower bound
    warning_high:   float | None  # None = open upper bound
    unit:           str = ""
    note:           str = ""

    def classify(self, value: float) -> str:
        """Return 'excellent' | 'acceptable' | 'warning' | 'critical'."""
        if self.excellent_low <= value <= self.excellent_high:
            return "excellent"
        if self.acceptable_low <= value <= self.acceptable_high:
            return "acceptable"
        lo = self.warning_low if self.warning_low is not None else float("-inf")
        hi = self.warning_high if self.warning_high is not None else float("inf")
        if lo <= value <= hi:
            return "warning"
        return "critical"

    def score(self, value: float) -> float:
        """Return a [0,1] quality score — 1.0 = dead centre of excellent range."""
        centre = (self.excellent_low + self.excellent_high) / 2
        half_exc = (self.excellent_high - self.excellent_low) / 2 + 1e-9
        half_acc = (self.acceptable_high - self.acceptable_low) / 2 + 1e-9

        dist = abs(value - centre)
        if dist <= half_exc:
            return 1.0 - 0.2 * (dist / half_exc)   # [0.8, 1.0]
        if dist <= half_acc:
            return 0.8 - 0.3 * ((dist - half_exc) / (half_acc - half_exc + 1e-9))
        return max(0.0, 0.5 - 0.5 * ((dist - half_acc) / (half_acc + 1e-9)))


@dataclass
class ContextProfile:
    name:         str
    description:  str
    rmsd_nm:      SoftRange
    rg_nm:        SoftRange
    temperature_k: SoftRange
    pressure_bar:  SoftRange
    energy_stability_kj: SoftRange   # std of potential energy / atom (kJ/mol)
    # Optional context-specific notes
    notes:        list[str] = field(default_factory=list)


# ── Context definitions ───────────────────────────────────────────────────────

SYSTEM_CONTEXTS: dict[SystemContext, ContextProfile] = {

    SystemContext.GLOBULAR_PROTEIN: ContextProfile(
        name="Globular Protein",
        description="Well-folded soluble protein with defined tertiary structure.",
        rmsd_nm=SoftRange(
            excellent_low=0.0, excellent_high=0.20,
            acceptable_low=0.0, acceptable_high=0.35,
            warning_low=0.0,   warning_high=0.50,
            unit="nm",
            note="RMSD plateau < 0.2 nm is excellent; > 0.5 nm suggests unfolding",
        ),
        rg_nm=SoftRange(
            excellent_low=1.0, excellent_high=2.5,
            acceptable_low=0.8, acceptable_high=3.0,
            warning_low=0.5,   warning_high=4.0,
            unit="nm",
        ),
        temperature_k=SoftRange(
            excellent_low=295.0, excellent_high=305.0,
            acceptable_low=290.0, acceptable_high=310.0,
            warning_low=280.0,   warning_high=320.0,
            unit="K",
        ),
        pressure_bar=SoftRange(
            excellent_low=-50.0, excellent_high=50.0,
            acceptable_low=-200.0, acceptable_high=200.0,
            warning_low=-500.0,  warning_high=500.0,
            unit="bar",
        ),
        energy_stability_kj=SoftRange(
            excellent_low=0.0, excellent_high=0.5,
            acceptable_low=0.0, acceptable_high=1.5,
            warning_low=0.0,   warning_high=3.0,
            unit="kJ/mol/atom",
        ),
        notes=[
            "Backbone RMSD is the primary convergence indicator.",
            "RMSD > 0.35 nm at plateau warrants visual inspection for partial unfolding.",
        ],
    ),

    SystemContext.MEMBRANE_PROTEIN: ContextProfile(
        name="Membrane Protein",
        description="Integral or peripheral membrane protein in lipid bilayer.",
        rmsd_nm=SoftRange(
            excellent_low=0.0, excellent_high=0.30,
            acceptable_low=0.0, acceptable_high=0.50,
            warning_low=0.0,   warning_high=0.80,
            unit="nm",
            note="Higher RMSD acceptable due to loop flexibility and bilayer coupling",
        ),
        rg_nm=SoftRange(
            excellent_low=1.5, excellent_high=4.0,
            acceptable_low=1.0, acceptable_high=5.0,
            warning_low=0.5,   warning_high=7.0,
            unit="nm",
        ),
        temperature_k=SoftRange(
            excellent_low=295.0, excellent_high=305.0,
            acceptable_low=290.0, acceptable_high=310.0,
            warning_low=280.0,   warning_high=320.0,
            unit="K",
        ),
        pressure_bar=SoftRange(
            excellent_low=-50.0, excellent_high=50.0,
            acceptable_low=-200.0, acceptable_high=200.0,
            warning_low=-500.0,  warning_high=500.0,
            unit="bar",
        ),
        energy_stability_kj=SoftRange(
            excellent_low=0.0, excellent_high=0.8,
            acceptable_low=0.0, acceptable_high=2.0,
            warning_low=0.0,   warning_high=4.0,
            unit="kJ/mol/atom",
        ),
        notes=[
            "Use TM-helix RMSD separately from loop RMSD for convergence assessment.",
            "Bilayer APL and thickness are important additional convergence metrics.",
            "Semiisotropic pressure coupling required for correct bilayer dynamics.",
        ],
    ),

    SystemContext.IDR: ContextProfile(
        name="Intrinsically Disordered Region/Protein",
        description="Protein lacking stable tertiary structure; high RMSD is physically correct.",
        rmsd_nm=SoftRange(
            excellent_low=0.5, excellent_high=3.0,
            acceptable_low=0.3, acceptable_high=5.0,
            warning_low=0.1,   warning_high=8.0,
            unit="nm",
            note="High RMSD is expected; convergence assessed via Rg distribution width",
        ),
        rg_nm=SoftRange(
            excellent_low=0.8, excellent_high=4.0,
            acceptable_low=0.5, acceptable_high=6.0,
            warning_low=0.3,   warning_high=8.0,
            unit="nm",
        ),
        temperature_k=SoftRange(
            excellent_low=295.0, excellent_high=305.0,
            acceptable_low=290.0, acceptable_high=310.0,
            warning_low=280.0,   warning_high=320.0,
            unit="K",
        ),
        pressure_bar=SoftRange(
            excellent_low=-50.0, excellent_high=50.0,
            acceptable_low=-200.0, acceptable_high=200.0,
            warning_low=-500.0,  warning_high=500.0,
            unit="bar",
        ),
        energy_stability_kj=SoftRange(
            excellent_low=0.0, excellent_high=1.0,
            acceptable_low=0.0, acceptable_high=2.5,
            warning_low=0.0,   warning_high=5.0,
            unit="kJ/mol/atom",
        ),
        notes=[
            "Do NOT use RMSD as primary convergence metric for IDRs.",
            "Rg distribution, SASA, and secondary structure content are more informative.",
            "Enhanced sampling methods (REST2, metadynamics) are often necessary.",
        ],
    ),

    SystemContext.PROTEIN_LIGAND_COMPLEX: ContextProfile(
        name="Protein–Ligand Complex",
        description="Folded protein with a small-molecule ligand in the binding pocket.",
        rmsd_nm=SoftRange(
            excellent_low=0.0, excellent_high=0.20,
            acceptable_low=0.0, acceptable_high=0.30,
            warning_low=0.0,   warning_high=0.45,
            unit="nm",
            note="Protein backbone RMSD; ligand RMSD should be assessed separately",
        ),
        rg_nm=SoftRange(
            excellent_low=1.0, excellent_high=2.5,
            acceptable_low=0.8, acceptable_high=3.0,
            warning_low=0.5,   warning_high=4.0,
            unit="nm",
        ),
        temperature_k=SoftRange(
            excellent_low=295.0, excellent_high=305.0,
            acceptable_low=290.0, acceptable_high=310.0,
            warning_low=280.0,   warning_high=320.0,
            unit="K",
        ),
        pressure_bar=SoftRange(
            excellent_low=-50.0, excellent_high=50.0,
            acceptable_low=-200.0, acceptable_high=200.0,
            warning_low=-500.0,  warning_high=500.0,
            unit="bar",
        ),
        energy_stability_kj=SoftRange(
            excellent_low=0.0, excellent_high=0.5,
            acceptable_low=0.0, acceptable_high=1.5,
            warning_low=0.0,   warning_high=3.0,
            unit="kJ/mol/atom",
        ),
        notes=[
            "Monitor ligand-pocket distance as primary stability indicator.",
            "Ligand RMSD > 0.3 nm from starting pose may indicate pocket exit.",
            "H-bond occupancy between ligand and key residues is critical.",
        ],
    ),

    SystemContext.MULTIMERIC_COMPLEX: ContextProfile(
        name="Multimeric Complex",
        description="Oligomeric protein assembly (dimer, tetramer, etc.).",
        rmsd_nm=SoftRange(
            excellent_low=0.0, excellent_high=0.25,
            acceptable_low=0.0, acceptable_high=0.40,
            warning_low=0.0,   warning_high=0.60,
            unit="nm",
            note="Global RMSD; inter-subunit RMSD often more informative",
        ),
        rg_nm=SoftRange(
            excellent_low=2.0, excellent_high=6.0,
            acceptable_low=1.5, acceptable_high=8.0,
            warning_low=1.0,   warning_high=12.0,
            unit="nm",
        ),
        temperature_k=SoftRange(
            excellent_low=295.0, excellent_high=305.0,
            acceptable_low=290.0, acceptable_high=310.0,
            warning_low=280.0,   warning_high=320.0,
            unit="K",
        ),
        pressure_bar=SoftRange(
            excellent_low=-50.0, excellent_high=50.0,
            acceptable_low=-200.0, acceptable_high=200.0,
            warning_low=-500.0,  warning_high=500.0,
            unit="bar",
        ),
        energy_stability_kj=SoftRange(
            excellent_low=0.0, excellent_high=0.6,
            acceptable_low=0.0, acceptable_high=1.8,
            warning_low=0.0,   warning_high=4.0,
            unit="kJ/mol/atom",
        ),
        notes=[
            "Interface RMSD and buried surface area are key stability metrics.",
            "Rg drift may indicate dissociation — monitor subunit separation distances.",
        ],
    ),

    SystemContext.ENZYME: ContextProfile(
        name="Enzyme",
        description="Catalytic protein; active site geometry is primary quantity of interest.",
        rmsd_nm=SoftRange(
            excellent_low=0.0, excellent_high=0.18,
            acceptable_low=0.0, acceptable_high=0.28,
            warning_low=0.0,   warning_high=0.40,
            unit="nm",
            note="Active-site residue RMSD is more critical than global RMSD",
        ),
        rg_nm=SoftRange(
            excellent_low=1.0, excellent_high=2.5,
            acceptable_low=0.8, acceptable_high=3.0,
            warning_low=0.5,   warning_high=4.0,
            unit="nm",
        ),
        temperature_k=SoftRange(
            excellent_low=295.0, excellent_high=305.0,
            acceptable_low=290.0, acceptable_high=310.0,
            warning_low=280.0,   warning_high=320.0,
            unit="K",
        ),
        pressure_bar=SoftRange(
            excellent_low=-50.0, excellent_high=50.0,
            acceptable_low=-200.0, acceptable_high=200.0,
            warning_low=-500.0,  warning_high=500.0,
            unit="bar",
        ),
        energy_stability_kj=SoftRange(
            excellent_low=0.0, excellent_high=0.5,
            acceptable_low=0.0, acceptable_high=1.5,
            warning_low=0.0,   warning_high=3.0,
            unit="kJ/mol/atom",
        ),
        notes=[
            "Active site RMSD and catalytic residue distances are primary metrics.",
            "H-bond network in active site should be monitored for occupancy.",
        ],
    ),

    SystemContext.PEPTIDE: ContextProfile(
        name="Peptide",
        description="Short peptide (< 50 residues); may be structured or unstructured.",
        rmsd_nm=SoftRange(
            excellent_low=0.0, excellent_high=0.30,
            acceptable_low=0.0, acceptable_high=0.60,
            warning_low=0.0,   warning_high=1.20,
            unit="nm",
            note="Short peptides are highly flexible; absolute RMSD less meaningful",
        ),
        rg_nm=SoftRange(
            excellent_low=0.3, excellent_high=1.5,
            acceptable_low=0.2, acceptable_high=2.0,
            warning_low=0.1,   warning_high=3.0,
            unit="nm",
        ),
        temperature_k=SoftRange(
            excellent_low=295.0, excellent_high=305.0,
            acceptable_low=290.0, acceptable_high=310.0,
            warning_low=280.0,   warning_high=320.0,
            unit="K",
        ),
        pressure_bar=SoftRange(
            excellent_low=-50.0, excellent_high=50.0,
            acceptable_low=-200.0, acceptable_high=200.0,
            warning_low=-500.0,  warning_high=500.0,
            unit="bar",
        ),
        energy_stability_kj=SoftRange(
            excellent_low=0.0, excellent_high=1.0,
            acceptable_low=0.0, acceptable_high=3.0,
            warning_low=0.0,   warning_high=6.0,
            unit="kJ/mol/atom",
        ),
        notes=[
            "Enhanced sampling recommended for small peptides in explicit solvent.",
            "Secondary structure content (helix fraction) often more informative than RMSD.",
        ],
    ),

    SystemContext.MEMBRANE_SYSTEM: ContextProfile(
        name="Pure Membrane System",
        description="Lipid bilayer without embedded protein; bilayer properties are primary.",
        rmsd_nm=SoftRange(
            excellent_low=0.0, excellent_high=1.0,
            acceptable_low=0.0, acceptable_high=2.0,
            warning_low=0.0,   warning_high=4.0,
            unit="nm",
            note="Lipid RMSD is not a useful metric; use APL, bilayer thickness instead",
        ),
        rg_nm=SoftRange(
            excellent_low=2.0, excellent_high=8.0,
            acceptable_low=1.0, acceptable_high=12.0,
            warning_low=0.5,   warning_high=20.0,
            unit="nm",
        ),
        temperature_k=SoftRange(
            excellent_low=295.0, excellent_high=310.0,
            acceptable_low=285.0, acceptable_high=325.0,
            warning_low=270.0,   warning_high=340.0,
            unit="K",
        ),
        pressure_bar=SoftRange(
            excellent_low=-50.0, excellent_high=50.0,
            acceptable_low=-200.0, acceptable_high=200.0,
            warning_low=-500.0,  warning_high=500.0,
            unit="bar",
        ),
        energy_stability_kj=SoftRange(
            excellent_low=0.0, excellent_high=1.0,
            acceptable_low=0.0, acceptable_high=2.5,
            warning_low=0.0,   warning_high=5.0,
            unit="kJ/mol/atom",
        ),
        notes=[
            "Area per lipid (APL) convergence is the primary equilibration indicator.",
            "DPPC at 323 K: APL ~64 Å², bilayer thickness ~40 Å.",
            "Equilibration for bilayer systems typically requires 20–100 ns.",
        ],
    ),

    SystemContext.FLEXIBLE_DOMAIN_PROTEIN: ContextProfile(
        name="Flexible Domain Protein",
        description="Protein with rigid core and one or more flexible loops/domains.",
        rmsd_nm=SoftRange(
            excellent_low=0.0, excellent_high=0.35,
            acceptable_low=0.0, acceptable_high=0.55,
            warning_low=0.0,   warning_high=0.80,
            unit="nm",
            note="Separate core and flexible domain RMSD for accurate convergence assessment",
        ),
        rg_nm=SoftRange(
            excellent_low=1.0, excellent_high=3.0,
            acceptable_low=0.8, acceptable_high=4.0,
            warning_low=0.5,   warning_high=6.0,
            unit="nm",
        ),
        temperature_k=SoftRange(
            excellent_low=295.0, excellent_high=305.0,
            acceptable_low=290.0, acceptable_high=310.0,
            warning_low=280.0,   warning_high=320.0,
            unit="K",
        ),
        pressure_bar=SoftRange(
            excellent_low=-50.0, excellent_high=50.0,
            acceptable_low=-200.0, acceptable_high=200.0,
            warning_low=-500.0,  warning_high=500.0,
            unit="bar",
        ),
        energy_stability_kj=SoftRange(
            excellent_low=0.0, excellent_high=0.7,
            acceptable_low=0.0, acceptable_high=2.0,
            warning_low=0.0,   warning_high=4.0,
            unit="kJ/mol/atom",
        ),
        notes=[
            "High RMSF in flexible regions is expected and physically meaningful.",
            "Core RMSD convergence is sufficient for thermodynamic analysis.",
        ],
    ),

    SystemContext.UNKNOWN: ContextProfile(
        name="Unknown System",
        description="System context not specified; use conservative generic ranges.",
        rmsd_nm=SoftRange(
            excellent_low=0.0, excellent_high=0.25,
            acceptable_low=0.0, acceptable_high=0.45,
            warning_low=0.0,   warning_high=0.70,
            unit="nm",
        ),
        rg_nm=SoftRange(
            excellent_low=0.5, excellent_high=4.0,
            acceptable_low=0.3, acceptable_high=6.0,
            warning_low=0.1,   warning_high=10.0,
            unit="nm",
        ),
        temperature_k=SoftRange(
            excellent_low=295.0, excellent_high=305.0,
            acceptable_low=290.0, acceptable_high=310.0,
            warning_low=280.0,   warning_high=320.0,
            unit="K",
        ),
        pressure_bar=SoftRange(
            excellent_low=-50.0, excellent_high=50.0,
            acceptable_low=-200.0, acceptable_high=200.0,
            warning_low=-500.0,  warning_high=500.0,
            unit="bar",
        ),
        energy_stability_kj=SoftRange(
            excellent_low=0.0, excellent_high=0.8,
            acceptable_low=0.0, acceptable_high=2.0,
            warning_low=0.0,   warning_high=4.0,
            unit="kJ/mol/atom",
        ),
    ),
}
