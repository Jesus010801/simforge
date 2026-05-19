SYSTEM_TYPES = [

    "protein",
    "protein-ligand",
    "protein-membrane",
    "protein-membrane-ligand",

    "enzyme-substrate",
    "competitive-inhibition",
    "allosteric-modulation",

    "peptide",
    "peptide-membrane",

    "small-molecule",

    "multicomponent-system",
    "essential-oil-system",

    "coarse-grained",

]

SIMULATION_GOALS = [

    "stability",
    "binding",
    "competitive_binding",

    "membrane_insertion",
    "membrane_perturbation",

    "conformational_sampling",

    "aggregation",

    "permeability",

    "allosteric_effect",

    "active_site_dynamics",

    "active_site_stability",
]

COMPONENT_ROLES = [

    "protein",

    "substrate",

    "competitive_ligand",

    "allosteric_ligand",

    "cofactor",

    "membrane",

    "solvent",

    "ion",

    "lipid",

    "peptide",

    "essential_oil_component",

]

FORCEFIELDS = {

    "atomistic": [

        "charmm36",
        "amber99sb",
        "amber14",
        "opls-aa",

    ],

    "ligands": [

        "gaff",
        "cgenff",
        "openff",

    ],

    "coarse_grained": [

        "martini",

    ]
}

WATER_MODELS = [

    "tip3p",
    "tip4p",
    "spce",

]

MEMBRANE_TYPES = [

    "POPC",
    "POPE",
    "DLPC",
    "mixed",
    "custom",

]

BOX_TYPES = [

    "cubic",
    "triclinic",
    "dodecahedron",
    "rectangular",

]

BIOLOGICAL_CONTEXT = [

    "soluble",

    "membrane_associated",

    "transmembrane",

    "partially_truncated",

    "oligomeric",

    "crystal_structure",

    "experimental_mutation",

]

RESTRAINT_TYPES = [

    "position_restraints",

    "terminal_restraints",

    "backbone_restraints",

    "distance_restraints",

]

ANALYSIS_TYPES = [

    "rmsd",

    "rmsf",

    "radius_of_gyration",

    "sasa",

    "hydrogen_bonds",

    "distance_analysis",

    "contact_map",

    "secondary_structure",

    "binding_energy",

    "membrane_thickness",

    "density_profile",

]

OUTPUT_FORMATS = [

    "png",
    "svg",
    "pdf",
    "csv",

]
