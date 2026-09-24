# Policy Engine Specification

## Purpose

Policies convert evidence-bearing candidates into deterministic workflow-scoped decisions. They do not rewrite global scientific knowledge.

## Decision states

### RESOLVED
Preferred acceptable solution found.

### RESOLVED_WITH_WARNING
Valid fallback exists but limitations must propagate into DecisionContext.

### DEFERRED
Resolution is delegated to an explicitly defined local scientific computation at a later preparation stage.

`DEFERRED` is not "unknown". It must define target stage, action, method/policy, required inputs, success criteria and failure behavior.

### ABORT
No scientifically acceptable path exists under active policy.

## Structure selection example

```yaml
policy_name: high_fidelity_structure_selection
version: 1
rules:
  - action: REQUIRE
    predicates:
      structure_type: experimental
      resolution_angstrom_lte: 2.5
  - action: FALLBACK
    predicates:
      structure_type: experimental
      resolution_angstrom_lte: 3.0
  - action: FALLBACK
    predicates:
      structure_type: predicted
      plddt_gte: 90.0
  - action: ABORT
    message: No suitable structure found.
```

Attributes such as conformational state should carry evidence/confidence requirements rather than being treated as unquestionable booleans.

## Protonation policy

```yaml
policy_name: physiological_protonation
version: 1
preserve:
  heavy_atom_connectivity: true
  stereochemistry: true
determine_locally:
  hydrogens: true
  formal_protonation_state: true
  tautomer_if_ambiguous: true
environment:
  pH: 7.4
  ionic_strength_M: 0.15
on_failure: ABORT
```

The local calculation is recorded in the RunDerivationLedger.

## Membrane policy

A generic template is a fallback model, not biological fact. Policy may prefer:

1. source-supported composition
2. organism/tissue-specific model
3. broader compartment template
4. generic fallback
5. abort if specificity is scientifically required

Warnings/confidence propagate.

## Implementation guidance

MVP policies should initially be typed Python objects. Do not build a general-purpose policy DSL until semantics are stable and tested.
