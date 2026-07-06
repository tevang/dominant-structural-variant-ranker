# Plausible Variant Workflow

DSVR's default workflow is a LigPrep-like ligand-preparation protocol for docking,
ligand-based modeling, and batch-library preparation. It is designed to produce a
bounded set of plausible pH, tautomer, stereochemical, and 3D ligand variants
quickly. It is not an exhaustive conformational free-energy workflow.

```text
Input SMILES/SDF
-> standardization and validity checks
-> plausible pH/protomer generation at target pH, default 7.0
-> early protomer filtering
-> Auto3D tautomer enumeration/ranking/filtering using RDKit tautomer engine and ANI2xt/AIMNet2
-> RDKit stereoisomer enumeration with timeouts/caps after tautomer filtering
-> Auto3D one-conformer optimization/ranking/filtering of stereoisomers
-> final SDF/CSV/JSON report with one optimized 3D conformer per surviving structural variant
-> optional CREST/xTB validation only if explicitly enabled
```

## Why Tautomers Are Filtered Before Stereoisomers

RDKit can enumerate many tautomer candidates, but RDKit tautomer enumeration is
candidate generation, not abundance ranking. Expanding stereoisomers for every
tautomer multiplies the candidate count and can create a combinatorial explosion
before any energy triage happens.

The default workflow therefore ranks and filters tautomer candidates before
stereoisomer enumeration. This keeps the downstream stereochemistry and 3D
optimization steps focused on plausible low-energy tautomer candidates.

## Energy Ranking Scope

Auto3D ranking is approximate potential-energy triage. It uses optimized
conformer energies to prioritize low-energy tautomers and stereoisomers, but it
does not predict true solution abundance, pKa, solvent speciation, or rigorous
free energies.

If Auto3D thermodynamic outputs are used, they still are not substitutes for
validated solvated free energies. Treat them as screening signals unless an
explicit validation protocol is run.

## Optional Validation

CREST/xTB conformer searches, xTB thermo, CREST entropy estimates, CENSO, and
Psi4/PySCF rescoring are optional validation/refinement paths. They are useful
for selected small candidate sets, but they are expensive and are not part of the
default ligand-preparation workflow.
