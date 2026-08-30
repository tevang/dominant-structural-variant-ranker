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

## Duplicate Elimination and Fill-to-Cap Semantics

Each expansion step (protomers per input, tautomers per protomer,
stereoisomers per tautomer, one conformer per variant) is deduplicated
branch-locally, and additionally every step is followed by a cross-branch
exact-duplicate guard per input molecule: identical structures that arise from
different protomer branches (for example protomers that are tautomers of one
another, such as neutral azoles protonated on different ring nitrogens) are
collapsed to one representative.

"Same structure" is decided by a toolkit-local exact-duplicate key: molecular
formula, net formal charge, and canonical isomeric SMILES computed after RDKit
`rdMolStandardize.Cleanup` (ChEMBL-style standardization; never uncharged and
never tautomer/charge-parent). Standard InChI/InChIKey is deliberately **not**
used, because its normalization collapses mobile-hydrogen tautomers into one
identifier and would over-merge distinct tautomers; raw unstandardized SMILES
is likewise not used, because Kekulé or charge-representation differences
would under-merge. The key is valid only within this RDKit-based pipeline and
this stated policy; it is not a global identifier.

When deduplication drops a branch below its configured cap, the branch is
refilled to the cap from its own next-best unused unique ranked candidates
(`tauto_k` per protomer branch for tautomers, `max_stereoisomers_per_tautomer`
for stereoisomers) — fill-to-cap semantics — and a final safety-net dedupe
(lowest final energy wins) runs before `final_variants.sdf/csv/json` are
written. Eliminated duplicates, refill promotions, and unfillable shortfalls
are recorded in `enumeration/tautomers/tautomer_dedupe.csv`,
`stereoisomer_filtering/stereo_dedupe.csv`, and `final_dedupe_audit.csv`, and
the retained record carries `merged_from` provenance.

Because deduplication changes the record sets that downstream steps consume,
resuming a run from before this policy was introduced re-executes the tautomer
and later stages automatically (their input hashes change). Within a current
run, resume re-derives the same dedupe decisions from the per-branch SDF
outputs; ranked refill pools only exist in memory, so resume-loaded branches
can dedupe but cannot refill (shortfalls are recorded instead).

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
