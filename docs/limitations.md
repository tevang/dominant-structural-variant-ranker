# Limitations

DSVR is an orchestration layer for practical structural-variant preparation and ranking. It can coordinate useful open-source tools, but it does not make their results more rigorous than the underlying enumeration coverage, physics, scoring model, and thermodynamic corrections support.

## Default Scope

The default workflow is designed for fast ligand preparation, docking, ligand-based modeling, and batch-library preparation. It is not an exhaustive conformational free-energy workflow.

## pH and Protonation

Default pH handling is candidate generation plus model predictions, not
rigorous constant-pH thermodynamics. The default Uni-Pka path enumerates
microspecies and predicts a free energy per microspecies; protomer selection
uses the predicted Boltzmann occupancies at the working pH, and six
per-molecule summary properties (ensemble entropy, charge populations,
microstate count, top-two occupancy gap, isoelectric point, distance to the
nearest macroscopic pKa transition) are stored for later uncertainty
analysis.

Uni-Pka predictions come with sharp caveats: the model was trained on
drug-like chemistry (Uni-Mol/Dwars) and its absolute occupancies are
estimates; the pH-distribution grid for pKa/pI/interpolation is coarse for
publication values (default step 0.25, analysis refined to ≤0.05); pure
acids and pure bases have a null isoelectric point by construction (net
charge never crosses zero); and the downstream ranking in this release does
not re-weight variants by Uni-Pka occupancies — they are recorded features,
not corrections. molscrub (`protonation.tool: molscrub`) remains available
and uses the original heuristic plausibility scoring with no occupancy
predictions.

Because Uni-Pka selection differs from molscrub's, default runs are not
comparable to molscrub-era runs; provenance records which tool produced
each run.

DSVR does not perform rigorous pH-dependent population calculations unless a micro-pKa/proton chemical-potential correction provider is added. Without that correction, cross-protomer and cross-charge populations are approximate.

## Tautomers

RDKit tautomer canonicalization and enumeration are not stability or abundance ranking. RDKit alone can enumerate too many tautomers and does not identify the dominant tautomer in solution.

The default workflow uses Auto3D to rank low-energy tautomer candidates by optimized conformer energies, but this is approximate potential-energy triage, not true solution abundance. Tautomer filtering happens before stereoisomer enumeration to avoid multiplying every tautomer by every possible stereoisomer.

## Stereoisomers

RDKit stereoisomer enumeration is explicit and controlled by timeout and caps. The workflow can expand undefined stereocenters, but generated candidates are only as complete as input chemistry and configured limits allow.

In achiral solvent, DSVR may collapse enantiomer-equivalent work and run energy calculations for representative enantiomers only. This is not valid for chiral binding pockets, chiral chromatography, chiral catalysts, or other chiral environments.

## Auto3D Ranking and Thermodynamics

Auto3D is useful for tautomer/stereoisomer energy triage and one-conformer final 3D optimization. Auto3D energies are optimized-conformer potential-energy signals. They are not pKa predictions, not solvent speciation, and not rigorous abundance estimates.

Auto3D thermodynamics, if enabled, are not substitutes for validated solvated free energies. Treat them as screening signals unless an explicit validation protocol is run.

## Optional Physics-Heavy Validation

CREST/xTB conformer searches, xTB thermo, CREST entropy estimates, and CENSO are optional validation/refinement steps for selected candidates. They are expensive and are not part of the default ligand-preparation path.

Psi4/PySCF rescoring is outside the default workflow and should remain an optional legacy/advanced module unless explicitly enabled.

## Molecular Identity and Deduplication

"Same structure" is decided by a toolkit-local exact-duplicate key — molecular
formula, net formal charge, and canonical isomeric SMILES after RDKit
`rdMolStandardize.Cleanup` (ChEMBL-style normalization; never uncharged, never
a tautomer or charge parent). The key is valid only inside this RDKit-based
pipeline under this policy; it is not a global identifier and does not
interoperate with other toolkits' canonicalization.

Standard InChI/InChIKey is deliberately not used as a dedupe key: its
normalization collapses mobile-hydrogen tautomers into one identifier, which
would over-merge distinct tautomers at exactly the stages that must keep them
apart. IUPAC names and raw, unstandardized SMILES strings are likewise never
used for identity.

Exact-duplicate elimination runs after every expansion step (tautomer
selection, stereoisomer selection, and a final safety net before
`final_variants` is written), because dedupe that is only branch-local leaves
duplicates when protomer branches converge on the same structures (for
example protomers that are tautomers of one another). After dedupe, branches
are refilled to their configured caps only from already-ranked unused unique
candidates: the refill pool is bounded (stereoisomer enumeration is
over-enumerated by `stereoisomer_filtering.enumeration_ceiling_multiplier`,
default 4), so a branch can remain below cap when its candidate space is
exhausted; that shortfall is recorded rather than padded with duplicates.
Resumed runs re-execute the tautomer-and-later stages because the deduped
record sets change downstream input hashes; refill pools are in-memory only,
so resume-loaded branches can dedupe but cannot refill. Rank metadata is not
persisted to the per-branch SDFs, so on a tautomers-step resume the
representative of a duplicate group falls back to record-ID order instead of
best-ranked provenance; the final structures still converge because any
representative has the same identity key and the final safety net holds.

Per-branch artifacts (`enumeration/tautomers/*_tautomers.sdf`,
`enumeration/stereoisomers/*_stereoisomers.sdf`) are branch-local snapshots
written before cross-branch deduplication; the deduplicated view lives in
`all_tautomers.sdf`, `all_stereoisomers.sdf`, the `*_dedupe.csv` audit files,
and downstream artifacts.

## Enumeration Bias

Missing candidate states cannot be recovered by downstream ranking. If molscrub, RDKit, or Auto3D omits a relevant protomer, tautomer, stereoisomer, or conformer, DSVR can only rank the candidates it receives.

## Dependency and Version Sensitivity

Results may change with tool versions, default parameters, CPU/GPU settings, threading, random seeds, and parser behavior. DSVR should record versions, command lines, configuration, input hashes, and logs for every run.
