## Why

The `prepare-ligands` workflow emits duplicate structures in `final_variants.sdf`. In run `Vaclav_set1_p8_t10_s24`, 13 of 81 final structures are redundant copies (11 duplicate groups across 7 of 10 inputs): `duplicate_structures_report.txt` shows entries with identical atom types, 3D coordinates, and bond connectivity, differing only in the protomer index of the variant ID.

Root cause: molscrub protomer candidates can be tautomers of one another (e.g. neutral azoles protonated on different ring nitrogens pass the protomer dedupe key because their canonical SMILES differ legitimately). Each protomer branch then enumerates its own tautomer set, and these sets overlap — for `mol_000007`, 21 selected tautomer records covered only 6 unique structures. Deduplication is branch-local at every stage (protomers per input, tautomers per protomer, stereoisomers per tautomer, one conformer per parent), so nothing collapses duplicates across protomer branches. The redundancy also multiplies downstream Auto3D compute (~3.5× for the affected molecule).

Duplicate detection must not rely on raw string comparison: SMILES equivalence is only defined relative to a toolkit and a standardization policy, and Standard InChI deliberately collapses mobile-hydrogen tautomers (over-merge) while unstandardized SMILES risks under-merge from Kekulé/charge-representation differences.

## What Changes

- Introduce a documented molecular-identity policy: an exact-duplicate key derived from RDKit-standardized structure objects (formula, net formal charge, canonical isomeric SMILES after `rdMolStandardize.Cleanup` without uncharging). Standard InChI is explicitly not used as a dedupe key (tautomer-collapsing); keys are documented as toolkit-local.
- Eliminate cross-protomer exact-duplicate tautomers per input molecule after tautomer selection (covers both the Auto3D tautomer filter and the RDKit enumeration path), keeping the best-ranked representative and recording merged provenance plus an audit trail.
- Add a cross-tautomer duplicate guard at stereoisomer stage and a final content dedupe before writing `final_variants.sdf/csv/json`, keeping the lowest-energy record with merged provenance.
- Adopt fill-to-cap semantics at enumeration stages: after deduplication, refill each protomer branch's tautomers to `tauto_k` from its next-best unique unused candidates, and refill stereoisomer sets to `max_stereoisomers_per_tautomer` by enumerating against a bounded internal ceiling, until the cap is reached or the candidate pool is exhausted (exhaustion is recorded).
- Add regression tests reproducing the mol_000007 duplicate case and update workflow/limitations docs to describe the identity policy.

## Capabilities

### New Capabilities

- `molecular-identity`: the standardized exact-duplicate key policy — what identifies two structures as the same structure, which representations are forbidden as keys and why, and how merged provenance is recorded.
- `enumeration-deduplication`: cross-branch duplicate elimination after each enumeration expansion step, and fill-to-cap semantics for protomer, tautomer, and stereoisomer stages.

### Modified Capabilities

## Impact

- Affected code: new `dsvr/chemistry/identity.py`; `dsvr/workflow/engine.py` (post-tautomer cross-branch dedupe and refill wiring); `dsvr/chemistry/tautomer_auto3d_filter.py` and `dsvr/chemistry/tautomers.py` (ranked candidate pool exposure for refill); `dsvr/chemistry/stereochemistry.py` (bounded re-enumeration and refill); `dsvr/chemistry/stereo_auto3d_filter.py` (cross-tautomer guard); `dsvr/chemistry/final3d.py` (final safety-net dedupe); tests and docs.
- Outputs: `final_variants.sdf/csv/json` contain no exact duplicates per input molecule; new audit artifacts record merged/eliminated duplicates and cap-refill events; resumed runs re-execute tautomer-and-later stages because stage inputs change.
- No breaking changes to CLI or config schema; optional config knobs may be added for the stereo enumeration ceiling and for disabling the final safety-net dedupe.
