## Purpose

Ensures enumeration stages produce unique candidates per input molecule: duplicates that arise across protomer branches are eliminated after each expansion step, and stage caps are refilled from remaining unique candidates so variant counts stay at the configured limits.

## ADDED Requirements

### Requirement: Cross-branch tautomer deduplication

After tautomer selection for an input molecule, the workflow SHALL eliminate exact-duplicate tautomers across all protomer branches of that input molecule (using the molecular-identity exact key), retaining one representative per unique structure. This applies to both tautomer generation paths (Auto3D-filtered and plain RDKit enumeration).

#### Scenario: Protomers converge on the same tautomer

- **WHEN** two or more protomer branches of one input molecule select the same exact tautomer structure
- **THEN** exactly one tautomer record per unique structure proceeds to stereoisomer enumeration and the others are recorded as merged duplicates

#### Scenario: Representative keeps best ranking evidence

- **WHEN** a duplicate group contains records with different ranking evidence (e.g. Auto3D relative energies)
- **THEN** the retained representative is the record with the best available rank evidence, with deterministic tie-breaking

### Requirement: Tautomer refill to cap

When cross-branch deduplication leaves a protomer branch with fewer selected tautomers than `tauto_k`, the workflow SHALL promote that branch's next-best ranked, still-unique candidates until the branch reaches `tauto_k` or its ranked candidate pool is exhausted; exhaustion SHALL be recorded.

#### Scenario: Refill after duplicate elimination

- **WHEN** deduplication drops a branch below `tauto_k` and the branch has ranked candidates not already selected anywhere for that input molecule
- **THEN** the next-best unique candidates are promoted into the branch's selected set up to `tauto_k`

#### Scenario: Refill pool exhausted

- **WHEN** deduplication drops a branch below `tauto_k` and no unused unique ranked candidates remain
- **THEN** the branch proceeds with fewer than `tauto_k` tautomers and the shortfall is recorded in run artifacts

### Requirement: Stereoisomer deduplication and refill to cap

Stereoisomer enumeration SHALL deduplicate within each tautomer and SHALL guard against exact duplicates across tautomers of one input molecule. If deduplication leaves a tautomer below `max_stereoisomers_per_tautomer`, the workflow SHALL enumerate additional candidates up to a bounded internal ceiling and select the first unique candidates until the cap is reached or the enumeration space is exhausted.

#### Scenario: Unique stereoisomers per tautomer up to cap

- **WHEN** a tautomer's enumeration yields duplicates and additional distinct stereoisomers exist within the bounded ceiling
- **THEN** the tautomer proceeds with exactly `max_stereoisomers_per_tautomer` unique stereoisomers, or fewer with the shortfall recorded

#### Scenario: Cross-tautomer stereoisomer guard

- **WHEN** two tautomers of one input molecule would produce the same exact stereoisomer structure
- **THEN** only one copy proceeds downstream and the elimination is audited

### Requirement: Final output contains no exact duplicates

Before writing `final_variants.sdf`, `final_variants.csv`, and `final_variants.json`, the workflow SHALL remove exact-duplicate final records per input molecule, retaining the record with the lowest available final energy and deterministic tie-breaking.

#### Scenario: Clean final variant files

- **WHEN** final 3D generation converges to identical structures from different branches (identical standardized structure and geometry)
- **THEN** each `final_variants` artifact contains exactly one record per unique structure per input molecule

#### Scenario: Protomers as tautomers of each other (regression)

- **WHEN** an input molecule's protomer candidates are tautomers of one another (e.g. a neutral molecule with multiple azole nitrogens) and their tautomer sets overlap
- **THEN** the final outputs contain no duplicate structures and every eliminated duplicate appears in the run's dedupe audit artifacts

### Requirement: Protomer selection keeps fill-to-cap behavior

Protomer selection SHALL continue to fill up to `max_protomers_per_molecule` from the unique candidate pool after protomer deduplication; if the enumeration backend returns fewer unique candidates than the cap, the shortfall SHALL be recorded rather than padded with duplicates.

#### Scenario: Protomer pool smaller than cap

- **WHEN** the protonation backend yields fewer unique candidates than `max_protomers_per_molecule`
- **THEN** all unique candidates are retained and the shortfall is visible in run artifacts
