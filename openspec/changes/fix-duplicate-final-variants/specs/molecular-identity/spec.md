## Purpose

Defines how the workflow decides that two structures are the same structure, using keys derived from standardized structure objects rather than raw strings, so duplicate handling is reproducible and does not silently over- or under-merge chemical variants.

## ADDED Requirements

### Requirement: Exact-duplicate key from standardized structures

The workflow SHALL derive exact-duplicate keys from RDKit structure objects after a documented, reproducible standardization pass, and the key SHALL consist of molecular formula, net formal charge, and canonical isomeric SMILES of the standardized structure. Keys SHALL be documented as toolkit-local (valid within this RDKit-based pipeline and its stated policy), not as global identifiers.

#### Scenario: Same structure, different input representation

- **WHEN** two records represent the same exact structure via different valid representations (e.g. different Kekulé forms or equivalent charge-separated drawings) from different enumeration branches
- **THEN** both records yield the same exact-duplicate key

#### Scenario: Distinct structures are not merged

- **WHEN** two records differ in any of molecular formula, net formal charge, or standardized structure (including stereochemistry)
- **THEN** they yield different exact-duplicate keys and both are retained

### Requirement: Standardization preserves chemically meaningful state

The standardization used for key derivation MUST NOT neutralize charges, strip fragments, remove stereochemistry, or collapse tautomers: distinct protonation states, distinct tautomers, and distinct stereoisomers SHALL remain distinct keys.

#### Scenario: Tautomers stay distinct

- **WHEN** two records are different tautomers of the same input molecule (same formula and charge, different proton placement)
- **THEN** they yield different exact-duplicate keys

#### Scenario: Ionization states stay distinct

- **WHEN** two records differ only in formal charge state
- **THEN** they yield different exact-duplicate keys

### Requirement: Forbidden dedupe keys

The workflow MUST NOT use Standard InChI or Standard InChIKey as an exact-duplicate key (they normalize mobile-hydrogen tautomers into one identifier), MUST NOT compare raw, unstandardized SMILES strings for identity, and MUST NOT use IUPAC names for identity. Full structure records SHALL remain the provenance of record for every retained and merged variant.

#### Scenario: Tautomer-specific comparison survives InChI collapsing

- **WHEN** deduplication is applied at a stage whose candidates are distinct tautomers
- **THEN** the comparison uses the standardized exact key, not Standard InChI, so distinct tautomers are never merged by the identity layer

### Requirement: Merged provenance for deduplicated records

Whenever duplicates are eliminated, the retained record SHALL carry metadata identifying every merged contributing record (parent chain and record IDs), and the elimination SHALL be written to a machine-readable audit artifact of the run.

#### Scenario: Duplicate collapse is auditable

- **WHEN** two or more records from different enumeration branches are merged as exact duplicates
- **THEN** the retained record lists all contributing variant IDs and the run contains an audit entry naming retained and eliminated IDs
