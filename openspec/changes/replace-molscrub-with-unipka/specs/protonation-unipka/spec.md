## Purpose

Generates protomer candidates with Uni-Pka: batched container execution, occupancy-based candidate selection at the working pH, storage of predicted protonation properties, and a stored pH-distribution artifact.

## ADDED Requirements

### Requirement: Batched container execution

The system SHALL execute Uni-Pka through the configured container runtime (Docker or Apptainer) using the configured image or `.sif` path, submitting all input molecules of the protonation stage in a single invocation per run, at the configured working pH.

#### Scenario: All molecules in one call

- **WHEN** the protonation stage processes the run's input molecules with `protonation.tool: unipka`
- **THEN** Uni-Pka is invoked at most once with a tab-separated SMILES+name input covering all accepted molecules, and results are matched back by molecule name

#### Scenario: Container invocation failure

- **WHEN** the container exits non-zero or produces no parsable output
- **THEN** the protonation stage fails with the tool error surfaced, and no protomer records claim Uni-Pka provenance

#### Scenario: Molecule Uni-Pka cannot protonate

- **WHEN** Uni-Pka returns the input SMILES with no occupancy (`NA`) for a molecule
- **THEN** the input molecule's own state is retained as its sole protomer with a warning recorded, consistent with the existing molscrub fallback behavior

### Requirement: Occupancy-based protomer selection

The system SHALL request multiple protonation forms per molecule (up to a configurable maximum with a configurable minimum-occupancy threshold at the working pH) and SHALL select protomers by predicted occupancy: forms meeting the occupancy threshold are retained, ordered by decreasing occupancy, and trimmed to `max_protomers_per_molecule`; the `keep_best_per_charge` rule SHALL apply using each form's formal charge. The heuristic plausibility scoring used for molscrub candidates SHALL NOT be applied to Uni-Pka candidates.

#### Scenario: Dominant form selected

- **WHEN** a molecule's protonation forms have occupancies 0.97 and 0.03 with a threshold of 0.05
- **THEN** only the 0.97 form proceeds as a protomer

#### Scenario: Several forms above threshold within cap

- **WHEN** a molecule has forms with occupancies 0.60, 0.35, and 0.05, threshold 0.05, and `max_protomers_per_molecule` is 4
- **THEN** all three forms proceed as protomers ordered by decreasing occupancy

#### Scenario: Forms exceed cap

- **WHEN** more forms meet the threshold than `max_protomers_per_molecule`
- **THEN** the highest-occupancy forms up to the cap are retained and the trimming is recorded

#### Scenario: Keep input state

- **WHEN** `protonation.keep_input_state` is true and the input protonation state is not among the selected forms
- **THEN** the input state is retained in addition, subject to the cap, as in current behavior

### Requirement: Per-protomer predicted properties in metadata

Each Uni-Pka protomer record SHALL store in its metadata the predicted occupancy at the working pH and the predicted microstate free energy (dG).

#### Scenario: Metadata fields present

- **WHEN** a Uni-Pka protomer record is written to stage outputs
- **THEN** its metadata contains `unipka_occupancy` (fraction in [0,1]) and `unipka_dg` (free energy value)

### Requirement: Per-molecule protonation summary properties

For each input molecule processed by Uni-Pka, the system SHALL compute and store a protonation summary containing: the top-two occupancy gap at the working pH, the occupancy-distribution entropy at the working pH, the net-charge population distribution at the working pH, the number of microstates enumerated, the distance from the working pH to the nearest macroscopic pKa transition, and the isoelectric point.

#### Scenario: Summary written for every Uni-Pka molecule

- **WHEN** the protonation stage completes for a molecule via Uni-Pka
- **THEN** the molecule's stage output contains all six summary properties, with values derived only from the returned ensemble free energies and occupancies

#### Scenario: Degenerate ensembles

- **WHEN** a molecule retains only its input state (no predicted forms)
- **THEN** the summary is written with null values and a warning rather than failing the run

### Requirement: pH-distribution artifact

The system SHALL store Uni-Pka's microspecies pH-distribution output (per microspecies: input SMILES, microspecies SMILES, dG, occupancy, pH over the configured pH range) as a run artifact under the protonation stage directory, SHALL cover the working pH in the range, and SHALL reference the artifact from the run report.

#### Scenario: Artifact written and referenced

- **WHEN** a Uni-Pka run completes
- **THEN** the distribution file exists in the protonation stage directory for every molecule with predictions and the run report links it

### Requirement: Protonation property reporting

The run report SHALL present per-molecule protonation properties (selected-form occupancies and the summary properties) so they are usable for downstream uncertainty analysis.

#### Scenario: Report shows protonation evidence

- **WHEN** a Uni-Pka run's markdown report is generated
- **THEN** each molecule's section lists its protomers with occupancies and the summary values including pI and nearest-pKa distance
