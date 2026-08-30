# auto3d-engine-selection

## Purpose

Selects an Auto3D engine compatible with each molecule's elemental composition, splits batches when molecules need different engines, and never retries an engine already proven incompatible for a molecule.

## ADDED Requirements

### Requirement: Compatibility-checked engine selection
Before invoking Auto3D, the system SHALL verify that the selected engine supports the elements present in every molecule of the batch. Molecules unsupported by the primary engine SHALL be routed to a configured engine that supports them, or the batch SHALL be split so each sub-batch uses a supporting engine.

#### Scenario: Mixed batch needs AIMNET
- **WHEN** a final-3D batch configured for ANI2xt contains molecules with elements only AIMNET supports
- **THEN** those molecules are computed with AIMNET and the rest with ANI2xt, without "Only AIMNET can handle" retry loops

#### Scenario: No configured engine supports a molecule
- **WHEN** neither the primary nor fallback engine supports a molecule's elements
- **THEN** the molecule is recorded as an engine-incompatibility failure and the RDKit fallback is used for it, without invoking Auto3D with an incompatible engine

### Requirement: No retry of proven-incompatible engines
An engine rejected by Auto3D's own validation as incapable of handling given structures SHALL NOT be retried for those structures — not via alternative command candidates, not at smaller batch sizes.

#### Scenario: Auto3D validation rejects the engine
- **WHEN** Auto3D output states the engine cannot handle specific structures
- **THEN** subsequent attempts for those structures use a different, supporting engine or the declared fallback path, and the incompatible engine is invoked at most once for them
