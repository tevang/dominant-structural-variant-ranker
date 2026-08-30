# auto3d-execution-environment

## Purpose

Invokes Auto3D in a multiprocessing-context-safe way with honest CPU/GPU mode handling and explicit stereochemistry treatment, so all engines can actually execute and CPU runs stay CUDA-free.

## ADDED Requirements

### Requirement: Compatible multiprocessing contexts
All multiprocessing primitives created by or injected into the Auto3D invocation SHALL be compatible with the start method of the worker processes that use them; a run MUST NOT fail with the fork/spawn SemLock sharing error for any engine.

#### Scenario: Tautomer ranking with each engine
- **WHEN** the Auto3D tautomer filter runs with ANI2xt, AIMNET, or AIMNet2
- **THEN** no SemLock fork/spawn RuntimeError occurs and results or a genuine chemistry/engine failure are returned

#### Scenario: Smaller-batch retry
- **WHEN** a batch fails at runtime and is retried in smaller batches
- **THEN** the retries run under the same fixed multiprocessing setup and do not reproduce the SemLock error

### Requirement: Genuine CPU execution
When GPU use is disabled or no usable GPU exists, the system SHALL run Auto3D without passing GPU flags and without initializing CUDA-dependent components; CPU execution MUST NOT emit CUDA initialization errors.

#### Scenario: CPU-configured run
- **WHEN** the run is configured with GPU disabled
- **THEN** Auto3D is invoked without GPU arguments and no "no CUDA-capable device" errors appear in the run output

### Requirement: Verified GPU usability
GPU use SHALL require a verified working CUDA setup, not merely the presence of device nodes; when verification fails, the system SHALL degrade to CPU mode and record that degradation.

#### Scenario: Device node present, driver broken
- **WHEN** `/dev/nvidia*` nodes exist but CUDA initialization fails during the availability probe
- **THEN** the run proceeds in CPU mode and a warning records that GPU was unavailable despite device nodes

### Requirement: Explicit handling of unspecified stereochemistry
Molecules carrying unspecified stereochemical elements SHALL NOT enter Auto3D with isomer enumeration disabled. The system SHALL follow one explicit, configurable policy: enumerate unspecified stereochemistry before the Auto3D call, or enable Auto3D isomer enumeration for those molecules. The policy and its application SHALL be recorded in variant provenance.

#### Scenario: Molecule with undefined chiral center reaches an Auto3D stage
- **WHEN** a molecule with unspecified stereochemistry is passed to an Auto3D stage (tautomer filter, stereo filter, or final 3D)
- **THEN** the configured policy is applied, Auto3D does not receive unspecified-stereo structures with enumeration disabled, and the variant's provenance records the treatment

#### Scenario: Default policy documented
- **WHEN** the user has not configured a policy
- **THEN** the documented default policy is applied consistently across all Auto3D stages
