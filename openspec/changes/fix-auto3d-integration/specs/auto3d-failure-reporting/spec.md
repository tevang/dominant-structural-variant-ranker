# auto3d-failure-reporting

## Purpose

Records Auto3D execution failures once with their root error, links affected molecules/variants to that root record, and keeps per-candidate notes short so warning stores stay readable.

## ADDED Requirements

### Requirement: Single root-cause record
An Auto3D invocation failure affecting multiple candidates SHALL be stored once as a root-cause record capturing the error class, a bounded excerpt of the diagnostic output, the stage, and the engine(s) attempted; individual candidates affected by it SHALL reference the root record instead of embedding its text.

#### Scenario: Global Auto3D outage during tautomer filtering
- **WHEN** Auto3D fails for infrastructure reasons across many protomers in one stage
- **THEN** the run contains one root-cause record for that stage's failure and each affected variant carries only a short reference/status, not the full traceback

### Requirement: Bounded per-candidate notes
Per-candidate warning text SHALL be bounded in length; raw tracebacks and tool output MUST NOT be duplicated into every rejected candidate's record.

#### Scenario: Inspecting a rejected tautomer
- **WHEN** a tautomer is rejected because an Auto3D failure forced the RDKit fallback
- **THEN** its record contains a short note (failure class + reference to the root record) traceable to the full diagnostic

### Requirement: Failure memory within a run
Within a run, the system SHALL remember per-stage Auto3D invocation failures that have a global cause and SHALL NOT re-attempt the same invocation pattern for every remaining molecule/protomer.

#### Scenario: Auto3D binary broken for the whole stage
- **WHEN** the first Auto3D invocation of a stage fails with an environment-level error (e.g. multiprocessing context crash)
- **THEN** remaining units in that stage use the declared fallback path directly, noting the earlier root failure, instead of each paying a full failing invocation
