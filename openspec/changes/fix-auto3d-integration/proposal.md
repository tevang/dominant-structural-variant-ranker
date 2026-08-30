# Proposal: fix-auto3d-integration

## Why

Every Auto3D stage currently fails and silently falls back to RDKit. Root causes: the Auto3D v3 wrapper injects a fake `multiprocessing.Manager` whose queue is created in the fork context and then shared with Auto3D/torch spawn-context workers (deterministic `SemLock ... fork ... spawn` RuntimeError for ANI2xt, AIMNET, and AIMNet2 alike); engine selection ignores element compatibility, so ANI2xt is offered molecules Auto3D itself says "Only AIMNET can handle", and the incompatible engine is retried through every command candidate and again at smaller batch sizes; GPU "detection" only checks for `/dev/nvidia*` device nodes, so CUDA components (warp) initialize even on CPU runs; and molecules reaching Auto3D with unspecified stereochemistry are passed with isomer enumeration disabled, which Auto3D explicitly warns against. Tautomer filtering additionally stores the full concatenated traceback of all engine attempts on every rejected candidate of every molecule, duplicating kilobytes of identical warning text, and per-molecule smaller-batch retries deterministically reproduce the same SemLock failure.

## What Changes

- Fix the multiprocessing context conflict in Auto3D invocation so that all multiprocessing primitives and worker processes use compatible contexts; Auto3D stages must be able to complete without the fork/spawn SemLock error.
- Engine capability awareness: molecule–engine compatibility (allowed elements) is checked before invocation; molecules are routed to a supporting engine or the batch is split by required engine; an engine proven incompatible for a molecule is never retried for it (not in later candidates, not in smaller batches).
- Failure memory: within a run, an engine invocation pattern that fails for infrastructure reasons is not blindly re-attempted per molecule/protomer.
- Genuine CPU mode: with GPU disabled or unavailable, Auto3D runs without initializing CUDA-dependent components and without passing GPU flags; GPU usability is verified (driver/runtime), not inferred from device-node presence.
- Explicit stereochemistry policy: molecules with unspecified stereochemistry are either enumerated before Auto3D receives them or sent with isomer enumeration enabled; the chosen treatment is explicit, consistent across tautomer/stereo/final-3D stages, and recorded in provenance.
- Root-cause warning consolidation: an Auto3D execution failure is recorded once with its root error; rejected/affected candidates reference it, and each candidate carries only a short status note instead of duplicated tracebacks.
- Fallbacks must still guarantee that ranking continues (RDKit fallback retained), while clearly recording that a fallback occurred (accounting handled in change `unify-run-status-and-gui-consistency`).

## Capabilities

### New Capabilities
- `auto3d-engine-selection`: compatibility-aware selection of Auto3D engines, batch splitting by required engine, and no-retry guarantees for proven-incompatible engines.
- `auto3d-execution-environment`: multiprocessing-context-safe invocation of Auto3D with genuine CPU/GPU mode handling and verified GPU usability.
- `auto3d-failure-reporting`: consolidated, deduplicated recording of Auto3D execution failures with affected-candidate linkage.

### Modified Capabilities
(none — the existing stereochemistry behavior has no archived spec; the stereo policy lands in `auto3d-execution-environment` requirements.)

## Impact

- Code: `src/dsvr/runners/auto3d_runner.py` (generated wrappers, GPU/CPU env, candidate loop), `src/dsvr/chemistry/tautomer_auto3d_filter.py`, `stereo_auto3d_filter.py`, `final3d.py`, `conformers_auto3d.py`, `src/dsvr/chemistry/stereochemistry.py`, seed-retry logic in `workflow/engine.py`, config schema for any new policy knobs (`src/dsvr/config.py`).
- Behavior: Auto3D stages expected to actually run instead of falling back; fallback paths retained. Results may change versus previous runs (real Auto3D rankings instead of RDKit fallback rankings) — scientifically intended.
- Warnings artifacts get smaller and deduplicated; per-record warning text shape changes (consumed by GUI/errors panel).
- Removed: nothing user-facing; internal fake-manager shim replaced.
