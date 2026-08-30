# Tasks: fix-auto3d-integration

## 1. Multiprocessing context

- [x] 1.1 DONE-DIFFERENTLY (shipped in PR #3): `_FakeManager` was retained but extracted to `runners/_auto3d_mp_shim.py`; its `Queue()` is built from the spawn context, so no fork-context primitives are shared with Auto3D/torch spawn workers. The AF_UNIX SyncManager redesign proved unnecessary. Verify: regression test `test_install_fake_manager_queue_survives_spawn_worker` passes; real runs after 2026-08-28 show zero SemLock occurrences and real Auto3D rankings (previously 100% RDKit fallback)
- [x] 1.2 DSVR-injected multiprocessing primitives are spawn-context (the only ones DSVR injects); verified: no "SemLock created in a fork context" text appears in post-fix run logs, guarded by the regression test above
- [x] 1.3 Smaller-batch retries (`final3d.py::_run_final_auto3d_smaller_batches`) and seed retries (`engine.py::_generate_auto3d_seeds_with_retries`) both route through `run_auto3d`, hence through the same wrapped invocation; no separate fix needed

## 2. Engine selection

- [x] 2.1 Add engine→element capability table with config override in `auto3d_runner.py`; implement `partition_by_engine` batch splitter; verify: unit tests cover an ANI2xt-configured batch containing AIMNET-only elements, and an unsupported-by-all molecule going to recorded RDKit fallback
- [x] 2.2 Extend `_is_terminal_auto3d_selection_failure` with engine-incompatibility and CUDA-absence signatures; stop candidate loop on terminal signatures and never retry that engine for those structures (including smaller batches); verify: replayed "Only AIMNET can handle" output triggers exactly one invocation of the incompatible engine
- [x] 2.3 Apply `partition_by_engine` in final-3D and tautomer/stereo Auto3D callers; verify: mixed-element fixture produces per-engine sub-batches in recorded commands

## 3. CPU/GPU handling

- [x] 3.1 Replace device-node check with cached runtime GPU probe (nvidia-smi + CUDA init in subprocess); on failure degrade to CPU with one recorded warning; verify: probe unit tests mock both working and broken-driver cases
- [x] 3.2 CPU mode: omit all GPU flags, remove legacy `--gpu_idx` injection, keep CUDA hidden; verify: recorded command lines for a CPU run contain no GPU arguments and no warp CUDA error appears in captured output (with warp behavior documented/suppressed)
- [x] 3.3 Pass `TautomerFilteringConfig.timeout_seconds_per_protomer` through to the tautomer Auto3D invocation; verify: config test shows the timeout on the recorded invocation

## 4. Stereochemistry policy

- [x] 4.1 Add `on_unspecified_stereo` config (default `enumerate`) and preflight unspecified-stereo detection before every Auto3D stage; policy `enumerate` enumerates via existing stereo enumeration, `auto3d_enumerate` enables `--enumerate-isomer` for affected sub-batches; verify: tests for both policies and for the provenance note on treated variants
- [x] 4.2 Record policy application in variant provenance; verify: run artifacts show the treatment for a fixture molecule with an undefined chiral center

## 5. Failure reporting

- [x] 5.1 Implement root-cause records (`root_cause_id`, ≤1 KB excerpt, stage, engines) plus short per-candidate references in tautomer/stereo/final-3D fallback paths; verify: a multi-molecule failing fixture produces exactly one full diagnostic record and short notes elsewhere
- [x] 5.2 Implement per-stage terminal failure memory (infra signatures only) that skips subsequent identical invocations; verify: stage-level test with all units failing shows one Auto3D invocation total
- [x] 5.3 Bound per-candidate warning text; verify: generated CSVs/JSONL contain no field longer than the bound for failure notes

## 6. Verification

- [x] 6.1 `python -m pytest tests -q` passes including new Auto3D tests (non-external; mark real-Auto3D tests `external`)
- [x] 6.2 End-to-end run on the examples input with Auto3D available: tautomer filtering completes via Auto3D or a cleanly recorded fallback, no SemLock/CUDA-on-CPU errors, warnings deduplicated
- [x] 6.3 `ruff check` and `mypy` pass on touched modules
