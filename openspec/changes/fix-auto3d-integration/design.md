# Design: fix-auto3d-integration

## Context

All Auto3D calls go through `runners/auto3d_runner.py::run_auto3d` as subprocesses. For Auto3D v3, DSVR generates a wrapper script (`_auto3d_v3_wrapper.py`) that monkey-patches `multiprocessing.Manager` with a `_FakeManager` whose `Queue()` returns a fork-context `mp.Queue` (auto3d_runner.py:476-488); Auto3D/torch spawn workers then crash with the SemLock error. The shim exists because some sandboxes block `socket.setsockopt`, which a real `SyncManager` needs at startup. Engine candidates loop without capability knowledge; `_is_terminal_auto3d_selection_failure` only recognizes two selection-failure signatures, so "Only AIMNET can handle: ..." does not stop retries. `_should_use_gpu` only tests device-node existence; `_auto3d_env` sets `CUDA_VISIBLE_DEVICES=""` for CPU, but warp still initializes CUDA at import; the legacy path injects `--gpu_idx 0` whenever `cpu_workers>1` (auto3d_runner.py:244-245). All filtering/final-3D callers pass `internal_tautomer_stereo_enum=False`, while `stereochemistry.py` timeout/error fallbacks can retain mols with unspecified stereo. Failure text balloons via `_run_auto3d_for_tautomers` concatenating per-engine 4 KB output tails into one string copied into every `TautomerRecord`.

## Goals / Non-Goals

**Goals:**
- Make Auto3D stages execute for real in the supported environments; keep RDKit fallbacks as genuine degradation.
- Deterministic no-retry rules for proven-incompatible engines and globally-broken invocations.
- One explicit stereo policy across all Auto3D stages.

**Non-Goals:**
- Replacing subprocess invocation with an in-process Auto3D API.
- Auto-tuning engine choice for accuracy (only compatibility, not quality).
- Status/GUI presentation of failures (change `unify-run-status-and-gui-consistency`).

## Decisions

1. **Remove the `_FakeManager` patch; make the v3 wrapper spawn-consistent.** *(As implemented in PR #3: the fake manager was NOT removed — `_FakeManager.Queue()` now creates queues from the spawn context in `src/dsvr/runners/_auto3d_mp_shim.py`, which eliminates the fork/spawn SemLock conflict with far less risk than redesigning manager startup. The AF_UNIX `SyncManager` redesign and its env-flag fallback were not needed. Evidence: spawn-worker regression test passes, and runs after 2026-08-28 contain no SemLock errors with real Auto3D rankings.)* Replace the shim with `multiprocessing.get_context("spawn")` usage inside the wrapper for any queue Auto3D needs, and call `torch.multiprocessing`/Auto3D entry points under a single start method. Where the real `SyncManager` is needed but blocked by `setsockopt`, bind the manager to a Unix-domain address (`AF_UNIX`) instead of a TCP socket, which avoids `setsockopt` entirely — this removes the original motivation for the fake without reintroducing fork-context primitives. Acceptance test: the SemLock error signature must not appear in a run that exercises any engine.
   - Alternative considered: force `fork` everywhere (`multiprocessing.set_start_method("fork", force=True)` + torch). Rejected: torch/CUDA explicitly discourages fork with CUDA; spawn-consistent is the documented-safe route.
2. **Static engine-element capability table + preflight check.** Encode element support per engine (ANI2x/ANI2xt: H,C,N,O,F,S,Cl; AIMNET/AIMNet2: broader main-group set per Auto3D's own validation messages) as data in `auto3d_runner.py` derived from Auto3D's documented/observed validation, with a config override map. A `partition_by_engine(molecules, primary, fallback)` helper groups a batch into per-engine sub-batches before invocation; molecules with no supporting configured engine go straight to the recorded RDKit fallback with failure class `ENGINE_INCOMPATIBLE`.
   - Alternative: query Auto3D per batch — rejected; the validation only runs after expensive model load and only in-process.
3. **Extend terminal-failure detection.** Add the engine-incompatibility signature ("Only X can handle") and CUDA-absence signature to the terminal detector so the candidate loop stops instead of trying every command form; after an incompatibility hit, retry those structures only with a supporting engine. Smaller-batch retry (`final3d.py::_run_final_auto3d_smaller_batches`, engine.py seed retries) consumes the same terminal detection, so a proven-global failure can't recur per molecule.
4. **Runtime GPU probe.** Replace `_should_use_gpu`'s node-existence check with a probe subprocess (`nvidia-smi` plus a minimal `torch.cuda.is_available()`/warp import check inside the Auto3D env); cache the verdict per run. CPU mode: omit all GPU flags (including the legacy `--gpu_idx` injection — delete that dead-ish path), keep `CUDA_VISIBLE_DEVICES=""`, and set env vars that keep warp off the driver (`WARP_CACHE_PATH` kept; document the residual import-time noise and suppress it via `warp` config env if available). CPU/GPU verdict and mode recorded into the run's progress warnings once, not per molecule.
5. **Stereo policy `on_unspecified_stereo: enumerate | auto3d_enumerate`, default `enumerate`.** Before any Auto3D stage call, count unassigned stereo (`includeUnassigned=True`); with `enumerate`, run the existing stereo enumeration for those mols up-front (science-preferred: the pipeline already enumerates stereo deterministically); with `auto3d_enumerate`, pass `--enumerate-isomer` for the affected sub-batch only. Record the policy per variant in provenance. This also fixes the contradictory `--no-enumerate-isomer` + unspecified-stereo combination Auto3D warns about.
6. **Root-cause failure records.** Introduce a `root_cause_id` (hash of stage + error class + normalized diagnostic excerpt): the first occurrence writes one warnings record containing class, ≤1 KB excerpt, stage, engines; subsequent occurrences write only `{ref: root_cause_id}` on the affected candidate with a short note (`auto3d_failed:ENGINE_INCOMPATIBLE`). Per-run failure-memory: a per-stage set of terminal root causes; when a new unit matches, skip invocation and apply fallback immediately. Diagnostic stores (`progress.jsonl`, `warnings.jsonl`, record CSVs) all use the same reference so the GUI can link them.

## Risks / Trade-offs

- Capability table may drift from actual Auto3D behavior → config override + terminal-signature fallback ensures behavior stays correct even if the table is wrong; verified against Auto3D's own messages in tests with recorded outputs.
- Removing `_FakeManager` may break the sandbox that motivated it → CI/manual check in that sandbox; AF_UNIX manager kept as fallback path behind env flag.
- Up-front stereo enumeration changes variant counts/provenance versus previous runs → intended and recorded; release note in change summary.
- Failure-memory could mask a transient per-molecule error → memory keys on *terminal/infrastructure* signatures only; molecule-level chemistry errors still invoked individually.
