# Design: replace-molscrub-with-unipka

## Context

The protonation stage currently loops per molecule in `engine.py` and calls `generate_protomer_candidates()` → `generate_molscrub_candidates()` (Python API, CLI fallback), with per-molecule recovery checkpoints and `_load_existing_protomer_outputs()` reuse. The Uni-Pka implementation to integrate is EasyDock's `containers/unipka/unipka.py` (error-corrected fork), executed via its container entrypoint `protonate`. Its CLI accepts tab-separated `SMILES<TAB>name` input and writes `form_smi<TAB>name<TAB>occupancy` rows (one per form, `NA` occupancy for unpredictable molecules, input SMILES echoed); `--distribution-file` writes `name, input_smi, microstate_smi, dG, occupancy, pH` rows for the whole `--ph-range`/`--ph-step` grid. The container is based on a pinned Uni-Mol/torch image; DSVR's env has no torch.

## Goals / Non-Goals

**Goals:**
- One container invocation per run for the whole protonation stage, with per-molecule recovery semantics preserved.
- Occupancy/dG parsed once and reused for selection, metadata, summary properties, and artifact — no second parse of model internals.
- Molscrub path behavior bit-identical when selected.

**Non-Goals:**
- Importing `unipka.py` as a Python module in DSVR's env (torch 1.11 pin conflicts).
- Caching container results across runs beyond existing stage-skip/checkpoint mechanics.
- Streaming/incremental input to Uni-Pka (batch size is one run).

## Decisions

**D1 — Container-only execution, runtime auto-detected.** `unipka_runner.py` builds `docker run -i --rm <image> protonate ...` or `apptainer exec <sif> unipka.py protonate ...` (runscript `protonate` for `.sif` via `apptainer run`). Config: `protonation.unipka.container` (image name or `.sif` path); optional `protonation.unipka.runtime: auto|docker|apptainer` (`shutil.which` autodetect: apptainer preferred when a `.sif` is configured). Alternative considered: local pip install — rejected (heavy pinned torch deps, no benefit; container recipe is the upstream-supported path and forward-compatible with HPC Apptainer).

**D2 — Batch pre-pass in the engine.** For `tool: unipka` and a non-skipped protonation step, the engine first collects molecules needing work (not checkpointed, no existing outputs), calls `generate_unipka_batch(molecules, config)` once, and consumes the per-molecule results inside the existing loop so recovery, fallback (`keep_fallback_parent_state`), progress, and `_load_existing_protomer_outputs` semantics are unchanged. When every molecule is checkpointed, no call happens. `protonation.enabled: false` continues to bypass the tool entirely. Alternative (lazy cache on first per-molecule call): hidden ordering coupling, rejected.

**D3 — Single source of truth is the distribution TSV.** Run Uni-Pka always with `--distribution-file` (range `ph_range = [protonation.unipka.ph_range_low, ph_range_high]`, default `2–12`, step default `0.25`, always containing `chemistry.ph`). The main output gives selection (forms + occupancies); the distribution TSV is parsed to attach per-microstate dG to each selected protomer (match on canonicalized microstate SMILES) and to compute the per-molecule summary. Main output kept in memory only (not a run artifact); distribution TSV stored at `enumeration/protomers/unipka_distribution.tsv`. Alternative (parse occupancies only, skip dG): loses dG and pKa-derived features; rejected — features 1–7 are an explicit goal.

**D4 — Summary properties computed in DSVR from the grid (not by Uni-Pka).** From per-microstate dG the script re-derives occupancies at any pH (same Boltzmann reweighting Uni-Pka uses: shift by max exponent). Derived in `chemistry/protonation_summary.py`:
- `occupancy_entropy` = Shannon entropy of microspecies occupancies at working pH (nats, pH-independent constant basis);
- `charge_population` = microspecies occupancies binned by formal charge (RDKit `GetFormalCharge`) at working pH;
- `microstate_count` = number of distinct microstates in the ensemble;
- `pI` = smallest pH in the grid where the occupancy-weighted net-charge crosses zero (linear interpolation between grid points; null if charge never changes sign — for always-positive molecules, null with warning);
- `pka_nearest_distance` = min over adjacent-charge transitions of |working pH − transition pH|, where a transition pH is where the occupancy-weighted net charge crosses an integer value;
- `top_two_occupancy_gap` = occupancy difference between the first and second selected forms (null/0 for single-form molecules).
All six stored in the protomer records' molecule-level metadata block written to the protonation stage outputs, and per-protomer `unipka_occupancy` / `unipka_dg` in each record's metadata (spec fields).

**D5 — Occupancy selection replaces plausibility scoring only on the Uni-Pka path.** `_select_plausible_protomers` stays for molscrub. Uni-Pka path: forms already arrive occupancy-filtered (`--occupancy`) and `-n`-capped by the tool; DSVR additionally applies `keep_best_per_charge` and `max_protomers_per_molecule` trimming by occupancy, `keep_input_state` prioritization (the input state is prioritized only when Uni-Pka predicted that form — unlike molscrub, it is never force-injected below the occupancy threshold; documented in `docs/limitations.md`), and the same exact-duplicate dedupe used cross-molscrub-candidates. Config: `protonation.unipka.min_occupancy` (default 0.05), `protonation.unipka.max_forms` (default = `max_protomers_per_molecule`, must be ≥ it — otherwise trimming by occupancy silently discards forms the cap would have kept).

**D6 — Config shape.** `ProtonationConfig.tool` default becomes `"unipka"`; new nested `unipka` model (`container`, `runtime`, `min_occupancy`, `max_forms`, `ph_range_low`, `ph_range_high`, `ph_step`, `timeout_seconds` default 3600 for the batch). `skip_gen3d_in_molscrub` renamed `molscrub_skip_gen3d` is deferred (out of scope churn): keep current key, document it as molscrub-only. Timeout semantics: the single batch call gets `timeout_seconds`; per-molecule timeout key remains molscrub-only.

**D7 — Failure classification.** New `UnipkaUnavailableError` (runtime or image missing — analogous to `MolscrubUnavailableError`, fails all items in stage) vs `UnipkaExecutionError` (non-zero exit/empty output — stage-level failure, no partial provenance). Molecule-level `NA` rows map to the existing input-state fallback warning path. Engine catches both new errors alongside `MolscrubUnavailableError` in the protonation stage's existing handler.

## Risks / Trade-offs

- [Name-based result matching breaks on duplicate input names] → DSVR input ids are already unique; runner writes its own SMILES+id TSV keyed by `input_id` as name and maps back, never user display names.
- [Batch call makes a single bad molecule fail the whole invocation] → upstream script is per-molecule fault-tolerant (NA rows, stderr warnings); add integration test with a poison molecule; if a torch-level crash still occurs, fall back per spec to input state with warning — verify during Vaclav run.
- [Grid resolution (0.25 pH) limits pI/pKa-transition precision] → documented interpolation to sub-grid precision; adequate for uncertainty features, not for publication pKa values (docs note).
- [Rankings shift vs molscrub-era runs] → provenance records tool + version (image digest recorded); changelog + limitations updated; molscrub selectable for A/B.
- [Container runtime absent on some workstations] → doctor check with acquisition hint (Zenodo `unipka.sif`, EasyDock Dockerfile); explicit failure, no silent molscrub fallback.
- [Long first-run install iteration] → acceptance protocol (tasks §6) iterates container build/pull until Vaclav run completes; environment fixes are apply-phase work, not code.

## Migration Plan

1. Implement + unit-test with recorded fixture outputs (no container needed in CI).
2. Locally build/pull the Uni-Pka image; run doctor; run smoke dataset; iterate install issues.
3. Full Vaclav run per acceptance gate; multi-agent review; PR.
4. Rollback: `protonation.tool: molscrub` restores prior behavior without code revert.

## Open Questions

- Exact Apptainer bind-mount flags needed for workdir visibility on this workstation vs clusters (resolve during §6 iteration; runner keeps a small, documented flag set).
- Whether `--distribution-min-occupancy` should be user-configurable now (leave at default 0.01; revisit if artifact size is problematic on large runs).
