# Tasks: replace-molscrub-with-unipka

## 1. Config and dispatch

- [x] 1.1 Add `UnipkaConfig` (container, runtime auto/docker/apptainer, min_occupancy 0.05, max_forms, ph_range_low 2.0, ph_range_high 12.0, ph_step 0.25, timeout_seconds 3600) to `ProtonationConfig` in `src/dsvr/config.py`; validator requiring `max_forms >= max_protomers_per_molecule` and pH range containing `chemistry.ph`; change `protonation.tool` default to `unipka` with allowed values `unipka|molscrub`. Verify: config unit tests for defaults, validation errors, and molscrub still accepted.
- [x] 1.2 Update `configs/ligprep_like_*.yaml`, `configs/physics_validation_optional.yaml`, and `examples/example_config.yaml` with the `unipka` block; keep one example with `tool: molscrub`. Verify: `verify_configs.py` passes; `dsvr` loads each config.

## 2. Uni-Pka runner

- [x] 2.1 Create `src/dsvr/runners/unipka_runner.py`: build docker/apptainer command per design D1, write `SMILES<TAB>input_id` TSV, invoke once per batch with `-n`, `--occupancy`, `--pH`, `--distribution-file`, `--ph-range`, `--ph-step` via `run_command`; parse main output (forms+occupancy, `NA` rows) and distribution TSV into per-`input_id` results; raise `UnipkaUnavailableError`/`UnipkaExecutionError` per D7. Verify: unit tests with fake executable and recorded fixture files covering multi-form, NA, empty-output, and non-zero-exit cases.
- [x] 2.2 Record container provenance: image/sif path, runtime, and full command into runner results. Verify: provenance assertion in runner unit test.

## 3. Selection, metadata, and summary

- [x] 3.1 Implement `src/dsvr/chemistry/protonation_summary.py`: occupancy reweighting from dG grid, occupancy entropy, charge populations, microstate count, top-two gap, pI with zero-crossing interpolation, nearest pKa-transition distance. Verify: unit tests with hand-computable ensembles (single form, 50/50 tie, triprotic pattern, never-changing charge sign → null pI + warning).
- [x] 3.2 Rework `src/dsvr/chemistry/protonation.py`: dispatch on `config.protonation.tool`; Uni-Pka path does dedupe, `keep_best_per_charge` via `GetFormalCharge`, occupancy-ordered trim to `max_protomers_per_molecule`, `keep_input_state`; attach `unipka_occupancy` + `unipka_dg` per protomer and the six summary fields per molecule; molscrub path untouched. Verify: unit tests for trim-by-occupancy, charge rule, input-state retention, NA fallback, and that molscrub output is unchanged (golden test).
- [x] 3.3 Batch pre-pass in `src/dsvr/workflow/engine.py` per design D2 (collect uncheckpointed molecules, single runner call, per-molecule consumption; skip call when nothing to do); catch new Unipka errors in the existing protonation-stage handler; update `recovery.py` if stage-skip semantics need it. Verify: workflow smoke test (fake runner) asserts exactly one batch call, resume skips call when all checkpointed, and failure fallback paths fire.

## 4. Artifacts, doctor, CLI, reporting

- [x] 4.1 Write `enumeration/protomers/unipka_distribution.tsv` (copy of runner distribution output) and surface it in stage outputs listing. Verify: smoke test asserts file presence and coverage of working pH.
- [x] 4.2 Replace the molscrub group in `src/dsvr/utils/tool_check.py`: required Uni-Pka group (runtime executable + configured image/sif exists + `protonate` entrypoint responds to help), optional molscrub group. Verify: `tests/test_doctor.py` updated and passing for both tool selections.
- [x] 4.3 Update `src/dsvr/workflow/steps.py` (step description + required tools follow `protonation.tool`) and `src/dsvr/cli.py` help strings that name molscrub. Verify: step-listing test for both tools.
- [x] 4.4 Extend `src/dsvr/reporting/markdown.py`: per-molecule protomer occupancies, summary properties (pI, nearest pKa distance, entropy, charge populations), distribution-file link; stage summary counts reflect Uni-Pka warnings. Verify: report snapshot/regex test from fixture run.

## 5. Provenance, tests, docs

- [x] 5.1 Update provenance model/tests (`tests/test_provenance_models.py`) for Uni-Pka `source_software` values and metadata fields. Verify: tests pass.
- [x] 5.2 Update `tests/test_protonation.py` and `tests/test_workflow_smoke.py` for the new default tool (fixtures, no real container in CI). Verify: full `pytest`, `ruff check src tests`, `mypy src` pass.
- [x] 5.3 Documentation: `docs/external_tools.md` (Uni-Pka entry: container acquisition — Zenodo unipka.sif / EasyDock build — molscrub demoted to optional alternative), `docs/installation.md` (container runtime install + image acquisition, bootstrap script flags), `docs/workflow.md` (occupancy-based protomer stage), `docs/architecture.md` (runner dispatch), `docs/limitations.md` (grid-resolution pKa/pI precision, occupancy semantics, ranking shift vs molscrub), `docs/changelog.md` entry; check `README.md` and `ONBOARDING.md` for molscrub-as-default statements and update; update `scripts/bootstrap_*.sh` flags (`--with-unipka` doc/pull hint, molscrub optional). Verify: `grep -rn molscrub docs/ README.md ONBOARDING.md` shows no stale default claims; use the `docs` skill for consistency.

## 6. End-to-end verification and review protocol

- [x] 6.1 Acquire Uni-Pka container locally (build easydock Dockerfile or convert Zenodo SIF); iteratively resolve install/runtime issues until `dsvr doctor` reports Uni-Pka available. Verify: doctor green. (Resolved: conda-forge Apptainer + Zenodo record 19627026 SIF, md5-verified; all Zenodo SIF builds ship a pre-occupancy `unipka.py`, fixed by vendoring the current EasyDock script to `containers/unipka.py` and bind-mounting it via new `protonation.unipka.script_path`; `container: unipka` now auto-maps to `containers/unipka.sif`.)
- [x] 6.2 Run the smoke dataset end-to-end with Uni-Pka; fix integration bugs (parsing, bind mounts, name mapping) until clean completion with distribution artifact and metadata fields 1–7 present. Verify: inspect `final_variants.sdf` + protomer stage outputs. (`runs/smoke_unipka_192941`: 8/8 molecules in distribution TSV, provenance metadata carries occupancy/dG + all summary fields, report shows Protonation Properties section, 58 final variants with 0 duplicate groups.)
- [x] 6.3 Run `dsvr prepare-ligands examples/Vaclavs_heterocycles_set1.smi --max-protomers 8 --tauto-k 10 --max-stereoisomers 24 --out runs/Vaclav_set1_p8_t10_s24_unipka` to completion, iterating on failures (background task; log triage may be delegated to a subagent). Verify: run completes; `duplicate_structures_report.txt` shows 0 duplicate groups; every Uni-Pka property present in outputs; no stage-summary errors. (`runs/Vaclav_set1_p8_t10_s24_unipka`: 11/11 molecules, 110 final variants, 0 duplicate groups, distribution TSV + occupancy/dG/summary fields in provenance, report section complete, failures.jsonl empty. Fixed: `--max-protomers` CLI override now re-defaults a cap-equal `unipka.max_forms` instead of hard-failing validation.)
- [ ] 6.4 Multi-agent pre-PR review: Workflow fan-out review of the diff (correctness; molscrub-path regression; config/provenance schema; docs completeness) plus built-in `/code-review`, adversarial verification of findings, fix all confirmed findings; re-run tests. Verify: review reports clean; pytest/ruff/mypy green.
- [ ] 6.5 Open PR, trigger the agentic PR reviewer, iterate on its findings until clean. Verify: reviewer reports no outstanding issues; CI green.
