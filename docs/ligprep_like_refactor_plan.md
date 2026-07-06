# LigPrep-like Plausible-Variant Workflow Refactor Plan

## Current Repository Structure Summary

- `src/dsvr/cli.py` defines the Typer CLI, including `dsvr run`, `doctor`, manual step commands, and global config loading.
- `src/dsvr/config.py` defines the strict Pydantic configuration model and current default values.
- `src/dsvr/workflow/engine.py` orchestrates the full workflow and contains both the current default path and the `auto3d_entropy` protocol branch.
- `src/dsvr/workflow/steps.py` defines planned workflow stages, dry-run plans, resume markers, and expected external tools.
- `src/dsvr/chemistry/standardize.py`, `src/dsvr/io/read_inputs.py`, `src/dsvr/io/smiles.py`, and `src/dsvr/io/sdf.py` handle input reading, validation, and standardization.
- `src/dsvr/chemistry/protonation.py` and `src/dsvr/runners/molscrub_runner.py` handle pH/protomer generation through molscrub with fallback behavior.
- `src/dsvr/chemistry/tautomers.py` implements RDKit tautomer enumeration with caps and timeouts, but it is candidate generation only.
- `src/dsvr/chemistry/stereochemistry.py` implements RDKit stereoisomer enumeration with caps and `tryEmbedding` heuristics.
- `src/dsvr/chemistry/conformers_auto3d.py` and `src/dsvr/runners/auto3d_runner.py` integrate Auto3D for representative conformer generation and approximate NNP energy triage.
- `src/dsvr/chemistry/conformers_rdkit.py` provides ETKDG conformer seeding for the current default path.
- `src/dsvr/filtering/selection.py`, `src/dsvr/filtering/stereo_reduce.py`, `src/dsvr/filtering/variant_score.py`, and `src/dsvr/filtering/xtb_prefilter.py` implement cheap filtering, enantiomer collapse, and optional xTB prefiltering.
- `src/dsvr/runners/crest_runner.py`, `src/dsvr/runners/xtb_runner.py`, `src/dsvr/runners/xtb_prefilter_runner.py`, `src/dsvr/runners/censo_runner.py`, `src/dsvr/runners/psi4_runner.py`, and `src/dsvr/runners/pyscf_runner.py` provide physics-heavy optional execution.
- `configs/default.yaml` is currently described as the default physics-heavy workflow.
- `configs/physics_heavy.yaml` already provides a named high-accuracy/expensive profile.
- `configs/auto3d_entropy_protocol.yaml`, `configs/auto3d_entropy_smoke.yaml`, `configs/auto3d_internal_enum.yaml`, and related Auto3D configs provide partial precedent for a cheaper Auto3D-centered protocol.

## Current Workflow Entry Points

- Main CLI entry point: `dsvr.cli:app` via the `dsvr` console script in `pyproject.toml`.
- Normal orchestration: `dsvr run` loads `RunConfig` through `load_config()` / `merge_cli_overrides()` and calls `run_workflow()` in `src/dsvr/workflow/engine.py`.
- Planned steps are surfaced by `planned_steps()` in `src/dsvr/workflow/steps.py`.
- The default `run_workflow()` path is:
  input reading -> standardization -> molscrub protomers -> RDKit tautomers -> RDKit stereoisomers -> cheap filtering -> RDKit/Auto3D 3D seeds -> seed filtering -> optional xTB prefilter -> CREST/xTB -> xTB thermo -> ranking -> optional CENSO -> optional Psi4/PySCF -> reports.
- The `protocol: auto3d_entropy` branch in `run_workflow()` already bypasses the RDKit tautomer/stereo stages, CREST/xTB thermo, CENSO, and QM rescoring in favor of Auto3D representative generation and approximate ranking.

## Current Config Files and Defaults

- `src/dsvr/config.py` sets `RunConfig.protocol = "default"`, `CrestConfig.enabled = True`, `ThermoConfig.enabled = True`, `SeedingConfig.method = "etkdg"`, `SeedingConfig.auto3d_k = 3`, and broad enumeration caps such as 32 protomers, 32 tautomers per protomer, and 64 stereoisomers per tautomer.
- `configs/default.yaml` describes itself as `Default configuration for the physics-heavy workflow (molscrub -> RDKit -> CREST/xTB)` and enables CREST and xTB thermo by default.
- `configs/physics_heavy.yaml` enables CENSO and uses larger tautomer/stereo caps; it should become the opt-in profile for expensive validation/refinement.
- `configs/auto3d_entropy_protocol.yaml` disables CREST, xTB prefilter, thermo, CENSO, and QM, uses `protocol: auto3d_entropy`, `seeding.method: auto3d`, `auto3d_k: 1`, and Auto3D internal tautomer/stereo enumeration.
- `configs/fast_smoke.yaml` and smoke tests should remain small CI-oriented checks.

## Risky Defaults Causing Combinatorial Explosion

- The current default multiplies `max_protomers_per_molecule` x `max_tautomers_per_protomer` x `max_stereoisomers_per_tautomer` before 3D generation. With defaults from `src/dsvr/config.py` and `configs/default.yaml`, one molecule can expose up to 32 x 32 x 64 candidate states before later filtering.
- `stereo_try_embedding: true` can be expensive because RDKit embedding is attempted during stereoisomer enumeration.
- `rdkit_num_conformers: 30` and `variant_filtering.max_seeds_per_variant: 2` still leave many conformers for downstream physics-heavy stages.
- `crest.enabled: true` and `thermo.enabled: true` in the default profile send selected variants into CREST/xTB conformer search and xTB Hessian/thermo.
- The optional Psi4/PySCF machinery is disabled by default, but the CLI and workflow still frame QM rescoring as part of the final protocol path rather than a clearly separate validation/refinement mode.
- `configs/physics_heavy.yaml` and the old default both make the expensive workflow easy to select accidentally.

## Files Needing Changes

- `src/dsvr/config.py`: add a first-class LigPrep-like protocol name/default, add explicit Auto3D triage and local-agent config sections, reduce default caps, and make CREST/xTB thermo/QM opt-in.
- `src/dsvr/workflow/engine.py`: make the default path use early protomer pruning, Auto3D tautomer energy triage, tautomer pruning before stereochemistry, one-conformer stereoisomer triage, and one final optimized 3D conformer per surviving variant.
- `src/dsvr/workflow/steps.py`: update planned default stages and external-tool requirements so CREST/xTB/QM are no longer shown as default requirements.
- `src/dsvr/chemistry/protonation.py`: expose bounded protomer filtering metadata and ensure molscrub/fallback states are capped before downstream enumeration.
- `src/dsvr/chemistry/tautomers.py`: keep RDKit tautomer enumeration available, but use it as a candidate generator feeding Auto3D energy triage rather than as an abundance predictor.
- `src/dsvr/chemistry/stereochemistry.py`: enumerate stereoisomers after tautomer pruning, preserve assigned stereo by default, and keep enantiomer-collapse behavior explicit for achiral environments.
- `src/dsvr/chemistry/conformers_auto3d.py`: add reusable helpers for tautomer ranking, stereoisomer one-conformer ranking, final one-conformer generation, strict caps, timeouts, and output provenance.
- `src/dsvr/runners/auto3d_runner.py`: ensure all Auto3D calls can be bounded by timeout, model, workers, GPU use, and failure/fallback policy.
- `src/dsvr/filtering/selection.py`, `src/dsvr/filtering/stereo_reduce.py`, and `src/dsvr/filtering/variant_score.py`: align filtering decisions with the new early-pruning protocol and preserve audit tables.
- `src/dsvr/runners/crest_runner.py`, `src/dsvr/runners/xtb_runner.py`, `src/dsvr/runners/psi4_runner.py`, and `src/dsvr/runners/pyscf_runner.py`: keep as optional validation/refinement tools only.
- `src/dsvr/utils/tool_check.py`: update required-vs-optional tool reporting for the new default.
- `src/dsvr/models.py`: add metadata fields if needed for Auto3D triage energies, pruning reasons, and local-agent diagnostics.
- `configs/default.yaml`: replace the physics-heavy default with the bounded LigPrep-like protocol.
- `configs/physics_heavy.yaml`: keep the old expensive workflow under an explicit opt-in name.
- `configs/auto3d_entropy_protocol.yaml`: either migrate into the new default or keep as a backward-compatible alias.
- `README.md`, `docs/workflow.md`, `docs/architecture.md`, `docs/limitations.md`, `docs/external_tools.md`, and `docs/installation.md`: document the new default, limits, optional physics-heavy validation, and interpretation of Auto3D energies.
- `tests/test_config.py`, `tests/test_workflow_smoke.py`, `tests/test_fast_smoke.py`, `tests/test_auto3d.py`, `tests/test_auto3d_runner.py`, `tests/test_tautomer_enumeration.py`, `tests/test_stereochemistry.py`, `tests/test_filtering.py`, `tests/test_doctor.py`, and CLI/config tests: update and extend coverage for the new default.

## Proposed New Default Workflow

1. Read molecules from SMILES or SDF using existing input readers and write invalid-input diagnostics.
2. Standardize and validate molecules with RDKit logic in `src/dsvr/chemistry/standardize.py`.
3. Generate plausible pH/protomer states at target pH, default 7.0, using molscrub when available and current fallback pH-normalization behavior when not.
4. Filter protomers immediately with bounded per-molecule caps, preserving the original state and best representative charge/protomer families.
5. Generate tautomer candidates with RDKit or Auto3D-supported enumeration under strict caps and timeouts.
6. Rank tautomer candidates with Auto3D NNP potential-energy triage using an ANI2xt/AIMNet2-style model configuration.
7. Keep only low-energy tautomer candidates before stereoisomer enumeration.
8. Enumerate stereoisomers only for surviving tautomers, preserving assigned stereochemistry by default, applying timeout/caps, and collapsing enantiomer-equivalent work in achiral environments.
9. Rank stereoisomers with a fast Auto3D/NNP one-conformer optimization and filter high-energy stereoisomers.
10. Generate exactly one final optimized 3D conformer per surviving structural variant by default.
11. Rank final variants using approximate Auto3D potential-energy deltas as a practical plausibility ordering.
12. Keep CREST/xTB, xTB thermo, CENSO, Psi4, and PySCF off by default and available only through explicit validation/refinement configs.
13. Add an experimental local-agent layer using `codex --oss -m qwen3.6:35b`, disabled by default, limited to bounded diagnostic/retry tasks such as summarizing failed external-tool logs, suggesting parameter reductions, or classifying known failure modes. This layer must not launch unbounded chemistry jobs or modify chemistry outputs without explicit workflow support.

Auto3D energy ranking must be documented as approximate potential-energy triage for pruning and prioritization, not rigorous abundance prediction, pKa prediction, solvent thermodynamics, or a replacement for validated ensemble free energies.

## Backward Compatibility Strategy

- The old physics-heavy workflow should become optional, not default.
- Keep `protocol: default` accepted, but migrate its meaning to the LigPrep-like bounded workflow in the next implementation prompts.
- Preserve the old physics-heavy behavior in `configs/physics_heavy.yaml` and consider an explicit protocol value such as `physics_heavy` if needed.
- Keep `protocol: auto3d_entropy` as an alias or compatibility mode until configs and docs are fully migrated.
- Preserve existing output filenames where practical: `all_protomers.sdf`, `all_tautomers.sdf`, `all_stereoisomers.sdf`, `all_3d_conformers.sdf`, `ranked_variants.*`, `manifest.json`, and `report.md`.
- Add provenance fields and warnings rather than silently changing the scientific interpretation of rankings.
- Keep CREST/xTB/QM runner APIs stable so existing users can opt into the previous expensive validation path.

## Test Plan

- Config tests: assert the new default disables CREST, xTB thermo, CENSO, Psi4, and PySCF; uses tight enumeration caps; and generates one final conformer per surviving variant.
- CLI dry-run tests: assert planned default steps do not list CREST/xTB/QM as required external tools.
- Workflow smoke tests: run the LigPrep-like default on `examples/test_molecules_minimal.smi` with mocked/fallback Auto3D behavior and verify outputs, manifest, report, and resolved config.
- Auto3D runner tests: verify timeout/cap propagation, model selection, fallback handling, and one-conformer output.
- Tautomer pruning tests: verify low-energy tautomer filtering happens before stereoisomer enumeration.
- Stereochemistry tests: verify caps, timeout behavior, assigned-stereo preservation, and enantiomer-collapse metadata in achiral environments.
- Optional physics-heavy tests: assert `configs/physics_heavy.yaml` still plans/runs CREST/xTB/CENSO/QM branches only when explicitly enabled.
- Regression tests: ensure Auto3D energy deltas are labeled approximate and not reported as rigorous abundance predictions.
