# LigPrep-like Refactor Report

## 1. Summary of implemented changes

DSVR now has a bounded LigPrep-like ligand-preparation workflow intended for docking, ligand-based modeling, and batch library preparation. The workflow prioritizes plausible structural variants and early pruning instead of expanding every protomer, tautomer, stereoisomer, and conformer into a physics-heavy CREST/xTB workflow.

Implemented capabilities include:

- First-class LigPrep-like configs: `configs/ligprep_like_default.yaml`, `configs/ligprep_like_conservative.yaml`, and `configs/ligprep_like_aggressive.yaml`.
- `dsvr prepare-ligands` CLI entry point with overrides for pH, solvent, protomer cap, tautomer top-K/window, stereoisomer cap, optional CREST/xTB validation, and the local diagnostic agent.
- Bounded protomer generation through molscrub with per-molecule caps and rejection CSVs.
- RDKit tautomer enumeration followed by Auto3D energy triage before stereoisomer expansion.
- RDKit stereoisomer enumeration after tautomer filtering, with caps, timeouts, embedding checks, enantiomer collapse options, and rejection outputs.
- Final 3D generation that writes one selected conformer per surviving structural variant.
- Auto3D fallback behavior for GPU-to-CPU retry, optimizing-engine retry, smaller-batch retry, and RDKit fallback when Auto3D is unavailable or fails.
- Optional CREST/xTB, thermochemistry, CENSO, and QM refinement paths remain available but are not default.
- Experimental local qwen diagnostic-agent wiring is present but disabled by default and constrained to diagnostic/retry-menu tasks.

## 2. New default workflow

Recommended command:

```bash
dsvr prepare-ligands examples/test_molecules_minimal.smi \
  --config configs/ligprep_like_default.yaml \
  --ph 7.0 \
  --solvent water \
  --out runs/smoke
```

Default workflow:

```text
Input SMILES/SDF
-> standardization and validity checks
-> plausible pH/protomer generation at target pH
-> early protomer filtering
-> RDKit tautomer enumeration
-> Auto3D tautomer ranking/filtering
-> RDKit stereoisomer enumeration for selected tautomers
-> stereoisomer filtering
-> final Auto3D one-conformer 3D generation
-> final SDF/CSV/JSON outputs and summary reports
-> optional CREST/xTB validation only when explicitly enabled
```

The new default config exists at `configs/ligprep_like_default.yaml`, and `README.md` now recommends `dsvr prepare-ligands` with this config.

## 3. What remains optional

The following remain opt-in validation or refinement features:

- CREST conformer searches.
- xTB thermochemistry.
- CREST entropy estimates.
- CENSO refinement.
- Psi4/PySCF rescoring.
- Local qwen diagnostic agent.
- Physics-heavy configs such as `configs/physics_validation_optional.yaml`, `configs/physics_heavy.yaml`, and `configs/exhaustive_debug.yaml`.

These are useful for selected small candidate sets, not default large-library ligand preparation.

## 4. How protomers are filtered

The default LigPrep-like config uses molscrub in plausible mode:

- `protonation.enabled: true`
- `protonation.tool: molscrub`
- `protonation.mode: plausible`
- `protonation.max_protomers_per_molecule: 4`
- `protonation.keep_input_state: true`
- `protonation.keep_best_per_charge: true`
- `protonation.skip_gen3d_in_molscrub: true`

The conservative config lowers the protomer cap to 3. Rejected protomers are written to `protomers_rejected.csv` with rejection reasons so users can audit pruning decisions.

## 5. How tautomers are ranked and filtered with Auto3D

RDKit generates tautomer candidates, then Auto3D ranks them by approximate optimized-conformer energy before stereoisomer expansion. This avoids multiplying every tautomer by every possible stereoisomer before any energy signal is available.

Default tautomer controls:

- `tautomer_filtering.enabled: true`
- `tautomer_filtering.tool: auto3d`
- `tautomer_filtering.tauto_engine: rdkit`
- `tautomer_filtering.optimizing_engine: ANI2xt`
- `tautomer_filtering.fallback_optimizing_engine: AIMNET`
- `tautomer_filtering.max_rdkit_tautomers_before_auto3d: 64`
- `tautomer_filtering.rdkit_tautomer_timeout_seconds: 30`
- `tautomer_filtering.tauto_k: 3`
- `tautomer_filtering.tauto_window_kcal_mol: 5.0`
- `tautomer_filtering.keep_input_tautomer: true`

If RDKit tautomer enumeration times out, the fallback keeps the input tautomer and canonical tautomer rather than expanding unbounded candidates. Rejected and ranked tautomer decisions are written to files such as `tautomers_auto3d_ranked.csv`, `tautomers_selected.csv`, and `tautomers_rejected.csv`.

Auto3D energies are screening signals. They are not pKa predictions, solvated free energies, or rigorous tautomer abundance estimates.

## 6. How stereoisomers are enumerated and filtered

Stereoisomers are enumerated only after protomer and tautomer pruning.

Default stereoisomer controls:

- `stereoisomer_filtering.enabled: true`
- `stereoisomer_filtering.enumerator: rdkit`
- `stereoisomer_filtering.max_stereoisomers_per_tautomer: 16`
- `stereoisomer_filtering.timeout_seconds_per_tautomer: 300`
- `stereoisomer_filtering.try_embedding: true`
- `stereoisomer_filtering.only_unassigned: true`
- `stereoisomer_filtering.collapse_enantiomers_in_achiral_solvent: true`
- `stereoisomer_filtering.run_energy_for_enantiomer_representatives_only: true`
- `stereoisomer_filtering.stereo_energy_window_kcal_mol: 7.0`
- `stereoisomer_filtering.keep_top_n_diastereomers: 8`

The conservative config lowers `max_stereoisomers_per_tautomer` to 8 and `keep_top_n_diastereomers` to 4. Rejected stereoisomers are written to `stereoisomers_rejected.csv`.

## 7. How final one-conformer output is generated

Final 3D generation uses Auto3D by default:

- `final_3d.tool: auto3d`
- `final_3d.optimizing_engine: AIMNET`
- `final_3d.fallback_optimizing_engine: ANI2xt`
- `final_3d.use_gpu: true`
- `final_3d.k: 1`
- `final_3d.max_confs: 10`
- `final_3d.one_conformer_per_variant: true`

The final outputs include:

- `final_variants.sdf`
- `final_variants.csv`
- `final_variants.json`
- `ranked_variants.sdf`
- `ranked_variants.csv`
- `ranked_variants.json`

When Auto3D is unavailable or fails, final generation retries GPU then CPU, retries the fallback optimizing engine, retries smaller batches, and then falls back to an RDKit one-conformer output with explicit warnings. The conservative smoke run on this machine used the RDKit fallback because Auto3D is not installed, but still completed and produced five final variants.

## 8. Why CREST/xTB is no longer default

CREST/xTB is more expensive than needed for routine ligand-preparation batches. Running CREST/xTB before early pruning increases runtime, disk use, and failure surface area, especially when protomer, tautomer, and stereoisomer counts multiply.

The default workflow now uses Auto3D as fast triage and final one-conformer generation. CREST/xTB remains appropriate for targeted validation of small selected candidate sets where the extra cost is justified.

## 9. How to enable optional CREST/xTB validation

Use the CLI flag:

```bash
dsvr prepare-ligands examples/test_molecules_minimal.smi \
  --config configs/ligprep_like_default.yaml \
  --enable-crest-validation \
  --out runs/ligprep_with_crest_validation
```

Or set this in YAML:

```yaml
optional_validation:
  crest_xtb_enabled: true
  xtb_thermo_enabled: true
  crest_entropy_enabled: true
```

Run `dsvr doctor` first to verify that optional external tools are installed and visible on PATH.

## 10. Local qwen agent limitations and safety policy

The local qwen agent is disabled by default:

```yaml
agent:
  enabled: false
  backend: ollama_codex_cli
  command: "codex --oss -m qwen3.6:35b"
```

Its intended scope is constrained diagnostics:

- Classify external-tool failures.
- Summarize compact bug packages.
- Suggest bounded retry actions from a deterministic menu.

It must not silently patch code, change scientific thresholds, delete outputs, or launch large reruns. The config requires explicit user approval for code patches, science-threshold changes, output deletion, and large reruns. The agent receives bounded context only and should not receive secrets, credentials, private keys, tokens, browser cookies, `.env` contents, or credential stores.

## 11. Test results

Commands run for this final prompt:

```bash
python -m pip install -e .
```

Result: passed.

```bash
dsvr --help
```

Result: passed after rerunning outside the sandbox wrapper.

```bash
dsvr prepare-ligands --help
```

Result: passed after rerunning outside the sandbox wrapper.

```bash
pytest
```

Result: passed: `160 passed, 2 skipped, 43 warnings in 12.95s`.

```bash
ruff check .
```

Result: not run because `ruff` is not installed on PATH: `/bin/bash: line 1: ruff: command not found`.

```bash
mypy src
```

Result: not run because `mypy` is not installed on PATH: `/bin/bash: line 1: mypy: command not found`.

Integration command:

```bash
dsvr prepare-ligands examples/test_molecules_minimal.smi \
  --config configs/ligprep_like_conservative.yaml \
  --ph 7.0 \
  --solvent water \
  --out runs/ligprep_conservative_smoke \
  --overwrite
```

Result: passed after rerunning outside the sandbox wrapper. The run completed all nine workflow stages and wrote `runs/ligprep_conservative_smoke`. Stage summary reported two input molecules and five final variants. Auto3D was unavailable in this environment, so final 3D generation used the explicit RDKit fallback and wrote warnings in `final_variants.csv`.

## 12. Exact commands for users

Install and check the project:

```bash
python -m pip install -e .
dsvr doctor
dsvr --help
dsvr prepare-ligands --help
```

Run the default LigPrep-like workflow:

```bash
dsvr prepare-ligands examples/test_molecules_minimal.smi \
  --config configs/ligprep_like_default.yaml \
  --ph 7.0 \
  --solvent water \
  --out runs/ligprep_like_default_smoke \
  --overwrite
```

Run the conservative workflow:

```bash
dsvr prepare-ligands examples/test_molecules_minimal.smi \
  --config configs/ligprep_like_conservative.yaml \
  --ph 7.0 \
  --solvent water \
  --out runs/ligprep_conservative_smoke \
  --overwrite
```

Override key caps from the CLI:

```bash
dsvr prepare-ligands examples/test_molecules_minimal.smi \
  --config configs/ligprep_like_default.yaml \
  --max-protomers 4 \
  --tauto-k 3 \
  --tauto-window 5.0 \
  --max-stereoisomers 16 \
  --out runs/ligprep_custom_caps \
  --overwrite
```

Enable optional CREST/xTB validation:

```bash
dsvr prepare-ligands examples/test_molecules_minimal.smi \
  --config configs/ligprep_like_default.yaml \
  --enable-crest-validation \
  --out runs/ligprep_crest_validation \
  --overwrite
```

Run developer checks:

```bash
pytest
ruff check .
mypy src
```

Install missing developer tools before running `ruff` or `mypy` in environments where they are not already available.
