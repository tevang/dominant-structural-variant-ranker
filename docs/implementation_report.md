# Implementation Report

Generated during final local integration in this checkout.

## 1. What Was Implemented

- Python package scaffold for `dsvr`, a wrapper/orchestrator for dominant structural variant ranking.
- Typer CLI with commands for doctor checks, input validation, staged enumeration/seeding/ranking, full workflow execution, and summary reporting.
- YAML configuration system using Pydantic models with defaults for pH `7.0`, solvent `water`, temperature `298.15 K`, ETKDG seeding, and `same_formula` population scope.
- Input handling for `.smi`, `.smiles`, `.txt`, `.sdf`, and `.sd`, including two-column SMILES, multi-molecule SDF, invalid input reporting, and deduplication.
- Lineage/provenance models for input molecules, protomers, tautomers, stereoisomers, seed conformers, CREST conformers, thermochemistry records, and ranked variants.
- Candidate generation and enumeration modules for molscrub protomer candidates, RDKit tautomers, RDKit stereoisomers, RDKit ETKDG seeds, and optional Auto3D seeds.
- External runner wrappers and parsers for molscrub, Auto3D, CREST, xTB, CENSO, Psi4, and PySCF.
- Ranking logic for relative free energies and Boltzmann populations with explicit population scope and approximate cross-protonation warnings.
- Future `ProtonationCorrectionProvider` interface for microstate pH corrections.
- Final output writers for CSV, SDF, JSON, provenance JSONL, manifest, logs, and Markdown reports.
- Installation helper scripts and documentation for local, conda/mamba, and HPC-style installs.

## 2. What Was Mocked In Tests

- molscrub output/API behavior is mocked so CI does not require molscrub.
- Auto3D output SDF parsing and missing-tool behavior are mocked.
- CREST output parsing and missing-tool behavior are mocked; real CREST execution is marked `external`.
- xTB subprocess execution and sample energy/thermo logs are mocked for parser/ranking tests.
- CENSO output parsing and optional execution paths are mocked.
- Psi4 and PySCF rescoring output paths are mocked; imports remain optional.
- Subprocess log monitoring and command failure behavior are tested with synthetic subprocesses.

## 3. External Tools Detected

Detected by `PYTHONPATH=src python -m dsvr.cli doctor --json`:

- `rdkit` Python module: available/importable.
- Output directory writability: available.
- CPU count: available, `32`.
- Disk space: available, approximately `377 GiB` free at `runs/dsvr`.

## 4. External Tools Missing

Detected by doctor in this environment:

- Required Python runtime check: missing because active `python` is `3.10.12`; the package requires Python `>=3.11`.
- `molscrub` Python module: missing.
- `scrub.py` / `molscrub` CLI: missing.
- `xtb` executable: missing.
- `crest` executable: missing.
- Optional `Auto3D` Python module and CLI: missing.
- Optional `censo` executable: missing.
- Optional `psi4` Python module and executable: missing.
- Optional `pyscf` Python module: missing.

`PYTHONPATH=src python -m dsvr.cli doctor --strict` exits nonzero with an actionable missing-tool list instead of crashing.

## 5. Workflow Steps Fully Functional Locally

The following steps are functional in the current checkout without external CREST/xTB binaries:

- Input validation for all eight supplied SMILES in `examples/test_molecules.smi`.
- SMILES and SDF parsing, name/property preservation, and invalid input reporting.
- RDKit standardization.
- Fallback single-protomer smoke candidate generation when molscrub is unavailable.
- RDKit tautomer enumeration.
- RDKit stereoisomer enumeration.
- RDKit ETKDG seed generation.
- Smoke-mode ranking from generated seed/CREST-like records when external CREST/xTB is disabled or unavailable.
- Final CSV, JSON, SDF, manifest, logs, provenance JSONL, summary tables, and report generation.

Validated command:

```bash
PYTHONPATH=src python -m dsvr.cli run examples/test_molecules_minimal.smi \
  --config configs/fast_smoke.yaml \
  --ph 7.0 \
  --solvent water \
  --out runs/fast_smoke_test \
  --overwrite
```

The command completed and wrote `runs/fast_smoke_test`.

## 6. Steps Requiring Installed External Binaries

The following production workflow steps require external tools:

- molscrub pH/protomer candidate generation requires molscrub Python API or CLI.
- Auto3D seed generation requires Auto3D Python package or CLI.
- CREST conformer search and ensemble reduction requires `crest` and `xtb`.
- xTB optimization, Hessian, and thermochemistry require `xtb`.
- CENSO refinement requires `censo`.
- Psi4 final rescoring requires Psi4.
- PySCF final rescoring requires PySCF.

External CREST/xTB tests were invoked with:

```bash
DSVR_RUN_EXTERNAL=1 pytest -m external
```

Result in this environment: `1 skipped, 89 deselected`, because `crest` and `xtb` are not on `PATH`.

## 7. Scientific Limitations

- DSVR does not provide rigorous pH-dependent populations across protonation states unless a micro-pKa or proton chemical-potential correction provider is explicitly available.
- In the default implementation, pH affects molscrub candidate generation only.
- CREST/xTB-derived free energies provide approximate ranking over generated candidates in the configured solvent model.
- RDKit tautomer enumeration is candidate generation, not tautomer stability ranking.
- RDKit stereoisomer enumeration is explicit and controlled by configuration; `tryEmbedding` is a heuristic.
- Auto3D should not double-enumerate tautomers or stereoisomers unless `auto3d_internal_tautomer_stereo_enum=true`.
- Boltzmann populations are comparable by default only within the same formula/proton-count group.
- `all_approximate` population scope deliberately mixes generated records and marks populations approximate.

## 8. Exact Commands For Full Physics-Heavy Workflow

Create an environment with Python 3.11+ and core Python dependencies:

```bash
conda env create -f environment.yml
conda activate dsvr
python -m pip install -e ".[dev]"
```

Install external tools through conda-forge, official binaries, or HPC modules. For example:

```bash
scripts/bootstrap_mamba.sh --with-molscrub --with-xtb --with-crest
```

Verify the environment:

```bash
dsvr doctor --strict --json
```

Run the full physics-heavy workflow:

```bash
dsvr run examples/test_molecules.smi \
  --config configs/physics_heavy.yaml \
  --ph 7.0 \
  --solvent water \
  --out runs/physics_heavy_water_pH7 \
  --overwrite
```

Optional refinement commands after preliminary ranking:

```bash
dsvr rank runs/physics_heavy_water_pH7 --population-scope same_formula
dsvr summarize runs/physics_heavy_water_pH7
```

## Final Integration Command Results

- `python -m pip install -e .`: failed in this shell because active Python is `3.10.12`, while the package requires `>=3.11`.
- `dsvr doctor --json`: failed as a console command because editable install did not complete; source-tree equivalent passed.
- `PYTHONPATH=src python -m dsvr.cli doctor --json`: passed and wrote `doctor_report.json`.
- `PYTHONPATH=src python -m dsvr.cli validate-input examples/test_molecules.smi --out runs/validate_test`: passed; `8` valid molecules, `0` invalid records.
- `PYTHONPATH=src python -m dsvr.cli run examples/test_molecules_minimal.smi --config configs/fast_smoke.yaml --ph 7.0 --solvent water --out runs/fast_smoke_test --overwrite`: passed.
- `pytest`: passed after configuring checkout `pythonpath`; `88 passed, 2 skipped`.
- `ruff check .`: passed when `ruff` was available on `PATH`.
- `mypy src/dsvr`: passed when `mypy` was available on `PATH`.

## Output Files

The fast smoke run produced the expected top-level outputs under `runs/fast_smoke_test`, including:

- `manifest.json`
- `resolved_config.yaml`
- `logs/workflow.log`
- `invalid_inputs.csv`
- `inputs.csv`
- `protomers.csv`
- `tautomers.csv`
- `stereoisomers.csv`
- `seeds.csv`
- `crest_conformers.csv`
- `thermo.csv`
- `ranked_variants.csv`
- `ranked_variants.json`
- `ranked_variants.sdf`
- `report.md`
- `summary.md`
- provenance JSONL files

The input validation run produced:

- `runs/validate_test/validation_report.json`
