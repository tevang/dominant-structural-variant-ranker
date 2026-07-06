# Dominant Structural Variant Ranker

`dominant-structural-variant-ranker` (`dsvr`) is a Python orchestration package for preparing and ranking pH- and solvent-dependent small-molecule structural variants with maintained open-source tools.

This repository is a wrapper/orchestrator. It does **not** vendor, mirror, or clone third-party repositories. RDKit, molscrub, Auto3D, xTB, CREST, CENSO, Psi4, and PySCF remain external tools installed through conda-forge, pip, official binary distributions, or user-managed software modules.

## Default LigPrep-like Workflow

The recommended default is a bounded plausible-variant ligand-preparation workflow for docking, ligand-based modeling, and batch-library preparation:

```text
Input SMILES/SDF -> standardization and validity checks -> plausible pH/protomer generation at target pH, default pH 7.0 -> early protomer filtering -> Auto3D tautomer enumeration/ranking/filtering using RDKit tautomer engine and ANI2xt/AIMNet2 -> RDKit stereoisomer enumeration with timeout/caps after tautomer filtering -> Auto3D one-conformer optimization/ranking/filtering of stereoisomers -> final SDF/CSV/JSON report with one optimized 3D conformer per surviving structural variant -> optional CREST/xTB validation only if explicitly enabled
```

Start with:

```bash
dsvr run examples/test_molecules.smi   --config configs/ligprep_like_default.yaml   --outdir runs/ligprep_like_water_pH7
```

The old CREST/xTB-centered workflow is expensive and optional. Use `configs/physics_validation_optional.yaml` or `configs/physics_heavy.yaml` only for selected validation/refinement runs after the candidate set is small. `configs/exhaustive_debug.yaml` remains useful for small-molecule debugging, but it is intentionally expensive.

## Auto3D Energy Triage

RDKit alone can enumerate too many tautomers and does not rank tautomer abundance. The default workflow filters tautomers before stereoisomer enumeration because expanding stereoisomers for every tautomer multiplies the candidate count before any energy signal is available.

Auto3D ranking is approximate potential-energy triage. It ranks low-energy tautomer and stereoisomer candidates by optimized conformer energies, not by true solution abundance. Auto3D thermodynamics, when used, are not substitutes for validated solvated free energies.

## Scientific Warning

The default pipeline is fast ligand preparation, not an exhaustive conformational free-energy workflow. It does not perform rigorous pH-dependent population calculations, pKa prediction, or solution speciation.

CREST/xTB, xTB thermo, CREST entropy estimates, CENSO, and Psi4/PySCF rescoring are optional validation/refinement steps. Psi4/PySCF rescoring is outside the default workflow and should be treated as an advanced legacy module unless explicitly enabled.

## Quick Start

```bash
conda env create -f environment.yml
conda activate dsvr
python -m pip install -e ".[dev]"
dsvr doctor
dsvr run examples/test_molecules_minimal.smi --config configs/ligprep_like_default.yaml --outdir runs/smoke
```

For direct source-tree smoke checks:

```bash
PYTHONPATH=src python -m dsvr.cli --help
PYTHONPATH=src python -m pytest
```

## Dependency Strategy

- Do not vendor third-party repositories.
- Install Python packages via conda or pip.
- Install external binaries via conda, official binaries, or user-managed modules.
- Use `dsvr doctor` to verify the environment before running optional physics-heavy workflows.

Optional Python tools:

```bash
scripts/bootstrap_mamba.sh --with-auto3d --with-molscrub
```

## CLI

```bash
python -m dsvr.cli --help
dsvr --help
dsvr doctor
dsvr inspect examples/test_molecules.smi
dsvr run examples/test_molecules_minimal.smi --config configs/ligprep_like_default.yaml --outdir runs/smoke
```

## Documentation

- [Architecture](docs/architecture.md)
- [Workflow](docs/workflow.md)
- [Plausible variant workflow](docs/plausible_variant_workflow.md)
- [Limitations](docs/limitations.md)
- [External tools](docs/external_tools.md)
- [File formats](docs/file_formats.md)
- [Installation](docs/installation.md)

## Development

```bash
pytest
ruff check src tests
mypy src
```
