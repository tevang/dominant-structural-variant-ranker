# Dominant Structural Variant Ranker

`dominant-structural-variant-ranker` (`dsvr`) is a Python orchestration package
for ranking pH- and solvent-dependent structural variants of small molecules
using maintained open-source tools.

This repository is a wrapper/orchestrator. It does **not** vendor, mirror, or
clone third-party repositories. Tools such as RDKit, molscrub, Auto3D, xTB,
CREST, CENSO, Psi4, and PySCF remain external dependencies installed through
conda-forge, pip, or their upstream installation instructions.

## Scope

The package coordinates:

1. Input parsing from SMILES and SDF.
2. Standardization and identifier generation.
3. pH/protomer/protonation candidate generation through optional external tools.
4. Tautomer and stereoisomer enumeration.
5. RDKit or Auto3D seed conformer generation.
6. Optional xTB/CREST conformer search and thermochemistry parsing.
7. Optional CENSO, Psi4, or PySCF refinement/rescoring.
8. Approximate population-oriented ranking and reporting.

## Scientific Limitation

Default pH handling controls candidate generation. It does not by itself produce
rigorous pH-dependent populations across protonation states. Cross-protomer or
cross-charge population estimates must be treated as approximate unless explicit
micro-pKa or proton chemical-potential corrections are supplied by a future
extension.

## Install

```bash
conda env create -f environment.yml
conda activate dominant-structural-variant-ranker
python -m pip install -e ".[dev]"
```

Or use the bootstrap script:

```bash
scripts/bootstrap_conda.sh
```

Optional Python tools:

```bash
scripts/bootstrap_mamba.sh --with-auto3d --with-molscrub
```

## CLI

```bash
python -m dsvr.cli --help
dsvr --help
dsvr doctor
dsvr run examples/test_molecules_minimal.smi --config configs/fast_smoke.yaml --outdir runs/smoke
```

## Development

```bash
pytest
ruff check src tests
mypy src
```

