# Dominant Structural Variant Ranker

`dominant-structural-variant-ranker` (`dsvr`) is a Python orchestration package
for ranking pH- and solvent-dependent structural variants of small molecules
using maintained open-source tools.

This repository is a wrapper/orchestrator. It does **not** vendor, mirror, or
clone third-party repositories. RDKit, molscrub, Auto3D, xTB, CREST, CENSO,
Psi4, and PySCF remain external tools installed through conda-forge, pip,
official binary distributions, or user-managed software modules.

## Default Physics-Heavy Workflow

The package implements this default workflow:

```text
molscrub protonation/protomer generation at target pH
-> RDKit tautomer enumeration
-> RDKit stereoisomer enumeration
-> RDKit ETKDG or Auto3D conformer seeding
-> CREST/xTB conformer search and ensemble reduction
-> xTB thermo / CREST entropy Delta G ranking
-> optional CENSO
-> optional Psi4/PySCF final rescoring
```

## Auto3D Entropy Protocol

The repository also includes an Auto3D-owned enumeration and ranking protocol:

```text
molscrub protonation/protomer generation at target pH
-> Auto3D tautomer enumeration
-> Auto3D stereoisomer enumeration
-> Auto3D conformer generation
-> Auto3D geometry optimization ensemble reduction
-> Auto3D configurational entropy Delta G ranking
```

Run it with:

```bash
dsvr run examples/test_molecules.smi \
  --config configs/auto3d_entropy_protocol.yaml \
  --outdir runs/auto3d_entropy_water_pH7
```

For a fast sanity-check, use:

```bash
dsvr run examples/test_molecules_minimal.smi \
  --config configs/auto3d_entropy_smoke.yaml \
  --outdir runs/auto3d_entropy_smoke
```

The corresponding Auto3D-native parameter template is
`configs/auto3d_entropy.auto3d.yaml`. DSVR writes lineage and ranking outputs
under `seeding/auto3d_protocol`, `auto3d_entropy`, and `ranking`.

Default implementation notes:

- Default pH: `7.0`
- Default solvent: `water`
- Default temperature: `298.15 K`
- Default initial seeder: RDKit ETKDG
- Optional seeder/prefilter: Auto3D
- Main decision engine: CREST/xTB
- Optional high-confidence refinement: CENSO
- Optional final QM rescoring: Psi4 or PySCF

## Scientific Warning

DSVR does not perform rigorous pH-dependent population calculations unless a
micro-pKa/proton chemical potential correction plugin is added. By default,
molscrub is used for practical pH/protomer candidate generation, then
CREST/xTB-derived free energies rank the generated candidates in the configured
solvent model.

Boltzmann populations are derived from relative free energies and must be
labeled with their scope:

- Comparable within the same formula/proton count.
- Approximate across different protonation/protomer states unless corrections
  are available.

RDKit tautomer canonicalization is not stability ranking. RDKit stereoisomer
enumeration is explicit and controlled by configuration. Auto3D can be useful
for seed generation or prefiltering, but it must not double-enumerate tautomers
or stereoisomers unless the user explicitly enables Auto3D internal enumeration.

## Quick Start

```bash
conda env create -f environment.yml
conda activate dsvr
python -m pip install -e ".[dev]"
dsvr doctor
dsvr run examples/test_molecules_minimal.smi --config configs/fast_smoke.yaml --outdir runs/smoke
```

For production local runs, start with the bounded balanced profile:

```bash
dsvr run examples/test_molecules.smi \
  --config configs/production_balanced.yaml \
  --outdir runs/production_balanced_water_pH7
```

`configs/exhaustive_debug.yaml` is intentionally expensive and may generate
very large variant and XYZ counts. Use it only for small molecules or debugging.

For a direct source-tree smoke check:

```bash
PYTHONPATH=src python -m dsvr.cli --help
PYTHONPATH=src python -m pytest
```

For GitHub Actions debugging, use:

```bash
python scripts/inspect_ci_run.py https://github.com/tevang/dominant-structural-variant-ranker/actions/runs/<run_id>
```

If `GH_TOKEN` is not set, the script will reuse `GITHUB_TOKEN` when available.

Short form:

```bash
make ci-log RUN=https://github.com/tevang/dominant-structural-variant-ranker/actions/runs/<run_id>
```

## Dependency Strategy

- Do not vendor third-party repositories.
- Install Python packages via conda or pip.
- Install external binaries via conda, official binaries, or user-managed
  modules.
- Use `dsvr doctor` to verify the environment before running physics-heavy
  workflows.

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
dsvr run examples/test_molecules_minimal.smi --config configs/fast_smoke.yaml --outdir runs/smoke
```

## Documentation

- [Architecture](docs/architecture.md)
- [Workflow](docs/workflow.md)
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
