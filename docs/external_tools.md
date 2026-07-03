# External Tools

DSVR coordinates external chemistry tools but does not vendor them. Install
Python packages with conda or pip. Install external binaries with conda,
official binary distributions, source builds, or user-managed modules on HPC
systems.

Use:

```bash
dsvr doctor
```

to check whether expected Python modules and executables are available.

## Dependency Strategy

- Do not vendor third-party repos.
- Install Python packages via conda/pip.
- Install external binaries via conda, official binaries, or user-managed
  modules.
- Use `dsvr doctor` to verify the environment.
- Keep tool versions in run provenance because rankings can change across
  versions.

## Tool Notes and References

| Tool | Role in DSVR | Install notes | URLs |
| --- | --- | --- | --- |
| RDKit | Core cheminformatics toolkit for reading, standardization hooks, tautomer enumeration, stereoisomer enumeration, and ETKDG seeding. | Prefer conda-forge: `conda install -c conda-forge rdkit`. | Docs: https://www.rdkit.org/docs/ |
| molscrub | Practical pH/protomer/protonation candidate generation. | Upstream currently documents pip-from-GitHub style installs; keep it optional and check with `dsvr doctor`. | GitHub: https://github.com/forlilab/molscrub ; Docs: https://molscrub.readthedocs.io/ |
| Auto3D | Optional conformer seeder or prefilter using neural-network potentials. | Install only when selected. The upstream project documents pip and conda-forge options. Disable internal tautomer/stereoisomer enumeration unless explicitly requested. | GitHub: https://github.com/isayevlab/Auto3D_pkg ; Docs: https://auto3d.readthedocs.io/ |
| AIMNet / aimnetcentral | Neural-network potential ecosystem relevant to Auto3D engines. | Usually pulled through the selected Auto3D configuration or installed as required by Auto3D. | GitHub: https://github.com/isayevlab/aimnetcentral |
| xTB | Semiempirical quantum engine for optimization, thermo, solvation, and CREST-backed workflows. | Prefer conda-forge where available, or official upstream binaries/source builds. Ensure `xtb` is on `PATH`. | GitHub: https://github.com/grimme-lab/xtb ; Docs: https://xtb-docs.readthedocs.io/en/latest/ |
| CREST | Main conformer search and ensemble reduction layer for the physics-heavy workflow. | Prefer conda-forge or official releases. Ensure `crest` is on `PATH`; CREST workflows often require xTB availability. | GitHub: https://github.com/crest-lab/crest ; Docs: https://crest-lab.github.io/crest-docs/ |
| CENSO | Optional high-confidence ensemble refinement and energetic sorting. | Install only for refinement workflows. Follow current CENSO documentation and ensure its own backend requirements are available. | Docs: https://xtb-docs.readthedocs.io/en/latest/CENSO_docs/censo.html |
| Psi4 | Optional final quantum-chemistry rescoring. | Prefer conda-forge or official Psi4 installation instructions. Ensure the Python module or executable is visible to the selected workflow. | Site: https://psicode.org/ ; Manual: https://psicode.org/psi4manual/master/index.html |
| PySCF | Optional Python-native final quantum-chemistry rescoring. | Install as an optional Python dependency only for PySCF workflows. | Site/docs: https://pyscf.org/ |

## Installation Examples

Core conda environment:

```bash
conda env create -f environment.yml
conda activate dominant-structural-variant-ranker
python -m pip install -e ".[dev]"
```

Optional Python packages through bootstrap flags:

```bash
scripts/bootstrap_conda.sh --with-molscrub --with-auto3d
scripts/bootstrap_mamba.sh --with-molscrub --with-auto3d
```

External binaries through conda-forge where available:

```bash
conda install -c conda-forge xtb crest
```

For cluster environments, prefer site-managed modules when available:

```bash
module load xtb
module load crest
module load psi4
dsvr doctor
```

## Runtime Behavior

Missing optional tools should not break import of `dsvr`. A missing tool should
produce a clear runtime error only when the selected workflow step requires that
tool.

Examples:

- A fast RDKit-only smoke workflow should not require Auto3D, xTB, CREST, CENSO,
  Psi4, or PySCF.
- A physics-heavy workflow must require xTB and CREST.
- A CENSO refinement workflow must require CENSO plus its configured backend
  requirements.
- A final Psi4/PySCF rescoring workflow must require the selected QM backend.
