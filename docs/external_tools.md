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
| Uni-Pka | Default protomer generator for `prepare-ligands`: per-microspecies free energies and Boltzmann occupancies at the working pH, driving occupancy-based protomer selection. Executed through a container (Docker or Apptainer); see acquisition below. | No Python dependencies in the DSVR env. Acquire either the pre-built `unipka.sif` from Zenodo or build from the EasyDock recipe; set `protonation.unipka.container` to the `.sif` path or Docker image name and check with `dsvr doctor`. | EasyDock container & script: https://github.com/ci-lab-cz/easydock/tree/master/containers/unipka ; docs: https://easydock.readthedocs.io/en/latest/usage/#standalone-uni-pka-usage |
| molscrub | Optional alternative protomer generator, selected with `protonation.tool: molscrub` (kept for A/B comparison; selection uses the legacy heuristic plausibility score). | Upstream documents pip-from-GitHub style installs; install with the `molscrub` extra and check with `dsvr doctor`. | GitHub: https://github.com/forlilab/molscrub ; Docs: https://molscrub.readthedocs.io/ |
| Auto3D | Neural-network-potential conformer generation and energy triage. Required for the default `prepare-ligands` workflow (tautomer energy triage and final 3D conformer generation). Optional as a seeding method when `seeding.method` is `auto3d` or `both`. | Install with the `auto3d` extra. The upstream project documents pip and conda-forge options. Disable internal tautomer/stereoisomer enumeration unless explicitly requested. | GitHub: https://github.com/isayevlab/Auto3D_pkg ; Docs: https://auto3d.readthedocs.io/ |
| AIMNet / aimnetcentral | Neural-network potential ecosystem relevant to Auto3D engines. | Usually pulled through the selected Auto3D configuration or installed as required by Auto3D. | GitHub: https://github.com/isayevlab/aimnetcentral |
| xTB | Semiempirical quantum engine for optimization, thermo, solvation, and CREST-backed workflows. | Prefer conda-forge where available, or official upstream binaries/source builds. Ensure `xtb` is on `PATH`. | GitHub: https://github.com/grimme-lab/xtb ; Docs: https://xtb-docs.readthedocs.io/en/latest/ |
| CREST | Opt-in conformer search and ensemble validation. | Prefer conda-forge or official releases. Ensure `crest` is on `PATH`; CREST workflows often require xTB availability. | GitHub: https://github.com/crest-lab/crest ; Docs: https://crest-lab.github.io/crest-docs/ |
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

Acquiring the Uni-Pka container (required for the default workflow):

```bash
# Apptainer: download the pre-built image (see the Zenodo "Uni-Pka protonation
# container and model files" record 19627026) into containers/unipka.sif —
# the default `container: unipka` finds it there automatically —
# or point protonation.unipka.container at any path / Docker image name.
wget -O containers/unipka.sif "https://zenodo.org/records/19627026/files/unipka.sif?download=1"
# or build from the EasyDock recipe
#   git clone https://github.com/ci-lab-cz/easydock.git && cd easydock/containers/unipka
#   apptainer build unipka.sif unipka.def        # or: docker build -t unipka .
```

Note: every published Zenodo `unipka.sif` build (through 2026-04-19) bakes an
older `unipka.py` that lacks the occupancy flags DSVR needs (`-n`,
`--occupancy`, `--distribution-file`). DSVR therefore bind-mounts the current
EasyDock script — vendored at `containers/unipka.py` (BSD-3-Clause, see
`containers/README.md`) — over the image copy. The vendored copy carries one
DSVR patch on top of upstream (`unipka-shared-microstate-fix.patch`): without
it, molecules sharing microstates within one batch (e.g. a conjugate acid/base
pair) are silently dropped from the output. Set
`protonation.unipka.script_path` to another checkout of the EasyDock script, or
to `""` to disable the override (only correct for freshly built images; note
that an unpatched script still has the shared-microstate bug until upstream
fixes it).

`dsvr doctor` (with a global `-c/--config`) checks the runtime and image of the
configured container reference; without a config it probes the default
selection. The `unipka` row is required only when the configured protonation
stage uses Uni-Pka (`enabled` + `tool: unipka`); molscrub-selected or disabled
protonation demotes it to informational, so `dsvr doctor --strict` reflects the
actual workflow.

For cluster environments, prefer site-managed modules when available:

```bash
module load xtb
module load crest
module load psi4
module load apptainer
dsvr doctor
```

## Runtime Behavior

Missing optional tools should not break import of `dsvr`. A missing tool should
produce a clear runtime error only when the selected workflow step requires that
tool.

Examples:

- An RDKit-only `--dry-run` or CI-safe smoke workflow should not require Auto3D,
  Uni-Pka, molscrub, xTB, CREST, CENSO, Psi4, or PySCF.
- The default `prepare-ligands` workflow requires the Uni-Pka container
  (protomer generation) and Auto3D (tautomer energy triage and final 3D
  conformer generation). With `protonation.tool: molscrub`, molscrub replaces
  the Uni-Pka requirement.
- A physics-heavy workflow must require xTB and CREST.
- A CENSO refinement workflow must require CENSO plus its configured backend
  requirements.
- A final Psi4/PySCF rescoring workflow must require the selected QM backend.
