# Installation

DSVR is a wrapper/orchestrator. It does not vendor RDKit, molscrub, Auto3D, xTB, CREST, CENSO, Psi4, PySCF, or their source trees.

## Minimal Install

Use conda-forge for the core Python stack:

```bash
conda env create -f environment.yml
conda activate dsvr
python -m pip install -e ".[dev]"
dsvr doctor
```

This minimal environment supports input parsing, RDKit enumeration, ETKDG seeding, dry runs, reports, and CI-safe smoke tests. External physics steps require additional tools.

## Default `prepare-ligands` Workflow

The default `prepare-ligands` (LigPrep-like) workflow uses molscrub for pH/protomer
generation, Auto3D for tautomer energy triage, and Auto3D for final 3D conformer
generation. Both are required for any non-dry-run `prepare-ligands` invocation.

With `uv` (recommended):

```bash
uv sync --extra dev --extra auto3d --extra molscrub
source .venv/bin/activate
dsvr doctor
```

With conda/pip:

```bash
conda activate dsvr
python -m pip install -e ".[dev]"
python -m pip install Auto3D
python -m pip install git+https://github.com/forlilab/molscrub.git
dsvr doctor
```

For RDKit-only `--dry-run` checks or CI-safe smoke tests, the `auto3d` and
`molscrub` extras are not required.

## Bootstrap Helpers

The helper scripts create or update an environment named `dsvr`, install the package editable, and run `dsvr doctor`.

```bash
scripts/bootstrap_mamba.sh
scripts/bootstrap_conda.sh
```

Optional flags:

```bash
scripts/bootstrap_mamba.sh --with-molscrub --with-auto3d --with-xtb --with-crest
scripts/bootstrap_conda.sh --with-pyscf --with-psi4 --strict
```

Optional installs are best-effort by default. Use `--strict` if optional install failures should stop the bootstrap. The scripts do not modify shell profiles.

## Full Local Workstation Install

A fuller local setup may look like:

```bash
scripts/bootstrap_mamba.sh \
  --with-molscrub \
  --with-auto3d \
  --with-pyscf \
  --with-xtb \
  --with-crest
```

If `crest` is unavailable from conda-forge on your platform, install CREST from official release binaries or your package manager, then ensure `crest` and `xtb` are on `PATH`.

For the local `dsvr` conda environment used in this repository, we also install a conda activation hook that sets up both tools automatically. After `conda activate dsvr`, `xtb` and `crest` should already be on `PATH` and `XTBHOME` / `CRESTHOME` should be set for the session.

Check the result:

```bash
dsvr doctor
scripts/check_external_tools.sh --json --json-out doctor_report.json
```

## External Binaries

xTB and CREST may be installed through any of these routes:

- conda-forge, when packages are available for your platform;
- official upstream release binaries;
- system package managers;
- HPC module systems.

Manual binary fallback:

```bash
which xtb
xtb --version
which crest
crest --version
dsvr doctor
```

If commands are not found, add their installation directory to `PATH` in your current shell or load the relevant module. DSVR does not edit shell startup files.

## HPC Install With Modules

On clusters, prefer the site-provided modules for heavy external tools:

```bash
module load xtb
module load crest
module load psi4  # if available

conda env create -f environment.yml
conda activate dsvr
python -m pip install -e ".[dev]"
dsvr doctor
```

Use workflow configuration to avoid CPU oversubscription. For example, set `crest.nproc` to the scheduler allocation per CREST job and keep global workflow workers conservative.

## Auto3D / AIMNet2 GPU Note

Auto3D can use ML models such as AIMNet2. GPU support depends on the Auto3D/PyTorch/AIMNet installation and the local CUDA driver stack. Install GPU-enabled dependencies according to Auto3D and PyTorch guidance for your machine. CPU-only Auto3D may work but can be slower.

Auto3D is required for the default `prepare-ligands` workflow: it performs
tautomer energy triage (`--tauto-k`, `--tauto-window`) and final 3D conformer
generation. It is only optional for RDKit-only `--dry-run` checks and CI-safe
smoke tests. Install it with the `auto3d` extra:

```bash
uv sync --extra auto3d
# or
python -m pip install Auto3D
```

## molscrub

molscrub is optional at install time but required for the default pH/protomer candidate-generation step:

```bash
python -m pip install git+https://github.com/forlilab/molscrub.git
dsvr doctor
```

If molscrub is unavailable, fast smoke workflows can use fallback original-molecule candidates, but production pH/protomer enumeration should install molscrub.

## Verify

Always run:

```bash
dsvr doctor
```

Use strict mode when preparing production runs:

```bash
dsvr doctor --strict
```

Strict mode fails if required default-workflow checks are missing. Optional tools such as Auto3D, CENSO, Psi4, and PySCF are reported but not required unless the corresponding workflow step is enabled.
