#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="dsvr"
STRICT=0
WITH_MOLSCRUB=0
WITH_AUTO3D=0
WITH_PYSCF=0
WITH_PSI4=0
WITH_XTB=0
WITH_CREST=0

usage() {
  cat <<'EOF'
Usage: scripts/bootstrap_conda.sh [options]

Create/update the dsvr conda environment and install this package editable.

Options:
  --with-molscrub   pip install git+https://github.com/forlilab/molscrub.git
  --with-auto3d     pip install Auto3D
  --with-pyscf      pip install pyscf
  --with-psi4       conda install psi4 -c conda-forge
  --with-xtb        conda install xtb -c conda-forge
  --with-crest      conda install crest -c conda-forge; print manual instructions if unavailable
  --strict          fail if an optional install fails
  -h, --help        show this help
EOF
}

for arg in "$@"; do
  case "$arg" in
    --with-molscrub) WITH_MOLSCRUB=1 ;;
    --with-auto3d) WITH_AUTO3D=1 ;;
    --with-pyscf) WITH_PYSCF=1 ;;
    --with-psi4) WITH_PSI4=1 ;;
    --with-xtb) WITH_XTB=1 ;;
    --with-crest) WITH_CREST=1 ;;
    --strict) STRICT=1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $arg" >&2; usage >&2; exit 2 ;;
  esac
done

run_optional() {
  local description=$1
  shift
  echo "Optional install: ${description}"
  if "$@"; then
    echo "Optional install succeeded: ${description}"
  else
    echo "Optional install failed: ${description}" >&2
    if [[ "$STRICT" == "1" ]]; then
      exit 1
    fi
    return 0
  fi
}

if ! command -v conda >/dev/null 2>&1; then
  echo "conda is not on PATH. Install conda/mambaforge or use scripts/bootstrap_mamba.sh." >&2
  exit 1
fi

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "Environment ${ENV_NAME} exists; updating from environment.yml."
  conda env update -n "$ENV_NAME" -f environment.yml --prune
else
  echo "Creating environment ${ENV_NAME} from environment.yml."
  conda env create -n "$ENV_NAME" -f environment.yml
fi

conda run -n "$ENV_NAME" python -m pip install -e ".[dev]"

if [[ "$WITH_MOLSCRUB" == "1" ]]; then
  run_optional "molscrub" \
    conda run -n "$ENV_NAME" python -m pip install \
    "git+https://github.com/forlilab/molscrub.git"
fi

if [[ "$WITH_AUTO3D" == "1" ]]; then
  run_optional "Auto3D" conda run -n "$ENV_NAME" python -m pip install Auto3D
fi

if [[ "$WITH_PYSCF" == "1" ]]; then
  run_optional "PySCF" conda run -n "$ENV_NAME" python -m pip install pyscf
fi

if [[ "$WITH_PSI4" == "1" ]]; then
  run_optional "Psi4 from conda-forge" conda install -n "$ENV_NAME" -c conda-forge -y psi4
fi

if [[ "$WITH_XTB" == "1" ]]; then
  run_optional "xTB from conda-forge" conda install -n "$ENV_NAME" -c conda-forge -y xtb
fi

if [[ "$WITH_CREST" == "1" ]]; then
  run_optional "CREST from conda-forge" conda install -n "$ENV_NAME" -c conda-forge -y crest
  if ! conda run -n "$ENV_NAME" bash -lc 'command -v crest >/dev/null 2>&1'; then
    cat >&2 <<'EOF'
CREST was not found after the optional install attempt.
Manual options:
  - install CREST from conda-forge if available for your platform;
  - install official release binaries and add them to PATH;
  - load an HPC module that provides crest and xtb.
EOF
  fi
fi

conda run -n "$ENV_NAME" dsvr doctor

cat <<EOF

Bootstrap complete.
Activate with:
  conda activate ${ENV_NAME}

This script did not modify shell profiles.
EOF
