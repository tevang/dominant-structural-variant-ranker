#!/usr/bin/env bash
set -euo pipefail

WITH_AUTO3D=0
WITH_MOLSCRUB=0

for arg in "$@"; do
  case "$arg" in
    --with-auto3d) WITH_AUTO3D=1 ;;
    --with-molscrub) WITH_MOLSCRUB=1 ;;
    *) echo "Unknown option: $arg" >&2; exit 2 ;;
  esac
done

conda env create -f environment.yml
conda run -n dominant-structural-variant-ranker python -m pip install -e ".[dev]"

if [[ "$WITH_AUTO3D" == "1" ]]; then
  conda run -n dominant-structural-variant-ranker python -m pip install Auto3D
fi

if [[ "$WITH_MOLSCRUB" == "1" ]]; then
  conda run -n dominant-structural-variant-ranker python -m pip install molscrub
fi

