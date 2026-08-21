#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$ROOT_DIR"

export PATH="$ROOT_DIR/.venv/bin:$PATH"
export XDG_CACHE_HOME="$ROOT_DIR/.uv/xdg-cache"
export AIMNET_CACHE_DIR="$ROOT_DIR/.uv/aimnet-cache"
export WARP_CACHE_PATH="$ROOT_DIR/.uv/warp-cache"
export CUDA_VISIBLE_DEVICES=""
mkdir -p "$XDG_CACHE_HOME" "$AIMNET_CACHE_DIR" "$WARP_CACHE_PATH"

LOG_DIR="$ROOT_DIR/logs/screen-tests"
mkdir -p "$LOG_DIR"
RUN_TAG="${RUN_TAG:-$(date +%Y%m%d-%H%M%S)}"
RUN_DIR="$ROOT_DIR/runs/standard-protocol-screen-$RUN_TAG"
LOG_FILE="$LOG_DIR/standard-protocol-$RUN_TAG.log"

exec > >(tee "$LOG_FILE") 2>&1
echo "Repository: $ROOT_DIR"
echo "Run directory: $RUN_DIR"
echo "Started: $(date --iso-8601=seconds)"
"$ROOT_DIR/.venv/bin/dsvr" prepare-ligands \
  "$ROOT_DIR/examples/test_molecules_minimal.smi" \
  --config "$ROOT_DIR/configs/ligprep_like_default.yaml" \
  --out "$RUN_DIR" --overwrite --no-resume
echo "Completed: $(date --iso-8601=seconds)"
echo "Status summary:"
"$ROOT_DIR/.venv/bin/dsvr" status "$RUN_DIR"
