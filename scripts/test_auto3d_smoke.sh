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
LOG_FILE="$LOG_DIR/auto3d-smoke-$(date +%Y%m%d-%H%M%S).log"
JOB_NAME="auto3d-smoke-${RUN_TAG:-$(date +%Y%m%d-%H%M%S)}"

exec > >(tee "$LOG_FILE") 2>&1
echo "Repository: $ROOT_DIR"
echo "Started: $(date --iso-8601=seconds)"
"$ROOT_DIR/.venv/bin/auto3d" run "$ROOT_DIR/examples/test_molecules_minimal.smi" \
  --job-name "$JOB_NAME" --engine AIMNET --no-gpu --k 1
echo "Completed: $(date --iso-8601=seconds)"
