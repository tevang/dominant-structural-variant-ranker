#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$ROOT_DIR"

export PATH="$ROOT_DIR/.venv/bin:$PATH"
LOG_DIR="$ROOT_DIR/logs/screen-tests"
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/doctor-$(date +%Y%m%d-%H%M%S).log"

exec > >(tee "$LOG_FILE") 2>&1
echo "Repository: $ROOT_DIR"
echo "Started: $(date --iso-8601=seconds)"
"$ROOT_DIR/.venv/bin/dsvr" doctor --strict
echo "Completed: $(date --iso-8601=seconds)"
