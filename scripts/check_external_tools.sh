#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: scripts/check_external_tools.sh [--strict] [--json] [--json-out PATH]

Run DSVR environment/tool diagnostics through `dsvr doctor`.
EOF
}

ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --strict|--json)
      ARGS+=("$1")
      shift
      ;;
    --json-out)
      if [[ $# -lt 2 ]]; then
        echo "--json-out requires a path" >&2
        exit 2
      fi
      ARGS+=("$1" "$2")
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if command -v dsvr >/dev/null 2>&1; then
  dsvr doctor "${ARGS[@]}"
else
  python -m dsvr.cli doctor "${ARGS[@]}"
fi
