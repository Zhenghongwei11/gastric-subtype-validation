#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$PWD}"
RUN_ID="${1:-full-$(date -u +%Y%m%dT%H%M%SZ)}"

python3 "$SCRIPT_DIR/workflow.py" --mode full --project-root "$PROJECT_ROOT" --run-id "$RUN_ID"