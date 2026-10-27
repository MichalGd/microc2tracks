#!/usr/bin/env bash
set -Eeuo pipefail

# Both assays use the same integrated implementation and resolved JSON run plan.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONDONTWRITEBYTECODE=1
exec "${PYTHON:-python}" -B "${SCRIPT_DIR}/pipeline.py" run "$@"
