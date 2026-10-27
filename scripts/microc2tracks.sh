#!/usr/bin/env bash
set -Eeuo pipefail

# Both assays use the same integrated implementation and resolved JSON run plan.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "${PYTHON:-python}" "${SCRIPT_DIR}/pipeline.py" run "$@"
