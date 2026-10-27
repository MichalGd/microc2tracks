#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Generic restriction Hi-C. Declare the actual enzyme in the CSV/config; do not
# infer it from the assembly. my_reference in the example must be registered.
bash "${REPO_DIR}/scripts/microc2tracks.sh" \
  -c "${1:-${REPO_DIR}/config/config.conf}" \
  -s "${2:-${REPO_DIR}/config/samplesheet_hic_example.csv}"
