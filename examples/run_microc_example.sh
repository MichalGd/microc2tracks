#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

bash "${REPO_DIR}/scripts/preflight_check.sh" \
  -c "${REPO_DIR}/config/config.conf" \
  -s "${REPO_DIR}/config/samplesheet.csv"

bash "${REPO_DIR}/scripts/microc2tracks.sh" \
  -c "${REPO_DIR}/config/config.conf" \
  -s "${REPO_DIR}/config/samplesheet.csv"
