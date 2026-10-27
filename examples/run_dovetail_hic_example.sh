#!/usr/bin/env bash
set -euo pipefail
# Original restriction-based Dovetail Hi-C (DpnII), not Omni-C/LinkPrep.
# Edit the reference registry, config and CSV paths before running.
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="${1:-${REPO_DIR}/config/config.conf}"
SHEET="${2:-${REPO_DIR}/config/samplesheet_dovetail_example.csv}"
bash "${REPO_DIR}/scripts/preflight_check.sh" -c "${CONFIG}" -s "${SHEET}"
bash "${REPO_DIR}/scripts/microc2tracks.sh" -c "${CONFIG}" -s "${SHEET}"
