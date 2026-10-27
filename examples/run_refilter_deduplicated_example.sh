#!/usr/bin/env bash
set -euo pipefail
# Existing deduplicated .pairs/.pairs.gz only. Original files/indexes are read-only.
# Review the reference, enzymes, distance cutoff and NEW output roots first.
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="${1:-${REPO_DIR}/config/config.conf}"
SHEET="${2:-${REPO_DIR}/config/samplesheet_refilter_example.csv}"
bash "${REPO_DIR}/scripts/preflight_check.sh" --mode refilter-deduplicated -c "${CONFIG}" -s "${SHEET}"
bash "${REPO_DIR}/scripts/microc2tracks.sh" --mode refilter-deduplicated -c "${CONFIG}" -s "${SHEET}"
