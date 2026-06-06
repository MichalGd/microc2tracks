#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Use assay=hic in the sample sheet. The first version uses the shared
# pairs-to-matrix path; add restriction-fragment filtering once enzyme metadata
# and fragment BED files are available.
bash "${REPO_DIR}/scripts/microc2tracks.sh" \
  -c "${REPO_DIR}/config/config.conf" \
  -s "${REPO_DIR}/config/samplesheet.csv"
