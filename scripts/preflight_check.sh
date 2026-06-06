#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  preflight_check.sh -c config/config.conf -s config/samplesheet.csv

Checks required tools, config values, reference files, and sample-sheet format.
USAGE
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

CONFIG=""
SAMPLESHEET=""

while getopts ":c:s:h" opt; do
  case "${opt}" in
    c) CONFIG="${OPTARG}" ;;
    s) SAMPLESHEET="${OPTARG}" ;;
    h) usage; exit 0 ;;
    :) die "Option -${OPTARG} requires an argument" ;;
    \?) die "Unknown option: -${OPTARG}" ;;
  esac
done

[ -n "${CONFIG}" ] || die "Missing -c config file"
[ -n "${SAMPLESHEET}" ] || die "Missing -s sample sheet"
[ -f "${CONFIG}" ] || die "Config not found: ${CONFIG}"
[ -f "${SAMPLESHEET}" ] || die "Sample sheet not found: ${SAMPLESHEET}"

# shellcheck source=/dev/null
source "${CONFIG}"

required_vars=(
  OUTDIR TMPDIR GENOME_ASSEMBLY REFERENCE_FASTA CHROM_SIZES
  MAPQ_THRESHOLD BASE_RESOLUTION RESOLUTIONS
  THREADS_ALIGN THREADS_SORT THREADS_MATRIX
)

for var in "${required_vars[@]}"; do
  [ -n "${!var:-}" ] || die "Config variable is empty or missing: ${var}"
done

required_tools=(
  awk
  bgzip
  bwa-mem2
  cooler
  fastp
  java
  multiqc
  pairix
  pairtools
  python
  samtools
  zcat
)

echo "Checking required tools..."
missing=0
for tool in "${required_tools[@]}"; do
  if command -v "${tool}" >/dev/null 2>&1; then
    printf '  OK   %s\n' "${tool}"
  else
    printf '  MISS %s\n' "${tool}"
    missing=1
  fi
done

optional_tools=(cooltools hictk mustache chromosight coolpup.py hicFindTADs)

echo "Checking optional downstream tools..."
for tool in "${optional_tools[@]}"; do
  if command -v "${tool}" >/dev/null 2>&1; then
    printf '  OK   %s\n' "${tool}"
  else
    printf '  MISS %s\n' "${tool}"
  fi
done

[ "${missing}" -eq 0 ] || die "Install missing required tools before running"

[ -f "${REFERENCE_FASTA}" ] || die "REFERENCE_FASTA not found: ${REFERENCE_FASTA}"
[ -f "${CHROM_SIZES}" ] || die "CHROM_SIZES not found: ${CHROM_SIZES}"

if [ "${RUN_HIC:-true}" = "true" ]; then
  [ -f "${JUICER_TOOLS_JAR:-}" ] || die "RUN_HIC=true but JUICER_TOOLS_JAR not found: ${JUICER_TOOLS_JAR:-unset}"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python "${SCRIPT_DIR}/validate_samplesheet.py" "${SAMPLESHEET}"

echo "Checking chromosome-size first entry..."
head -n 1 "${CHROM_SIZES}" || true

echo "Preflight finished successfully."

