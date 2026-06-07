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

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python "${SCRIPT_DIR}/sanitize_text_inputs.py" --kind config "${CONFIG}"
python "${SCRIPT_DIR}/sanitize_text_inputs.py" --kind samplesheet "${SAMPLESHEET}"

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
BWA_INDEX_PREFIX="${BWA_INDEX_PREFIX:-$REFERENCE_FASTA}"

if [ ! -f "${BWA_INDEX_PREFIX}.bwt.2bit.64" ]; then
  echo "WARNING: BWA-MEM2 index file not found: ${BWA_INDEX_PREFIX}.bwt.2bit.64" >&2
  echo "         Create it with: bwa-mem2 index \"${BWA_INDEX_PREFIX}\"" >&2
  echo "         If the index prefix differs from REFERENCE_FASTA, set BWA_INDEX_PREFIX in config.conf." >&2
fi

if [ "${FILTER_CANONICAL_CHROMS:-true}" = "true" ]; then
  retained_chroms="$(
    awk -v regex="${CANONICAL_CHROMS_REGEX:-^(chr)?([1-9][0-9]?|X|Y|M|MT)$}" '$1 ~ regex {count++} END{print count+0}' "${CHROM_SIZES}"
  )"
  [ "${retained_chroms}" -gt 0 ] || die "FILTER_CANONICAL_CHROMS=true but CANONICAL_CHROMS_REGEX retained no chromosomes from CHROM_SIZES"
  echo "Canonical chromosome filter retains ${retained_chroms} chromosome(s)."
fi

if [ "${RUN_HIC:-true}" = "true" ]; then
  [ -f "${JUICER_TOOLS_JAR:-}" ] || die "RUN_HIC=true but JUICER_TOOLS_JAR not found: ${JUICER_TOOLS_JAR:-unset}"
fi

python "${SCRIPT_DIR}/validate_samplesheet.py" --require-files "${SAMPLESHEET}"

echo "Checking chromosome-size first entry..."
head -n 1 "${CHROM_SIZES}" || true

echo "Preflight finished successfully."
