#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  merge_replicates.sh -c config/config.conf -n merge_name -p pairs1.gz,pairs2.gz[,pairs3.gz...]

Merges filtered valid pairs with pairtools merge, indexes the merged pairs,
computes pairtools stats, and rebuilds .cool/.mcool/.hic matrices.
USAGE
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

log_msg() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

bool_true() {
  [ "${1:-false}" = "true" ] || [ "${1:-false}" = "1" ] || [ "${1:-false}" = "yes" ]
}

CONFIG=""
MERGE_NAME=""
PAIRS_CSV=""

while getopts ":c:n:p:h" opt; do
  case "${opt}" in
    c) CONFIG="${OPTARG}" ;;
    n) MERGE_NAME="${OPTARG}" ;;
    p) PAIRS_CSV="${OPTARG}" ;;
    h) usage; exit 0 ;;
    :) die "Option -${OPTARG} requires an argument" ;;
    \?) die "Unknown option: -${OPTARG}" ;;
  esac
done

[ -n "${CONFIG}" ] || die "Missing -c config file"
[ -n "${MERGE_NAME}" ] || die "Missing -n merge name"
[ -n "${PAIRS_CSV}" ] || die "Missing -p comma-separated pair files"
[ -f "${CONFIG}" ] || die "Config not found: ${CONFIG}"

# shellcheck source=/dev/null
source "${CONFIG}"

IFS=, read -r -a PAIR_FILES <<< "${PAIRS_CSV}"
[ "${#PAIR_FILES[@]}" -ge 2 ] || die "At least two pair files are required"

for pair_file in "${PAIR_FILES[@]}"; do
  [ -f "${pair_file}" ] || die "Pair file not found: ${pair_file}"
done

MERGE_DIR="${OUTDIR}/merged/${MERGE_NAME}"
PAIRS_DIR="${MERGE_DIR}/03_pairs"
MATRIX_DIR="${MERGE_DIR}/04_matrices"
LOG_DIR="${MERGE_DIR}/logs"
MERGE_TMP="${TMPDIR}/merged_${MERGE_NAME}"

mkdir -p "${PAIRS_DIR}" "${MATRIX_DIR}" "${LOG_DIR}" "${MERGE_TMP}"

MERGED_PAIRS="${PAIRS_DIR}/${MERGE_NAME}.valid.mapq${MAPQ_THRESHOLD}.pairs.gz"

if [ ! -s "${MERGED_PAIRS}" ]; then
  log_msg "${MERGE_NAME}: merging ${#PAIR_FILES[@]} pair files"
  pairtools merge \
    --nproc "${THREADS_MERGE}" \
    --memory "${PAIRTOOLS_MERGE_MEMORY}" \
    --tmpdir "${MERGE_TMP}" \
    -o "${MERGED_PAIRS}" \
    "${PAIR_FILES[@]}" \
    > "${LOG_DIR}/${MERGE_NAME}.pairtools_merge.log" 2>&1
else
  log_msg "${MERGE_NAME}: merged pairs exist, skipping merge"
fi

if [ ! -s "${MERGED_PAIRS}.px2" ]; then
  pairix "${MERGED_PAIRS}"
fi

pairtools stats "${MERGED_PAIRS}" > "${PAIRS_DIR}/${MERGE_NAME}.pairs.stats.txt"

RAW_COOL="${MATRIX_DIR}/${MERGE_NAME}.raw.${BASE_RESOLUTION}.cool"
NORM_COOL="${MATRIX_DIR}/${MERGE_NAME}.norm.${BASE_RESOLUTION}.cool"
RAW_MCOOL="${MATRIX_DIR}/${MERGE_NAME}.raw.mcool"
NORM_MCOOL="${MATRIX_DIR}/${MERGE_NAME}.norm.mcool"
JUICER_PAIRS="${PAIRS_DIR}/${MERGE_NAME}.valid.mapq${MAPQ_THRESHOLD}.juicer.pairs.gz"
RAW_HIC="${MATRIX_DIR}/${MERGE_NAME}.raw.hic"
NORM_HIC="${MATRIX_DIR}/${MERGE_NAME}.norm.hic"

if [ ! -s "${RAW_COOL}" ]; then
  log_msg "${MERGE_NAME}: building raw cooler"
  cooler cload pairix --assembly "${GENOME_ASSEMBLY}" \
    -p "${THREADS_MATRIX}" \
    "${CHROM_SIZES}:${BASE_RESOLUTION}" \
    "${MERGED_PAIRS}" \
    "${RAW_COOL}" \
    > "${LOG_DIR}/${MERGE_NAME}.cooler_cload.log" 2>&1
fi

if [ ! -s "${NORM_COOL}" ]; then
  log_msg "${MERGE_NAME}: balancing cooler"
  cp "${RAW_COOL}" "${NORM_COOL}"
  cooler balance -p "${THREADS_MATRIX}" -f "${NORM_COOL}" \
    > "${LOG_DIR}/${MERGE_NAME}.cooler_balance.log" 2>&1
fi

if bool_true "${RUN_MCOOL:-true}"; then
  if [ ! -s "${RAW_MCOOL}" ]; then
    cooler zoomify -p "${THREADS_MATRIX}" \
      -r "${RESOLUTIONS}" \
      -o "${RAW_MCOOL}" \
      "${RAW_COOL}" \
      > "${LOG_DIR}/${MERGE_NAME}.cooler_zoomify_raw.log" 2>&1
  fi

  if [ ! -s "${NORM_MCOOL}" ]; then
    cooler zoomify -p "${THREADS_MATRIX}" \
      -r "${RESOLUTIONS}" \
      --balance \
      -o "${NORM_MCOOL}" \
      "${RAW_COOL}" \
      > "${LOG_DIR}/${MERGE_NAME}.cooler_zoomify_norm.log" 2>&1
  fi
fi

cooler info "${RAW_COOL}" > "${MATRIX_DIR}/${MERGE_NAME}.raw.${BASE_RESOLUTION}.cool.info.json"

if bool_true "${RUN_HIC:-true}"; then
  [ -f "${JUICER_TOOLS_JAR}" ] || die "JUICER_TOOLS_JAR not found: ${JUICER_TOOLS_JAR}"

  if [ ! -s "${JUICER_PAIRS}" ]; then
    zcat "${MERGED_PAIRS}" \
      | awk 'BEGIN{OFS="\t"} /^## pairs format/ {print; next} /^#columns:/ {print; next} /^#/ {next} {print}' \
      | bgzip -@ "${THREADS_MATRIX}" \
      > "${JUICER_PAIRS}"
  fi

  if [ ! -s "${RAW_HIC}" ]; then
    java -Xmx"${JAVA_HEAP}" -jar "${JUICER_TOOLS_JAR}" pre \
      -n \
      -r "${RESOLUTIONS}" \
      "${JUICER_PAIRS}" \
      "${RAW_HIC}" \
      "${CHROM_SIZES}" \
      > "${LOG_DIR}/${MERGE_NAME}.juicer_pre.log" 2>&1
  fi

  if [ ! -s "${NORM_HIC}" ]; then
    cp "${RAW_HIC}" "${NORM_HIC}"
    java -Xmx"${JAVA_HEAP}" -jar "${JUICER_TOOLS_JAR}" addNorm \
      -k VC,VC_SQRT,KR,SCALE \
      "${NORM_HIC}" \
      > "${LOG_DIR}/${MERGE_NAME}.juicer_addNorm.log" 2>&1
  fi
fi

log_msg "${MERGE_NAME}: finished"

