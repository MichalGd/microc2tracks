#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  merge_replicates.sh -c config/config.conf -n merge_name -r reference_id -p pairs1.gz,pairs2.gz[,pairs3.gz...]

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

outputs_complete() {
  local path

  for path in "$@"; do
    [ -s "${path}" ] || return 1
  done

  return 0
}

mark_done() {
  local done_file="$1"
  local step="$2"
  shift 2

  mkdir -p "$(dirname "${done_file}")"
  {
    printf 'step\t%s\n' "${step}"
    printf 'completed_at\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S')"
    printf 'config\t%s\n' "${CONFIG}"
    printf 'reference_id\t%s\n' "${REFERENCE_ID:-unknown}"
    printf 'assembly\t%s\n' "${GENOME_ASSEMBLY:-unknown}"
    printf 'outputs\t%s\n' "$*"
  } > "${done_file}.tmp.$$"
  mv -f "${done_file}.tmp.$$" "${done_file}"
}

step_done() {
  local done_file="$1"
  local step="$2"
  shift 2

  if outputs_complete "$@"; then
    if [ -s "${done_file}" ]; then
      if sentinel_matches_reference "${done_file}" "${REFERENCE_ID:-unknown}" "${GENOME_ASSEMBLY:-unknown}"; then
        return 0
      fi
      log_msg "${step}: completion metadata has a different or missing reference; rebuilding"
      return 1
    fi
    if bool_true "${BOOTSTRAP_SENTINELS:-true}" \
      && bool_true "${ALLOW_LEGACY_MM39_RESUME:-true}" \
      && [ "${REFERENCE_ID:-}" = "mm39" ]; then
      log_msg "${step}: adopting legacy mouse output and recording reference_id=mm39"
      mark_done "${done_file}" "${step}" "$@"
      return 0
    fi
  fi

  return 1
}

write_ucsc_hic_track() {
  local sample="$1"
  local hic_path="$2"
  local matrix_dir="$3"
  local hic_for_url
  local public_url

  if [ -n "${PUBLIC_HIC_DIR:-}" ]; then
    mkdir -p "${PUBLIC_HIC_DIR}"
    cp "${hic_path}" "${PUBLIC_HIC_DIR}/"
    hic_for_url="$(basename "${hic_path}")"
  else
    hic_for_url="$(basename "${hic_path}")"
  fi

  if [ -n "${PUBLIC_HIC_BASE_URL:-}" ]; then
    public_url="${PUBLIC_HIC_BASE_URL%/}/${hic_for_url}"
    cat > "${matrix_dir}/${sample}.ucsc.hic.track.txt" <<TRACK
# reference_id=${REFERENCE_ID} assembly=${GENOME_ASSEMBLY} browser_preset=${BROWSER_PRESET}
track type=hic name="${sample}" description="${sample} normalized Hi-C/Micro-C contacts" bigDataUrl=${public_url}
TRACK
  fi
}

prepare_matrix_chrom_sizes() {
  local output_chrom_sizes="$1"

  if bool_true "${FILTER_CANONICAL_CHROMS:-true}"; then
    awk -v regex="${CANONICAL_CHROMS_REGEX:-^(chr)?([1-9][0-9]?|X|Y|M|MT)$}" \
      'BEGIN{OFS="\t"} $1 ~ regex {print $1, $2}' \
      "${CHROM_SIZES}" > "${output_chrom_sizes}"
  else
    cp "${CHROM_SIZES}" "${output_chrom_sizes}"
  fi

  [ -s "${output_chrom_sizes}" ] || die "No chromosomes retained in ${output_chrom_sizes}; check CHROM_SIZES and CANONICAL_CHROMS_REGEX"
}

filter_pairs_to_chrom_sizes() {
  local chrom_sizes="$1"

  awk -v chrom_sizes="${chrom_sizes}" '
    BEGIN {
      FS = OFS = "\t"
      while ((getline line < chrom_sizes) > 0) {
        split(line, fields, "\t")
        keep[fields[1]] = 1
      }
      close(chrom_sizes)
    }
    /^#/ { print; next }
    (($2 in keep) && ($4 in keep)) { print }
  '
}

cleanup_merge_intermediates() {
  if ! bool_true "${KEEP_RAW_MERGED_PAIRS:-false}"; then
    rm -f "${RAW_MERGED_PAIRS}"
  fi

  if ! bool_true "${KEEP_JUICER_PAIRS:-false}"; then
    rm -f "${JUICER_PAIRS}"
  fi

  if ! bool_true "${KEEP_SINGLE_RES_COOL:-false}"; then
    rm -f "${RAW_COOL}" "${NORM_COOL}"
  fi

  if ! bool_true "${KEEP_RAW_MCOOL:-false}"; then
    rm -f "${RAW_MCOOL}"
  fi

  if ! bool_true "${KEEP_RAW_HIC:-false}"; then
    rm -f "${RAW_HIC}"
  fi

  if bool_true "${CLEAN_TMP_ON_SUCCESS:-true}"; then
    rm -rf "${MERGE_TMP}"
  fi
}

CONFIG=""
MERGE_NAME=""
PAIRS_CSV=""
REFERENCE_REQUEST=""

while getopts ":c:n:p:r:h" opt; do
  case "${opt}" in
    c) CONFIG="${OPTARG}" ;;
    n) MERGE_NAME="${OPTARG}" ;;
    p) PAIRS_CSV="${OPTARG}" ;;
    r) REFERENCE_REQUEST="${OPTARG}" ;;
    h) usage; exit 0 ;;
    :) die "Option -${OPTARG} requires an argument" ;;
    \?) die "Unknown option: -${OPTARG}" ;;
  esac
done

[ -n "${CONFIG}" ] || die "Missing -c config file"
[ -n "${MERGE_NAME}" ] || die "Missing -n merge name"
[ -n "${PAIRS_CSV}" ] || die "Missing -p comma-separated pair files"
[ -f "${CONFIG}" ] || die "Config not found: ${CONFIG}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=reference_sentinel.sh
source "${SCRIPT_DIR}/reference_sentinel.sh"
python "${SCRIPT_DIR}/sanitize_text_inputs.py" --kind config "${CONFIG}"

# shellcheck source=/dev/null
source "${CONFIG}"

TAB=$'\t'
REFERENCE_REGISTRY="${REFERENCE_REGISTRY:-${SCRIPT_DIR}/../config/references.tsv}"
if [[ "${REFERENCE_REGISTRY}" != /* ]] && [ ! -f "${REFERENCE_REGISTRY}" ]; then
  if [ -f "$(dirname "${CONFIG}")/${REFERENCE_REGISTRY}" ]; then
    REFERENCE_REGISTRY="$(dirname "${CONFIG}")/${REFERENCE_REGISTRY}"
  elif [ -f "${SCRIPT_DIR}/../config/$(basename "${REFERENCE_REGISTRY}")" ]; then
    REFERENCE_REGISTRY="${SCRIPT_DIR}/../config/$(basename "${REFERENCE_REGISTRY}")"
  fi
fi
[ -f "${REFERENCE_REGISTRY}" ] || die "REFERENCE_REGISTRY not found: ${REFERENCE_REGISTRY}"
DEFAULT_REFERENCE_ID="${DEFAULT_REFERENCE_ID:-${GENOME_ASSEMBLY:-mm39}}"
REFERENCE_REQUEST="${REFERENCE_REQUEST:-${DEFAULT_REFERENCE_ID}}"
mapfile -t reference_fields < <(
  python "${SCRIPT_DIR}/reference_registry.py" resolve \
    --registry "${REFERENCE_REGISTRY}" "${REFERENCE_REQUEST}"
)
[ "${#reference_fields[@]}" -eq 11 ] || die "Could not resolve reference ${REFERENCE_REQUEST}"
REFERENCE_ID="${reference_fields[0]}"
REFERENCE_SPECIES="${reference_fields[1]}"
GENOME_ASSEMBLY="${reference_fields[2]}"
REFERENCE_FASTA_RESOLVED="${reference_fields[3]}"
BWA_INDEX_PREFIX_RESOLVED="${reference_fields[4]}"
CHROM_SIZES_RESOLVED="${reference_fields[5]}"
CANONICAL_CHROMS_REGEX_RESOLVED="${reference_fields[6]}"
BROWSER_PRESET="${reference_fields[7]}"
PHASING_TRACK_RESOLVED="${reference_fields[8]}"
ANNOTATION_METADATA="${reference_fields[9]}"
BLACKLIST_METADATA="${reference_fields[10]}"

mapfile -t default_fields < <(
  python "${SCRIPT_DIR}/reference_registry.py" resolve \
    --registry "${REFERENCE_REGISTRY}" "${DEFAULT_REFERENCE_ID}"
)
DEFAULT_REFERENCE_ID="${default_fields[0]}"
if [ "${REFERENCE_ID}" = "${DEFAULT_REFERENCE_ID}" ] \
  && [ -n "${REFERENCE_FASTA:-}" ] && [ -n "${CHROM_SIZES:-}" ]; then
  REFERENCE_FASTA_RESOLVED="${REFERENCE_FASTA}"
  BWA_INDEX_PREFIX_RESOLVED="${BWA_INDEX_PREFIX:-${REFERENCE_FASTA}}"
  CHROM_SIZES_RESOLVED="${CHROM_SIZES}"
  [ -z "${CANONICAL_CHROMS_REGEX:-}" ] || CANONICAL_CHROMS_REGEX_RESOLVED="${CANONICAL_CHROMS_REGEX}"
  [ -z "${PHASING_TRACK:-}" ] || PHASING_TRACK_RESOLVED="${PHASING_TRACK}"
fi
REFERENCE_FASTA="${REFERENCE_FASTA_RESOLVED}"
BWA_INDEX_PREFIX="${BWA_INDEX_PREFIX_RESOLVED}"
CHROM_SIZES="${CHROM_SIZES_RESOLVED}"
CANONICAL_CHROMS_REGEX="${CANONICAL_CHROMS_REGEX_RESOLVED}"
PHASING_TRACK="${PHASING_TRACK_RESOLVED}"
[ -f "${REFERENCE_FASTA}" ] || die "Reference FASTA not found for ${REFERENCE_ID}: ${REFERENCE_FASTA}"
[ -f "${CHROM_SIZES}" ] || die "Chromosome sizes not found for ${REFERENCE_ID}: ${CHROM_SIZES}"

MCOOL_RESOLUTIONS="${MCOOL_RESOLUTIONS:-${RESOLUTIONS}}"
HIC_RESOLUTIONS="${HIC_RESOLUTIONS:-${RESOLUTIONS}}"
THREADS_HIC_NORM="${THREADS_HIC_NORM:-24}"
HIC_NORMALIZATIONS="${HIC_NORMALIZATIONS:-VC,VC_SQRT,KR,SCALE}"
BOOTSTRAP_SENTINELS="${BOOTSTRAP_SENTINELS:-true}"
ALLOW_LEGACY_MM39_RESUME="${ALLOW_LEGACY_MM39_RESUME:-true}"

IFS=, read -r -a PAIR_FILES <<< "${PAIRS_CSV}"
[ "${#PAIR_FILES[@]}" -ge 2 ] || die "At least two pair files are required"

for pair_file in "${PAIR_FILES[@]}"; do
  [ -f "${pair_file}" ] || die "Pair file not found: ${pair_file}"
done

python "${SCRIPT_DIR}/validate_pairs_reference.py" \
  --assembly "${GENOME_ASSEMBLY}" "${PAIR_FILES[@]}"

MERGE_DIR="${OUTDIR}/merged/${MERGE_NAME}"
PAIRS_DIR="${MERGE_DIR}/03_pairs"
MATRIX_DIR="${MERGE_DIR}/04_matrices"
LOG_DIR="${MERGE_DIR}/logs"
MERGE_TMP="${TMPDIR}/merged_${MERGE_NAME}"
DONE_DIR="${LOG_DIR}/done"

mkdir -p "${PAIRS_DIR}" "${MATRIX_DIR}" "${LOG_DIR}" "${DONE_DIR}" "${MERGE_TMP}"
MERGE_START_EPOCH="$(date '+%s')"
MERGE_STATUS_FILE="${LOG_DIR}/${MERGE_NAME}.status.tsv"
printf 'merge_group\treference_id\tassembly\tstatus\tstart_epoch\tend_epoch\truntime_seconds\n' > "${MERGE_STATUS_FILE}"
printf '%s\t%s\t%s\trunning\t%s\t\t\n' "${MERGE_NAME}" "${REFERENCE_ID}" "${GENOME_ASSEMBLY}" "${MERGE_START_EPOCH}" >> "${MERGE_STATUS_FILE}"
{
  printf 'reference_id\t%s\n' "${REFERENCE_ID}"
  printf 'species\t%s\n' "${REFERENCE_SPECIES}"
  printf 'assembly\t%s\n' "${GENOME_ASSEMBLY}"
  printf 'fasta\t%s\n' "${REFERENCE_FASTA}"
  printf 'chrom_sizes\t%s\n' "${CHROM_SIZES}"
  printf 'browser_preset\t%s\n' "${BROWSER_PRESET}"
  printf 'phasing_track\t%s\n' "${PHASING_TRACK}"
  printf 'annotation_metadata\t%s\n' "${ANNOTATION_METADATA}"
  printf 'blacklist_metadata\t%s\n' "${BLACKLIST_METADATA}"
} > "${LOG_DIR}/${MERGE_NAME}.reference.tsv"

MATRIX_CHROM_SIZES="${PAIRS_DIR}/${MERGE_NAME}.matrix.chrom.sizes"
prepare_matrix_chrom_sizes "${MATRIX_CHROM_SIZES}"

RAW_MERGED_PAIRS="${PAIRS_DIR}/${MERGE_NAME}.merged.valid.mapq${MAPQ_THRESHOLD}.pairs.gz"
MERGED_PAIRS="${PAIRS_DIR}/${MERGE_NAME}.valid.mapq${MAPQ_THRESHOLD}.pairs.gz"
RAW_MERGED_DONE="${DONE_DIR}/raw_merged_pairs.done"
MERGED_DONE="${DONE_DIR}/merged_pairs.done"
MERGED_INDEX_DONE="${DONE_DIR}/merged_pairix.done"
MERGED_STATS_DONE="${DONE_DIR}/merged_pairs_stats.done"

if step_done "${RAW_MERGED_DONE}" "${MERGE_NAME}:raw_merged_pairs" "${RAW_MERGED_PAIRS}"; then
  log_msg "${MERGE_NAME}: raw merged pairs exist, skipping merge"
else
  log_msg "${MERGE_NAME}: merging ${#PAIR_FILES[@]} pair files"
  RAW_MERGED_TMP="${RAW_MERGED_PAIRS%.pairs.gz}.tmp.$$.pairs.gz"
  rm -f "${RAW_MERGED_TMP}"
  pairtools merge \
    --nproc "${THREADS_MERGE}" \
    --memory "${PAIRTOOLS_MERGE_MEMORY}" \
    --tmpdir "${MERGE_TMP}" \
    -o "${RAW_MERGED_TMP}" \
    "${PAIR_FILES[@]}" \
    > "${LOG_DIR}/${MERGE_NAME}.pairtools_merge.log" 2>&1
  mv -f "${RAW_MERGED_TMP}" "${RAW_MERGED_PAIRS}"
  mark_done "${RAW_MERGED_DONE}" "${MERGE_NAME}:raw_merged_pairs" "${RAW_MERGED_PAIRS}"
fi

if step_done "${MERGED_DONE}" "${MERGE_NAME}:merged_pairs" "${MERGED_PAIRS}"; then
  log_msg "${MERGE_NAME}: filtered merged pairs exist, skipping chromosome filter"
else
  log_msg "${MERGE_NAME}: filtering merged pairs to matrix chromosomes"
  MERGED_TMP="${MERGED_PAIRS%.pairs.gz}.tmp.$$.pairs.gz"
  rm -f "${MERGED_TMP}"
  zcat "${RAW_MERGED_PAIRS}" \
    | filter_pairs_to_chrom_sizes "${MATRIX_CHROM_SIZES}" \
    | bgzip -@ "${THREADS_MATRIX}" \
    > "${MERGED_TMP}"
  mv -f "${MERGED_TMP}" "${MERGED_PAIRS}"
  mark_done "${MERGED_DONE}" "${MERGE_NAME}:merged_pairs" "${MERGED_PAIRS}"
fi

if step_done "${MERGED_INDEX_DONE}" "${MERGE_NAME}:merged_pairix" "${MERGED_PAIRS}.px2"; then
  :
else
  pairix "${MERGED_PAIRS}"
  mark_done "${MERGED_INDEX_DONE}" "${MERGE_NAME}:merged_pairix" "${MERGED_PAIRS}.px2"
fi

if step_done "${MERGED_STATS_DONE}" "${MERGE_NAME}:merged_pairs_stats" "${PAIRS_DIR}/${MERGE_NAME}.pairs.stats.txt"; then
  :
else
  MERGED_STATS_TMP="${PAIRS_DIR}/${MERGE_NAME}.pairs.stats.txt.tmp.$$"
  pairtools stats "${MERGED_PAIRS}" > "${MERGED_STATS_TMP}"
  mv -f "${MERGED_STATS_TMP}" "${PAIRS_DIR}/${MERGE_NAME}.pairs.stats.txt"
  mark_done "${MERGED_STATS_DONE}" "${MERGE_NAME}:merged_pairs_stats" "${PAIRS_DIR}/${MERGE_NAME}.pairs.stats.txt"
fi

RAW_COOL="${MATRIX_DIR}/${MERGE_NAME}.raw.${BASE_RESOLUTION}.cool"
NORM_COOL="${MATRIX_DIR}/${MERGE_NAME}.norm.${BASE_RESOLUTION}.cool"
RAW_MCOOL="${MATRIX_DIR}/${MERGE_NAME}.raw.mcool"
NORM_MCOOL="${MATRIX_DIR}/${MERGE_NAME}.norm.mcool"
JUICER_PAIRS="${PAIRS_DIR}/${MERGE_NAME}.valid.mapq${MAPQ_THRESHOLD}.juicer.pairs.gz"
RAW_HIC="${MATRIX_DIR}/${MERGE_NAME}.raw.hic"
NORM_HIC="${MATRIX_DIR}/${MERGE_NAME}.norm.hic"
RAW_COOL_DONE="${DONE_DIR}/raw_cool.done"
NORM_COOL_DONE="${DONE_DIR}/norm_cool.done"
RAW_MCOOL_DONE="${DONE_DIR}/raw_mcool.done"
NORM_MCOOL_DONE="${DONE_DIR}/norm_mcool.done"
JUICER_PAIRS_DONE="${DONE_DIR}/juicer_pairs.done"
RAW_HIC_DONE="${DONE_DIR}/raw_hic.done"
NORM_HIC_DONE="${DONE_DIR}/norm_hic.done"

if bool_true "${RUN_MCOOL:-true}" \
  && ! bool_true "${KEEP_SINGLE_RES_COOL:-false}" \
  && ! bool_true "${KEEP_RAW_MCOOL:-false}" \
  && step_done "${NORM_MCOOL_DONE}" "${MERGE_NAME}:norm_mcool" "${NORM_MCOOL}"; then
  log_msg "${MERGE_NAME}: normalized mcool exists, skipping cooler/mcool rebuild"
else

if step_done "${RAW_COOL_DONE}" "${MERGE_NAME}:raw_cool" "${RAW_COOL}"; then
  log_msg "${MERGE_NAME}: raw cooler exists, skipping"
else
  log_msg "${MERGE_NAME}: building raw cooler"
  RAW_COOL_TMP="${RAW_COOL}.tmp.$$"
  rm -f "${RAW_COOL_TMP}"
  cooler cload pairix --assembly "${GENOME_ASSEMBLY}" \
    -p "${THREADS_MATRIX}" \
    "${MATRIX_CHROM_SIZES}:${BASE_RESOLUTION}" \
    "${MERGED_PAIRS}" \
    "${RAW_COOL_TMP}" \
    > "${LOG_DIR}/${MERGE_NAME}.cooler_cload.log" 2>&1
  mv -f "${RAW_COOL_TMP}" "${RAW_COOL}"
  mark_done "${RAW_COOL_DONE}" "${MERGE_NAME}:raw_cool" "${RAW_COOL}"
fi

if step_done "${NORM_COOL_DONE}" "${MERGE_NAME}:norm_cool" "${NORM_COOL}"; then
  log_msg "${MERGE_NAME}: normalized cooler exists, skipping"
else
  log_msg "${MERGE_NAME}: balancing cooler"
  NORM_COOL_TMP="${NORM_COOL}.tmp.$$"
  rm -f "${NORM_COOL_TMP}"
  cp "${RAW_COOL}" "${NORM_COOL_TMP}"
  cooler balance -p "${THREADS_MATRIX}" -f "${NORM_COOL_TMP}" \
    > "${LOG_DIR}/${MERGE_NAME}.cooler_balance.log" 2>&1
  mv -f "${NORM_COOL_TMP}" "${NORM_COOL}"
  mark_done "${NORM_COOL_DONE}" "${MERGE_NAME}:norm_cool" "${NORM_COOL}"
fi

if bool_true "${RUN_MCOOL:-true}"; then
  if step_done "${RAW_MCOOL_DONE}" "${MERGE_NAME}:raw_mcool" "${RAW_MCOOL}"; then
    log_msg "${MERGE_NAME}: raw mcool exists, skipping"
  else
    log_msg "${MERGE_NAME}: building raw mcool"
    RAW_MCOOL_TMP="${RAW_MCOOL}.tmp.$$"
    rm -f "${RAW_MCOOL_TMP}"
    cooler zoomify -p "${THREADS_MATRIX}" \
      -r "${MCOOL_RESOLUTIONS}" \
      -o "${RAW_MCOOL_TMP}" \
      "${RAW_COOL}" \
      > "${LOG_DIR}/${MERGE_NAME}.cooler_zoomify_raw.log" 2>&1
    mv -f "${RAW_MCOOL_TMP}" "${RAW_MCOOL}"
    mark_done "${RAW_MCOOL_DONE}" "${MERGE_NAME}:raw_mcool" "${RAW_MCOOL}"
  fi

  if step_done "${NORM_MCOOL_DONE}" "${MERGE_NAME}:norm_mcool" "${NORM_MCOOL}"; then
    log_msg "${MERGE_NAME}: normalized mcool exists, skipping"
  else
    log_msg "${MERGE_NAME}: building balanced mcool"
    NORM_MCOOL_TMP="${NORM_MCOOL}.tmp.$$"
    rm -f "${NORM_MCOOL_TMP}"
    cooler zoomify -p "${THREADS_MATRIX}" \
      -r "${MCOOL_RESOLUTIONS}" \
      --balance \
      -o "${NORM_MCOOL_TMP}" \
      "${RAW_COOL}" \
      > "${LOG_DIR}/${MERGE_NAME}.cooler_zoomify_norm.log" 2>&1
    mv -f "${NORM_MCOOL_TMP}" "${NORM_MCOOL}"
    mark_done "${NORM_MCOOL_DONE}" "${MERGE_NAME}:norm_mcool" "${NORM_MCOOL}"
  fi
fi

if [ -s "${RAW_COOL}" ]; then
  cooler info "${RAW_COOL}" > "${MATRIX_DIR}/${MERGE_NAME}.raw.${BASE_RESOLUTION}.cool.info.json"
fi

fi

if bool_true "${RUN_HIC:-true}"; then
  [ -f "${JUICER_TOOLS_JAR}" ] || die "JUICER_TOOLS_JAR not found: ${JUICER_TOOLS_JAR}"

  if step_done "${NORM_HIC_DONE}" "${MERGE_NAME}:norm_hic" "${NORM_HIC}"; then
    log_msg "${MERGE_NAME}: normalized hic exists, skipping hic generation"
    write_ucsc_hic_track "${MERGE_NAME}" "${NORM_HIC}" "${MATRIX_DIR}"
  else

  if step_done "${JUICER_PAIRS_DONE}" "${MERGE_NAME}:juicer_pairs" "${JUICER_PAIRS}"; then
    log_msg "${MERGE_NAME}: Juicer-compatible pairs exist, skipping"
  else
    log_msg "${MERGE_NAME}: writing Juicer-compatible pairs"
    JUICER_PAIRS_TMP="${JUICER_PAIRS%.pairs.gz}.tmp.$$.pairs.gz"
    rm -f "${JUICER_PAIRS_TMP}"
    zcat "${MERGED_PAIRS}" \
      | awk 'BEGIN{OFS="\t"} /^## pairs format/ {print; next} /^#columns:/ {print; next} /^#/ {next} {print}' \
      | bgzip -@ "${THREADS_MATRIX}" \
      > "${JUICER_PAIRS_TMP}"
    mv -f "${JUICER_PAIRS_TMP}" "${JUICER_PAIRS}"
    mark_done "${JUICER_PAIRS_DONE}" "${MERGE_NAME}:juicer_pairs" "${JUICER_PAIRS}"
  fi

  if step_done "${RAW_HIC_DONE}" "${MERGE_NAME}:raw_hic" "${RAW_HIC}"; then
    log_msg "${MERGE_NAME}: raw hic exists, skipping"
  else
    log_msg "${MERGE_NAME}: building raw hic"
    RAW_HIC_TMP="${RAW_HIC}.tmp.$$"
    rm -f "${RAW_HIC_TMP}"
    java -Xmx"${JAVA_HEAP}" -jar "${JUICER_TOOLS_JAR}" pre \
      -n \
      -r "${HIC_RESOLUTIONS}" \
      "${JUICER_PAIRS}" \
      "${RAW_HIC_TMP}" \
      "${MATRIX_CHROM_SIZES}" \
      > "${LOG_DIR}/${MERGE_NAME}.juicer_pre.log" 2>&1
    mv -f "${RAW_HIC_TMP}" "${RAW_HIC}"
    mark_done "${RAW_HIC_DONE}" "${MERGE_NAME}:raw_hic" "${RAW_HIC}"
  fi

    log_msg "${MERGE_NAME}: adding Juicer normalizations"
    NORM_HIC_TMP="${NORM_HIC}.tmp.$$"
    rm -f "${NORM_HIC_TMP}"
    cp "${RAW_HIC}" "${NORM_HIC_TMP}"
    java -Xmx"${JAVA_HEAP}" -jar "${JUICER_TOOLS_JAR}" addNorm \
      -j "${THREADS_HIC_NORM}" \
      -k "${HIC_NORMALIZATIONS}" \
      "${NORM_HIC_TMP}" \
      > "${LOG_DIR}/${MERGE_NAME}.juicer_addNorm.log" 2>&1
    mv -f "${NORM_HIC_TMP}" "${NORM_HIC}"
    mark_done "${NORM_HIC_DONE}" "${MERGE_NAME}:norm_hic" "${NORM_HIC}"

  write_ucsc_hic_track "${MERGE_NAME}" "${NORM_HIC}" "${MATRIX_DIR}"
  fi
fi

PRELIM_DONE="${DONE_DIR}/prelim_downstream.done"
if bool_true "${RUN_PRELIM_DOWNSTREAM:-true}"; then
  if [ -s "${PRELIM_DONE}" ]; then
    log_msg "${MERGE_NAME}: preliminary downstream already marked complete, skipping"
  elif [ -s "${NORM_MCOOL}" ]; then
    log_msg "${MERGE_NAME}: running preliminary downstream analyses"
    MICROC2TRACKS_PHASING_TRACK_OVERRIDE="${PHASING_TRACK}" \
    MICROC2TRACKS_REFERENCE_ID="${REFERENCE_ID}" \
    bash "${SCRIPT_DIR}/run_downstream.sh" \
      -l \
      -c "${CONFIG}" \
      -s "${MERGE_NAME}" \
      -m "${NORM_MCOOL}" \
      -o "${MERGE_DIR}" \
      > "${LOG_DIR}/${MERGE_NAME}.prelim_downstream.log" 2>&1 || true
    mark_done "${PRELIM_DONE}" "${MERGE_NAME}:prelim_downstream" "${NORM_MCOOL}"
  else
    log_msg "${MERGE_NAME}: normalized mcool not found; skipping preliminary downstream"
  fi
fi

MERGE_END_EPOCH="$(date '+%s')"
printf 'merge_group\treference_id\tassembly\tstatus\tstart_epoch\tend_epoch\truntime_seconds\n' > "${MERGE_STATUS_FILE}"
printf '%s\t%s\t%s\tsuccess\t%s\t%s\t%s\n' \
  "${MERGE_NAME}" "${REFERENCE_ID}" "${GENOME_ASSEMBLY}" "${MERGE_START_EPOCH}" "${MERGE_END_EPOCH}" "$((MERGE_END_EPOCH - MERGE_START_EPOCH))" \
  >> "${MERGE_STATUS_FILE}"

cleanup_merge_intermediates

log_msg "${MERGE_NAME}: finished"
