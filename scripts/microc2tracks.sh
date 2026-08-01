#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'USAGE'
Usage:
  microc2tracks.sh -c config/config.conf -s config/samplesheet.csv

Runs paired-end FASTQ through:
  fastp -> bwa-mem2 -> pairtools -> cooler/mcool -> Juicer .hic

Sample sheet columns:
  sample,assay,reference_genome,condition,biological_replicate,technical_replicate,fastq_r1,fastq_r2

reference_genome accepts controlled aliases from config/references.tsv. Legacy
sheets without the column use DEFAULT_REFERENCE_ID (mm39 by default) with a warning.
USAGE
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

log_msg() {
  local message="$1"
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${message}"
}

on_error() {
  local exit_code="$1"
  local line_no="$2"
  local command="$3"
  echo "ERROR: microc2tracks failed with exit code ${exit_code} at line ${line_no}" >&2
  echo "ERROR: failed command: ${command}" >&2
  echo "ERROR: check the newest log under ${OUTDIR:-results}/*/logs/ or ${OUTDIR:-results}/merged/*/logs/" >&2
  exit "${exit_code}"
}

trap 'on_error "$?" "$LINENO" "$BASH_COMMAND"' ERR

bool_true() {
  [ "${1:-false}" = "true" ] || [ "${1:-false}" = "1" ] || [ "${1:-false}" = "yes" ]
}

sanitize_id() {
  printf '%s' "$1" | tr -c 'A-Za-z0-9_.-' '_'
}

clean_csv_field() {
  printf '%s' "$1" \
    | tr -d '\r\n' \
    | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' -e 's/^"//' -e 's/"$//'
}

is_positive_int() {
  [[ "${1:-}" =~ ^[1-9][0-9]*$ ]]
}

is_nonnegative_int() {
  [[ "${1:-}" =~ ^[0-9]+$ ]]
}

require_positive_int() {
  local name="$1"
  local value="$2"
  is_positive_int "${value}" || die "${name} must be a positive integer, got '${value}'"
}

acquire_slot() {
  local name="$1"
  local limit="$2"
  local lock_dir="${LOCK_ROOT}/${name}"
  local slot

  mkdir -p "${lock_dir}"

  while true; do
    for slot in $(seq 1 "${limit}"); do
      if mkdir "${lock_dir}/slot_${slot}" 2>/dev/null; then
        printf '%s\n' "${lock_dir}/slot_${slot}"
        return 0
      fi
    done
    sleep "${LOCK_POLL_SECONDS:-5}"
  done
}

run_with_slot() {
  local name="$1"
  local limit="$2"
  shift 2

  require_positive_int "MAX_PARALLEL_${name}" "${limit}"

  local slot_dir
  slot_dir="$(acquire_slot "${name}" "${limit}")"

  (
    trap 'rmdir "${slot_dir}" 2>/dev/null || true' EXIT
    "$@"
  )
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

load_reference() {
  local requested="$1"
  local -a fields
  mapfile -t fields < <(
    python "${SCRIPT_DIR}/reference_registry.py" resolve \
      --registry "${REFERENCE_REGISTRY}" "${requested}"
  )
  [ "${#fields[@]}" -eq 11 ] || die "Could not resolve complete registry entry for ${requested}"

  RESOLVED_REFERENCE_ID="${fields[0]}"
  RESOLVED_SPECIES="${fields[1]}"
  RESOLVED_ASSEMBLY="${fields[2]}"
  RESOLVED_FASTA="${fields[3]}"
  RESOLVED_BWA_INDEX_PREFIX="${fields[4]}"
  RESOLVED_CHROM_SIZES="${fields[5]}"
  RESOLVED_CANONICAL_REGEX="${fields[6]}"
  RESOLVED_BROWSER_PRESET="${fields[7]}"
  RESOLVED_PHASING_TRACK="${fields[8]}"
  RESOLVED_ANNOTATION_METADATA="${fields[9]}"
  RESOLVED_BLACKLIST_METADATA="${fields[10]}"

  if [ "${RESOLVED_REFERENCE_ID}" = "${DEFAULT_REFERENCE_ID}" ] \
    && [ -n "${LEGACY_REFERENCE_FASTA}" ] && [ -n "${LEGACY_CHROM_SIZES}" ]; then
    RESOLVED_FASTA="${LEGACY_REFERENCE_FASTA}"
    RESOLVED_BWA_INDEX_PREFIX="${LEGACY_BWA_INDEX_PREFIX:-${LEGACY_REFERENCE_FASTA}}"
    RESOLVED_CHROM_SIZES="${LEGACY_CHROM_SIZES}"
    [ -z "${LEGACY_CANONICAL_REGEX}" ] || RESOLVED_CANONICAL_REGEX="${LEGACY_CANONICAL_REGEX}"
    [ -z "${LEGACY_PHASING_TRACK}" ] || RESOLVED_PHASING_TRACK="${LEGACY_PHASING_TRACK}"
  fi
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
# shellcheck source=reference_sentinel.sh
source "${SCRIPT_DIR}/reference_sentinel.sh"
python "${SCRIPT_DIR}/sanitize_text_inputs.py" --kind config "${CONFIG}"
python "${SCRIPT_DIR}/sanitize_text_inputs.py" --kind samplesheet "${SAMPLESHEET}"

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
LEGACY_REFERENCE_FASTA="${REFERENCE_FASTA:-}"
LEGACY_BWA_INDEX_PREFIX="${BWA_INDEX_PREFIX:-}"
LEGACY_CHROM_SIZES="${CHROM_SIZES:-}"
LEGACY_CANONICAL_REGEX="${CANONICAL_CHROMS_REGEX:-}"
LEGACY_PHASING_TRACK="${PHASING_TRACK:-}"
DEFAULT_REFERENCE_ID="${DEFAULT_REFERENCE_ID:-${GENOME_ASSEMBLY:-mm39}}"
mapfile -t default_reference_fields < <(
  python "${SCRIPT_DIR}/reference_registry.py" resolve \
    --registry "${REFERENCE_REGISTRY}" "${DEFAULT_REFERENCE_ID}"
)
[ "${#default_reference_fields[@]}" -eq 11 ] || die "Could not resolve DEFAULT_REFERENCE_ID=${DEFAULT_REFERENCE_ID}"
DEFAULT_REFERENCE_ID="${default_reference_fields[0]}"
THREADS_FASTP="${THREADS_FASTP:-${THREADS_ALIGN:-1}}"
THREADS_HIC_NORM="${THREADS_HIC_NORM:-24}"
MCOOL_RESOLUTIONS="${MCOOL_RESOLUTIONS:-${RESOLUTIONS}}"
HIC_RESOLUTIONS="${HIC_RESOLUTIONS:-${RESOLUTIONS}}"
HIC_NORMALIZATIONS="${HIC_NORMALIZATIONS:-VC,VC_SQRT,KR,SCALE}"
BOOTSTRAP_SENTINELS="${BOOTSTRAP_SENTINELS:-true}"
ALLOW_LEGACY_MM39_RESUME="${ALLOW_LEGACY_MM39_RESUME:-true}"

MAX_PARALLEL_SAMPLES="${MAX_PARALLEL_SAMPLES:-1}"
MAX_PARALLEL_MATRIX="${MAX_PARALLEL_MATRIX:-1}"
MAX_PARALLEL_HIC="${MAX_PARALLEL_HIC:-1}"
LOCK_POLL_SECONDS="${LOCK_POLL_SECONDS:-5}"

require_positive_int "MAX_PARALLEL_SAMPLES" "${MAX_PARALLEL_SAMPLES}"
require_positive_int "MAX_PARALLEL_MATRIX" "${MAX_PARALLEL_MATRIX}"
require_positive_int "MAX_PARALLEL_HIC" "${MAX_PARALLEL_HIC}"
require_positive_int "LOCK_POLL_SECONDS" "${LOCK_POLL_SECONDS}"
require_positive_int "THREADS_FASTP" "${THREADS_FASTP}"
require_positive_int "THREADS_ALIGN" "${THREADS_ALIGN}"
require_positive_int "THREADS_SORT" "${THREADS_SORT}"
require_positive_int "THREADS_MATRIX" "${THREADS_MATRIX}"
require_positive_int "THREADS_HIC_NORM" "${THREADS_HIC_NORM}"

FASTP_TIMEOUT_SECONDS="${FASTP_TIMEOUT_SECONDS:-0}"
is_nonnegative_int "${FASTP_TIMEOUT_SECONDS}" || die "FASTP_TIMEOUT_SECONDS must be a non-negative integer, got '${FASTP_TIMEOUT_SECONDS}'"

if [ "${FASTP_TIMEOUT_SECONDS}" -gt 0 ] && ! command -v timeout >/dev/null 2>&1; then
  die "FASTP_TIMEOUT_SECONDS is set but GNU/coreutils timeout was not found"
fi

PIPELINE_RUN_ID="$(date '+%Y%m%d_%H%M%S')_$$"
LOCK_ROOT="${TMPDIR}/microc2tracks_locks/${PIPELINE_RUN_ID}"
mkdir -p "${LOCK_ROOT}"
NORMALIZED_SAMPLESHEET="${LOCK_ROOT}/resolved_samplesheet.tsv"

python "${SCRIPT_DIR}/validate_samplesheet.py" \
  --require-files \
  --reference-registry "${REFERENCE_REGISTRY}" \
  --default-reference "${DEFAULT_REFERENCE_ID}" \
  --normalized-output "${NORMALIZED_SAMPLESHEET}" \
  --output-format tsv \
  "${SAMPLESHEET}"

build_select_expr() {
  local assay="$1"
  local min_cis_dist="0"
  local expr

  expr='((pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU"))'
  expr="${expr} and (mapq1 >= ${MAPQ_THRESHOLD}) and (mapq2 >= ${MAPQ_THRESHOLD})"

  case "${assay}" in
    microc) min_cis_dist="${MICROC_MIN_CIS_DIST:-1000}" ;;
    hic) min_cis_dist="${HIC_MIN_CIS_DIST:-0}" ;;
    *) die "Unsupported assay: ${assay}" ;;
  esac

  if [ "${min_cis_dist}" != "0" ]; then
    expr="${expr} and ((chrom1 != chrom2) or (abs(pos2 - pos1) >= ${min_cis_dist}))"
  fi

  printf '%s\n' "${expr}"
}

merge_technical_replicate_groups() {
  local manifest="$1"

  if ! bool_true "${RUN_MERGE_TECHNICAL_REPLICATES:-true}"; then
    log_msg "Technical replicate merging disabled"
    return 0
  fi

  [ -s "${manifest}" ] || return 0

  tail -n +2 "${manifest}" | cut -f1 | sort -u | while read -r group_id; do
    [ -n "${group_id}" ] || continue

    local count
    count="$(awk -F'\t' -v group="${group_id}" 'NR > 1 && $1 == group {count++} END{print count+0}' "${manifest}")"

    if [ "${count}" -le 1 ]; then
      log_msg "${group_id}: only one technical replicate, skipping merge"
      continue
    fi

    local pairs_csv
    pairs_csv="$(
      awk -F'\t' -v group="${group_id}" 'NR > 1 && $1 == group {print $3}' "${manifest}" \
        | paste -sd, -
    )"

    local reference_id
    reference_id="$(awk -F'\t' -v group="${group_id}" 'NR > 1 && $1 == group {print $4}' "${manifest}" | sort -u)"
    [ "$(printf '%s\n' "${reference_id}" | sed '/^$/d' | wc -l)" -eq 1 ] \
      || die "${group_id}: prohibited merge across multiple references: ${reference_id}"

    log_msg "${group_id}: merging ${count} technical replicates"
    bash "${SCRIPT_DIR}/merge_replicates.sh" \
      -c "${CONFIG}" \
      -n "${group_id}" \
      -r "${reference_id}" \
      -p "${pairs_csv}"
  done
}

build_cool_mcool_products() {
  local sample="$1"
  local pairs="$2"
  local log_dir="$3"
  local matrix_chrom_sizes="$4"
  local raw_cool="$5"
  local norm_cool="$6"
  local raw_mcool="$7"
  local norm_mcool="$8"
  local done_dir="$9"
  local matrix_dir
  local raw_cool_done="${done_dir}/raw_cool.done"
  local norm_cool_done="${done_dir}/norm_cool.done"
  local raw_mcool_done="${done_dir}/raw_mcool.done"
  local norm_mcool_done="${done_dir}/norm_mcool.done"
  local tmp

  matrix_dir="$(dirname "${raw_cool}")"

  if bool_true "${RUN_MCOOL:-true}" \
    && ! bool_true "${KEEP_SINGLE_RES_COOL:-false}" \
    && ! bool_true "${KEEP_RAW_MCOOL:-false}" \
    && step_done "${norm_mcool_done}" "${sample}:norm_mcool" "${norm_mcool}"; then
    log_msg "${sample}: normalized mcool exists, skipping cooler/mcool rebuild"
    return 0
  fi

  if step_done "${raw_cool_done}" "${sample}:raw_cool" "${raw_cool}"; then
    log_msg "${sample}: raw cooler exists, skipping"
  else
    log_msg "${sample}: building raw cooler"
    tmp="${raw_cool}.tmp.$$"
    rm -f "${tmp}"
    cooler cload pairix --assembly "${GENOME_ASSEMBLY}" \
      -p "${THREADS_MATRIX}" \
      "${matrix_chrom_sizes}:${BASE_RESOLUTION}" \
      "${pairs}" \
      "${tmp}" \
      > "${log_dir}/${sample}.cooler_cload.log" 2>&1
    mv -f "${tmp}" "${raw_cool}"
    mark_done "${raw_cool_done}" "${sample}:raw_cool" "${raw_cool}"
  fi

  if step_done "${norm_cool_done}" "${sample}:norm_cool" "${norm_cool}"; then
    log_msg "${sample}: normalized cooler exists, skipping"
  else
    log_msg "${sample}: balancing single-resolution cooler"
    tmp="${norm_cool}.tmp.$$"
    rm -f "${tmp}"
    cp "${raw_cool}" "${tmp}"
    cooler balance -p "${THREADS_MATRIX}" -f "${tmp}" \
      > "${log_dir}/${sample}.cooler_balance.log" 2>&1
    mv -f "${tmp}" "${norm_cool}"
    mark_done "${norm_cool_done}" "${sample}:norm_cool" "${norm_cool}"
  fi

  if bool_true "${RUN_MCOOL:-true}"; then
    if step_done "${raw_mcool_done}" "${sample}:raw_mcool" "${raw_mcool}"; then
      log_msg "${sample}: raw mcool exists, skipping"
    else
      log_msg "${sample}: building raw mcool"
      tmp="${raw_mcool}.tmp.$$"
      rm -f "${tmp}"
      cooler zoomify -p "${THREADS_MATRIX}" \
        -r "${MCOOL_RESOLUTIONS}" \
        -o "${tmp}" \
        "${raw_cool}" \
        > "${log_dir}/${sample}.cooler_zoomify_raw.log" 2>&1
      mv -f "${tmp}" "${raw_mcool}"
      mark_done "${raw_mcool_done}" "${sample}:raw_mcool" "${raw_mcool}"
    fi

    if step_done "${norm_mcool_done}" "${sample}:norm_mcool" "${norm_mcool}"; then
      log_msg "${sample}: normalized mcool exists, skipping"
    else
      log_msg "${sample}: building balanced mcool"
      tmp="${norm_mcool}.tmp.$$"
      rm -f "${tmp}"
      cooler zoomify -p "${THREADS_MATRIX}" \
        -r "${MCOOL_RESOLUTIONS}" \
        --balance \
        -o "${tmp}" \
        "${raw_cool}" \
        > "${log_dir}/${sample}.cooler_zoomify_norm.log" 2>&1
      mv -f "${tmp}" "${norm_mcool}"
      mark_done "${norm_mcool_done}" "${sample}:norm_mcool" "${norm_mcool}"
    fi
  fi

  if [ -s "${raw_cool}" ]; then
    cooler info "${raw_cool}" > "${matrix_dir}/${sample}.raw.${BASE_RESOLUTION}.cool.info.json"
  fi
}

build_hic_products() {
  local sample="$1"
  local pairs="$2"
  local matrix_dir="$3"
  local pairs_dir="$4"
  local log_dir="$5"
  local matrix_chrom_sizes="$6"
  local juicer_pairs="$7"
  local raw_hic="$8"
  local norm_hic="$9"
  local done_dir="${10}"
  local juicer_pairs_done="${done_dir}/juicer_pairs.done"
  local raw_hic_done="${done_dir}/raw_hic.done"
  local norm_hic_done="${done_dir}/norm_hic.done"
  local tmp

  if bool_true "${RUN_HIC:-true}"; then
    [ -f "${JUICER_TOOLS_JAR}" ] || die "JUICER_TOOLS_JAR not found: ${JUICER_TOOLS_JAR}"

    if step_done "${norm_hic_done}" "${sample}:norm_hic" "${norm_hic}"; then
      log_msg "${sample}: normalized hic exists, skipping hic generation"
      write_ucsc_hic_track "${sample}" "${norm_hic}" "${matrix_dir}"
      return 0
    fi

    if step_done "${juicer_pairs_done}" "${sample}:juicer_pairs" "${juicer_pairs}"; then
      log_msg "${sample}: Juicer-compatible pairs exist, skipping"
    else
      log_msg "${sample}: writing Juicer-compatible pairs"
      tmp="${juicer_pairs%.pairs.gz}.tmp.$$.pairs.gz"
      rm -f "${tmp}"
      zcat "${pairs}" \
        | awk 'BEGIN{OFS="\t"} /^## pairs format/ {print; next} /^#columns:/ {print; next} /^#/ {next} {print}' \
        | bgzip -@ "${THREADS_MATRIX}" \
        > "${tmp}"
      mv -f "${tmp}" "${juicer_pairs}"
      mark_done "${juicer_pairs_done}" "${sample}:juicer_pairs" "${juicer_pairs}"
    fi

    if step_done "${raw_hic_done}" "${sample}:raw_hic" "${raw_hic}"; then
      log_msg "${sample}: raw hic exists, skipping"
    else
      log_msg "${sample}: building raw hic"
      tmp="${raw_hic}.tmp.$$"
      rm -f "${tmp}"
      java -Xmx"${JAVA_HEAP}" -jar "${JUICER_TOOLS_JAR}" pre \
        -n \
        -r "${HIC_RESOLUTIONS}" \
        "${juicer_pairs}" \
        "${tmp}" \
        "${matrix_chrom_sizes}" \
        > "${log_dir}/${sample}.juicer_pre.log" 2>&1
      mv -f "${tmp}" "${raw_hic}"
      mark_done "${raw_hic_done}" "${sample}:raw_hic" "${raw_hic}"
    fi

    log_msg "${sample}: adding Juicer normalizations"
    tmp="${norm_hic}.tmp.$$"
    rm -f "${tmp}"
    cp "${raw_hic}" "${tmp}"
    java -Xmx"${JAVA_HEAP}" -jar "${JUICER_TOOLS_JAR}" addNorm \
      -j "${THREADS_HIC_NORM}" \
      -k "${HIC_NORMALIZATIONS}" \
      "${tmp}" \
      > "${log_dir}/${sample}.juicer_addNorm.log" 2>&1
    mv -f "${tmp}" "${norm_hic}"
    mark_done "${norm_hic_done}" "${sample}:norm_hic" "${norm_hic}"

    write_ucsc_hic_track "${sample}" "${norm_hic}" "${matrix_dir}"
  fi
}

build_matrices_from_pairs() {
  local sample="$1"
  local pairs="$2"
  local sample_dir="$3"
  local log_dir="$4"
  local matrix_chrom_sizes="$5"

  local matrix_dir="${sample_dir}/04_matrices"
  local pairs_dir="${sample_dir}/03_pairs"
  local raw_cool="${matrix_dir}/${sample}.raw.${BASE_RESOLUTION}.cool"
  local norm_cool="${matrix_dir}/${sample}.norm.${BASE_RESOLUTION}.cool"
  local raw_mcool="${matrix_dir}/${sample}.raw.mcool"
  local norm_mcool="${matrix_dir}/${sample}.norm.mcool"
  local juicer_pairs="${pairs_dir}/${sample}.valid.mapq${MAPQ_THRESHOLD}.juicer.pairs.gz"
  local raw_hic="${matrix_dir}/${sample}.raw.hic"
  local norm_hic="${matrix_dir}/${sample}.norm.hic"
  local done_dir="${log_dir}/done"

  mkdir -p "${matrix_dir}"

  run_with_slot matrix "${MAX_PARALLEL_MATRIX}" \
    build_cool_mcool_products \
      "${sample}" "${pairs}" "${log_dir}" "${matrix_chrom_sizes}" \
      "${raw_cool}" "${norm_cool}" "${raw_mcool}" "${norm_mcool}" "${done_dir}"

  run_with_slot hic "${MAX_PARALLEL_HIC}" \
    build_hic_products \
      "${sample}" "${pairs}" "${matrix_dir}" "${pairs_dir}" "${log_dir}" \
      "${matrix_chrom_sizes}" "${juicer_pairs}" "${raw_hic}" "${norm_hic}" "${done_dir}"
}

cleanup_sample_intermediates() {
  local sample="$1"
  local trim_r1="$2"
  local trim_r2="$3"
  local dedup_pairs="$4"
  local juicer_pairs="$5"
  local raw_cool="$6"
  local norm_cool="$7"
  local raw_mcool="$8"
  local raw_hic="$9"
  local sample_tmp="${10}"

  if bool_true "${RUN_FASTP:-true}" && ! bool_true "${KEEP_TRIMMED_FASTQ:-false}"; then
    log_msg "${sample}: removing trimmed FASTQ intermediates"
    rm -f "${trim_r1}" "${trim_r2}"
  fi

  if ! bool_true "${KEEP_DEDUP_PAIRS:-true}"; then
    log_msg "${sample}: removing deduplicated pair intermediates"
    rm -f "${dedup_pairs}" "${dedup_pairs}.px2"
  fi

  if ! bool_true "${KEEP_JUICER_PAIRS:-false}"; then
    rm -f "${juicer_pairs}"
  fi

  if ! bool_true "${KEEP_SINGLE_RES_COOL:-false}"; then
    rm -f "${raw_cool}" "${norm_cool}"
  fi

  if ! bool_true "${KEEP_RAW_MCOOL:-false}"; then
    rm -f "${raw_mcool}"
  fi

  if ! bool_true "${KEEP_RAW_HIC:-false}"; then
    rm -f "${raw_hic}"
  fi

  if bool_true "${CLEAN_TMP_ON_SUCCESS:-true}"; then
    rm -rf "${sample_tmp}"
  fi
}

generate_final_report() {
  local report_dir="${OUTDIR}/final_report"

  mkdir -p "${report_dir}"

  if bool_true "${RUN_GLOBAL_MULTIQC:-true}"; then
    log_msg "Generating global MultiQC report"
    multiqc "${OUTDIR}" -n "microc2tracks_multiqc.html" -o "${report_dir}" \
      > "${report_dir}/microc2tracks_multiqc.log" 2>&1 || true
  fi

  if bool_true "${RUN_FINAL_REPORT:-true}"; then
    log_msg "Generating microc2tracks final summary report"
    python "${SCRIPT_DIR}/summarize_run.py" \
      -c "${CONFIG}" \
      -s "${SAMPLESHEET}" \
      -o "${report_dir}" \
      --sample-manifest "${SAMPLE_MANIFEST}" \
      --technical-manifest "${TECH_MANIFEST}" \
      --merge-manifest "${MERGE_MANIFEST}"
  fi
}

run_prelim_downstream() {
  local sample="$1"
  local matrix="$2"
  local log_dir="$3"
  local done_dir="$4"
  local done_file="${done_dir}/prelim_downstream.done"

  if ! bool_true "${RUN_PRELIM_DOWNSTREAM:-true}"; then
    return 0
  fi

  if step_done "${done_file}" "${sample}:prelim_downstream" "${matrix}"; then
    log_msg "${sample}: preliminary downstream already marked complete, skipping"
    return 0
  fi

  if [ ! -s "${matrix}" ]; then
    log_msg "${sample}: normalized mcool not found; skipping preliminary downstream"
    return 0
  fi

  log_msg "${sample}: running preliminary downstream analyses"
  MICROC2TRACKS_PHASING_TRACK_OVERRIDE="${PHASING_TRACK}" \
  MICROC2TRACKS_REFERENCE_ID="${REFERENCE_ID}" \
  bash "${SCRIPT_DIR}/run_downstream.sh" \
    -l \
    -c "${CONFIG}" \
    -s "${sample}" \
    -m "${matrix}" \
    > "${log_dir}/${sample}.prelim_downstream.log" 2>&1 || true
  mark_done "${done_file}" "${sample}:prelim_downstream" "${matrix}"
}

fastq_file_summary() {
  local path="$1"

  if stat -c '%s bytes, mtime=%y' "${path}" 2>/dev/null; then
    return 0
  fi

  ls -lh "${path}" 2>/dev/null || printf 'metadata unavailable'
}

run_fastp_step() {
  local sample="$1"
  local fastq_r1="$2"
  local fastq_r2="$3"
  local trim_r1="$4"
  local trim_r2="$5"
  local qc_dir="$6"
  local log_dir="$7"
  local fastp_log="${log_dir}/${sample}.fastp.log"
  local tmp_trim_r1="${trim_r1%.fastq.gz}.tmp.$$.fastq.gz"
  local tmp_trim_r2="${trim_r2%.fastq.gz}.tmp.$$.fastq.gz"
  local fastp_html="${qc_dir}/${sample}.fastp.html"
  local fastp_json="${qc_dir}/${sample}.fastp.json"
  local tmp_html="${fastp_html}.tmp.$$"
  local tmp_json="${fastp_json}.tmp.$$"
  local status

  log_msg "${sample}: fastq_r1 ${fastq_r1} ($(fastq_file_summary "${fastq_r1}"))"
  log_msg "${sample}: fastq_r2 ${fastq_r2} ($(fastq_file_summary "${fastq_r2}"))"
  rm -f "${tmp_trim_r1}" "${tmp_trim_r2}" "${tmp_html}" "${tmp_json}"

  set +e
  if [ "${FASTP_TIMEOUT_SECONDS}" -gt 0 ]; then
    timeout "${FASTP_TIMEOUT_SECONDS}" fastp \
      -i "${fastq_r1}" \
      -I "${fastq_r2}" \
      -o "${tmp_trim_r1}" \
      -O "${tmp_trim_r2}" \
      --thread "${THREADS_FASTP}" \
      --html "${tmp_html}" \
      --json "${tmp_json}" \
      ${FASTP_EXTRA_ARGS:-} \
      > "${fastp_log}" 2>&1
  else
    fastp \
      -i "${fastq_r1}" \
      -I "${fastq_r2}" \
      -o "${tmp_trim_r1}" \
      -O "${tmp_trim_r2}" \
      --thread "${THREADS_FASTP}" \
      --html "${tmp_html}" \
      --json "${tmp_json}" \
      ${FASTP_EXTRA_ARGS:-} \
      > "${fastp_log}" 2>&1
  fi
  status="$?"
  set -e

  if [ "${status}" -ne 0 ]; then
    echo "ERROR: ${sample}: fastp failed with exit code ${status}" >&2
    if [ "${status}" -eq 124 ]; then
      echo "ERROR: ${sample}: fastp exceeded FASTP_TIMEOUT_SECONDS=${FASTP_TIMEOUT_SECONDS}" >&2
    fi
    echo "ERROR: fastp log: ${fastp_log}" >&2
    tail -n 40 "${fastp_log}" >&2 || true
    rm -f "${tmp_trim_r1}" "${tmp_trim_r2}" "${tmp_html}" "${tmp_json}"
    return "${status}"
  fi

  mv -f "${tmp_trim_r1}" "${trim_r1}"
  mv -f "${tmp_trim_r2}" "${trim_r2}"
  mv -f "${tmp_html}" "${fastp_html}"
  mv -f "${tmp_json}" "${fastp_json}"

  log_msg "${sample}: fastp finished; trimmed R1 $(fastq_file_summary "${trim_r1}")"
  log_msg "${sample}: fastp finished; trimmed R2 $(fastq_file_summary "${trim_r2}")"

}

process_sample() {
  local sample="$1"
  local assay="$2"
  local fastq_r1="$3"
  local fastq_r2="$4"
  local REFERENCE_ID="$5"
  local REFERENCE_SPECIES="$6"
  local GENOME_ASSEMBLY="$7"
  local REFERENCE_FASTA="$8"
  local BWA_INDEX_PREFIX="$9"
  local CHROM_SIZES="${10}"
  local CANONICAL_CHROMS_REGEX="${11}"
  local BROWSER_PRESET="${12}"
  local PHASING_TRACK="${13}"
  local ANNOTATION_METADATA="${14}"
  local BLACKLIST_METADATA="${15}"

  assay="$(printf '%s' "${assay}" | tr '[:upper:]' '[:lower:]')"

  local sample_dir="${OUTDIR}/${sample}"
  local qc_dir="${sample_dir}/01_qc"
  local trimmed_dir="${sample_dir}/02_trimmed"
  local pairs_dir="${sample_dir}/03_pairs"
  local log_dir="${sample_dir}/logs"
  local sample_tmp="${TMPDIR}/${sample}"
  local trim_r1="${trimmed_dir}/${sample}.trim.R1.fastq.gz"
  local trim_r2="${trimmed_dir}/${sample}.trim.R2.fastq.gz"
  local dedup_pairs="${pairs_dir}/${sample}.dedup.pairs.gz"
  local valid_pairs="${pairs_dir}/${sample}.valid.mapq${MAPQ_THRESHOLD}.pairs.gz"
  local matrix_chrom_sizes="${pairs_dir}/${sample}.matrix.chrom.sizes"
  local raw_cool="${sample_dir}/04_matrices/${sample}.raw.${BASE_RESOLUTION}.cool"
  local norm_cool="${sample_dir}/04_matrices/${sample}.norm.${BASE_RESOLUTION}.cool"
  local raw_mcool="${sample_dir}/04_matrices/${sample}.raw.mcool"
  local norm_mcool="${sample_dir}/04_matrices/${sample}.norm.mcool"
  local juicer_pairs="${pairs_dir}/${sample}.valid.mapq${MAPQ_THRESHOLD}.juicer.pairs.gz"
  local raw_hic="${sample_dir}/04_matrices/${sample}.raw.hic"
  local status_file="${log_dir}/${sample}.status.tsv"
  local done_dir="${log_dir}/done"
  local fastp_done="${done_dir}/fastp.done"
  local dedup_done="${done_dir}/dedup_pairs.done"
  local dedup_index_done="${done_dir}/dedup_pairix.done"
  local valid_done="${done_dir}/valid_pairs.done"
  local valid_index_done="${done_dir}/valid_pairix.done"
  local stats_done="${done_dir}/pairs_stats.done"
  local sample_start_epoch
  local sample_end_epoch
  local select_expr
  local tmp_pairs
  local tmp_stats

  mkdir -p "${qc_dir}" "${trimmed_dir}" "${pairs_dir}" "${log_dir}" "${done_dir}" "${sample_tmp}"
  [ -f "${REFERENCE_FASTA}" ] || die "${sample}: reference FASTA not found for ${REFERENCE_ID}: ${REFERENCE_FASTA}"
  [ -f "${CHROM_SIZES}" ] || die "${sample}: chromosome sizes not found for ${REFERENCE_ID}: ${CHROM_SIZES}"
  sample_start_epoch="$(date '+%s')"
  printf 'sample\treference_id\tassembly\tstatus\tstart_epoch\tend_epoch\truntime_seconds\n' > "${status_file}"
  printf '%s\t%s\t%s\trunning\t%s\t\t\n' "${sample}" "${REFERENCE_ID}" "${GENOME_ASSEMBLY}" "${sample_start_epoch}" >> "${status_file}"
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
  } > "${log_dir}/${sample}.reference.tsv"
  prepare_matrix_chrom_sizes "${matrix_chrom_sizes}"

  log_msg "${sample}: assay=${assay}, reference_id=${REFERENCE_ID}, assembly=${GENOME_ASSEMBLY}, browser_preset=${BROWSER_PRESET}"

  if ! step_done "${dedup_done}" "${sample}:dedup_pairs" "${dedup_pairs}"; then
    if bool_true "${RUN_FASTP:-true}"; then
      if step_done "${fastp_done}" "${sample}:fastp" "${trim_r1}" "${trim_r2}" "${qc_dir}/${sample}.fastp.json"; then
        log_msg "${sample}: trimmed FASTQ exists, skipping fastp"
      else
        log_msg "${sample}: running fastp"
        run_fastp_step \
          "${sample}" "${fastq_r1}" "${fastq_r2}" \
          "${trim_r1}" "${trim_r2}" "${qc_dir}" "${log_dir}" \
          || exit 1
        mark_done "${fastp_done}" "${sample}:fastp" "${trim_r1}" "${trim_r2}" "${qc_dir}/${sample}.fastp.json"
      fi
    else
      trim_r1="${fastq_r1}"
      trim_r2="${fastq_r2}"
    fi

    if bool_true "${CHECK_TRIMMED_FASTQ_SYNC:-false}"; then
      log_msg "${sample}: streaming alignment FASTQs to verify R1/R2 names and record counts"
      python "${SCRIPT_DIR}/check_fastq_pair_names.py" "${trim_r1}" "${trim_r2}" \
        > "${log_dir}/${sample}.fastq_sync_check.log"
    fi

    log_msg "${sample}: aligning, parsing, sorting, and deduplicating pairs"
    tmp_pairs="${dedup_pairs%.pairs.gz}.tmp.$$.pairs.gz"
    tmp_stats="${pairs_dir}/${sample}.dedup.stats.txt.tmp.$$"
    rm -f "${tmp_pairs}" "${tmp_stats}"
    bwa-mem2 mem ${BWA_EXTRA_ARGS:-} -t "${THREADS_ALIGN}" "${BWA_INDEX_PREFIX}" "${trim_r1}" "${trim_r2}" \
      2> "${log_dir}/${sample}.bwa_mem2.log" \
      | pairtools parse \
          -c "${CHROM_SIZES}" \
          --assembly "${GENOME_ASSEMBLY}" \
          --walks-policy 5unique \
          --drop-sam \
          --drop-seq \
          --add-columns mapq \
      | pairtools sort \
          --nproc "${THREADS_SORT}" \
          --memory "${PAIRTOOLS_SORT_MEMORY}" \
          --tmpdir "${sample_tmp}" \
      | pairtools dedup \
          --output-stats "${tmp_stats}" \
      | bgzip -@ "${THREADS_SORT}" \
      > "${tmp_pairs}"
    mv -f "${tmp_pairs}" "${dedup_pairs}"
    mv -f "${tmp_stats}" "${pairs_dir}/${sample}.dedup.stats.txt"
    mark_done "${dedup_done}" "${sample}:dedup_pairs" "${dedup_pairs}"
  else
    log_msg "${sample}: deduplicated pairs exist, skipping alignment"
  fi

  if step_done "${dedup_index_done}" "${sample}:dedup_pairix" "${dedup_pairs}.px2"; then
    :
  else
    pairix "${dedup_pairs}"
    mark_done "${dedup_index_done}" "${sample}:dedup_pairix" "${dedup_pairs}.px2"
  fi

  if step_done "${valid_done}" "${sample}:valid_pairs" "${valid_pairs}"; then
    log_msg "${sample}: valid pairs exist, skipping selection"
  else
    select_expr="$(build_select_expr "${assay}")"
    log_msg "${sample}: selecting valid pairs with expression: ${select_expr}"
    tmp_pairs="${valid_pairs%.pairs.gz}.tmp.$$.pairs.gz"
    rm -f "${tmp_pairs}"
    pairtools select "${select_expr}" "${dedup_pairs}" \
      | filter_pairs_to_chrom_sizes "${matrix_chrom_sizes}" \
      | bgzip -@ "${THREADS_SORT}" \
      > "${tmp_pairs}"
    mv -f "${tmp_pairs}" "${valid_pairs}"
    mark_done "${valid_done}" "${sample}:valid_pairs" "${valid_pairs}"
  fi

  if step_done "${valid_index_done}" "${sample}:valid_pairix" "${valid_pairs}.px2"; then
    :
  else
    pairix "${valid_pairs}"
    mark_done "${valid_index_done}" "${sample}:valid_pairix" "${valid_pairs}.px2"
  fi

  if step_done "${stats_done}" "${sample}:pairs_stats" "${pairs_dir}/${sample}.pairs.stats.txt"; then
    :
  else
    tmp_stats="${pairs_dir}/${sample}.pairs.stats.txt.tmp.$$"
    pairtools stats "${valid_pairs}" > "${tmp_stats}"
    mv -f "${tmp_stats}" "${pairs_dir}/${sample}.pairs.stats.txt"
    mark_done "${stats_done}" "${sample}:pairs_stats" "${pairs_dir}/${sample}.pairs.stats.txt"
  fi

  build_matrices_from_pairs "${sample}" "${valid_pairs}" "${sample_dir}" "${log_dir}" "${matrix_chrom_sizes}"
  run_prelim_downstream "${sample}" "${norm_mcool}" "${log_dir}" "${done_dir}"

  if bool_true "${RUN_MULTIQC:-true}"; then
    log_msg "${sample}: running MultiQC"
    multiqc "${sample_dir}" -n "${sample}.multiqc.html" -o "${qc_dir}" \
      > "${log_dir}/${sample}.multiqc.log" 2>&1 || true
  fi

  log_msg "${sample}: finished"
  sample_end_epoch="$(date '+%s')"
  printf 'sample\treference_id\tassembly\tstatus\tstart_epoch\tend_epoch\truntime_seconds\n' > "${status_file}"
  printf '%s\t%s\t%s\tsuccess\t%s\t%s\t%s\n' \
    "${sample}" "${REFERENCE_ID}" "${GENOME_ASSEMBLY}" "${sample_start_epoch}" "${sample_end_epoch}" "$((sample_end_epoch - sample_start_epoch))" \
    >> "${status_file}"

  cleanup_sample_intermediates \
    "${sample}" "${trim_r1}" "${trim_r2}" "${dedup_pairs}" "${juicer_pairs}" \
    "${raw_cool}" "${norm_cool}" "${raw_mcool}" "${raw_hic}" "${sample_tmp}"
}

mkdir -p "${OUTDIR}" "${TMPDIR}"
RUN_METADATA_DIR="${OUTDIR}/run_metadata"
SAMPLE_MANIFEST="${RUN_METADATA_DIR}/sample_manifest.tsv"
TECH_MANIFEST="${RUN_METADATA_DIR}/technical_replicates.tsv"
MERGE_MANIFEST="${RUN_METADATA_DIR}/merge_manifest.tsv"
mkdir -p "${RUN_METADATA_DIR}"

printf 'sample\tassay\treference_id\tspecies\tassembly\tbrowser_preset\tcondition\tbiological_replicate\ttechnical_replicate\tmerge_group\tfastq_r1\tfastq_r2\tsample_dir\n' > "${SAMPLE_MANIFEST}"
printf 'merge_group\tsample\tpairs\treference_id\tassay\tcondition\tbiological_replicate\ttechnical_replicate\n' > "${TECH_MANIFEST}"
printf 'merge_group\treference_id\tassembly\tbrowser_preset\n' > "${MERGE_MANIFEST}"

active_jobs=0
sample_failures=0
sample_count=0

while IFS=$'\t' read -r sample assay reference_id condition biological_replicate technical_replicate fastq_r1 fastq_r2 requested_merge_group; do

  [ -n "${sample// }" ] || continue

  assay_group="$(printf '%s' "${assay}" | tr '[:upper:]' '[:lower:]')"
  if [ -n "${requested_merge_group}" ]; then
    group_id="${requested_merge_group}"
  else
    group_id="$(sanitize_id "${condition}_${assay_group}_B${biological_replicate}_tech_merged")"
  fi
  valid_pairs="${OUTDIR}/${sample}/03_pairs/${sample}.valid.mapq${MAPQ_THRESHOLD}.pairs.gz"

  load_reference "${reference_id}"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${sample}" "${assay_group}" "${RESOLVED_REFERENCE_ID}" "${RESOLVED_SPECIES}" "${RESOLVED_ASSEMBLY}" "${RESOLVED_BROWSER_PRESET}" \
    "${condition}" "${biological_replicate}" "${technical_replicate}" "${group_id}" "${fastq_r1}" "${fastq_r2}" "${OUTDIR}/${sample}" \
    >> "${SAMPLE_MANIFEST}"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${group_id}" "${sample}" "${valid_pairs}" "${RESOLVED_REFERENCE_ID}" "${assay_group}" "${condition}" "${biological_replicate}" "${technical_replicate}" \
    >> "${TECH_MANIFEST}"
  if ! awk -F'\t' -v group="${group_id}" 'NR > 1 && $1 == group {found=1} END{exit !found}' "${MERGE_MANIFEST}"; then
    printf '%s\t%s\t%s\t%s\n' "${group_id}" "${RESOLVED_REFERENCE_ID}" "${RESOLVED_ASSEMBLY}" "${RESOLVED_BROWSER_PRESET}" >> "${MERGE_MANIFEST}"
  fi

  log_msg "${sample}: queueing sample worker for ${RESOLVED_REFERENCE_ID}"
  process_sample \
    "${sample}" "${assay_group}" "${fastq_r1}" "${fastq_r2}" \
    "${RESOLVED_REFERENCE_ID}" "${RESOLVED_SPECIES}" "${RESOLVED_ASSEMBLY}" \
    "${RESOLVED_FASTA}" "${RESOLVED_BWA_INDEX_PREFIX}" "${RESOLVED_CHROM_SIZES}" \
    "${RESOLVED_CANONICAL_REGEX}" "${RESOLVED_BROWSER_PRESET}" "${RESOLVED_PHASING_TRACK}" \
    "${RESOLVED_ANNOTATION_METADATA}" "${RESOLVED_BLACKLIST_METADATA}" &
  active_jobs=$((active_jobs + 1))
  sample_count=$((sample_count + 1))

  if [ "${active_jobs}" -ge "${MAX_PARALLEL_SAMPLES}" ]; then
    if ! wait -n; then
      sample_failures=1
    fi
    active_jobs=$((active_jobs - 1))
  fi
done < <(tail -n +2 "${NORMALIZED_SAMPLESHEET}")

while [ "${active_jobs}" -gt 0 ]; do
  if ! wait -n; then
    sample_failures=1
  fi
  active_jobs=$((active_jobs - 1))
done

[ "${sample_count}" -gt 0 ] || die "No samples found in ${SAMPLESHEET}"
[ "${sample_failures}" -eq 0 ] || die "One or more sample workers failed; check per-sample logs under ${OUTDIR}/*/logs"

merge_technical_replicate_groups "${TECH_MANIFEST}"
generate_final_report

if bool_true "${CLEAN_TMP_ON_SUCCESS:-true}"; then
  rm -rf "${LOCK_ROOT}"
fi

log_msg "All samples finished"
