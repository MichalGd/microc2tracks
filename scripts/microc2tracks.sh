#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  microc2tracks.sh -c config/config.conf -s config/samplesheet.csv

Runs paired-end FASTQ through:
  fastp -> bwa-mem2 -> pairtools -> cooler/mcool -> Juicer .hic

Sample sheet columns:
  sample,assay,replicate_group,condition,fastq_r1,fastq_r2
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

bool_true() {
  [ "${1:-false}" = "true" ] || [ "${1:-false}" = "1" ] || [ "${1:-false}" = "yes" ]
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

[ -f "${REFERENCE_FASTA}" ] || die "REFERENCE_FASTA not found: ${REFERENCE_FASTA}"
[ -f "${CHROM_SIZES}" ] || die "CHROM_SIZES not found: ${CHROM_SIZES}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python "${SCRIPT_DIR}/validate_samplesheet.py" "${SAMPLESHEET}"

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

build_matrices_from_pairs() {
  local sample="$1"
  local pairs="$2"
  local sample_dir="$3"
  local log_dir="$4"

  local matrix_dir="${sample_dir}/04_matrices"
  local raw_cool="${matrix_dir}/${sample}.raw.${BASE_RESOLUTION}.cool"
  local norm_cool="${matrix_dir}/${sample}.norm.${BASE_RESOLUTION}.cool"
  local raw_mcool="${matrix_dir}/${sample}.raw.mcool"
  local norm_mcool="${matrix_dir}/${sample}.norm.mcool"
  local juicer_pairs="${sample_dir}/03_pairs/${sample}.valid.mapq${MAPQ_THRESHOLD}.juicer.pairs.gz"
  local raw_hic="${matrix_dir}/${sample}.raw.hic"
  local norm_hic="${matrix_dir}/${sample}.norm.hic"

  mkdir -p "${matrix_dir}"

  if [ ! -s "${raw_cool}" ]; then
    log_msg "${sample}: building raw cooler"
    cooler cload pairix --assembly "${GENOME_ASSEMBLY}" \
      -p "${THREADS_MATRIX}" \
      "${CHROM_SIZES}:${BASE_RESOLUTION}" \
      "${pairs}" \
      "${raw_cool}" \
      > "${log_dir}/${sample}.cooler_cload.log" 2>&1
  else
    log_msg "${sample}: raw cooler exists, skipping"
  fi

  if [ ! -s "${norm_cool}" ]; then
    log_msg "${sample}: balancing single-resolution cooler"
    cp "${raw_cool}" "${norm_cool}"
    cooler balance -p "${THREADS_MATRIX}" -f "${norm_cool}" \
      > "${log_dir}/${sample}.cooler_balance.log" 2>&1
  else
    log_msg "${sample}: normalized cooler exists, skipping"
  fi

  if bool_true "${RUN_MCOOL:-true}"; then
    if [ ! -s "${raw_mcool}" ]; then
      log_msg "${sample}: building raw mcool"
      cooler zoomify -p "${THREADS_MATRIX}" \
        -r "${RESOLUTIONS}" \
        -o "${raw_mcool}" \
        "${raw_cool}" \
        > "${log_dir}/${sample}.cooler_zoomify_raw.log" 2>&1
    else
      log_msg "${sample}: raw mcool exists, skipping"
    fi

    if [ ! -s "${norm_mcool}" ]; then
      log_msg "${sample}: building balanced mcool"
      cooler zoomify -p "${THREADS_MATRIX}" \
        -r "${RESOLUTIONS}" \
        --balance \
        -o "${norm_mcool}" \
        "${raw_cool}" \
        > "${log_dir}/${sample}.cooler_zoomify_norm.log" 2>&1
    else
      log_msg "${sample}: normalized mcool exists, skipping"
    fi
  fi

  cooler info "${raw_cool}" > "${matrix_dir}/${sample}.raw.${BASE_RESOLUTION}.cool.info.json"

  if bool_true "${RUN_HIC:-true}"; then
    [ -f "${JUICER_TOOLS_JAR}" ] || die "JUICER_TOOLS_JAR not found: ${JUICER_TOOLS_JAR}"

    if [ ! -s "${juicer_pairs}" ]; then
      log_msg "${sample}: writing Juicer-compatible pairs"
      zcat "${pairs}" \
        | awk 'BEGIN{OFS="\t"} /^## pairs format/ {print; next} /^#columns:/ {print; next} /^#/ {next} {print}' \
        | bgzip -@ "${THREADS_MATRIX}" \
        > "${juicer_pairs}"
    fi

    if [ ! -s "${raw_hic}" ]; then
      log_msg "${sample}: building raw hic"
      java -Xmx"${JAVA_HEAP}" -jar "${JUICER_TOOLS_JAR}" pre \
        -n \
        -r "${RESOLUTIONS}" \
        "${juicer_pairs}" \
        "${raw_hic}" \
        "${CHROM_SIZES}" \
        > "${log_dir}/${sample}.juicer_pre.log" 2>&1
    else
      log_msg "${sample}: raw hic exists, skipping"
    fi

    if [ ! -s "${norm_hic}" ]; then
      log_msg "${sample}: adding Juicer normalizations"
      cp "${raw_hic}" "${norm_hic}"
      java -Xmx"${JAVA_HEAP}" -jar "${JUICER_TOOLS_JAR}" addNorm \
        -k VC,VC_SQRT,KR,SCALE \
        "${norm_hic}" \
        > "${log_dir}/${sample}.juicer_addNorm.log" 2>&1
    else
      log_msg "${sample}: normalized hic exists, skipping"
    fi
  fi
}

process_sample() {
  local sample="$1"
  local assay="$2"
  local fastq_r1="$3"
  local fastq_r2="$4"

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
  local select_expr

  mkdir -p "${qc_dir}" "${trimmed_dir}" "${pairs_dir}" "${log_dir}" "${sample_tmp}"

  log_msg "${sample}: assay=${assay}"

  if bool_true "${RUN_FASTP:-true}"; then
    if [ ! -s "${trim_r1}" ] || [ ! -s "${trim_r2}" ]; then
      log_msg "${sample}: running fastp"
      fastp \
        -i "${fastq_r1}" \
        -I "${fastq_r2}" \
        -o "${trim_r1}" \
        -O "${trim_r2}" \
        --thread "${THREADS_ALIGN}" \
        --html "${qc_dir}/${sample}.fastp.html" \
        --json "${qc_dir}/${sample}.fastp.json" \
        ${FASTP_EXTRA_ARGS:-} \
        > "${log_dir}/${sample}.fastp.log" 2>&1
    else
      log_msg "${sample}: trimmed FASTQ exists, skipping fastp"
    fi
  else
    trim_r1="${fastq_r1}"
    trim_r2="${fastq_r2}"
  fi

  if [ ! -s "${dedup_pairs}" ]; then
    log_msg "${sample}: aligning, parsing, sorting, and deduplicating pairs"
    bwa-mem2 mem ${BWA_EXTRA_ARGS:-} -t "${THREADS_ALIGN}" "${REFERENCE_FASTA}" "${trim_r1}" "${trim_r2}" \
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
          --output-stats "${pairs_dir}/${sample}.dedup.stats.txt" \
      | bgzip -@ "${THREADS_SORT}" \
      > "${dedup_pairs}"
  else
    log_msg "${sample}: deduplicated pairs exist, skipping alignment"
  fi

  if [ ! -s "${dedup_pairs}.px2" ]; then
    pairix "${dedup_pairs}"
  fi

  if [ ! -s "${valid_pairs}" ]; then
    select_expr="$(build_select_expr "${assay}")"
    log_msg "${sample}: selecting valid pairs with expression: ${select_expr}"
    pairtools select "${select_expr}" "${dedup_pairs}" \
      | bgzip -@ "${THREADS_SORT}" \
      > "${valid_pairs}"
  else
    log_msg "${sample}: valid pairs exist, skipping selection"
  fi

  if [ ! -s "${valid_pairs}.px2" ]; then
    pairix "${valid_pairs}"
  fi

  pairtools stats "${valid_pairs}" > "${pairs_dir}/${sample}.pairs.stats.txt"

  build_matrices_from_pairs "${sample}" "${valid_pairs}" "${sample_dir}" "${log_dir}"

  if bool_true "${RUN_MULTIQC:-true}"; then
    log_msg "${sample}: running MultiQC"
    multiqc "${sample_dir}" -n "${sample}.multiqc.html" -o "${qc_dir}" \
      > "${log_dir}/${sample}.multiqc.log" 2>&1 || true
  fi

  log_msg "${sample}: finished"
}

mkdir -p "${OUTDIR}" "${TMPDIR}"

tail -n +2 "${SAMPLESHEET}" | while IFS=, read -r sample assay replicate_group condition fastq_r1 fastq_r2 rest; do
  [ -n "${sample// }" ] || continue
  process_sample "${sample}" "${assay}" "${fastq_r1}" "${fastq_r2}"
done

log_msg "All samples finished"
