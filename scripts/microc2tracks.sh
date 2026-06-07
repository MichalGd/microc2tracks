#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat <<'USAGE'
Usage:
  microc2tracks.sh -c config/config.conf -s config/samplesheet.csv

Runs paired-end FASTQ through:
  fastp -> bwa-mem2 -> pairtools -> cooler/mcool -> Juicer .hic

Sample sheet columns:
  sample,assay,condition,biological_replicate,technical_replicate,fastq_r1,fastq_r2
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

[ -f "${REFERENCE_FASTA}" ] || die "REFERENCE_FASTA not found: ${REFERENCE_FASTA}"
[ -f "${CHROM_SIZES}" ] || die "CHROM_SIZES not found: ${CHROM_SIZES}"
BWA_INDEX_PREFIX="${BWA_INDEX_PREFIX:-$REFERENCE_FASTA}"

python "${SCRIPT_DIR}/validate_samplesheet.py" --require-files "${SAMPLESHEET}"

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

  cut -f1 "${manifest}" | sort -u | while read -r group_id; do
    [ -n "${group_id}" ] || continue

    local count
    count="$(awk -F'\t' -v group="${group_id}" '$1 == group {count++} END{print count+0}' "${manifest}")"

    if [ "${count}" -le 1 ]; then
      log_msg "${group_id}: only one technical replicate, skipping merge"
      continue
    fi

    local pairs_csv
    pairs_csv="$(
      awk -F'\t' -v group="${group_id}" '$1 == group {print $3}' "${manifest}" \
        | paste -sd, -
    )"

    log_msg "${group_id}: merging ${count} technical replicates"
    bash "${SCRIPT_DIR}/merge_replicates.sh" \
      -c "${CONFIG}" \
      -n "${group_id}" \
      -p "${pairs_csv}"
  done
}

build_matrices_from_pairs() {
  local sample="$1"
  local pairs="$2"
  local sample_dir="$3"
  local log_dir="$4"
  local matrix_chrom_sizes="$5"

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
      "${matrix_chrom_sizes}:${BASE_RESOLUTION}" \
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
        "${matrix_chrom_sizes}" \
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

    write_ucsc_hic_track "${sample}" "${norm_hic}" "${matrix_dir}"
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
  local matrix_chrom_sizes="${pairs_dir}/${sample}.matrix.chrom.sizes"
  local select_expr

  mkdir -p "${qc_dir}" "${trimmed_dir}" "${pairs_dir}" "${log_dir}" "${sample_tmp}"
  prepare_matrix_chrom_sizes "${matrix_chrom_sizes}"

  log_msg "${sample}: assay=${assay}"

  if bool_true "${RUN_FASTP:-true}"; then
    if [ ! -s "${trim_r1}" ] || [ ! -s "${trim_r2}" ]; then
      log_msg "${sample}: running fastp"
      if ! fastp \
        -i "${fastq_r1}" \
        -I "${fastq_r2}" \
        -o "${trim_r1}" \
        -O "${trim_r2}" \
        --thread "${THREADS_ALIGN}" \
        --html "${qc_dir}/${sample}.fastp.html" \
        --json "${qc_dir}/${sample}.fastp.json" \
        ${FASTP_EXTRA_ARGS:-} \
        > "${log_dir}/${sample}.fastp.log" 2>&1; then
        echo "ERROR: ${sample}: fastp failed" >&2
        echo "ERROR: fastp log: ${log_dir}/${sample}.fastp.log" >&2
        tail -n 40 "${log_dir}/${sample}.fastp.log" >&2 || true
        exit 1
      fi
    else
      log_msg "${sample}: trimmed FASTQ exists, skipping fastp"
    fi
  else
    trim_r1="${fastq_r1}"
    trim_r2="${fastq_r2}"
  fi

  if [ ! -s "${dedup_pairs}" ]; then
    log_msg "${sample}: aligning, parsing, sorting, and deduplicating pairs"
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
      | filter_pairs_to_chrom_sizes "${matrix_chrom_sizes}" \
      | bgzip -@ "${THREADS_SORT}" \
      > "${valid_pairs}"
  else
    log_msg "${sample}: valid pairs exist, skipping selection"
  fi

  if [ ! -s "${valid_pairs}.px2" ]; then
    pairix "${valid_pairs}"
  fi

  pairtools stats "${valid_pairs}" > "${pairs_dir}/${sample}.pairs.stats.txt"

  build_matrices_from_pairs "${sample}" "${valid_pairs}" "${sample_dir}" "${log_dir}" "${matrix_chrom_sizes}"

  if bool_true "${RUN_MULTIQC:-true}"; then
    log_msg "${sample}: running MultiQC"
    multiqc "${sample_dir}" -n "${sample}.multiqc.html" -o "${qc_dir}" \
      > "${log_dir}/${sample}.multiqc.log" 2>&1 || true
  fi

  log_msg "${sample}: finished"
}

mkdir -p "${OUTDIR}" "${TMPDIR}"
TECH_MANIFEST="${TMPDIR}/microc2tracks.technical_replicates.$$.tsv"
: > "${TECH_MANIFEST}"

tail -n +2 "${SAMPLESHEET}" | while IFS=, read -r sample assay condition biological_replicate technical_replicate fastq_r1 fastq_r2 rest; do
  sample="$(clean_csv_field "${sample}")"
  assay="$(clean_csv_field "${assay}")"
  condition="$(clean_csv_field "${condition}")"
  biological_replicate="$(clean_csv_field "${biological_replicate}")"
  technical_replicate="$(clean_csv_field "${technical_replicate}")"
  fastq_r1="$(clean_csv_field "${fastq_r1}")"
  fastq_r2="$(clean_csv_field "${fastq_r2}")"

  [ -n "${sample// }" ] || continue
  process_sample "${sample}" "${assay}" "${fastq_r1}" "${fastq_r2}"
  assay_group="$(printf '%s' "${assay}" | tr '[:upper:]' '[:lower:]')"
  group_id="$(sanitize_id "${condition}_${assay_group}_B${biological_replicate}_tech_merged")"
  valid_pairs="${OUTDIR}/${sample}/03_pairs/${sample}.valid.mapq${MAPQ_THRESHOLD}.pairs.gz"
  printf '%s\t%s\t%s\n' "${group_id}" "${sample}" "${valid_pairs}" >> "${TECH_MANIFEST}"
done

merge_technical_replicate_groups "${TECH_MANIFEST}"
rm -f "${TECH_MANIFEST}"

log_msg "All samples finished"
