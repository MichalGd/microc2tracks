#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  run_downstream.sh -c config/config.conf -s sample_name -m sample.norm.mcool

Runs available downstream tools from a balanced .mcool file.
USAGE
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

log_msg() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

CONFIG=""
SAMPLE=""
MATRIX=""

while getopts ":c:s:m:h" opt; do
  case "${opt}" in
    c) CONFIG="${OPTARG}" ;;
    s) SAMPLE="${OPTARG}" ;;
    m) MATRIX="${OPTARG}" ;;
    h) usage; exit 0 ;;
    :) die "Option -${OPTARG} requires an argument" ;;
    \?) die "Unknown option: -${OPTARG}" ;;
  esac
done

[ -n "${CONFIG}" ] || die "Missing -c config file"
[ -n "${SAMPLE}" ] || die "Missing -s sample name"
[ -n "${MATRIX}" ] || die "Missing -m matrix"
[ -f "${CONFIG}" ] || die "Config not found: ${CONFIG}"
[ -f "${MATRIX}" ] || die "Matrix not found: ${MATRIX}"

# shellcheck source=/dev/null
source "${CONFIG}"

DOWNSTREAM_DIR="${OUTDIR}/${SAMPLE}/05_downstream"
LOG_DIR="${OUTDIR}/${SAMPLE}/logs"
mkdir -p \
  "${DOWNSTREAM_DIR}/expected" \
  "${DOWNSTREAM_DIR}/insulation" \
  "${DOWNSTREAM_DIR}/compartments" \
  "${DOWNSTREAM_DIR}/loops" \
  "${DOWNSTREAM_DIR}/saddle" \
  "${LOG_DIR}"

if ! command -v cooltools >/dev/null 2>&1; then
  log_msg "cooltools not found; skipping cooltools downstream analyses"
else
  TAD_URI="${MATRIX}::resolutions/${TAD_RESOLUTION}"
  COMP_URI="${MATRIX}::resolutions/${COMPARTMENT_RESOLUTION}"
  LOOP_URI="${MATRIX}::resolutions/${LOOP_RESOLUTION}"

  EXPECTED_COMP="${DOWNSTREAM_DIR}/expected/${SAMPLE}.expected.${COMPARTMENT_RESOLUTION}.tsv"
  EXPECTED_LOOP="${DOWNSTREAM_DIR}/expected/${SAMPLE}.expected.${LOOP_RESOLUTION}.tsv"
  INSULATION_OUT="${DOWNSTREAM_DIR}/insulation/${SAMPLE}.insulation.${TAD_RESOLUTION}.tsv"
  EIGS_PREFIX="${DOWNSTREAM_DIR}/compartments/${SAMPLE}.eigs.${COMPARTMENT_RESOLUTION}"
  DOTS_OUT="${DOWNSTREAM_DIR}/loops/${SAMPLE}.cooltools_dots.${LOOP_RESOLUTION}.bedpe"

  if [ ! -s "${EXPECTED_COMP}" ]; then
    log_msg "${SAMPLE}: cooltools expected-cis at ${COMPARTMENT_RESOLUTION}"
    cooltools expected-cis \
      -p "${COOLTOOLS_THREADS}" \
      -o "${EXPECTED_COMP}" \
      "${COMP_URI}" \
      > "${LOG_DIR}/${SAMPLE}.expected_compartments.log" 2>&1
  fi

  if [ ! -s "${INSULATION_OUT}" ]; then
    log_msg "${SAMPLE}: cooltools insulation at ${TAD_RESOLUTION}"
    cooltools insulation \
      -p "${COOLTOOLS_THREADS}" \
      -o "${INSULATION_OUT}" \
      "${TAD_URI}" \
      ${INSULATION_WINDOWS_BP} \
      > "${LOG_DIR}/${SAMPLE}.insulation.log" 2>&1
  fi

  if [ -n "${PHASING_TRACK:-}" ] && [ -f "${PHASING_TRACK}" ]; then
    log_msg "${SAMPLE}: cooltools eigs-cis at ${COMPARTMENT_RESOLUTION}"
    cooltools eigs-cis \
      -p "${COOLTOOLS_THREADS}" \
      -o "${EIGS_PREFIX}" \
      --phasing-track "${PHASING_TRACK}" \
      "${COMP_URI}" \
      "${EXPECTED_COMP}" \
      > "${LOG_DIR}/${SAMPLE}.eigs_cis.log" 2>&1 || true
  else
    log_msg "${SAMPLE}: PHASING_TRACK not set; skipping compartment phasing"
  fi

  if [ ! -s "${EXPECTED_LOOP}" ]; then
    log_msg "${SAMPLE}: cooltools expected-cis at ${LOOP_RESOLUTION}"
    cooltools expected-cis \
      -p "${COOLTOOLS_THREADS}" \
      -o "${EXPECTED_LOOP}" \
      "${LOOP_URI}" \
      > "${LOG_DIR}/${SAMPLE}.expected_loops.log" 2>&1
  fi

  if [ ! -s "${DOTS_OUT}" ]; then
    log_msg "${SAMPLE}: cooltools dots at ${LOOP_RESOLUTION}"
    cooltools dots \
      -p "${COOLTOOLS_THREADS}" \
      -o "${DOTS_OUT}" \
      "${LOOP_URI}" \
      "${EXPECTED_LOOP}" \
      > "${LOG_DIR}/${SAMPLE}.cooltools_dots.log" 2>&1 || true
  fi
fi

if command -v mustache >/dev/null 2>&1; then
  MUSTACHE_OUT="${DOWNSTREAM_DIR}/loops/${SAMPLE}.mustache.${LOOP_RESOLUTION}.tsv"
  if [ ! -s "${MUSTACHE_OUT}" ]; then
    log_msg "${SAMPLE}: Mustache loop calling at ${LOOP_RESOLUTION}"
    mustache \
      -f "${MATRIX}" \
      -r "${LOOP_RESOLUTION}" \
      -p "${COOLTOOLS_THREADS}" \
      -o "${MUSTACHE_OUT}" \
      > "${LOG_DIR}/${SAMPLE}.mustache.log" 2>&1 || true
  fi
else
  log_msg "mustache not found; skipping Mustache loop calling"
fi

if command -v chromosight >/dev/null 2>&1; then
  log_msg "chromosight found; stripe/loop modules can be added here when parameters are finalized"
fi

log_msg "${SAMPLE}: downstream finished"

