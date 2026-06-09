#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  run_downstream.sh [-l] -c config/config.conf -s sample_name -m sample.norm.mcool [-o output_base_dir]

Runs available downstream tools from a balanced .mcool file.

Options:
  -l  Light mode: expected/insulation/TADs/compartments only; skip loop callers.
  -o  Output base directory. Defaults to OUTDIR/sample_name.
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

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG=""
SAMPLE=""
MATRIX=""
OUTPUT_BASE=""
LIGHT_MODE="false"

while getopts ":c:s:m:o:lh" opt; do
  case "${opt}" in
    c) CONFIG="${OPTARG}" ;;
    s) SAMPLE="${OPTARG}" ;;
    m) MATRIX="${OPTARG}" ;;
    o) OUTPUT_BASE="${OPTARG}" ;;
    l) LIGHT_MODE="true" ;;
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

OUTPUT_BASE="${OUTPUT_BASE:-${OUTDIR}/${SAMPLE}}"
DOWNSTREAM_DIR="${OUTPUT_BASE}/05_downstream"
LOG_DIR="${OUTPUT_BASE}/logs"
mkdir -p \
  "${DOWNSTREAM_DIR}/expected" \
  "${DOWNSTREAM_DIR}/insulation" \
  "${DOWNSTREAM_DIR}/insulation/bedgraph" \
  "${DOWNSTREAM_DIR}/tads" \
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
  INSULATION_OUT="${DOWNSTREAM_DIR}/insulation/${SAMPLE}.insulation.${TAD_RESOLUTION}.tsv"
  INSULATION_BEDGRAPH_DIR="${DOWNSTREAM_DIR}/insulation/bedgraph"
  INSULATION_BEDGRAPH_SENTINEL="${INSULATION_BEDGRAPH_DIR}/${SAMPLE}.insulation_bedgraph.done"
  TADS_OUT="${DOWNSTREAM_DIR}/tads/${SAMPLE}.tads.${TAD_BOUNDARY_WINDOW_BP:-100000}.bed"
  EIGS_PREFIX="${DOWNSTREAM_DIR}/compartments/${SAMPLE}.eigs.${COMPARTMENT_RESOLUTION}"
  EXPECTED_LOOP="${DOWNSTREAM_DIR}/expected/${SAMPLE}.expected.${LOOP_RESOLUTION}.tsv"
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
      --threshold "${TAD_BOUNDARY_THRESHOLD:-Li}" \
      -p "${COOLTOOLS_THREADS}" \
      -o "${INSULATION_OUT}" \
      "${TAD_URI}" \
      ${INSULATION_WINDOWS_BP} \
      > "${LOG_DIR}/${SAMPLE}.insulation.log" 2>&1
  fi

  if [ ! -s "${INSULATION_BEDGRAPH_SENTINEL}" ]; then
    log_msg "${SAMPLE}: exporting insulation bedGraph tracks"
    BEDGRAPH_ARGS=(
      --insulation "${INSULATION_OUT}"
      --out-dir "${INSULATION_BEDGRAPH_DIR}"
      --sample "${SAMPLE}"
      --windows "${INSULATION_WINDOWS_BP}"
    )
    if bool_true "${INSULATION_BEDGRAPH_TRACKLINE:-true}"; then
      BEDGRAPH_ARGS+=(--track-line)
    fi
    python "${SCRIPT_DIR}/export_insulation_bedgraph.py" "${BEDGRAPH_ARGS[@]}" \
      > "${LOG_DIR}/${SAMPLE}.insulation_bedgraph.log" 2>&1 || true
    if ls "${INSULATION_BEDGRAPH_DIR}/${SAMPLE}.insulation_score."*.bedGraph >/dev/null 2>&1; then
      date '+%Y-%m-%d %H:%M:%S' > "${INSULATION_BEDGRAPH_SENTINEL}"
    fi
  fi

  if [ ! -s "${TADS_OUT}" ]; then
    log_msg "${SAMPLE}: calling TAD-like intervals from insulation boundaries"
    TAD_ARGS=(
      --insulation "${INSULATION_OUT}"
      --out "${TADS_OUT}"
      --window-bp "${TAD_BOUNDARY_WINDOW_BP:-100000}"
      --min-size-bp "${TAD_MIN_SIZE_BP:-40000}"
    )
    if [ -n "${TAD_MIN_BOUNDARY_STRENGTH:-}" ]; then
      TAD_ARGS+=(--min-boundary-strength "${TAD_MIN_BOUNDARY_STRENGTH}")
    fi
    python "${SCRIPT_DIR}/call_tads_from_insulation.py" "${TAD_ARGS[@]}" \
      > "${LOG_DIR}/${SAMPLE}.tads.log" 2>&1 || true
  fi

  if [ -n "${PHASING_TRACK:-}" ] && [ -f "${PHASING_TRACK}" ]; then
    log_msg "${SAMPLE}: cooltools eigs-cis at ${COMPARTMENT_RESOLUTION} with phasing track"
    cooltools eigs-cis \
      -p "${COOLTOOLS_THREADS}" \
      -o "${EIGS_PREFIX}" \
      --phasing-track "${PHASING_TRACK}" \
      "${COMP_URI}" \
      > "${LOG_DIR}/${SAMPLE}.eigs_cis.log" 2>&1 || true
  else
    log_msg "${SAMPLE}: cooltools eigs-cis at ${COMPARTMENT_RESOLUTION} without phasing track; PC1 sign is arbitrary"
    cooltools eigs-cis \
      -p "${COOLTOOLS_THREADS}" \
      -o "${EIGS_PREFIX}" \
      "${COMP_URI}" \
      > "${LOG_DIR}/${SAMPLE}.eigs_cis.log" 2>&1 || true
  fi

  if bool_true "${LIGHT_MODE:-false}"; then
    log_msg "${SAMPLE}: light downstream mode; skipping loop callers"
  elif [ ! -s "${EXPECTED_LOOP}" ]; then
    log_msg "${SAMPLE}: cooltools expected-cis at ${LOOP_RESOLUTION}"
    cooltools expected-cis \
      -p "${COOLTOOLS_THREADS}" \
      -o "${EXPECTED_LOOP}" \
      "${LOOP_URI}" \
      > "${LOG_DIR}/${SAMPLE}.expected_loops.log" 2>&1
  fi

  if ! bool_true "${LIGHT_MODE:-false}" && [ ! -s "${DOTS_OUT}" ]; then
    log_msg "${SAMPLE}: cooltools dots at ${LOOP_RESOLUTION}"
    cooltools dots \
      -p "${COOLTOOLS_THREADS}" \
      -o "${DOTS_OUT}" \
      "${LOOP_URI}" \
      "${EXPECTED_LOOP}" \
      > "${LOG_DIR}/${SAMPLE}.cooltools_dots.log" 2>&1 || true
  fi
fi

if bool_true "${LIGHT_MODE:-false}"; then
  :
else
  MUSTACHE_CMD=()
  if command -v mustache >/dev/null 2>&1; then
    MUSTACHE_CMD=(mustache)
  elif python -m mustache --help >/dev/null 2>&1; then
    MUSTACHE_CMD=(python -m mustache)
  fi

  if [ "${#MUSTACHE_CMD[@]}" -gt 0 ]; then
    MUSTACHE_OUT="${DOWNSTREAM_DIR}/loops/${SAMPLE}.mustache.${LOOP_RESOLUTION}.tsv"
    if [ ! -s "${MUSTACHE_OUT}" ]; then
      log_msg "${SAMPLE}: Mustache loop calling at ${LOOP_RESOLUTION}"
      "${MUSTACHE_CMD[@]}" \
        -f "${MATRIX}" \
        -r "${LOOP_RESOLUTION}" \
        -p "${COOLTOOLS_THREADS}" \
        -o "${MUSTACHE_OUT}" \
        > "${LOG_DIR}/${SAMPLE}.mustache.log" 2>&1 || true
    fi
  else
    log_msg "mustache not found; skipping Mustache loop calling"
  fi
fi

if ! bool_true "${LIGHT_MODE:-false}" && command -v chromosight >/dev/null 2>&1; then
  log_msg "chromosight found; stripe/loop modules can be added here when parameters are finalized"
fi

log_msg "${SAMPLE}: downstream finished"
