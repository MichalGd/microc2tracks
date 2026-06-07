#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  prepare_reference.sh -f reference.fa -a assembly_name -o output_dir

Creates FASTA index, chrom.sizes, and bwa-mem2 index files for the exact FASTA
used by the pipeline.
USAGE
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

FASTA=""
ASSEMBLY=""
OUTDIR=""

while getopts ":f:a:o:h" opt; do
  case "${opt}" in
    f) FASTA="${OPTARG}" ;;
    a) ASSEMBLY="${OPTARG}" ;;
    o) OUTDIR="${OPTARG}" ;;
    h) usage; exit 0 ;;
    :) die "Option -${OPTARG} requires an argument" ;;
    \?) die "Unknown option: -${OPTARG}" ;;
  esac
done

[ -n "${FASTA}" ] || die "Missing -f reference FASTA"
[ -n "${ASSEMBLY}" ] || die "Missing -a assembly name"
[ -n "${OUTDIR}" ] || die "Missing -o output directory"
[ -f "${FASTA}" ] || die "Reference FASTA not found: ${FASTA}"

mkdir -p "${OUTDIR}"

if [ "$(cd "$(dirname "${FASTA}")" && pwd)" != "$(cd "${OUTDIR}" && pwd)" ]; then
  echo "Using FASTA in place: ${FASTA}"
fi

samtools faidx "${FASTA}"
cut -f1,2 "${FASTA}.fai" > "${OUTDIR}/${ASSEMBLY}.chrom.sizes"
bwa-mem2 index "${FASTA}"

echo "Reference prepared:"
echo "  FASTA: ${FASTA}"
echo "  FAI: ${FASTA}.fai"
echo "  chrom sizes: ${OUTDIR}/${ASSEMBLY}.chrom.sizes"
echo "  BWA_INDEX_PREFIX for config.conf: ${FASTA}"
