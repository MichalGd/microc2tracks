#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=../scripts/reference_sentinel.sh
source "${ROOT}/scripts/reference_sentinel.sh"

tmp="$(mktemp)"
trap 'rm -f "${tmp}"' EXIT
printf 'reference_id\tmm39\nassembly\tmm39\n' > "${tmp}"

sentinel_matches_reference "${tmp}" mm39 mm39
if sentinel_matches_reference "${tmp}" hg38 hg38; then
  echo "ERROR: stale mm39 sentinel was accepted for hg38" >&2
  exit 1
fi

printf 'reference_id\thg38\nassembly\thg38\n' > "${tmp}"
sentinel_matches_reference "${tmp}" hg38 hg38
if sentinel_matches_reference "${tmp}" mm39 mm39; then
  echo "ERROR: stale hg38 sentinel was accepted for mm39" >&2
  exit 1
fi

echo "Reference sentinel tests OK"
