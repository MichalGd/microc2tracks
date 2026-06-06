#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

for script in "${ROOT_DIR}"/scripts/*.sh "${ROOT_DIR}"/examples/*.sh; do
  echo "Checking ${script}"
  bash -n "${script}"
done

echo "Bash syntax OK"

