#!/usr/bin/env bash
# Shared reference-aware completion-sentinel predicate.

sentinel_matches_reference() {
  local done_file="$1"
  local reference_id="$2"
  local assembly="$3"
  local tab=$'\t'

  [ -s "${done_file}" ] \
    && grep -Fqx "reference_id${tab}${reference_id}" "${done_file}" \
    && grep -Fqx "assembly${tab}${assembly}" "${done_file}"
}
