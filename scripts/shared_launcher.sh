#!/usr/bin/env bash
# Symlink this file from /usr/local/bin; do not copy it away from the release.
# Run as the submitting user, not sudo. Only installation requires an admin.
set -Eeuo pipefail

LAUNCHER_PATH="$(readlink -f -- "${BASH_SOURCE[0]}")"
RELEASE_DIR="$(cd "$(dirname "$LAUNCHER_PATH")/.." && pwd)"

# This is trusted, administrator-owned site configuration, never a project config.
ENV_PREFIX=/opt/conda/envs/microc2tracks
EXTRA_TOOL_PATH=""
SITE_CONFIG="${MICROC2TRACKS_LAUNCHER_CONFIG:-/etc/microc2tracks/launcher.conf}"
if [[ -e "$SITE_CONFIG" ]]; then
    if [[ ! -f "$SITE_CONFIG" || ! -r "$SITE_CONFIG" ]]; then
        printf 'ERROR: site launcher config is not a readable file: %s\n' "$SITE_CONFIG" >&2
        exit 1
    fi
    source "$SITE_CONFIG"
fi
if [[ "$ENV_PREFIX" != /* ]]; then
    printf 'ERROR: ENV_PREFIX in the site launcher config must be absolute.\n' >&2
    exit 1
fi
if [[ ! -x "$ENV_PREFIX/bin/python" ]]; then
    printf 'ERROR: shared Python is inaccessible: %s/bin/python\n' "$ENV_PREFIX" >&2
    exit 1
fi

# No conda executable, activation script or maintainer home is needed.
# Do not inherit an unrelated active environment's Python import overrides.
export PATH="$ENV_PREFIX/bin${EXTRA_TOOL_PATH:+:$EXTRA_TOOL_PATH}:$PATH"
export CONDA_PREFIX="$ENV_PREFIX"
export PYTHON="$ENV_PREFIX/bin/python"
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
unset PYTHONHOME PYTHONPATH

COMMAND=run
case "$(basename "$0")" in
    microc2tracks-preflight) COMMAND=preflight ;;
    microc2tracks-merge) COMMAND=merge ;;
    *)
        case "${1:-}" in
            run|preflight|merge|check-access) COMMAND="$1"; shift ;;
        esac
        ;;
esac
if [[ "$COMMAND" == check-access ]]; then
    exec "$PYTHON" -B "$RELEASE_DIR/scripts/check_shared_access.py" --env-prefix "$ENV_PREFIX" "$@"
fi
exec "$PYTHON" -B "$RELEASE_DIR/scripts/pipeline.py" "$COMMAND" "$@"
