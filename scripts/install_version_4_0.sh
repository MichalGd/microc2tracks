#!/usr/bin/env bash
# Run from the transferred package. Installs only a new version-specific release.
set -euo pipefail

V4_CHECK=false
case "${1:-}" in
    --check) V4_CHECK=true ;;
    '') ;;
    *) printf 'Usage: bash install_version_4_0.sh [--check]\n' >&2; exit 2 ;;
esac
V4_PACKAGE_DIR="${V4_PACKAGE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
V4_TARGET="${V4_TARGET:-/opt/microc2tracks/versions/4.0}"
V4_BIN_DIR="${V4_BIN_DIR:-/usr/local/bin}"
V4_LINK="$V4_BIN_DIR/microc2tracks-4.0"
case "$V4_TARGET" in /*/4.0) ;; *) printf 'Release path must be absolute and end in /4.0\n' >&2; exit 1;; esac
case "$V4_BIN_DIR" in /*) ;; *) printf 'Binary directory must be absolute\n' >&2; exit 1;; esac

# Match the existing shared launcher's environment without changing it.
ENV_PREFIX=/opt/conda/envs/microc2tracks
EXTRA_TOOL_PATH=""
V4_SITE_CONFIG="${MICROC2TRACKS_LAUNCHER_CONFIG:-/etc/microc2tracks/launcher.conf}"
if [ -f "$V4_SITE_CONFIG" ]; then source "$V4_SITE_CONFIG"; fi
test -x "$ENV_PREFIX/bin/python" || { printf 'Existing shared Python missing: %s/bin/python\n' "$ENV_PREFIX" >&2; exit 1; }
export PATH="$ENV_PREFIX/bin${EXTRA_TOOL_PATH:+:$EXTRA_TOOL_PATH}:$PATH"
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
unset PYTHONHOME PYTHONPATH

cd "$V4_PACKAGE_DIR"
sha256sum -c microc2tracks_4.0_20261108.sha256
if [ -e "$V4_TARGET" ] || [ -L "$V4_TARGET" ] || [ -e "$V4_LINK" ] || [ -L "$V4_LINK" ]; then
    printf 'STOP: release/launcher already exists; no existing installation will be overwritten.\n' >&2
    exit 1
fi
command -v g++ >/dev/null || { printf 'Building the optional linker needs g++, OpenMP and zlib headers.\n' >&2; exit 1; }
printf '#include <zlib.h>\n#include <omp.h>\nint main() { return 0; }\n' | g++ -fopenmp -x c++ -fsyntax-only -
printf 'Existing environment: %s\nNew release: %s\nNew command: %s\n' "$ENV_PREFIX" "$V4_TARGET" "$V4_LINK"
if "$V4_CHECK"; then
    printf 'Checks passed. No installation, mapping, data download or package-manager operation performed.\n'
    exit 0
fi

# Build/test as the submitting user, before copying into the new shared path.
V4_STAGE=$(mktemp -d "$V4_PACKAGE_DIR/stage_4.0.XXXXXXXX")
printf 'Retained staging directory: %s\n' "$V4_STAGE"
tar -xzf microc2tracks_4.0_20261108.tar.gz -C "$V4_STAGE" --no-same-owner
V4_SOURCE="$V4_STAGE/microc2tracks"
test "$("$ENV_PREFIX/bin/python" -B "$V4_SOURCE/scripts/pipeline.py" run --version)" = 'microc2tracks 4.0'
"$ENV_PREFIX/bin/python" -B "$V4_SOURCE/scripts/build_trimlinker.py"
(
    cd "$V4_SOURCE"
    PYTHONPATH=tests:scripts "$ENV_PREFIX/bin/python" -B -m unittest test_blhic -v
)

sudo -v
# Recheck immediately before creating the new installation and launcher.
test ! -e "$V4_TARGET" && test ! -L "$V4_TARGET"
test ! -e "$V4_LINK" && test ! -L "$V4_LINK"
sudo mkdir -p -- "$(dirname "$V4_TARGET")" "$V4_BIN_DIR"
sudo mkdir -- "$V4_TARGET"
sudo cp -a -- "$V4_SOURCE/." "$V4_TARGET/"
sudo chown -R root:root -- "$V4_TARGET"
sudo chmod -R a+rX,go-w -- "$V4_TARGET"
sudo ln -s -- "$V4_TARGET/scripts/shared_launcher.sh" "$V4_LINK"
cd /tmp
test "$("$V4_LINK" --version)" = 'microc2tracks 4.0'
printf 'Installed: %s\n' "$V4_LINK"
printf 'The existing default commands, environments, references and running releases were not changed.\n'
