#!/usr/bin/env bash
# Run the transferred copy beside the verified source archive/checksum.
# Install only the new version; existing releases and default commands are retained.
set -euo pipefail

V41_CHECK=false
if [ "$#" -gt 1 ]; then
    printf 'Usage: bash install_version_4_1.sh [--check]\n' >&2
    exit 2
fi
case "${1:-}" in
    --check) V41_CHECK=true ;;
    '') ;;
    *) printf 'Usage: bash install_version_4_1.sh [--check]\n' >&2; exit 2 ;;
esac
V41_PACKAGE_DIR="${V41_PACKAGE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
V41_TARGET="${V41_TARGET:-/opt/microc2tracks/versions/4.1}"
V41_BIN_DIR="${V41_BIN_DIR:-/usr/local/bin}"
V41_LINK="$V41_BIN_DIR/microc2tracks-4.1"
case "$V41_TARGET" in /*/4.1) ;; *) printf 'Release path must be absolute and end in /4.1\n' >&2; exit 1;; esac
case "$V41_BIN_DIR" in /*) ;; *) printf 'Binary directory must be absolute\n' >&2; exit 1;; esac

# Reuse the administrator's existing environment without changing packages.
ENV_PREFIX=/opt/conda/envs/microc2tracks
EXTRA_TOOL_PATH=""
V41_SITE_CONFIG="${MICROC2TRACKS_LAUNCHER_CONFIG:-/etc/microc2tracks/launcher.conf}"
if [ -e "$V41_SITE_CONFIG" ]; then
    test -f "$V41_SITE_CONFIG" && test -r "$V41_SITE_CONFIG" || {
        printf 'Site launcher config is not a readable file: %s\n' "$V41_SITE_CONFIG" >&2
        exit 1
    }
    source "$V41_SITE_CONFIG"
fi
case "$ENV_PREFIX" in /*) ;; *) printf 'ENV_PREFIX must be absolute\n' >&2; exit 1;; esac
test -x "$ENV_PREFIX/bin/python" || {
    printf 'Existing shared Python missing: %s/bin/python\n' "$ENV_PREFIX" >&2
    exit 1
}
export PATH="$ENV_PREFIX/bin${EXTRA_TOOL_PATH:+:$EXTRA_TOOL_PATH}:$PATH"
export CONDA_PREFIX="$ENV_PREFIX"
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
unset PYTHONHOME PYTHONPATH

cd "$V41_PACKAGE_DIR"
V41_PACKAGE_DIR="$(pwd -P)"
sha256sum --strict -c microc2tracks_4.1_20261111.sha256
if [ -e "$V41_TARGET" ] || [ -L "$V41_TARGET" ] || [ -e "$V41_LINK" ] || [ -L "$V41_LINK" ]; then
    printf 'STOP: release/launcher already exists; no existing installation will be overwritten.\n' >&2
    exit 1
fi
command -v g++ >/dev/null || {
    printf 'Building the optional linker needs g++, OpenMP and zlib headers.\n' >&2
    exit 1
}
printf '#include <zlib.h>\n#include <omp.h>\nint main() { return 0; }\n' | g++ -fopenmp -x c++ -fsyntax-only -
printf 'Existing environment: %s\nNew release: %s\nNew command: %s\n' "$ENV_PREFIX" "$V41_TARGET" "$V41_LINK"
if "$V41_CHECK"; then
    printf 'Checks passed. No installation, staging directory, mapping, data download or package-manager operation performed.\n'
    exit 0
fi

# Build and run tiny local checks as the submitting user, before requesting sudo.
V41_STAGE=$(mktemp -d "$V41_PACKAGE_DIR/stage_4.1.XXXXXXXX")
printf 'Retained staging directory: %s\n' "$V41_STAGE"
tar -xzf microc2tracks_4.1_20261111.tar.gz -C "$V41_STAGE" --no-same-owner
V41_SOURCE="$V41_STAGE/microc2tracks"
test "$("$ENV_PREFIX/bin/python" -B "$V41_SOURCE/scripts/pipeline.py" run --version)" = 'microc2tracks 4.1'
"$ENV_PREFIX/bin/python" -B "$V41_SOURCE/scripts/build_trimlinker.py"
(
    cd "$V41_SOURCE/tests"
    PYTHONPATH="$V41_SOURCE/tests:$V41_SOURCE/scripts" "$ENV_PREFIX/bin/python" -B -m unittest -v \
        test_arima test_blhic \
        test_rc6_safety.SafetyTests.test_resume_with_removed_fastqs_uses_recorded_lineage_and_rejects_science_changes \
        test_rc6_safety.SafetyTests.test_completed_library_resumes_after_raw_and_sorted_inputs_removed \
        test_shared_installation.SharedInstallationTests.test_actual_version_runs_without_writes_in_release
)

sudo -v
# Recheck immediately before creating the new directory and versioned launcher.
test ! -e "$V41_TARGET" && test ! -L "$V41_TARGET"
test ! -e "$V41_LINK" && test ! -L "$V41_LINK"
# -m applies only to newly created directories; existing parent modes stay intact.
sudo mkdir -p -m 0755 -- "$(dirname "$V41_TARGET")" "$V41_BIN_DIR"
sudo mkdir -m 0755 -- "$V41_TARGET"
sudo cp -a -- "$V41_SOURCE/." "$V41_TARGET/"
sudo chown -R root:root -- "$V41_TARGET"
sudo chmod -R a+rX,go-w -- "$V41_TARGET"
sudo ln -s -- "$V41_TARGET/scripts/shared_launcher.sh" "$V41_LINK"
cd /tmp
test "$("$V41_LINK" --version)" = 'microc2tracks 4.1'
printf 'Installed: %s\n' "$V41_LINK"
printf 'The existing default commands, environments, references and running releases were not changed.\n'
