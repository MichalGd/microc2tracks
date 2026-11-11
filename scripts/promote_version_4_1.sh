#!/usr/bin/env bash
# Explicitly promote standard launcher entries; keep installed releases intact.
set -euo pipefail
V41_LAUNCHER="${V41_LAUNCHER:-/opt/microc2tracks/versions/4.1/scripts/shared_launcher.sh}"
V41_BIN_DIR="${V41_BIN_DIR:-/usr/local/bin}"
V41_BACKUP_ROOT="${V41_BACKUP_ROOT:-/opt/microc2tracks}"
case "$V41_LAUNCHER:$V41_BIN_DIR:$V41_BACKUP_ROOT" in
    /*:/*:/*) ;;
    *) printf 'Launcher, binary and backup paths must be absolute.\n' >&2; exit 1 ;;
esac
test -x "$V41_LAUNCHER"
test "$("$V41_BIN_DIR/microc2tracks-4.1" --version)" = 'microc2tracks 4.1'
test "$(readlink -f -- "$V41_BIN_DIR/microc2tracks-4.1")" = "$(readlink -f -- "$V41_LAUNCHER")"
for V41_NAME in microc2tracks microc2tracks-preflight microc2tracks-merge; do
    V41_DEST="$V41_BIN_DIR/$V41_NAME"
    if [ -e "$V41_DEST" ] || [ -L "$V41_DEST" ]; then
        test -f "$V41_DEST" || test -L "$V41_DEST"
    fi
done
sudo -v
V41_BACKUP=$(sudo mktemp -d "$V41_BACKUP_ROOT/launcher-backup-4.1.XXXXXXXX")
printf 'Keep this launcher backup: %s\n' "$V41_BACKUP"
for V41_NAME in microc2tracks microc2tracks-preflight microc2tracks-merge; do
    V41_DEST="$V41_BIN_DIR/$V41_NAME"
    if [ -e "$V41_DEST" ] || [ -L "$V41_DEST" ]; then
        sudo cp -a -- "$V41_DEST" "$V41_BACKUP/$V41_NAME"
    fi
done
for V41_NAME in microc2tracks microc2tracks-preflight microc2tracks-merge; do
    V41_TEMP="$V41_BIN_DIR/.${V41_NAME}.v41.$$"
    sudo ln -s -T -- "$V41_LAUNCHER" "$V41_TEMP"
    sudo mv -Tf -- "$V41_TEMP" "$V41_BIN_DIR/$V41_NAME"
done
cd /tmp
for V41_NAME in microc2tracks microc2tracks-preflight microc2tracks-merge; do
    V41_ACTUAL=$("$V41_BIN_DIR/$V41_NAME" --version)
    printf '%s: %s\n' "$V41_NAME" "$V41_ACTUAL"
    test "$V41_ACTUAL" = 'microc2tracks 4.1'
done
printf 'Default launchers now select 4.1. Run hash -r in your interactive shell.\n'
