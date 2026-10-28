# Make the tested shared installation the default

The installed `3.0.0-rc4` runtime is ready for default-command promotion after
the user-reported Linux checks recorded on 2026-10-28: 63 regression tests and
four opt-in scientific-tool tests passed. This is the same runtime as commit
`22d10efa98a36234ab798db9096a071aff9a7c4b`; publishing documentation does not
require reinstalling software, rebuilding indexes or rerunning those checks.

Default-command promotion is an administrator action on the server. Publishing
GitHub does not perform it. Keep the version-specific command and old workflow
directories for reproducibility/rollback. The `rc4` identifier is retained so
the published code remains traceable to the version actually tested.

## Promote only when workflow jobs are idle

The following assumes RC4 is already installed and verified at the shared path.
It backs up existing launcher entries (including symlinks), then atomically
replaces only the three explicitly named standard commands. It does not modify
the old home-directory workflow, environment, references, indexes or outputs.
Do not execute this while a job may resolve these launchers again.

```bash
(
set -euo pipefail
RC4_LAUNCHER=/opt/microc2tracks/versions/3.0.0-rc4/scripts/shared_launcher.sh
test -x "$RC4_LAUNCHER"
/usr/local/bin/microc2tracks-rc4 --version

RC4_BACKUP=$(sudo mktemp -d /opt/microc2tracks/launcher-backup.XXXXXXXX)
printf 'Launcher backup: %s\n' "$RC4_BACKUP"

# Finish every backup before changing any standard command.
for RC4_NAME in microc2tracks microc2tracks-preflight microc2tracks-merge
do
    RC4_DEST="/usr/local/bin/$RC4_NAME"
    if [ -e "$RC4_DEST" ] || [ -L "$RC4_DEST" ]; then
        if [ ! -f "$RC4_DEST" ] && [ ! -L "$RC4_DEST" ]; then
            printf 'STOP: not a regular file/symlink: %s\n' "$RC4_DEST"
            exit 1
        fi
        sudo cp -a -- "$RC4_DEST" "$RC4_BACKUP/$RC4_NAME"
    fi
done

for RC4_NAME in microc2tracks microc2tracks-preflight microc2tracks-merge
do
    RC4_DEST="/usr/local/bin/$RC4_NAME"
    RC4_STAGE="/usr/local/bin/.${RC4_NAME}.rc4.$$"
    sudo ln -s -T -- "$RC4_LAUNCHER" "$RC4_STAGE"
    sudo mv -Tf -- "$RC4_STAGE" "$RC4_DEST"
done

# Avoid testing from a directory private to the administrator's home.
cd /tmp
/usr/local/bin/microc2tracks --version
/usr/local/bin/microc2tracks-preflight --version
/usr/local/bin/microc2tracks-merge --version
sudo -u nobody -- /usr/local/bin/microc2tracks --version
)
```

Each version check should print `microc2tracks 3.0.0-rc4`. Then, in the normal
submitting shell, refresh Bash's command lookup and inspect what will run:

```bash
hash -r
type -a microc2tracks microc2tracks-preflight microc2tracks-merge
```

Aliases/functions or an earlier executable directory on an individual user's
PATH may still override the standard commands. Inspect those rather than
rewriting everyone's shell profiles automatically; `/usr/local/bin/microc2tracks`
and `microc2tracks-rc4` remain explicit entry points. Check from another real
login account where site group/ACL rules differ; testing `nobody` alone cannot
prove every possible account policy.

## Future jobs and rollback

Users can run `microc2tracks preflight ...`, `microc2tracks run ...`, or the
`microc2tracks-preflight`/`microc2tracks-merge` aliases without conda activation.
Use a reviewed project config/CSV and new writable output roots. Shared defaults
do not infer actual Hi-C chemistry or repair historical deduplication.

The archived launchers are root-readable backups. Record the printed backup
directory. To roll back, an administrator should inspect the saved entries and
restore the appropriate file/symlink atomically while jobs are idle; do not
delete old workflow directories or broadly reset `/usr/local/bin`.

The existing shared environment/indexes are reused, not copied. Software updates
on GitHub do not automatically change this archive-based installation, which has
no Git checkout requirement. Future versions should be staged separately, tested
and explicitly promoted in the same way.
