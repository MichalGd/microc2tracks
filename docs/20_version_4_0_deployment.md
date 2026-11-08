# Version 4.0 deployment

Version **4.0**, Git tag **v4.0**, names the native BL-Hi-C implementation
prepared as RC7. Its preset, filtering, matrix behavior and scientific stage
contracts are unchanged. Checkpoint recognition now includes version 4.0.

## Versioned shared installation

The transfer package contains the source tarball, checksum, an installer and
Windows PowerShell transfer instructions. Use the installer copied beside the
tarball, not the copy still inside an unextracted source archive.

On biolserv, from the transferred package directory:

```bash
sha256sum -c SHA256SUMS
bash install_version_4_0.sh --check
bash install_version_4_0.sh
/usr/local/bin/microc2tracks-4.0 --version
```

The installer verifies the archive, reuses the existing
`/opt/conda/envs/microc2tracks` environment (or the existing root-owned launcher
configuration), checks build prerequisites, builds the pinned standalone linker
and runs the 14 tiny BL tests before copying into:

```text
/opt/microc2tracks/versions/4.0
/usr/local/bin/microc2tracks-4.0
```

It requires sudo for the new shared directory/command, refuses to overwrite an
existing destination, grants ordinary users read/execute access to the new
release, and retains its staging directory for inspection. No apt/conda update,
reference/index generation, raw-data download or mapping is performed. Existing
installed releases, shared environment and standard command symlinks are left
in place. Concurrent existing analyses can continue using their original release.

If g++/OpenMP/zlib headers are unavailable, the prerequisite check stops before
installation. Resolve those build prerequisites separately; the script never
installs system packages automatically. The optional trimmer source/license are
bundled and the build makes no network request.

The version check establishes executable availability for the current user;
check another real login account if site ACL/group policies differ. Existing
project references/data keep their existing permissions.

## Use version 4.0 for future runs

Copy/edit normal project config and samplesheet files, choosing fresh output
and bulk scratch paths where needed, then run:

```bash
/usr/local/bin/microc2tracks-4.0 preflight -c /path/config.conf -s /path/samplesheet.csv
nohup /usr/local/bin/microc2tracks-4.0 run -c /path/config.conf -s /path/samplesheet.csv > mapping_v4.log 2>&1 < /dev/null &
echo $! > mapping_v4.pid
tail -f mapping_v4.log
```

The `/path/...` entries are placeholders. After the nohup launch succeeds,
disconnecting SSH does not stop the workflow. BL rows use `assay=hic` and
`hic_protocol=blhic`, with HaeIII and raw paired inputs as described in the
[BL guide](19_blhic_preprocessing.md). Arima N/GANTC support remains separate.

## Optional standard-command promotion

The versioned command already uses 4.0. To make the usual `microc2tracks`,
`microc2tracks-preflight` and `microc2tracks-merge` commands select it, wait until
active workflow jobs finish, then back up and replace only those launcher entries.
Do not use a broad recursive update of old installations or environments.

```bash
(
set -euo pipefail
V4_LAUNCHER=/opt/microc2tracks/versions/4.0/scripts/shared_launcher.sh
test -x "$V4_LAUNCHER"
test "$(/usr/local/bin/microc2tracks-4.0 --version)" = 'microc2tracks 4.0'
V4_BACKUP=$(sudo mktemp -d /opt/microc2tracks/launcher-backup-4.0.XXXXXXXX)
printf 'Keep this backup directory: %s\n' "$V4_BACKUP"
for V4_NAME in microc2tracks microc2tracks-preflight microc2tracks-merge; do
    V4_DEST="/usr/local/bin/$V4_NAME"
    if [ -e "$V4_DEST" ] || [ -L "$V4_DEST" ]; then
        test -f "$V4_DEST" || test -L "$V4_DEST"
        sudo cp -a -- "$V4_DEST" "$V4_BACKUP/$V4_NAME"
    fi
done
for V4_NAME in microc2tracks microc2tracks-preflight microc2tracks-merge; do
    V4_TEMP="/usr/local/bin/.${V4_NAME}.v4.$$"
    sudo ln -s -- "$V4_LAUNCHER" "$V4_TEMP"
    sudo mv -Tf -- "$V4_TEMP" "/usr/local/bin/$V4_NAME"
done
cd /tmp
/usr/local/bin/microc2tracks --version
/usr/local/bin/microc2tracks-preflight --version
/usr/local/bin/microc2tracks-merge --version
)
hash -r
type -a microc2tracks
```

Keep older version directories and the printed launcher backup for rollback.
No server deployment is performed merely by publishing the GitHub tag.
