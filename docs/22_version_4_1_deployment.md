# Version 4.1 server update

Version **4.1**, Git tag **v4.1**, adds native `arima_2enzyme` Hi-C and full
DNA IUPAC cut patterns in the existing config/samplesheet enzyme fields.
See [the protocol guide](21_arima_hic.md) for syntax and scientific limits.

## Install a separate shared release

The prepared deployment bundle contains `microc2tracks_4.1_20261111.tar.gz`,
its `.sha256`, `install_version_4_1.sh`, `promote_version_4_1.sh`, a checksum
manifest, a PowerShell transfer helper and complete transfer instructions.
Run the installer **beside the transferred archive**, not from inside the
source checkout without its corresponding archive/checksum.

On biolserv in your Linux SSH terminal:

```bash
cd /home/micgdu/Analysis/workflows/releases/microc2tracks_4.1_2026-11-11
sha256sum -c SHA256SUMS
bash install_version_4_1.sh --check
bash install_version_4_1.sh
/usr/local/bin/microc2tracks-4.1 --version
```

Expected final version: `microc2tracks 4.1`. Installation creates:

```text
/opt/microc2tracks/versions/4.1
/usr/local/bin/microc2tracks-4.1
```

The installer reuses `/opt/conda/envs/microc2tracks` or your existing trusted
site launcher configuration. It verifies the source checksum and prerequisites,
builds the bundled pinned linker, and runs 29 tiny Arima/BL/checkpoint/launcher
tests as the submitting user. Only after these pass does it request sudo for
the new shared directory and command. All-user source read/execute permissions
are applied to that new release. Existing data/reference permissions are not
changed; use `microc2tracks-4.1 check-access` with the documented project inputs
when validating another user's actual access.

No apt/conda operations, reference building, sequencing downloads or mapping
are performed. Existing releases, environments and default commands remain
in place. Existing jobs using their resolved old release can continue while
4.1 is installed. The small build/tests use some CPU and local disk but do not
replace files used by those jobs. A future job uses 4.1 only when its launcher
selects 4.1.

`--check` verifies prerequisites without installing or creating staging. A real
installation retains its staging folder for inspection. Any existing target or
versioned launcher, including a broken symlink, is refused. If installation is
interrupted, inspect the printed staging/target before retrying; the installer
will not delete or overwrite a partial installed release automatically.

The optional linker build needs the already available g++, OpenMP and zlib
headers. The check stops if missing; no package manager is invoked. The 29 tests
exercise small fixtures, not a production BWA/cooler/Java mapping pilot.

## Select 4.1 in a project

Use the new versioned command in the normal project working directory:

```bash
/usr/local/bin/microc2tracks-4.1 preflight -c /absolute/project/config.conf -s /absolute/project/samplesheet.csv
nohup /usr/local/bin/microc2tracks-4.1 run -c /absolute/project/config.conf -s /absolute/project/samplesheet.csv > mapping_v41.log 2>&1 < /dev/null &
echo $! > mapping_v41.pid
tail -f mapping_v41.log
```

The `/absolute/project/...` entries are placeholders. Choose real project paths,
and proceed to mapping only after preflight succeeds. The update procedure
above does not start this analysis. Outputs remain controlled by the project's
config; installing 4.1 does not relocate reads or results.

For Chen's confirmed two-pattern Arima dataset, use `assay=hic` and:

```bash
HIC_PROTOCOL="arima_2enzyme"
HIC_RESTRICTION_ENZYMES=""
```

Use the new bundled protocol registry: remove an old
`HIC_PROTOCOL_REGISTRY="hic_protocols_chen.json"` override, or explicitly select
`/opt/microc2tracks/versions/4.1/config/hic_protocols.json`. Existing row-specific
settings override global settings; update `hic_protocol`/`restriction_enzymes`
in rows if they contain old values. Explicit `DpnII+HinfI` or `^GATC+G^ANTC`
is also valid for the preset. Keep the established mm39 references and separate
WT/DKO sample selection.

The earlier private Chen launcher still invokes its private runtime. It is not
converted by installing or promoting 4.1. Use the native versioned command and
a reviewed project config; do not run the old compatibility preparer over the
new configuration. Do not launch a second process into an active OUTDIR.
Checkpoint migration from RC4, `4.1-dev`, or private `4.0+chen-arimaN1` is not
newly certified. Existing provenance checks remain active.

## Optional: make 4.1 the default

After the versioned installation passes, run this separate command if future
standard invocations should use 4.1:

```bash
cd /home/micgdu/Analysis/workflows/releases/microc2tracks_4.1_2026-11-11
bash promote_version_4_1.sh
hash -r
type -a microc2tracks
microc2tracks --version
```

This backs up and atomically replaces only the three default launchers
`microc2tracks`, `microc2tracks-preflight`, and `microc2tracks-merge`. Both
versioned commands and old source trees are retained. Already running workflow
processes keep their resolved release paths. A custom controller that later
launches another bare `microc2tracks` command would select the new default;
schedule promotion between jobs when that distinction matters. Inspect `type -a`
for a shell alias or earlier PATH entry hiding `/usr/local/bin/microc2tracks`.

For rollback of the defaults, retain the printed backup directory. Set
`V41_RESTORE` below to that **exact** path and restore each saved launcher:

```bash
(
set -euo pipefail
V41_RESTORE='/opt/microc2tracks/launcher-backup-4.1.REPLACE_WITH_ACTUAL_SUFFIX'
sudo test -d "$V41_RESTORE"
for V41_NAME in microc2tracks microc2tracks-preflight microc2tracks-merge; do
    if sudo test -e "$V41_RESTORE/$V41_NAME" || sudo test -L "$V41_RESTORE/$V41_NAME"; then
        sudo cp -a --remove-destination -- "$V41_RESTORE/$V41_NAME" "/usr/local/bin/$V41_NAME"
    fi
done
)
hash -r
/usr/local/bin/microc2tracks --version
```

That restores entries that existed before promotion. Keep the versioned 4.0
and 4.1 directories; rollback does not require deleting either release.
The backup directory is root-owned; the existence checks use sudo as well as
the copy so an ordinary account does not silently skip protected backup entries.
