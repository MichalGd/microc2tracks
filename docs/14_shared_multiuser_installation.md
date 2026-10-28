# Shared, versioned Linux installation (3.0.0-rc4)

Install code outside every user's home. All ordinary local accounts must be able
to read code/references and execute tools; only administrators/maintainers may
change shared software. Each run has user/project-owned writable data. Making
the software available does **not** make private sample data readable to others.

## Layout and unchanged running work

```text
/opt/microc2tracks/versions/3.0.0-rc4/       root-owned release, dirs 0755/files 0644
/usr/local/bin/microc2tracks-rc4            symlink to that release's shared_launcher.sh
/etc/microc2tracks/launcher.conf            optional root-owned site config, 0644
/opt/conda/envs/microc2tracks/               EXISTING environment, not copied/updated
existing reference/index/JAR paths          reused read-only, no genome reindexing
user/project work directory                own config, sample sheet and unique OUTDIR
```

The shared launcher resolves its symlink to a specific release and directly runs
that environment's Python with its `bin` first on PATH. It does not source a
maintainer's conda initialization or require `conda activate`, `conda` on PATH,
or software in a maintainer's home. User-site packages/PYTHONPATH are disabled for
reproducibility; Python bytecode is not written into the read-only installation.
Activation-hook-dependent external tools still need the real-tool checks below.

Leave the running workflow path, environment, configs, references and outputs
untouched. Do not replace an existing launcher that an active job may call. Add
only the version-specific command first; promotion of a default command can wait.
Nothing here launches a production run, changes packages or copies genome indexes.

## Administrator staging: verified archive only

Upload the RC4 archive and checksum, e.g. to the submitting user's existing
`/home/micgdu/Analysis/workflows/incoming`. This is only a temporary transfer
location, not the installation. The following example uses that upload location;
other installations can use any administrator-accessible upload directory.

```bash
(
set -e
UPLOAD=/home/micgdu/Analysis/workflows/incoming
RELEASE=/opt/microc2tracks/versions/3.0.0-rc4
cd "$UPLOAD"
sha256sum -c microc2tracks_3.0.0-rc4_20261028.sha256

if sudo test -e "$RELEASE" || sudo test -L "$RELEASE"; then
    printf 'STOP: release path already exists: %s\n' "$RELEASE"
    exit 1
fi
sudo install -d -o root -g root -m 0755 /opt/microc2tracks/versions
sudo mkdir -- "$RELEASE"
sudo tar --no-same-owner -xzf microc2tracks_3.0.0-rc4_20261028.tar.gz \
    --strip-components=1 -C "$RELEASE"

# These permissions apply ONLY to the newly extracted release, not old software.
sudo chown -R root:root -- "$RELEASE"
sudo find "$RELEASE" -type d -exec chmod 0755 {} +
sudo find "$RELEASE" -type f -exec chmod 0644 {} +
sudo find "$RELEASE/scripts" "$RELEASE/examples" -type f -name '*.sh' -exec chmod 0755 {} +

# No -f: refuse to replace an existing command/symlink.
sudo ln -s -- "$RELEASE/scripts/shared_launcher.sh" /usr/local/bin/microc2tracks-rc4
)
```

First inspect existing `/opt/microc2tracks` components if they already exist.
They must be searchable/readable, and not writable by ordinary users. Do not
recursively change permissions on an existing installation or shared environment
while analyses are running. Never use `chmod 777` on software, environments,
indexes or private data.

The launcher defaults to `/opt/conda/envs/microc2tracks`. A different existing
environment or necessary **shared** extra executable directories can be specified
in `/etc/microc2tracks/launcher.conf`, using `config/launcher.conf.example`. Copy
it only if that site config does not already exist; review existing settings
instead of overwriting them. Keep the file and its directory administrator-owned,
not writable by ordinary users. No change is required for the default prefix.

## Check from an unrelated account, not just the maintainer

Every parent directory needs search (`x`) permission, data/scripts need read (`r`),
and executable tools need execute (`x`). A symlink cannot bypass missing access to
its target/ancestors. An environment in `/opt` can still contain symlinks into a
private home, so the pathname alone does not prove all-user availability.

```bash
/usr/local/bin/microc2tracks-rc4 --version
/usr/local/bin/microc2tracks-rc4 check-access
sudo -u nobody -- /usr/local/bin/microc2tracks-rc4 --version
sudo -u nobody -- /usr/local/bin/microc2tracks-rc4 check-access
```

The `nobody` checks exercise an unrelated, unprivileged identity without starting
an analysis. They do not replace a check in another real login account, especially
where ACLs/groups differ. `check-access` reports the effective identity and checks
read/search/execute access; it does not install, execute scientific tools or hash
whole genomes. Missing reference/JAR arguments are explicitly **not checked**.

Also run `check-access --registry ... --reference ... --juicer-jar ...` with the
actual shared registry and JAR paths, both as the submitting user and another
account. `--reference` is repeatable; without it all registry entries are checked.
`--refilter` skips FASTQ/BWA checks; `--no-hic` skips Java/JAR checks if `.hic` is
disabled. It reads existing `.fai`/chromosome sizes, never builds indexes. Any
legacy reference overrides in a project config need their own project preflight.

If access fails, inspect the exact failing path with `namei -l` (and ACLs if in
use). Stop and have the administrator decide the narrow correction. Do **not**
open up `/home/micgdu`, chmod the active environment recursively, move references,
or copy large resources just to hide a permission problem. An administrator may
choose a narrowly scoped traversal/read ACL or an explicitly managed read-only
bind mount for an existing reference directory; neither is applied automatically.

Passing as one identity does not prove every server account has access. Site
permissions/ACL policy and another normal-account test are required for that
claim. Shared installed code and the environment must remain non-writable to
ordinary users; do not execute untrusted project configs or analyses with sudo.

## Each user's project

From a user-owned project directory, copy only small config/sample-sheet templates
and edit those copies. In the user's config set absolute existing reference/JAR
paths, a new unique `OUTDIR`, and optional alternative roots. For example:

```bash
microc2tracks-rc4 preflight -c config.conf -s samplesheet.csv
microc2tracks-rc4 run -c config.conf -s samplesheet.csv
```

Refiltering uses the same command with `--mode refilter-deduplicated` and its
separate CSV. Direct Bash entry points remain available for compatibility.
Do not run preflight until the project config/paths/resources have been reviewed.

Do not make every user's output world-writable. Different users/runs need distinct
OUTDIR and category roots, or deliberate project-group access. Empty optional
roots preserve the standard layout under OUTDIR. Leave the restriction cache at
its per-run default (`OUTDIR/reference_cache`) unless a trusted writable project
cache is deliberately configured. The current cache locks need write access; a
root-owned read-only shared cache is not transparently supported. Do not use a
world-writable shared cache. Scratch is created in unique private run directories;
its configured parent must be writable/searchable by the submitting user.

Sharing software does not share Java heaps or coordinate CPU/RAM/I/O between
users. Concurrency limits are per run, not server-wide. Use site scheduling/quotas
or coordinate jobs. Separate folders on the same filesystem do not create extra
capacity or protect the system disk; configure actual data mounts/reserves.

## Minimum validation and promotion

Reuse the current environment without changing packages if imports and tiny
functional tests pass. The proposed `environment.yml` is a fresh-install target,
not an instruction to replace this server's active environment. Defer tests if
their shared resource use could affect active work. In the staged release, run
the unit suite and the four opt-in toy scientific tests in TEST_REPORT.md with
the existing environment/JAR. These tests create a tiny temporary synthetic BWA
index, not a production index. No production FASTQs or existing outputs are used.

Only after access checks, functional validation and review, add an unversioned
`/usr/local/bin/microc2tracks` symlink to this release's launcher if it is absent.
If a default launcher already exists, inspect it and agree a safe promotion time;
do not overwrite it during an active run. Optional aliases
`microc2tracks-preflight`/`microc2tracks-merge` may point to the same launcher.
Retain version-specific commands and old installations for reproducibility and
rollback. No server-side validation is claimed by preparing this package.
