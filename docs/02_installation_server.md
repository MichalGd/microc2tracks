# Shared Server Installation

These notes target a shared Linux server with approximately 500 GB RAM, 70 physical CPU cores, 140 logical CPUs, and no/limited GPU.

## Recommended Directory Layout

```text
/opt/microc2tracks/                # cloned pipeline repo, readable by all users
/opt/conda/envs/microc2tracks/     # shared conda environment
/shared/references/
  mm39/
    GRCm39.primary_assembly.genome.fa
    GRCm39.primary_assembly.genome.fa.fai
    mm39.chrom.sizes
    bwa-mem2-index-files
/shared/software/juicer_tools/
  juicer_tools.2.18.00.jar
/data/projects/
  project_name/
    raw_fastq/
    results/
    work/
```

## Install Environment

Use a shared conda prefix, not a user-private environment in one person's home directory. This is the most important point for all-user access.

Recommended one-time setup by an admin or by a user with permission to write under `/opt`:

```bash
sudo mkdir -p /opt
sudo git clone https://github.com/MichalGd/microc2tracks.git /opt/microc2tracks
cd /opt/microc2tracks

# For future updates:
# cd /opt/microc2tracks
# sudo git pull

# Load conda from wherever it is installed on this server.
# This is safer than assuming /opt/miniconda3.
source "$(conda info --base)/etc/profile.d/conda.sh"

# Recommended before creating bioinformatics environments with conda.
conda config --set channel_priority strict
conda config --set solver libmamba

# Create the environment at an explicit shared path.
conda env create \
  -p /opt/conda/envs/microc2tracks \
  -f environment.yml

# Activate by path, not by user-local environment name.
conda activate /opt/conda/envs/microc2tracks
```

## Tested `biolserv` Installation Recipe

This is the known-good recipe from the first successful installation on `biolserv`. It assumes:

- repository path: `/opt/microc2tracks`
- shared conda environment path: `/opt/conda/envs/microc2tracks`
- conda base path: `/home/micgdu/miniconda3`
- no dedicated `bioinfo` group exists
- `micgdu` maintains the installation

Clone or update the repository:

```bash
sudo mkdir -p /opt
sudo git clone https://github.com/MichalGd/microc2tracks.git /opt/microc2tracks
cd /opt/microc2tracks

# For later updates:
# cd /opt/microc2tracks
# sudo git pull
```

Prepare the shared conda environment directory:

```bash
sudo mkdir -p /opt/conda/envs
sudo chown -R micgdu:micgdu /opt/conda
sudo chmod -R a+rX /opt/conda
sudo chmod -R u+rwX /opt/conda
```

Load conda and use the faster conda solver:

```bash
source /home/micgdu/miniconda3/etc/profile.d/conda.sh

conda config --set channel_priority strict
conda config --set solver libmamba
```

Create and activate the shared environment:

```bash
cd /opt/microc2tracks

conda env create \
  -p /opt/conda/envs/microc2tracks \
  -f environment.yml

conda activate /opt/conda/envs/microc2tracks
```

Verify the main tools and scripts:

```bash
which fastp
which bwa-mem2
which pairtools
which cooler
which multiqc

bash tests/check_bash_syntax.sh
bash scripts/microc2tracks.sh -h
```

Expected tool paths should look like:

```text
/opt/conda/envs/microc2tracks/bin/fastp
/opt/conda/envs/microc2tracks/bin/bwa-mem2
/opt/conda/envs/microc2tracks/bin/pairtools
/opt/conda/envs/microc2tracks/bin/cooler
/opt/conda/envs/microc2tracks/bin/multiqc
```

Make the workflow and environment readable/executable for all users, while preventing accidental edits to the environment:

```bash
sudo chmod -R a+rX /opt/conda/envs/microc2tracks
sudo chmod -R go-w /opt/conda/envs/microc2tracks
sudo chmod -R a+rX /opt/microc2tracks
```

Because this server uses conda from `/home/micgdu/miniconda3`, other users also need read/execute access to that conda installation:

```bash
chmod -R a+rX /home/micgdu/miniconda3
```

Create system-wide launcher commands:

```bash
sudo tee /usr/local/bin/microc2tracks >/dev/null <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
source /home/micgdu/miniconda3/etc/profile.d/conda.sh
conda activate /opt/conda/envs/microc2tracks
exec bash /opt/microc2tracks/scripts/microc2tracks.sh "$@"
EOF

sudo tee /usr/local/bin/microc2tracks-preflight >/dev/null <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
source /home/micgdu/miniconda3/etc/profile.d/conda.sh
conda activate /opt/conda/envs/microc2tracks
exec bash /opt/microc2tracks/scripts/preflight_check.sh "$@"
EOF

sudo chmod 755 /usr/local/bin/microc2tracks /usr/local/bin/microc2tracks-preflight
```

Test the launchers:

```bash
which microc2tracks
which microc2tracks-preflight

microc2tracks -h
microc2tracks-preflight -h
```

From another user account, run:

```bash
which microc2tracks
microc2tracks -h
microc2tracks-preflight -h
```

If those commands work from another account, the shared installation is ready.

For real project use:

```bash
cd /data/projects/my_project
cp /opt/microc2tracks/config/config_template.conf config.conf
cp /opt/microc2tracks/config/samplesheet_template.csv samplesheet.csv

# Edit config.conf and samplesheet.csv first.
microc2tracks-preflight -c config.conf -s samplesheet.csv
microc2tracks -c config.conf -s samplesheet.csv
```

If `environment.yml` is too large to solve cleanly, use a two-environment strategy. Conda can do this directly; it may simply be slower than mamba:

```bash
# Upstream matrix environment
conda create -p /opt/conda/envs/microc2tracks-core -c conda-forge -c bioconda \
  python fastp bwa-mem2 samtools htslib pairtools pairix cooler multiqc openjdk \
  pandas numpy matplotlib bioframe

# Downstream environment
conda env create \
  -p /opt/conda/envs/microc2tracks-downstream \
  -f envs/downstream_optional.yml
```

If `mamba` is available, prefer it for the optional downstream environment:

```bash
mamba env create \
  -p /opt/conda/envs/microc2tracks-downstream \
  -f envs/downstream_optional.yml
```

The downstream YAML installs Mustache through `pip` because `mustache-hic` is
not a conda package on the standard `conda-forge`/`bioconda` channels. It also
pins `setuptools<81` for `coolpuppy`/`h5sparse` compatibility.

## Permissions

The pipeline repo and shared environment should be readable/executable by all users. Only admins or the bioinformatics-maintainer group should need write access.

```bash
sudo chgrp -R bioinfo /opt/microc2tracks /opt/conda/envs/microc2tracks
sudo chmod -R a+rX /opt/microc2tracks /opt/conda/envs/microc2tracks
sudo chmod -R a+rX /shared/references /shared/software/juicer_tools
```

If the `bioinfo` group should maintain the installation, allow group write on the repository but keep the environment itself read-only for normal users:

```bash
sudo chmod -R g+rwX /opt/microc2tracks
sudo find /opt/microc2tracks -type d -exec chmod g+s {} \;
```

Project-specific `results/` and `work/` directories should normally be group-writable:

```bash
mkdir -p /data/projects/my_project/results /data/projects/my_project/work
chgrp -R bioinfo /data/projects/my_project
chmod -R g+rwX /data/projects/my_project
```

## Optional System-Wide Launchers

To let users run the workflow without typing long paths, create small wrapper scripts in `/usr/local/bin`.

```bash
sudo tee /usr/local/bin/microc2tracks >/dev/null <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
CONDA_BASE="$(conda info --base)"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate /opt/conda/envs/microc2tracks
exec bash /opt/microc2tracks/scripts/microc2tracks.sh "$@"
EOF

sudo tee /usr/local/bin/microc2tracks-preflight >/dev/null <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
CONDA_BASE="$(conda info --base)"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate /opt/conda/envs/microc2tracks
exec bash /opt/microc2tracks/scripts/preflight_check.sh "$@"
EOF

sudo chmod 755 /usr/local/bin/microc2tracks /usr/local/bin/microc2tracks-preflight
```

Then any user can run:

```bash
microc2tracks-preflight -c config.conf -s samplesheet.csv
microc2tracks -c config.conf -s samplesheet.csv
```

## Prepare References

```bash
bash scripts/prepare_reference.sh \
  -f /shared/references/mm39/GRCm39.primary_assembly.genome.fa \
  -a mm39 \
  -o /shared/references/mm39

bash scripts/prepare_reference.sh \
  -f /shared/references/hg38/GRCh38.primary_assembly.genome.fa \
  -a hg38 \
  -o /shared/references/hg38
```

Copy `config/references.tsv` to a writable site configuration, edit the mm39 and
hg38 paths, and keep the exact UCSC canonical names described in
`09_multi_reference.md`. Preflight checks only references used by the current
sample sheet, so a mouse-only run does not require hg38 to be installed.

Always confirm that the chromosome names in the pair files match the chromosome sizes file. A mismatch can create an apparently valid but empty cooler.

## Project Setup

```bash
cd /data/projects/my_project
cp /opt/microc2tracks/config/config_template.conf config.conf
cp /opt/microc2tracks/config/samplesheet_template.csv samplesheet.csv
```

Edit `config.conf` so `OUTDIR`, `TMPDIR`, `REFERENCE_REGISTRY`, and
`JUICER_TOOLS_JAR` match the server. Existing mouse-only global reference values
may remain as default-mm39 overrides.

If the BWA-MEM2 index was created with the same FASTA path, keep:

```bash
BWA_INDEX_PREFIX="$REFERENCE_FASTA"
```

If the index was created with another prefix, set `BWA_INDEX_PREFIX` explicitly.

Then run:

```bash
bash /opt/microc2tracks/scripts/preflight_check.sh -c config.conf -s samplesheet.csv
bash /opt/microc2tracks/scripts/microc2tracks.sh -c config.conf -s samplesheet.csv
```

## Resource Profiles

Small test dataset:

- alignment: 4 to 8 threads
- pairtools sort: 4 threads, 8 GB
- matrix: 4 to 8 threads
- Java heap: 32 GB

Typical single mammalian sample:

- alignment: 24 to 32 threads
- pairtools sort: 12 to 16 threads, 24 to 32 GB
- matrix: 16 to 32 threads
- Java heap: 128 to 256 GB

Large merged Micro-C dataset:

- merge: 48 to 64 threads, 128 GB
- matrix: 32 to 64 threads
- Java heap: 256 GB

On a shared server, avoid using all 140 logical CPUs unless this has been coordinated with other users.
