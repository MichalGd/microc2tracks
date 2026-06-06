# Shared Server Installation

These notes target a shared Linux server with approximately 500 GB RAM, 70 physical CPU cores, 140 logical CPUs, and no/limited GPU.

## Recommended Directory Layout

```text
/opt/microc2tracks/                # cloned pipeline repo, readable by all users
/opt/conda/envs/microc2tracks/     # shared conda/mamba environment
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

```bash
cd /opt/microc2tracks
mamba env create -f environment.yml
conda activate microc2tracks
```

If `environment.yml` is too large to solve cleanly, use a two-environment strategy:

```bash
# Upstream matrix environment
mamba create -n microc2tracks-core -c conda-forge -c bioconda \
  python fastp bwa-mem2 samtools htslib pairtools pairix cooler multiqc openjdk \
  pandas numpy matplotlib bioframe

# Downstream environment
mamba create -n microc2tracks-downstream -c conda-forge -c bioconda \
  python cooltools hictk mustache-hic chromosight coolpuppy hicexplorer \
  pandas numpy matplotlib bioframe
```

The optional downstream environment can also be created from the repository:

```bash
mamba env create -f envs/downstream_optional.yml
```

## Permissions

The pipeline repo and shared environment should be readable/executable by all users:

```bash
sudo chgrp -R bioinfo /opt/microc2tracks /opt/conda/envs/microc2tracks
sudo chmod -R g+rX /opt/microc2tracks /opt/conda/envs/microc2tracks
sudo chmod -R a+rX /shared/references /shared/software/juicer_tools
```

Project-specific `results/` and `work/` directories should normally be group-writable:

```bash
mkdir -p /data/projects/my_project/results /data/projects/my_project/work
chgrp -R bioinfo /data/projects/my_project
chmod -R g+rwX /data/projects/my_project
```

## Prepare References

```bash
bash scripts/prepare_reference.sh \
  -f /shared/references/mm39/GRCm39.primary_assembly.genome.fa \
  -a mm39 \
  -o /shared/references/mm39
```

Always confirm that the chromosome names in the pair files match the chromosome sizes file. A mismatch can create an apparently valid but empty cooler.

## Project Setup

```bash
cd /data/projects/my_project
cp /opt/microc2tracks/config/config_template.conf config.conf
cp /opt/microc2tracks/config/samplesheet_template.csv samplesheet.csv
```

Edit `config.conf` so `OUTDIR`, `TMPDIR`, `REFERENCE_FASTA`, `CHROM_SIZES`, and `JUICER_TOOLS_JAR` match the server.

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
