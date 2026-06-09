# Troubleshooting

## `fastp` Says `Failed to open file` But `ls` Shows The FASTQ Exists

Most likely cause: the sample sheet has hidden Windows CRLF line endings, so the last CSV field may be passed to Bash with an invisible carriage return. The file name looks correct on screen, but `fastp` receives a different path.

Check for hidden carriage returns:

```bash
cd /path/to/project
sed -n '1,5l' config/samplesheet.csv
```

If lines end with `\r$`, current versions of `microc2tracks-preflight` and `microc2tracks` fix this automatically at startup. You can also convert the sample sheet manually:

```bash
sed -i 's/\r$//' config/samplesheet.csv
```

Then rerun preflight and the workflow:

```bash
microc2tracks-preflight -c config/config.conf -s config/samplesheet.csv
microc2tracks -c config/config.conf -s config/samplesheet.csv
```

Current versions also normalize the config and sample sheet in place before parsing. The automatic cleanup removes CRLF/CR line endings, UTF-8 byte-order marks, non-breaking spaces, copied smart quotes, and surrounding whitespace around sample-sheet fields.

## `fastp` Appears Stalled On One Sample

If earlier samples finished and the next sample sits at `running fastp` with almost no CPU, suspect that sample's input FASTQ files or the storage path first. A previously downloaded `.hic` file is not used by `fastp`; it only matters if it filled the output/tmp filesystem or stressed the same storage mount.

Check disk space and whether the trimmed outputs are still growing:

```bash
df -h "$OUTDIR" "$TMPDIR" /dysk3
ls -lh "$OUTDIR"/PMKbas_MC1/02_trimmed "$OUTDIR"/PMKbas_MC1/logs/PMKbas_MC1.fastp.log
ps -o pid,stat,etime,pcpu,pmem,cmd -p FASTP_PID
```

Check the suspicious FASTQs directly:

```bash
ls -lh /path/to/PMKbas_MC1/*_1.fq.gz /path/to/PMKbas_MC1/*_2.fq.gz
gzip -t /path/to/PMKbas_MC1/*_1.fq.gz
gzip -t /path/to/PMKbas_MC1/*_2.fq.gz
```

For a stronger preflight before rerunning, enable full gzip checks and an optional `fastp` timeout in `config.conf`:

```bash
CHECK_FASTQ_GZIP="true"
FASTP_TIMEOUT_SECONDS="14400"
```

Then rerun:

```bash
bash scripts/preflight_check.sh -c config/config.conf -s config/samplesheet.csv
bash scripts/microc2tracks.sh -c config/config.conf -s config/samplesheet.csv
```

The main runner logs each input FASTQ size/mtime before `fastp`; compare the stalled sample against samples that already finished.

## Empty Cooler: `nnz: 0`, `sum: 0`

Most likely cause: chromosome names in `.pairs.gz` do not match `CHROM_SIZES`.

Check:

```bash
zcat sample.valid.mapq30.pairs.gz | awk 'BEGIN{FS="\t"} $1 !~ /^#/ {print $2; exit}'
head "$CHROM_SIZES"
```

Fix by regenerating `CHROM_SIZES` from the exact FASTA used for alignment:

```bash
samtools faidx reference.fa
cut -f1,2 reference.fa.fai > genome.chrom.sizes
```

## Juicer Tools Fails Or Runs Out Of Memory

Increase `JAVA_HEAP`, for example:

```bash
JAVA_HEAP="256g"
```

For very large datasets, build `.mcool` first, then run `.hic` generation separately when the server is quiet.

If sample-level parallelism is enabled, keep `.hic` work serialized:

```bash
MAX_PARALLEL_SAMPLES="2"
MAX_PARALLEL_HIC="1"
```

## Merged `.hic` Normalization Takes Too Long

If the main log appears to sit near matrix balancing but `top` shows a Java
process, check the Juicer logs:

```bash
ps -fp JAVA_PID
tail -n 40 results/merged/GROUP/logs/GROUP.juicer_addNorm.log
```

`juicer_tools addNorm` can be the long step for deep merged Micro-C datasets,
especially at 1000 or 2000 bp. Use separate `.hic` resolutions or fewer
normalizations for exploratory runs:

```bash
MCOOL_RESOLUTIONS="$RESOLUTIONS"
HIC_RESOLUTIONS="5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000"
THREADS_HIC_NORM="24"
HIC_NORMALIZATIONS="VC,VC_SQRT,KR,SCALE"
```

Set `RUN_HIC="false"` to skip `.hic` creation entirely while keeping `.mcool`
outputs for open2c analyses.

## Pairtools Sort Is Slow

Use a fast local temporary directory:

```bash
TMPDIR="/scratch/$USER/microc2tracks_tmp"
```

Increase `THREADS_SORT` and `PAIRTOOLS_SORT_MEMORY` if RAM is available.

## Temporary Disk Space Fills Up

Large Micro-C/Hi-C samples create large pairtools sort and merge temporary files. Put `TMPDIR` on fast, spacious scratch storage:

```bash
TMPDIR="/scratch/$USER/microc2tracks_tmp"
```

Keep cleanup enabled:

```bash
CLEAN_TMP_ON_SUCCESS="true"
KEEP_TRIMMED_FASTQ="false"
KEEP_JUICER_PAIRS="false"
```

For very large samples, avoid increasing `MAX_PARALLEL_SAMPLES` until one complete serial run has been checked.

## Conda Environment Does Not Solve

First try conda's faster solver:

```bash
conda config --set channel_priority strict
conda config --set solver libmamba
```

Then retry `conda env create`. This still uses `conda`; it only changes the solver.

If it still does not solve, install the core matrix environment first, then install downstream tools separately.

Core:

```bash
conda create -p /opt/conda/envs/microc2tracks-core -c conda-forge -c bioconda \
  python fastp bwa-mem2 samtools htslib pairtools pairix cooler multiqc openjdk
```

Downstream:

```bash
conda create -p /opt/conda/envs/microc2tracks-downstream -c conda-forge -c bioconda \
  python cooltools hictk mustache-hic chromosight coolpuppy hicexplorer
```

## Micro-C And Hi-C Look Different At The Same Resolution

This is expected. Compare matrices at matched resolution and matched valid-contact depth when making quantitative claims. Use observed/expected or distance-aware summaries for cross-assay comparisons.
