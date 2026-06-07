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

## Pairtools Sort Is Slow

Use a fast local temporary directory:

```bash
TMPDIR="/scratch/$USER/microc2tracks_tmp"
```

Increase `THREADS_SORT` and `PAIRTOOLS_SORT_MEMORY` if RAM is available.

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
