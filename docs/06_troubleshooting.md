# Troubleshooting

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

Install the core matrix environment first, then install downstream tools separately.

Core:

```bash
conda create -n microc2tracks-core -c conda-forge -c bioconda \
  python fastp bwa-mem2 samtools htslib pairtools pairix cooler multiqc openjdk
```

Downstream:

```bash
conda create -n microc2tracks-downstream -c conda-forge -c bioconda \
  python cooltools hictk mustache-hic chromosight coolpuppy hicexplorer
```

## Micro-C And Hi-C Look Different At The Same Resolution

This is expected. Compare matrices at matched resolution and matched valid-contact depth when making quantitative claims. Use observed/expected or distance-aware summaries for cross-assay comparisons.
