# Analysis And Design

## Legacy Script Summary

The legacy Micro-C script is a successful manual prototype. It downloads vendor archives, checks MD5 sums, unpacks FASTQ files, prepares an mm39 reference, trims reads with `fastp`, aligns with `bwa-mem2 mem -SP5M`, parses alignments with `pairtools parse`, sorts and deduplicates contacts with `pairtools sort` and `pairtools dedup`, filters to valid MAPQ30 pairs, builds `.cool` and `.mcool` matrices with `cooler`, and builds `.hic` files with Juicer Tools.

It also merges four Micro-C replicates with `pairtools merge` and rebuilds matrices from the merged pair file. The later notes in the script document an important failure mode: chromosome names in the pair file and chromosome sizes file did not match, producing an empty cooler with `nnz: 0` and `sum: 0`. This is now treated as a first-class validation issue.

## Tools Already Present In The Existing `microc` Conda Environment

The current environment already contains the upstream matrix backbone:

- `fastp`
- `bwa-mem2`
- `samtools`
- `pairtools`
- `pairix`
- `bgzip` via `htslib`
- `cooler`
- `multiqc`
- `openjdk`
- Python packages such as `numpy`, `pandas`, `matplotlib`, and `bioframe`

Main missing or optional downstream tools:

- `cooltools`
- `hictk`
- `mustache-hic`
- `chromosight`
- `coolpuppy`
- `hicexplorer`

For stability, keep the upstream FASTQ-to-matrix environment small if conda solving becomes difficult. Downstream tools can live in a second environment that is activated after `.mcool` files exist.

## Main Limitations In The Legacy Script

- Hard-coded paths to raw data, references, Juicer Tools, and web export folders.
- Repeated code blocks for each replicate.
- No sample sheet or structured config.
- Sample naming mistakes are easy, for example repeated `SAMPLE=microC3`.
- No automated preflight check for reference/FASTA/chromosome-size consistency.
- No consistent output tree.
- `.dedup.all.pairs.gz` is a confusing name because the file is already deduplicated.
- No reusable downstream modules for TADs, compartments, loops, stripes, saddle plots, or comparisons.
- No shared-server installation notes or permission model.

## Method Review Synthesis

The review notes strongly support one unified workflow with assay-specific configuration rather than two independent pipelines.

Shared for Micro-C and Hi-C:

- read QC and trimming
- `bwa mem` or `bwa-mem2 mem` with Hi-C style chimeric-read flags
- pairs parsing, sorting, deduplication, and stats
- `.pairs.gz` as the central intermediate
- `.cool` and `.mcool` generation
- `.hic` generation
- balancing and normalization
- downstream analyses from contact matrices

Micro-C-specific:

- MNase/no-restriction mode
- stronger attention to short insert and near-diagonal contacts
- default near-diagonal cis filter of 1000 bp
- finer matrix resolutions, often 1 to 5 kb when sequencing depth allows
- optional read stitching in future versions

Hi-C-specific:

- restriction enzyme metadata
- restriction fragment BED generation
- fragment-level filtering of dangling ends, self-circles, and intra-fragment artifacts
- usually coarser default downstream resolutions than deep Micro-C

## Recommended Architecture

This first version uses a practical shell layout:

```text
config/
docs/
examples/
scripts/
tests/
environment.yml
README.md
```

The central abstraction is a filtered, deduplicated, indexed `.pairs.gz` file. Once both Micro-C and Hi-C are converted into comparable pairs and matrices, downstream modules can operate on the same `.mcool` and `.hic` products.

The sample sheet separates biological and technical replicates. Each row is one technical replicate and is processed to its own `.mcool` and `.hic`; after all rows finish, rows with matching assay, condition, and biological replicate are merged as technical replicates and matrices are rebuilt from the merged pair file.

## Default Pipeline Path

```text
FASTQ
  -> fastp
  -> bwa-mem2 mem -SP5M
  -> pairtools parse
  -> pairtools sort
  -> pairtools dedup
  -> pairtools select valid MAPQ-filtered contacts
  -> pairix
  -> pairtools stats
  -> cooler cload pairix
  -> cooler balance
  -> cooler zoomify
  -> Juicer Tools pre/addNorm
  -> MultiQC
```

## CPU/RAM Strategy

The target server has 500 GB RAM, 70 physical cores, and 140 logical CPUs. Defaults intentionally reserve capacity for other users:

- alignment: 32 threads
- pair sorting: 16 threads and 32 GB memory per sample
- matrix generation: 32 threads
- replicate merge: 64 threads and 128 GB memory
- Juicer Tools: 256 GB Java heap for large mammalian matrices

Use more threads for a single large run only when the server is otherwise quiet. Pair sorting, replicate merging, matrix loading, balancing, and `.hic` generation are the main CPU/RAM bottlenecks.

## Downstream Defaults

- insulation/TADs: `cooltools insulation`
- compartments: `cooltools eigs-cis`
- saddle plots: `cooltools saddle` when phasing data are available
- loops: `cooltools dots` and/or Mustache
- loop anchors and pileups: future `coolpuppy` module
- stripes: future Chromosight stripe module
- side-by-side Micro-C vs Hi-C plots: `scripts/compare_matrices.py`

## First-Version Boundary

This repository implements the shared backbone and Micro-C-ready defaults now. Hi-C is supported through the same contact-pair path, but strict restriction-fragment filtering is marked as the main next module because it requires reliable enzyme and restriction fragment metadata per sample.
