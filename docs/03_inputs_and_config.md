# Inputs And Configuration

## Sample Sheet

Required columns:

```csv
sample,assay,condition,biological_replicate,technical_replicate,fastq_r1,fastq_r2
```

Column meanings:

- `sample`: unique technical replicate ID, used in output file names
- `assay`: `microc` or `hic`
- `condition`: free text metadata, for example `basal` or `treated`
- `biological_replicate`: biological replicate number within an assay and condition
- `technical_replicate`: technical replicate number within one biological replicate; use `1` when there is no technical replicate
- `fastq_r1`: absolute or project-relative path to R1 FASTQ
- `fastq_r2`: absolute or project-relative path to R2 FASTQ

Each row is processed independently first and produces its own `.mcool` and `.hic` outputs. After all rows finish, rows with the same `assay`, `condition`, and `biological_replicate` are treated as technical replicates. If a group has two or more technical replicates, the workflow merges the filtered pair files and rebuilds merged `.mcool` and `.hic` outputs.

Example:

```csv
sample,assay,condition,biological_replicate,technical_replicate,fastq_r1,fastq_r2
microC_B1_T1,microc,basal,1,1,/data/project/PMK_Basal_MicroC_B1_T1_R1.fq.gz,/data/project/PMK_Basal_MicroC_B1_T1_R2.fq.gz
microC_B1_T2,microc,basal,1,2,/data/project/PMK_Basal_MicroC_B1_T2_R1.fq.gz,/data/project/PMK_Basal_MicroC_B1_T2_R2.fq.gz
microC_B2_T1,microc,basal,2,1,/data/project/PMK_Basal_MicroC_B2_T1_R1.fq.gz,/data/project/PMK_Basal_MicroC_B2_T1_R2.fq.gz
```

In this example, `microC_B1_T1` and `microC_B1_T2` are processed separately and then merged as `basal_microc_B1_tech_merged`. `microC_B2_T1` is processed separately and is not copied into a redundant merged output because it has only one technical replicate.

Validate the file before running:

```bash
python scripts/validate_samplesheet.py config/samplesheet.csv
```

The preflight check and main runner first run `scripts/sanitize_text_inputs.py` on `config.conf` and the sample sheet. This automatically fixes common Windows/copy-paste artifacts before Bash or CSV parsing:

- CRLF or stray carriage-return line endings
- UTF-8 byte-order marks
- non-breaking spaces
- copied smart quotes
- surrounding whitespace around sample-sheet fields

After this cleanup, preflight and the main runner use stricter validation and require FASTQ files to exist and be readable before the run starts.

## Main Config Values

`REFERENCE_FASTA` and `CHROM_SIZES` must use matching chromosome names.

Good:

```text
pair file: chr1
chrom sizes: chr1
```

Bad:

```text
pair file: chr1
chrom sizes: mchr1
```

`BASE_RESOLUTION` controls the first `.cool` file. For deep Micro-C, 1000 is practical. For lower-depth Hi-C, 5000 or 10000 may be safer.

`BWA_INDEX_PREFIX` tells the pipeline where the BWA-MEM2 index is. In the usual case, it is the same as `REFERENCE_FASTA`:

```bash
REFERENCE_FASTA="/shared/references/mm39/GRCm39.primary_assembly.genome.fa"
BWA_INDEX_PREFIX="$REFERENCE_FASTA"
```

This assumes the index was created with:

```bash
bwa-mem2 index /shared/references/mm39/GRCm39.primary_assembly.genome.fa
```

If the index was created with a different prefix, set `BWA_INDEX_PREFIX` to that prefix.

`RESOLUTIONS` is the shared fallback resolution list. By default, `.mcool` and
Juicer `.hic` use the same values:

```bash
RESOLUTIONS="1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000"
MCOOL_RESOLUTIONS="$RESOLUTIONS"
HIC_RESOLUTIONS="$RESOLUTIONS"
```

For large merged Micro-C datasets, keep `.mcool` high resolution and make `.hic`
coarser if Juicer normalization is too slow:

```bash
MCOOL_RESOLUTIONS="$RESOLUTIONS"
HIC_RESOLUTIONS="5000,10000,25000,50000,100000,250000,500000,1000000,2500000,5000000,10000000"
```

`MICROC_MIN_CIS_DIST=1000` removes very short-range cis contacts in Micro-C mode. Do not blindly apply this to Hi-C if the goal is restriction-fragment-aware processing.

`RUN_MERGE_TECHNICAL_REPLICATES=true` enables automatic post-processing merges for sample-sheet rows that share `assay`, `condition`, and `biological_replicate`. Set it to `false` if you want only per-row outputs and plan to merge manually later.

For suspicious FASTQ inputs, preflight can summarize file sizes and optionally read every `.gz` file fully:

```bash
CHECK_FASTQ_GZIP="true"
```

The main runner can also put a wall-clock guard around `fastp`:

```bash
FASTP_TIMEOUT_SECONDS="14400"
```

The default `FASTP_TIMEOUT_SECONDS="0"` disables the timeout. Use a large value rather than a tight one, because real Micro-C FASTQs can be large.

Thread settings are split by step:

```bash
THREADS_FASTP="12"
THREADS_ALIGN="32"
THREADS_SORT="16"
THREADS_MATRIX="32"
THREADS_HIC_NORM="24"
```

`THREADS_FASTP` is intentionally lower than `THREADS_ALIGN`; compressed FASTQ
I/O and `fastp` can stall or scale poorly at very high thread counts. BWA-MEM2
can still use more threads.

Juicer normalizations can be selected with:

```bash
HIC_NORMALIZATIONS="VC,VC_SQRT,KR,SCALE"
```

For exploratory runs, fewer `.hic` normalizations or coarser `HIC_RESOLUTIONS`
can save substantial time.

The runner writes per-step sentinel files under each `logs/done/` directory.
If valid outputs already exist, missing sentinels are bootstrapped by default:

```bash
BOOTSTRAP_SENTINELS="true"
```

This improves resumability without forcing old completed runs to recompute.

Light preliminary downstream analysis is enabled by default:

```bash
RUN_PRELIM_DOWNSTREAM="true"
```

The main workflow calls `run_downstream.sh -l` from each technical-replicate
balanced `.mcool`, and the merge workflow does the same for each merged
technical-replicate `.mcool`. Light mode runs expected contacts, insulation
scores/boundaries, TAD-like BED intervals, and compartment eigenvectors, and
skips heavier loop callers. Set this to `false` to finish matrix production
first and run downstream analyses later.

The default TAD/domain calls are derived from cooltools insulation boundaries:

```bash
TAD_RESOLUTION="10000"
INSULATION_WINDOWS_BP="50000 100000 250000"
INSULATION_BEDGRAPH_TRACKLINE="true"
TAD_BOUNDARY_THRESHOLD="Li"
TAD_BOUNDARY_WINDOW_BP="100000"
TAD_MIN_SIZE_BP="40000"
TAD_MIN_BOUNDARY_STRENGTH=""
```

`TAD_BOUNDARY_WINDOW_BP` should be one of the values in
`INSULATION_WINDOWS_BP`. The TAD BED file contains intervals between adjacent
strong boundaries at that window size.

The pipeline also exports one log2 insulation score bedGraph per insulation
window under `05_downstream/insulation/bedgraph/`. With
`INSULATION_BEDGRAPH_TRACKLINE="true"`, each file starts with a UCSC
`track type=bedGraph` line for direct browser upload. Set it to `false` if a
strict downstream parser needs pure four-column bedGraph.

Compartment eigenvectors are calculated even when `PHASING_TRACK` is empty. In
that case the PC1 sign is arbitrary: positive PC1 does not automatically mean A
or B compartment. To orient PC1 during the run, provide a BED-like phasing track
with a biological signal such as GC content, gene density, ATAC-seq, H3K27ac, or
RNA-seq.

See `docs/08_light_downstream_analysis.md` for output interpretation and links
to the tool documentation.

## Parallel Execution

By default, the workflow processes one sample-sheet row at a time:

```bash
MAX_PARALLEL_SAMPLES="1"
MAX_PARALLEL_MATRIX="1"
MAX_PARALLEL_HIC="1"
```

On a quiet large-memory server, a practical first parallel setting is:

```bash
MAX_PARALLEL_SAMPLES="2"
MAX_PARALLEL_MATRIX="1"
MAX_PARALLEL_HIC="1"
```

This allows two sample workers to run through FASTQ/alignment/pairs work while keeping cooler/mcool and Juicer `.hic` work serialized. This is safer than running two large `juicer_tools addNorm` jobs at the same time. Values of `3` or `4` are possible, but increase disk I/O, temporary storage pressure, and RAM risk.

## Retention Policy

Raw FASTQ files are never deleted by the pipeline because they live outside `OUTDIR` and are the primary data. Reproducible intermediates can be controlled with:

```bash
KEEP_TRIMMED_FASTQ="false"
KEEP_DEDUP_PAIRS="true"
KEEP_JUICER_PAIRS="false"
KEEP_SINGLE_RES_COOL="false"
KEEP_RAW_MCOOL="false"
KEEP_RAW_HIC="false"
KEEP_RAW_MERGED_PAIRS="false"
CLEAN_TMP_ON_SUCCESS="true"
```

Recommended during method development: keep `KEEP_DEDUP_PAIRS="true"` so you can re-run MAPQ or distance filters without realigning. After the workflow is stable and storage is tight, switch it to `false` and keep raw FASTQ, valid pairs, final `.norm.mcool`, final `.norm.hic`, reports, config, and logs.

## Final Report

These switches control end-of-run reporting:

```bash
RUN_GLOBAL_MULTIQC="true"
RUN_FINAL_REPORT="true"
```

The final report combines parsed `fastp` JSON, `pairtools` stats, cooler info JSON, sample status files, merge metadata, and output file sizes into TSV/HTML summaries under `results/final_report/`.

## Canonical Chromosome Filtering

By default, the workflow filters valid pairs and matrix chromosome sizes to canonical chromosomes before creating `.cool`, `.mcool`, and `.hic` files:

```bash
FILTER_CANONICAL_CHROMS="true"
CANONICAL_CHROMS_REGEX="^(chr)?([1-9][0-9]?|X|Y|M|MT)$"
```

This removes alt, random, unplaced, and other non-canonical contigs from the matrix and `.hic` outputs. This is recommended for UCSC Genome Browser compatibility.

Important: this filter removes non-canonical contigs but does not rename chromosomes. For UCSC `mm39`, use a reference and chromosome sizes file with UCSC-style names such as `chr1`, `chr2`, and `chrX`. Ensembl-style names such as `1`, `2`, and `X` may still be incompatible with the UCSC `mm39` browser even after non-canonical contigs are removed.

## Assay Behavior

For both assays the current script uses:

```text
fastp -> bwa-mem2 -> pairtools -> cooler -> Juicer Tools
```

For Micro-C, the valid-pair filter also removes cis contacts closer than `MICROC_MIN_CIS_DIST`.

For Hi-C, the first version uses the shared path and `HIC_MIN_CIS_DIST`. A future restriction-fragment module should add enzyme-specific fragment annotation and filtering when reliable restriction enzyme metadata are available.
