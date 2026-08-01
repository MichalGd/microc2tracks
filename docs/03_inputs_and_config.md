# Inputs And Configuration

## Sample Sheet

Required columns:

```csv
sample,assay,reference_genome,condition,biological_replicate,technical_replicate,fastq_r1,fastq_r2
```

Column meanings:

- `sample`: unique technical replicate ID, used in output file names
- `assay`: `microc` or `hic`
- `reference_genome`: controlled alias `mouse`, `mm39`, `human`, or `hg38`
- `condition`: free text metadata, for example `basal` or `treated`
- `biological_replicate`: biological replicate number within an assay and condition
- `technical_replicate`: technical replicate number within one biological replicate; use `1` when there is no technical replicate
- `fastq_r1`: absolute or project-relative path to R1 FASTQ
- `fastq_r2`: absolute or project-relative path to R2 FASTQ

Each row is processed independently first and produces its own `.mcool` and `.hic` outputs. After all rows finish, rows with the same reference, `assay`, `condition`, and `biological_replicate` are treated as technical replicates. Different references are never merged. An optional safe `merge_group` column can name a group explicitly.

Example:

```csv
sample,assay,reference_genome,condition,biological_replicate,technical_replicate,fastq_r1,fastq_r2
microC_B1_T1,microc,mouse,basal,1,1,/data/project/PMK_Basal_MicroC_B1_T1_R1.fq.gz,/data/project/PMK_Basal_MicroC_B1_T1_R2.fq.gz
microC_B1_T2,microc,mm39,basal,1,2,/data/project/PMK_Basal_MicroC_B1_T2_R1.fq.gz,/data/project/PMK_Basal_MicroC_B1_T2_R2.fq.gz
human_B1_T1,microc,human,human_control,1,1,/data/project/Human_MicroC_B1_T1_R1.fq.gz,/data/project/Human_MicroC_B1_T1_R2.fq.gz
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

`REFERENCE_REGISTRY` points to the centralized TSV. Each reference row supplies
its FASTA, BWA-MEM2 index prefix, chromosome sizes, canonical regex, browser
preset, phasing track, and optional annotation/blacklist metadata. The FASTA
`.fai` and chromosome sizes must contain exactly matching names and lengths.

The legacy `REFERENCE_FASTA`, `BWA_INDEX_PREFIX`, `CHROM_SIZES`,
`CANONICAL_CHROMS_REGEX`, and `PHASING_TRACK` values override only the default
reference row, preserving old mouse configs. See `09_multi_reference.md`.

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

Within each registry row, `bwa_index_prefix` tells the pipeline where the
BWA-MEM2 index is. It usually equals that row's FASTA:

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
THREADS_FASTP="8"
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
Sentinels record canonical reference identity. A mismatched or legacy human
sentinel is never reused. Missing sentinels can be bootstrapped only for legacy
mm39 outputs by default:

```bash
BOOTSTRAP_SENTINELS="true"
ALLOW_LEGACY_MM39_RESUME="true"
```

Changing reference invalidates reuse. Prefer a new output directory; otherwise
archive and remove the affected sample and merged-group outputs before rerunning.

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

By default, the workflow processes two sample-sheet rows at a time while
serializing matrix and `.hic` work:

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
# Registry values:
# mm39 ^chr([1-9]|1[0-9]|X|Y|M)$
# hg38 ^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$
```

This removes alt, random, unplaced, and other non-canonical contigs from the matrix and `.hic` outputs. This is recommended for UCSC Genome Browser compatibility.

Important: this filter removes non-canonical contigs but does not rename chromosomes. The shipped mm39 and hg38 definitions require UCSC-style `chr1` names and deliberate `chrM`. Names such as `1`, `MT`, `chrMT`, or `NC_000001.11` require a distinct explicit registry/custom-browser design and are not silently rewritten. `chrY` stays in both matrix schemas, including for expected XX samples.

## Assay Behavior

For both assays the current script uses:

```text
fastp -> bwa-mem2 -> pairtools -> cooler -> Juicer Tools
```

For Micro-C, the valid-pair filter also removes cis contacts closer than `MICROC_MIN_CIS_DIST`.

For Hi-C, the first version uses the shared path and `HIC_MIN_CIS_DIST`. A future restriction-fragment module should add enzyme-specific fragment annotation and filtering when reliable restriction enzyme metadata are available.
