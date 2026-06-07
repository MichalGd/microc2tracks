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

`RESOLUTIONS` controls the `.mcool` zoom levels and Juicer `.hic` resolutions.

`MICROC_MIN_CIS_DIST=1000` removes very short-range cis contacts in Micro-C mode. Do not blindly apply this to Hi-C if the goal is restriction-fragment-aware processing.

`RUN_MERGE_TECHNICAL_REPLICATES=true` enables automatic post-processing merges for sample-sheet rows that share `assay`, `condition`, and `biological_replicate`. Set it to `false` if you want only per-row outputs and plan to merge manually later.

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
