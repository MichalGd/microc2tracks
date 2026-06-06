# Inputs And Configuration

## Sample Sheet

Required columns:

```csv
sample,assay,replicate_group,condition,fastq_r1,fastq_r2
```

Column meanings:

- `sample`: unique sample ID, used in output file names
- `assay`: `microc` or `hic`
- `replicate_group`: biological group for later merging, for example `PMK_Basal`
- `condition`: free text metadata
- `fastq_r1`: absolute or project-relative path to R1 FASTQ
- `fastq_r2`: absolute or project-relative path to R2 FASTQ

Validate the file before running:

```bash
python scripts/validate_samplesheet.py config/samplesheet.csv
```

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

`RESOLUTIONS` controls the `.mcool` zoom levels and Juicer `.hic` resolutions.

`MICROC_MIN_CIS_DIST=1000` removes very short-range cis contacts in Micro-C mode. Do not blindly apply this to Hi-C if the goal is restriction-fragment-aware processing.

## Assay Behavior

For both assays the current script uses:

```text
fastp -> bwa-mem2 -> pairtools -> cooler -> Juicer Tools
```

For Micro-C, the valid-pair filter also removes cis contacts closer than `MICROC_MIN_CIS_DIST`.

For Hi-C, the first version uses the shared path and `HIC_MIN_CIS_DIST`. A future restriction-fragment module should add enzyme-specific fragment annotation and filtering when reliable restriction enzyme metadata are available.

