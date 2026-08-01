# Server Deployment Notes

This candidate is prepared for manual deployment only. No server was changed.

## Files To Deploy

Deploy the complete file set listed in `UPLOAD_MANIFEST.txt`. In particular,
keep `scripts/reference_registry.py`, `scripts/reference_sentinel.sh`,
`scripts/validate_pairs_reference.py`, `scripts/check_fastq_pair_names.py`, and
`config/references.tsv` together with the updated runner, merge, preflight,
reporting, downstream, and validator scripts.

## Reference Preparation

For each assembly, use one FASTA with UCSC-style chromosome names and derive the
`.fai`, chromosome sizes, and BWA-MEM2 index from that exact FASTA:

```bash
bash scripts/prepare_reference.sh \
  -f /shared/references/mm39/GRCm39.primary_assembly.genome.fa \
  -a mm39 -o /shared/references/mm39

bash scripts/prepare_reference.sh \
  -f /shared/references/hg38/GRCh38.primary_assembly.genome.fa \
  -a hg38 -o /shared/references/hg38
```

Edit the two registry rows to the deployed paths. Keep mm39's exact canonical
set at `chr1`-`chr19`, `chrX`, `chrY`, `chrM`; keep hg38's at `chr1`-`chr22`,
`chrX`, `chrY`, `chrM`. Optional phasing tracks must match their row's assembly.

## Configuration Migration

Add to existing configs:

```bash
REFERENCE_REGISTRY="/opt/microc2tracks/config/references.tsv"
DEFAULT_REFERENCE_ID="mm39"
THREADS_FASTP="8"
MAX_PARALLEL_SAMPLES="2"
CHECK_TRIMMED_FASTQ_SYNC="false"
```

Old mouse global reference values may remain; they override the default mm39
row. Add `reference_genome` to new sample sheets. Old sheets remain valid and
warn while defaulting to mm39.

## Staged Deployment And Validation

1. Extract to a new versioned directory, not over the active installation.
2. Edit only the staged registry/config paths.
3. Run the syntax, compilation, unit, and preflight commands from
   `TEST_REPORT.md`.
4. Run a small existing mouse sample into a new output directory and compare its
   filtering parameters, pairs counts, resolutions, and matrices with production.
5. Run a small synthetic or pilot hg38 sample into a separate output directory.
6. Switch a launcher/symlink only after review; retain the previous version.

Reference changes invalidate completion metadata. Prefer a new `OUTDIR`. For a
legacy mm39 resume, `ALLOW_LEGACY_MM39_RESUME=true` permits adoption of complete
outputs lacking sentinels and immediately records mm39 metadata. Set it to
`false` for a forced clean rebuild.

The optional `CHECK_TRIMMED_FASTQ_SYNC=true` safeguard fully streams both
alignment FASTQs. Enable it when validating a new installation or investigating
read-name synchronization; it adds I/O and is not evidence of a universal
`fastp` defect.

## Rollback

Stop new submissions, restore the previous launcher/symlink or versioned path,
and point users back to the previous config. Do not reuse candidate outputs in
the old workflow if they contain hg38 or mixed-reference data. Preserve failed
candidate outputs/logs for diagnosis; rollback does not require deleting them.
