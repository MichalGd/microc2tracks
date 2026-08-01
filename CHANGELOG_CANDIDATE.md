# Candidate Changelog

## Human And Per-Sample Reference Support

- Added controlled per-row `reference_genome` selection with `mouse`/`mm39`
  and `human`/`hg38` aliases normalized to canonical assembly IDs.
- Added one validated TSV registry for species, assembly, FASTA, BWA-MEM2
  prefix, chromosome sizes, canonical regex, browser preset, optional phasing
  track/assembly, and annotation/blacklist metadata.
- Propagated canonical reference identity through BWA-MEM2, Pairtools assembly
  metadata, chromosome filtering, Cooler assembly metadata, `.mcool`, `.hic`,
  downstream phasing, browser packaging, status files, manifests, sentinels,
  and final reports.
- Added used-reference preflight for exact `.fai`/chromosome-size agreement,
  BWA-MEM2 components, expected canonical chromosomes, mitochondrial naming,
  browser preset, and optional phasing track.
- Prohibited cross-reference technical-replicate merges in both sample-sheet
  grouping and pairs-header validation.
- Bound resumability sentinels to `reference_id` and assembly.
- Added assembly-explicit browser packaging; hg38 cannot silently select the
  mouse preset.

## Backward Compatibility

- Sample sheets without `reference_genome` warn and default to mm39.
- Existing global reference config values override the default mm39 registry
  row, preserving installed mouse paths.
- Existing output names/layout remain unchanged.
- Missing-sentinel adoption remains available only for legacy mm39 outputs and
  can be disabled with `ALLOW_LEGACY_MM39_RESUME=false`.
- Scientific defaults for MAPQ, Micro-C cis distance, BWA arguments, matrix
  resolutions, Juicer normalization, retention, and downstream analyses were
  preserved.

## Production Behavior Incorporated

- Kept the production gzip-suffix-safe atomic temporary names in the runner and
  merge workflow.
- Added production-only saddle comparison and UCSC/HiGlass packaging utilities.
- Set the public template to the requested production profile:
  `THREADS_FASTP=8`, `MAX_PARALLEL_SAMPLES=2`.
- Added optional full alignment-FASTQ name/count synchronization checking after
  trimming (`CHECK_TRIMMED_FASTQ_SYNC=false` by default).

## Minor Corrections Separate From Human Support

- Made merge manifests self-describing with headers and updated report parsing,
  avoiding a header being interpreted as a data group.
- Made browser reference/preset contradictions explicit errors.
- Preserved compression suffixes on atomic temporary outputs, incorporating the
  low-risk production fix.
- Expanded diagnostics for chromosome mismatch, mitochondrial naming, missing
  index files, stale sentinels, and cross-reference merges.
- Restricted sample and explicit merge-group identifiers before using them in
  output paths.
