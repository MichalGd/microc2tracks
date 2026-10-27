# Changelog

## 2026-10-27 - Installed-server source synchronization

Source update: [`8abf4cd`](https://github.com/MichalGd/microc2tracks/commit/8abf4cd06e27d6d7b609ff318d30e79bff2f041d).

### Fixed

- Normalized sample-sheet CSV/TSV now uses LF line endings, preventing an empty
  trailing `merge_group` from becoming a carriage-return directory name in Bash.
- HTML summary cells now pass through the existing `as_text` helper before HTML
  escaping, preventing crashes on missing/`None` values while preserving escaping.

### Restored

- The installed server's independent fastp concurrency guard in the runner and
  its validation/reporting in preflight. `MAX_PARALLEL_FASTP` defaults to `1`
  when omitted; overall sample parallelism is unchanged.

### Verification and scope

- Added `SERVER_BASELINE_SHA256SUMS.txt`: all 18 enumerated workflow scripts,
  `environment.yml`, and the config template match the 20 supplied server hashes.
- Added three regression tests for normalized line endings and HTML cell values.
  Recorded checks passed: 18 Python tests, Bash syntax, reference-sentinel
  regression, and Git whitespace validation.
- Source identity is based on supplied fingerprints, not a live server inspection
  or a production-data rerun. No server files or analysis outputs were changed.
- No scientific filtering defaults, merge behavior, or matrix-generation changes.
  Restriction-aware Hi-C processing, cross-input deduplication, optional output
  roots, and the later all-sample pooling extension are not part of this update.

See [source synchronization details and test limitations](docs/11_server_baseline_sync.md).

## 2026-08-01 - Published baseline

### Human And Per-Sample Reference Support

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

### Backward Compatibility

- Sample sheets without `reference_genome` warn and default to mm39.
- Existing global reference config values override the default mm39 registry
  row, preserving installed mouse paths.
- Existing output names/layout remain unchanged.
- Missing-sentinel adoption remains available only for legacy mm39 outputs and
  can be disabled with `ALLOW_LEGACY_MM39_RESUME=false`.
- Scientific defaults for MAPQ, Micro-C cis distance, BWA arguments, matrix
  resolutions, Juicer normalization, retention, and downstream analyses were
  preserved.

### Production Behavior Incorporated

- Kept the production gzip-suffix-safe atomic temporary names in the runner and
  merge workflow.
- Added production-only saddle comparison and UCSC/HiGlass packaging utilities.
- Set the public template to the requested production profile:
  `THREADS_FASTP=8`, `MAX_PARALLEL_SAMPLES=2`.
- Added optional full alignment-FASTQ name/count synchronization checking after
  trimming (`CHECK_TRIMMED_FASTQ_SYNC=false` by default).

### Minor Corrections Separate From Human Support

- Made merge manifests self-describing with headers and updated report parsing,
  avoiding a header being interpreted as a data group.
- Made browser reference/preset contradictions explicit errors.
- Preserved compression suffixes on atomic temporary outputs, incorporating the
  low-risk production fix.
- Expanded diagnostics for chromosome mismatch, mitochondrial naming, missing
  index files, stale sentinels, and cross-reference merges.
- Restricted sample and explicit merge-group identifiers before using them in
  output paths.
