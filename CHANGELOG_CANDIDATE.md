# Changelog

## Unreleased - 3.0.0-rc4 (2026-10-28)

- Added `HIC_OUTPUT_SCOPE="merged"` as the default for both assays: no individual
  technical-row Juicer exports, `.hic` generation/normalization or public tracks;
  requested technical/all-sample pools still get `.hic` after their merges.
- Scope `all` restores individual exports; `RUN_HIC=false` remains the master
  off switch. Unpooled singleton/imported rows are intentionally skipped under
  `merged`, with no automatic extra pooling. Explicit standalone merged matrix
  requests remain supported; internal dedup groups do not acquire extra matrices.
- Left contact filtering/deduplication and cooler processing unchanged. Java/JAR
  are required only for actual targets, which preflight now lists explicitly.
- Reports distinguish intentional skips, requested outputs and preserved older
  files. Scope changes reuse completed pair/cooler stages within this release
  and do not remove previously generated `.hic` or retained Juicer/raw products.
- Added six focused mocked regressions, including merged export/publication,
  no-target dependency checks, scope-change reuse and failed-normalization resume.
  Extended the unexecuted tiny real-tool suite's merged-only assertions.
- Updated the shared installation guide to stage RC4 outside home directories;
  existing server code/jobs/environments/references are not modified.

## Unreleased - 3.0.0-rc3 (2026-10-27)

- Added a shared Linux launcher for a read-only versioned installation outside
  user homes, using an existing environment without conda activation. It resolves
  symlinks to the selected release, preserves arguments/current directory and
  avoids user-site/PYTHONPATH contamination and source-tree bytecode writes.
- Added `check-access` for current-account tool, code, reference/index and JAR
  permissions, with explicit warnings for unchecked resources and home paths.
  It makes no changes and does not claim every account or scientific validation.
- Replaced outdated installation recipes with administrator-only code staging,
  cross-account checks and user/project-owned output/scratch/cache guidance.
  Existing environments, running jobs and genome indexes are not changed.
- Added five focused launcher/access regressions. No scientific filtering,
  deduplication or matrix algorithms were changed by this shared-install update.

## Unreleased - 3.0.0-rc2 (2026-10-27)

- Added `--mode refilter-deduplicated` to the existing runner/preflight, also
  selectable with `WORKFLOW_MODE`. Default FASTQ mode is unchanged.
- Added an explicit CSV schema for existing deduplicated mapped pairs, with no
  new-version sidecar requirement and no FASTQ/BWA-index requirement.
- Skips trimming, alignment, parsing and deduplication. Streams validated input
  through canonicalization, restriction/distance selection, sorting/indexing and
  the same matrix builders; no extra full-size imported copy is retained.
- Keeps source pairs/indexes read-only, requires new output roots, records prior
  deduplication as a user assertion, and detects declared/physical input overlaps.
- Adds per-input import QC, source/provenance identity checks, policy-specific
  restart, examples and ten focused tests. The generic real-tool Hi-C smoke test
  now also covers refiltering a joint dedup stream; real-tool execution remains
  an outstanding release gate, not a claimed pass.

## Unreleased - 3.0.0-rc1 (2026-10-27)

One complete integration candidate; not yet promoted to production.

- Kept Hi-C and Micro-C behind the same existing shell commands, using shared
  planning, pairs, matrix and artifact-state helpers.
- Added explicit restriction profiles, DpnII/MboI/NlaIII/HindIII digests and
  combinations, optional original-kit Dovetail preset, fragment artifact
  classification and consistent assay-specific distance filtering.
- Added optional `dedup_group` metadata and a single joint pass over original
  sorted streams, while preserving independent sample products. No UMIs.
- Added technical/all-sample pools of disjoint units and overlap/partial-group
  rejection, plus explicit standalone pool/joint-dedup/legacy-pool modes.
- Added optional output roots, isolated SSD scratch, HDD destination staging,
  verified cross-filesystem publication, reserves/mount checks and owned cleanup.
- Added parameter/input/tool fingerprints, stage-specific restarts, safe managed
  relocation without recomputation, explicit failure states and live tool logs.
- Added plan-aware summaries/browser packaging and custom-reference support
  without mandatory human/mouse chromosome names.
- Added generic examples, target environment pins, small regressions and tiny
  opt-in real-tool tests. Linux real-tool tests and an environment solve remain
  release gates; local mocked contracts do not establish scientific-tool success.
- Preserved the earlier LF and missing-HTML-value fixes and fastp guard. Historical
  `SERVER_BASELINE_SHA256SUMS.txt` describes the pre-integration server source,
  not this candidate's changed files.

See [the integration guide](docs/12_integrated_hic_and_storage.md) and
[verification results](TEST_REPORT.md).

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
