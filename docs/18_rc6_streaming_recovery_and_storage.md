# RC6: bounded streaming, recovery and large-file efficiency

RC6 implements the retained items from the 2026-11-04 revised RC5 plan as one
release candidate. It does not change duplicate matching, scientific contact
filters, references, resolution settings, normalization methods, or library
membership. It does not execute deployment or production analyses automatically.

## Processing and checkpoints

```text
per-run fastp -> alignment/parse/sort -> retained original sorted pairs
    -> checked sorted merge | one library dedup pass -> retained pre-filter dedup
    -> restriction/distance filtering -> accepted pairs + index/QC
    -> biological mcool
    -> explicit final accepted pool + compatible raw-count pooling
    -> newly normalized final mcool + final hic
```

`dedup_group` explicitly declares one physical-library deduplication scope.
Biological labels never implicitly create a PCR-deduplication group. Independent
libraries retain identical-coordinate contacts when pooled.

`MERGE_DEDUP_MODE=stream` removes the complete joint sorted intermediate unless
`KEEP_RAW_MERGED_PAIRS=true` or a valid completed merge already exists. Ordinary
bounded pipes apply backpressure; the whole library is not held in RAM.
`materialized` keeps a separate completed-merge checkpoint. A failed fused stage
restarts merge/dedup from retained sorted inputs, not halfway through a stream.

Every input reader is waited, verifies gzip completion, and writes an input
identity/count/logical-digest receipt. The deduplication total must agree with
the merged count. The same GNU comparator/locale/full-line tie fallback is used.
`MERGE_READERS=2` limits active decoding through short-lived cross-process locks;
it does not limit which input files participate or hold locks while pipes block.
GNU merge fan-in is separately controlled by `MERGE_FAN_IN=8`. Above that fan-in,
uncompressed internal merge runs can be large even though output is streamed.
`MAX_MERGE_INPUTS=256` bounds total reader processes. Larger libraries require
bounded sorted pre-merges followed by one joint dedup pass; do not deduplicate
those batches separately.

Keep `DEDUP_BACKEND`, `DEDUP_METHOD`, `DEDUP_MAX_MISMATCH`, `DEDUP_CHUNKSIZE`, and
`DEDUP_CARRYOVER` unchanged. Existing tests demonstrate chunk-sensitive duplicate
representatives. Changing those values or sharding dedup is not a neutral speedup.

## Recovery with FASTQs removed

Use the same reviewed config/samplesheet and managed OUTDIR:

```bash
microc2tracks-rc6 preflight --resume-checkpoints -c config.conf -s samplesheet.csv
microc2tracks-rc6 run --resume-checkpoints -c config.conf -s samplesheet.csv
```

Only run the second command when intentionally ready to resume. The option also
exists as `RESUME_FROM_CHECKPOINTS=true`. It does not authorize legacy-file adoption.
Startup verifies prior source paths, membership, reference identity and upstream
processing/tool contracts before using recorded sorted-stage lineage. Normal
stage checks must still find a complete reusable sorted, deduplicated or accepted
checkpoint. If none remains, missing FASTQs cannot be recreated automatically.

RC5 reuse requires the exact reviewed source fingerprint. RC6 records upstream
keys/input identities and saves the preceding run plan under
`run_metadata/history/<run_id>/previous_run_plan.json`. Metadata/provenance of
reused scientific artifacts are not rewritten merely to mark them as reused.
Changes to duplicate criteria invalidate deduplication, not alignment; changes
to restriction/distance filters invalidate accepted contacts, not alignment.

For older unmanifested outputs, use the existing declared `joint-dedup` merge or
`refilter-deduplicated` input modes with appropriate provenance. Already accepted
or deduplicated inputs cannot restore contacts removed upstream. Do not bypass
manifest validation or relabel accepted files as pre-deduplication sorted input.

## Products and QC

- Filtering uses NumPy batches bounded by both `FILTER_BATCH_RECORDS=50000` and
  `FILTER_BATCH_BYTES=16777216`. A single larger record fails explicitly rather
  than exhausting an unbounded queue. Output order and all QC counters are kept.
  `CONTACT_FILTER_ENGINE=scalar` selects the bounded reference loop instead.
  The small synthetic check was approximately 1.44 s scalar versus 1.34 s batched,
  with roughly 56 versus 97 MiB peak process RSS and identical output/QC. An
  earlier smaller-digest fixture was slightly slower when batched. These are
  kernel observations, not evidence for a large or universal workflow speedup.
  See [the recorded synthetic result](rc6_filter_benchmark.json).
- Juicer projection runs in a monitored subprocess, with identical zero-based
  short-format coordinates and separate checked bgzip compression. Configure
  `JUICER_COMPRESSION_LEVEL=3`; compression differences are not contact changes.
- `POOL_STATS_MODE=merge` aggregates compatible accepted counters/histograms and
  recalculates derived summaries. Standalone external pools without managed
  constituent statistics fall back to scanning. `scan` restores the old path.
- `POOL_COOLER_MODE=counts` uses complete constituent raw count matrices (or raw
  pixels in a normalized mcool) for the final pool when their provenance and
  exact membership match. Weight columns are never pooled. Missing constituent
  products fall back to pooled-pairs binning; stale products fail validation.
  `pairs` always uses the original binning path.
- Integer totals/axes are checked and pooled counts use int64. Cooler merge
  batches can exceed their target for a single dense bin row; its buffer is not
  a hard memory limit. Final counts are normalized anew using unchanged methods.
- Essential accepted-output validation finishes before matrices/export. Full
  descriptive statistics may overlap them, but remain a completion prerequisite.
  Final pooling still waits for selected biological outputs to finish.

## Resource and storage policy

Empty alternative output settings still use the standard structure under
`OUTDIR`. Use HDD paths for durable sorted/dedup/accepted pairs and final products.
Use SSD only for deliberately budgeted active scratch and a small trim backlog.
`LARGE_OUTPUT_STAGING=destination` avoids copying complete new large outputs from
SSD to HDD; `scratch` can temporarily hold both copies during publication.

The optional `config/server_140cpu_500gb.conf` supplies two 50-thread aligners,
140 total accounted CPU slots, 400 GiB requested RAM admission, 32 GiB available
host-memory headroom, two bounded-decoder joint stages with 96 GiB reservations,
and a 128-GiB estimated prepared-input backlog. Lower the total budget when
sharing resources with other users. Budgets apply to this run, not all server jobs.

RAM reservations are admission estimates. `MAX_PROCESS_TREE_RSS_GB` optionally
stops monitored child work when summed Linux RSS exceeds its polling limit.
Shared pages are counted repeatedly; this is conservative and not PSS/cgroups.
`run_metadata/resource_usage.json` records recent measurements. Scientific
Python child libraries retain one native BLAS/OpenMP thread per admitted process.

Disk admission groups roots by underlying filesystem device. Concurrent tasks
reserve estimated growth plus the configured free-space floor; cross-filesystem
staging/publication requests both locations. Estimates initially use compressed
input sizes with `SCRATCH_ESTIMATE_FACTOR=6` and `OUTPUT_ESTIMATE_FACTOR=2`.
Override/calibrate factors or minimum `STAGE_*_RESERVATION_GB` from measurements;
these are not claims of a universal compression ratio. Conservative outstanding
reservations can reduce concurrency. An insufficient-space queue is logged.

`SCRATCH_BUDGET_GB` checks this run's scratch tree, not unrelated contents of a
shared `/tmp`. Old partials elsewhere still consume filesystem free space.
For a system SSD with approximately 1 TB currently free, explicitly consider:

```bash
SCRATCH_BUDGET_GB="250"
MIN_FREE_SCRATCH_GB="250"
```

These illustrative settings are not enough for every sort. A stage whose
predicted scratch exceeds the budget fails before starting; choose explicit HDD
scratch, a larger budget, or a justified estimate. Do not silently write onto `/`
when an expected HDD mount is missing: declare `REQUIRED_STORAGE_MOUNTS`.
Concurrent run/global-user enforcement requires administrator quotas/cgroups.

Publication journals identify generation-owned final/transfer files before any
rename. A completion record is written last. Interrupted bundles are never
silently adopted as complete, and retries never append to earlier partials.
Checksums collected during sorted writes can be checked by later merge readers;
legacy path/size/mtime records cannot retrospectively prove absence of same-size
corruption. No automatic removal of previous-run outputs is added.

## Validation and deployment

See `TEST_REPORT.md` for actual results and unexecuted gates. There are three small
test levels, not three implementation releases:

```bash
MICROC2TRACKS_REAL_TOOLS=0 python -B -m unittest discover -s tests -p 'test_*.py' -q
MICROC2TRACKS_PAIR_TOOLS_TESTS=1 python -B -m unittest discover -s tests -p 'test_rc6_real_science.py' -v
MICROC2TRACKS_REAL_TOOLS=1 JUICER_TOOLS_JAR=/path/to/existing.jar \
  python -B -m unittest discover -s tests -p 'test_real_tools.py' -v
```

Stage RC6 alongside RC5 at `/opt/microc2tracks/versions/3.0.0-rc6`, root-owned and
readable/executable by other accounts. Add a version-specific launcher first;
reuse the existing conda environment, references, indexes and JAR read-only.
Verify `--version`, `check-access`, imports, the small server real-tool gate and
preflight before promoting the default launcher. Do not overwrite an active
release, change software permissions, rebuild production indexes or launch
production data as part of installation. Reverting the launcher to RC5 is a
deployment rollback, not authorization to overwrite RC6 scientific artifacts.

Direct native-pairs/FIFO input into Juicer, changed deduplication algorithms or
chunk sizes, and earlier final-pool admission remain deferred by design.
