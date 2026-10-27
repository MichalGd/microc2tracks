# Integrated Hi-C, joint deduplication and optional storage roots (3.0.0-rc2)

This is one update to **microc2tracks**, not a second Hi-C pipeline. The existing
commands remain:

```bash
bash scripts/preflight_check.sh -c config/config.conf -s config/samplesheet.csv
bash scripts/microc2tracks.sh -c config/config.conf -s config/samplesheet.csv
```

The shell entry points share Python helpers for planning, pair processing,
matrix generation, artifact state and publication. No workflow engine is added.
To start from existing mapped, already-deduplicated pairs, use the explicit
[refilter-deduplicated run mode](13_refilter_deduplicated_inputs.md). It skips
alignment/deduplication and writes a separate output tree; it does not adopt
legacy files as new-stage completion proof.
Example paths, references, sample names and resource values are illustrative.
No server installation or production output is automatically changed.

## Compatibility and release gate

Existing Micro-C configs/sheets retain the assay's pair-type, MAPQ and inclusive
minimum-distance rules and normal output names. Restriction processing never
runs for Micro-C. New columns are optional. The processing target is pairtools
1.1.3 with explicit duplicate settings; preflight refuses a different version.
This does **not** claim equivalence to an unknown older server dependency stack.
Create a separate environment from `environment.yml`, not an in-place update to
the environment of a running job. Its version pins are target specifications,
not a solved Linux build lock. Actual tool/library versions are saved per run.

The candidate's small Python/contract tests and shell checks can run without
genomics tools. Tiny real-tool tests are provided separately. They must pass in
the target Linux environment before production promotion; mocked orchestration
tests do not establish real-tool compatibility or matrix correctness.

Reference-only `.done` markers are not adopted as enhanced-stage proof. Use a
**new OUTDIR** for old runs. Nonempty unowned output files are never overwritten.
Old configs' `BOOTSTRAP_SENTINELS`/legacy-adoption switches do not bypass this
safety rule. Existing originals and prior run metadata remain available.

## Protocol and sample-sheet metadata

| Optional column | Meaning when nonempty | Blank/absent |
|---|---|---|
| `hic_protocol` | `legacy`, `restriction_hic`, `dovetail_hic` | Inherit config; old configs without HIC_PROTOCOL fall back to legacy |
| `restriction_enzymes` | Explicit `+`-separated supported enzyme names | Inherit config or selected preset |
| `dedup_group` | User-selected rows to deduplicate together once | Independent singleton row |

Existing `merge_group` defines an output pool, **not** duplicate scope. Existing
technical/biological labels do not silently enable joint deduplication. No
`library_id` questionnaire or UMI fields are required. Only assign the same
dedup_group where joint coordinate/PCR deduplication is intended; independent
biological libraries must not be grouped indiscriminately.

The generic template uses:

```bash
HIC_PROTOCOL="restriction_hic"
HIC_RESTRICTION_ENZYMES=""  # must be resolved explicitly for restriction Hi-C
HIC_MIN_CIS_DIST="0"        # disables distance filtering, not fragment filtering
HIC_REMOVE_SAME_FRAGMENT="true"
HIC_REMOVE_ADJACENT_INWARD="false"
```

Micro-C rows leave the Hi-C columns blank. `assay=hic` selects assay processing;
`RUN_HIC` controls **.hic file generation for either assay**. These are different
settings. CSV is the supported sample-sheet format; export Excel workbooks first.
CRLF CSV and empty optional fields are parsed structurally, not by Bash `read`.
Configs are trusted Bash code; sample metadata is never evaluated as shell code.

### Enzymes and filtering

`config/hic_protocols.json` defines DpnII (`^GATC`), MboI (`^GATC`), NlaIII
(`CATG^`) and HindIII (`A^AGCTT`). Use one or an explicit combination such as
`DpnII+NlaIII`. Combined digests union cut positions. Recognition motifs are not
inferred from an assembly. The initial extensible registry validates palindromic
literal A/C/G/T motifs and their plus-strand cut offsets.

The optional `dovetail_hic` profile denotes the original DpnII-based kit, not
Omni-C or LinkPrep. It rejects contradictory enzyme settings. It adds no
unverified kit-specific aligner changes, bridge trimming or UMI processing.

Digests are built from the exact FASTA and validated chromosome sizes. Cache
identity includes reference content, enzyme definitions and generator version;
concurrent generation is locked and partial entries are not treated as complete.
Intervals are zero-based half-open. Pair positions stay one-based and become
`pos-1` at fragment assignment. Restriction pairs receive `rfrag1/rfrag2` columns.

Selection after deduplication is shared by samples and joint groups:

1. Eligible pair types (`UU`, `UR`, `RU`) and both MAPQ values at/above threshold.
2. Declared chromosome policy; eligible mapped positions must fit the reference.
3. Restriction profiles: classify fragments/orientations, then apply rejection rules.
4. Assay-specific minimum cis distance: positive values retain equality; trans
   pairs bypass this cutoff. Zero disables it.

Same-fragment inward pairs are dangling-end candidates; outward pairs are
self-circle candidates; equal-position/same-strand configurations are other
same-fragment contacts. These are rejected by default in restriction profiles.
Adjacent-fragment inward pairs are counted but retained unless explicitly
disabled. Different chromosomes never share fragment identity, even when their
local fragment ordinals match. Inconsistent digest/reference data is an error.

These are operational artifact candidates, not proof of molecular origin.
Same-fragment or distance rules can discard genuine short-range contacts.
Adjacent-fragment filtering can discard genuine contacts as well; hence its
conservative default. No arbitrary distance-to-recognition-motif gate is added.
Legacy Hi-C explicitly warns that it has no fragment filtering.

## Joint deduplication and pooling

Each row retains sorted original **pre-dedup** pairs and independent products.
Rows in a declared dedup_group are merged from those original sorted streams,
then deduplicated once and filtered. Internal joint products live below
`dedup_groups/<id>`, separate from individual samples and `merged/<pool>` outputs.

```bash
DEDUP_BACKEND="scipy"       # or explicit legacy-compatible cython
DEDUP_METHOD="max"          # or sum
DEDUP_MAX_MISMATCH="3"
DEDUP_CHUNKSIZE="10000"
DEDUP_CARRYOVER="100"
KEEP_SORTED_PAIRS="true"
BUILD_SAMPLE_MATRICES="true"
RUN_MERGE_TECHNICAL_REPLICATES="true"
RUN_MERGE_ALL_SAMPLES="false"
ALL_SAMPLES_MERGE_NAME="all_samples_merged"
```

Chromosomes, positions and strands determine coordinate matches, not read IDs.
Backend/tolerance/chunk/carryover values and dependency versions are recorded.
The pinned scipy implementation is transitive within a chunk, but carries
representatives rather than discarded bridge records across chunks. Therefore
fuzzy results can depend on chunk boundaries; this release matches the explicit
pinned single-pass policy, not a claim of global transitive closure.
[Pinned implementation](https://github.com/open2c/pairtools/blob/v1.1.3/pairtools/lib/dedup.py#L130-L201).

Technical and optional all-sample pools consume disjoint dedup units once.
They cannot contain both a joint group and its individual rows, and cannot
select only part of a group. Incompatible references, assays, profiles, digests
or filter policies fail before processing. Independently deduplicated units
are pooled **without another coordinate dedup pass**. Individual counts need
not sum to joint counts; their difference is not an exact PCR-molecule count.

Standalone merging keeps the old command flags, with explicit input modes:

```bash
# Verified accepted disjoint units; default mode is pool.
bash scripts/merge_replicates.sh -c config/config.conf -n my_pool -r my_reference \
  --mode pool -p /data/unit1.pairs.gz,/data/unit2.pairs.gz

# Original per-row sorted pre-dedup streams from this version, with sidecars.
bash scripts/merge_replicates.sh -c config/config.conf -n my_joint -r my_reference \
  --mode joint-dedup -p /data/row1.sorted.pairs.gz,/data/row2.sorted.pairs.gz

# Explicit compatibility mode for old accepted pairs lacking stage provenance.
bash scripts/merge_replicates.sh -c config/config.conf -n old_pool -r my_reference \
  --mode legacy-pool -p /data/old1.pairs.gz,/data/old2.pairs.gz
```

Strict modes validate stage sidecars and payload identities. Legacy pool mode
checks assembly and selected chromosomes only; upstream filtering and duplicate
scope remain unverified and are reported as such. It cannot certify enhanced
pool inputs. Manual merges keep their plan/report under
`run_metadata/manual_merges/<name>` and do not replace the full-run plan/report.
Standalone joint-dedup resolves the desired filtering profile from its config;
standalone strict pooling inherits its inputs' verified filtering policy and
uses the config for matrix/storage settings, not for recovering discarded pairs.

Already independently deduplicated files can be reused for simple independent
pooling, but cannot reconstruct a correct original-stream joint fuzzy pass.
For joint reprocessing, use retained original sorted pairs; otherwise regenerate
them from earlier alignments/FASTQs. No automatic legacy-file adoption is added.

## Optional disks, safe publication and restart

`PAIRS_OUTDIR`, `MATRIX_OUTDIR`, `TRIMMED_OUTDIR` and `REFERENCE_CACHE_DIR` are
optional. Empty category roots retain the standard relative layout under OUTDIR;
the cache defaults to `OUTDIR/reference_cache`. Empty TMPDIR uses `OUTDIR/tmp`.

For example, accepted pairs remain `sample/03_pairs/sample.valid.mapq30.pairs.gz`
relative to either OUTDIR or PAIRS_OUTDIR. Pools keep `merged/name/03_pairs` or
`04_matrices`; joint units keep `dedup_groups/id`. Main metadata, logs and reports
stay under OUTDIR. All consumers use absolute paths from `run_metadata/run_plan.json`
and `artifacts.tsv`. Browser packaging accepts `--run-plan` plus optional `--entity`.

A practical layout is OUTDIR on a large HDD and TMPDIR on SSD scratch. Final
matrices/reports can stay on HDD. Retained sorted/accepted pairs can also stay
there, with slower merging, filtering and matrix input scans. Sort spills and
temporary merge files benefit most from SSD. Active matrix staging can benefit
from SSD but requires space for the large output and its copy. Repeated digest
loading and future matrix browsing/downstream access can be slower on HDD.

`LARGE_OUTPUT_STAGING=destination` writes unique temporary siblings on the output
disk, minimizing SSD use. `scratch` builds in SSD scratch, then verifies size and
SHA-256 after copying to a destination sibling before the final atomic rename.
Success is recorded only after all bundle outputs/indexes publish. Interrupted
transfers/checksum failures never receive success; valid source data is retained.

Set `REQUIRED_STORAGE_MOUNTS` to the expected mount points to detect absent HDDs
instead of accidentally writing to the system filesystem beneath a mount.
Free-space reserves are checked before stages and throughout external commands.
Filesystem identity changes fail the run. Reserves are safety margins, not exact
expansion predictions or hard quotas. Use a dedicated scratch partition/quota
for a strict system-disk cap, and account for per-worker memory plus concurrency.
Concurrency limits are per run; simultaneous independent runs still share RAM/I/O.

Scratch has a unique run/entity/stage namespace. OS-released target locks prevent
concurrent writes. Cleanup deletes only recorded, unchanged, owned files within
allowed roots; it never recursively deletes a storage root or follows a convenience
symlink. Sorted pairs are removed, if requested, only after all consumers finish.
PUBLIC_HIC_DIR remains empty by default to avoid an unexpected large extra copy.

Stages record scientific signatures, output identities and state. Input FASTQs
use path/size/mtime identity; the reference uses content hashes. Changing distance
or enzymes reuses suitable earlier stages; changing group membership rebuilds the
joint unit and dependents; changing resolutions rebuilds matrices only. A managed
storage-root change republishes verified outputs without recomputation and keeps
old copies. Migrating untracked legacy outputs is not automatic.

Run states distinguish success, partial_success, failed and interrupted. Optional
downstream failures are recorded and give exit code 2/partial_success rather than
being marked complete. Long tool stdout/stderr remains live in per-stage logs.

## Minimum verification

```bash
python -m unittest discover -s tests -p 'test_*.py' -v
bash tests/check_bash_syntax.sh
bash tests/test_reference_sentinel.sh

# In the isolated Linux target environment, with the actual configured JAR:
MICROC2TRACKS_REAL_TOOLS=1 JUICER_TOOLS_JAR=/path/to/juicer_tools.2.18.00.jar \
  python -m unittest discover -s tests -p test_real_tools.py -v
```

The opt-in suite creates only tiny synthetic references/FASTQs in temporary
directories. It runs Micro-C, generic restriction Hi-C and the Dovetail preset,
joint deduplication, split roots, restart and raw `.mcool`/`.hic` bin-count
comparisons. A small pinned chunk/bridge case records fuzzy matching limitations.
No production-data reruns or extensive benchmarks are required.
