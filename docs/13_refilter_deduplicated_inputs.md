# Refilter existing deduplicated mapped pairs

Available in candidate **3.0.0-rc2** as a mode of the existing microc2tracks
workflow. No separate pipeline, remapping or new deduplication pass is required.

## Run mode

```bash
bash scripts/preflight_check.sh --mode refilter-deduplicated \
  -c config/config.conf -s config/samplesheet_refilter.csv

bash scripts/microc2tracks.sh --mode refilter-deduplicated \
  -c config/config.conf -s config/samplesheet_refilter.csv
```

Alternatively set `WORKFLOW_MODE="refilter-deduplicated"` in the config. The CLI
flag overrides that setting. The default `fastq` mode is unchanged.

This mode skips fastp, BWA, alignment parsing and coordinate/PCR deduplication,
even if their ordinary config switches are enabled. Neither FASTQs nor a BWA
index/tool are required. Pairtools is used for sorting/statistics, **not dedup**.
The same pinned environment, shared filtering and matrix builders are used.

Since RC4 the default `HIC_OUTPUT_SCOPE="merged"` skips individual-row `.hic`,
including imported rows that already represent a merged library. To regenerate
`.hic` for such a final input, explicitly select `HIC_OUTPUT_SCOPE="all"`.
Requested new technical/all-sample pools get `.hic` under the default scope.
Other matrix formats and the no-alignment/no-dedup guarantees are unchanged; see
[the output-scope guide](15_merged_hic_outputs.md).

## Input CSV

Required columns:

```csv
sample,assay,reference_genome,deduplicated_pairs
existing_library,hic,my_reference,/data/old_run/library.dedup.pairs.gz
```

Each row is one already-deduplicated input unit, which may represent a previously
jointly deduplicated set of technical inputs. Declare the actual registered
reference and assay; no enzyme or chemistry is inferred from them.

Optional columns:

- `hic_protocol` and `restriction_enzymes`: same profile/enzyme resolution as the
  FASTQ workflow. Otherwise inherit the config/profile values.
- `input_members`: semicolon-separated original member IDs, e.g. `part1;part2`.
  Defaults to known upstream member metadata, or otherwise the row's sample ID.
  Explicit values must agree with known upstream membership. Declare these when inputs are previously
  pooled/joint products so overlapping pools can be rejected.
- `condition`, `biological_replicate`, `technical_replicate`, `merge_group`:
  optional output-pooling metadata. Without a merge_group or a complete set of
  replicate labels, no automatic technical pool is inferred for that row.

Do not supply nonempty `dedup_group` or FASTQ columns in this mode. It cannot
repair incorrect prior deduplication or reconstruct records already discarded.
Existing files need no new-version stage sidecars: selecting this mode is the
user's explicit assertion that inputs are already appropriately deduplicated.
That assertion is recorded as **unverified input history**, not a new dedup pass.

See [the example CSV](../config/samplesheet_refilter_example.csv) and
[example command](../examples/run_refilter_deduplicated_example.sh). All names,
paths, reference IDs and enzyme selections there are illustrative.

## Supported data and validation

Accepts standard text `.pairs` or gzip/BGZF `.pairs.gz` with the first seven
columns `readID chrom1 pos1 chrom2 pos2 strand1 strand2`, plus `pair_type`,
`mapq1`, and `mapq2` in the `#columns:` schema. Mapping qualities and pair types
are never invented. The usual old microc2tracks `.dedup.pairs.gz` files satisfy
this schema.

The `#genome_assembly:` declaration must match the selected registry reference.
Any declared chromosome lengths must match it. All mapped coordinates, strands,
row widths and numeric mapping fields are validated while streaming. Missing,
incompatible or malformed data fails rather than silently changing the reference.
An input sidecar explicitly identifying original pre-dedup pairs is rejected.
If a known FASTA hash is available, it must agree; otherwise assembly/length
checks do not prove which exact FASTA originally produced the alignments.

Direct BAM/SAM/CRAM or matrix import is **not** implemented in this mode. Such
alignments would first need a separately reviewed conversion retaining the
required pair/mapping information. Do not rename them to `.pairs`.

## Processing and configuration

The mode streams existing pairs through:

```text
read-only deduplicated input
  -> validate / canonicalize ends / remove old fragment annotations
  -> shared MAPQ / chromosome / restriction / cis-distance selection
  -> sort / compress / index new accepted pairs
  -> shared matrix builders, QC and optional downstream analyses
```

End swapping also swaps corresponding mapping fields and reverses pair_type.
Old rfrag fields are discarded and regenerated from the currently selected digest.
Sorting never removes matching coordinates: input multiplicity is preserved
except for the explicitly selected filters. No extra full-size imported copy is
saved; streams and sorting scratch feed the new accepted output directly.

For restriction Hi-C, declare the actual enzyme/profile, for example:

```bash
HIC_PROTOCOL="restriction_hic"
HIC_RESTRICTION_ENZYMES="DpnII"  # example only: use the library's actual digest
HIC_MIN_CIS_DIST="0"
```

Zero means no additional distance exclusion; fragment filtering remains enabled.
Choose a positive protocol-appropriate cutoff if wanted. Equality is retained,
and trans contacts bypass the cis cutoff. No universal 1,000-bp Hi-C cutoff is
imposed. Micro-C rows use MICROC_MIN_CIS_DIST and bypass restriction classification.

`BUILD_SAMPLE_MATRICES=false` can be used for filtering/QC only. Disable requested
pools too if no pooled matrices are wanted. Output formats, resolutions, optional
disk roots and retention otherwise follow the normal workflow.

Optional technical/all-sample pools concatenate disjoint accepted units without
another coordinate dedup pass. Duplicate physical input files/aliases are
rejected, as are overlapping declared input_members within a requested pool.
Unreported molecular overlap cannot be inferred from old files; input membership
and prior duplicate scope remain the user's responsibility.

## Preservation, provenance and restart

Use a **new OUTDIR and new category output roots**, outside the input files.
An input located under an output/cache root is rejected. Existing source files
and their `.px2` indexes are never rewritten, indexed in place or cleaned up.

New outputs retain the standard `sample/03_pairs` and `sample/04_matrices`
structure, including `.filter_qc.json`, `.import_qc.json`, accepted pairs/indexes,
matrix products and provenance. Reports identify the input path, mode and
user-asserted deduplication history. No current DEDUP_* or parse settings are
misrepresented as having been applied to the imported data.

Source identities are checked before and after processing. Restart signatures
include the source, reference and actual selection policy; changing unused BWA,
parse or DEDUP_* settings does not trigger refiltering. Changing an enzyme or
distance cutoff rebuilds accepted pairs and dependent matrices, not alignment.
Previously removed contacts cannot be recovered by lowering a cutoff. Any known
upstream filter policy is retained in the import provenance.

The candidate's local scientific/unit and mocked orchestration tests are listed
in [TEST_REPORT.md](../TEST_REPORT.md). The real-tool Linux gate is still required;
the generic Hi-C smoke test now also refilters its joint pre-selection dedup pairs.
