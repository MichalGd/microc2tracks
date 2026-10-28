# Generate .hic files only for merged outputs

Introduced in candidate **3.0.0-rc4**, 2026-10-28. This is a matrix-output choice
for **both Micro-C and Hi-C assays**, not a change to restriction/distance
filtering, deduplication scope or matrix counts.

## Default configuration

```bash
RUN_HIC="true"
HIC_OUTPUT_SCOPE="merged"
```

Individual sample-sheet rows retain their accepted pairs, QC and cooler outputs
according to the existing matrix/retention switches. They do not run Juicer pair
export, `java ... pre` or `java ... addNorm`, create raw/normalized `.hic`, or
publish new `.hic` copies/UCSC links. This removes the expensive per-technical-row
`.hic` work; it does not disable `assay=hic` processing.

Technical-replicate pools are processed after the individual/joint units. Each
pool receives its `.hic` after its merged accepted pairs are ready. If optional
all-sample pooling is enabled, that pooled product also gets `.hic`. Existing
matrix/merge concurrency limits still apply; this does not introduce a separate
barrier waiting for every pool's cooler to finish before any pool's `.hic` starts.

Normal internal dedup-group products do not get matrices just because they combine
inputs. An explicit standalone `merge_replicates.sh --mode joint-dedup` request
does build its merged matrices, including `.hic` when enabled.

No extra merging is inferred to create a `.hic`. If a row has no merge partner,
technical merging is disabled, or no pools are requested, there is no `.hic` for
that row in `merged` scope. Preflight lists exact `hic_targets` and prints a
warning if `RUN_HIC=true` resolves to no targets. Intentional skips are not errors.

## Explicit alternatives

```bash
# Previous behavior: individual sample matrices AND requested merged matrices.
HIC_OUTPUT_SCOPE="all"

# Master switch: no .hic at all, regardless of scope.
RUN_HIC="false"
```

`BUILD_SAMPLE_MATRICES=false` continues to suppress all individual matrix formats,
even with scope `all`. `.cool`/`.mcool` behavior is otherwise unchanged. The default
applies to old configs that omit the new setting as well as the new template.

Refilter mode follows the same rule. To regenerate `.hic` from one input row that
already represents a final merged/deduplicated library, set scope `all` explicitly.
An `input_members` annotation alone does not turn that row into a new pool.

## Restart, reporting and safety

Changing only the scope does not change pair filtering, dedup or cooler stage
signatures within this release. Enabling individual exports later reuses existing
accepted pairs/coolers; it only adds the missing `.hic` work. Skipping individual
exports leaves previously generated `.hic` and retained raw/Juicer products
unchanged, rather than deleting them to match the new scope.

This does not mean arbitrary legacy outputs can be adopted: existing stage
ownership/provenance and new-OUTDIR migration safeguards remain in effect.

Reports include `hic_output_scope`, `hic_requested` and `hic_output_status`,
distinguishing an intentional omission from a requested missing/failed product
and a preserved older file. Java/JAR are required/hashed only if at least one
actual `.hic` target exists. A merged `.hic` failure preserves completed sample
outputs; restart can reuse successful pairs/cooler/raw-HIC stages.

The shared outside-home installation model remains unchanged. RC4 is staged as
a new version, never by overwriting code, settings or jobs currently running.
The tiny real-tool Linux release tests still need to be run before production
promotion; local mocked output tests do not replace them.
