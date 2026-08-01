# Mouse And Human References

## Design

The sample sheet keeps the user-facing `reference_genome` column. Values are
resolved through `config/references.tsv` before any paths are used. The shipped
aliases are controlled and case-insensitive:

| Input | Canonical ID | Species | Assembly | Browser preset |
|---|---|---|---|---|
| `mouse`, `mm39` | `mm39` | *Mus musculus* | mm39 | `mouse` |
| `human`, `hg38` | `hg38` | *Homo sapiens* | hg38 | `human` |

Unknown values are errors. Sample-sheet text is never evaluated as a shell
expression and is never interpolated into a reference path.

The registry contains the canonical ID, aliases, species, assembly, FASTA,
BWA-MEM2 index prefix, chromosome sizes, canonical regex, browser preset,
optional phasing track/assembly, and annotation/blacklist metadata. Annotation
and blacklist metadata are recorded only; the core workflow does not currently
apply gene annotations or blacklist filtering.

## Sample Sheets

Mouse only:

```csv
sample,assay,reference_genome,condition,biological_replicate,technical_replicate,fastq_r1,fastq_r2
mouse_1,microc,mouse,control,1,1,/data/mouse_1_R1.fastq.gz,/data/mouse_1_R2.fastq.gz
```

Human only:

```csv
sample,assay,reference_genome,condition,biological_replicate,technical_replicate,fastq_r1,fastq_r2
human_1,microc,hg38,control,1,1,/data/human_1_R1.fastq.gz,/data/human_1_R2.fastq.gz
```

Mixed-reference runs are supported. Use distinct conditions/biological groups
where appropriate because otherwise-compatible technical replicate rows may not
cross references:

```csv
sample,assay,reference_genome,condition,biological_replicate,technical_replicate,fastq_r1,fastq_r2
mouse_1,microc,mm39,mouse_control,1,1,/data/mouse_R1.fastq.gz,/data/mouse_R2.fastq.gz
human_1,microc,human,human_control,1,1,/data/human_R1.fastq.gz,/data/human_R2.fastq.gz
```

An optional `merge_group` column can explicitly name a technical-replicate
group. It accepts only letters, digits, dot, underscore, and hyphen. A group
cannot contain different references.

## Backward Compatibility

If `reference_genome` is absent, every row uses `DEFAULT_REFERENCE_ID`, which
defaults to `mm39`, and validation emits a warning. A blank field behaves the
same way and warns. Explicit `human`/`hg38` is never changed to mouse.

For the default reference only, old `GENOME_ASSEMBLY`, `REFERENCE_FASTA`,
`BWA_INDEX_PREFIX`, `CHROM_SIZES`, `CANONICAL_CHROMS_REGEX`, and
`PHASING_TRACK` config values override the matching registry entry. This lets an
existing mouse server config continue to use its installed files. New references
must be configured in the registry.

Existing output directory and file names are retained. New manifests/status
files add reference columns. Completion sentinels now contain `reference_id` and
`assembly`. A sentinel from another reference is rejected and the step rebuilds.
Legacy outputs without sentinels may be adopted only for `mm39` when both
`BOOTSTRAP_SENTINELS=true` and `ALLOW_LEGACY_MM39_RESUME=true`; this is the
documented migration bridge. Human outputs are never adopted without matching
reference metadata.

When a sample's reference changes, use a new `OUTDIR` if possible. Otherwise
remove that sample's output directory (and its corresponding merged group) after
archiving it. Do not copy old `.done`, pairs, matrices, or downstream outputs
into the new reference run.

## Chromosome Policy

The matrix schema is never silently renamed. The supported registry entries
expect UCSC names and retain:

- mm39: `chr1`-`chr19`, `chrX`, `chrY`, `chrM`;
- hg38: `chr1`-`chr22`, `chrX`, `chrY`, `chrM`.

`chrY` remains in both schemas even for expected XX samples. `1`, `MT`,
`chrMT`, NCBI accessions such as `NC_000001.11`, alt contigs, random contigs,
and unplaced contigs are not renamed. Preflight rejects a supported reference
whose canonical set differs from the expected UCSC set. If non-UCSC naming is
scientifically required, create a separate, explicitly named registry entry and
use a matching custom browser assembly; do not relabel an existing matrix.

## Used-Reference Preflight

Preflight validates only references selected by the current sheet. For each it
requires a readable FASTA, `.fai`, chromosome sizes, all expected BWA-MEM2
index components, exact FASTA-index/chromosome-size name and length agreement,
the expected canonical chromosomes, deliberate `chrM` naming, the correct
browser preset, and a matching readable phasing track when configured.

```bash
python scripts/reference_registry.py validate --registry config/references.tsv
bash scripts/preflight_check.sh -c config/config.conf -s config/samplesheet.csv
```

The first command checks registry structure. The second performs installed-file
checks for the references actually used.

## Downstream And Browser Behavior

The resolved phasing track is passed per sample/merge. If it is empty, the
workflow computes compartment eigenvectors without phasing and explicitly logs
that PC1 sign is arbitrary. Insulation, TAD-like calls, loops, BED, bedGraph,
`.mcool`, and `.hic` retain the matrix chromosome names.

Browser packaging requires `--reference-id`; `mm39` selects `mouse` and `hg38`
selects `human`. A contradictory explicit preset is rejected.
