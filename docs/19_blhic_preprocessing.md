# Integrated BL-Hi-C: version 4.0 (implemented in RC7)

The RC7 implementation is published under version 4.0. The linker preset,
contact processing and scientific stage contracts are unchanged. Use the
[4.0 deployment guide](20_version_4_0_deployment.md) for installation.

Select `assay=hic` and `hic_protocol=blhic` in the samplesheet, or set
`HIC_PROTOCOL="blhic"` for a homogeneous run. Point `fastq_r1`/`fastq_r2` at
the original verified FASTQs. The existing workflow command now performs
bridge selection, optional fastp, mapping, contact selection and matrices.

This preset supports **HaeIII, two-step bridge-ligation BL-Hi-C**. It is a
common reanalysis, not exact historical HiCpipe/HiC-Pro reproduction. It does
not implement other BL chemistries or reproduce paper-specific loop callers.
Arima GANTC/N support remains separate work.

## Set up the optional tool once

Ordinary Hi-C/Micro-C and refilter-only runs need no new executable. Raw BL
runs need only the standalone pinned `trimLinker` component, not HiC-Pro or
the complete historical ChIA-PET2 stack. Its source/license are bundled.
An administrator can build it before making a versioned release read-only:

```bash
cd /path/to/new/microc2tracks
python3 -B scripts/build_trimlinker.py
```

This uses an existing C++ compiler, OpenMP and zlib headers; it makes no
network request or system installation. The default destination is
`tools/blhic/trimLinker`, accompanied by a required `.build.json` receipt.
Alternatively build with `--output /path/to/versioned/tools/trimLinker`
and set `BLHIC_TRIMLINKER` to that absolute executable. Keep its receipt.
Do not update the environment or launcher of an active job in place.

Source revision: ChIA-PET2 `e120726d6440b24034f70bc3c51c17f351fef496`.
The source manifest checks the bundled files. Preflight verifies the build
receipt and executable SHA256; no unsupported `trimLinker --version` is used.
If the optional tool is missing or unrecognized, the BL preflight fails
clearly rather than installing it or substituting another algorithm.

## Configuration and samplesheet

Copy the normal template and choose:

```bash
HIC_PROTOCOL="blhic"
HIC_RESTRICTION_ENZYMES="HaeIII"
BLHIC_INPUT_STATE="raw"
RUN_FASTP="true"             # independent optional quality/adapter cleanup
MAX_PARALLEL_FASTP="1"       # includes the bridge stage
THREADS_FASTP="8"            # bridge and fastp run sequentially
MEMORY_BLHIC_GB="16"         # enforced minimum for raw BL tasks
TMPDIR="/absolute/bulk/scratch/project"
TRIMMED_OUTDIR="/absolute/bulk/processed_reads/project"
```

See [the example config](../config/config_blhic_example.conf) and
[samplesheet](../config/samplesheet_blhic_example.csv). Replace their paths.
The normal `microc2tracks.sh -c CONFIG -s CSV` and preflight commands apply.

The profile defaults to HaeIII if no enzyme is supplied. Explicit row values
override global values. A copied global DpnII setting is rejected for BL;
override it with row `restriction_enzymes=HaeIII` or change the global setting.
The fragment map uses GGCC, cut offset 2. Keep `assay=hic`; `assay=blhic` is
not an accepted assay.

For mixed inputs, optional row `blhic_input_state` overrides
`BLHIC_INPUT_STATE`. `raw` runs bridge selection. `pretrimmed` explicitly
declares prior bridge processing and skips only that step. The workflow
validates supplied pairs and records user-declared provenance; without a
historical trusted receipt, original archive counts/bridge retention remain
unknown. No automatic bypass is based on low linker frequency.

`RUN_FASTP=false` **still performs required BL bridge processing**, validates
the retained pair and compresses it internally. With fastp enabled, bridge
outputs feed fastp first, and the final paired mapping inputs are validated.
Read merging is rejected. fastp may apply additional quality/length filters;
its losses are separate from bridge retention.

## Fixed scientific preset

Internal ID: `haeiii_two_step_v1`.

- A: `ACGCGATATCTTATC`; B: `AGTCAGATAAGATAT` (15-base search motifs).
- Bridge mode `-m 1`; keep-empty `-k 2`; substitution allowance `-e 1`;
  minimum retained length `-l 15` in both mates.
- Exact match preferred over fuzzy match for each motif. Matching is
  case-sensitive and substitution-only; no additional reverse-complement or
  indel search is silently added.
- Retain genomic sequence/quality before the linker. Opposite A/B is eligible;
  same A/A or B/B is rejected. One/both mates without a linker are eligible
  under `-k 2`; ambiguous assignments and too-short pairs are rejected.
- No post-linker segment is emitted as an additional contact.

The motifs are not the entire modified 20-base bridge. They and these options
are supported by the [original main methods](https://www.nature.com/articles/s41467-017-01754-3)
and [HiCpipe invocation](https://github.com/ChenFengling/HiCpipe/blob/af5e9089df91de4eb9e9da5cdcf37b2806a04e4e/main.sh).
No custom matching knobs are exposed in this first preset.

## Integrity, resources and storage

The unchanged upstream reader has weak error checking and assumes a particular
FASTQ format. The wrapper therefore validates the entire input pair, including
gzip CRC/EOF, equal counts/normalized identifiers, four complete LF lines,
sequence/quality lengths and uppercase A/C/G/T/N. It rejects lines exceeding
999 bytes including LF, CRLF, NULs and malformed records before invoking the
upstream reader. Ordinary assay validation retains its previous behavior.

Input validation receipts are reused only for matching file path/size/mtime
identities, following the workflow's existing input identity convention.
Original archive MD5 verification remains the downloader's responsibility.
Do not change inputs during a run. Output records and linker statistics must
reconcile; exit zero or filename existence alone cannot mark a stage complete.

The trimmer holds two million records per mate per batch. Raw BL tasks reserve
at least **16 GiB**, conservatively based on bounded line/record allocations;
this is not a measured peak for the user's libraries. Ordinary jobs in a mixed
run keep their own reservation. A production resource pilot must fill that
batch or account for it analytically; a tiny semantic fixture is insufficient.

The trimmer emits **plain FASTQ**. `TMPDIR`, not `TRIMMED_OUTDIR` alone, controls
its owned attempt directory. Use bulk scratch with adequate space. After the
strict scan measures uncompressed input bytes, the stage grows its filesystem
reservation for the plain pair plus conservative compressed-output/publication
copies before creating large outputs. It fails clearly if that bound cannot be
admitted or exceeds `SCRATCH_BUDGET_GB`. A failure can be retried after selecting
larger storage; no original input is removed. Existing minimum-free-space and
process-tree/scratch guards remain active. Resource reservations are not quotas
and do not coordinate independent jobs submitted by other users.

The wrapper removes only its successful attempt's named plain intermediates.
Final processed reads follow `KEEP_TRIMMED_FASTQ`; compact linker QC and receipts
remain. Incomplete attempts/logs are retained for diagnosis.

## Outputs, QC and restart

Final mapping inputs use the existing `02_trimmed/*.trim.R[12].fastq.gz` paths.
The existing trim transaction also owns:

- `01_qc/SAMPLE.blhic.trim.stat`: original linker statistics, or an explicit
  pretrimmed/unknown-history note;
- `01_qc/SAMPLE.blhic.json`: supplied-input validation, preset, tool/build identity,
  command, retention/length statistics and final mapping-input identities;
- normal fastp JSON/HTML when enabled.

The raw-input count equals `Total PETs`. `Expect` equals opposite-linker plus
one-empty plus both-empty classifications. `Valid` is the retained count after
both-mate length filtering. Ambiguous and too-short counts are derived from
reconciled totals; an unused racy upstream counter is not reported. Final
summaries distinguish original pairs, supplied pairs, bridge-retained pairs and
mapping-input pairs. Pretrimmed inputs do not receive invented raw statistics.

BL trim signatures include the fixed preset, declared state and pinned binary
identity. Changing chemistry invalidates affected read preparation/mapping;
changing only contact MAPQ/distance does not. Moving an identical tool preserves
the scientific signature. Ordinary Hi-C/Micro-C stage contracts remain RC6
compatible. Missing-raw checkpoint recovery additionally checks BL lineage;
an old untrimmed mapping cannot become BL-processed by changing its label.

Sorted/deduplicated pairs retain the preprocessing receipt reference. Refilter
mode preserves available upstream provenance but neither requires the trimmer
nor repairs historical missing bridge processing.

## Validation and production scope

`tests/test_blhic.py` exercises the real pinned tool on tiny synchronized reads
at one and two threads, plus configuration, orchestration, failure/publication,
QC, storage and restart cases. Global output order can differ with OpenMP;
compare paired content, not byte-identical compressed files across thread counts.

Before full GSE157339 processing, use one bounded normal-workflow pilot to confirm
archive read state, retention, mapping/contact QC and resources. This release
does not claim a real GSE157339 production pilot was run. For comparisons with
Bonev, align raw data to the same exact mm39 reference, retain independent
libraries and compatible matrix masks/bins/normalization. BL chemistry and
study/biology remain potential confounders. HiCCUPS/Arrowhead/HiCDB and differential
analyses are separate downstream choices.
