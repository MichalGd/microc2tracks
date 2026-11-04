# Stage scheduling and library/biological outputs (RC5)

RC5 is one integrated Micro-C/Hi-C performance update. No scientific package or
production genome index needs to be installed or duplicated. It preserves the
existing assay/reference/filter/resolution/normalization configuration.

## Data and scheduling model

Each FASTQ run has file-based fastp followed by the existing BWA/parse/sort/bgzip
pipeline and an original sorted-pairs checkpoint. Physical FASTQ concatenation
and raw-stdin fastp are not introduced. This preserves run-specific automatic
preprocessing and restart boundaries.

When `BUILD_SAMPLE_MATRICES=false` and `RUN_GROUPED_ROW_CONTACT_QC=false`, rows
assigned to `dedup_group` stop at those original sorted pairs. One complete joint
dedup invocation per library consumes the merged original streams, then applies
the existing shared contact filter. Standalone ungrouped rows retain normal
dedup/filter behavior. Enabling individual matrices or the contact-QC switch
explicitly restores individual contact processing.

Tasks are admitted only when dependencies and CPU/RAM/stage slots are available.
They do not wait inside a small whole-sample worker pool. A ready library need
not wait for unrelated FASTQs. All scientific subprocess errors are checked;
failure/SIGTERM sets the shared stop event and stops only that run's process
groups. Retained complete checkpoints remain reusable; incomplete stages restart.

## One CSV and one project config

Keep using `dedup_group` as an explicit physical-library identity. Never infer it
from a biological label. A library cannot span independent biological groups.
`merge_group` identifies a biological output and can combine several independent
libraries without another dedup pass. Distinguish biological labels using the
existing reference/assay/condition context and explicit safe output IDs.

Optional `final_pool` is the desired final output ID, repeated on every run of
each selected biological group. Leave it empty to exclude that group. Partial
biological/library selection and conflicting explicit/legacy pooling requests
fail. Use either explicit membership or `RUN_MERGE_ALL_SAMPLES=true`, not both.
Never automatically pool unrelated conditions. Membership is also supported by
the separate refilter input CSV without changing its read-only scientific path.

Recommended output switches:

```bash
BUILD_SAMPLE_MATRICES="false"
BUILD_BIOLOGICAL_MATRICES="true"
RUN_GROUPED_ROW_CONTACT_QC="false"
RUN_MERGE_TECHNICAL_REPLICATES="true"
RUN_MERGE_ALL_SAMPLES="false"
RUN_MCOOL="true"
RUN_HIC="true"
HIC_OUTPUT_SCOPE="final"
```

This produces biological `.mcool` outputs, including singleton biological groups,
and final `.mcool` plus `.hic`. Final pooling waits for selected biological output
completion (including requested downstream work, retaining the existing optional
failure policy). It merges disjoint accepted library contacts without dedup or
new filtering and normalizes each format anew. Do not sum balanced matrices or
their weights. Raw Cooler counts and Juicer observed/NONE counts should match
at shared resolutions; their normalized values need not match.

A biological output with one library references its canonical accepted pairs,
index, chromosome sizes and QC instead of copying/reindexing them. Manifests
record that source. Do not manually move/delete shared artifacts. Older hic or
individual outputs already present are preserved and marked not requested, not
silently deleted because output scope changed.

## Server resource profile

Copy settings from `config/server_140cpu_500gb.conf` into the project config, or
source that absolute installed profile path after ordinary project settings.
It does not supply/change project, FASTQ, reference, scratch or JAR paths.
Keep using a writable project config, not edits to shared read-only software.

The profile requests two 32-thread BWA alignments. Each alignment reservation
also accounts for parsing, four sort threads, four compression threads and
helper allowance: 44 CPU slots for the whole pipeline. Two concurrent fastp
jobs use eight threads each. The 120-CPU / 400-GiB admission budget prevents
all stage maxima from being active together. Available RAM is also checked
before admission, leaving 16 GiB headroom in this profile.

`MAX_PARALLEL_SAMPLES` is retained only as an accepted legacy setting; it no
longer defines alignment slots. Use `MAX_PARALLEL_ALIGNMENTS`. Other independent
caps include FASTP, MERGES, DEDUP, FILTER, MATRIX, HIC, DOWNSTREAM and REPORTS.
`THREADS_BALANCE` is passed separately to Cooler balancing. More BWA threads do
not make the current scipy dedup kernel a 32-thread operation.
Scientific Python/downstream subprocesses cap hidden BLAS/OpenMP/Numba worker
threads at one per process, so process-level parallelism is not multiplied by
another automatic native thread pool. BWA/fastp/Java retain their configured
threading environment and explicit command-line allocations.

The effective RAM admission budget is capped at 80% of currently available RAM
at startup, even if the requested ceiling is higher; preflight reports both.
This avoids reserving the nominal 400 GiB when other users already occupy it.
Budgets/reservations are submitting-user admission estimates, not hard RSS or
host-wide limits. Separate runs/users need coordinated allocations. Affinity
and cgroup-v2 CPU quota are respected; a requested CPU budget exceeding the
available allowance fails rather than silently assuming all host CPUs are free.
`TOTAL_*` zero defaults use affinity/quota CPUs minus one and 80% of available
RAM. Portable configurations must fit the actual machine; the large profile
is not automatically applied to every server.

`MAX_PREPARED_SAMPLES` bounds running fastp, waiting prepared runs and active
alignments together. Sort scratch belongs on the chosen fast device with free
space reserves; retained contacts/final files can remain on larger storage.
Existing optional roots and OUTDIR fallbacks remain unchanged. This bound and
space floors are not exact predictions of temporary bytes for a very large run.

New disposable trims are released after their sorted checkpoint succeeds.
`KEEP_SORTED_PAIRS` controls original per-run streams. `KEEP_RAW_MERGED_PAIRS`
separately controls newly generated library merged originals, which can be
released after dedup succeeds. Cleanup only removes explicitly owned files
generated in this run; prior runs' retained artifacts are preserved.

## Visibility

Preflight prints effective global budgets, independent limits, reservations,
final hic targets and grouped rows that skip individual contact processing.
Run metadata includes `resource_plan.json`, `stage_status.json` and append-only
`stage_events.jsonl`. Waiting reasons distinguish dependencies, stage limits,
CPU budget, RAM budget, currently available RAM and prepared-input backlog.
Reports mark individual duplicate/contact metrics not computed when omitted.
An idle CPU interval can still be legitimate when input readiness or storage
throughput limits progress; no fixed speedup is promised.

## RC4 checkpoint reuse and handover

Do not replace an installation used by an active run. Stage RC5 as a separate
root-owned `/opt/microc2tracks/versions/3.0.0-rc5` release with a version-specific
launcher. Keep the same existing environment/indexes/JAR. Ordinary user accounts
need read/execute access to code/tools/reference paths and writable private
project/output/scratch/cache paths; never use chmod 777 or run analyses as root.

`REUSE_RC4_OUTPUTS=true` enables only the exact source fingerprint allowlisted
in `rc4_compatibility.py`. RC4 source metadata is preserved separately before
the new run metadata is written. Candidate old completion keys are reconstructed
from the current scientific inputs, settings and tool versions plus that known
old contract. Stage reuse still validates successful records and all physical
payload identities. Original RC4 data, provenance and completion records are not
rewritten just to claim reuse. Unknown source versions are not certified.

Test a staged release first. A later controlled handover must stop the old
controller, wait for its children and release locks before starting RC5 against
the same output tree. Do not run both versions concurrently on the same outputs.
Keep FASTQ/reference/assay/dedup/filter settings unchanged for reuse. Partial
stages may restart; preserved sorted checkpoints avoid repeating valid mapping.
Previously independently deduplicated representatives cannot replace original
sorted streams when joint library dedup has not yet been done. Previously
accepted contacts cannot restore pairs removed by old filters.

## Minimum verification before promotion

Run the ordinary regression suite, then on Linux with the existing shared tools:

```bash
cd /opt/microc2tracks/versions/3.0.0-rc5
env -u PYTHONPATH -u PYTHONHOME \
  PATH="/opt/conda/envs/microc2tracks/bin:$PATH" \
  PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 \
  /opt/conda/envs/microc2tracks/bin/python -B \
  -m unittest discover -s tests -p 'test_*.py' -q

env -u PYTHONPATH -u PYTHONHOME \
  PATH="/opt/conda/envs/microc2tracks/bin:$PATH" \
  PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 \
  MICROC2TRACKS_REAL_TOOLS=1 JUICER_TOOLS_JAR="/actual/readable/juicer_tools.jar" \
  /opt/conda/envs/microc2tracks/bin/python -B \
  -m unittest discover -s tests -p 'test_real_tools.py' -v
```

The five opt-in tests use tiny temporary references/reads, now include fastp and
the two-biological-library/final-only-hic contract, and do not touch production
indexes/data. No production reruns or extensive benchmark are required. Existing
RC4 server test results do not certify changed RC5 orchestration. Promote the
default launcher only after these checks; archive preparation is not promotion.
