# Test Report

## RC6 bounded-streaming candidate: 2026-11-04

Runtime: `3.0.0-rc6`; reviewed RC5 base:
`2ab2a8bcb27f82f8342ec60ab3c631433c962463`. Historical RC4/RC5 server results
below do not certify the changed RC6 runtime.

### Executed locally

- Windows regression discovery: **93 tests, 83 passed, 10 skipped**. The skipped
  tests were five full-tool tests, three opt-in pairtools/Cooler checks and two
  POSIX-only merge tests. Git Bash syntax checks passed for all shipped shell
  runners/examples. These are not server production-throughput measurements.
- Linux/WSL regression discovery with `MICROC2TRACKS_PAIR_TOOLS_TESTS=1` and
  `MICROC2TRACKS_REAL_TOOLS=0`: **93 tests, 88 passed, five skipped**. The five
  skips are the full FASTQ/BWA/fastp/Juicer tests, not failures or passes.
- Isolated Linux scientific environment: Python 3.12.3, pairtools 1.1.3,
  Cooler 0.10.4, NumPy 1.26.4, SciPy 1.17.1, pandas 3.0.3, h5py 3.16.0,
  pysam 0.23.3. Pairtools was built from its pinned v1.1.3 source tag using
  Cython 0.29.37. This local test environment did not change biolserv software.
- Real checked merge/GNU sort tests include eleven inputs (above the fan-in),
  full-line tie order, reader receipts, truncated gzip rejection and FIFO cleanup.
- Real pairtools streamed versus materialized merging/deduplication retained
  identical ordered data records, including the existing 10,000-row chunk bridge.
  Its two representatives remained unchanged; no dedup backend/chunk change.
- Real accepted-statistics aggregation matched a direct pooled pass.
- Real Cooler raw-count pooling matched direct integer pixels, excluded weight
  columns, reconciled accepted totals, checked convergence/NaN masks and finite
  weights (`rtol=1e-6`, `atol=1e-10`), and exercised dense-bin buffer overshoot.
  Extended checks cover 1/2/5-kb coarsening and int64-overflow rejection.
- Scalar/batched restriction and Micro-C filtering compared exact records and
  complete QC across seeded contacts, orientations, cutoffs, batch boundaries
  and reordered columns. Juicer projection matched the prior decompressed short
  format. These export tests do not certify a new end-to-end `.hic` build.
- Recovery tests cover removed FASTQs, removed sorted prerequisites with retained
  completed library outputs, scientific-contract rejection, exact RC5 key-formula
  migration without rewriting original payload/provenance, interrupted publication
  journals and shared-filesystem reservations. Existing ownership, relocation,
  cancellation and unchanged Micro-C regressions remain included.

Expected test diagnostics include a deliberately corrupted gzip reader, synthetic
addNorm failure, sparse-fixture pairtools division warnings and a deliberately
small Cooler buffer warning. They do not indicate a failed final test run.

### Bounded performance check and limits

The [recorded synthetic filter comparison](docs/rc6_filter_benchmark.json) uses
200,000 generated contacts (9,651,135 input bytes), two 100-Mb synthetic chromosome
axes and periodic restriction fragments. Scalar and batched output SHA256/QC are
identical. One run measured **1.435 s / 56.1 MiB peak process RSS** for scalar and
**1.344 s / 96.6 MiB** for batched. A smaller-digest fixture was slightly slower
when batched. These small shared-WSL kernel observations establish neither a
large filtering acceleration nor a whole-workflow speedup. Imports were warmed
equally and outputs were hashed rather than written as large compressed files.
Use `CONTACT_FILTER_ENGINE=scalar` when site measurements favor it.

The principal expected gains are eliminating redundant complete merged-file I/O,
pooled-pairs QC rereads/re-binning, and serial gzip export compression. Their
production wall time, peak scratch and process-tree RSS have not been benchmarked
on the 140-thread/500-GB server. Resource reservations and polling are estimates,
not OS-enforced CPU/RAM/disk quotas or proof of safety for every dataset.

### Required server gate before default promotion

Use the existing environment/JAR/references read-only, stage the complete new
version alongside RC5, check access as an unprivileged account, then run:

```bash
MICROC2TRACKS_PAIR_TOOLS_TESTS=1 python -B -m unittest discover -s tests -p test_rc6_real_science.py -v
MICROC2TRACKS_REAL_TOOLS=1 JUICER_TOOLS_JAR=/path/to/existing/juicer_tools.2.18.00.jar \
  python -B -m unittest discover -s tests -p test_real_tools.py -v
```

Those five tiny full-tool tests remain unexecuted for RC6 in this local setup.
No production FASTQs or live server outputs were used, removed or modified.
The GitHub tag denotes a release candidate, not automatic server installation.

## RC5 performance update preparation: 2026-10-29

Runtime candidate: `3.0.0-rc5`, implemented on the clean RC4 default commit
`d55761095d545353d9bc099a4c17e1ff2456ae95`. This section supersedes older
verification claims for the new candidate; the historical RC4 evidence below is
not being presented as a scientific-tool pass for changed RC5 orchestration.

- Full local regression discovery: **79 tests; 74 passed, five opt-in real-tool
  tests skipped** on Windows using
  the existing local Python and Git Bash. Skips are not passes.
- New focused checks cover independent ready alignment admission while two
  downstream jobs are blocked, CPU/RAM caps, prepared-input backlog, dependency
  readiness, oversized/cyclic task rejection, failure cancellation/credit release,
  subprocess cancellation preserving an unrelated process, one joint dedup per
  library, independent biological/library multiplicity, canonical accepted-file
  reuse, exact output inventory, final biological completion barrier and resume.
- Checked exact-fingerprint RC4 sorted checkpoint reuse without changing its
  original data/provenance/completion-record bytes or mtimes. Changed mapping
  settings, tool versions, input identities and failed markers reject reuse.
- Native-thread environment checks bound scientific Python worker pools while
  retaining BWA's own threading environment.
- Tested the large-server profile and singleton biological outputs. Existing
  filter, enzyme, reference, storage, refilter, output-scope, ownership and restart
  tests remain included. In-process mocked refilter stdin is serialized because
  actual production subprocesses have independent stdin streams.
- All shipped shell runners/examples pass Bash syntax checks. These tests do
  not establish production throughput or peak RSS.
- Five opt-in scientific tests now include fastp and a new two-biological-library
  final-only-hic fixture, alongside Micro-C, generic restriction Hi-C, original
  Dovetail Hi-C, and the pinned dedup chunk-boundary behavior. They compare raw
  Cooler pixels with Juicer observed/NONE contacts and test restart.

**Not yet verified:** RC5 with the actual Linux scientific binaries. This local
Windows preparation does not have the working biolserv environment. Run the
small opt-in suite after staging RC5 separately and before switching the default
launcher. No production FASTQs, shared environment, genomic indices, installed
release, server permissions or active analyses were changed during preparation.
See [the performance guide](docs/17_stage_scheduling_and_library_outputs.md).

## Linux server validation and default promotion: 2026-10-28

Evidence: user-supplied terminal output from `biolserv`, not direct remote
execution by the preparation agent. Tested runtime: `3.0.0-rc4`, commit
`22d10efa98a36234ab798db9096a071aff9a7c4b`, archive SHA-256
`de0f83a045457189d2f3581e6593ed518c670f139e9ba850342eb31b4e7f66c4`.
The following documentation-only promotion does not change that runtime.

- The server regression run discovered 67 tests: **63 passed, four opt-in tests
  skipped**, in 2.110 seconds. The separate opt-in run then completed **all four
  real-tool tests successfully** in 167.086 seconds. All 67 distinct checks have
  therefore passed on that server; skipped tests are not being counted as passes.
- Both maintainer and unrelated unprivileged (`nobody`, euid 65534) version routes
  worked. Code, shared executables, scientific imports, actual reference/index
  paths and the configured Juicer JAR were accessible in the checks provided.
  The configured reference-size table matched its existing `.fai` names/lengths.
- The scientific fixtures exercise Micro-C, generic restriction Hi-C, the original
  DpnII Dovetail profile, one joint dedup pass, independent-library retention and
  chunk-boundary behavior, merged-only `.hic`, split output roots, restart, and
  refiltering an existing jointly deduplicated stream without mapping/deduplication.
  Raw `.mcool` and `.hic` bin counts matched on those fixtures.
- No additional environment, package installations or production genome indexes
  were needed. Tests created their own tiny synthetic reference/index and temporary
  outputs, not production analyses.

Reported shared environment: Python 3.12.13, pairtools 1.1.3, cooler 0.10.4,
cooltools 0.7.1, bwa-mem2 2.3, pairix 0.3.9, htslib 1.23.1, numpy 1.26.4,
scipy 1.17.1, pandas 3.0.3, bioframe 0.8.0, h5py 3.16.0, OpenJDK 25.0.2 and
Juicer Tools 2.18.00. This is the observed working server environment, not a
claim that the fresh-install target in `environment.yml` was solved or tested.

### Warnings and limits

The bridge fixture is written in coordinate order but lacks a `#sorted` header.
Pairtools 1.1.3 emits the displayed sorting warning on a missing header marker;
it is not a detected record-order violation. See the
[upstream condition](https://github.com/open2c/pairtools/blob/v1.1.3/pairtools/cli/dedup.py#L464-L468).
The normal pipeline explicitly sorts before deduplication. The test also pins
pairtools' representative-only carryover: chunk-sensitive fuzzy clusters are
not claimed to have exact global transitive closure.

Java reflection/development-mode warnings were nonfatal in these tests; they
do not establish a measured performance penalty. The earlier intentional
synthetic addNorm failure was a regression test, not a production error.

The toy scientific runs disable fastp, preliminary downstream analyses and
MultiQC, use 1 kb matrices and VC normalization, and do not benchmark large
production runs or validate every enzyme/protocol/normalization combination.
The unit suite covers additional filtering/configuration contracts. No real
production library chemistry or dedup history is established by toy tests.
Project-specific configuration/preflight is still required, as is site access
review where another real login account has different ACL/group policies.

The software-validation gate has passed for this environment. The server's
default-launcher switch is a separate administrator action described in
[the promotion guide](docs/16_default_promotion.md); a GitHub push does not perform it.
The earlier sections below retain the historical preparation-stage results.

## Merged-only .hic candidate 3.0.0-rc4: 2026-10-28

- Full local suite: **67 discovered, 63 passed, four real-tool tests skipped**.
- Six new tests check default/override/master-off policies, invalid scope rejection,
  all requested technical/all-sample exports with no individual Juicer/export/link
  work, unchanged per-row cooler/pairs, no-Java singleton runs and explicit report
  skips, scope-change reuse/preservation, and failed merged normalization restart.
- The synthetic addNorm-failure test deliberately exercises failure reporting;
  its expected diagnostic is not a real scientific-tool failure.
- Bash syntax, Python AST syntax, manifest completeness and Git whitespace checks
  passed. No filtering, deduplication or cooler algorithms were changed.
- The existing opt-in Linux tests now assert that individual/internal-dedup `.hic`
  stages are absent by default; the final already-merged refilter fixture explicitly
  selects `all`. Those four real-tool checks remain unexecuted here, not passed.

The shared outside-home installation model and actual cross-account access gate
remain in force. No server files/jobs/packages/indexes or production outputs were
modified. See [the output-scope guide](docs/15_merged_hic_outputs.md).

## Shared-installation candidate 3.0.0-rc3: 2026-10-27

- Full local suite: **61 discovered, 57 passed, four real-tool tests skipped**.
- Five new tests cover shared-prefix selection without conda activation,
  run/preflight/merge/check-access routing, argument/working-directory preservation,
  user-site/PYTHONPATH isolation, missing-Python failure, no installed-source writes
  on the actual version route, and read-only reference/JAR permission checks.
- Bash syntax, reference-sentinel regressions and Python AST syntax checks passed.
  Source-tree identity snapshots confirm no bytecode/new source artifacts from the
  launcher. Tests used Windows/Git Bash, not another actual Linux login identity.
- No scientific algorithms were changed by the shared-installation update.
  The four existing toy scientific tests remain unexecuted here. Cross-account
  access on the real server has not been tested or claimed.

For this server upgrade, **reuse the existing environment** after checking it;
do not follow historical environment-creation commands below against an active
environment. Existing imports/version checks are encouraging but not end-to-end
validation. See [shared deployment](docs/14_shared_multiuser_installation.md)
for root-owned code staging outside user homes and other-account access checks.

After reviewing resource availability, run the opt-in real-tool tests with the
existing environment and actual shared JAR path, not production data. They create
only a tiny temporary synthetic reference/index. Any failure is a release blocker
to investigate, not permission to update packages or modify active jobs.

## Refilter mode candidate 3.0.0-rc2: 2026-10-27

Extends local candidate `d158d3f` with the agreed already-deduplicated input mode.
No production analysis or server deployment was performed.

- Local suite: **56 discovered, 52 passed, four real-tool tests skipped**.
- Ten new tests cover no FASTQ/BWA requirement, skipped mapping/dedup commands,
  matching reference/schema checks, end swapping and stale-fragment removal,
  preserved multiplicity, correct restriction/cis boundaries, new-root safeguards,
  source/index immutability, overlapping pools, restart and ignored unused settings.
- The mocked orchestration contract rebuilds both output branches and proves no
  dedup command is invoked; this is not real scientific-tool validation.
- The generic opt-in Hi-C smoke test now additionally imports its jointly
  deduplicated pre-selection pairs, refilters/rebuilds products and compares raw
  bin counts while checking source preservation. It has not been run here.

The existing Linux environment/real-tool release gate below still applies. Use a
separate `microc2tracks-3.0.0-rc2` environment when running those checks. Source
input deduplication is an explicit user assertion, not verified or repaired by
this mode. See [the mode guide](docs/13_refilter_deduplicated_inputs.md).

## Integrated update candidate 3.0.0-rc1: 2026-10-27

Baseline: published `main` commit `0ce6761` (the verified installed-server source
sync and its changelog). All agreed feature areas are implemented in one local
candidate; no server installation, environment or production outputs were changed.

### Executed checks

- `python -m unittest discover -s tests -p 'test_*.py' -q`: **46 discovered,
  42 passed, 4 explicitly skipped**; about seven seconds in this workspace.
- Python compilation of `scripts/` and `tests/`: **passed**.
- `bash tests/check_bash_syntax.sh`: **passed**.
- `bash tests/test_reference_sentinel.sh`: **passed** (historical helper regression).
- Main, merge and preflight shell-entry help routes: **passed**.
- Browser-package CLI help and example CSV validation: **passed**. Missing
  example FASTQs warn as expected; examples are not production run instructions.
- Git whitespace validation: **passed**.

The 42 passing tests include exact enzyme cut maps/offsets and chunk overlaps,
fragment/orientation classes, chromosome/coordinate checks, disabled/inclusive
distance behavior, Micro-C legacy-predicate contact equivalence, explicit profiles,
disjoint/partial dedup units, arbitrary custom contigs, CRLF/empty metadata,
path fallback/overrides, selective restart signatures, locks, reserves/mounts,
interrupted bundles, failed checksums, owned-only cleanup and legacy-output refusal.

A **mocked external-tool orchestration contract**, explicitly not a scientific
tool test, exercises row/joint/pool scheduling, independent-unit retention,
resume without new commands, verified disk relocation without recomputation,
and standalone pool/joint modes without replacing full-run metadata.

Environment used here: Windows, Thonny Python 3.14.4, Git for Windows Bash 5.3.9.
The target environment remains separate Python 3.11 with the scientific version
targets in `environment.yml`. That Linux environment has not been solved/tested here.

### Required remaining release gate

The four real-tool tests in `tests/test_real_tools.py` were **not run**. This
workspace has no pairtools/cooler/BWA/SciPy installation; WSL has no usable Linux
distribution and Docker's Linux engine is not running, confirmed outside the
sandbox as well. Do not interpret skipped tests as passing real-tool validation.

The small opt-in suite covers Micro-C, generic restriction Hi-C and original-kit
Dovetail FASTQ-to-`.mcool`/`.hic` runs, raw binned-count equality, joint scope,
split output roots, restart, and pinned fuzzy bridge/chunk behavior. Run it in an
isolated target environment before production promotion:

```bash
conda env create -n microc2tracks-3.0.0-rc1 -f environment.yml
conda activate microc2tracks-3.0.0-rc1
MICROC2TRACKS_REAL_TOOLS=1 JUICER_TOOLS_JAR=/path/to/juicer_tools.2.18.00.jar \
  python -m unittest discover -s tests -p test_real_tools.py -v
```

No production-data rerun or extensive benchmark is required. The historical
server fingerprint manifest below describes the earlier synchronized commit,
not the changed integration-candidate files.

## Installed-server source sync: 2026-10-27

The [source synchronization note](docs/11_server_baseline_sync.md) records scope,
commands, and limitations. Verification in the Windows preparation workspace:

- Supplied installed-server fingerprints: **20/20 exact matches**.
- Python synthetic/regression tests: **18/18 passed**, including the three new
  normalized-line-ending and HTML-cell regression tests.
- Bash syntax checks and reference-sentinel regression: **passed**.
- Git whitespace check: **passed**.
- No production-data rerun or server-side environment test was performed.

The report below is the retained historical 2026-08-01 validation, not a new
production validation.

Test date: 2026-08-01 (Europe/Berlin)

## Environment

Candidate tests ran in the Windows preparation workspace. WSL's system `bash`
service was disabled, so Bash checks used Git for Windows Bash. The Windows
Store `python` alias was unavailable in the sandbox, so Python tests used Thonny
Python 3.14.4. The deployed conda environment remains pinned to Python 3.11.

## Commands And Results

Repository/baseline inspection (pass):

```text
git -c safe.directory=<local-baseline> -C microc2tracks remote -v
git -c safe.directory=<local-baseline> -C microc2tracks status --short
git -c safe.directory=<local-baseline> -C microc2tracks branch --show-current
git -c safe.directory=<local-baseline> -C microc2tracks log -1 --oneline
git ls-remote https://github.com/MichalGd/microc2tracks.git HEAD refs/heads/main refs/tags/v2.0
```

The local tree was clean on `main`; local and public `main`/HEAD both resolved
to `3ac7607079fd45dc589ddf6d1fc56dace97ae54f`.

Bash syntax and reference sentinel regression (pass):

```bash
bash tests/check_bash_syntax.sh
bash tests/test_reference_sentinel.sh
```

All runner, merge, preflight, reference helper, downstream, preparation, and
example shell scripts passed `bash -n`. A sentinel written for mm39 was rejected
for hg38 and vice versa.

Python compilation (pass):

```bash
python -m py_compile scripts/*.py tests/test_multireference.py
```

Registry/sample-sheet smoke tests (pass; missing example FASTQs warn as
expected):

```bash
python scripts/reference_registry.py validate --registry config/references.tsv
python scripts/validate_samplesheet.py \
  --reference-registry config/references.tsv \
  config/samplesheet_example.csv
python scripts/package_browser_tracks.py --help
```

The example resolved both `mm39` and `hg38`. A separate negative smoke test
confirmed that `--reference-id hg38 --canonical-preset mouse` exits with an
error.

Mouse/backward-compatibility regression first (4/4 pass):

```bash
python -m unittest -v \
  tests.test_multireference.MultiReferenceTests.test_alias_normalization_and_browser_presets \
  tests.test_multireference.MultiReferenceTests.test_missing_reference_column_defaults_to_mouse_with_warning \
  tests.test_multireference.MultiReferenceTests.test_expected_canonical_mouse_and_human_chromosomes \
  tests.test_multireference.MultiReferenceTests.test_registry_file_preflight_validates_both_references
```

Human/mixed/merge regression (5/5 pass):

```bash
python -m unittest -v \
  tests.test_multireference.MultiReferenceTests.test_valid_mouse_human_and_mixed_sheets \
  tests.test_multireference.MultiReferenceTests.test_cross_reference_replicate_merge_is_rejected \
  tests.test_multireference.MultiReferenceTests.test_explicit_merge_group_requires_compatible_metadata \
  tests.test_multireference.MultiReferenceTests.test_mitochondrial_mt_and_chrmt_are_not_silently_renamed \
  tests.test_multireference.MultiReferenceTests.test_reference_is_carried_into_final_report_rows
```

Complete synthetic suite (15/15 pass):

```bash
python -m unittest discover -s tests -p 'test_*.py' -v
```

Requested resource-default check (pass): the final template contains exactly
`THREADS_FASTP="8"` and `MAX_PARALLEL_SAMPLES="2"`; Bash syntax and the complete
15-test suite were rerun after this update.

Coverage includes `mouse`/`human` aliases, explicit `mm39`/`hg38`, unknown
references, missing reference column/default warning, mixed rows, cross-reference
and incompatible explicit merge groups, exact mouse/human canonical sets,
`chrM`/`MT`/`chrMT`, FASTA/chromosome-size mismatch, missing BWA-MEM2 component,
human browser preset, pairs assembly header, report reference fields, stale
sentinels, and alignment-FASTQ read-name mismatch.

Diff hygiene (pass):

```bash
git diff --no-index --check -- <exported-baseline> microc2tracks_updated
```

No whitespace errors were reported. Candidate text files were also checked as
LF-only, and generated `__pycache__`/`.pyc` files were removed.

## Skipped Tests

- Full `preflight_check.sh` against installed mm39/hg38 references: skipped in
  the preparation workspace because complete references and the Linux
  bioinformatics toolchain are not present.
- End-to-end fastp/BWA-MEM2/Pairtools/Cooler/Juicer run with real FASTQs: skipped
  for the same reason and because the production server was explicitly out of
  scope.
- Real-data mouse matrix comparison and human pilot run: deployment-gate tests,
  documented in `SERVER_DEPLOYMENT_NOTES.md`; not replaced by synthetic claims.

## Failures And Risks

No routine test failures remain. Synthetic tests validate control flow and
metadata but cannot establish biological equivalence, tool/runtime performance,
or correctness of server-installed reference files. Those remain explicit
pre-deployment gates.
