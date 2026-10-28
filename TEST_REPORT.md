# Test Report

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
