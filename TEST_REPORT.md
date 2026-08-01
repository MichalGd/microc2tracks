# Test Report

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
