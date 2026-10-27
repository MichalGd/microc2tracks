# Installed-server source synchronization

Historical note: the fingerprints below identify commit `8abf4cd`, before the
3.0 integration changes. They are not checksums of the newer candidate source.

Date: 2026-10-27

This small maintenance update synchronizes the public workflow source with the
installed-server fingerprints supplied by the maintainer. It is not the planned
Hi-C integration release.

## Source evidence and scope

The installation has no Git metadata. Its source was identified using supplied
SHA-256 hashes for all 18 enumerated Bash/Python scripts, `environment.yml`, and
`config/config_template.conf`. These 20 files are recorded in
[`SERVER_BASELINE_SHA256SUMS.txt`](../SERVER_BASELINE_SHA256SUMS.txt).

The saved 2026-08-01 package provides the matching source for 16 scripts and the
environment. Two one-line fixes reproduce the remaining script fingerprints.
The public configuration template already matches the supplied server hash.
This establishes byte-identical source for the supplied files, not live SSH
verification of the server, its current dependency versions, reference files,
run-specific configs, or analysis outputs.

The update is based on public `main` at
`85b95b64aecd7c93bc68b40b317b37f4d9779862`. Relative to that commit, it:

- Restores the installed runner's independent fastp slot limit and preflight
  validation/reporting. `MAX_PARALLEL_FASTP` defaults to `1` when absent from a
  config; explicit positive values can override it. Sample parallelism is
  unchanged.
- Writes normalized CSV/TSV using LF line endings. An empty trailing
  `merge_group` therefore cannot become a carriage-return directory name when
  Bash consumes the normalized TSV.
- Converts HTML cell values with the existing `as_text` helper before escaping,
  preventing the final-summary crash on `None` cells while retaining HTML
  escaping.
- Adds small regression tests and the source fingerprint manifest.

The merge script, scientific filtering defaults, sample schema, and matrix
generation are unchanged. This does not add the later all-sample pooling
extension, restriction-aware Hi-C classification, joint technical-library
deduplication, or alternative output roots. Those remain separate future work.

No server files, sample outputs, or running analyses are modified by publishing
this repository update. Existing incorrectly named output directories are not
automatically renamed or removed.

## Lightweight verification

```bash
sha256sum -c SERVER_BASELINE_SHA256SUMS.txt
python -m unittest discover -s tests -p 'test_*.py' -v
bash tests/check_bash_syntax.sh
bash tests/test_reference_sentinel.sh
```

These checks validate source identity and small regressions; they are not a
production-data rerun or proof of the installed runtime environment.

Results in the preparation workspace on 2026-10-27: all 20 source fingerprints
matched, all 18 Python tests passed, Bash syntax checks passed, and the reference
sentinel regression passed. Python checks used Thonny Python 3.14.4 on Windows;
shell checks used Git for Windows Bash 5.3.9. The server conda environment was
not exercised.
