# Integrated update handoff: 3.0.0-rc2

Prepared: 2026-10-27. Base: `0ce6761` (published server-source synchronization).
Local branch: `feature/integrated-hic-storage-20261027`.

## Implemented together

- One integrated Micro-C/Hi-C workflow with the existing shell commands.
- Four supported restriction enzymes and combinations; explicit generic,
  original-kit Dovetail, and legacy profiles; fragment and independent distance
  filtering shared by individual/joint products and both matrix formats.
- User-selected joint deduplication of original sorted streams; individual
  products retained; technical/all-sample pools use disjoint units once.
- Optional disk roots with OUTDIR fallback, isolated scratch, checked publication,
  parameter-aware restart, safe managed relocation, locks and owned cleanup.
- Path-aware reports/browser packaging, generic examples, target version pins,
  and small scientific/contract tests.
- Explicit refilter-deduplicated run mode for existing mapped pairs: no alignment,
  parsing or deduplication, read-only sources/indexes, new output tree, import QC
  and the shared restriction/distance/matrix path. See [the refilter guide](docs/13_refilter_deduplicated_inputs.md).

## Verification and release gate

52 local tests passed; four opt-in real-tool tests were skipped. Python
compilation, Bash syntax/sentinel checks, shell help routes, example CSV validation,
browser help and Git whitespace checks passed. The runner contract includes
mocked external commands, not proof of real BWA/pairtools/cooler/Java behavior.

Before promotion, solve the isolated Linux target environment and run the four
tiny real-tool tests documented in [TEST_REPORT.md](TEST_REPORT.md). This checks
FASTQ-to-matrix behavior and raw `.mcool`/`.hic` equality without production data
or extensive benchmarking. Local Docker's Linux engine was unavailable, so those
checks could not be truthfully marked passed here.

No server deployment, running-job changes, production-output deletion, or
production-data reruns were performed. GitHub main was not changed by preparing
this candidate. Original planning documents and prior baseline directories were
preserved. Do not overwrite the active server installation with an unvalidated
candidate.

See [the user guide](docs/12_integrated_hic_and_storage.md),
[configuration template](config/config_template.conf), and
[changelog](CHANGELOG_CANDIDATE.md).
