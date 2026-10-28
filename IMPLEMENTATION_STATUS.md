# Integrated workflow status: 3.0.0-rc4

Prepared: 2026-10-28. Base: `0ce6761` (published server-source synchronization).
Runtime baseline: `22d10efa98a36234ab798db9096a071aff9a7c4b`.
The published default retains the exact RC4 runtime tested on the server.

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
- Read-only shared installation outside user homes, a version-resolving launcher
  using the existing environment without conda activation, and explicit
  current-account access checks. See [the multi-user guide](docs/14_shared_multiuser_installation.md).
- Merged-only `.hic` generation by default, retaining individual cooler matrices;
  explicit all-output override, conditional Java/JAR requirements, exact preflight
  targets and report/restart safeguards. See [the scope guide](docs/15_merged_hic_outputs.md).

## Verification and default promotion

Preparation checks passed locally, and the user subsequently supplied Linux
server results: 63 regressions passed, followed by all four opt-in real-tool
checks (167.086 seconds). The exact scope, working environment, warnings and
limitations are recorded in [TEST_REPORT.md](TEST_REPORT.md). Unlike the earlier
mocked-only status, these results include real BWA/pairtools/cooler/Java runs.

The server software-validation gate has passed using the existing environment.
No duplicate software installation or production reindexing was required. These
small tests are not exhaustive production validation: actual chemistry, dedup
scope, reference policy, resource limits and a new output tree still need project
configuration/preflight. The fresh-install target environment has not been solved
here and is not an instruction to replace the working server environment.

The user installed RC4 outside the home directory and tested it from maintainer
and unrelated unprivileged identities. The preparation agent did not remotely
deploy or run tests itself. No production-output deletion or production rerun
was performed. Original planning documents and archived candidates are preserved.

The user has confirmed workflow jobs are idle and requested default promotion.
Follow [the promotion guide](docs/16_default_promotion.md) to back up existing
launchers and switch the three standard commands. The GitHub update is separate
from this server action. Keep the versioned launcher, old code and environment;
do not recursively change existing resource permissions. Other real login
accounts may still need checks where site ACLs differ.

See [the user guide](docs/12_integrated_hic_and_storage.md),
[configuration template](config/config_template.conf), and
[changelog](CHANGELOG_CANDIDATE.md).
