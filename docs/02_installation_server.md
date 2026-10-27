# Shared Server Installation

Follow [the versioned multi-user installation guide](14_shared_multiuser_installation.md)
for the current candidate. It supersedes the earlier recipes that sourced conda
initialization from a maintainer's home, changed active permissions recursively,
or updated a shared checkout in place.

Code belongs under `/opt/microc2tracks/versions/<version>`, with a system-wide
launcher in `/usr/local/bin`. Code, environments and shared reference resources
must be readable/executable by ordinary users but maintained only by
administrators. Configs, run data and scratch belong to the user or project.

For this upgrade, reuse the existing `/opt/conda/envs/microc2tracks` and genome
indexes. Do not clone/update the environment, rebuild production indexes, change
active configs, or overwrite old code while an analysis runs. Stage new code
separately, verify access from another account, and run the tiny functional tests
when resources permit. Imports alone are not functional validation.

A fresh server without an environment needs a separate administrator-reviewed
dependency installation using `environment.yml`; that is not part of this
side-by-side upgrade. Reference paths in shipped examples are not instructions
to download/index genomes that are already installed.

See [inputs/configuration](03_inputs_and_config.md),
[references](09_multi_reference.md), [Hi-C/storage](12_integrated_hic_and_storage.md)
and [refilter mode](13_refilter_deduplicated_inputs.md).

Resource limits are per run, not a server-wide scheduler. Coordinate CPU/RAM/I/O
with other users. Keep system-disk reserves and use actual separate data/scratch
filesystems where available.
