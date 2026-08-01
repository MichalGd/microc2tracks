# Baseline Comparison

## Repository State

The local public-repository baseline was clean on branch `main` at commit
`3ac7607079fd45dc589ddf6d1fc56dace97ae54f` (`v2.0 snetinels,
paralelization, downstream analysis`, 2026-06-09). A read-only GitHub
`ls-remote` check confirmed that the public repository's `main`/HEAD was the
same commit at analysis time.

The configured local origin is
`https://github.com/MichalGd/micorc2tracks.git`, which contains the typo
`micorc2tracks`. The intended public URL is
`https://github.com/MichalGd/microc2tracks`. The baseline remote was not edited.

## File Comparison And Candidate Seeds

All published tracked files seeded the candidate. Matching production scripts
were then classified as follows:

- Identical production/repository files: `call_tads_from_insulation.py`,
  `compare_matrices.py`, `export_insulation_bedgraph.py`, `plot_ps_curves.py`,
  `preflight_check.sh`, `prepare_reference.sh`, `run_downstream.sh`,
  `sanitize_text_inputs.py`, `summarize_run.py`, and
  `validate_samplesheet.py`. Repository copies seeded these files.
- Newer production fixes: `microc2tracks.sh` and `merge_replicates.sh` differed
  only in four temporary gzip filenames that retain `.fastq.gz`/`.pairs.gz`
  suffixes so compression-aware tools recognize the format. Production copies
  seeded those two files before multi-reference edits.
- Production-only maintained utilities: `package_browser_tracks.py` and
  `compare_saddle_plots.py` were added to the candidate and documented.
- Timestamped `*.before_tmpfix_*` server backups are obsolete/generated
  operational artifacts and were not copied.

The repository supplied documentation, examples, tests, environment files,
packaging attributes, and public templates. Production configs supplied
behavioral evidence only: all were mm39; canonical filtering and preliminary
downstream analysis were enabled; production examples used `THREADS_FASTP=8`
and sample parallelism of one or two. Their server-specific paths and sample
FASTQs were not copied. The candidate defaults to `THREADS_FASTP=8` and
`MAX_PARALLEL_SAMPLES=2` as explicitly requested for this update.

## Difference Categories

- Functional workflow differences: gzip-suffix-safe temporary files and two
  production-only downstream/browser utilities.
- Production-specific configuration: installed reference, scratch, output,
  software, publication, and FASTQ paths; excluded from public candidate files.
- Documentation/packaging differences: repository-only and retained as the
  public baseline, then updated for multi-reference support.
- Potentially newer server fixes: the four temporary-name fixes, incorporated.
- Potentially newer public features: none beyond the clean commit; public and
  local HEAD were identical.
- Local uncommitted changes: none.
- Generated/obsolete files: `.git`, Git objects/logs, timestamped server
  backups, caches, staged distributions, and run outputs; excluded.
