# microc2tracks

`microc2tracks` is a shell-script-first Micro-C / Hi-C workflow for producing common contact matrix outputs and downstream analyses from paired-end FASTQ files.

The first implementation follows the style of the existing Shell repositories `fastq2tracks` and `rnaseq2tracksP`: a top-level runner, `config/`, `docs/`, `examples/`, `scripts/`, `environment.yml`, and simple tests.

Version **3.0.0-rc3** implements integrated restriction Hi-C, user-selected joint
deduplication, disjoint pooling and optional output disks. Read the
[integration/storage guide](docs/12_integrated_hic_and_storage.md) before upgrading.
It preserves the [earlier server fixes](docs/11_server_baseline_sync.md).
This candidate still requires the tiny real-tool Linux smoke tests before
production promotion; see [the test report](TEST_REPORT.md).

Existing mapped, deduplicated `.pairs`/`.pairs.gz` can be reprocessed with
`--mode refilter-deduplicated`, skipping alignment and deduplication. See the
[refilter mode guide](docs/13_refilter_deduplicated_inputs.md) for its separate
CSV schema, read-only input safeguards and restriction/distance settings.

For all-user Linux access, install a read-only version under `/opt/microc2tracks`
and use the [shared launcher and permission checks](docs/14_shared_multiuser_installation.md).
It reuses an existing environment without conda activation or maintainer-home
initialization. Configs, outputs, scratch and writable caches belong to each
user/project, not the shared installation. Verify actual shared resources from
another account before claiming server-wide access.

## Workflow Schematic

```mermaid
flowchart LR
  R["Validated reference registry<br/>including custom assemblies"] --> A
  A["Sample sheet<br/>paired FASTQ + reference_genome"] --> B["fastp<br/>QC + trimming"]
  B --> C["bwa-mem2<br/>Hi-C/Micro-C alignment"]
  C --> D["pairtools parse/sort<br/>retained original sorted pairs"]
  D --> Q["Independent per-row dedup"]
  Q --> E{"Assay profile"}
  E -->|"Micro-C"| F["MAPQ filter<br/>near-diagonal cis filter"]
  E -->|"Hi-C"| G["MAPQ + optional cis filter<br/>profile-specific restriction classification"]
  F --> H["indexed valid .pairs.gz"]
  G --> H
  H --> I["per-row matrices<br/>.cool + .mcool + .hic"]
  D --> U["User-selected dedup_group<br/>merge originals, dedup once"]
  U --> V["Same shared assay filtering"]
  V --> W["Accepted joint-unit pairs"]
  H --> T["Select disjoint units once<br/>joint unit replaces its rows"]
  W --> T
  T --> N["Technical / optional all-sample pool<br/>no extra global dedup"]
  N --> O["merged matrices<br/>.cool + .mcool + .hic"]
  I --> K["cooltools / Mustache / Chromosight<br/>TADs, compartments, loops, stripes"]
  O --> K
  I --> L["comparison plots<br/>Micro-C vs Hi-C"]
  O --> L
  I --> M["Juicebox / UCSC hic tracks"]
  O --> M
  I --> P["final report<br/>summary TSV/HTML + global MultiQC"]
  O --> P
```

## What It Produces

For each sample-sheet row, interpreted as one technical replicate:

- `fastp` QC, with trimmed FASTQ kept only if configured
- deduplicated `.pairs.gz` contacts
- MAPQ-filtered valid `.pairs.gz`
- raw and balanced `.cool` intermediates if configured
- final balanced `.mcool`
- final normalized `.hic`
- pairtools stats and MultiQC report

When multiple rows share the same `assay`, `condition`, `biological_replicate`, and canonical reference, the workflow also merges those technical replicates after per-row processing and rebuilds merged `.pairs.gz`, `.cool`, `.mcool`, and `.hic` files. Cross-reference merges are rejected. Samples without technical replicates should use `technical_replicate=1`; no redundant merged output is created for a one-row technical replicate group.

At the end of a successful run, `results/final_report/` contains compact all-sample summary tables, a small HTML report, a file manifest, and a global MultiQC report. Large reproducible intermediates such as trimmed FASTQ, temporary Juicer pairs, raw `.hic`, raw `.mcool`, and single-resolution `.cool` files are controlled by retention switches in `config.conf`.

Default light downstream analysis generates:

- expected contacts / contact-decay summaries
- insulation scores, browser-ready bedGraph tracks, boundary tables, and TAD-like BED intervals
- compartment eigenvectors, oriented when a phasing track is supplied

Optional heavier downstream scripts can generate:

- loop calls with `cooltools` and/or Mustache
- saddle plot inputs
- Micro-C vs Hi-C side-by-side matrix plots
- P(s) curve comparison plots

## Functionality Overview

| Area | Current functionality | Main outputs | Status |
|---|---|---|---|
| Input hygiene | Checks and normalizes config/sample-sheet text artifacts before parsing | Unix line endings, cleaned CSV fields | implemented |
| Per-sample references | Validated registry, including custom references and contig names | assembly-correct alignment, matrices, tracks, manifests | implemented |
| Resumability | Parameter/input/tool fingerprints, owned artifacts and verified publication | `logs/done/*.json`, safer restarts | implemented |
| Parallel samples | Optional sample-level parallelism with matrix and `.hic` semaphores | 1-4 concurrent sample workers, conservative defaults | implemented |
| FASTQ QC and trimming | Paired-end adapter detection and trimming | `fastp.html`, `fastp.json`, trimmed FASTQ | implemented |
| Alignment | Hi-C/Micro-C style chimeric-read alignment with `bwa-mem2 mem -SP5M` | streamed SAM into pairtools | implemented |
| Contact parsing | Parse, sort, deduplicate, and index contact pairs | `.dedup.pairs.gz`, `.valid.mapq*.pairs.gz`, `.px2` | implemented |
| Micro-C filtering | MAPQ filtering plus configurable short cis-distance filter | Micro-C-ready valid pairs | implemented |
| Canonical chromosomes | Optional canonical chromosome filter before `.cool`, `.mcool`, and `.hic` generation | UCSC-friendlier matrices and `.hic` files | enabled by default |
| Hi-C support | Shared workflow with explicit protocol, digest and artifact policy | annotated accepted pairs, QC, `.mcool`, `.hic` | implemented in candidate |
| Joint deduplication | User-selected dedup_group, original sorted streams, one joint pass | independent row outputs and joint-unit pairs/QC | implemented in candidate |
| Refilter old mapped pairs | Import already-deduplicated pairs, skip mapping/dedup, apply current filters | new accepted pairs, QC and rebuilt matrices | implemented in candidate |
| Biological pooling | Optional all-sample aggregate of disjoint units, without global dedup | explicitly named pooled matrices | implemented in candidate |
| Multi-disk storage | Optional pairs/matrix/trim/cache roots; SSD scratch; HDD final output | resolved paths and verified transfers | implemented in candidate |
| Matrix generation | Raw and balanced single-resolution and multiresolution matrices | `.cool`, `.mcool` | implemented |
| `.hic` export | Juicer-compatible `.hic` generation and normalization | `.raw.hic`, `.norm.hic` | implemented |
| Technical replicate merging | Merge filtered pair files within assay/condition/biological replicate groups and rebuild matrices | merged `.pairs.gz`, `.mcool`, `.hic` | implemented |
| Final run reporting | Summarize samples, merged groups, important files, and global QC | TSV, HTML, global MultiQC | implemented |
| Storage cleanup | Optional cleanup of large reproducible intermediates | smaller run directories | implemented |
| QC aggregation | Collect QC reports where available | MultiQC report | implemented |
| Preliminary downstream | Default light downstream pass from each technical-replicate and merged `.mcool` | expected contacts, insulation bedGraph, TAD BED, compartments | enabled; switchable |
| TADs / insulation | Run insulation score, bedGraph export, boundary calling, and TAD interval extraction from `.mcool` | insulation tables, bedGraph tracks, boundary calls, TAD BED | implemented |
| Compartments / saddle | Compartment eigenvectors and saddle-ready inputs | eigenvector tables, expected TSV | compartments enabled; phasing track optional for PC1 orientation |
| Loops | CPU-friendly loop calling | `cooltools dots`, Mustache outputs | implemented when tools are installed |
| Stripes | Stripe detection pathway | Chromosight outputs | planned/optional |
| Micro-C vs Hi-C comparison | Side-by-side matched-resolution matrix plots | PNG plots | implemented |

## Software Stack

| Step | Default software | Link | Notes |
|---|---|---|---|
| Environment management | Conda | [conda](https://docs.conda.io/) | `environment.yml` is conda-first; mamba is optional, not required. |
| FASTQ trimming/QC | fastp | [OpenGene/fastp](https://github.com/OpenGene/fastp) | Fast paired-end QC and adapter trimming. |
| Alignment | bwa-mem2 | [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) | Uses Hi-C/Micro-C-friendly `-SP5M` flags. |
| Reference indexing | samtools | [samtools](https://www.htslib.org/) | Creates `.fai` and chromosome sizes. |
| Contact parsing/filtering | pairtools | [open2c/pairtools](https://github.com/open2c/pairtools) | Central `.pairs.gz` workflow backbone. |
| Pair indexing | pairix | [4dn-dcic/pairix](https://github.com/4dn-dcic/pairix) | Indexes pairs for random access and cooler loading. |
| Matrix storage | cooler | [open2c/cooler](https://github.com/open2c/cooler) | Creates `.cool` and `.mcool` matrices. |
| `.hic` generation | Juicer Tools | [aidenlab/juicer](https://github.com/aidenlab/juicer) | Produces Juicebox/UCSC-compatible `.hic` files. |
| QC reporting | MultiQC | [MultiQC](https://multiqc.info/) | Aggregates available QC outputs. |
| Insulation/TADs | cooltools | [open2c/cooltools](https://github.com/open2c/cooltools) | Default downstream matrix analysis toolkit. |
| Loop calling | cooltools dots, Mustache | [cooltools](https://github.com/open2c/cooltools), [Mustache](https://github.com/ay-lab/mustache) | CPU-friendly defaults for no/limited GPU servers. |
| Stripes | Chromosight | [Chromosight](https://github.com/koszullab/chromosight) | Optional downstream module. |
| Pileups | coolpuppy | [coolpuppy](https://github.com/open2c/coolpuppy) | Optional aggregate loop/anchor analysis. |
| Format conversion | hictk | [hictk](https://github.com/paulsengroup/hictk) | Optional fast `.hic`/`.cool` toolkit. |
| Alternative downstream suite | HiCExplorer | [HiCExplorer](https://github.com/deeptools/HiCExplorer) | Optional TAD/visualization tools. |

## Recommended Design

Use one unified workflow with assay-specific configuration:

- Micro-C: no restriction fragment logic, default near-diagonal cis filter of 1000 bp.
- Hi-C: the shared path supports DpnII, MboI, NlaIII and HindIII digests, restriction-fragment classification and an independent optional distance cutoff. Legacy mode remains explicit.
- Both assays produce `.pairs.gz`, `.mcool`, and `.hic`, which makes downstream comparison consistent.
- By default, matrices and `.hic` files are filtered to canonical chromosomes with `FILTER_CANONICAL_CHROMS=true`; for UCSC tracks, the retained chromosome names must also match the UCSC assembly, e.g. `chr1`, `chr2`, `chrX`.

## Reference Selection

Each new sample-sheet row selects a registered reference. The supplied examples
include `mouse`/`mm39` and `human`/`hg38`; custom entries are supported. One TSV
registry supplies the corresponding FASTA, BWA-MEM2 index, chromosome sizes,
canonical regex, browser preset, and optional metadata. Unknown values fail.

Legacy sample sheets without the column default to `mm39` with a warning. Old
global mouse paths remain usable for that default row. Human selection is never
silently changed to mouse. See `docs/09_multi_reference.md` for complete mouse,
human, mixed-run, migration, chromosome-policy, and resumability examples.

## Quick Start

The following is a fresh/development installation example, not a shared-server
upgrade recipe. Existing shared installations should follow
[the versioned deployment guide](docs/14_shared_multiuser_installation.md) and
reuse working environments/indexes without changing active analyses.

```bash
git clone https://github.com/MichalGd/microc2tracks.git
cd microc2tracks

conda env create -f environment.yml
conda activate microc2tracks

cp config/config_template.conf config/config.conf
cp config/samplesheet_template.csv config/samplesheet.csv

# Edit config/config.conf and config/samplesheet.csv first.
# Required sample-sheet columns:
# sample,assay,reference_genome,condition,biological_replicate,technical_replicate,fastq_r1,fastq_r2
# Edit config/references.tsv for the exact references used by this run.
# Hi-C restriction rows require actual enzyme metadata; examples are not inference.
# Resource values are configurable examples; do not assume a particular server.
bash scripts/preflight_check.sh -c config/config.conf -s config/samplesheet.csv
bash scripts/microc2tracks.sh -c config/config.conf -s config/samplesheet.csv
```

Optional heavier downstream tools, beyond the default `cooltools` light pass,
can be installed in a separate environment:

```bash
conda env create -f envs/downstream_optional.yml
```

If classic conda spends a long time solving, enable conda's faster solver and retry:

```bash
conda config --set channel_priority strict
conda config --set solver libmamba
```

This still uses `conda`; it only changes the dependency solver.

## Shared Server Environment

Install versioned code outside user homes and use the existing shared environment
when compatible. The [shared installation guide](docs/14_shared_multiuser_installation.md)
provides the administrator-owned `/opt` layout, `/usr/local/bin` launcher,
cross-account access checks and project-owned outputs. Do not create duplicate
environments/indexes, recursively change active permissions or overwrite old
launchers while jobs are running. All-user server permissions and real scientific
tool behavior have not been validated merely by preparing this candidate.

For downstream analysis of one matrix:

```bash
bash scripts/run_downstream.sh \
  -c config/config.conf \
  -s microC1 \
  -m results/microC1/04_matrices/microC1.norm.mcool
```

Compare a Micro-C and Hi-C region:

```bash
python scripts/compare_matrices.py \
  --matrix-a results/microC1/04_matrices/microC1.norm.mcool \
  --matrix-b results/hic1/04_matrices/hic1.norm.mcool \
  --label-a Micro-C \
  --label-b Hi-C \
  --resolution 10000 \
  --region chr1:30000000-35000000 \
  --out results/comparison/microc_vs_hic_chr1.png
```

## Main Files

- `config/config_template.conf`: server, legacy-default-reference, tool, and resource defaults
- `config/references.tsv`: centralized mm39/hg38 registry with controlled aliases and browser/chromosome policy
- `config/samplesheet_template.csv`: per-row reference plus biological and technical replicate metadata
- `scripts/sanitize_text_inputs.py`: optional explicit in-place cleanup utility; the integrated runner reads configs without rewriting their originals
- `scripts/microc2tracks.sh`: FASTQ to pairs, `.cool`, `.mcool`, and `.hic`
- `scripts/merge_replicates.sh`: merge filtered pairs and rebuild matrices
- `scripts/run_plan.py`: resolved profiles, dedup units, pools and physical paths
- `scripts/filter_contact_pairs.py`: shared quality, distance and restriction filtering
- `scripts/artifact_state.py`: fingerprinted stages, safe publication and owned cleanup
- `scripts/run_downstream.sh`: common downstream analyses from `.mcool`
- `scripts/call_tads_from_insulation.py`: TAD-like BED intervals from cooltools insulation boundaries
- `scripts/export_insulation_bedgraph.py`: browser-ready insulation score bedGraph tracks
- `scripts/package_browser_tracks.py`: assembly-explicit UCSC/HiGlass track packages
- `scripts/reference_registry.py`: registry resolution and installed-reference preflight
- `scripts/check_fastq_pair_names.py`: optional full post-fastp R1/R2 synchronization safeguard
- `scripts/compare_matrices.py`: side-by-side matrix plotting
- `scripts/summarize_run.py`: final sample/merge summary reports
- `environment.yml`: core upstream matrix environment
- `envs/downstream_optional.yml`: optional heavier downstream tools
- `docs/`: analysis notes, installation, inputs, outputs, light downstream analysis, and troubleshooting

## Server Fit

The template provides conservative example resources; adjust threads, memory,
concurrency and storage reserves to the actual machine. GPU tools are not required.

## Status

The complete agreed feature set is implemented in one release candidate. Local
scientific/unit and orchestration-contract checks are distinct from the opt-in
real-tool tests. Do not deploy the candidate over an active installation or use
old reference-only completion markers to certify enhanced outputs.
