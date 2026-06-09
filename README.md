# microc2tracks

`microc2tracks` is a shell-script-first Micro-C / Hi-C workflow for producing common contact matrix outputs and downstream analyses from paired-end FASTQ files.

The first implementation follows the style of the existing Shell repositories `fastq2tracks` and `rnaseq2tracksP`: a top-level runner, `config/`, `docs/`, `examples/`, `scripts/`, `environment.yml`, and simple tests.

## Workflow Schematic

```mermaid
flowchart LR
  A["Sample sheet<br/>paired FASTQ"] --> B["fastp<br/>QC + trimming"]
  B --> C["bwa-mem2<br/>Hi-C/Micro-C alignment"]
  C --> D["pairtools parse/sort/dedup<br/>deduplicated pairs"]
  D --> E{"Assay profile"}
  E -->|"Micro-C"| F["MAPQ filter<br/>near-diagonal cis filter"]
  E -->|"Hi-C"| G["MAPQ filter<br/>restriction-aware module planned"]
  F --> H["indexed valid .pairs.gz"]
  G --> H
  H --> I["per-row matrices<br/>.cool + .mcool + .hic"]
  H --> N["technical replicate merge<br/>same assay + condition + biological replicate"]
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

When multiple rows share the same `assay`, `condition`, and `biological_replicate`, the workflow also merges those technical replicates after per-row processing and rebuilds merged `.pairs.gz`, `.cool`, `.mcool`, and `.hic` files. Samples without technical replicates should use `technical_replicate=1`; no redundant merged output is created for a one-row technical replicate group.

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
| Resumability | Per-step sentinels plus atomic temp outputs for major files | `logs/done/*.done`, safer restarts | implemented |
| Parallel samples | Optional sample-level parallelism with matrix and `.hic` semaphores | 1-4 concurrent sample workers, conservative defaults | implemented |
| FASTQ QC and trimming | Paired-end adapter detection and trimming | `fastp.html`, `fastp.json`, trimmed FASTQ | implemented |
| Alignment | Hi-C/Micro-C style chimeric-read alignment with `bwa-mem2 mem -SP5M` | streamed SAM into pairtools | implemented |
| Contact parsing | Parse, sort, deduplicate, and index contact pairs | `.dedup.pairs.gz`, `.valid.mapq*.pairs.gz`, `.px2` | implemented |
| Micro-C filtering | MAPQ filtering plus configurable short cis-distance filter | Micro-C-ready valid pairs | implemented |
| Canonical chromosomes | Optional canonical chromosome filter before `.cool`, `.mcool`, and `.hic` generation | UCSC-friendlier matrices and `.hic` files | enabled by default |
| Hi-C support | Shared FASTQ-to-pairs-to-matrices path | Hi-C `.pairs.gz`, `.mcool`, `.hic` | implemented; restriction-fragment filtering planned |
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
- Hi-C: shared alignment, pairs, matrix, and downstream path; restriction-enzyme-aware filtering can be added as a focused module when enzyme metadata are available.
- Both assays produce `.pairs.gz`, `.mcool`, and `.hic`, which makes downstream comparison consistent.
- By default, matrices and `.hic` files are filtered to canonical chromosomes with `FILTER_CANONICAL_CHROMS=true`; for UCSC tracks, the retained chromosome names must also match the UCSC assembly, e.g. `chr1`, `chr2`, `chrX`.

## Quick Start

```bash
git clone https://github.com/MichalGd/microc2tracks.git
cd microc2tracks

conda env create -f environment.yml
conda activate microc2tracks

cp config/config_template.conf config/config.conf
cp config/samplesheet_template.csv config/samplesheet.csv

# Edit config/config.conf and config/samplesheet.csv first.
# Required sample-sheet columns:
# sample,assay,condition,biological_replicate,technical_replicate,fastq_r1,fastq_r2
# For a quiet large-memory server, try MAX_PARALLEL_SAMPLES=2 after one serial test.
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

For an installation available to all users, create the conda environment at a shared prefix instead of using a user-local named environment:

```bash
sudo mkdir -p /opt
sudo git clone https://github.com/MichalGd/microc2tracks.git /opt/microc2tracks
cd /opt/microc2tracks

# Load conda from wherever it is installed on this server.
source "$(conda info --base)/etc/profile.d/conda.sh"

conda env create \
  -p /opt/conda/envs/microc2tracks \
  -f environment.yml

conda activate /opt/conda/envs/microc2tracks
```

Then make the pipeline, environment, references, and Juicer Tools readable/executable:

```bash
sudo chgrp -R bioinfo /opt/microc2tracks /opt/conda/envs/microc2tracks
sudo chmod -R a+rX /opt/microc2tracks /opt/conda/envs/microc2tracks
sudo chmod -R a+rX /shared/references /shared/software/juicer_tools
```

See `docs/02_installation_server.md` for optional `/usr/local/bin` launchers.

The server guide also includes a tested `biolserv` installation recipe covering `/opt` cloning, `/opt/conda/envs/microc2tracks`, `libmamba` solver setup, permissions without a `bioinfo` group, and shared launchers.

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

- `config/config_template.conf`: server, reference, tool, and resource defaults
- `config/samplesheet_template.csv`: sample metadata template with biological and technical replicate columns
- `scripts/sanitize_text_inputs.py`: automatic line-ending/BOM/CSV-field cleanup used by preflight and the main runner
- `scripts/microc2tracks.sh`: FASTQ to pairs, `.cool`, `.mcool`, and `.hic`
- `scripts/merge_replicates.sh`: merge filtered pairs and rebuild matrices
- `scripts/run_downstream.sh`: common downstream analyses from `.mcool`
- `scripts/call_tads_from_insulation.py`: TAD-like BED intervals from cooltools insulation boundaries
- `scripts/export_insulation_bedgraph.py`: browser-ready insulation score bedGraph tracks
- `scripts/compare_matrices.py`: side-by-side matrix plotting
- `scripts/summarize_run.py`: final sample/merge summary reports
- `environment.yml`: core upstream matrix environment
- `envs/downstream_optional.yml`: optional heavier downstream tools
- `docs/`: analysis notes, installation, inputs, outputs, light downstream analysis, and troubleshooting

## Server Fit

The defaults target a shared CPU-heavy server with 500 GB RAM, 70 physical CPU cores, 140 logical CPUs, and no/limited GPU. GPU-dependent tools such as HiCCUPS are not required for the default path.

## Status

This is a first usable script-based workflow. The most important next extension is a stricter Hi-C restriction-fragment module when enzyme metadata and fragment BED files are available.
