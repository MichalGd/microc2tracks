# microc2tracks

`microc2tracks` is a shell-script-first Micro-C / Hi-C workflow for producing common contact matrix outputs and downstream analyses from paired-end FASTQ files.

The first implementation follows the style of the existing Shell repositories `fastq2tracks` and `rnaseq2tracksP`: a top-level runner, `config/`, `docs/`, `examples/`, `scripts/`, `environment.yml`, and simple tests.

## What It Produces

For each sample:

- trimmed FASTQ files and `fastp` QC
- deduplicated `.pairs.gz` contacts
- MAPQ-filtered valid `.pairs.gz`
- raw and balanced `.cool`
- raw and balanced `.mcool`
- raw and normalized `.hic`
- pairtools stats and MultiQC report

Optional downstream scripts can generate:

- insulation scores and TAD boundary tables
- compartment eigenvectors
- saddle plot inputs
- loop calls with `cooltools` and/or Mustache
- Micro-C vs Hi-C side-by-side matrix plots
- P(s) curve comparison plots

## Recommended Design

Use one unified workflow with assay-specific configuration:

- Micro-C: no restriction fragment logic, default near-diagonal cis filter of 1000 bp.
- Hi-C: shared alignment, pairs, matrix, and downstream path; restriction-enzyme-aware filtering can be added as a focused module when enzyme metadata are available.
- Both assays produce `.pairs.gz`, `.mcool`, and `.hic`, which makes downstream comparison consistent.

## Quick Start

```bash
conda env create -f environment.yml
conda activate microc2tracks

cp config/config_template.conf config/config.conf
cp config/samplesheet_template.csv config/samplesheet.csv

# Edit config/config.conf and config/samplesheet.csv first.
bash scripts/preflight_check.sh -c config/config.conf -s config/samplesheet.csv
bash scripts/microc2tracks.sh -c config/config.conf -s config/samplesheet.csv
```

Optional heavier downstream tools can be installed in a separate environment:

```bash
conda env create -f envs/downstream_optional.yml
```

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
- `config/samplesheet_template.csv`: sample metadata template
- `scripts/microc2tracks.sh`: FASTQ to pairs, `.cool`, `.mcool`, and `.hic`
- `scripts/merge_replicates.sh`: merge filtered pairs and rebuild matrices
- `scripts/run_downstream.sh`: common downstream analyses from `.mcool`
- `scripts/compare_matrices.py`: side-by-side matrix plotting
- `environment.yml`: core upstream matrix environment
- `envs/downstream_optional.yml`: optional heavier downstream tools
- `docs/`: analysis notes, installation, inputs, outputs, and troubleshooting

## Server Fit

The defaults target a shared CPU-heavy server with 500 GB RAM, 70 physical CPU cores, 140 logical CPUs, and no/limited GPU. GPU-dependent tools such as HiCCUPS are not required for the default path.

## Status

This is a first usable script-based workflow. The most important next extension is a stricter Hi-C restriction-fragment module when enzyme metadata and fragment BED files are available.
