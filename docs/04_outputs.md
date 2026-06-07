# Outputs

For each sample-sheet row, outputs are written below. A row represents one technical replicate, so every technical replicate gets its own `.mcool` and `.hic` files before any merging happens:

```text
results/{sample}/
  01_qc/
    {sample}.fastp.html
    {sample}.fastp.json
    {sample}.multiqc.html
  02_trimmed/
    {sample}.trim.R1.fastq.gz
    {sample}.trim.R2.fastq.gz
  03_pairs/
    {sample}.dedup.pairs.gz
    {sample}.dedup.pairs.gz.px2
    {sample}.dedup.stats.txt
    {sample}.matrix.chrom.sizes
    {sample}.valid.mapq30.pairs.gz
    {sample}.valid.mapq30.pairs.gz.px2
    {sample}.pairs.stats.txt
    {sample}.valid.mapq30.juicer.pairs.gz
  04_matrices/
    {sample}.raw.1000.cool
    {sample}.norm.1000.cool
    {sample}.raw.mcool
    {sample}.norm.mcool
    {sample}.raw.hic
    {sample}.norm.hic
    {sample}.ucsc.hic.track.txt
  05_downstream/
    expected/
    insulation/
    compartments/
    loops/
    saddle/
  logs/
```

Merged technical replicate outputs use the same structure under:

```text
results/merged/{condition}_{assay}_B{biological_replicate}_tech_merged/
```

These merged outputs are created only when two or more rows share the same `assay`, `condition`, and `biological_replicate`. A biological replicate with a single technical replicate still has complete per-row outputs under `results/{sample}/`, but the workflow skips the redundant merge.

## Primary Outputs

Use `.norm.mcool` for most open2c downstream analyses and browser-style plots.

Use `.norm.hic` for Juicebox/Juicer ecosystem tools and UCSC `track type=hic` publication.

Keep `.valid.mapq30.pairs.gz` files because they are the most reusable intermediate for technical replicate merging and rebuilding matrices at new resolutions.

## UCSC Genome Browser Compatibility

The workflow creates `.hic` files with Juicer Tools:

```text
juicer_tools pre -> raw .hic
juicer_tools addNorm -> normalized .hic
```

This is the `.hic` route documented by UCSC for `track type=hic` custom tracks and track hubs. The file can be displayed in UCSC when these conditions are met:

- The `.hic` file was produced successfully by Juicer Tools.
- The matrix-specific `{sample}.matrix.chrom.sizes` file used by `juicer_tools pre` matches the chromosome names in the contact pairs.
- `FILTER_CANONICAL_CHROMS=true` is enabled, or `CHROM_SIZES` already contains only UCSC-compatible canonical chromosomes.
- The chromosome names match the UCSC genome database you want to view. For UCSC `mm39`, this usually means `chr1`, `chr2`, etc.; Ensembl-style names such as `1`, `2`, etc. may not display on the UCSC `mm39` browser unless you use a matching custom assembly hub.
- The `.hic` file is hosted at a public or intranet-accessible `http`, `https`, or `ftp` URL.
- UCSC can reach that URL through `bigDataUrl`.

If `PUBLIC_HIC_BASE_URL` is set in `config.conf`, the workflow writes:

```text
results/{sample}/04_matrices/{sample}.ucsc.hic.track.txt
```

Example:

```text
track type=hic name="microC1" description="microC1 normalized Hi-C/Micro-C contacts" bigDataUrl=https://server.example.org/hic/microC1.norm.hic
```

If `PUBLIC_HIC_DIR` is also set, the normalized `.hic` file is copied there so it can be served by a web server.
