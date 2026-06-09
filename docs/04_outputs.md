# Outputs

For each sample-sheet row, outputs are written below. A row represents one technical replicate, so every technical replicate gets its own `.mcool` and `.hic` files before any merging happens:

```text
results/{sample}/
  01_qc/
    {sample}.fastp.html
    {sample}.fastp.json
    {sample}.multiqc.html
  02_trimmed/
    {sample}.trim.R1.fastq.gz              # optional; removed by default after success
    {sample}.trim.R2.fastq.gz              # optional; removed by default after success
  03_pairs/
    {sample}.dedup.pairs.gz                # kept by default during method development
    {sample}.dedup.pairs.gz.px2
    {sample}.dedup.stats.txt
    {sample}.matrix.chrom.sizes
    {sample}.valid.mapq30.pairs.gz
    {sample}.valid.mapq30.pairs.gz.px2
    {sample}.pairs.stats.txt
    {sample}.valid.mapq30.juicer.pairs.gz  # optional; removed by default after .hic creation
  04_matrices/
    {sample}.raw.1000.cool                 # optional; removed by default after .mcool creation
    {sample}.norm.1000.cool                # optional; removed by default after .mcool creation
    {sample}.raw.mcool                     # optional; removed by default
    {sample}.norm.mcool
    {sample}.raw.hic                       # optional; removed by default after addNorm
    {sample}.norm.hic
    {sample}.ucsc.hic.track.txt
  05_downstream/
    expected/
    insulation/
      bedgraph/
    tads/
    compartments/
    loops/
    saddle/
  logs/
    {sample}.status.tsv
    done/
      fastp.done
      dedup_pairs.done
      valid_pairs.done
      norm_mcool.done
      norm_hic.done
```

Merged technical replicate outputs use the same structure under:

```text
results/merged/{condition}_{assay}_B{biological_replicate}_tech_merged/
```

These merged outputs are created only when two or more rows share the same `assay`, `condition`, and `biological_replicate`. A biological replicate with a single technical replicate still has complete per-row outputs under `results/{sample}/`, but the workflow skips the redundant merge.

Run-level metadata and final reports are written under:

```text
results/run_metadata/
  sample_manifest.tsv
  technical_replicates.tsv

results/final_report/
  microc2tracks_sample_summary.tsv
  microc2tracks_merge_summary.tsv
  microc2tracks_file_manifest.tsv
  microc2tracks_summary.html
  microc2tracks_multiqc.html
```

## Primary Outputs

Use `.norm.mcool` for most open2c downstream analyses and browser-style plots.

Use `.norm.hic` for Juicebox/Juicer ecosystem tools and UCSC `track type=hic` publication.

Keep `.valid.mapq30.pairs.gz` files because they are the most reusable intermediate for technical replicate merging and rebuilding matrices at new resolutions.

Trimmed FASTQ, temporary Juicer pair files, raw `.hic`, raw `.mcool`, and single-resolution `.cool` files are reproducible from retained upstream files. Their retention is controlled in `config.conf`.

`logs/done/*.done` files are resumability sentinels. They record successful
step completion after outputs exist. Current runners can bootstrap missing
sentinels from valid existing outputs, which lets older completed runs resume
without recomputation.

By default, `RUN_PRELIM_DOWNSTREAM=true`, so light downstream outputs are written
under `05_downstream/` for each technical replicate and for each merged
technical-replicate matrix. This mode is intended as a quick first pass and
skips heavier loop callers. It includes expected contacts, insulation/TAD
boundary-style tables, TAD-like BED intervals, and compartment eigenvectors. If
`PHASING_TRACK` is not set, compartment PC1 signs are arbitrary and should be
oriented later. Set `RUN_PRELIM_DOWNSTREAM=false` to omit this step. See
`docs/08_light_downstream_analysis.md` for details.

The default TAD BED file is:

```text
05_downstream/tads/{sample}.tads.{TAD_BOUNDARY_WINDOW_BP}.bed
```

Insulation score bedGraph tracks are written for each `INSULATION_WINDOWS_BP`
value:

```text
05_downstream/insulation/bedgraph/{sample}.insulation_score.{window}.bedGraph
```

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
