# Outputs

For each sample, outputs are written below:

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
  05_downstream/
    expected/
    insulation/
    compartments/
    loops/
    saddle/
  logs/
```

Merged replicate outputs use the same structure under:

```text
results/merged/{merge_name}/
```

## Primary Outputs

Use `.norm.mcool` for most open2c downstream analyses and browser-style plots.

Use `.norm.hic` for Juicebox/Juicer ecosystem tools and UCSC `track type=hic` publication.

Keep `.valid.mapq30.pairs.gz` files because they are the most reusable intermediate for merging and rebuilding matrices at new resolutions.

