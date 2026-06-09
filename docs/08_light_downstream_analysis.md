# Light Downstream Analysis

`RUN_PRELIM_DOWNSTREAM=true` enables a lightweight, default downstream pass from
each balanced `.norm.mcool` matrix. It runs for every technical replicate and
for every merged technical-replicate matrix. The goal is to produce the common
first-pass 3D genome summaries without delaying the run with heavier loop or
stripe callers.

Light mode runs:

- cis expected contacts / distance-decay summaries
- insulation scores and boundary calls
- TAD-like intervals derived from strong insulation boundaries
- compartment eigenvectors

Light mode skips:

- `cooltools dots`
- Mustache loop calling
- Chromosight stripe/loop calling

## Configuration

The relevant config block is:

```bash
RUN_PRELIM_DOWNSTREAM="true"
COOLTOOLS_THREADS="16"

TAD_RESOLUTION="10000"
INSULATION_WINDOWS_BP="50000 100000 250000"
TAD_BOUNDARY_THRESHOLD="Li"
TAD_BOUNDARY_WINDOW_BP="100000"
TAD_MIN_SIZE_BP="40000"
TAD_MIN_BOUNDARY_STRENGTH=""

COMPARTMENT_RESOLUTION="100000"
PHASING_TRACK=""
```

`TAD_RESOLUTION` and `COMPARTMENT_RESOLUTION` must exist in
`MCOOL_RESOLUTIONS`. If the selected resolution is absent, downstream analysis
will fail or skip that matrix.

## Outputs

For each technical replicate:

```text
results/{sample}/05_downstream/
  expected/
    {sample}.expected.{COMPARTMENT_RESOLUTION}.tsv
  insulation/
    {sample}.insulation.{TAD_RESOLUTION}.tsv
    bedgraph/
      {sample}.insulation_score.{window}.bedGraph
  tads/
    {sample}.tads.{TAD_BOUNDARY_WINDOW_BP}.bed
  compartments/
    {sample}.eigs.{COMPARTMENT_RESOLUTION}*
```

Merged technical replicates use the same structure under:

```text
results/merged/{merge_name}/05_downstream/
```

## Expected Contacts

Command:

```bash
cooltools expected-cis
```

Output:

```text
05_downstream/expected/{sample}.expected.{COMPARTMENT_RESOLUTION}.tsv
```

This table summarizes how contact frequency changes with genomic separation.
It is useful for global quality control and for comparing contact decay between
samples, assays, or conditions. In common 3D-genome language, this is closely
related to a P(s) contact-probability curve.

Further reading:

- Cooltools CLI reference: https://cooltools.readthedocs.io/en/latest/cli.html
- Cooltools paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC11098495/

## Insulation And Boundaries

Command:

```bash
cooltools insulation --threshold "${TAD_BOUNDARY_THRESHOLD}"
```

Output:

```text
05_downstream/insulation/{sample}.insulation.{TAD_RESOLUTION}.tsv
05_downstream/insulation/bedgraph/{sample}.insulation_score.{window}.bedGraph
```

Insulation measures local separation between upstream and downstream contacts.
Strong valleys in the insulation profile are commonly interpreted as candidate
domain boundaries. The output table contains window-specific columns such as:

- `log2_insulation_score_<window>`
- `boundary_strength_<window>`
- `is_boundary_<window>`

`TAD_BOUNDARY_THRESHOLD="Li"` asks cooltools to threshold boundary strengths
with the Li method. This can be changed if a project needs a different threshold
rule.

The bedGraph files contain `log2_insulation_score_<window>` as a genome-browser
track. By default they include a UCSC-style first line:

```text
track type=bedGraph name="sample_insulation_100000" visibility=full autoScale=on color=180,40,40
```

The data lines are four-column bedGraph:

```text
chrom  start  end  log2_insulation_score
```

Set `INSULATION_BEDGRAPH_TRACKLINE="false"` if another program requires a pure
four-column bedGraph without the track definition line.

Further reading:

- Cooltools insulation and boundaries notebook:
  https://cooltools.readthedocs.io/en/latest/notebooks/insulation_and_boundaries.html
- Cooltools CLI reference:
  https://cooltools.readthedocs.io/en/latest/cli.html
- UCSC bedGraph format:
  https://genome.ucsc.edu/goldenPath/help/bedgraph

## TAD-Like Domain Calls

Command:

```bash
python scripts/call_tads_from_insulation.py
```

Output:

```text
05_downstream/tads/{sample}.tads.{TAD_BOUNDARY_WINDOW_BP}.bed
```

The pipeline converts strong insulation boundaries into TAD-like intervals by
taking adjacent `is_boundary_<TAD_BOUNDARY_WINDOW_BP>` calls from the cooltools
insulation table. For example, with:

```bash
TAD_BOUNDARY_WINDOW_BP="100000"
```

the helper uses:

```text
is_boundary_100000
boundary_strength_100000
```

The output is BED-like:

```text
chrom  start  end  name  score  strand  left_boundary  right_boundary  left_strength  right_strength  window_bp
```

Interpretation:

- these are practical TAD/domain intervals for quick inspection;
- boundaries come directly from cooltools insulation;
- domains are intervals between adjacent strong boundaries;
- very small domains are removed with `TAD_MIN_SIZE_BP`;
- an optional extra cutoff can be set with `TAD_MIN_BOUNDARY_STRENGTH`.

Important caveat: TAD calls depend strongly on resolution, window size, library
depth, balancing, and boundary thresholding. Treat this as a first-pass,
consistent domain set, not as a final biological ground truth.

## Compartments

Command:

```bash
cooltools eigs-cis
```

Output prefix:

```text
05_downstream/compartments/{sample}.eigs.{COMPARTMENT_RESOLUTION}*
```

The pipeline calculates compartment eigenvectors by default. If
`PHASING_TRACK` is empty, eigenvectors are still produced, but the sign of PC1 is
arbitrary. That means positive PC1 does not automatically mean A compartment and
negative PC1 does not automatically mean B compartment.

For oriented A/B compartments, provide a phasing track:

```bash
PHASING_TRACK="/path/to/gc_or_gene_density_or_atac_or_rnaseq_track.tsv::track_column"
```

Cooltools expects a BedGraph-like phasing track with chromosome, start, end, and
a numeric track column. Common choices include GC content, gene density, ATAC-seq
signal, active histone-mark signal, or RNA-seq signal.

Further reading:

- Cooltools compartments and saddleplots notebook:
  https://cooltools.readthedocs.io/en/latest/notebooks/compartments_and_saddles.html
- Cooltools `eigs-cis` CLI documentation:
  https://cooltools.readthedocs.io/en/latest/cli.html#eigs-cis
- Cooltools paper:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC11098495/

## Runtime Expectations

For large Micro-C matrices, light downstream analysis is usually much cheaper
than alignment, pairtools sorting, matrix construction, and Juicer `.hic`
normalization. The most relevant cost drivers are:

- number of matrices: technical replicates plus merged matrices;
- `TAD_RESOLUTION`: 10 kb is more expensive than 25 kb or 50 kb;
- `INSULATION_WINDOWS_BP`: more windows means more boundary columns;
- `COMPARTMENT_RESOLUTION`: 100 kb is usually a reasonable default;
- `COOLTOOLS_THREADS`: more threads can help, but disk and matrix size still
  matter.

For urgent large-data runs, the recommended default is to keep light downstream
enabled but leave loop callers disabled. If the downstream step becomes a
bottleneck, set:

```bash
RUN_PRELIM_DOWNSTREAM="false"
```

and run `scripts/run_downstream.sh` later on selected `.norm.mcool` files.

## Practical Interpretation

The light analysis gives a compact first view of common 3D chromatin structure:

| Output | Biological readout | Best use |
|---|---|---|
| expected contacts | distance decay / P(s)-like behavior | global QC and sample comparison |
| insulation scores | local domain boundary signal | boundary strength and sample QC |
| TAD BED | TAD-like intervals between strong boundaries | genome browser tracks and quick domain comparisons |
| compartment eigenvectors | A/B compartment structure | broad active/inactive chromatin organization |

For final publication-grade calls, inspect representative loci, compare
replicates, and consider re-running downstream analysis with project-specific
resolutions, thresholds, and phasing tracks.
