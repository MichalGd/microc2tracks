# Tool Strategy

## Upstream Defaults

| Step | Default | Notes |
|---|---|---|
| QC/trimming | `fastp` | Already present in the existing environment. Fast and simple paired-end QC. |
| Alignment | `bwa-mem2 mem -SP5M` | Matches the legacy script and common Hi-C/Micro-C chimeric-read practice. |
| Pair parsing | `pairtools parse` | Produces the canonical `.pairs` intermediate. |
| Sorting/dedup | `pairtools sort`, `pairtools dedup` | CPU/RAM-heavy for large datasets. |
| Pair indexing | `pairix` | Required for `cooler cload pairix`. |
| Matrix generation | `cooler cload pairix`, `cooler zoomify` | Produces `.cool` and `.mcool`. |
| `.hic` generation | Juicer Tools `pre` and `addNorm` | Matches the legacy script and creates Juicebox-compatible output. |
| QC aggregation | `multiqc` | Collects fastp and other reports. |

## Downstream Defaults

| Output | Default | Notes |
|---|---|---|
| Insulation/TADs | `cooltools insulation` | Recommended first-line TAD/boundary method. |
| Compartments | `cooltools eigs-cis` | Requires a phasing track for confident A/B orientation. |
| Saddle plots | `cooltools saddle` | Uses compartment eigenvectors and expected contacts. |
| Loops | `cooltools dots`, Mustache | CPU-friendly alternatives to GPU-heavy HiCCUPS. |
| Loop anchors | derive from loop BEDPE/TSV | Anchor extraction can be added after choosing the loop caller. |
| Stripes | Chromosight | Optional downstream module. |
| Matrix comparison | `scripts/compare_matrices.py` | Plots matched-resolution regions from two `.mcool` files. |

## Missing From The Existing Conda Environment

The old `microc` environment already has the upstream tools. Add optional expanded downstream tools in a separate environment:

```bash
conda env create -f envs/downstream_optional.yml
```

This keeps the upstream matrix environment smaller and reduces the risk of dependency conflicts.

## Bottlenecks

- `bwa-mem2`: CPU-heavy, scales well to 32 threads.
- `pairtools sort`: CPU and disk I/O heavy; needs fast temporary storage.
- `pairtools merge`: heavy for multi-replicate Micro-C; allocate 64 threads and 128 GB RAM when available.
- `cooler cload` at 1 kb: memory and I/O heavy.
- `cooler balance`: can be slow at fine resolution.
- Juicer Tools `.hic` normalization: high Java heap recommended for large mammalian datasets.
- loop calling at 1 to 2 kb: the heaviest downstream task on CPU-only servers.
