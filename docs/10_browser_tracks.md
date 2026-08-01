# Browser Track Packaging

The workflow writes `.hic` track text beside each normalized `.hic`. Its comment
records `reference_id`, assembly, and browser preset. When publishing the file,
use the matching UCSC database or custom assembly; do not use an mm39 URL/track
for hg38 data.

The production-derived `scripts/package_browser_tracks.py` packages insulation
and compartment bedGraphs for UCSC and HiGlass. Reference identity is required:

```bash
python scripts/package_browser_tracks.py \
  --reference-id hg38 \
  --sample human_control=/results/merged/human_control_microc_B1_tech_merged \
  --out-dir /results/browser/human_control
```

`--reference-id mm39` selects the mouse canonical set; `--reference-id hg38`
selects the human set. `--canonical-preset` is optional, but if supplied it must
agree with the reference. `--canonical-chroms` is an explicit advanced override
and does not rename source chromosomes.

Do not combine mouse and human samples in one browser package. Create separate
packages so chromosome-size and ingestion metadata remain assembly-specific.
The generated `track_manifest.tsv` records the canonical reference and preset.

bedGraph exports keep the source matrix names. Thus `chr1` is accepted for the
shipped UCSC presets, while `1`, `MT`, `chrMT`, alt/random/unplaced contigs, and
NCBI accessions are neither guessed nor silently rewritten.
