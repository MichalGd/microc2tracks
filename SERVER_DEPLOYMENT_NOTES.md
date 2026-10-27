# Server Deployment Notes: 3.0.0-rc1

This is a complete feature candidate, prepared for staged validation. No server
installation or running analysis was changed. The real-tool release gate in
[TEST_REPORT.md](TEST_REPORT.md) remains outstanding.

## Safe staging

1. Keep the active installation and its environment untouched.
2. Extract the complete candidate archive to a new versioned staging directory.
   Deploy the entire set in UPLOAD_MANIFEST.txt, not selected shell scripts.
3. Create a separate environment, for example:

   ```bash
   conda env create -n microc2tracks-3.0.0-rc1 -f environment.yml
   conda activate microc2tracks-3.0.0-rc1
   ```

4. Run the small unit tests and opt-in synthetic real-tool tests from TEST_REPORT.
   No production FASTQs are needed. An environment solve failure or smoke-test
   failure must be resolved before promotion.
5. Configure the references actually used, actual library chemistry/enzymes,
   intentional dedup groups, resource limits and output/storage roots. Examples
   are not server-specific deployment instructions.
6. Run preflight against a reviewed CSV/config and a new OUTDIR. No automatic
   adoption or deletion of old sample outputs is performed.
7. Switch a launcher/symlink only after the real-tool gate and review succeed.
   Retain the previous installation and environment for rollback.

## References

Derive FASTA index, chromosome sizes and BWA-MEM2 components from the same FASTA.
Use the validated registry for arbitrary assemblies and contig names. The named
mm39/hg38 browser presets keep their explicit naming policies; those policies are
not applied to custom reference entries. Generated restriction maps live in a
run-owned or explicitly configured cache, never in the source reference directory.

## Storage and restart

Prefer durable outputs on HDD and active scratch on SSD, with reserves and
REQUIRED_STORAGE_MOUNTS configured for the actual mount points. Alternative
category roots are optional; empty settings keep the standard OUTDIR structure.
Publication/cleanup uses stage ownership, never recursive deletion of roots.

Legacy reference-only .done files do not certify the enhanced workflow. Use a
new OUTDIR for old-run migration/reference changes and preserve originals.
Already accepted independent outputs can be pooled only using the documented
input modes; joint deduplication needs original sorted pre-dedup streams.

## Rollback

Stop only new submissions to the candidate, restore the earlier launcher/config
and keep candidate output trees separate. Do not move candidate outputs into an
old run or remove original unmerged outputs. Preserve failed stage logs and
retained scratch for diagnosis.

See [the integration/storage guide](docs/12_integrated_hic_and_storage.md).
