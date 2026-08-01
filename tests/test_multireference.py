#!/usr/bin/env python3
"""Small synthetic regression tests for multi-reference behavior."""

from __future__ import annotations

import csv
import gzip
import re
import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

from check_fastq_pair_names import check_pairs  # noqa: E402
from reference_registry import (  # noqa: E402
    BWA_MEM2_INDEX_SUFFIXES,
    EXPECTED_CANONICAL,
    REGISTRY_COLUMNS,
    Reference,
    Registry,
    RegistryError,
    validate_reference_files,
)
from summarize_run import sample_summary_rows  # noqa: E402
from validate_pairs_reference import declared_assembly  # noqa: E402
from validate_samplesheet import NORMALIZED_COLUMNS, validate_sheet  # noqa: E402


class Fixture:
    def __init__(self, root: Path):
        self.root = root
        self.registry_path = root / "references.tsv"
        self.refs: dict[str, dict[str, str]] = {}
        self.add_reference("mm39", "mouse", "Mus musculus", EXPECTED_CANONICAL["mm39"])
        self.add_reference("hg38", "human", "Homo sapiens", EXPECTED_CANONICAL["hg38"])
        self.write_registry()

    def add_reference(self, reference_id: str, alias: str, species: str, chroms: list[str]) -> None:
        ref_dir = self.root / reference_id
        ref_dir.mkdir()
        fasta = ref_dir / f"{reference_id}.fa"
        fasta.write_text(">synthetic\nA\n", encoding="utf-8")
        sizes = {chrom: 1000 + index for index, chrom in enumerate(chroms)}
        (ref_dir / f"{reference_id}.fa.fai").write_text(
            "".join(f"{name}\t{length}\t0\t0\t0\n" for name, length in sizes.items()),
            encoding="utf-8",
        )
        chrom_sizes = ref_dir / f"{reference_id}.chrom.sizes"
        chrom_sizes.write_text(
            "".join(f"{name}\t{length}\n" for name, length in sizes.items()),
            encoding="utf-8",
        )
        for suffix in BWA_MEM2_INDEX_SUFFIXES:
            Path(f"{fasta}{suffix}").write_text("synthetic\n", encoding="utf-8")
        regex = (
            r"^chr([1-9]|1[0-9]|X|Y|M)$"
            if reference_id == "mm39"
            else r"^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$"
        )
        self.refs[reference_id] = {
            "reference_id": reference_id,
            "aliases": alias,
            "species": species,
            "assembly": reference_id,
            "fasta": str(fasta),
            "bwa_index_prefix": str(fasta),
            "chrom_sizes": str(chrom_sizes),
            "canonical_regex": regex,
            "browser_preset": alias,
            "phasing_track": "",
            "phasing_assembly": "",
            "annotation_metadata": "not consumed",
            "blacklist_metadata": "not consumed",
        }

    def write_registry(self) -> None:
        with self.registry_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=REGISTRY_COLUMNS, delimiter="\t")
            writer.writeheader()
            writer.writerows(self.refs.values())

    def write_sheet(self, name: str, rows: list[dict[str, str]], include_reference: bool = True) -> Path:
        path = self.root / name
        fields = NORMALIZED_COLUMNS if include_reference else [
            "sample",
            "assay",
            "condition",
            "biological_replicate",
            "technical_replicate",
            "fastq_r1",
            "fastq_r2",
        ]
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
        return path


def row(sample: str, reference: str, condition: str = "control", bio: str = "1", tech: str = "1") -> dict[str, str]:
    return {
        "sample": sample,
        "assay": "microc",
        "reference_genome": reference,
        "condition": condition,
        "biological_replicate": bio,
        "technical_replicate": tech,
        "fastq_r1": "/synthetic/R1.fastq.gz",
        "fastq_r2": "/synthetic/R2.fastq.gz",
        "merge_group": "",
    }


class MultiReferenceTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.fixture = Fixture(self.root)
        self.registry = Registry(self.fixture.registry_path)

    def tearDown(self) -> None:
        self.temp.cleanup()

    def test_alias_normalization_and_browser_presets(self) -> None:
        self.assertEqual(self.registry.resolve("mouse").reference_id, "mm39")
        self.assertEqual(self.registry.resolve("mm39").reference_id, "mm39")
        self.assertEqual(self.registry.resolve("human").reference_id, "hg38")
        self.assertEqual(self.registry.resolve("hg38").reference_id, "hg38")
        self.assertEqual(self.registry.resolve("human").browser_preset, "human")

    def test_unknown_reference_is_rejected(self) -> None:
        with self.assertRaisesRegex(RegistryError, "unknown reference"):
            self.registry.resolve("GRCh37")

    def test_valid_mouse_human_and_mixed_sheets(self) -> None:
        sheet = self.fixture.write_sheet(
            "mixed.csv",
            [row("mouse_sample", "mouse", "mouse_control"), row("human_sample", "human", "human_control")],
        )
        rows, _, _ = validate_sheet(sheet, self.registry, "mm39")
        self.assertEqual([item["reference_genome"] for item in rows], ["mm39", "hg38"])

    def test_missing_reference_column_defaults_to_mouse_with_warning(self) -> None:
        sheet = self.fixture.write_sheet("legacy.csv", [row("legacy_mouse", "")], include_reference=False)
        rows, _, warnings = validate_sheet(sheet, self.registry, "mm39")
        self.assertEqual(rows[0]["reference_genome"], "mm39")
        self.assertTrue(any("absent" in warning and "mm39" in warning for warning in warnings))

    def test_explicit_unknown_reference_in_sheet_is_rejected(self) -> None:
        sheet = self.fixture.write_sheet("unknown.csv", [row("bad", "unknown")])
        with self.assertRaisesRegex(RegistryError, "unknown reference"):
            validate_sheet(sheet, self.registry, "mm39")

    def test_cross_reference_replicate_merge_is_rejected(self) -> None:
        sheet = self.fixture.write_sheet(
            "cross.csv",
            [row("mouse_t1", "mm39", tech="1"), row("human_t2", "hg38", tech="2")],
        )
        with self.assertRaisesRegex(RegistryError, "prohibited cross-reference replicate merge"):
            validate_sheet(sheet, self.registry, "mm39")

    def test_explicit_merge_group_requires_compatible_metadata(self) -> None:
        first = row("mouse_t1", "mm39", condition="control", tech="1")
        second = row("mouse_t2", "mm39", condition="treated", tech="2")
        first["merge_group"] = second["merge_group"] = "forced_group"
        sheet = self.fixture.write_sheet("incompatible.csv", [first, second])
        with self.assertRaisesRegex(RegistryError, "incompatible rows in merge_group"):
            validate_sheet(sheet, self.registry, "mm39")

    def test_expected_canonical_mouse_and_human_chromosomes(self) -> None:
        for reference_id, expected in EXPECTED_CANONICAL.items():
            ref = self.registry.resolve(reference_id)
            retained = [name for name in expected if re.fullmatch(ref.canonical_regex, name)]
            self.assertEqual(retained, expected)
            self.assertNotRegex("chrUn_random", re.compile(ref.canonical_regex))
        self.assertIn("chrY", EXPECTED_CANONICAL["mm39"])
        self.assertIn("chrY", EXPECTED_CANONICAL["hg38"])

    def test_registry_file_preflight_validates_both_references(self) -> None:
        for reference_id in ("mm39", "hg38"):
            self.assertEqual(validate_reference_files(self.registry.resolve(reference_id)), [])

    def test_fasta_chrom_size_mismatch_is_rejected(self) -> None:
        ref = self.registry.resolve("mm39")
        path = Path(ref.chrom_sizes)
        text = path.read_text(encoding="utf-8").replace("chr1\t1000", "chr1\t999")
        path.write_text(text, encoding="utf-8")
        errors = validate_reference_files(ref)
        self.assertTrue(any("mismatch" in error for error in errors))

    def test_missing_bwa_index_is_rejected(self) -> None:
        ref = self.registry.resolve("hg38")
        Path(f"{ref.bwa_index_prefix}.0123").unlink()
        errors = validate_reference_files(ref)
        self.assertTrue(any("BWA-MEM2 index component" in error for error in errors))

    def test_mitochondrial_mt_and_chrmt_are_not_silently_renamed(self) -> None:
        ref = self.registry.resolve("mm39")
        for mitochondrial in ("MT", "chrMT"):
            source = Path(ref.chrom_sizes).read_text(encoding="utf-8").replace("chrM\t1021", f"{mitochondrial}\t1021")
            chrom_path = self.root / f"{mitochondrial}.chrom.sizes"
            fai_path = self.root / f"{mitochondrial}.fa.fai"
            fasta_path = self.root / f"{mitochondrial}.fa"
            chrom_path.write_text(source, encoding="utf-8")
            fai_path.write_text(source, encoding="utf-8")
            fasta_path.write_text(">synthetic\nA\n", encoding="utf-8")
            custom = replace(ref, fasta=str(fasta_path), chrom_sizes=str(chrom_path))
            errors = validate_reference_files(custom)
            self.assertTrue(any("mitochondrial names" in error for error in errors))

    def test_pairs_header_assembly_is_detected(self) -> None:
        pairs = self.root / "pairs.gz"
        with gzip.open(pairs, "wt", encoding="utf-8") as handle:
            handle.write("## pairs format v1.0\n#genome_assembly: hg38\nread\tchr1\t1\tchr1\t2\n")
        self.assertEqual(declared_assembly(pairs), "hg38")

    def test_fastq_name_mismatch_is_detected(self) -> None:
        r1 = self.root / "R1.fastq"
        r2 = self.root / "R2.fastq"
        r1.write_text("@read1/1\nAC\n+\nII\n", encoding="utf-8")
        r2.write_text("@read2/2\nTG\n+\nII\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "read-name mismatch"):
            check_pairs(r1, r2)

    def test_reference_is_carried_into_final_report_rows(self) -> None:
        manifest = [{
            "sample": "human_sample",
            "assay": "microc",
            "reference_id": "hg38",
            "species": "Homo sapiens",
            "assembly": "hg38",
            "browser_preset": "human",
            "sample_dir": str(self.root / "human_sample"),
        }]
        rows = sample_summary_rows(manifest, self.root, "30", "1000")
        self.assertEqual(rows[0]["reference_id"], "hg38")
        self.assertEqual(rows[0]["browser_preset"], "human")


if __name__ == "__main__":
    unittest.main(verbosity=2)
