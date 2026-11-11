"""Native two-motif Arima selection and fragment semantics; tiny local fixtures."""
import copy
import itertools
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import test_integration_features as fixtures
from pairs_processing import PairProcessor
from pipeline import Runtime
from prepare_restriction_fragments import load_digest, prepare, scan_fasta
from run_plan import build_plan
from workflow_common import registry_data, resolve_policy, resolve_preprocessing


ARIMA_DEFINITIONS = {
    "DpnII": {"motif": "GATC", "cut_offset": 0},
    "HinfI": {"motif": "GANTC", "cut_offset": 1},
}


class ArimaTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.IntegrationFeatures(
            "test_profiles_generic_dovetail_legacy_and_explicit_enzyme_requirement")
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.addCleanup(self.fixture.tearDown)
        self.root = self.fixture.root
        self.c = {**self.fixture.c, "HIC_PROTOCOL": "arima_2enzyme",
                  "HIC_RESTRICTION_ENZYMES": "", "RUN_MERGE_TECHNICAL_REPLICATES": "false",
                  "RUN_MERGE_ALL_SAMPLES": "false", "RUN_HIC": "false"}
        self.fixture.rows = [
            {**row, "assay": "hic", "dedup_group": "", "condition": condition,
             "hic_protocol": "", "restriction_enzymes": ""}
            for row, condition in zip(self.fixture.rows[:2], ("WT", "DKO"))]
        self.fixture.write_sheet()

    def registry_file(self, data):
        path = self.root / "protocols.json"
        path.write_text(json.dumps(data))
        return path

    def test_profile_defaults_explicit_precedence_and_conflicts(self):
        policy = resolve_policy({"assay": "hic"}, self.c)
        self.assertEqual(policy["profile"], "arima_2enzyme")
        self.assertEqual(policy["enzyme_definitions"], ARIMA_DEFINITIONS)
        self.assertEqual(policy["enzymes"], ["DpnII", "HinfI"])
        self.assertEqual(resolve_preprocessing({}, self.c, policy), {})
        for incorrect in ("DpnII", "MboI+HinfI", "HaeIII", "DpnII+HinfI+NlaIII"):
            with self.subTest(enzymes=incorrect), self.assertRaises(ValueError):
                resolve_policy({"assay": "hic"}, {**self.c, "HIC_RESTRICTION_ENZYMES": incorrect})
        explicit = resolve_policy({"assay": "hic", "restriction_enzymes": "hinfi + dpnii"},
                                  {**self.c, "HIC_RESTRICTION_ENZYMES": "HaeIII"})
        self.assertEqual(explicit, policy)
        generic = resolve_policy({"assay": "hic", "hic_protocol": "restriction_hic",
                                  "restriction_enzymes": "DpnII+HinfI"}, self.c)
        self.assertEqual(generic["enzyme_definitions"], ARIMA_DEFINITIONS)

    def test_named_profile_checks_definitions_not_only_enzyme_names(self):
        original = registry_data(ROOT / "config/hic_protocols.json")
        for enzyme, definition in (("HinfI", {"motif": "GANTC", "cut_offset": 0}),
                                   ("HinfI", {"motif": "GAATTC", "cut_offset": 1}),
                                   ("DpnII", {"motif": "GATC", "cut_offset": 1})):
            data = copy.deepcopy(original)
            data["enzymes"][enzyme] = definition
            config = {**self.c, "HIC_PROTOCOL_REGISTRY": str(self.registry_file(data))}
            with self.subTest(enzyme=enzyme, definition=definition), self.assertRaises(ValueError):
                resolve_policy({"assay": "hic"}, config)

    def test_registry_accepts_degenerate_palindrome_but_not_regex_or_nonpalindromes(self):
        original = registry_data(ROOT / "config/hic_protocols.json")
        self.assertEqual(original["enzymes"]["HinfI"], ARIMA_DEFINITIONS["HinfI"])
        for motif in ("G.NTC", "GA[ACGT]TC", "^GANTC", "G^ANTC", "GACTC", "GARTC", ""):
            data = copy.deepcopy(original)
            data["enzymes"]["HinfI"]["motif"] = motif
            with self.subTest(motif=motif), self.assertRaises(ValueError):
                registry_data(self.registry_file(data))

    def test_all_four_n_variants_chunk_edges_and_ambiguous_reference(self):
        # Motif starts: 0, 6, 13, 20, 27, 34. The final GANTC in the
        # reference is unknown DNA and must not be accepted as a fifth site.
        sequence = "GATCNNGAATCNNGACTCNNGAGTCNNGATTCNNGANTC".lower()
        fasta = self.root / "variants.fa"
        # One base per line exercises real chunk boundaries, not just a small
        # chunk_size argument applied to a single long input line.
        fasta.write_text(">variants\n" + "\n".join(sequence) + "\n>second\nG\nA\nA\nT\nC\n")
        for chunk in (1, 3, 4, 5, 8):
            with self.subTest(chunk=chunk):
                observed = list(scan_fasta(fasta, ARIMA_DEFINITIONS, chunk_size=chunk))
                self.assertEqual(observed, [("variants", 39, [0, 7, 14, 21, 28, 39]),
                                            ("second", 5, [0, 1, 5])])
        duplicated = {**ARIMA_DEFINITIONS, "MboI": ARIMA_DEFINITIONS["DpnII"]}
        self.assertEqual(list(scan_fasta(fasta, duplicated, chunk_size=3)), observed)

    def test_existing_literal_enzymes_and_overlapping_hits_are_unchanged(self):
        fasta = self.root / "literal.fa"
        fasta.write_text(">literal\n" + "\n".join("ttGATCaaCATGaaAAGCTTnnGGCC") + "\n")
        data = registry_data(ROOT / "config/hic_protocols.json")
        for enzyme, position in (("DpnII", 2), ("MboI", 2), ("NlaIII", 12),
                                 ("HindIII", 15), ("HaeIII", 24)):
            with self.subTest(enzyme=enzyme):
                self.assertEqual(list(scan_fasta(fasta, {enzyme: data["enzymes"][enzyme]}, 3)),
                                 [("literal", 26, [0, position, 26])])
        fasta.write_text(">overlap\nA\nT\nA\nT\nA\nT\n")
        self.assertEqual(list(scan_fasta(fasta, {"toy": {"motif": "ATAT", "cut_offset": 2}}, 3)),
                         [("overlap", 6, [0, 2, 4, 6])])

    def test_combined_digest_cache_and_native_fragment_filter(self):
        length = len(self.fixture.sequence)
        self.fixture.fasta.write_text(
            ">scaffoldA\n" + "A"*8 + "GATC" + "AA" + "GAATC" + "A"*(length-19)
            + "\n>plasmidZ\n" + "A"*length + "\n")
        policy = resolve_policy({"assay": "hic"}, self.c)
        bed, metadata = prepare(self.fixture.fasta, self.fixture.sizes,
                                policy["enzyme_definitions"], self.root / "cache")
        digest = load_digest(bed)
        self.assertEqual(list(digest["scaffoldA"]), [0, 8, 15, length])
        self.assertEqual(metadata["identity"]["generator_version"], 2)
        modified = bed.stat().st_mtime_ns
        self.assertEqual(prepare(self.fixture.fasta, self.fixture.sizes,
                                 policy["enzyme_definitions"], self.root / "cache"), (bed, metadata))
        self.assertEqual(bed.stat().st_mtime_ns, modified)
        records = [fixtures.contact("same_fragment", 10, 14),
                   fixtures.contact("across_DpnII", 8, 9),
                   fixtures.contact("across_HinfI", 15, 16),
                   fixtures.contact("self_circle", 16, 19, s1="-", s2="+"),
                   fixtures.contact("trans", 10, 14, c2="plasmidZ")]
        accepted, _ = self.fixture.filter(records, policy, digest)
        self.assertEqual([row[0] for row in accepted], ["across_DpnII", "across_HinfI", "trans"])
        self.assertEqual([row[-2:] for row in accepted[:2]], [["0", "1"], ["1", "2"]])
        policy["remove_adjacent_inward"] = True
        self.assertEqual([row[0] for row in self.fixture.filter(records, policy, digest)[0]], ["trans"])

    def test_native_planner_accepts_global_and_per_row_protocol(self):
        config = {**self.c, "RUN_HIC": "true", "HIC_OUTPUT_SCOPE": "all"}
        global_plan = build_plan(config, self.fixture.sheet)
        self.assertEqual(len(global_plan["samples"]), 2)
        self.assertEqual(set(global_plan["entities"]), {"sample:row1", "sample:row2"})
        self.assertEqual(set(global_plan["hic_targets"]), {"sample:row1", "sample:row2"})
        for row in global_plan["samples"].values():
            self.assertEqual(row["policy"]["enzyme_definitions"], ARIMA_DEFINITIONS)
            self.assertNotIn("preprocessing", row)
        for row in self.fixture.rows:
            row["hic_protocol"] = "arima_2enzyme"
        self.fixture.write_sheet()
        row_plan = build_plan({**config, "HIC_PROTOCOL": "legacy"}, self.fixture.sheet)
        self.assertEqual([row["policy"] for row in global_plan["samples"].values()],
                         [row["policy"] for row in row_plan["samples"].values()])

    def test_protocol_change_reuses_mapping_signature_but_changes_accepted_signature(self):
        versions = {"bwa-mem2": "fixture", "pairtools": "1.1.3", "fastp": "fixture"}
        old_config = {**self.c, "HIC_PROTOCOL": "restriction_hic", "HIC_RESTRICTION_ENZYMES": "DpnII"}
        old = Runtime(build_plan(old_config, self.fixture.sheet), self.root / "config", versions=versions)
        new = Runtime(build_plan(self.c, self.fixture.sheet), self.root / "config", versions=versions)
        old_keys = PairProcessor(old).keys["sample:row1"]
        new_keys = PairProcessor(new).keys["sample:row1"]
        for stage in ("trim", "sorted", "dedup"):
            self.assertEqual(new_keys[stage], old_keys[stage], stage)
        self.assertNotEqual(new_keys["accepted"], old_keys["accepted"])

    def test_full_iupac_alphabet_matches_concrete_dna_and_not_reference_ambiguity(self):
        # Independent standard IUPAC base sets and complementary ambiguity
        # symbols: exercise every allowed concrete choice, not only N.
        bases = {"R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
                 "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
        complement = dict(zip("RYSWKMBDHVN", "YRSWMKVHDBN"))
        original = registry_data(ROOT / "config/hic_protocols.json")
        fasta = self.root / "iupac.fa"
        for symbol, choices in bases.items():
            partner = complement[symbol]
            motif = "G" + symbol + partner + "C"
            definition = {"motif": motif, "cut_offset": 1}
            data = copy.deepcopy(original)
            data["enzymes"]["Degenerate"] = definition
            registry_data(self.registry_file(data))
            concrete = ["G" + left + right + "C"
                        for left, right in itertools.product(choices, bases[partner])]
            # The final word contains reference ambiguity and must never count
            # as an additional recognition site, even when it spells the motif.
            sequence = "NNNNN".join(concrete + [motif]).lower()
            fasta.write_text(">iupac\n" + "\n".join(sequence) + "\n")
            expected = [0, *[1 + index*9 for index in range(len(concrete))], len(sequence)]
            with self.subTest(symbol=symbol, motif=motif):
                self.assertEqual(list(scan_fasta(fasta, {"Degenerate": definition}, chunk_size=3)),
                                 [("iupac", len(sequence), expected)])

    def test_direct_multi_motif_config_and_sheet_build_combined_digest(self):
        # GCWGC is the ApeKI recognition motif. This is a scanner fixture, not
        # an assertion that the Chen or Arima library contains this enzyme.
        config = {**self.c, "HIC_PROTOCOL": "restriction_hic",
                  "HIC_RESTRICTION_ENZYMES": "^gatc + g^cwgc"}
        plan = build_plan(config, self.fixture.sheet)
        expected = {"DpnII": {"motif": "GATC", "cut_offset": 0},
                    "motif_GCWGC_cut1": {"motif": "GCWGC", "cut_offset": 1}}
        self.assertEqual(plan["samples"]["row1"]["policy"]["enzyme_definitions"], expected)
        for row in self.fixture.rows:
            row["restriction_enzymes"] = "g^cwgc+DpnII"
        self.fixture.write_sheet()
        row_plan = build_plan({**config, "HIC_RESTRICTION_ENZYMES": "HaeIII"}, self.fixture.sheet)
        self.assertEqual(row_plan["samples"]["row1"]["policy"], plan["samples"]["row1"]["policy"])
        fasta = self.root / "apeki.fa"
        # GCTGC and GCAGC match W; GCSGC is ambiguous reference DNA and does not.
        sequence = "GCTGCNNGCAGCNNGCSGCNNGATC"
        fasta.write_text(">apeki\n" + "\n".join(sequence) + "\n")
        self.assertEqual(list(scan_fasta(fasta, expected, chunk_size=2)),
                         [("apeki", 25, [0, 1, 8, 21, 25])])

    def test_explicit_arima_patterns_canonicalize_to_names_and_reuse_stage_keys(self):
        named_policy = resolve_policy({"assay": "hic"}, self.c)
        for value in ("^GATC+G^ANTC", "g^antc+^gatc", "DpnII+G^ANTC+^GATC"):
            config = {**self.c, "HIC_RESTRICTION_ENZYMES": value}
            with self.subTest(value=value):
                self.assertEqual(resolve_policy({"assay": "hic"}, config), named_policy)
        versions = {"bwa-mem2": "fixture", "pairtools": "1.1.3", "fastp": "fixture"}
        named = Runtime(build_plan(self.c, self.fixture.sheet), self.root / "config", versions=versions)
        config = {**self.c, "HIC_RESTRICTION_ENZYMES": "^GATC+G^ANTC"}
        explicit = Runtime(build_plan(config, self.fixture.sheet), self.root / "config", versions=versions)
        self.assertEqual(PairProcessor(explicit).keys, PairProcessor(named).keys)
        # Profile-specific canonicalization must honor its preset even when
        # another registry enzyme has exactly the same recognition definition.
        data = registry_data(ROOT / "config/hic_protocols.json")
        data["enzymes"]["AaaGatc"] = {"motif": "GATC", "cut_offset": 0}
        self.assertEqual(resolve_policy({"assay": "hic"}, {
            **config, "HIC_PROTOCOL_REGISTRY": str(self.registry_file(data))}), named_policy)

    def test_direct_patterns_require_explicit_valid_within_motif_cut(self):
        config = {**self.c, "HIC_PROTOCOL": "restriction_hic"}
        invalid = ("GANTC", "GCWGC", "G^^ANTC", "G^[ACGT]TC", "G^AT.C", "G^ACTC",
                   "G^ARTC", "G^AUT C", "^", "^GATC,G^ANTC", "GANTC(1/5)")
        for value in invalid:
            with self.subTest(value=value), self.assertRaises(ValueError):
                resolve_policy({"assay": "hic", "restriction_enzymes": value}, config)


if __name__ == "__main__":
    unittest.main()
