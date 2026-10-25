import gzip
import json
import tempfile
import unittest
from pathlib import Path

from mva.track1 import GeneModel, run
from mva.vcf import VCFReader


class Track1Tests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.addCleanup(self.directory.cleanup)

    def gene(self, **changes):
        gene = {"symbol": "TEST", "chrom": "1", "start": 100, "end": 130, "strand": 1,
                "transcript_id": "TEST_TRANSCRIPT", "exons": [{"start": 100, "end": 111}],
                "translation_start": 100, "translation_end": 111, "cds_sequence": "ATGGAATTTTAA"}
        gene.update(changes)
        return gene

    def write_vcf(self, rows, samples=("sample",)):
        path = self.root / "test.vcf.gz"
        header = ("##fileformat=VCFv4.2\n##contig=<ID=1,length=248956422>\n"
                  "##contig=<ID=2,length=242193529>\n"
                  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
        with gzip.open(path, "wt") as stream:
            stream.write(header + "\n".join(rows) + "\n")
        return path

    def analyze(self, rows, genes=None, **options):
        reference = self.root / "reference.json"
        reference.write_text(json.dumps({"assembly": "GRCh38", "genes": genes or [self.gene()]}))
        summary = run(self.write_vcf(rows), reference, self.root / "out", **options)
        candidates = json.loads((self.root / "out/candidates.json").read_text())
        hypotheses = json.loads((self.root / "out/hypotheses.json").read_text())
        return summary, candidates, hypotheses

    def test_multiallelic_alternate_depth_and_no_cohort_af_misuse(self):
        summary, candidates, hypotheses = self.analyze([
            "1\t103\t.\tG\tA,T\t100\tPASS\tAF=0.5,0.5\tGT:AD:DP:GQ\t1/2:0,12,18:30:99"])
        self.assertEqual(summary["records_seen"], 1)
        by_alt = {c["alt"]: c for c in candidates}
        self.assertEqual(by_alt["A"]["alternate_depth"], 12)
        self.assertEqual(by_alt["T"]["alternate_depth"], 18)
        self.assertEqual(by_alt["A"]["allele_fraction"], 0.4)
        self.assertEqual(by_alt["T"]["consequence"], "stop_gained")
        self.assertIsNone(by_alt["T"]["population_af"])
        self.assertEqual(by_alt["T"]["score_components"]["external_population_frequency"], 0)
        self.assertEqual(hypotheses[0]["phase"], "unknown")

    def test_missing_quality_retained_and_flagged(self):
        _, candidates, _ = self.analyze(["1\t103\t.\tG\tT\t.\t.\t.\tGT\t0/1"])
        candidate = candidates[0]
        self.assertIsNone(candidate["qual"])
        self.assertEqual(candidate["quality_status"], "review")
        self.assertIn("genotype_quality_missing", candidate["quality_flags"])
        self.assertIn("filter_not_applied_or_missing", candidate["quality_flags"])

    def test_negative_strand_codon(self):
        negative = self.gene(strand=-1)
        _, candidates, _ = self.analyze(["1\t108\t.\tC\tA\t100\tPASS\t.\tGT:AD:DP:GQ\t0/1:15,15:30:99"], [negative])
        candidate = candidates[0]
        self.assertEqual(candidate["cds_position"], 4)
        self.assertEqual(candidate["consequence"], "stop_gained")
        self.assertEqual(candidate["reference_codon"], "GAA")
        self.assertEqual(candidate["alternate_codon"], "TAA")

    def test_phased_gt_different_or_same_blocks_remain_unconfirmed(self):
        _, _, hypotheses = self.analyze([
            "1\t103\t.\tG\tT\t100\tPASS\t.\tGT:AD:DP:GQ:PS\t0|1:15,15:30:99:1",
            "1\t106\t.\tT\tC\t100\tPASS\t.\tGT:AD:DP:GQ:PS\t1|0:15,15:30:99:2"])
        self.assertEqual(len(hypotheses), 1)
        self.assertEqual(hypotheses[0]["phase"], "unknown")
        self.assertTrue(hypotheses[0]["phase_unconfirmed"])

    def test_supplied_cis_pair_is_excluded_but_evidence_is_retained(self):
        summary, candidates, hypotheses = self.analyze([
            "1\t103\t.\tG\tT\t100\tPASS\t.\tGT:AD:DP:GQ:PS\t0|1:15,15:30:99:1",
            "1\t106\t.\tT\tC\t100\tPASS\t.\tGT:AD:DP:GQ:PS\t0|1:15,15:30:99:1"])
        self.assertEqual(hypotheses, [])
        self.assertEqual(summary["qc"]["supplied_cis_pair_exclusions"], 1)
        self.assertEqual(candidates[0]["supplied_cis_partner_ids"], [candidates[1]["variant_id"]])

    def test_supplied_trans_phase_is_reported_with_validation_needed(self):
        _, _, hypotheses = self.analyze([
            "1\t103\t.\tG\tT\t100\tPASS\t.\tGT:AD:DP:GQ:PS\t0|1:15,15:30:99:1",
            "1\t106\t.\tT\tC\t100\tPASS\t.\tGT:AD:DP:GQ:PS\t1|0:15,15:30:99:1"])
        self.assertEqual(hypotheses[0]["phase"], "supplied_in_trans")
        self.assertTrue(hypotheses[0]["phase_unconfirmed"])

    def test_autosomal_haploid_call_is_not_assumed_hemizygous(self):
        _, candidates, hypotheses = self.analyze([
            "1\t103\t.\tG\tT\t100\tPASS\t.\tGT:AD:DP:GQ\t1:0,30:30:99"])
        self.assertEqual(candidates[0]["zygosity"], "haploid_alternate")
        self.assertIn("haploid_call_requires_ploidy_review", candidates[0]["quality_flags"])
        self.assertEqual(hypotheses, [])

    def test_scan_qc_and_prefix_are_explicit(self):
        rows = ["1\t103\t.\tG\tT\t100\tPASS\t.\tGT:AD:DP:GQ\t0/1:15,15:30:99",
                "2\t500\t.\tA\tG\t.\tLowQual\t.\tGT\t./."]
        summary, _, _ = self.analyze(rows)
        self.assertTrue(summary["scan_complete"])
        self.assertEqual(summary["qc"]["chromosome_counts"], {"1": 1, "2": 1})
        self.assertEqual(summary["qc"]["counts"]["missing_or_partial_genotype_records"], 1)
        prefix, _, _ = self.analyze(rows, max_records=1)
        self.assertFalse(prefix["scan_complete"])
        self.assertEqual(prefix["records_seen"], 1)

    def test_requires_single_sample(self):
        path = self.write_vcf([], samples=("one", "two"))
        with self.assertRaisesRegex(ValueError, "exactly one"):
            list(VCFReader(path).records())

    def test_codon_mapping_crosses_exon_on_reverse_strand(self):
        gene = GeneModel(self.gene(strand=-1, exons=[{"start": 100, "end": 105}, {"start": 120, "end": 125}],
                                   translation_start=100, translation_end=125))
        self.assertEqual(gene.cds_offset(125), 0)
        self.assertEqual(gene.cds_offset(120), 5)
        self.assertEqual(gene.cds_offset(105), 6)
        self.assertIsNone(gene.cds_offset(115))

    def test_wrong_reference_base_is_not_a_precise_protein_prediction(self):
        _, candidates, hypotheses = self.analyze(["1\t103\t.\tA\tT\t100\tPASS\t.\tGT:AD:DP:GQ\t1/1:0,30:30:99"])
        self.assertEqual(candidates[0]["consequence"], "coding_sequence_variant")
        self.assertIsNone(candidates[0]["protein_change"])
        self.assertIn("reference_sequence_mismatch", candidates[0]["annotation_flags"])
        self.assertEqual(hypotheses, [])

    def test_indel_reference_mismatch_cannot_be_called_frameshift_or_hypothesis(self):
        _, candidates, hypotheses = self.analyze([
            "1\t103\t.\tAC\tA\t100\tPASS\t.\tGT:AD:DP:GQ\t1/1:0,30:30:99"])
        self.assertEqual(candidates[0]["consequence"], "coding_sequence_variant")
        self.assertIn("reference_sequence_mismatch", candidates[0]["annotation_flags"])
        self.assertEqual(hypotheses, [])

    def test_missing_cds_does_not_promote_indel_to_hypothesis(self):
        _, candidates, hypotheses = self.analyze([
            "1\t103\t.\tGA\tG\t100\tPASS\t.\tGT:AD:DP:GQ\t1/1:0,30:30:99"],
            [self.gene(cds_sequence="")])
        self.assertEqual(candidates[0]["consequence"], "coding_sequence_variant")
        self.assertEqual(hypotheses, [])

    def test_low_quality_not_promoted_to_recessive_hypothesis(self):
        _, candidates, hypotheses = self.analyze(["1\t103\t.\tG\tT\t5\tLowQual\t.\tGT:AD:DP:GQ\t1/1:0,2:2:5"])
        self.assertEqual(candidates[0]["quality_status"], "fail")
        self.assertEqual(hypotheses, [])

    def test_indel_anchor_does_not_create_false_coding_precision(self):
        _, candidates, _ = self.analyze([
            "1\t111\t.\tA\tAC\t100\tPASS\t.\tGT:AD:DP:GQ\t0/1:15,15:30:99"])
        self.assertEqual(candidates[0]["consequence"], "coding_boundary_variant")
        self.assertIsNone(candidates[0]["protein_change"])

    def test_unaffected_anchor_outside_cds_does_not_hide_coding_deletion(self):
        gene = self.gene(start=99)
        _, candidates, _ = self.analyze([
            "1\t99\t.\tCA\tC\t100\tPASS\t.\tGT:AD:DP:GQ\t0/1:15,15:30:99"], [gene])
        self.assertEqual(candidates[0]["consequence"], "frameshift_variant")

    def test_duplicate_rows_do_not_form_a_compound_heterozygous_pair(self):
        row = "1\t103\t.\tG\tT\t100\tPASS\t.\tGT:AD:DP:GQ\t0/1:15,15:30:99"
        _, _, hypotheses = self.analyze([row, row])
        self.assertEqual(hypotheses, [])


if __name__ == "__main__":
    unittest.main()
