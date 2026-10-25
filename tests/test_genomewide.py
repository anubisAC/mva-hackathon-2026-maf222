import gzip
import json
import tempfile
import unittest
from pathlib import Path

from mva.genomewide import (ExonIndex, PhenotypeIndex, _clinical_score, build_hypotheses,
                           join_clinvar, load_reference, rank_genes, remove_clinvar_evidence, run)
from mva.track1 import GeneModel


class GenomeWideTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.addCleanup(self.directory.cleanup)

    @staticmethod
    def gene(symbol="TEST", **changes):
        result = {"symbol": symbol, "chrom": "1", "start": 100, "end": 300, "strand": 1,
                  "transcript_id": symbol + "_TRANSCRIPT", "exons": [{"start": 100, "end": 111}],
                  "translation_start": 100, "translation_end": 111, "cds_sequence": "ATGGAATTTTAA"}
        result.update(changes)
        return result

    @staticmethod
    def variant(gene="TEST", variant="1:103:G:T", **changes):
        result = {"gene": gene, "variant_id": variant, "chrom": "1", "consequence": "stop_gained",
                  "protein_change": "p.E2*", "quality_status": "pass", "annotation_flags": [],
                  "population_af": None, "population_af_source": None, "zygosity": "heterozygous",
                  "phase_set": None, "genotype": "0/1", "alt_index": 1, "triage_score": 9,
                  "score_components": {"predicted_consequence": 8, "quality": 1, "external_population_frequency": 0},
                  "clinical_score": 0, "supplied_cis_partner_ids": []}
        result.update(changes)
        return result

    def analyze(self, rows, genes=None, phenotype=None, clinvar=None, **options):
        reference = self.root / "reference.jsonl.gz"
        with gzip.open(reference, "wt") as stream:
            stream.write(json.dumps({"assembly": "GRCh38"}) + "\n")
            for gene in genes or [self.gene()]:
                stream.write(json.dumps(gene) + "\n")
        phenotype_path = self.root / "phenotype.json"
        phenotype_path.write_text(json.dumps(phenotype or {"genes": {}, "patient_features": []}))
        vcf = self.root / "test.vcf"
        vcf.write_text("##fileformat=VCFv4.2\n##reference=GRCh38\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n" + "\n".join(rows) + "\n")
        summary = run(vcf, reference, phenotype_path, self.root / "out", clinvar_path=clinvar, **options)
        with gzip.open(self.root / "out/candidates.json.gz", "rt") as stream:
            candidates = json.load(stream)
        return summary, candidates, json.loads((self.root / "out/genes.json").read_text())

    def test_interval_index_is_precise_with_overlapping_genes_and_long_deletion(self):
        genes = [GeneModel(self.gene()), GeneModel(self.gene("OTHER", start=109, exons=[{"start": 109, "end": 120}]))]
        index = ExonIndex(genes, bin_size=10)
        self.assertEqual(index.query("chr1", 110, 110), [0, 1])
        self.assertEqual(index.query("1", 200, 201), [])
        self.assertEqual(index.query("1", 1, 300000000), [0, 1])

    def test_whole_stream_scanned_but_deep_intronic_not_annotated(self):
        summary, candidates, ranked = self.analyze([
            "1\t103\t.\tG\tT\t100\tPASS\tAF=0.5\tGT:AD:DP:GQ\t0/1:15,15:30:99",
            "1\t200\t.\tA\tG\t100\tPASS\t.\tGT:AD:DP:GQ\t0/1:15,15:30:99",
            "2\t100\t.\tA\tG\t100\tPASS\t.\tGT:AD:DP:GQ\t0/1:15,15:30:99"])
        self.assertEqual(summary["records_seen"], 3)
        self.assertTrue(summary["scan_complete"])
        self.assertEqual(len(candidates), 1)
        self.assertIsNone(candidates[0]["population_af"])
        self.assertEqual(ranked[0]["type"], "single_heterozygous_exploratory")
        self.assertIn("single_allele_does_not_establish_recessive_disease", ranked[0]["review_flags"])

    def test_unselected_gene_can_rank_above_initial_gene(self):
        phenotype = {"patient_features": [{"hpo_id": "HP:1"}], "family_features": [],
                     "genes": {"NEW": ["HP:1"], "OLD": []}}
        variants = [self.variant("OLD"), self.variant("NEW", "1:203:G:T")]
        index = PhenotypeIndex(phenotype)
        hypotheses, _ = build_hypotheses(variants, index)
        self.assertEqual(rank_genes(hypotheses, index)[0]["gene"], "NEW")
        no_features = rank_genes(hypotheses, index, "no_phenotype")
        self.assertEqual(no_features[0]["score"], no_features[1]["score"])
        self.assertEqual(no_features[0]["score_rank"], no_features[1]["score_rank"])

    def test_feature_weights_use_all_genes_and_separate_family(self):
        index = PhenotypeIndex({"patient_features": [{"hpo_id": "HP:1"}, {"hpo_id": "HP:2"}],
                                "family_features": [{"hpo_id": "HP:F"}],
                                "genes": {"A": ["HP:1"], "B": ["HP:1", "HP:2"], "F": ["HP:F"]}})
        self.assertGreater(index.weights["HP:2"], index.weights["HP:1"])
        self.assertEqual(index.coverage("F"), (0, []))
        self.assertEqual(index.coverage("F", family=True), (1, ["HP:F"]))
        hypotheses, _ = build_hypotheses([self.variant("F")], index)
        self.assertEqual(hypotheses[0]["score_components"]["phenotype"], 0)
        self.assertEqual(hypotheses[0]["score_components"]["family_history"], .25)

    def test_clinvar_exact_allele_and_conflicts_no_positive_score(self):
        path = self.root / "clinvar.vcf"
        path.write_text("##reference=GRCh38\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                        "1\t103\t123\tG\tA\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel\n"
                        "1\t103\t124\tG\tT\t.\t.\tCLNSIG=Conflicting_classifications_of_pathogenicity;CLNREVSTAT=criteria_provided,_conflicting_classifications;CLNSIGCONF=Benign(1)|Pathogenic(1);AF_TGP=0.2\n")
        candidates = [self.variant()]
        qc = join_clinvar(path, candidates)
        self.assertEqual(qc["matched_unique_alleles"], 1)
        self.assertEqual(candidates[0]["clinvar"][0]["variation_id"], "124")
        self.assertEqual(candidates[0]["clinical_score"], 0)
        self.assertEqual(candidates[0]["population_af"], .2)
        self.assertIn("AF_TGP", candidates[0]["population_af_source"])

    def test_clinvar_review_status_and_benign_are_respected(self):
        def evidence(sig, review):
            return [{"clinical_significance": sig, "review_status": review, "allele_assignment_unambiguous": True}]
        self.assertEqual(_clinical_score(evidence("Pathogenic", "no_assertion_criteria_provided")), 0)
        self.assertEqual(_clinical_score(evidence("Pathogenic", "reviewed_by_expert_panel")), 4)
        self.assertEqual(_clinical_score(evidence("Benign/Likely_benign", "criteria_provided,_multiple_submitters,_no_conflicts")), -6)
        self.assertEqual(_clinical_score(evidence("Pathogenic", "reviewed_by_expert_panel") + evidence("Benign", "reviewed_by_expert_panel")), 0)

    def test_cis_pair_excluded_and_trans_still_unvalidated(self):
        index = PhenotypeIndex({})
        one = self.variant(phase_set="block", genotype="0|1")
        two = self.variant(variant="1:106:T:C", phase_set="block", genotype="0|1")
        hypotheses, qc = build_hypotheses([one, two], index)
        self.assertEqual(qc["supplied_cis_pair_exclusions"], 1)
        self.assertFalse(any(row["type"] == "possible_compound_heterozygous" for row in hypotheses))
        two["genotype"] = "1|0"
        hypotheses, _ = build_hypotheses([one, two], index)
        pair = next(row for row in hypotheses if row["type"] == "possible_compound_heterozygous")
        self.assertEqual(pair["phase"], "supplied_in_trans")
        self.assertTrue(pair["phase_unconfirmed"])

    def test_clinvar_ablation_preserves_alternative_within_gene(self):
        index = PhenotypeIndex({})
        variants = [self.variant(variant=f"1:{100+i}:G:T", triage_score=i, clinical_score=-6 if i > 6 else 4)
                    for i in range(1, 12)]
        hypotheses, qc = build_hypotheses(variants, index)
        self.assertGreater(qc["generated_hypothesis_count"], len(hypotheses))
        stripped = rank_genes(hypotheses, index, "no_clinvar")
        self.assertIn("1:111:G:T", stripped[0]["variant_ids"])

    def test_pair_effect_fields_follow_sorted_variant_ids(self):
        variants = [self.variant(variant="1:106:T:C", consequence="missense_variant", protein_change="p.F3L"),
                    self.variant(variant="1:103:G:T", consequence="stop_gained", protein_change="p.E2*")]
        hypotheses, _ = build_hypotheses(variants, PhenotypeIndex({}))
        pair = next(row for row in hypotheses if row["type"] == "possible_compound_heterozygous")
        self.assertEqual(pair["variant_ids"], ["1:103:G:T", "1:106:T:C"])
        self.assertEqual(pair["consequences"], ["stop_gained", "missense_variant"])
        self.assertEqual(pair["protein_changes"], ["p.E2*", "p.F3L"])

    def test_clinvar_ablation_rebuilds_common_frequency_eligibility(self):
        variant = self.variant(population_af=.3, population_af_source="ClinVar INFO/AF_TGP", triage_score=6,
                               population_af_without_clinvar=None, population_af_source_without_clinvar=None,
                               genomic_score_without_clinvar=9, frequency_score_without_clinvar=0)
        index = PhenotypeIndex({})
        self.assertEqual(build_hypotheses([variant], index)[0], [])
        stripped = remove_clinvar_evidence([variant])
        hypotheses, _ = build_hypotheses(stripped, index)
        self.assertEqual(hypotheses[0]["population_frequency_status"], "unknown")
        self.assertEqual(hypotheses[0]["score_components"]["genomic"], 9)
        self.assertEqual(variant["population_af"], .3)

    def test_known_recessive_gene_pair_not_hidden_by_single_stop(self):
        index = PhenotypeIndex({"gene_inheritance": {"TEST": ["autosomal_recessive"]}})
        hypotheses, _ = build_hypotheses([self.variant(), self.variant(variant="1:106:T:C", triage_score=5)], index)
        ranking = rank_genes(hypotheses, index)
        self.assertEqual(ranking[0]["type"], "possible_compound_heterozygous")
        self.assertLess(ranking[0]["score"], max(row["score"] for row in hypotheses))
        single = next(row for row in hypotheses if row["type"] == "single_heterozygous_exploratory")
        self.assertEqual(single["inheritance_compatibility"], "incomplete_known_recessive")

    def test_selenocysteine_model_cannot_generate_false_stop_hypothesis(self):
        _, candidates, ranked = self.analyze(["1\t103\t.\tG\tT\t100\tPASS\t.\tGT:AD:DP:GQ\t0/1:15,15:30:99"],
                                             genes=[self.gene(annotation_flags=["selenocysteine_requires_specialist_annotation"])])
        self.assertEqual(ranked, [])
        self.assertEqual(candidates[0]["hypothesis_exclusion_reason"], "selenocysteine_requires_specialist_annotation")

    def test_leave_one_out_recalculates_denominator(self):
        index = PhenotypeIndex({"patient_features": [{"hpo_id": "HP:1"}, {"hpo_id": "HP:2"}], "genes": {"TEST": ["HP:1"]}})
        hypotheses, _ = build_hypotheses([self.variant()], index)
        self.assertEqual(rank_genes(hypotheses, index, "omit:HP:1")[0]["score_components"]["phenotype"], 0)
        self.assertEqual(rank_genes(hypotheses, index, "omit:HP:2")[0]["score_components"]["phenotype"], 6)

    def test_qc_gating_unknown_af_and_genotype(self):
        variants = [self.variant("UNKNOWN"), self.variant("COMMON", population_af=.2),
                    self.variant("BAD", quality_status="fail"), self.variant("HAPLOID", zygosity="haploid_alternate"),
                    self.variant("MISMATCH", annotation_flags=["reference_sequence_mismatch"])]
        hypotheses, qc = build_hypotheses(variants, PhenotypeIndex({}))
        self.assertEqual({row["gene"] for row in hypotheses}, {"UNKNOWN"})
        self.assertEqual(sum(qc["allele_exclusions"].values()), 4)
        self.assertEqual(hypotheses[0]["population_frequency_status"], "unknown")

    def test_splice_only_variant_indexed_and_protein_not_invented(self):
        gene = self.gene(exons=[{"start": 100, "end": 105}, {"start": 120, "end": 125}], translation_end=125)
        _, candidates, _ = self.analyze(["1\t106\t.\tG\tA\t100\tPASS\t.\tGT:AD:DP:GQ\t0/1:15,15:30:99"], genes=[gene])
        self.assertEqual(candidates[0]["consequence"], "splice_donor_variant")
        self.assertIsNone(candidates[0]["protein_change"])


if __name__ == "__main__":
    unittest.main()
