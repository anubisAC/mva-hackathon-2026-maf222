"""Exercise the boundary between strategy B modules using synthetic alleles."""
import gzip
import json
import tempfile
import unittest
from pathlib import Path

from mva.genomewide import run as genomewide_run
from mva.intake import sha256
from mva.rescue import run as rescue_run
from mva.strategy_b import baseline_fingerprints, comparison, create_rescue_context, read, write


EVIDENCE = Path(__file__).resolve().parents[1] / "config" / "rescue_evidence.json"


class StrategyBIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.addCleanup(self.temporary.cleanup)

    def test_genome_context_rescue_and_comparison_use_compatible_schemas(self):
        baseline, output = self.root / "baseline", self.root / "alternative"
        models = self.root / "genes.jsonl.gz"
        gene = {"symbol": "BUB1B", "chrom": "1", "start": 100, "end": 1199,
                "strand": 1, "transcript_id": "SYNTHETIC", "transcript_version": "1",
                "exons": [{"start": 100, "end": 399}, {"start": 500, "end": 799},
                          {"start": 900, "end": 1199}],
                "translation_start": 100, "translation_end": 1199,
                "cds_sequence": "ATG" + "GCT" * 59 + "TTA" + "GCT" * 238 + "TAA"}
        with gzip.open(models, "wt", encoding="utf-8") as handle:
            handle.write(json.dumps({"assembly": "GRCh38"}) + "\n" + json.dumps(gene) + "\n")
        phenotype = {"features": [
            {"hpo_id": "HP:0000121", "label": "Synthetic kidney feature", "subject": "proband"},
            {"hpo_id": "HP:0200067", "label": "Synthetic family feature", "subject": "family_history"}]}
        write(baseline / "intake/phenotype.json", phenotype)
        write(baseline / "track1/qc.json", {"panel_gene_count": 7})
        write(baseline / "track1/hypotheses.json", [{"gene": "BUB1B", "type": "possible_compound_heterozygous",
                                                   "variant_ids": ["1:281:T:G", "1:500:G:T"]}])
        write(baseline / "track2/track2_candidates.json", {"counts": {"deferred": 3, "rejected": 1}})
        baseline_before = {str(path): path.read_bytes() for path in baseline.rglob("*") if path.is_file()}
        phenotype_reference = self.root / "phenotype_reference.json"
        write(phenotype_reference, {"genes": {"BUB1B": ["HP:0000121"]},
                                    "patient_features": [{"hpo_id": "HP:0000121"}], "family_features": []})
        vcf = self.root / "synthetic.vcf"
        vcf.write_text("##fileformat=VCFv4.2\n##reference=GRCh38\n"
                       "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsynthetic\n"
                       "1\t281\t.\tT\tG\t100\tPASS\t.\tGT:AD:DP:GQ\t0/1:15,15:30:99\n"
                       "1\t500\t.\tG\tT\t100\tPASS\t.\tGT:AD:DP:GQ\t0/1:15,15:30:99\n")
        genomewide_run(vcf, models, phenotype_reference, output / "track1")
        # With unknown inheritance, the single stop can have the higher mean score.
        self.assertEqual(read(output / "track1/genes.json")[0]["type"], "single_heterozygous_exploratory")
        spec_path = create_rescue_context(output, baseline, models, "BUB1B")
        spec = read(spec_path)
        self.assertEqual(spec["selection_context"]["selected_genotype_model"], "possible_compound_heterozygous")
        self.assertEqual(len(read(output / "context/candidates.json")), 2)
        summary = rescue_run(EVIDENCE, output / "track2", spec_path)
        context = read(output / "track2/variant_context.json")
        self.assertTrue(context["stops"][0]["verified_stop"])
        self.assertEqual(context["stops"][0]["stop_codon_rna"], "UGA")
        self.assertEqual(context["stops"][0]["bases_between_stop_end_and_final_junction"], 417)
        self.assertEqual(len(context["missense_variants"]), 1)
        self.assertEqual(summary["validated_rescue_count"], 0)
        result = comparison(output, baseline)
        self.assertEqual(result["records"], 2)
        self.assertEqual(result["baseline_gene_count"], 7)
        self.assertEqual(result["alternative_gene_count"], 1)
        self.assertEqual(result["approved_lab_probes"], summary["approved_lab_shortlist_count"])
        self.assertEqual(result["validated_rescues"], 0)
        self.assertIn("not blinded replication", (output / "COMPARISON.md").read_text())
        self.assertEqual({str(path): path.read_bytes() for path in baseline.rglob("*") if path.is_file()}, baseline_before)

    def test_baseline_fingerprint_rejects_changed_artifacts(self):
        baseline = self.root / "baseline"
        write(baseline / "track1/qc.json", {"synthetic": 1})
        write(baseline / "run_manifest.json", {"outputs_sha256": {
            "track1/qc.json": sha256(baseline / "track1/qc.json")}})
        before = baseline_fingerprints(baseline)
        self.assertEqual(len(before), 2)
        write(baseline / "track1/qc.json", {"synthetic": 2})
        with self.assertRaisesRegex(ValueError, "Baseline artifact changed"):
            baseline_fingerprints(baseline)


if __name__ == "__main__":
    unittest.main()
