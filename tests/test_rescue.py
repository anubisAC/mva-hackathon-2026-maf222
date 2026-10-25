"""Sequence geometry and anti-overclaiming tests for allele-specific rescue."""
import copy
import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from mva.rescue import build_variant_context, rank_candidates, run, stop_context, validate_config


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "config" / "rescue_evidence.json"


def fixture(strand=1, codon_offset=180, exon_length=300):
    """Three exons, synthetic leucine-to-stop SNV, both genomic strands."""
    exons = [{"start": 100 + i * (exon_length + 100),
              "end": 99 + exon_length + i * (exon_length + 100)} for i in range(3)]
    sequence = "GCT" * exon_length
    sequence = sequence[:codon_offset] + "TTA" + sequence[codon_offset + 3:]
    transcript_exons = exons if strand == 1 else exons[::-1]
    offset = codon_offset + 1
    exon = transcript_exons[offset // exon_length]
    pos = exon["start"] + offset % exon_length if strand == 1 else exon["end"] - offset % exon_length
    model = {"symbol": "BUB1B", "chrom": "1", "start": 100, "end": exons[-1]["end"],
             "strand": strand, "transcript_id": "SYNTHETIC", "exons": exons,
             "translation_start": 100, "translation_end": exons[-1]["end"], "cds_sequence": sequence}
    variant = {"variant_id": "synthetic_stop", "gene": "BUB1B", "transcript_id": "SYNTHETIC",
               "chrom": "1", "pos": pos, "ref": "T" if strand == 1 else "A",
               "alt": "G" if strand == 1 else "C", "consequence": "stop_gained",
               "alternate_codon": "TGA", "quality_status": "pass", "protein_change": "synthetic"}
    return model, variant


class RescueTests(unittest.TestCase):
    def setUp(self):
        self.config = json.loads(EVIDENCE.read_text(encoding="utf-8"))

    def context(self, flags=True):
        model, variant = fixture()
        phenotype = {"features": [{"hpo_id": hpo, "subject": "proband"} for hpo in
                                  (["HP:0000121", "HP:0001508", "HP:0002859"] if flags else [])]}
        return build_variant_context([variant, {"gene": "BUB1B", "consequence": "missense_variant"}],
                                     {"genes": [model]}, phenotype, [], self.config)

    def test_recomputes_plus_and_minus_strand_stop_geometry(self):
        for strand in (1, -1):
            with self.subTest(strand=strand):
                model, variant = fixture(strand)
                result = stop_context(variant, model, self.config["readthrough_amino_acids"])
                self.assertTrue(result["verified_stop"])
                self.assertEqual(result["stop_codon_rna"], "UGA")
                self.assertEqual(result["bases_between_stop_end_and_final_junction"], 417)
                self.assertEqual(result["nmd_prediction"], "likely_nmd")
                self.assertFalse(result["native_amino_acid_in_reported_set"])
                self.assertEqual(result["predicted_truncated_length_if_translated"], 60)
                self.assertEqual(result["c_terminal_residues_missing_if_translated"], 240)

    def test_last_exon_stop_predicts_positional_escape(self):
        model, variant = fixture(codon_offset=660)
        result = stop_context(variant, model, self.config["readthrough_amino_acids"])
        self.assertEqual(result["exon_number"], 3)
        self.assertEqual(result["nmd_prediction"], "positional_escape_predicted")
        self.assertEqual(result["bases_between_stop_end_and_final_junction"], -63)

    def test_boundary_and_start_proximal_are_uncertain(self):
        for offset, expected in ((543, "boundary_uncertain"), (60, "start_proximal_uncertain")):
            model, variant = fixture(codon_offset=offset)
            result = stop_context(variant, model, self.config["readthrough_amino_acids"])
            self.assertEqual(result["nmd_prediction"], expected)

    def test_long_exon_is_not_unqualified_nmd(self):
        model, variant = fixture(exon_length=450)
        result = stop_context(variant, model, self.config["readthrough_amino_acids"])
        self.assertEqual(result["nmd_prediction"], "long_exon_uncertain")

    def test_intronless_model_has_no_final_junction_rule(self):
        model, variant = fixture()
        model["exons"] = [{"start": 100, "end": 999}]
        model["end"] = model["translation_end"] = 999
        result = stop_context(variant, model, self.config["readthrough_amino_acids"])
        self.assertEqual(result["nmd_prediction"], "intronless_rule_not_applicable")
        self.assertIsNone(result["bases_between_stop_end_and_final_junction"])

    def test_reference_mismatch_and_failed_quality_block_stop_substrate(self):
        model, variant = fixture()
        for mutation in ({"ref": "A"}, {"quality_status": "fail"}, {"transcript_id": "OTHER"}):
            result = stop_context({**variant, **mutation}, model, self.config["readthrough_amino_acids"])
            self.assertFalse(result["verified_stop"])
            self.assertEqual(result["nmd_prediction"], "unknown")

    def test_nmd_and_non_native_residue_change_priority(self):
        context = self.context()
        original = {r["id"]: r for r in rank_candidates(self.config, context)}
        context["stops"][0]["nmd_prediction"] = "positional_escape_predicted"
        context["stops"][0]["native_amino_acid_in_reported_set"] = True
        edited = {r["id"]: r for r in rank_candidates(self.config, context)}
        self.assertEqual(edited["gentamicin"]["experimental_priority_score"] - original["gentamicin"]["experimental_priority_score"], 5)
        self.assertNotIn("likely_nmd_substrate_loss", edited["gentamicin"]["score_components"])

    def test_local_flags_change_rank_without_inventing_renal_failure(self):
        flagged = rank_candidates(self.config, self.context(True))
        unflagged = rank_candidates(self.config, self.context(False))
        self.assertEqual(flagged[0]["id"], "azithromycin")
        self.assertEqual(unflagged[0]["id"], "gentamicin")
        self.assertIn("does not establish impaired renal function", self.context()["phenotype_note"])

    def test_family_history_does_not_trigger_proband_flags(self):
        model, variant = fixture()
        context = build_variant_context([variant], {"genes": [model]},
                                        {"features": [{"hpo_id": "HP:0000121", "subject": "family_history"}]}, [], self.config)
        self.assertFalse(context["phenotype_flags"]["kidney"])

    def test_higher_score_never_overrides_expired_or_discontinued_approval(self):
        self.config["candidates"][3]["base_priority"] = 10000
        rows = rank_candidates(self.config, self.context())
        eligible = [r for r in rows if r["approved_shortlist_eligible"]]
        self.assertEqual({r["id"] for r in eligible}, {"azithromycin", "gentamicin", "amikacin"})
        for row in rows:
            self.assertFalse(row["clinical_recommendation"])
            self.assertNotEqual(row["state"], "accepted")
            if row["id"] in {"amlexanox", "ataluren", "nmn_sirt2"}:
                self.assertIsNone(row["lab_priority_rank"])

    def test_missing_or_other_gene_context_cannot_activate_bubr1_shortlist(self):
        context = build_variant_context([{"gene": "OTHER", "consequence": "stop_gained"}],
                                        {"genes": []}, {}, [], self.config)
        self.assertFalse(any(r["approved_shortlist_eligible"] for r in rank_candidates(self.config, context)))

    def test_approval_requires_regulator_not_primary_paper(self):
        validate_config(self.config)
        self.config["candidates"][0]["approval"]["source_ids"] = ["pcd_context"]
        with self.assertRaisesRegex(ValueError, "regulatory"):
            validate_config(self.config)

    def test_non_tga_and_flanking_context_are_not_ignored(self):
        context = self.context()
        context["stops"][0]["stop_codon_dna"] = "TAG"
        context["stops"][0]["downstream_base"] = "A"
        row = next(r for r in rank_candidates(self.config, context) if r["id"] == "azithromycin")
        self.assertEqual(row["score_components"]["outside_reported_codon_context"], -2)
        self.assertEqual(row["score_components"]["outside_reported_flanking_context"], -1)

    def test_complete_offline_run_preserves_provenance_and_no_success_claim(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            model, variant = fixture()
            inputs = {"reference": {"genes": [model]}, "candidates": [variant],
                      "phenotype": {"features": [{"hpo_id": "HP:0000121", "subject": "proband"}]},
                      "hypotheses": [{"gene": "BUB1B", "phase": "unknown"}]}
            spec = {"selection_context": {"track1_top_gene": "OTHER", "mechanism_gene_rank": 2}}
            for key, value in inputs.items():
                (root / (key + ".json")).write_text(json.dumps(value))
                spec[key + "_path"] = key + ".json"
            context_path = root / "context.json"
            context_path.write_text(json.dumps(spec))
            with patch("socket.socket", side_effect=AssertionError("Network forbidden")):
                summary = run(EVIDENCE, root / "out", context_path)
            self.assertEqual(summary["approved_lab_shortlist_count"], 3)
            self.assertEqual(summary["validated_rescue_count"], 0)
            self.assertEqual(summary["clinical_recommendation_count"], 0)
            context = json.loads((root / "out" / "variant_context.json").read_text())
            self.assertEqual(context["selection_context"]["track1_top_gene"], "OTHER")
            self.assertEqual(len(summary["inputs"]["candidates"]["sha256"]), 64)
            with (root / "out" / "rescue_candidates.csv").open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 6)
            self.assertTrue(all("https://" in row["source_urls"] for row in rows))
            report = (root / "out" / "RESCUE_REVIEW.md").read_text(encoding="utf-8")
            self.assertIn("zero demonstrated case-model rescues", report)
            self.assertIn("in_cis", report)
            self.assertIn("Sensitivity", report)

    def test_run_without_context_has_empty_shortlist(self):
        with tempfile.TemporaryDirectory() as directory:
            result = run(EVIDENCE, Path(directory))
            self.assertEqual(result["approved_lab_shortlist_count"], 0)


if __name__ == "__main__":
    unittest.main()
