"""Protect direction, background, context aggregation, and catalog identity."""
import gzip
import json
import math
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from mva.drug_signatures import (DATASETS, aggregate_contexts, catalog_name_index, download_sources,
                                 compare_signed, match_catalog, read_pairs, run)


class SignatureTests(unittest.TestCase):
    def test_signed_score_matches_independent_pearson_and_reverses(self):
        x = [1, 1, -1, 0, 0, 0, -1, 0]
        y = [-1, 0, 1, 1, -1, 0, 0, 0]
        masks = lambda v: {"up": sum(1 << i for i, a in enumerate(v) if a == 1),
                           "down": sum(1 << i for i, a in enumerate(v) if a == -1)}
        mx, my = sum(x) / len(x), sum(y) / len(y)
        expected = -sum((a-mx)*(b-my) for a, b in zip(x, y)) / math.sqrt(
            sum((a-mx)**2 for a in x) * sum((b-my)**2 for b in y))
        score = compare_signed(masks(x), masks(y), 255)
        self.assertAlmostEqual(score["reversal_score"], expected)
        self.assertAlmostEqual(compare_signed(masks(x), masks(x), 255)["reversal_score"], -1)
        self.assertAlmostEqual(compare_signed(masks(x), masks([-a for a in x]), 255)["reversal_score"], 1)
        self.assertEqual(score["mimic_score"], -score["reversal_score"])

    def test_genes_outside_shared_universe_do_not_change_score(self):
        q = {"up": 3, "down": 12}
        d = {"up": 12, "down": 3}
        expected = compare_signed(q, d, 255)
        d["up"] |= 1 << 90
        self.assertEqual(compare_signed(q, d, 255), expected)
        self.assertIsNone(compare_signed(q, d, 0)["reversal_score"])

    def test_replicated_cell_does_not_outvote_other_cells_and_no_best_dose(self):
        values = [{"cell": "A", "condition": "low", "score": .6},
                  {"cell": "A", "condition": "high", "score": -.6},
                  {"cell": "B", "condition": "low", "score": -.2}]
        expected = aggregate_contexts(values)
        # Repeating identical assays in cell A does not increase its weight.
        repeated = values[:2] * 50 + values[2:]
        self.assertAlmostEqual(expected["score"], -.1)
        self.assertAlmostEqual(aggregate_contexts(repeated)["score"], -.1)
        self.assertEqual(expected["context_count"], 3)

    def test_alias_matching_exposes_ambiguity_and_does_not_strip_salts(self):
        catalog = {"drugs": [{"id": "CHEMBL1", "name": "Drug A", "synonyms": [{"molecule_synonym": "Alias"}]},
                             {"id": "CHEMBL2", "name": "Drug B", "synonyms": ["Alias"]}]}
        index = catalog_name_index(catalog)
        self.assertEqual(match_catalog("  dRuG   A ", index)["phase4_ids"], ["CHEMBL1"])
        self.assertEqual(match_catalog("Alias", index)["mapping_status"], "ambiguous_exact")
        self.assertEqual(match_catalog("Drug A hydrochloride", index)["mapping_status"], "unmatched")

    def test_complete_library_background_not_only_selected_proxy(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fixture.gmt"
            path.write_text("BUB1B::A up\t\tA\tB\nBUB1B::A down\t\tB\tC\nOTHER::A up\t\tD\nOTHER::A down\t\tE\nORPHAN::A up\t\tF\n")
            pairs, background, audit = read_pairs(path, {}, "BUB1B")
            self.assertEqual(background.bit_count(), 6)
            self.assertEqual(len(pairs), 1)
            self.assertEqual(pairs[0]["up"].bit_count(), 1)
            self.assertEqual(pairs[0]["down"].bit_count(), 1)
            self.assertEqual(audit["ambiguous_direction_genes_removed"], 1)

    def test_conflicting_duplicate_directions_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fixture.gmt"
            path.write_text("Drug::A up\t\tA\nDrug::A up\t\tB\nDrug::A down\t\tC\n")
            with self.assertRaisesRegex(ValueError, "Conflicting duplicate"):
                read_pairs(path, {})

    def test_cached_acquisition_never_queries_and_detects_modified_source(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for filename in DATASETS.values():
                (root / filename).write_bytes(b"public fixture data")
            with patch("urllib.request.urlopen", side_effect=AssertionError("Network must not be used")):
                receipts = download_sources(root)
                self.assertEqual(len(receipts), 3)
                self.assertEqual(download_sources(root), receipts)
                (root / DATASETS["drug"]).write_bytes(b"changed after receipt")
                with self.assertRaisesRegex(ValueError, "differs from its receipt"):
                    download_sources(root)

    def test_end_to_end_preserves_all_compounds_and_separates_proxies(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            libraries = {
                "drug": "Drug A::A375::::24H::1uM::r1 up\tID1\tC\tD\nDrug A::A375::::24H::1uM::r1 down\tID1\tA\tB\nUnmapped::OTHER::::6H::2uM up\tID2\tA\tB\nUnmapped::OTHER::::6H::2uM down\tID2\tC\tD\n",
                "knockout": "BUB1B::A375::::96H::::KO::Cas9 up\t\tA\tB\nBUB1B::A375::::96H::::KO::Cas9 down\t\tC\tD\n",
                "knockdown": "BUB1B::KOLF2.1J::hiPSCs::144H::::KD::CRISPRi up\t\tA\tB\nBUB1B::KOLF2.1J::hiPSCs::144H::::KD::CRISPRi down\t\tC\tD\n"}
            for key, data in libraries.items():
                if key == "drug":
                    data += "Drug A::A375::::24H::1uM::r2 up\tID3\tC\tD\nDrug A::A375::::24H::1uM::r2 down\tID3\tA\tB\n"
                with gzip.open(root / DATASETS[key], "wt") as stream:
                    stream.write(data)
            catalog = root / "catalog.json"
            catalog.write_text(json.dumps({"drugs": [{"id": "CHEMBL1", "name": "Drug A"}]}))
            output = root / "output"
            summary = run(root, catalog, output)
            rows = json.loads((output / "drug_rankings.json").read_text())
            self.assertEqual(summary["unique_drugs"], 2)
            self.assertEqual(summary["signature_comparisons"], 5)
            self.assertEqual(summary["unique_exact_phase4_name_matches"], 1)
            self.assertAlmostEqual(rows[0]["same_cell_score"], 1)
            self.assertEqual(rows[0]["source_identifier_count"], 2)
            self.assertFalse(rows[0]["identity_resolved_for_candidate_review"])
            self.assertEqual(summary["drug_names_with_multiple_source_ids"], 1)
            self.assertTrue((output / "run_manifest.json").exists())
            self.assertIsNone(rows[1]["same_cell_score"])
            self.assertAlmostEqual(rows[1]["hipsc_cross_cell_score"], -1)
            with self.assertRaises(FileExistsError):
                run(root, catalog, output)


if __name__ == "__main__":
    unittest.main()
