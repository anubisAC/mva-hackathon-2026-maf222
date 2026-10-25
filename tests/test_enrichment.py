"""Statistical reference checks and input-integrity tests for local enrichment."""

import json
import math
from pathlib import Path
import tempfile
import unittest

from mva.enrichment import benjamini_hochberg, hypergeometric_tail, run


class HypergeometricTests(unittest.TestCase):
    def test_known_exact_tail(self):
        # 2 draws from 10 objects, with 3 marked: P(X >= 2) = 3/45.
        self.assertAlmostEqual(hypergeometric_tail(2, 10, 3, 2), 1 / 15, places=13)

    def test_against_combinatorial_reference(self):
        for population, marked, draws in [(20, 7, 5), (40, 33, 20), (60, 5, 40)]:
            low = max(0, draws - (population - marked))
            high = min(marked, draws)
            for overlap in range(low, high + 2):
                exact = sum(
                    math.comb(marked, hits) * math.comb(population - marked, draws - hits)
                    for hits in range(overlap, high + 1)
                ) / math.comb(population, draws)
                self.assertAlmostEqual(hypergeometric_tail(overlap, population, marked, draws), exact, places=12)

    def test_tiny_tail_is_not_lost_by_subtraction(self):
        expected = math.comb(100, 20) / math.comb(10000, 20)
        actual = hypergeometric_tail(20, 10000, 100, 20)
        self.assertGreater(actual, 0)
        self.assertAlmostEqual(actual / expected, 1, places=9)

    def test_edges_and_invalid_parameters(self):
        self.assertEqual(hypergeometric_tail(0, 10, 0, 0), 1)
        self.assertEqual(hypergeometric_tail(2, 10, 1, 3), 0)
        self.assertEqual(hypergeometric_tail(5, 10, 10, 5), 1)
        for arguments in [(1, 0, 0, 0), (1, 10, 11, 2), (1, 10, 3, 11), (-1, 10, 3, 2), (1, 10, True, 2)]:
            with self.assertRaises(ValueError):
                hypergeometric_tail(*arguments)


class CorrectionTests(unittest.TestCase):
    def test_known_bh_values_and_monotonicity(self):
        p_values = [0.04, 0.001, 0.03, 0.008, 0.05]
        adjusted = benjamini_hochberg(p_values)
        for observed, expected in zip(adjusted, [0.05, 0.005, 0.05, 0.02, 0.05]):
            self.assertAlmostEqual(observed, expected)
        sorted_adjusted = [adjusted[index] for index in sorted(range(len(p_values)), key=p_values.__getitem__)]
        self.assertEqual(sorted_adjusted, sorted(sorted_adjusted))
        self.assertTrue(all(q >= p for p, q in zip(p_values, adjusted)))

    def test_ties_empty_and_invalid(self):
        self.assertEqual(benjamini_hochberg([]), [])
        self.assertEqual(benjamini_hochberg([0.1, 0.1]), [0.1, 0.1])
        for values in [[-0.1], [float("nan")], [float("inf")], [1.1], [True]]:
            with self.assertRaises(ValueError):
                benjamini_hochberg(values)


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.signature = self.root / "signature.json"
        self.gmt = self.root / "sets.gmt"
        self.universe = self.root / "universe.txt"
        self.output = self.root / "results"
        self.metadata = {
            "name": "Test contrast", "species": "Homo sapiens", "context": "case versus control",
            "direction": "up", "genes": ["A", "B"], "source": "Local experiment",
            "evidence_type": "measured_expression",
        }
        self.write_signature()
        self.gmt.write_text("Hits\tReference\tA\tB\tC\nNoHits\tReference\tD\tE\nOutside\tReference\tZ\n", encoding="utf-8")
        self.universe.write_text("\n".join("ABCDEFGHIJ") + "\n", encoding="utf-8")

    def write_signature(self):
        self.signature.write_text(json.dumps(self.metadata), encoding="utf-8")

    def analyze(self):
        return run(self.signature, self.gmt, self.universe, self.output)

    def test_all_sets_enter_correction_and_outputs_reopen(self):
        report = self.analyze()
        hit = next(row for row in report["results"] if row["set_name"] == "Hits")
        self.assertAlmostEqual(hit["p_value"], 1 / 15)
        self.assertAlmostEqual(hit["fdr"], 0.2)
        self.assertEqual(hit["context"], "case versus control")
        self.assertEqual(hit["direction"], "up")
        self.assertEqual(report["library"]["tested_set_count"], 3)
        outside = next(row for row in report["results"] if row["set_name"] == "Outside")
        self.assertEqual((outside["set_size"], outside["overlap"], outside["p_value"], outside["fdr"]), (0, 0, 1, 1))
        self.assertIsNone(outside["fold_enrichment"])
        for path in report["outputs"].values():
            self.assertTrue(Path(path).is_file())
        saved = json.loads(Path(report["outputs"]["json"]).read_text(encoding="utf-8"))
        self.assertEqual(saved, report)

    def test_exact_duplicates_removed_without_case_conversion(self):
        self.metadata["genes"] = ["A", " A ", "a"]
        self.write_signature()
        self.universe.write_text("A\na\nA\nB\nC\nD\nE\n", encoding="utf-8")
        self.gmt.write_text("Mixed\tReference\tA\tA\tB\nLowercase\tReference\ta\n", encoding="utf-8")
        report = self.analyze()
        self.assertEqual(report["signature"]["genes"], ["A", "a"])
        self.assertEqual(report["signature_duplicates_removed"], 1)
        self.assertEqual(report["universe"]["duplicates_removed"], 1)
        mixed = next(row for row in report["results"] if row["set_name"] == "Mixed")
        self.assertEqual(mixed["overlap_genes"], ["A"])
        self.assertEqual(mixed["duplicates_removed"], 1)

    def test_invalid_or_missing_background(self):
        self.universe.write_text("A\nC\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "outside the explicit universe"):
            self.analyze()
        self.assertFalse(self.output.exists())
        self.universe.write_text("# no genes\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Universe is empty"):
            self.analyze()
        self.universe.unlink()
        with self.assertRaises(FileNotFoundError):
            self.analyze()

    def test_empty_signature_and_no_overlap(self):
        self.metadata["genes"] = []
        self.write_signature()
        report = self.analyze()
        self.assertTrue(all(row["p_value"] == row["fdr"] == 1 for row in report["results"]))
        self.metadata["genes"] = ["J"]
        self.write_signature()
        report = self.analyze()
        self.assertTrue(any("No selected gene overlaps" in warning for warning in report["warnings"]))

    def test_literature_direction_is_only_a_hypothesis(self):
        self.metadata["evidence_type"] = "literature_gene_set"
        self.write_signature()
        report = self.analyze()
        self.assertEqual(report["signature"]["direction_status"], "hypothesis_only_unvalidated_direction")
        self.assertTrue(any("not evidence of expression reversal" in warning for warning in report["warnings"]))
        self.assertTrue(any("hypothesis only" in warning for warning in report["warnings"]))

    def test_bad_signature_and_duplicate_set_names(self):
        self.metadata["direction"] = "reversed"
        self.write_signature()
        with self.assertRaisesRegex(ValueError, "direction"):
            self.analyze()
        self.metadata["direction"] = "up"
        self.metadata["source"] = ""
        self.write_signature()
        with self.assertRaisesRegex(ValueError, "source"):
            self.analyze()
        self.metadata["source"] = "Reference"
        self.write_signature()
        self.gmt.write_text("Repeated\tDescription\tA\nRepeated\tDescription\tB\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "Duplicate GMT"):
            self.analyze()


if __name__ == "__main__":
    unittest.main()
