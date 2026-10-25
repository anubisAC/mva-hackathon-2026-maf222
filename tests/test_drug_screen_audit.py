"""Independent checks of the drug-target statistical null and grouping invariants."""

import copy
import itertools
import random
import unittest

from mva.drug_screen import parent_families, score_families, stratified_tail


CONFIG = {"degree_bin_upper_limits": [1, 2, 4, 8, 1000000],
          "checkpoint_support_genes": [], "inhibitory_actions": []}


def mechanism(gene):
    return {"id": gene, "target_id": gene, "target_type": "SINGLE PROTEIN",
            "target_tax_id": 9606, "direct_interaction": True, "genes": [gene]}


def catalog():
    return {"drugs": [{"id": str(i), "name": str(i),
                       "mechanisms": [mechanism(g) for g in genes]}
                      for i, genes in enumerate(["AB", "AC", "AD", "BE", "CF", "DG", "HI", "JK"])]}


class DrugScreenAuditTests(unittest.TestCase):
    def test_stratified_null_matches_random_small_exhaustive_examples(self):
        rng = random.Random(771)
        for _ in range(50):
            strata = []
            combinations = []
            for _ in range(rng.randrange(1, 4)):
                n = rng.randrange(0, 7)
                k, draws = rng.randrange(n + 1), rng.randrange(n + 1)
                strata.append((n, k, draws))
                combinations.append(list(itertools.combinations(range(n), draws)))
            counts = [sum(sum(g < k for g in picked) for picked, (_, k, _) in zip(choice, strata))
                      for choice in itertools.product(*combinations)]
            for observed in range(sum(row[2] for row in strata) + 2):
                self.assertAlmostEqual(stratified_tail(observed, tuple(strata)),
                                       sum(count >= observed for count in counts) / len(counts), places=12)

    def test_single_popularity_bin_reduces_to_unstratified_null(self):
        config = {**CONFIG, "degree_bin_upper_limits": [1000000]}
        rows, _ = score_families(parent_families(catalog()), {"A", "B", "J"}, config)
        for row in rows:
            self.assertAlmostEqual(row["nominal_p"], row["degree_adjusted_p"], places=12)
            self.assertAlmostEqual(row["nominal_q"], row["degree_adjusted_q"], places=12)

    def test_duplicating_salt_form_does_not_change_annotation_degrees_or_scores(self):
        original = catalog()
        duplicate = copy.deepcopy(original["drugs"][0])
        duplicate.update({"id": "salt", "parent_id": "0", "name": "Salt of first molecule"})
        extended = copy.deepcopy(original)
        extended["drugs"].append(duplicate)
        before, summary_before = score_families(parent_families(original), {"A", "C"}, CONFIG)
        after, summary_after = score_families(parent_families(extended), {"A", "C"}, CONFIG)
        self.assertEqual(summary_before, summary_after)
        for first, second in zip(before, after):
            for key in ["id", "relevance_rank", "nominal_p", "degree_adjusted_p", "degree_adjusted_q"]:
                self.assertEqual(first[key], second[key])

    def test_catalog_and_mechanism_order_do_not_change_scores_or_tie_order(self):
        original = catalog()
        reordered = copy.deepcopy(original)
        reordered["drugs"].reverse()
        for drug in reordered["drugs"]:
            drug["mechanisms"].reverse()
        before, _ = score_families(parent_families(original), {"A", "B"}, CONFIG)
        after, _ = score_families(parent_families(reordered), {"A", "B"}, CONFIG)
        self.assertEqual(before, after)
        self.assertTrue(all(row["benefit_prediction"] == "not_established" for row in before))


if __name__ == "__main__":
    unittest.main()
