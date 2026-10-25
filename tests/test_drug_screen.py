import itertools
import math
import unittest

from mva.drug_screen import parent_families, primary_genes, score_families, stratified_tail


def mechanism(gene, **changes):
    return {"id": gene, "target_id": gene, "target_type": "SINGLE PROTEIN", "target_tax_id": 9606,
            "direct_interaction": True, "genes": [gene], "action_type": "INHIBITOR", **changes}


CONFIG = {"degree_bin_upper_limits": [1, 2, 4, 8, 1000000], "checkpoint_support_genes": ["A"], "inhibitory_actions": ["INHIBITOR"]}


class DrugScreenTests(unittest.TestCase):
    def test_stratified_tail_matches_exhaustive_sampling(self):
        strata = ((4, 2, 2), (5, 1, 2))
        possibilities = [list(itertools.combinations(range(n), draw)) for n, _, draw in strata]
        outcomes = [sum(sum(g < k for g in pick) for pick, (_, k, _) in zip(picks, strata)) for picks in itertools.product(*possibilities)]
        for hits in range(5):
            self.assertAlmostEqual(stratified_tail(hits, strata), sum(x >= hits for x in outcomes) / len(outcomes))

    def test_single_stratum_matches_known_hypergeometric(self):
        self.assertAlmostEqual(stratified_tail(2, ((10, 3, 2),)), math.comb(3, 2)/math.comb(10, 2))
        self.assertEqual(stratified_tail(0, ((0, 0, 0),)), 1)

    def test_complex_nonhuman_indirect_and_ambiguous_targets_not_direct(self):
        self.assertEqual(primary_genes(mechanism("A")), {"A"})
        for change in ({"target_tax_id": 10090}, {"target_type": "PROTEIN COMPLEX"}, {"direct_interaction": False}, {"genes": ["A", "B"]}):
            self.assertEqual(primary_genes(mechanism("A", **change)), set())

    def test_salt_forms_do_not_multiply_hypotheses(self):
        catalog = {"drugs": [{"id": "P", "name": "Parent", "mechanisms": [mechanism("A")]},
                              {"id": "S", "parent_id": "P", "name": "Salt", "mechanisms": [mechanism("A")], "withdrawn_flag": True},
                              {"id": "N", "name": "Unknown", "mechanisms": []}]}
        families = parent_families(catalog)
        self.assertEqual(len(families), 2)
        row = next(f for f in families if f["id"] == "P")
        self.assertEqual(row["primary_target_genes"], ["A"])
        self.assertEqual(len(row["mechanisms"]), 1)
        self.assertTrue(row["withdrawal_flag_present"])

    def test_missing_target_is_not_negative_evidence_and_all_tests_enter_bh(self):
        catalog = {"drugs": [{"id": g, "name": g, "mechanisms": [mechanism(g)]} for g in "ABCDE"] + [{"id": "none", "mechanisms": []}]}
        ranked, summary = score_families(parent_families(catalog), {"A"}, CONFIG)
        self.assertEqual(summary["assessable_families"], 5)
        self.assertEqual(summary["unassessed_families"], 1)
        hit = next(r for r in ranked if r["id"] == "A")
        self.assertAlmostEqual(hit["nominal_p"], .2)
        self.assertAlmostEqual(hit["nominal_q"], 1)
        self.assertTrue(hit["direction_flags"])
        self.assertEqual(hit["benefit_prediction"], "not_established")
        missing = next(r for r in ranked if r["id"] == "none")
        self.assertIsNone(missing["degree_adjusted_p"])
        self.assertIsNone(missing["relevance_rank"])

    def test_empty_testable_module_fails_instead_of_claiming_no_drugs(self):
        with self.assertRaisesRegex(ValueError, "no genes"):
            score_families(parent_families({"drugs": [{"id": "A", "mechanisms": [mechanism("A")]}]}), {"missing"}, CONFIG)


if __name__ == "__main__":
    unittest.main()
