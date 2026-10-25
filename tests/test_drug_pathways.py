"""Reactome schema and selection tests use public-format synthetic fixtures."""

import copy
import json
import tempfile
import unittest
from pathlib import Path

from mva.drug_pathways import DISEASE_ROOT, FILES, RELEASE, SOURCE_BASE, download_file, normalize_reference, select_modules


def reference_fixture():
    names = "\n".join([
        "R-HSA-1\tTop\tHomo sapiens", "R-HSA-2\tIntermediate\tHomo sapiens",
        "R-HSA-3\tSpecific\tHomo sapiens", "R-HSA-4\tOther branch\tHomo sapiens",
        "R-MMU-9\tMouse\tMus musculus",
        f"{DISEASE_ROOT}\tDisease\tHomo sapiens",
        "R-HSA-6\tDisease child\tHomo sapiens",
        "R-HSA-7\tCell Cycle\tHomo sapiens",
    ])
    gmt = "\n".join([
        "Top\tR-HSA-1\tSEED\tA\tB\tC\tD\tEXTRA",
        "Intermediate\tR-HSA-2\tSEED\tA\tB\tC\tD\tEXTRA",
        "Specific\tR-HSA-3\tSEED\tA\tB\tC\tD\tD",
        "Other branch\tR-HSA-4\tSEED\tE\tF\tG\tH",
        "Disease child\tR-HSA-6\tSEED\tW\tX\tY\tZ",
        "Cell Cycle\tR-HSA-7\tSEED\tA\tB\tC\tD",
        "Mouse\tR-MMU-9\tSeed\tA\tB\tC\tD",
    ])
    relations = "\n".join([
        "R-HSA-1\tR-HSA-2", "R-HSA-2\tR-HSA-3",
        "R-HSA-1\tR-HSA-4", f"{DISEASE_ROOT}\tR-HSA-6",
        "R-HSA-1\tR-HSA-7",
    ])
    descriptions = "Identifier\tName\tSummation\nR-HSA-3\tSpecific\tA mechanism.\n"
    return gmt, names, relations, descriptions


class DrugPathwayTests(unittest.TestCase):
    def reference(self):
        return normalize_reference(*reference_fixture())

    def test_human_pathway_filter_deduplication_and_transitive_hierarchy(self):
        value = self.reference()
        rows = {row["id"]: row for row in value["pathways"]}
        self.assertNotIn("R-MMU-9", rows)
        self.assertEqual(rows["R-HSA-3"]["gene_count"], 5)
        self.assertEqual(rows["R-HSA-3"]["ancestor_ids"], ["R-HSA-1", "R-HSA-2"])
        self.assertEqual(rows["R-HSA-3"]["description"], "A mechanism.")
        self.assertFalse(rows[DISEASE_ROOT]["gmt_available"])

    def test_leaf_most_keeps_independent_branches_and_exposes_exclusions(self):
        value = select_modules(self.reference(), ["SEED"])
        self.assertEqual({row["id"] for row in value["selected"]}, {"R-HSA-3", "R-HSA-4"})
        excluded = {row["id"]: row["exclusion_reasons"] for row in value["excluded"]}
        self.assertIn("broader_ancestor_of_eligible_seed_module", excluded["R-HSA-2"])
        self.assertIn("disease_category_or_descendant", excluded["R-HSA-6"])
        self.assertIn("broad_category", excluded["R-HSA-7"])
        self.assertIn("top_level_category", excluded["R-HSA-1"])
        self.assertEqual(value["union_gene_count"], 9)
        self.assertEqual(value["seed_coverage"]["SEED"]["selected_module_count"], 2)

    def test_size_thresholds_are_inclusive_and_too_small_descendant_does_not_hide_parent(self):
        reference = self.reference()
        selected = select_modules(reference, ["SEED"], min_genes=5, max_genes=5)
        self.assertEqual(selected["selected_module_count"], 2)
        selected = select_modules(reference, ["SEED"], min_genes=6, max_genes=6)
        self.assertEqual([row["id"] for row in selected["selected"]], ["R-HSA-2"])

    def test_each_seed_retains_its_own_closest_module(self):
        reference = self.reference()
        parent = next(row for row in reference["pathways"] if row["id"] == "R-HSA-2")
        parent["gene_symbols"].append("SECOND")
        parent["gene_count"] += 1
        selected = select_modules(reference, ["SEED", "SECOND"])
        rows = {row["id"]: row for row in selected["selected"]}
        self.assertEqual(rows["R-HSA-2"]["selected_for_seed_genes"], ["SECOND"])
        self.assertEqual(rows["R-HSA-3"]["selected_for_seed_genes"], ["SEED"])

    def test_missing_seed_is_reported_without_fallback_gene_injection(self):
        value = select_modules(self.reference(), ["UNKNOWN"])
        self.assertEqual(value["selected"], [])
        self.assertEqual(value["union_genes"], [])
        self.assertEqual(value["seed_coverage"]["UNKNOWN"]["reference_module_count"], 0)

    def test_symbol_case_is_not_changed(self):
        value = select_modules(self.reference(), ["seed"])
        self.assertEqual(value["selected_module_count"], 0)

    def test_selection_does_not_mutate_reference(self):
        value = self.reference()
        original = copy.deepcopy(value)
        select_modules(value, ["SEED"])
        self.assertEqual(value, original)

    def test_duplicate_gmt_and_unknown_hierarchy_ids_fail(self):
        gmt, names, relations, descriptions = reference_fixture()
        with self.assertRaisesRegex(ValueError, "Duplicate GMT"):
            normalize_reference(gmt + "\n" + gmt.splitlines()[0], names, relations, descriptions)
        with self.assertRaisesRegex(ValueError, "unknown pathway"):
            normalize_reference(gmt, names, relations + "\nR-HSA-1\tR-HSA-UNKNOWN", descriptions)

    def test_cycles_fail_without_hanging(self):
        gmt, names, relations, descriptions = reference_fixture()
        with self.assertRaisesRegex(ValueError, "Cycle"):
            normalize_reference(gmt, names, relations + "\nR-HSA-3\tR-HSA-1", descriptions)

    def test_disambiguated_gmt_names_and_embedded_summation_tabs(self):
        gmt, names, relations, descriptions = reference_fixture()
        gmt = gmt.replace("Specific\tR-HSA-3", "Specific_3\tR-HSA-3")
        names = names.replace("Specific\tHomo", "Specific \tHomo")
        descriptions += "R-HSA-3\tSpecific\tSecond\tsummation.\n"
        row = next(row for row in normalize_reference(gmt, names, relations, descriptions)["pathways"] if row["id"] == "R-HSA-3")
        self.assertEqual(row["gmt_name"], "Specific_3")
        self.assertEqual(row["summations"], ["A mechanism.", "Second\tsummation."])

    def test_cached_hash_mismatch_fails_before_network(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary) / f"reactome_v{RELEASE}"
            directory.mkdir()
            name = FILES[0]
            (directory / name).write_bytes(b"tampered")
            (directory / (name + ".receipt.json")).write_text(json.dumps({"url": SOURCE_BASE + name, "sha256": "wrong"}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "fingerprint mismatch"):
                download_file(Path(temporary), name)

    def test_invalid_selection_parameters_fail(self):
        for seeds, low, high in [([], 5, 200), (["SEED"], 0, 200), (["SEED"], 6, 5)]:
            with self.assertRaises(ValueError):
                select_modules(self.reference(), seeds, low, high)


if __name__ == "__main__":
    unittest.main()
