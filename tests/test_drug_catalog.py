"""Public-format fixtures guard broad catalog coverage and evidence semantics."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from mva.drug_catalog import collect, normalize, request_json


def molecule(identifier, **kwargs):
    return {"molecule_chembl_id": identifier, "pref_name": identifier, "max_phase": "4.0", **kwargs}


def target(identifier="TARGET", kind="SINGLE PROTEIN", tax_id=9606, symbols=("GENEA",)):
    return {"target_chembl_id": identifier, "tax_id": tax_id, "target_type": kind, "pref_name": kind,
            "target_components": [{"component_id": i, "accession": f"P{i}", "component_type": "PROTEIN",
                                   "target_component_synonyms": [{"syn_type": "GENE_SYMBOL", "component_synonym": symbol},
                                                                 {"syn_type": "GENE_SYMBOL_OTHER", "component_synonym": "NOT_CANONICAL"}]}
                                  for i, symbol in enumerate(symbols)]}


def mechanism(identifier=1, molecule_id="DRUG", **kwargs):
    return {"mec_id": identifier, "molecule_chembl_id": molecule_id, "target_chembl_id": "TARGET",
            "action_type": "INHIBITOR", "direct_interaction": 1, **kwargs}


class DrugCatalogTests(unittest.TestCase):
    def catalog(self, drugs, mechanisms=(), targets=(), warnings=()):
        return normalize(drugs, list(mechanisms), list(targets), list(warnings), {"test": True})

    def test_zero_target_molecules_remain_and_current_approval_not_inferred(self):
        catalog = self.catalog([molecule("UNKNOWN"), molecule("DRUG", withdrawn_flag=True)],
                               [mechanism()], [target()])
        self.assertEqual(catalog["summary"]["molecule_count"], 2)
        self.assertEqual(catalog["summary"]["without_human_target_genes"], 1)
        self.assertEqual(catalog["summary"]["withdrawn_flag_missing"], 1)
        self.assertEqual(catalog["summary"]["withdrawn_flag_true"], 1)
        self.assertEqual(catalog["summary"]["currently_approved_verified"], 0)
        self.assertTrue(all(d["current_regulatory_status"] == "unverified" for d in catalog["drugs"]))

    def test_nonhuman_target_is_not_converted_to_human_gene(self):
        catalog = self.catalog([molecule("DRUG")], [mechanism()], [target(tax_id=10090)])
        item = catalog["drugs"][0]["mechanisms"][0]
        self.assertEqual(item["genes"], [])
        self.assertIsNone(item["target_tax_id"])
        self.assertEqual(item["target_mapping_status"], "target_not_in_human_reference")

    def test_complex_gene_membership_preserves_ambiguity(self):
        catalog = self.catalog([molecule("DRUG")], [mechanism()], [target(kind="PROTEIN COMPLEX", symbols=("GENEA", "GENEB"))])
        item = catalog["drugs"][0]["mechanisms"][0]
        self.assertEqual(item["genes"], ["GENEA", "GENEB"])
        self.assertEqual(item["gene_specificity"], "component_membership_only")
        self.assertNotIn("NOT_CANONICAL", item["genes"])

    def test_parent_and_exact_mechanisms_merge_once_and_retain_join_origin(self):
        drugs = [molecule("SALT", molecule_hierarchy={"parent_chembl_id": "PARENT"})]
        records = [mechanism(1, "SALT", parent_molecule_chembl_id="PARENT"),
                   mechanism(2, "OTHER_SALT", parent_molecule_chembl_id="PARENT"),
                   mechanism(3, "UNRELATED")]
        items = self.catalog(drugs, records, [target()])["drugs"][0]["mechanisms"]
        self.assertEqual([item["id"] for item in items], [1, 2])
        self.assertEqual([item["match_basis"] for item in items], ["exact_molecule", "shared_parent"])

    def test_withdrawal_on_related_formulation_is_explicit(self):
        catalog = self.catalog([molecule("PARENT", withdrawn_flag=False)], warnings=[{
            "warning_id": 7, "molecule_chembl_id": "SALT", "parent_molecule_chembl_id": "PARENT",
            "warning_type": "Withdrawn", "warning_country": "Example"}])
        item = catalog["drugs"][0]
        self.assertFalse(item["withdrawn_flag"])
        self.assertTrue(item["withdrawal_warning_present"])
        self.assertEqual(item["warnings"][0]["molecule_chembl_id"], "SALT")

    def test_hgnc_fallback_and_missing_mechanism_gene_are_distinct(self):
        mapped = target()
        mapped["target_components"] = [{"target_component_synonyms": [], "target_component_xrefs": [
            {"xref_src_db": "HGNC", "xref_name": "GENEC"}]}]
        items = self.catalog([molecule("DRUG")], [mechanism(), mechanism(2, target_chembl_id=None)], [mapped])["drugs"][0]["mechanisms"]
        self.assertEqual(items[0]["genes"], ["GENEC"])
        self.assertEqual(items[0]["target_components"][0]["gene_mapping_basis"], "ChEMBL HGNC cross-reference")
        self.assertEqual(items[1]["target_mapping_status"], "no_target_record")

    def test_invalid_phase_or_duplicate_molecule_cannot_enter_universe(self):
        for rows in ([molecule("A"), molecule("A")], [molecule("A", max_phase="3")]):
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                self.catalog(rows)

    def pages(self, pages):
        calls = []
        def fetch(url, path, version):
            calls.append(url)
            return pages[len(calls) - 1], {"url": url, "sha256": "fixture", "chembl_version": version}
        return calls, fetch

    def page(self, identifiers, offset, total, next_page=None):
        return {"molecules": [molecule(i) for i in identifiers], "page_meta": {
            "offset": offset, "total_count": total, "next": next_page}}

    def test_pagination_completes_and_counts_every_id(self):
        pages = [self.page(["A", "B"], 0, 3, "/chembl/api/data/molecule.json?offset=2"), self.page(["C"], 2, 3)]
        calls, fetch = self.pages(pages)
        with tempfile.TemporaryDirectory() as temporary:
            result = collect("molecules", Path(temporary), "test", fetch=fetch)
        self.assertTrue(result["complete"])
        self.assertEqual(result["record_count"], 3)
        self.assertEqual(len(calls), 2)

    def test_truncation_duplicate_id_offset_change_and_cross_host_rejected(self):
        next_page = "/chembl/api/data/molecule.json?offset=1"
        cases = [[self.page(["A"], 0, 2)],
                 [self.page(["A"], 0, 2, next_page), self.page(["A"], 1, 2)],
                 [self.page(["A"], 0, 2, next_page), self.page(["B"], 2, 2)],
                 [self.page(["A"], 0, 2, next_page), self.page(["B"], 1, 3)],
                 [self.page(["A"], 0, 2, "https://example.com/next")]]
        for pages in cases:
            with self.subTest(pages=pages), tempfile.TemporaryDirectory() as temporary, self.assertRaises(ValueError):
                collect("molecules", Path(temporary), "test", fetch=self.pages(pages)[1])

    def test_valid_cache_requires_matching_url_release_and_hash(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "page.json"
            path.write_text('{"cached": true}', encoding="utf-8")
            receipt = {"url": "https://example.test/public.json", "chembl_version": "test",
                       "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            path.with_name("page.json.source.json").write_text(json.dumps(receipt), encoding="utf-8")
            data, returned = request_json(receipt["url"], path, "test")
            self.assertTrue(data["cached"])
            self.assertEqual(returned, receipt)

    def test_stale_or_corrupted_cache_is_not_returned(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "page.json"
            path.write_text('{"cached": true}', encoding="utf-8")
            correct = {"url": "https://example.test/public.json", "chembl_version": "test",
                       "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            for field, value in (("url", "https://other.test/"), ("chembl_version", "old"), ("sha256", "corrupted")):
                receipt = {**correct, field: value}
                path.with_name("page.json.source.json").write_text(json.dumps(receipt), encoding="utf-8")
                with self.subTest(field=field), patch("mva.drug_catalog.urllib.request.urlopen", side_effect=ValueError("fetch attempted")) as mocked:
                    with self.assertRaisesRegex(ValueError, "fetch attempted"):
                        request_json(correct["url"], path, "test")
                    mocked.assert_called_once()


if __name__ == "__main__":
    unittest.main()
