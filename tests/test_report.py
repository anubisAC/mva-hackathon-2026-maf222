import json
import tempfile
import unittest
from pathlib import Path

from mva.intake import sha256
from mva.report import build
from mva.submission import validate


class ReportTests(unittest.TestCase):
    def write_json(self, path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value), encoding="utf-8")

    def fixture(self, path):
        config = path / "config.json"
        self.write_json(config, {"mechanisms": [{"gene": "TESTGENE", "description": "Synthetic mechanism"}]})
        self.write_json(path / "track1/qc.json", {"scan_complete": True, "records_seen": 2, "candidate_count": 1, "hypothesis_count": 1, "panel_gene_count": 1, "assembly": "GRCh38"})
        variant = {"variant_id": "1:10:A:G", "chrom": "1", "pos": 10, "ref": "A", "alt": "G", "gene": "TESTGENE", "protein_change": "p.K2R", "quality_status": "pass"}
        self.write_json(path / "track1/candidates.json", [variant])
        self.write_json(path / "track1/hypotheses.json", [{"gene": "TESTGENE", "type": "homozygous_alternate", "variant_ids": [variant["variant_id"]], "phase": "not_applicable", "population_frequency_status": "unknown"}])
        self.write_json(path / "track2/track2_candidates.json", {"config_sha256": sha256(config), "counts": {"accepted": 0, "deferred": 1, "rejected": 0}, "candidates": [{"name": "Synthetic comparator", "state": "deferred", "decision_reasons": ["Unverified mechanism."]}]})
        self.write_json(path / "run_manifest.json", {"status": "completed"})
        return config

    def test_derived_summary_and_draft_fail_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)
            config = self.fixture(path)
            result = build(path, config)
            text = Path(result["review"]).read_text(encoding="utf-8")
            self.assertIn("1-gene starter panel", text)
            self.assertIn("leading gene hypothesis is **TESTGENE**", text)
            self.assertIn("Synthetic comparator", text)
            self.assertNotIn("N-acetylcysteine", text)
            self.assertFalse(validate(Path(result["draft_csv"]))["valid"])
            manifest = json.loads((path / "run_manifest.json").read_text())
            self.assertEqual(manifest["outputs_sha256"]["REVIEW.md"], sha256(path / "REVIEW.md"))
            self.assertEqual(len(manifest["report_revisions"]), 1)

    def test_modified_ledger_cannot_replace_prior_evidence(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)
            config = self.fixture(path)
            original = config.read_bytes()
            config.write_text('{"mechanisms": []}', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "evidence differs"):
                build(path, config)
            config.write_bytes(original)
            build(path, config)
            config.write_text('{"mechanisms": []}', encoding="utf-8")
            build(path, config)
            self.assertIn("Synthetic mechanism", (path / "REVIEW.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
