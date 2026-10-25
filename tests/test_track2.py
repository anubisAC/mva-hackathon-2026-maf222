"""Decision-rule and provenance tests using the checked-in literature ledger."""

import copy
import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from mva.track2 import REQUIRED_GATES, run, triage_candidate, validate_config


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "config" / "track2_evidence.json"


class Track2Tests(unittest.TestCase):
    def setUp(self):
        self.config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))

    def test_seed_does_not_conflate_literature_with_validated_treatment(self):
        validate_config(self.config)
        decisions = {c["id"]: triage_candidate(c) for c in self.config["candidates"]}
        self.assertEqual(decisions["colchicine"]["state"], "rejected")
        for candidate_id in ("acetylcysteine", "sirolimus", "nicotinamide_mononucleotide"):
            self.assertEqual(decisions[candidate_id]["state"], "deferred")
        self.assertEqual(sum(c["state"] == "accepted" for c in decisions.values()), 0)

    def test_approval_is_not_sufficient_for_acceptance(self):
        candidate = self.config["candidates"][0]
        self.assertEqual(candidate["approval"]["status"], "approved")
        result = triage_candidate(candidate)
        self.assertEqual(result["state"], "deferred")
        self.assertIn("Unresolved evidence gate: genome_safety.", result["decision_reasons"])

    def test_mechanistic_harm_cannot_be_overridden_by_passing_other_gates(self):
        candidate = copy.deepcopy(self.config["candidates"][2])
        for gate in candidate["gates"].values():
            gate["status"] = "pass"
        candidate.pop("exclusion_reason")
        result = triage_candidate(candidate)
        self.assertEqual(result["state"], "rejected")
        self.assertIn("conflicts", " ".join(result["decision_reasons"]))

    def test_failed_genome_gate_rejects_compatible_candidate(self):
        candidate = copy.deepcopy(self.config["candidates"][0])
        candidate["direction"] = "compatible"
        candidate["gates"]["genome_safety"]["status"] = "fail"
        self.assertEqual(triage_candidate(candidate)["state"], "rejected")

    def test_acceptance_requires_every_gate_and_compatible_direction(self):
        # Synthetic branch exercises rules; it is not written to the evidence ledger.
        candidate = copy.deepcopy(self.config["candidates"][0])
        candidate["direction"] = "compatible"
        for gate in REQUIRED_GATES:
            candidate["gates"][gate]["status"] = "pass"
        self.assertEqual(triage_candidate(candidate)["state"], "accepted")
        for gate in REQUIRED_GATES:
            trial = copy.deepcopy(candidate)
            trial["gates"][gate]["status"] = "pending"
            self.assertEqual(triage_candidate(trial)["state"], "deferred", gate)

    def test_unknown_approval_remains_deferred_even_if_gates_pass(self):
        candidate = copy.deepcopy(self.config["candidates"][3])
        candidate["direction"] = "compatible"
        for gate in candidate["gates"].values():
            gate["status"] = "pass"
        self.assertEqual(triage_candidate(candidate)["state"], "deferred")

    def test_evidence_source_must_exist(self):
        self.config["candidates"][0]["source_ids"] = ["invented_reference"]
        with self.assertRaisesRegex(ValueError, "unknown source"):
            validate_config(self.config)

    def test_approval_cannot_be_claimed_from_a_cell_study(self):
        self.config["candidates"][0]["approval"]["source_ids"] = ["nac_replication"]
        with self.assertRaisesRegex(ValueError, "regulatory provenance"):
            validate_config(self.config)

    def test_completed_gates_require_provenance(self):
        self.config["candidates"][0]["gates"]["functional_rescue"]["status"] = "pass"
        with self.assertRaisesRegex(ValueError, "evidence references required"):
            validate_config(self.config)

    def test_duplicate_candidate_ids_rejected(self):
        self.config["candidates"].append(copy.deepcopy(self.config["candidates"][0]))
        with self.assertRaisesRegex(ValueError, "unique"):
            validate_config(self.config)

    def test_run_is_offline_reproducible_and_preserves_uncertainty(self):
        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            with patch("socket.socket", side_effect=AssertionError("Network access forbidden")):
                result = run(CONFIG_PATH, Path(first))
                result2 = run(CONFIG_PATH, Path(second))
            self.assertEqual(result, result2)
            self.assertFalse(result["submission_ready"])
            self.assertEqual(result["counts"], {"accepted": 0, "deferred": 3, "rejected": 1})
            for filename in result["outputs"].values():
                self.assertEqual((Path(first) / filename).read_bytes(), (Path(second) / filename).read_bytes())
            report = (Path(first) / result["outputs"]["report"]).read_text(encoding="utf-8")
            self.assertIn("not efficacy scores", report)
            self.assertIn("No measured gene-expression signature", report)
            self.assertIn("UNCONFIRMED", report)
            self.assertIn("Falsification", report)
            self.assertIn("https://pubmed.ncbi.nlm.nih.gov/24477002/", report)
            with (Path(first) / result["outputs"]["csv"]).open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 4)
            self.assertEqual(rows[2]["state"], "rejected")
            saved = json.loads((Path(first) / result["outputs"]["json"]).read_text(encoding="utf-8"))
            self.assertEqual(saved["config_sha256"], result["config_sha256"])


if __name__ == "__main__":
    unittest.main()
