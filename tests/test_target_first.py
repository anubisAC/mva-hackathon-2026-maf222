import copy
import json
import math
import tempfile
import unittest
from pathlib import Path

from mva.target_first import (
    activity_fraction, classify_match, map_catalog, pools, protein_at_time,
    profile_status, read, run, run_drug_profiles, sample_parameters, sha, simulate,
)


def parameters(scenario="trans_unstable_residual"):
    return {
        "scenario": scenario, "sample": 0, "stop_rna_fraction": .1,
        "stop_degradation": 2, "stop_checkpoint_function": .1,
        "stop_attachment_function": .1, "readthrough_degradation": 2,
        "readthrough_checkpoint_function": .8, "readthrough_attachment_function": .7,
        "stabilizable_fraction": .8, "balance_width": 1,
        "tumor_synthesis_multiplier": 3, "missense_degradation": 4,
        "missense_checkpoint_function": .8, "missense_attachment_function": .6,
    }


class TargetFirstTests(unittest.TestCase):
    def test_turnover_half_life_null_and_equilibrium(self):
        self.assertEqual(protein_at_time(0, 0, .5, 100), 0)
        self.assertAlmostEqual(protein_at_time(1, 0, 2, math.log(2)/2), .5)
        self.assertAlmostEqual(protein_at_time(2, 6, 3, 999), 2)
        self.assertEqual(protein_at_time(7, 1, 1, 0), 7)
        with self.assertRaises(ValueError):
            protein_at_time(1, 1, 0, 1)

    def test_exact_solution_agrees_with_independent_euler_integration(self):
        value, step = .7, .0001
        for _ in range(10000):
            value += step * (.4 - 1.3 * value)
        self.assertAlmostEqual(protein_at_time(.7, .4, 1.3, 1), value, places=4)

    def test_zero_intervention_and_vehicle_are_identical(self):
        p = parameters()
        original = copy.deepcopy(p)
        baseline = simulate(p, "vehicle", 0, 5)
        for intervention in [
            "residual_stabilization", "residual_function_restore", "nmd_attenuation",
            "readthrough", "local_pp2a_bypass", "aurora_b_reduction",
            "apc_cdc20_reduction", "bubr1_suppression", "wild_type_addback",
        ]:
            self.assertEqual(simulate(p, intervention, 0, 5), baseline)
        self.assertEqual(p, original)

    def test_amount_rescue_cannot_supply_absent_function(self):
        p = parameters()
        for key in list(p):
            if key.endswith("_function"):
                p[key] = 0
        baseline = simulate(p, "vehicle", 0, 5)
        stabilized = simulate(p, "residual_stabilization", .9, 5)
        self.assertGreater(stabilized["total_bubr1"], baseline["total_bubr1"])
        self.assertEqual(stabilized["checkpoint_equivalent_bubr1"], 0)
        self.assertEqual(stabilized["attachment_equivalent_bubr1"], 0)
        self.assertGreater(simulate(p, "wild_type_addback", 1, 5)["checkpoint_equivalent_bubr1"], .9)

    def test_nonresponsive_turnover_gives_no_stabilization_gain(self):
        p = parameters()
        p["stabilizable_fraction"] = 0
        self.assertEqual(simulate(p, "residual_stabilization", 1, 5), simulate(p, "vehicle", 0, 5))

    def test_readthrough_conserves_translation_when_degradation_equal(self):
        p = parameters()
        for strength in (0, .2, .9, 1):
            self.assertAlmostEqual(simulate(p, "readthrough", strength, 5)["total_bubr1"],
                                   simulate(p, "vehicle", 0, 5)["total_bubr1"])

    def test_cis_readthrough_retains_missense_damage(self):
        p = parameters("cis_with_reference_allele")
        synthesis, degradation, checkpoint, attachment, _ = pools(p)
        self.assertEqual(synthesis[1], 0)
        self.assertEqual(synthesis[0], .5)
        self.assertAlmostEqual(checkpoint[3], .8*.8)
        self.assertAlmostEqual(attachment[3], .7*.6)
        self.assertEqual(degradation[3], 4)
        # An absent full-length missense pool cannot be repaired into existence.
        self.assertEqual(simulate(p, "residual_function_restore", .9, 5), simulate(p, "vehicle", 0, 5))

    def test_overinhibition_and_exit_arrest_are_not_called_rescue(self):
        p = parameters()
        baseline = simulate(p, "vehicle", 0, 5)
        full_aurora = simulate(p, "aurora_b_reduction", 1, 5)
        self.assertEqual(full_aurora["attachment_balance_proxy"], 0)
        self.assertGreater(baseline["attachment_balance_proxy"], 0)
        full_apc = simulate(p, "apc_cdc20_reduction", 1, 5)
        self.assertIsNone(full_apc["checkpoint_braking_proxy"])
        self.assertEqual(full_apc["exit_capacity_proxy"], 0)

    def test_local_attachment_bypass_does_not_repair_checkpoint(self):
        p = parameters()
        baseline = simulate(p, "vehicle", 0, 5)
        bypass = simulate(p, "local_pp2a_bypass", 1, 5)
        self.assertEqual(bypass["checkpoint_equivalent_bubr1"], baseline["checkpoint_equivalent_bubr1"])
        self.assertEqual(bypass["attachment_balance_proxy"], 1)

    def test_tumor_and_normal_are_not_assumed_selective(self):
        p = parameters()
        normal = simulate(p, "residual_stabilization", .75, 5)
        same_state = simulate(p, "residual_stabilization", .75, 5, "same_residual_state")
        self.assertEqual(normal, same_state)
        overexpressed = simulate(p, "residual_stabilization", .75, 5, "residual_overexpression_state")
        self.assertAlmostEqual(overexpressed["total_bubr1"], normal["total_bubr1"]*3)

    def test_missing_functional_pharmacology_is_not_assigned_a_dose_response(self):
        self.assertIsNone(activity_fraction(1, None, .8))
        self.assertIsNone(activity_fraction(1, 1, None))
        self.assertAlmostEqual(activity_fraction(1, 1, .8), .4)
        self.assertEqual(activity_fraction(0, 1, .8), 0)
        self.assertAlmostEqual(activity_fraction(1e300, 1, .8), .8)
        with self.assertRaises(ValueError):
            activity_fraction(1, 0, .8)
        with self.assertRaises(ValueError):
            activity_fraction(float("nan"), 1, .8)

    def test_mechanism_direction_complex_and_species_guards(self):
        t = {"genes": ["MTOR"], "desired_action": "decrease"}
        m = {"genes": ["MTOR"], "target_tax_id": 9606, "target_type": "SINGLE PROTEIN",
             "direct_interaction": 1, "action_type": "INHIBITOR"}
        self.assertEqual(classify_match(t, m), "direct_direction_compatible")
        self.assertEqual(classify_match({**t, "desired_action": "increase"}, m), "opposite_direction")
        self.assertEqual(classify_match({**t, "desired_action": "increase_local_recruitment"}, m),
                         "binding_does_not_establish_required_function")
        self.assertEqual(classify_match(t, {**m, "target_type": "PROTEIN COMPLEX", "genes": ["MTOR", "FKBP1A"]}),
                         "complex_or_indirect_context_only")
        self.assertIsNone(classify_match(t, {**m, "target_tax_id": 10090}))
        self.assertEqual(classify_match(t, {**m, "variant_sequence": "mutant"}),
                         "variant_specific_record_requires_review")

    def test_parent_family_dedup_preserves_formulation_flags(self):
        target = {"id": "m", "target": "mTOR", "genes": ["MTOR"], "desired_action": "decrease",
                  "priority": "conditional", "arm": "both", "condition": "high activity"}
        mechanism = {"id": 3, "genes": ["MTOR"], "target_tax_id": 9606, "target_type": "SINGLE PROTEIN",
                     "direct_interaction": 1, "action_type": "INHIBITOR"}
        drugs = [{"id": "A", "parent_id": "A", "name": "base", "withdrawn_flag": False, "mechanisms": [mechanism]},
                 {"id": "B", "parent_id": "A", "name": "salt", "withdrawn_flag": True, "mechanisms": [mechanism]}]
        rows, coverage = map_catalog([target], {"drugs": drugs})
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["member_ids"], ["A", "B"])
        self.assertEqual(coverage[0]["direct_compatible_parent_families"], 1)
        self.assertTrue(rows[0]["member_flags"][1]["withdrawn_flag"])

    def test_profile_adapter_distinguishes_synthetic_binding_and_missing_data(self):
        policy = read("config/target_first_framework.json")
        empty = read("config/target_first_drug_profiles.template.json")["profiles"][0]
        self.assertEqual(profile_status(empty, policy), "missing_functional_response_or_exposure")
        profile = {**empty, "drug_id": "SYNTHETIC_CONTROL", "chemical_identity": "not a medicine",
                   "source_url": "synthetic_control", "assay_context": "test only",
                   "evidence_status": "synthetic_control", "concentration_unit": "arbitrary",
                   "free_concentrations": [1], "ec50": 1, "emax": 1, "hill": 1}
        self.assertEqual(profile_status(profile, policy), "synthetic_adapter_control_only")
        self.assertEqual(profile_status({**profile, "assay_type": "binding_kd"}, policy),
                         "binding_or_other_assay_not_a_functional_response")
        self.assertEqual(profile_status({**profile, "exposure_basis": "total_plasma"}, policy),
                         "exposure_basis_unresolved")
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            result = run_drug_profiles([empty, profile], policy, [parameters()], output)
            self.assertEqual(result["measured_profiles_simulated"], 0)
            self.assertEqual(result["synthetic_profiles_simulated"], 1)
            self.assertEqual(result["unresolved_profiles"], 1)
            self.assertEqual(result["profile_predictions"], 3)

    def test_small_reproducible_end_to_end_and_manifest(self):
        policy = read("config/target_first_framework.json")
        policy["simulation"]["samples_per_allele_scenario"] = 2
        self.assertEqual(sample_parameters(policy["simulation"]), sample_parameters(policy["simulation"]))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for key, value in [("phenotype", {"features": []}), ("variants", []), ("allele_context", {})]:
                path = root / (key + ".json")
                path.write_text(json.dumps(value))
                policy["patient_inputs"][key] = str(path)
            catalog = root / "catalog.json"
            catalog.write_text('{"drugs": []}')
            policy["catalog_mapping"]["input"] = str(catalog)
            config = root / "policy.json"
            config.write_text(json.dumps(policy))
            output = root / "run"
            summary = run(config, output)
            self.assertEqual(summary["virtual_experiments"], 4*2*3*10*5)
            self.assertEqual(summary["measured_profiles_simulated"], 0)
            self.assertIsNone(summary["drug_winner"])
            self.assertFalse(summary["phenotype_outcomes_fitted"])
            for name, digest in read(output / "run_manifest.json")["outputs_sha256"].items():
                self.assertEqual(sha(output / name), digest)
            with self.assertRaises(FileExistsError):
                run(config, output)


if __name__ == "__main__":
    unittest.main()
