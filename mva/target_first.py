"""Target-first hypothesis experiments, with unfitted molecular proxies.

No clinical outcomes, affinity predictions, or drug efficacy scores are inferred.
The default run changes virtual target activities, not administered drug doses.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import random
import statistics
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


READOUTS = (
    "total_bubr1", "checkpoint_equivalent_bubr1", "attachment_equivalent_bubr1",
    "attachment_balance_proxy", "checkpoint_braking_proxy", "exit_capacity_proxy",
)
INTERVENTIONS = {
    "vehicle", "residual_stabilization", "residual_function_restore",
    "nmd_attenuation", "readthrough", "local_pp2a_bypass",
    "aurora_b_reduction", "apc_cdc20_reduction", "bubr1_suppression",
    "wild_type_addback",
}
NEGATIVE_ACTIONS = {"INHIBITOR", "ANTAGONIST", "NEGATIVE ALLOSTERIC MODULATOR",
                    "INVERSE AGONIST", "BLOCKER", "DEGRADER"}
POSITIVE_ACTIONS = {"ACTIVATOR", "AGONIST", "PARTIAL AGONIST",
                    "POSITIVE ALLOSTERIC MODULATOR"}


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write(path, data):
    Path(path).write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def csv_write(path, rows):
    if not rows:
        Path(path).write_text("", encoding="utf-8")
        return
    with Path(path).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, ensure_ascii=False) if isinstance(v, (dict, list)) else v for k, v in row.items()})


def protein_at_time(initial, synthesis, degradation, time):
    """Exact solution of dP/dt = synthesis - degradation*P."""
    if min(initial, synthesis, time) < 0 or degradation <= 0:
        raise ValueError("Nonnegative amounts/rates/time and positive degradation required")
    steady = synthesis / degradation
    return steady + (initial - steady) * math.exp(-degradation * time)


def activity_fraction(concentration, ec50, emax, hill=1):
    """Assay response adapter. Missing functional pharmacology stays missing."""
    if any(v is None for v in (concentration, ec50, emax, hill)):
        return None
    if not all(math.isfinite(v) for v in (concentration, ec50, emax, hill)) or concentration < 0 or ec50 <= 0 or hill <= 0 or not 0 <= emax <= 1:
        raise ValueError("Incompatible concentration-response values")
    if concentration == 0:
        return 0.0
    # Stable Hill equation; concentration and EC50 must use the same unit.
    z = hill * math.log(concentration / ec50)
    logistic = 1 / (1 + math.exp(-z)) if z >= 0 else math.exp(z) / (1 + math.exp(z))
    return emax * logistic


def sample_parameters(simulation):
    """Latin-hypercube design samples a range; it is not a probability model."""
    rng = random.Random(simulation["seed"])
    n = simulation["samples_per_allele_scenario"]
    if n < 1:
        raise ValueError("At least one scenario sample is required")
    result = []
    for scenario in simulation["allele_scenarios"]:
        ranges = {**simulation["shared_ranges"], **simulation["scenario_ranges"][scenario]}
        columns = {}
        for key, (low, high) in sorted(ranges.items()):
            if low > high:
                raise ValueError("Reversed range")
            values = [low + (high - low) * (i + rng.random()) / n for i in range(n)]
            rng.shuffle(values)
            columns[key] = values
        for i in range(n):
            result.append({"scenario": scenario, "sample": i, **{k: v[i] for k, v in columns.items()}})
    return result


def pools(parameters, tumor_state="same_residual_state"):
    p = parameters
    cis = p["scenario"] == "cis_with_reference_allele"
    if tumor_state not in {"same_residual_state", "residual_overexpression_state"}:
        raise ValueError("Unknown tumor state")
    multiplier = p["tumor_synthesis_multiplier"] if tumor_state == "residual_overexpression_state" else 1
    # 0.5 synthesis per allele, WT degradation=1: two reference alleles give P=1.
    synthesis = [0.5 if cis else 0, 0 if cis else 0.5, 0.5 * p["stop_rna_fraction"], 0]
    synthesis = [s * multiplier for s in synthesis]
    degradation = [1, p["missense_degradation"], p["stop_degradation"], p["readthrough_degradation"]]
    cp = [1, p["missense_checkpoint_function"], p["stop_checkpoint_function"], p["readthrough_checkpoint_function"]]
    at = [1, p["missense_attachment_function"], p["stop_attachment_function"], p["readthrough_attachment_function"]]
    if cis:
        # Readthrough in cis still contains the missense change.
        degradation[3] = max(degradation[3], p["missense_degradation"])
        cp[3] *= cp[1]
        at[3] *= at[1]
    return synthesis, degradation, cp, at, multiplier


def simulate(parameters, intervention, strength, time, tumor_state="same_residual_state"):
    if intervention not in INTERVENTIONS or not 0 <= strength <= 1 or time < 0:
        raise ValueError("Invalid intervention, strength, or time")
    synthesis, degradation, cp, at, multiplier = pools(parameters, tumor_state)
    initial = [s / d for s, d in zip(synthesis, degradation)]
    cp = cp.copy()
    at = at.copy()
    if intervention == "residual_stabilization":
        # This lever changes turnover of available protein; it does not repair function.
        factor = 1 - strength * parameters["stabilizable_fraction"]
        degradation = [d * max(0.01, factor) for d in degradation]
    elif intervention == "residual_function_restore":
        cp[1] += strength * (1 - cp[1])
        at[1] += strength * (1 - at[1])
        if parameters["scenario"] == "cis_with_reference_allele":
            cp[3] = parameters["readthrough_checkpoint_function"] * cp[1]
            at[3] = parameters["readthrough_attachment_function"] * at[1]
    elif intervention == "nmd_attenuation":
        rna = parameters["stop_rna_fraction"]
        synthesis[2] = 0.5 * multiplier * (rna + strength * (1 - rna))
    elif intervention == "readthrough":
        # Translation is redistributed, not counted twice; NMD feedback is not modeled.
        production = synthesis[2]
        synthesis[2], synthesis[3] = production * (1 - strength), production * strength
    elif intervention == "bubr1_suppression":
        synthesis = [s * (1 - strength) for s in synthesis]
    elif intervention == "wild_type_addback":
        synthesis[0] += strength
    amounts = [protein_at_time(p, s, d, time) for p, s, d in zip(initial, synthesis, degradation)]
    checkpoint = sum(p * f for p, f in zip(amounts, cp))
    attachment = sum(p * f for p, f in zip(amounts, at))
    local_pp2a = attachment
    if intervention == "local_pp2a_bypass":
        # An idealized localization intervention, not global phosphatase activation.
        local_pp2a += strength * max(0, 1 - local_pp2a)
    aurora = 1 - strength if intervention == "aurora_b_reduction" else 1
    apc = 1 - strength if intervention == "apc_cdc20_reduction" else 1
    if local_pp2a <= 0 or aurora <= 0:
        balance = 0.0
    else:
        balance = math.exp(-parameters["balance_width"] * math.log(aurora / local_pp2a) ** 2)
    # These are phenomenological molecular proxies, not mitotic timing or viability.
    checkpoint_signal = (checkpoint / (0.5 + checkpoint)) / (1 / 1.5)
    braking = checkpoint_signal / apc if apc > 0 else None
    return {
        "total_bubr1": sum(amounts),
        "checkpoint_equivalent_bubr1": checkpoint,
        "attachment_equivalent_bubr1": attachment,
        "attachment_balance_proxy": balance,
        "checkpoint_braking_proxy": braking,
        "exit_capacity_proxy": apc,
    }


def classify_match(target, mechanism):
    genes = set(mechanism.get("genes", []))
    overlap = genes & set(target["genes"])
    if not overlap or str(mechanism.get("target_tax_id")) != "9606":
        return None
    if mechanism.get("variant_sequence"):
        return "variant_specific_record_requires_review"
    is_direct = (mechanism.get("target_type") == "SINGLE PROTEIN" and len(genes) == 1
                 and mechanism.get("direct_interaction") in (1, True))
    if not is_direct:
        return "complex_or_indirect_context_only"
    desired = target["desired_action"]
    if desired not in {"increase", "decrease"}:
        return "binding_does_not_establish_required_function"
    action = str(mechanism.get("action_type") or "").upper()
    sign = -1 if action in NEGATIVE_ACTIONS else 1 if action in POSITIVE_ACTIONS else None
    if sign is None:
        return "functional_direction_unresolved"
    return "direct_direction_compatible" if sign == (1 if desired == "increase" else -1) else "opposite_direction"


def map_catalog(targets, catalog):
    """Keep source mechanisms, exact family members and warnings; never rank efficacy."""
    hits = {}
    for target in targets:
        for drug in catalog["drugs"]:
            family = drug.get("parent_id") or drug["id"]
            for mechanism in drug.get("mechanisms", []):
                status = classify_match(target, mechanism)
                if status is None:
                    continue
                key = (target["id"], family, mechanism["id"])
                row = hits.setdefault(key, {
                    "target_id": target["id"], "target": target["target"],
                    "desired_action": target["desired_action"], "arm": target["arm"],
                    "parent_id": family, "member_ids": [], "member_names": [],
                    "mapping_status": status, "mechanism_id": mechanism["id"],
                    "action_type": mechanism.get("action_type"),
                    "target_chembl_id": mechanism.get("target_id"),
                    "target_type": mechanism.get("target_type"), "genes": mechanism.get("genes", []),
                    "direct_interaction": mechanism.get("direct_interaction"),
                    "mechanism_of_action": mechanism.get("mechanism_of_action"),
                    "mechanism_refs": mechanism.get("refs", []),
                    "member_flags": [],
                    "compound_url": f"https://www.ebi.ac.uk/chembl/compound_report_card/{family}/",
                    "simulation_status": "no_drug_response_parameters",
                })
                if drug["id"] not in row["member_ids"]:
                    row["member_ids"].append(drug["id"])
                    row["member_names"].append(drug["name"])
                    row["member_flags"].append({
                        "id": drug["id"], "withdrawn_flag": drug.get("withdrawn_flag"),
                        "withdrawal_warning_present": drug.get("withdrawal_warning_present"),
                        "veterinary": drug.get("veterinary"),
                        "current_regulatory_status": drug.get("current_regulatory_status"),
                        "black_box_warning": drug.get("black_box_warning"),
                    })
    rows = [hits[k] for k in sorted(hits)]
    coverage = []
    for target in targets:
        found = [r for r in rows if r["target_id"] == target["id"]]
        coverage.append({
            "target_id": target["id"], "target": target["target"], "priority": target["priority"],
            "arm": target["arm"], "desired_action": target["desired_action"],
            "any_parent_families": len({r["parent_id"] for r in found}),
            "direct_compatible_parent_families": len({r["parent_id"] for r in found if r["mapping_status"] == "direct_direction_compatible"}),
            "context_only_parent_families": len({r["parent_id"] for r in found if r["mapping_status"] == "complex_or_indirect_context_only"}),
            "status": "mechanism_coverage_only" if found else "no_gene_mapping_for_process" if not target["genes"] else "no_catalog_mechanism_match",
            "condition": target["condition"],
        })
    return rows, coverage


def profile_status(profile, policy):
    """Require functional response data; never substitute affinity for efficacy."""
    target = next((t for t in policy["targets"] if t["id"] == profile.get("target_id")), None)
    if target is None or profile.get("intervention") not in target["interventions"]:
        return "unresolved_target_or_intervention"
    required = ("drug_id", "chemical_identity", "source_url", "assay_context",
                "concentration_unit", "ec50", "emax", "hill", "free_concentrations")
    if any(profile.get(k) is None or profile.get(k) == "" or profile.get(k) == [] for k in required):
        return "missing_functional_response_or_exposure"
    if profile.get("assay_type") != "functional_target_response":
        return "binding_or_other_assay_not_a_functional_response"
    if profile.get("exposure_basis") != "free_intracellular":
        return "exposure_basis_unresolved"
    if profile.get("evidence_status") not in {"measured", "synthetic_control"}:
        return "evidence_status_unresolved"
    # Validate before producing any profile-level output.
    for concentration in profile["free_concentrations"]:
        activity_fraction(concentration, profile["ec50"], profile["emax"], profile["hill"])
    return "ready_for_unfitted_model" if profile["evidence_status"] == "measured" else "synthetic_adapter_control_only"


def run_drug_profiles(profiles, policy, samples, output):
    """Each profile supplies one action; no compound-wide selectivity is implied."""
    audits, rows = [], []
    for profile in profiles:
        status = profile_status(profile, policy)
        audits.append({**profile, "adapter_status": status,
                       "scope": "One specified molecular action; other drug effects and phenotype outcomes are not modeled."})
        if status not in {"ready_for_unfitted_model", "synthetic_adapter_control_only"}:
            continue
        sim = policy["simulation"]
        contexts = [("healthy_tissue_rescue", "same_residual_state")] + [
            ("tumor_control", state) for state in sim["tumor_states"]]
        for concentration in profile["free_concentrations"]:
            effect = activity_fraction(concentration, profile["ec50"], profile["emax"], profile["hill"])
            for p in samples:
                for arm, state in contexts:
                    baseline = simulate(p, "vehicle", 0, sim["time"], state)
                    result = simulate(p, profile["intervention"], effect, sim["time"], state)
                    row = {"drug_id": profile["drug_id"], "target_id": profile["target_id"],
                           "adapter_status": status, "concentration": concentration,
                           "concentration_unit": profile["concentration_unit"], "modeled_target_fraction": effect,
                           "scenario": p["scenario"], "sample": p["sample"], "arm": arm,
                           "tumor_state": state if arm == "tumor_control" else ""}
                    for metric, value in result.items():
                        row["delta_" + metric] = value - baseline[metric] if value is not None and baseline[metric] is not None else None
                    rows.append(row)
    write(output / "drug_profile_audit.json", audits)
    csv_write(output / "drug_profile_predictions.csv", rows)
    return {
        "measured_profiles_simulated": sum(p["adapter_status"] == "ready_for_unfitted_model" for p in audits),
        "synthetic_profiles_simulated": sum(p["adapter_status"] == "synthetic_adapter_control_only" for p in audits),
        "unresolved_profiles": sum(p["adapter_status"] not in {"ready_for_unfitted_model", "synthetic_adapter_control_only"} for p in audits),
        "profile_predictions": len(rows),
    }


def run(config, output, drug_profiles=None):
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError("Choose a fresh output directory")
    policy = read(config)
    sim = policy["simulation"]
    unknown = set(sim["interventions"]) - INTERVENTIONS
    if unknown:
        raise ValueError(f"Unknown interventions: {unknown}")
    inputs = {"policy": str(config), "catalog": policy["catalog_mapping"]["input"], **policy["patient_inputs"]}
    if drug_profiles is not None:
        inputs["drug_profiles"] = str(drug_profiles)
    input_hashes = {key: {"path": str(Path(value).resolve()), "sha256": sha(value)} for key, value in inputs.items()}
    source_hash = sha(__file__)
    output.mkdir(parents=True, exist_ok=True)
    write(output / "policy.json", policy)
    # Local patient observations constrain scope; none are converted into measured kinetics.
    write(output / "patient_context.json", {
        "phenotype": read(inputs["phenotype"]), "allele_context": read(inputs["allele_context"]),
        "source_variant_count": len(read(inputs["variants"])),
        "parameter_estimation": "No kinetic parameter inferred from HPO presence or variant class.",
    })
    samples = sample_parameters(sim)
    write(output / "parameter_design.json", samples)
    grouped = defaultdict(list)
    count = 0
    contexts = [("healthy_tissue_rescue", "same_residual_state")] + [
        ("tumor_control", state) for state in sim["tumor_states"]]
    print("Running virtual interventions with separate normal and tumor contexts", flush=True)
    with gzip.open(output / "virtual_experiments.csv.gz", "wt", encoding="utf-8", newline="") as stream:
        writer = None
        for parameters in samples:
            for arm, state in contexts:
                baseline = simulate(parameters, "vehicle", 0, sim["time"], state)
                for intervention in sim["interventions"]:
                    for strength in sim["strengths"]:
                        result = simulate(parameters, intervention, strength, sim["time"], state)
                        meta = {"scenario": parameters["scenario"], "sample": parameters["sample"],
                                "arm": arm, "tumor_state": state if arm == "tumor_control" else "",
                                "intervention": intervention, "strength": strength}
                        row = {**meta, **result}
                        for metric, value in result.items():
                            row["delta_" + metric] = value - baseline[metric] if value is not None and baseline[metric] is not None else None
                        if writer is None:
                            writer = csv.DictWriter(stream, fieldnames=list(row))
                            writer.writeheader()
                        writer.writerow(row)
                        key = (parameters["scenario"], arm, meta["tumor_state"], intervention, strength)
                        grouped[key].append(row)
                        count += 1
    summaries = []
    for key, rows in sorted(grouped.items()):
        row = dict(zip(("scenario", "arm", "tumor_state", "intervention", "strength"), key))
        row["design_samples"] = len(rows)
        for metric in READOUTS:
            values = [r["delta_" + metric] for r in rows if r["delta_" + metric] is not None]
            row[metric + "_delta_min"] = min(values) if values else None
            row[metric + "_delta_median"] = statistics.median(values) if values else None
            row[metric + "_delta_max"] = max(values) if values else None
        summaries.append(row)
    write(output / "virtual_experiment_summary.json", summaries)
    csv_write(output / "virtual_experiment_summary.csv", summaries)
    print("Matching the frozen target directions against the complete cached medicine catalog", flush=True)
    catalog = read(inputs["catalog"])
    matches, coverage = map_catalog(policy["targets"], catalog)
    write(output / "target_drug_matches.json", matches)
    csv_write(output / "target_drug_matches.csv", matches)
    write(output / "target_coverage.json", coverage)
    csv_write(output / "target_coverage.csv", coverage)
    profile_result = run_drug_profiles(read(drug_profiles)["profiles"] if drug_profiles else [], policy, samples, output)
    write(output / "calibration_requirements.json", {
        "state": "unfitted",
        "protein": ["allele-specific transcript abundance", "protein synthesis and half-life", "checkpoint and attachment function per protein amount"],
        "attachment": ["local B56 recruitment", "Aurora/PP2A balance versus misalignment", "cytokinesis and checkpoint-silencing responses"],
        "healthy_tissue": ["segregation errors per completed division", "productive growth or lineage-specific function"],
        "tumor": ["RMS molecular subtype and somatic state", "target dependence", "clonogenic survival across target perturbation strength"],
        "pharmacology": ["exact chemical identity", "signed functional target response", "EC50/Emax and assay units", "free intracellular exposure", "off-target effects"],
        "held_out_checks": ["BUBR1 addback and loss", "local B56 rescue versus global phosphatase modulation", "excess Aurora inhibition", "APC/C inhibition and context-dependent slippage"],
        "fit_policy": "Fit biological responses on one subset of public experiments and test on held-out perturbations before interpreting a drug ranking.",
    })
    result = {
        "status": "completed_unfitted_framework_pilot",
        "target_hypotheses": len(policy["targets"]), "allele_scenarios": len(sim["allele_scenarios"]),
        "parameter_designs": len(samples), "virtual_experiments": count,
        "interventions": len(sim["interventions"]), "intervention_strengths": len(sim["strengths"]),
        "drug_catalog_records": len(catalog["drugs"]),
        "catalog_parent_families": len({d.get("parent_id") or d["id"] for d in catalog["drugs"]}),
        "matched_parent_families": len({r["parent_id"] for r in matches}),
        "direct_direction_compatible_parent_families": len({r["parent_id"] for r in matches if r["mapping_status"] == "direct_direction_compatible"}),
        **profile_result, "drug_winner": None, "phenotype_outcomes_fitted": False,
        "interpretation": "Virtual molecular interventions test model dependencies. Catalog matches establish coverage only. Neither output is a drug efficacy ranking.",
    }
    write(output / "summary.json", result)
    if sha(__file__) != source_hash:
        raise RuntimeError("Code changed during run")
    write(output / "run_manifest.json", {
        "status": "complete", "completed_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": input_hashes, "code_sha256": {str(Path(__file__).resolve()): source_hash},
        "outputs_sha256": {p.name: sha(p) for p in sorted(output.iterdir()) if p.is_file()},
    })
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("config/target_first_framework.json"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--drug-profiles", type=Path, help="Optional measured functional target-response profiles; missing data are retained as unresolved")
    args = parser.parse_args()
    print(json.dumps(run(args.config, args.output, args.drug_profiles), indent=2))


if __name__ == "__main__":
    main()
