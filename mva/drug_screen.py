"""Broad local drug-target/pathway relevance screen, without efficacy claims."""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path

from .enrichment import benjamini_hochberg, hypergeometric_tail
from .intake import sha256


def read(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def write(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def flag(value):
    return value is True or value == 1 or str(value).lower() == "true"


def primary_genes(mechanism: dict) -> set[str]:
    if str(mechanism.get("target_tax_id")) != "9606":
        return set()
    if mechanism.get("target_type") != "SINGLE PROTEIN" or not flag(mechanism.get("direct_interaction")):
        return set()
    genes = set(mechanism.get("genes", []))
    # Multiple symbols may be aliases or unresolved components, not distinct bindings.
    return genes if len(genes) == 1 else set()


def parent_families(catalog: dict) -> list[dict]:
    buckets = defaultdict(list)
    ids = set()
    for drug in catalog["drugs"]:
        if drug["id"] in ids:
            raise ValueError("Duplicate catalog molecule ID")
        ids.add(drug["id"])
        buckets[drug.get("parent_id") or drug["id"]].append(drug)
    families = []
    for identifier, members in sorted(buckets.items()):
        preferred = next((d for d in members if d["id"] == identifier), members[0])
        unique = {}
        for member in members:
            for mechanism in member.get("mechanisms", []):
                key = (mechanism.get("id"), mechanism.get("target_id"), mechanism.get("action_type"), tuple(sorted(mechanism.get("genes", []))))
                unique[key] = mechanism
        mechanisms = list(unique.values())
        targets = set().union(*(primary_genes(m) for m in mechanisms)) if mechanisms else set()
        all_genes = {g for m in mechanisms if str(m.get("target_tax_id")) == "9606" for g in m.get("genes", [])}
        families.append({"id": identifier, "name": preferred.get("name") or identifier,
                         "member_ids": sorted(d["id"] for d in members),
                         "member_names": sorted({d.get("name") for d in members if d.get("name")}),
                         "primary_target_genes": sorted(targets), "contextual_target_genes": sorted(all_genes),
                         "mechanisms": mechanisms,
                         "withdrawal_flag_present": any(flag(d.get("withdrawn_flag")) for d in members),
                         "black_box_flag_present": any(flag(d.get("black_box_warning")) for d in members),
                         "veterinary_flag_present": any(flag(d.get("veterinary")) for d in members),
                         "current_regulatory_status": "requires_jurisdiction_formulation_and_current_status_review",
                         "warnings": [w for d in members for w in d.get("warnings", [])]})
    return families


def hypergeom_mass(n: int, k: int, draws: int) -> dict[int, float]:
    if not 0 <= k <= n or not 0 <= draws <= n:
        raise ValueError("Invalid hypergeometric stratum")
    if n == 0:
        return {0: 1.0}
    denominator = math.comb(n, draws)
    return {x: math.comb(k, x) * math.comb(n-k, draws-x) / denominator
            for x in range(max(0, draws-(n-k)), min(k, draws)+1)}


@lru_cache(maxsize=10000)
def stratified_distribution(strata: tuple[tuple[int, int, int], ...]) -> tuple[tuple[int, float], ...]:
    distribution = {0: 1.0}
    for n, k, draws in strata:
        combined = defaultdict(float)
        for total, probability in distribution.items():
            for hits, mass in hypergeom_mass(n, k, draws).items():
                combined[total+hits] += probability * mass
        distribution = dict(combined)
    return tuple(sorted(distribution.items()))


def stratified_tail(hits: int, strata: tuple[tuple[int, int, int], ...]) -> float:
    return min(1.0, sum(prob for count, prob in stratified_distribution(strata) if count >= hits))


def direction_flags(family: dict, config: dict) -> list[dict]:
    result = []
    checkpoint = set(config["checkpoint_support_genes"])
    inhibitory = set(config["inhibitory_actions"])
    for m in family["mechanisms"]:
        action = (m.get("action_type") or "").upper().replace("_", " ")
        direct = primary_genes(m)
        for gene in sorted(direct & checkpoint):
            if action in inhibitory:
                result.append({"type": "possible_conflict_with_checkpoint_restoration", "gene": gene,
                               "target_id": m.get("target_id"), "action": action,
                               "interpretation": "Inhibiting this checkpoint component may oppose normal-cell restoration; source evidence is not a patient-specific harm estimate."})
        if "tubulin" in (m.get("target_name") or "").lower() and str(m.get("target_tax_id")) == "9606":
            result.append({"type": "microtubule_perturbation_requires_mitotic_safety_review",
                           "target_id": m.get("target_id"), "action": action,
                           "interpretation": "A microtubule-targeting mechanism can change cell division; cytostasis or selective death must not be mistaken for restored segregation."})
    return result


def score_families(families: list[dict], module_genes: set[str], config: dict) -> tuple[list[dict], dict]:
    frequency = Counter(g for family in families for g in family["primary_target_genes"])
    universe = set(frequency)
    if not universe:
        raise ValueError("No unambiguous direct human single-protein targets in catalog")
    module = module_genes & universe
    if not module:
        raise ValueError("Disease module has no genes in the eligible target universe")
    limits = config["degree_bin_upper_limits"]
    bins = {gene: next(i for i, limit in enumerate(limits) if degree <= limit) for gene, degree in frequency.items()}
    bin_n = Counter(bins.values())
    bin_k = Counter(bins[g] for g in module)
    rows = []
    for family in families:
        row = {k: v for k, v in family.items() if k not in {"mechanisms", "warnings"}}
        targets = set(family["primary_target_genes"])
        hits = targets & module
        profile = Counter(bins[g] for g in targets)
        strata = tuple((bin_n[i], bin_k[i], profile[i]) for i in sorted(bin_n))
        row.update({"target_count": len(targets), "module_overlap": sorted(hits), "overlap_count": len(hits),
                    "contextual_module_overlap": sorted(set(family["contextual_target_genes"]) & module_genes),
                    "assessment_status": "assessable_target_relevance" if targets else "no_eligible_target_data",
                    "nominal_p": hypergeometric_tail(len(hits), len(universe), len(module), len(targets)) if targets else None,
                    "degree_adjusted_p": stratified_tail(len(hits), strata) if targets else None,
                    "target_annotation_degree_profile": dict(sorted(profile.items())),
                    "direction_flags": direction_flags(family, config), "benefit_prediction": "not_established"})
        rows.append(row)
    assessed = [r for r in rows if r["target_count"]]
    for p_key, q_key in (("nominal_p", "nominal_q"), ("degree_adjusted_p", "degree_adjusted_q")):
        for row, q in zip(assessed, benjamini_hochberg(r[p_key] for r in assessed)):
            row[q_key] = q
    rows.sort(key=lambda r: (r["degree_adjusted_p"] is None, r["degree_adjusted_p"] if r["degree_adjusted_p"] is not None else 2,
                             -r["overlap_count"], r["id"]))
    previous, score_rank = None, None
    for rank, row in enumerate(rows, 1):
        row.setdefault("nominal_q", None)
        row.setdefault("degree_adjusted_q", None)
        row["relevance_rank"] = rank if row["target_count"] else None
        tie_key = (row["degree_adjusted_p"], row["overlap_count"])
        if tie_key != previous:
            score_rank = rank
        row["tied_relevance_rank"] = score_rank if row["target_count"] else None
        previous = tie_key
    summary = {"drug_families": len(families), "assessable_families": len(assessed),
               "unassessed_families": len(families)-len(assessed), "eligible_human_target_genes": len(universe),
               "module_genes_total": len(module_genes), "module_genes_in_target_universe": len(module),
               "families_with_primary_overlap": sum(bool(r["overlap_count"]) for r in rows),
               "families_with_contextual_overlap": sum(bool(r["contextual_module_overlap"]) for r in rows),
               "nominal_q_below_0_05": sum(r["nominal_q"] < .05 for r in assessed),
               "degree_adjusted_q_below_0_05": sum(r["degree_adjusted_q"] < .05 for r in assessed),
               "families_with_direction_flags": sum(bool(r["direction_flags"]) for r in rows),
               "universe": sorted(universe), "tested_module_genes": sorted(module),
               "annotation_degree_bins": [{"upper_limit": limit, "genes": bin_n[i], "module_genes": bin_k[i]} for i, limit in enumerate(limits)]}
    return rows, summary


def csv_rows(path: Path, rows: list[dict]):
    fields = ["relevance_rank", "tied_relevance_rank", "id", "name", "member_ids", "assessment_status", "target_count", "module_overlap", "contextual_module_overlap", "nominal_p", "nominal_q", "degree_adjusted_p", "degree_adjusted_q", "direction_flags", "withdrawal_flag_present", "black_box_flag_present", "veterinary_flag_present", "current_regulatory_status", "benefit_prediction"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(row.get(key)) if isinstance(row.get(key), (list, dict)) else row.get(key) for key in fields})


def run(catalog_path: Path, pathways_path: Path, config_path: Path, output: Path) -> dict:
    from .drug_pathways import select_modules
    if output.exists():
        raise ValueError("Choose a fresh output directory; completed screens are preserved")
    catalog, reference, config = read(catalog_path), read(pathways_path), read(config_path)
    selection = select_modules(reference, config["mechanism_seed_genes"], min_genes=config["pathway_min_genes"], max_genes=config["pathway_max_genes"])
    families = parent_families(catalog)
    ranked, summary = score_families(families, set(selection["union_genes"]), config)
    output.mkdir(parents=True)
    write(output / "policy.json", config)
    write(output / "module_selection.json", selection)
    write(output / "drug_families.json", families)
    write(output / "ranked_drugs.json", ranked)
    csv_rows(output / "ranked_drugs.csv", ranked)
    sensitivities = []
    base = {row["id"]: row for row in ranked}
    for maximum in config["sensitivity_max_genes"]:
        alternative = select_modules(reference, config["mechanism_seed_genes"], min_genes=config["pathway_min_genes"], max_genes=maximum)
        rows, alt_summary = score_families(families, set(alternative["union_genes"]), config)
        sensitivities.append({"scenario": f"maximum_pathway_size_{maximum}", "summary": alt_summary,
                              "selected_modules": alternative["selected"],
                              "ranks": [{"id": r["id"], "rank": r["relevance_rank"], "tied_rank": r["tied_relevance_rank"], "base_rank": base[r["id"]]["relevance_rank"], "degree_adjusted_q": r["degree_adjusted_q"]} for r in rows]})
    write(output / "sensitivity.json", sensitivities)
    summary.update({"catalog_molecules": len(catalog["drugs"]), "conditional_mechanism_genes": config["mechanism_seed_genes"],
                    "hypothesis": "BUB1B function loss proxy; mutation-specific drug interactions have not been simulated",
                    "statistical_interpretation": config["statistical_scope"], "clinical_candidates_validated": 0})
    write(output / "summary.json", summary)
    inputs = [catalog_path, pathways_path, config_path]
    code = [Path(__file__), Path(__file__).with_name("drug_pathways.py"), Path(__file__).with_name("enrichment.py")]
    write(output / "run_manifest.json", {"status": "completed", "created_at_utc": datetime.now(timezone.utc).isoformat(),
          "input_sha256": {str(p.resolve()): sha256(p) for p in inputs}, "code_sha256": {str(p.resolve()): sha256(p) for p in code},
          "outputs_sha256": {p.name: sha256(p) for p in sorted(output.iterdir()) if p.is_file()}})
    return {key: value for key, value in summary.items() if key not in {"universe", "tested_module_genes", "annotation_degree_bins"}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, default=Path("references/drug_screen/catalog.json"))
    parser.add_argument("--pathways", type=Path, default=Path("references/drug_screen/reactome_human_v97.json"))
    parser.add_argument("--config", type=Path, default=Path("config/drug_screen.json"))
    parser.add_argument("--output", type=Path, default=Path("outputs/drug_screen_c/targets"))
    args = parser.parse_args()
    print(json.dumps(run(args.catalog, args.pathways, args.config, args.output), indent=2))
