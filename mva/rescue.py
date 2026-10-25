"""Offline allele-specific rescue triage; rankings are laboratory priorities only.

Patient coordinates and phenotype are read locally, never submitted to a service.
The sequence/exon calculation is deliberately separate from literature evidence.
"""
from __future__ import annotations

import copy
import csv
import hashlib
import json
from pathlib import Path
from urllib.parse import urlparse

from .track1 import GENETIC_CODE, GeneModel, reverse_complement


DISCLAIMER = ("Laboratory hypotheses only: no candidate has demonstrated rescue in this case. "
              "Rank is an auditable experimental priority, not efficacy, pathogenicity, or a treatment recommendation.")
GATES = ("causal_variant", "allele_resolved_rna", "protein_identity_and_function",
         "division_normalized_rescue", "normal_vs_tumor_effect", "exposure_feasibility")


def validate_config(config: dict) -> None:
    if config.get("schema_version") != 1:
        raise ValueError("Unsupported rescue evidence schema")
    sources = config.get("sources", {})
    if not sources:
        raise ValueError("Source ledger required")
    if set(config.get("readthrough_amino_acids", {})) != {"TAA", "TAG", "TGA"}:
        raise ValueError("Three stop codon insertion-context records required")
    if set(config.get("gate_reasons", {})) != set(GATES):
        raise ValueError("All experimental gate reasons required")
    for key, source in sources.items():
        parsed = urlparse(source.get("url", ""))
        if source.get("kind") not in {"primary_study", "regulatory"}:
            raise ValueError(f"{key}: primary or regulatory source required")
        if parsed.scheme != "https" or not parsed.netloc or not source.get("title") or not source.get("accessed_on"):
            raise ValueError(f"{key}: complete source provenance required")
    seen = set()
    for candidate in config.get("candidates", []):
        cid = candidate.get("id")
        if not cid or cid in seen:
            raise ValueError("Unique candidate ids required")
        seen.add(cid)
        for field in ("name", "mechanism", "study_context", "counter_evidence", "exposure", "tumor_normal", "falsification"):
            if not candidate.get(field):
                raise ValueError(f"{cid}: {field} required")
        if candidate.get("branch") not in {"readthrough", "readthrough_nmd", "stabilization"}:
            raise ValueError(f"{cid}: unknown branch")
        approval = candidate.get("approval", {})
        if approval.get("status") not in {"approved_marketed", "historical_discontinued", "expired", "unverified"}:
            raise ValueError(f"{cid}: invalid approval status")
        refs = candidate.get("source_ids", [])
        if not refs or any(ref not in sources for ref in refs):
            raise ValueError(f"{cid}: source references missing or unknown")
        if not any(sources[ref]["kind"] == "primary_study" for ref in refs):
            raise ValueError(f"{cid}: primary evidence required")
        approval_refs = approval.get("source_ids", [])
        if any(ref not in sources for ref in approval_refs):
            raise ValueError(f"{cid}: unknown approval reference")
        if approval["status"] != "unverified":
            if not approval_refs or not all(sources[ref]["kind"] == "regulatory" for ref in approval_refs):
                raise ValueError(f"{cid}: regulatory approval provenance required")
        if not all(approval.get(k) for k in ("jurisdiction", "formulation", "indication", "qualification")):
            raise ValueError(f"{cid}: precise approval scope required")
        if not isinstance(candidate.get("base_priority"), int):
            raise ValueError(f"{cid}: integer experimental-priority baseline required")
        for flag, penalty in candidate.get("phenotype_penalties", {}).items():
            if flag not in config.get("phenotype_flag_hpo", {}) or not isinstance(penalty, int) or penalty > 0:
                raise ValueError(f"{cid}: invalid phenotype penalty")
            if not candidate.get("phenotype_reasons", {}).get(flag):
                raise ValueError(f"{cid}: phenotype penalty explanation required")


def _spliced_offset(position: int, exons: list[tuple[int, int]], strand: int) -> int | None:
    preceding = 0
    for start, end in exons:
        if start <= position <= end:
            return preceding + (position - start if strand == 1 else end - position)
        preceding += end - start + 1
    return None


def stop_context(variant: dict, raw_model: dict, readthrough_map: dict) -> dict:
    """Recompute an SNV stop and its spliced distance to the final exon junction.

    Distances count bases between the PTC's final nucleotide and the junction.
    >55 nt is NMD-likely, 50--55 nt is a boundary, and <=49 nt predicts escape
    from this positional rule only. Early stops and intronless models are caveated.
    """
    result = {"variant_id": variant.get("variant_id"), "gene": variant.get("gene"),
              "protein_change": variant.get("protein_change"), "verified_stop": False,
              "nmd_prediction": "unknown", "limitations": []}
    try:
        model = GeneModel(raw_model)
        if variant.get("transcript_id") != model.transcript:
            raise ValueError("Transcript does not match local reference")
        if variant.get("quality_status") == "fail":
            raise ValueError("Candidate failed local quality checks")
        if not model.sequence_valid or len(model.cds) % 3:
            raise ValueError("CDS sequence/model length mismatch")
        pos = int(variant["pos"])
        offset = model.cds_offset(pos)
        ref, alt = str(variant["ref"]).upper(), str(variant["alt"]).upper()
        if len(ref) != 1 or len(alt) != 1 or offset is None or alt not in "ACGT":
            raise ValueError("Only validated coding SNV stops supported")
        if model.strand == -1:
            ref, alt = reverse_complement(ref), reverse_complement(alt)
        if model.cds[offset] != ref:
            raise ValueError("Reference base mismatch")
        codon_offset = offset // 3 * 3
        old = model.cds[codon_offset:codon_offset + 3]
        new = old[:offset % 3] + alt + old[offset % 3 + 1:]
        if GENETIC_CODE.get(old) in {None, "*"} or GENETIC_CODE.get(new) != "*":
            raise ValueError("Local sequence does not verify a newly gained stop")
        if variant.get("alternate_codon") not in {None, new}:
            raise ValueError("Candidate codon disagrees with local sequence")
        exons = sorted(model.exons, reverse=model.strand == -1)
        if any(left[1] >= right[0] for left, right in zip(sorted(model.exons), sorted(model.exons)[1:])):
            raise ValueError("Overlapping exons in transcript model")
        translation_first = model.coding_start if model.strand == 1 else model.coding_end
        translation_offset = _spliced_offset(translation_first, exons, model.strand)
        if translation_offset is None:
            raise ValueError("Translation start outside exon model")
        stop_start = translation_offset + codon_offset
        stop_end = stop_start + 2
        final_junction = sum(end - start + 1 for start, end in exons[:-1])
        distance = final_junction - (stop_end + 1) if len(exons) > 1 else None
        exon_index = next(i + 1 for i, (start, end) in enumerate(exons) if start <= pos <= end)
        prediction = ("intronless_rule_not_applicable" if distance is None else
                      "likely_nmd" if distance > 55 else
                      "boundary_uncertain" if distance >= 50 else "positional_escape_predicted")
        if codon_offset < 150 and prediction == "likely_nmd":
            prediction = "start_proximal_uncertain"
        if exons[exon_index - 1][1] - exons[exon_index - 1][0] + 1 > 407 and prediction == "likely_nmd":
            prediction = "long_exon_uncertain"
        insertion = readthrough_map.get(new, [])
        reference_length = len(model.cds) // 3 - int(GENETIC_CODE.get(model.cds[-3:]) == "*")
        result.update({"verified_stop": True, "transcript_id": model.transcript,
                       "transcript_version": raw_model.get("transcript_version"),
                       "transcript_selection": raw_model.get("transcript_selection", "supplied local model"),
                       "cds_position": offset + 1, "protein_position": codon_offset // 3 + 1,
                       "reference_protein_length": reference_length,
                       "predicted_truncated_length_if_translated": codon_offset // 3,
                       "c_terminal_residues_missing_if_translated": reference_length - codon_offset // 3,
                       "reference_amino_acid": GENETIC_CODE[old], "reference_codon": old,
                       "stop_codon_dna": new, "stop_codon_rna": new.replace("T", "U"),
                       "downstream_base": model.cds[codon_offset + 3:codon_offset + 4] or None,
                       "reference_local_context": model.cds[max(0, codon_offset - 9):codon_offset + 12],
                       "exon_number": exon_index, "exon_count": len(exons),
                       "stop_end_spliced_position": stop_end + 1,
                       "final_junction_after_spliced_position": final_junction if len(exons) > 1 else None,
                       "bases_between_stop_end_and_final_junction": distance,
                       "nmd_prediction": prediction,
                       "reported_readthrough_amino_acids": insertion,
                       "native_amino_acid_in_reported_set": GENETIC_CODE[old] in insertion if insertion else None,
                       "source_ids": ["nmd_rule", "nmd_exceptions", "readthrough_identity"]})
        result["limitations"] = ["NMD is a positional prediction, not an RNA measurement; transcript use and reinitiation can alter it.",
                                 "Reported insertion identities come from other gene/reporter systems; no allele-specific yield or amino-acid distribution is predicted."]
    except (KeyError, TypeError, ValueError) as exc:
        result["limitations"].append(str(exc))
    return result


def build_variant_context(candidates: list[dict], reference: dict, phenotype: dict,
                          hypotheses: list[dict], evidence: dict) -> dict:
    models = {(g["symbol"], g.get("transcript_id")): g for g in reference.get("genes", [])}
    relevant = [v for v in candidates if v.get("gene") == "BUB1B"]
    stops = []
    for variant in relevant:
        if variant.get("consequence") == "stop_gained":
            model = models.get((variant["gene"], variant.get("transcript_id")))
            if model:
                stops.append(stop_context(variant, model, evidence["readthrough_amino_acids"]))
            else:
                stops.append({"variant_id": variant.get("variant_id"), "verified_stop": False,
                              "nmd_prediction": "unknown", "limitations": ["Matching transcript model unavailable"]})
    features = {x.get("hpo_id") for x in phenotype.get("features", []) if x.get("subject", "proband") == "proband"}
    flags = {name: bool(features.intersection(ids)) for name, ids in evidence["phenotype_flag_hpo"].items()}
    applicable_hypotheses = [h for h in hypotheses if h.get("gene") == "BUB1B"]
    phase = applicable_hypotheses[0].get("phase", "unknown") if applicable_hypotheses else "unknown"
    return {"gene": "BUB1B", "causal_status": "unvalidated", "phase": phase,
            "population_frequency_status": "unknown" if any(v.get("population_af") is None for v in relevant) else "available_requires_interpretation",
            "variant_count": len(relevant), "stops": stops,
            "missense_variants": [{"variant_id": v.get("variant_id"), "protein_change": v.get("protein_change"),
                                   "transcript_id": v.get("transcript_id"), "functional_status": "unmeasured"}
                                  for v in relevant if v.get("consequence") == "missense_variant"],
            "phenotype_flags": flags,
            "phenotype_note": "Flags describe locally documented features; a kidney feature does not establish impaired renal function. No phenotype was sent offsite.",
            "phase_scenarios": {
                "in_trans": "A readthrough product from the stop allele remains paired with the missense allele on the opposite homologue; its inserted residue still needs functional validation.",
                "in_cis": "Readthrough can produce a protein carrying both the inserted residue and the second missense change. The proposed recessive two-allele model would require revision.",
                "unknown": "Do not assume independent allele rescue until DNA/RNA phase is resolved."}}


def rank_candidates(config: dict, context: dict) -> list[dict]:
    """Uncalibrated, additive priorities plus hard approval and mechanism gates."""
    valid_stops = [s for s in context.get("stops", []) if s.get("verified_stop")]
    decisions = []
    for candidate in config["candidates"]:
        row = copy.deepcopy(candidate)
        parts = {"public_evidence_baseline": row["base_priority"]}
        reasons = []
        branch = row["branch"]
        applicable = bool(valid_stops) if branch.startswith("readthrough") else bool(context.get("missense_variants"))
        if not applicable:
            parts["missing_allele_substrate"] = -100
            reasons.append("No locally verified allele substrate for this branch.")
        if branch.startswith("readthrough") and valid_stops:
            parts["verified_stop_substrate"] = 3
            if any(s["nmd_prediction"] == "likely_nmd" for s in valid_stops):
                parts["likely_nmd_substrate_loss"] = -2 if branch == "readthrough" else 1
                reasons.append("Positional NMD prediction changes transcript availability; measure allele-specific RNA before screening.")
            if any(s.get("native_amino_acid_in_reported_set") is False for s in valid_stops):
                parts["non_native_insertion_risk"] = -3
                reasons.append("The native amino acid is outside the reported readthrough insertion set; full-length protein is insufficient evidence of rescue.")
            tested_codons = row.get("codon_evidence", {}).get("tested_favorable_codons", [])
            if tested_codons and all(s.get("stop_codon_dna") not in tested_codons for s in valid_stops):
                parts["outside_reported_codon_context"] = row["codon_evidence"]["outside_tested_codons_penalty"]
                reasons.append("No stop matches this candidate's favorable codon context in the cited reporter study; this is contextual counter-evidence, not proof of failure.")
            downstream = row.get("codon_evidence", {}).get("preferred_downstream_base")
            if downstream and all(s.get("downstream_base") != downstream for s in valid_stops):
                parts["outside_reported_flanking_context"] = -1
                reasons.append("Local +4 base differs from the strongest reporter contexts in the cited study; codon alone does not establish response.")
        if context.get("phase", "unknown") in {"unknown", "unresolved"}:
            parts["unknown_phase"] = -1
            reasons.append("Unknown phase changes the protein produced by stop bypass and the causal interpretation.")
        for flag, penalty in row.get("phenotype_penalties", {}).items():
            if context.get("phenotype_flags", {}).get(flag):
                parts[f"local_{flag}_concern"] = penalty
                reasons.append(row["phenotype_reasons"][flag])
        eligible = row["approval"]["status"] == "approved_marketed"
        if not eligible:
            reasons.append("Outside current approved-medication shortlist: " + row["approval"]["qualification"])
        state = "deferred_lab_hypothesis" if eligible and applicable else "mechanistic_comparator" if not eligible else "not_applicable"
        row.update({"state": state, "approved_shortlist_eligible": eligible and applicable,
                    "experimental_priority_score": sum(parts.values()), "score_components": parts,
                    "decision_reasons": reasons, "clinical_recommendation": False,
                    "gates": {gate: {"status": "pending", "reason": config["gate_reasons"][gate]} for gate in GATES}})
        if not applicable:
            row["gates"]["causal_variant"] = {"status": "fail", "reason": "Relevant allele substrate has not been locally verified."}
        decisions.append(row)
    decisions.sort(key=lambda row: (not row["approved_shortlist_eligible"], -row["experimental_priority_score"], row["id"]))
    rank = 0
    for row in decisions:
        if row["approved_shortlist_eligible"]:
            rank += 1
            row["lab_priority_rank"] = rank
        else:
            row["lab_priority_rank"] = None
    return decisions


def _read(path: Path) -> object:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def run(evidence_path: Path, output_dir: Path, variant_context_path: Path | None = None) -> dict:
    """Context JSON contains reference_path, candidates_path, phenotype_path,
    hypotheses_path. Relative paths resolve beside the context JSON, not cwd.
    Omitting context produces a comparator-only review with no active shortlist.
    """
    config = _read(evidence_path)
    validate_config(config)
    inputs = {"evidence": {"path": str(evidence_path.resolve()), "sha256": hashlib.sha256(evidence_path.read_bytes()).hexdigest()}}
    if variant_context_path:
        spec = _read(variant_context_path)
        loaded = {}
        for key in ("reference", "candidates", "phenotype", "hypotheses"):
            raw_path = spec.get(key + "_path")
            if not raw_path:
                raise ValueError(f"Context requires {key}_path")
            path = Path(raw_path)
            if not path.is_absolute():
                path = variant_context_path.parent / path
            loaded[key] = _read(path)
            inputs[key] = {"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        context = build_variant_context(loaded["candidates"], loaded["reference"], loaded["phenotype"], loaded["hypotheses"], config)
        context["selection_context"] = spec.get("selection_context", {})
    else:
        context = {"gene": "BUB1B", "causal_status": "unvalidated", "phase": "unknown", "stops": [],
                   "missense_variants": [], "phenotype_flags": {}, "context_note": "No local variant context supplied."}
    rows = rank_candidates(config, context)
    without_local_flags = {**context, "phenotype_flags": {}}
    unpenalized = rank_candidates(config, without_local_flags)
    ranked_ids = [r["id"] for r in rows if r["approved_shortlist_eligible"]]
    unpenalized_ids = [r["id"] for r in unpenalized if r["approved_shortlist_eligible"]]
    summary = {"strategy": "allele_specific_rescue", "reviewed_on": config["reviewed_on"],
               "approved_lab_shortlist_count": sum(r["approved_shortlist_eligible"] for r in rows),
               "validated_rescue_count": 0, "clinical_recommendation_count": 0,
               "comparator_count": sum(r["state"] == "mechanistic_comparator" for r in rows),
               "disclaimer": DISCLAIMER, "inputs": inputs,
               "ranking_sensitivity": {"with_local_flags": ranked_ids, "without_local_flags": unpenalized_ids,
                                       "order_changes": ranked_ids != unpenalized_ids,
                                       "interpretation": "Heuristic ranking depends on evidence/context weights and is not a validated estimate of benefit or safety."}}
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, payload in (("rescue_candidates.json", rows), ("variant_context.json", context),
                          ("rescue_summary.json", summary), ("rescue_evidence_snapshot.json", config)):
        (output_dir / name).write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    with (output_dir / "rescue_candidates.csv").open("w", encoding="utf-8", newline="") as handle:
        fields = ["lab_priority_rank", "id", "name", "state", "experimental_priority_score", "approval_status", "source_urls"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({**{key: row[key] for key in fields[:5]}, "approval_status": row["approval"]["status"],
                             "source_urls": " | ".join(config["sources"][sid]["url"] for sid in dict.fromkeys(row["source_ids"] + row["approval"]["source_ids"]))})
    _write_report(output_dir / "RESCUE_REVIEW.md", config, context, rows, summary)
    return summary


def _write_report(path: Path, config: dict, context: dict, rows: list[dict], summary: dict) -> None:
    sources = config["sources"]
    def links(ids):
        return "; ".join(f"[{sources[sid]['title']}]({sources[sid]['url']})" for sid in dict.fromkeys(ids))
    lines = ["# Alternative Track 2: restore the allele's missing function", "", DISCLAIMER, "",
             "This strategy asks whether a causal BUB1B allele can regain functional protein, rather than inferring a drug from downstream stress or transcription-factor enrichment. The causal model remains unvalidated.", "",
             f"{summary['approved_lab_shortlist_count']} approved-medication laboratory probes; zero demonstrated case-model rescues. Other compounds remain explicitly outside that shortlist.", "",
             "## Local allele context", ""]
    for stop in context.get("stops", []):
        if stop.get("verified_stop"):
            lines += [f"- {stop.get('protein_change')}: local {stop['transcript_id']} model confirms {stop['reference_codon']} to {stop['stop_codon_dna']} ({stop['stop_codon_rna']}); exon {stop['exon_number']}/{stop['exon_count']}; {stop['bases_between_stop_end_and_final_junction']} bases between stop end and final junction. Prediction: **{stop['nmd_prediction']}**.",
                      f"- If the RNA escapes decay, predicted truncation is {stop['predicted_truncated_length_if_translated']} amino acids versus a {stop['reference_protein_length']}-residue reference, removing {stop['c_terminal_residues_missing_if_translated']} C-terminal residues. Protein existence and function remain unmeasured.",
                      f"- Native residue {stop['reference_amino_acid']}; residues reported in other readthrough systems: {', '.join(stop['reported_readthrough_amino_acids'])}. Native residue present in that set: {stop['native_amino_acid_in_reported_set']}. This is an inference about a risk, not a measured product in this allele.",
                      "- " + " ".join(stop["limitations"]) + " " + links(stop["source_ids"])]
        else:
            lines += ["- Stop could not be independently recomputed: " + "; ".join(stop["limitations"])]
    lines += ["", f"Phase: {context.get('phase', 'unknown')}. Population-frequency interpretation and causal classification remain unresolved.", ""]
    if context.get("selection_context"):
        lines += ["Track 1 branch selection: " + json.dumps(context["selection_context"], ensure_ascii=False), ""]
    for phase, explanation in context.get("phase_scenarios", {}).items():
        lines += [f"- {phase}: {explanation}"]
    lines += ["", context.get("phenotype_note", context.get("context_note", "")), "",
              "## Ranked investigation and exclusions", "",
              "Scores order experiments only; point assignments are explicit heuristics, not fitted predictions. Approval is a hard eligibility gate and cannot be replaced by a high score.", "",
              "| Lab rank | Compound | Decision | Priority score | Approval scope |", "|---|---|---|---:|---|"]
    for row in rows:
        lines += [f"| {row['lab_priority_rank'] or 'Outside'} | {row['name']} | {row['state']} | {row['experimental_priority_score']} | {row['approval']['status']} |"]
    sensitivity = summary["ranking_sensitivity"]
    lines += ["", "Sensitivity to local phenotype penalties: " + " > ".join(sensitivity["with_local_flags"]) +
              "; with those penalties removed: " + " > ".join(sensitivity["without_local_flags"]) + ". " + sensitivity["interpretation"], ""]
    for row in rows:
        lines += ["", "### " + row["name"], "", row["mechanism"], "",
                  "**Observed evidence and setting.** " + row["study_context"] + " " + links(row["source_ids"]), "",
                  "**Counter-evidence.** " + row["counter_evidence"], "",
                  "**Approval.** " + row["approval"]["jurisdiction"] + "; " + row["approval"]["formulation"] + "; " + row["approval"]["indication"] + ". " + row["approval"]["qualification"] + " " + links(row["approval"]["source_ids"]), "",
                  "**Exposure feasibility.** " + row["exposure"], "",
                  "**Normal and tumor cells.** " + row["tumor_normal"], "",
                  "**Decision.** " + " ".join(row["decision_reasons"]), "",
                  "**Falsification.** " + row["falsification"], ""]
    lines += ["## Experiments that distinguish rescue from artifacts", ""]
    lines += [f"{i}. {step}" for i, step in enumerate(config["experiment_sequence"], 1)]
    lines += ["", "## Why this can improve the contest entry", "", config["comparison_claim"], "",
              "## Evidence boundary", "", "This is a reproducible computational and literature pass. No wet-laboratory results, expression measurements, exposure measurements, clinical benefit, or winning probability are claimed. All case-specific inputs remain local. Public literature was queried without patient coordinates or phenotype.", ""]
    path.write_text("\n".join(lines), encoding="utf-8")
