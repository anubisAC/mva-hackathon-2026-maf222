"""Run and compare an alternative strategy without changing baseline outputs."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

from .intake import sha256
from .report import ACKNOWLEDGEMENT


def read(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def write(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def fingerprints(paths):
    return {str(path.resolve()): sha256(path) for path in paths}


def baseline_fingerprints(directory: Path) -> dict:
    manifest = read(directory / "run_manifest.json")
    for name, expected in manifest["outputs_sha256"].items():
        if sha256(directory / name) != expected:
            raise ValueError(f"Baseline artifact changed since its run: {name}")
    return fingerprints([directory / "run_manifest.json", *[directory / name for name in manifest["outputs_sha256"]]])


def create_rescue_context(output: Path, baseline: Path, models: Path, gene: str) -> Path:
    candidates = read(output / "track1/candidates.json.gz")
    ranked = read(output / "track1/genes.json")
    biallelic = read(output / "track1/biallelic_genes.json")
    gene_row = next((r for r in ranked if r["gene"] == gene), None)
    hypothesis = next((r for r in biallelic if r["gene"] == gene), gene_row)
    variant_ids = set(hypothesis["variant_ids"]) if hypothesis else set()
    selected = [c for c in candidates if c["gene"] == gene and c["variant_id"] in variant_ids]
    with gzip.open(models, "rt", encoding="utf-8") as handle:
        metadata = json.loads(next(handle))
        reference_genes = [model for line in handle if (model := json.loads(line))["symbol"] == gene]
    context_dir = output / "context"
    write(context_dir / "candidates.json", selected)
    write(context_dir / "hypotheses.json", [hypothesis] if hypothesis else [])
    write(context_dir / "reference.json", {"assembly": metadata["assembly"], "genes": reference_genes})
    write(context_dir / "phenotype.json", read(baseline / "intake/phenotype.json"))
    spec = {key + "_path": str((context_dir / ("reference.json" if key == "reference" else key + ".json")).resolve()) for key in ("reference", "candidates", "hypotheses", "phenotype")}
    spec["selection_context"] = {"conditional_mechanism_gene": gene, "track1_top_gene": ranked[0]["gene"] if ranked else None, "mechanism_gene_rank": gene_row["rank"] if gene_row else None, "selected_genotype_model": hypothesis["type"] if hypothesis else None, "rationale": "Baseline-nominated mechanism branch; supported context is selected from the new genome-wide results, preferring its highest-ranked biallelic hypothesis. This does not establish causality."}
    write(context_dir / "spec.json", spec)
    return context_dir / "spec.json"


def comparison(output: Path, baseline: Path) -> dict:
    old_qc, new_qc = read(baseline / "track1/qc.json"), read(output / "track1/qc.json")
    old_hypotheses = read(baseline / "track1/hypotheses.json")
    new_genes = read(output / "track1/genes.json")
    new_biallelic = read(output / "track1/biallelic_genes.json")
    sensitivity = read(output / "track1/sensitivity.json")
    candidates = read(output / "track1/candidates.json.gz")
    by_allele = {(c["gene"], c["variant_id"]): c for c in candidates}
    old_drugs = read(baseline / "track2/track2_candidates.json")
    rescue = read(output / "track2/rescue_summary.json")
    rescue_rows = read(output / "track2/rescue_candidates.json")
    context = read(output / "track2/variant_context.json")
    lead = old_hypotheses[0]["gene"] if old_hypotheses else None
    lead_row = next((r for r in new_genes if r["gene"] == lead), None)
    lead_biallelic = next((r for r in new_biallelic if r["gene"] == lead), None)
    scenario_rows = [row for row in sensitivity if row["gene"] == lead]
    known = new_qc["population_af_known_candidates"]
    total = new_qc["candidate_count"]
    result = {"baseline": str(baseline.resolve()), "alternative": str(output.resolve()), "same_vcf": True,
              "records": new_qc["records_seen"], "baseline_gene_count": old_qc["panel_gene_count"],
              "alternative_gene_count": new_qc["reference_gene_count"], "ranked_gene_count": new_qc["ranked_gene_count"],
              "biallelic_ranked_gene_count": new_qc["biallelic_ranked_gene_count"], "baseline_lead_gene": lead,
              "alternative_lead_gene": new_genes[0]["gene"] if new_genes else None,
              "baseline_lead_new_rank": lead_row["rank"] if lead_row else None,
              "baseline_lead_biallelic_rank": lead_biallelic["rank"] if lead_biallelic else None,
              "lead_sensitivity": scenario_rows,
              "frequency_coverage": {"known_annotations": known, "total_annotations": total, "fraction": known / total if total else 0},
              "approved_lab_probes": rescue["approved_lab_shortlist_count"], "validated_rescues": 0,
              "top_lab_probe": next((r["name"] for r in rescue_rows if r["lab_priority_rank"] == 1), None),
              "interpretation": "Broader evidence and more specific experiments; no contest scoring or wet-lab validation has been performed."}
    write(output / "comparison.json", result)
    verification = []
    for hypothesis in old_hypotheses:
        for variant_id in hypothesis["variant_ids"]:
            candidate = by_allele.get((hypothesis["gene"], variant_id), {})
            verification.append({"baseline_gene": hypothesis["gene"], "variant": variant_id, "seen_in_alternative": bool(candidate), "protein_change": candidate.get("protein_change"), "population_af": candidate.get("population_af"), "population_af_source": candidate.get("population_af_source"), "clinvar": candidate.get("clinvar", []), "exclusion_reason": candidate.get("hypothesis_exclusion_reason")})
    write(output / "baseline_allele_comparison.json", verification)
    lines = ["# Comparing two MVA research strategies", "",
             "Strategy B replaces the small gene panel and regulator-enrichment starting point with genome-wide coding triage and allele-specific functional rescue. Both runs use the same supplied genome. This comparison assesses computational evidence and experimental priorities; it does not establish a diagnosis, drug benefit or higher contest score.", "",
             "## What changed", "", "| Dimension | Strategy A | Strategy B |", "| --- | --- | --- |",
             f"| Gene scope | {old_qc['panel_gene_count']} literature-selected genes | {new_qc['reference_gene_count']:,} eligible canonical coding genes |",
             "| Phenotypes | Retained for manual review | Explicit weighted HPO coverage; family history separate |",
             "| Clinical evidence | No local clinical database join | Exact local ClinVar matches with review/conflict labels |",
             f"| Population frequency | Unknown | Legacy database frequencies available for {known:,}/{total:,} annotations ({100*known/total:.1f}%) |",
             "| Inheritance | Recessive hypotheses | Separate biallelic and exploratory single-allele models; catalogued inheritance exposed |",
             "| Robustness | Not quantified | Remove phenotype, remove ClinVar and each individual patient feature |",
             "| Drug starting point | Downstream stress and regulator hypotheses | Stop codon, NMD prediction, amino-acid fidelity and residual protein function |",
             "| Drug decisions | Evidence gates | Approved-drug laboratory priorities with sequence/phenotype penalties and counter-evidence |", "",
             "## Track 1 comparison", "",
             f"The full alternative run processed **{new_qc['records_seen']:,} records**, retained **{total:,} exon/splice annotations**, and ranked **{new_qc['ranked_gene_count']:,} genes**. Reference QC exclusions are recorded; this covers canonical nuclear coding transcripts, not every possible WGS mechanism.", "",
             f"The baseline lead, **{lead}**, ranks **{result['baseline_lead_new_rank']}** in the alternative gene list and **{result['baseline_lead_biallelic_rank']}** in the separate biallelic list. The unrestricted alternative gene-list leader is **{result['alternative_lead_gene']}**.", "",
             "A high-ranking single allele can still be an incomplete recessive or carrier finding. The biallelic list is essential when evaluating a recessive explanation. The scores are not probabilities and do not combine independent likelihoods.", "",
             "### Leading genes", "", "| Rank | Gene | Genotype model | Inheritance context | Score |", "| --- | --- | --- | --- | --- |"]
    for row in new_genes[:10]:
        lines.append(f"| {row['rank']} | {row['gene']} | {row['type']} | {row['inheritance_compatibility']} | {row['score']:.3f} |")
    lines += ["", f"### Sensitivity of {lead}", "", "| Evidence scenario | Gene rank | Score |", "| --- | --- | --- |"]
    for row in scenario_rows:
        lines.append(f"| {row['scenario']} | {row['rank']} | {row['score']:.3f} |")
    lines += ["", "Rank changes identify which evidence supports prioritization. They are not independent validation, and the feature weights were not fitted to recover this case. Missing population-frequency data remains unknown. ClinVar annotations and phenotype catalogs have ascertainment bias, and HPO features can be correlated.", "",
              "[Full variant review](track1/review.md) · [Ranked genes](track1/genes.csv) · [Biallelic hypotheses](track1/biallelic_genes.csv) · [Sensitivity results](track1/sensitivity.csv) · [Recheck of baseline alleles](baseline_allele_comparison.json)", "",
              "## Track 2 comparison", "",
              f"Strategy A had {old_drugs['counts']['deferred']} deferred hypotheses/comparators and {old_drugs['counts']['rejected']} rejected candidate(s). Strategy B evaluates **{rescue['approved_lab_shortlist_count']} currently approved medicinal probes** for laboratory investigation and keeps other compounds outside that shortlist. Neither strategy has demonstrated rescue.", "",
              f"The first conditional laboratory priority is **{result['top_lab_probe']}**. This is an experiment-ordering decision, not a medication recommendation; point assignments are unvalidated heuristics.", ""]
    for stop in context.get("stops", []):
        if stop.get("verified_stop"):
            lines += [f"The local transcript calculation confirms **{stop['stop_codon_rna']}**, at exon **{stop['exon_number']}/{stop['exon_count']}**, with **{stop['bases_between_stop_end_and_final_junction']} nucleotides** from stop end to the final splice junction. This supports **{stop['nmd_prediction']}** as a positional prediction. Measure actual allele-specific RNA before assuming a drug has a substrate.", "",
                      f"The native amino acid is **{stop['reference_amino_acid']}**. Other readthrough systems reported **{', '.join(stop['reported_readthrough_amino_acids'])}** for this codon. Full-length protein therefore requires direct identity and checkpoint-function checks. This is a risk inferred from other systems, not a measured replacement residue in this allele.", ""]
    lines += ["Unknown phase is also consequential: readthrough of a cis double-mutant allele differs from readthrough of a stop allele in trans with a missense allele. Candidate approval status, prior evidence, unfavorable studies, normal-versus-tumor effects and exposure gaps are documented separately.", "",
              "[Allele-specific rescue review and source citations](track2/RESCUE_REVIEW.md) · [Drug comparison table](track2/rescue_candidates.csv)", "",
              "## What looks stronger for the entry", "",
              "Strategy B is the stronger computational foundation: it broadens coverage, can challenge the original hypothesis, exposes dependency on individual evidence sources, and yields allele-specific experiments that can disprove the drug rationale. That better addresses scientific rigor and reusability than a five-gene enrichment result alone. Strategy A remains a useful transparent comparator.", "",
              "The largest remaining gaps are current population-frequency coverage, independent transcript/variant confirmation, phase and read validation, alternate transcripts/CNV/SV analysis, and functional rescue at a plausible exposure. Drug-rank sensitivity must be retained in the final argument rather than hidden. A convincing contest proposal should lead with the discriminating experiment and its failure criteria, not an unsupported treatment claim.", "",
              "## Reproducibility and scope", "",
              "The baseline result was visible during design; this is not blinded replication. Strategy B uses public bulk reference downloads, matched locally. No patient data was uploaded. The baseline directory remains unchanged. The manifest fingerprints inputs, code at each execution stage, reference releases, outputs and baseline artifacts. No competition entry has been submitted.", "",
              "## Acknowledgement", "", ACKNOWLEDGEMENT, ""]
    (output / "COMPARISON.md").write_text("\n".join(lines), encoding="utf-8")
    return {key: value for key, value in result.items() if key != "lead_sensitivity"}


def run(baseline: Path, output: Path) -> dict:
    from .genomewide import run as genomewide_run, SCORING_POLICY
    if output.exists():
        raise ValueError("Output already exists; choose a fresh strategy B directory")
    baseline = baseline.resolve()
    original_baseline = baseline_fingerprints(baseline)
    config = read(Path("config/strategy_b.json"))
    refs = Path("references/alternative")
    model_path = refs / "genes_grch38.jsonl.gz"
    phenotype_path = Path("outputs/strategy_b_inputs/phenotype_reference.json")
    vcf = Path(read(baseline / "intake/inventory.json")["vcf"])
    input_paths = [vcf, model_path, phenotype_path, refs / "clinvar.vcf.gz", refs / "annotation.gtf.gz", refs / "cds.fa.gz", refs / "phenotype_to_genes.txt", refs / "phenotype.hpoa", Path("config/strategy_b.json")]
    inputs = fingerprints(input_paths)
    if inputs[str(vcf.resolve())] != read(baseline / "intake/inventory.json")["vcf_sha256"]:
        raise ValueError("VCF differs from the baseline input")
    output.mkdir(parents=True)
    write(output / "strategy_snapshot.json", config)
    manifest_path = output / "run_manifest.json"
    manifest = {"status": "running", "started_at_utc": datetime.now(timezone.utc).isoformat(), "python": platform.python_version(), "inputs_sha256": inputs, "baseline_sha256": original_baseline, "scoring_policy": SCORING_POLICY, "stages": {}}
    write(manifest_path, manifest)
    try:
        genome_code = [Path("mva") / name for name in ("genomewide.py", "track1.py", "vcf.py", "alternative_references.py", "strategy_b.py")]
        manifest["stages"]["track1"] = {"code_sha256": fingerprints(genome_code), "started_at_utc": datetime.now(timezone.utc).isoformat()}
        write(manifest_path, manifest)
        print("Running genome-wide coding triage and local clinical evidence matching...", flush=True)
        genome_summary = genomewide_run(vcf, model_path, phenotype_path, output / "track1", refs / "clinvar.vcf.gz")
        manifest["stages"]["track1"]["summary"] = genome_summary
        spec = create_rescue_context(output, baseline, model_path, config["mechanism_gene"])
        from .rescue import run as rescue_run
        rescue_code = [Path("mva/rescue.py"), Path("config/rescue_evidence.json")]
        manifest["stages"]["track2"] = {"code_sha256": fingerprints(rescue_code), "started_at_utc": datetime.now(timezone.utc).isoformat()}
        write(manifest_path, manifest)
        print("Running allele-specific rescue comparison...", flush=True)
        manifest["stages"]["track2"]["summary"] = rescue_run(Path("config/rescue_evidence.json"), output / "track2", spec)
        result = comparison(output, baseline)
        if baseline_fingerprints(baseline) != original_baseline:
            raise ValueError("Baseline changed during comparison")
        manifest.update(status="completed", completed_at_utc=datetime.now(timezone.utc).isoformat(), summary=result, baseline_unchanged=True)
        manifest["outputs_sha256"] = {p.relative_to(output).as_posix(): sha256(p) for p in sorted(output.rglob("*")) if p.is_file() and p != manifest_path}
    except Exception as error:
        manifest.update(status="failed", error=str(error))
        raise
    finally:
        write(manifest_path, manifest)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, default=Path("outputs/review_ready"))
    parser.add_argument("--output", type=Path, default=Path("outputs/strategy_b"))
    arguments = parser.parse_args()
    print(json.dumps(run(arguments.baseline, arguments.output), indent=2))
