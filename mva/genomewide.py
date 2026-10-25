"""Local, gene-panel-free canonical-transcript triage with explicit evidence scores.

Scores are frozen research priorities, not probabilities or pathogenicity calls.
Patient variants never leave the machine. Clinical databases are matched locally.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from .track1 import (GeneModel, IMPACT, QUALIFYING_CONSEQUENCES,
                     UNSUPPORTED_ANNOTATION_FLAGS, _candidate, _supplied_phase)
from .vcf import VCFReader, normalize_chrom, number, parse_info


SCORING_POLICY = {
    "genomic": "Mean of allele consequence (0-8), quality (-3/0/1), and source-labeled population AF (-3/0/1). Unknown AF contributes 0.",
    "inheritance": "Genotype-pattern support only: homozygous alternate 1; possible compound heterozygous .25, supplied trans 1; single heterozygous 0. Disease inheritance is unverified.",
    "phenotype": "6 times exact patient-HPO weighted coverage. Each HPO weight is log((number of annotated genes+1)/(genes with this HPO+1))+1, multiplied by any predefined feature weight.",
    "family_history": ".25 times independently calculated exact family-HPO weighted coverage. Family findings are never treated as patient findings.",
    "known_clinical": "Mean exact-allele ClinVar score: P/LP 4 at expert/guideline review, 3 at multiple submitters, 2 at one submitter; B/LB -6 with criteria; VUS/conflicting/unreviewed 0. Clinical disease context remains unverified.",
    "policy": "Scores are heuristic, frozen before the genome-wide run, not outcome-tuned, and not calibrated causal probabilities. Correlated HPO features are not independent evidence.",
}
LIMITATIONS = [
    "Canonical protein-coding transcript screen; alternative transcripts, deep intronic and regulatory variants are not exhaustively assessed.",
    "This comparison was designed with the earlier panel result visible; it is not a blinded recovery experiment.",
    "Exact HPO overlap is annotation-dependent; no ontology semantic similarity or phenotype independence is claimed.",
    "Unknown population frequency is not evidence of rarity. ClinVar frequency fields cover only a subset of alleles.",
    "Gene disease mechanism, inheritance and relevant ClinVar condition must be independently checked; single alleles do not establish recessive disease.",
    "Unphased pairs may be in cis. Same-block supplied cis pairs are excluded; supplied trans phase still requires validation.",
    "ClinVar exact REF/ALT joins do not normalize equivalent indels; a missing match does not imply benignity or novelty.",
    "Predictions use the supplied reference transcript and require independent annotation/read-level validation.",
    "No dedicated CNV, structural-variant, repeat, mitochondrial or mosaic analysis is performed.",
]


def _open_text(path: Path):
    path = Path(path)
    with path.open("rb") as stream:
        compressed = stream.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt", encoding="utf-8-sig") if compressed else path.open(encoding="utf-8-sig")


def load_reference(path: Path) -> tuple[dict, list[GeneModel]]:
    def model(raw):
        gene = GeneModel(raw)
        gene.reference_flags = list(raw.get("annotation_flags", []))
        return gene
    if str(path).endswith((".jsonl", ".jsonl.gz")):
        with _open_text(path) as stream:
            metadata = json.loads(next(stream))
            genes = [model(json.loads(line)) for line in stream if line.strip()]
    else:
        with _open_text(path) as stream:
            raw = json.load(stream)
        metadata = {key: value for key, value in raw.items() if key != "genes"}
        genes = [model(item) for item in raw["genes"]]
    if metadata.get("assembly") != "GRCh38":
        raise ValueError("Genome-wide triage requires a GRCh38 reference")
    if not genes:
        raise ValueError("Reference has no genes")
    symbols = [gene.symbol for gene in genes]
    if len(set(symbols)) != len(symbols):
        raise ValueError("Reference requires exactly one selected transcript per unique gene symbol")
    return metadata, genes


class ExonIndex:
    """Fixed-width exon bins followed by precise overlap; no per-record gene scan."""
    def __init__(self, genes: list[GeneModel], bin_size: int = 16384):
        self.bin_size = bin_size
        self.bins: dict[str, dict[int, list[tuple[int, int, int]]]] = defaultdict(lambda: defaultdict(list))
        for index, gene in enumerate(genes):
            for start, end in gene.exons:
                # Include both splice flanks; annotate() determines the actual effect.
                start, end = max(1, start - 8), end + 8
                for bin_id in range(start // bin_size, end // bin_size + 1):
                    self.bins[gene.chrom][bin_id].append((start, end, index))

    def query(self, chrom: str, start: int, end: int) -> list[int]:
        bins = self.bins.get(normalize_chrom(chrom), {})
        first, last = start // self.bin_size, end // self.bin_size
        keys = range(first, last + 1) if last - first < 10000 else (key for key in bins if first <= key <= last)
        return sorted({index for key in keys for low, high, index in bins.get(key, ()) if start <= high and end >= low})


class PhenotypeIndex:
    def __init__(self, raw: dict):
        self.genes = {symbol: set(values) for symbol, values in raw.get("genes", {}).items()}
        self.patient = raw.get("patient_features", [row for row in raw.get("features", []) if row.get("subject") != "family_history"])
        self.family = raw.get("family_features", [row for row in raw.get("features", []) if row.get("subject") == "family_history"])
        self.patient = self._unique(self.patient)
        self.family = self._unique(self.family)
        counts = Counter(term for terms in self.genes.values() for term in terms)
        total = len(self.genes)
        self.weights = {}
        for row in self.patient + self.family:
            multiplier = float(row.get("weight", 1))
            if not math.isfinite(multiplier) or multiplier <= 0:
                raise ValueError("Phenotype feature weights must be finite and positive")
            self.weights[row["hpo_id"]] = (math.log((total + 1) / (counts[row["hpo_id"]] + 1)) + 1) * multiplier
        self.source = raw.get("source", raw.get("sources", {}))
        self.inheritance = raw.get("gene_inheritance", {})
        self.inheritance_sources = raw.get("gene_inheritance_sources", {})

    def inheritance_modes(self, gene: str) -> set[str]:
        value = self.inheritance.get(gene, [])
        if isinstance(value, dict):
            value = value.get("modes", [])
        return set(value)

    @staticmethod
    def _unique(rows: list[dict]) -> list[dict]:
        return list({row["hpo_id"]: row for row in rows if row.get("hpo_id")}.values())

    def coverage(self, gene: str, family: bool = False, omitted: str | None = None) -> tuple[float, list[str]]:
        features = [row["hpo_id"] for row in (self.family if family else self.patient) if row["hpo_id"] != omitted]
        denominator = sum(self.weights[term] for term in features)
        hits = sorted(set(features) & self.genes.get(gene, set()))
        return (sum(self.weights[term] for term in hits) / denominator if denominator else 0), hits


def _clinical_score(evidence: list[dict]) -> float:
    if not evidence:
        return 0
    signatures = {row.get("clinical_significance", "").lower() for row in evidence}
    if any("conflict" in sig for sig in signatures) or any(row.get("conflicting_assertions") not in (None, "", ".") for row in evidence):
        return 0
    pathogenic = any("pathogenic" in sig and "conflict" not in sig for sig in signatures)
    benign = any("benign" in sig for sig in signatures)
    if pathogenic and benign:
        return 0
    scores = []
    for row in evidence:
        if not row.get("allele_assignment_unambiguous", False):
            continue
        review = row.get("review_status", "").lower()
        sig = row.get("clinical_significance", "").lower()
        supported = review.startswith("criteria_provided") or any(tag in review for tag in ("reviewed_by_expert_panel", "practice_guideline"))
        if not supported:
            continue
        if "benign" in sig:
            scores.append(-6)
        elif "pathogenic" in sig:
            scores.append(4 if "expert_panel" in review or "practice_guideline" in review else 3 if "multiple_submitters" in review else 2)
    return min(scores) if scores and max(scores) <= 0 else max(scores, default=0)


def join_clinvar(path: Path, candidates: list[dict]) -> dict:
    """Stream the public VCF and retain exact allele matches; never patient queries."""
    wanted = {row["variant_id"] for row in candidates}
    matched: dict[str, list[dict]] = defaultdict(list)
    frequencies = {}
    headers, rows, ambiguous = [], 0, 0
    with _open_text(path) as stream:
        for line in stream:
            if line.startswith("#"):
                if line.startswith(("##fileDate", "##reference", "##source", "##assembly")):
                    headers.append(line.strip())
                continue
            fields = line.rstrip("\r\n").split("\t", 8)
            if len(fields) < 8:
                continue
            rows += 1
            chrom, pos, accession, ref, alts, _, _, info_raw = fields[:8]
            for alt in alts.split(","):
                variant_id = f"{normalize_chrom(chrom)}:{pos}:{ref.upper()}:{alt.upper()}"
                if variant_id not in wanted:
                    continue
                info = parse_info(info_raw)
                unambiguous = "," not in alts and str(info.get("CLNALLE", "1")) in {"1", "1.0"}
                if not unambiguous:
                    ambiguous += 1
                matched[variant_id].append({
                    "variation_id": accession, "clinical_significance": str(info.get("CLNSIG", "")),
                    "review_status": str(info.get("CLNREVSTAT", "")), "condition": str(info.get("CLNDN", "")),
                    "condition_ids": str(info.get("CLNDISDB", "")), "conflicting_assertions": str(info.get("CLNSIGCONF", "")),
                    "allele_assignment_unambiguous": unambiguous, "match": "exact_chrom_pos_ref_alt_GRCh38",
                })
                if unambiguous:
                    observed = []
                    for tag in ("AF_EXAC", "AF_TGP", "AF_ESP"):
                        value = number(str(info.get(tag, ".")))
                        if value is not None and 0 <= value <= 1:
                            observed.append((float(value), tag))
                    if observed:
                        af, tag = max(observed)
                        frequencies[variant_id] = {"af": af, "source": f"ClinVar INFO/{tag} public database field; max of available ExAC/TGP/ESP fields"}
    if any("GRCh37" in header or "hg19" in header for header in headers):
        raise ValueError("ClinVar header indicates an incompatible GRCh37 assembly")
    for candidate in candidates:
        candidate["population_af_without_clinvar"] = candidate["population_af"]
        candidate["population_af_source_without_clinvar"] = candidate["population_af_source"]
        candidate["genomic_score_without_clinvar"] = candidate["triage_score"]
        candidate["frequency_score_without_clinvar"] = candidate["score_components"]["external_population_frequency"]
        candidate["clinvar"] = matched.get(candidate["variant_id"], [])
        candidate["clinical_score"] = _clinical_score(candidate["clinvar"])
        if candidate["population_af"] is None and candidate["variant_id"] in frequencies:
            population = frequencies[candidate["variant_id"]]
            candidate["population_af"], candidate["population_af_source"] = population["af"], population["source"]
            candidate["score_components"]["external_population_frequency"] = 1 if population["af"] < .001 else -3 if population["af"] > .01 else 0
            candidate["triage_score"] = sum(candidate["score_components"].values())
    return {"path": str(Path(path).resolve()), "records_scanned": rows, "matched_unique_alleles": len(matched),
            "ambiguous_matching_records": ambiguous, "population_frequency_matches": len(frequencies), "headers": headers}


def remove_clinvar_evidence(candidates: list[dict]) -> list[dict]:
    """Reset scores AND eligibility inputs before the database-removal ablation."""
    stripped = []
    for original in candidates:
        row = dict(original)
        row["clinvar"], row["clinical_score"], row["supplied_cis_partner_ids"] = [], 0, []
        row["population_af"] = row.get("population_af_without_clinvar", row["population_af"])
        row["population_af_source"] = row.get("population_af_source_without_clinvar", row["population_af_source"])
        row["triage_score"] = row.get("genomic_score_without_clinvar", row["triage_score"])
        row["score_components"] = dict(row["score_components"])
        row["score_components"]["external_population_frequency"] = row.get("frequency_score_without_clinvar", row["score_components"]["external_population_frequency"])
        stripped.append(row)
    return stripped


def _eligibility(candidate: dict) -> str | None:
    if candidate["quality_status"] == "fail":
        return "quality_failed"
    if UNSUPPORTED_ANNOTATION_FLAGS.intersection(candidate["annotation_flags"]):
        return "unsupported_consequence_annotation"
    if "selenocysteine_requires_specialist_annotation" in candidate["annotation_flags"]:
        return "selenocysteine_requires_specialist_annotation"
    if candidate["population_af"] is not None and candidate["population_af"] > .01:
        return "population_af_above_0.01"
    if candidate["consequence"] not in QUALIFYING_CONSEQUENCES and candidate.get("clinical_score", 0) <= 0:
        return "not_qualifying_consequence_or_reviewed_clinical_evidence"
    if candidate["zygosity"] not in {"heterozygous", "homozygous_alternate"}:
        return "unsupported_ploidy_or_genotype"
    return None


def build_hypotheses(candidates: list[dict], phenotype: PhenotypeIndex) -> tuple[list[dict], dict]:
    groups = defaultdict(list)
    exclusions, generated = Counter(), Counter()
    seen = set()
    for candidate in candidates:
        key = candidate["gene"], candidate["variant_id"]
        if key in seen:
            continue
        seen.add(key)
        reason = _eligibility(candidate)
        candidate["hypothesis_exclusion_reason"] = reason
        if reason:
            exclusions[reason] += 1
        else:
            groups[candidate["gene"]].append(candidate)
    retained = []
    cis_exclusions = 0
    for gene, variants in groups.items():
        coverage, hits = phenotype.coverage(gene)
        family_coverage, family_hits = phenotype.coverage(gene, family=True)
        buckets = defaultdict(list)

        def add(alleles: list[dict], model: str, phase: str, inheritance_score: float):
            alleles = sorted(alleles, key=lambda item: item["variant_id"])
            generated[model] += 1
            components = {"genomic": sum(row["triage_score"] for row in alleles) / len(alleles),
                          "inheritance": inheritance_score, "phenotype": 6 * coverage,
                          "family_history": .25 * family_coverage,
                          "known_clinical": sum(row.get("clinical_score", 0) for row in alleles) / len(alleles)}
            flags = ["disease_inheritance_and_mechanism_unverified", "requires_read_level_and_independent_annotation_validation"]
            if any(row["population_af"] is None for row in alleles):
                flags.append("population_frequency_unknown")
            if any(row.get("clinical_score", 0) < 0 for row in alleles):
                flags.append("contains_clinvar_benign_or_likely_benign_allele")
            if model == "possible_compound_heterozygous":
                flags.append("caller_phase_requires_validation" if phase == "supplied_in_trans" else "in_trans_not_established")
            if model == "single_heterozygous_exploratory":
                flags += ["dominant_mechanism_unverified", "single_allele_does_not_establish_recessive_disease"]
            modes = phenotype.inheritance_modes(gene)
            ar_only = bool(modes) and modes <= {"autosomal_recessive", "AR", "HP:0000007"}
            compatibility = "disease_inheritance_unverified"
            if ar_only and model == "single_heterozygous_exploratory":
                flags.append("incomplete_recessive_or_carrier_finding_in_current_disease_catalog")
                compatibility = "incomplete_known_recessive"
            elif ar_only:
                compatibility = "biallelic_pattern_compatible_with_catalogued_recessive_inheritance"
            row = {"gene": gene, "type": model, "variant_ids": sorted(item["variant_id"] for item in alleles),
                   "consequences": [item["consequence"] for item in alleles], "protein_changes": [item["protein_change"] for item in alleles],
                   "phase": phase, "phase_unconfirmed": model == "possible_compound_heterozygous",
                   "population_frequency_status": "unknown" if any(item["population_af"] is None for item in alleles) else "source_annotated",
                   "score_components": components, "score": sum(components.values()),
                   "matched_patient_hpo": hits, "matched_family_hpo": family_hits,
                   "patient_weighted_coverage": coverage, "family_weighted_coverage": family_coverage,
                   "catalogued_gene_inheritance": sorted(modes), "inheritance_compatibility": compatibility,
                   "status": "unvalidated_research_hypothesis", "review_flags": flags}
            # Preserve winners for both clinical-evidence scenarios. The phenotype
            # score is shared within a gene, so every phenotype ablation is exact.
            for mode in ("full", "no_clinvar"):
                bucket = buckets[(model, mode)]
                bucket.append(row)
                if len(bucket) > 6:
                    bucket.sort(key=lambda item: (-(item["score"] - (item["score_components"]["known_clinical"] if mode == "no_clinvar" else 0)), item["variant_ids"]))
                    del bucket[3:]

        heterozygotes = []
        for variant in variants:
            if variant["zygosity"] == "homozygous_alternate":
                add([variant], "homozygous_alternate", "not_applicable", 1)
            else:
                heterozygotes.append(variant)
                add([variant], "single_heterozygous_exploratory", "not_applicable", 0)
        for first, second in itertools.combinations(heterozygotes, 2):
            phase = _supplied_phase(first, second)
            if phase == "supplied_in_cis":
                cis_exclusions += 1
                first["supplied_cis_partner_ids"].append(second["variant_id"])
                second["supplied_cis_partner_ids"].append(first["variant_id"])
                continue
            add([first, second], "possible_compound_heterozygous", phase, 1 if phase == "supplied_in_trans" else .25)
        unique = {}
        for (model, mode), bucket in buckets.items():
            bucket.sort(key=lambda item: (-(item["score"] - (item["score_components"]["known_clinical"] if mode == "no_clinvar" else 0)), item["variant_ids"]))
            for row in bucket[:3]:
                unique[(model, tuple(row["variant_ids"]))] = row
        retained.extend(unique.values())
    retained.sort(key=lambda row: (-row["score"], row["gene"], row["type"], row["variant_ids"]))
    return retained, {"allele_exclusions": dict(exclusions), "generated_hypotheses_by_model": dict(generated),
                      "generated_hypothesis_count": sum(generated.values()), "retained_hypothesis_count": len(retained),
                      "supplied_cis_pair_exclusions": cis_exclusions,
                      "retention_policy": "Top 3 per gene/model under full and no-ClinVar scores; all eligible pairs evaluated. No fixed gene filter."}


def rank_genes(hypotheses: list[dict], phenotype: PhenotypeIndex, scenario: str = "full") -> list[dict]:
    best = {}
    for hypothesis in hypotheses:
        components = dict(hypothesis["score_components"])
        if scenario == "no_clinvar":
            components["known_clinical"] = 0
        elif scenario == "no_phenotype":
            components["phenotype"] = components["family_history"] = 0
        elif scenario.startswith("omit:"):
            components["phenotype"] = 6 * phenotype.coverage(hypothesis["gene"], omitted=scenario[5:])[0]
        row = {**hypothesis, "score_components": components, "score": sum(components.values()), "scenario": scenario}
        prior = best.get(row["gene"])
        def sort_key(item):
            return (item.get("inheritance_compatibility") == "incomplete_known_recessive", -item["score"], item["type"], item["variant_ids"])
        if prior is None or sort_key(row) < sort_key(prior):
            best[row["gene"]] = row
    ranked = sorted(best.values(), key=lambda row: (-row["score"], row["gene"]))
    previous_score, score_rank = None, 0
    for index, row in enumerate(ranked, 1):
        if previous_score is None or abs(row["score"] - previous_score) > 1e-10:
            score_rank = index
        row["rank"], row["score_rank"] = index, score_rank
        previous_score = row["score"]
    return ranked


def _write_csv(path: Path, rows: list[dict], columns: list[str]):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(value, ensure_ascii=False) if isinstance(value, (dict, list)) else value for key, value in row.items()})


def run(vcf_path: Path, reference_path: Path, phenotype_path: Path, output_dir: Path,
        clinvar_path: Path | None = None, population_path: Path | None = None, max_records: int | None = None) -> dict:
    vcf_path, reference_path, phenotype_path, output_dir = map(Path, (vcf_path, reference_path, phenotype_path, output_dir))
    metadata, genes = load_reference(reference_path)
    with _open_text(phenotype_path) as stream:
        phenotype = PhenotypeIndex(json.load(stream))
    frequencies = {}
    if population_path:
        frequencies = json.loads(Path(population_path).read_text(encoding="utf-8-sig"))
        frequencies = frequencies.get("population_frequencies", frequencies)
    index = ExonIndex(genes)
    reader = VCFReader(vcf_path)
    candidates, counts, chrom_counts = [], Counter(), Counter()
    for record in reader.records(max_records):
        if reader.assembly in {"GRCh37", "conflicting"}:
            raise ValueError(f"VCF/reference assembly mismatch: {reader.assembly}")
        counts["parsed_records"] += 1
        chrom_counts[record.chrom] += 1
        if not record.called_alt_indices:
            counts["no_called_alternate"] += 1
            continue
        overlaps = index.query(record.chrom, record.pos, record.end + 1)
        if overlaps:
            counts["exon_or_splice_window_overlapping_records"] += 1
        for gene_index in overlaps:
            for alt_index in record.called_alt_indices:
                candidate = _candidate(record, alt_index, genes[gene_index], frequencies)
                candidate["annotation_flags"].extend(genes[gene_index].reference_flags)
                if candidate["consequence"] == "intron_variant" and not candidate["additional_consequences"]:
                    continue
                candidate["clinvar"], candidate["clinical_score"] = [], 0
                candidates.append(candidate)
        if reader.records_seen % 1000000 == 0:
            print(f"Scanned {reader.records_seen:,} VCF records; retained {len(candidates):,} exon/splice annotations", flush=True)
    if reader.assembly in {"GRCh37", "conflicting"}:
        raise ValueError(f"VCF/reference assembly mismatch: {reader.assembly}")
    print(f"Genome scan finished: {reader.records_seen:,} records, {len(candidates):,} annotations", flush=True)
    clinical_qc = join_clinvar(Path(clinvar_path), candidates) if clinvar_path else {"status": "not_supplied"}
    hypotheses, hypothesis_qc = build_hypotheses(candidates, phenotype)
    no_clinvar_candidates = remove_clinvar_evidence(candidates)
    no_clinvar_hypotheses, no_clinvar_qc = build_hypotheses(no_clinvar_candidates, phenotype)
    del no_clinvar_candidates
    ranked = rank_genes(hypotheses, phenotype)
    biallelic_ranked = rank_genes([row for row in hypotheses if row["type"] != "single_heterozygous_exploratory"], phenotype)
    biallelic_by_gene = {row["gene"]: row for row in biallelic_ranked}
    for row in ranked:
        biallelic = biallelic_by_gene.get(row["gene"])
        row["best_biallelic_hypothesis"] = {key: biallelic[key] for key in ("rank", "type", "variant_ids", "score", "phase")} if biallelic else None
    scenarios = ["full", "no_phenotype", "no_clinvar"] + ["omit:" + row["hpo_id"] for row in phenotype.patient]
    sensitivity = []
    base_ranks = {row["gene"]: row["rank"] for row in ranked}
    for scenario in scenarios:
        scenario_hypotheses = no_clinvar_hypotheses if scenario == "no_clinvar" else hypotheses
        for row in rank_genes(scenario_hypotheses, phenotype, scenario):
            base_rank = base_ranks.get(row["gene"])
            sensitivity.append({"scenario": scenario, "gene": row["gene"], "rank": row["rank"], "score_rank": row["score_rank"],
                                "base_rank": base_rank, "rank_change": base_rank - row["rank"] if base_rank is not None else None,
                                "score": row["score"], "type": row["type"], "variant_ids": row["variant_ids"], "score_components": row["score_components"]})
    limitations = list(LIMITATIONS)
    if reader.assembly == "unknown":
        limitations.append("VCF assembly is unverified; coordinates provisionally interpreted as GRCh38.")
    if not reader.complete:
        limitations.append("The VCF scan is a prefix only because max_records was supplied.")
    qc = {"analysis": "canonical_transcript_genome_wide_genotype_first_triage", "created_at_utc": datetime.now(timezone.utc).isoformat(),
          "input_vcf": str(vcf_path.resolve()), "reference_path": str(reference_path.resolve()), "phenotype_path": str(phenotype_path.resolve()),
          "assembly": reader.assembly, "assembly_evidence": reader.assembly_evidence, "scan_complete": reader.complete,
          "records_seen": reader.records_seen, "malformed_records": reader.malformed_records, "malformed_examples": reader.malformed_examples,
          "counts": dict(counts), "chromosome_counts": dict(chrom_counts), "reference_metadata": metadata,
          "reference_gene_count": len(genes), "genes_with_candidate_annotations": len({row["gene"] for row in candidates}),
          "invalid_cds_models": sum(not gene.sequence_valid for gene in genes), "candidate_count": len(candidates),
          "unique_candidate_alleles": len({row["variant_id"] for row in candidates}), "ranked_gene_count": len(ranked),
          "consequence_counts": dict(Counter(row["consequence"] for row in candidates)),
          "candidate_quality_counts": dict(Counter(row["quality_status"] for row in candidates)),
          "population_af_known_candidates": sum(row["population_af"] is not None for row in candidates),
          "patient_features": phenotype.patient, "family_features": phenotype.family,
          "phenotype_weights": phenotype.weights, "phenotype_reference_genes": len(phenotype.genes), "phenotype_source": phenotype.source,
          "gene_inheritance_sources": phenotype.inheritance_sources, "biallelic_ranked_gene_count": len(biallelic_ranked),
          "clinvar": clinical_qc, "no_clinvar_ablation": {**no_clinvar_qc, "policy": "All ClinVar significance and ClinVar-supplied population AF removed; eligibility and pairs rebuilt."},
          **hypothesis_qc, "scoring_policy": SCORING_POLICY, "limitations": limitations}
    output_dir.mkdir(parents=True, exist_ok=True)
    for filename, value in (("hypotheses.json", hypotheses), ("no_clinvar_hypotheses.json", no_clinvar_hypotheses), ("genes.json", ranked), ("biallelic_genes.json", biallelic_ranked), ("qc.json", qc), ("sensitivity.json", sensitivity)):
        (output_dir / filename).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    with gzip.open(output_dir / "candidates.json.gz", "wt", encoding="utf-8") as stream:
        json.dump(candidates, stream, allow_nan=False)
    _write_csv(output_dir / "candidates.csv", candidates, ["variant_id", "gene", "transcript_id", "consequence", "additional_consequences", "protein_change", "genotype", "zygosity", "quality_status", "quality_flags", "population_af", "population_af_source", "phase_set", "supplied_cis_partner_ids", "triage_score", "clinical_score", "clinvar", "annotation_flags", "hypothesis_exclusion_reason"])
    gene_columns = ["rank", "score_rank", "gene", "score", "type", "variant_ids", "protein_changes", "phase", "population_frequency_status", "matched_patient_hpo", "matched_family_hpo", "score_components", "catalogued_gene_inheritance", "inheritance_compatibility", "best_biallelic_hypothesis", "review_flags"]
    _write_csv(output_dir / "genes.csv", ranked, gene_columns)
    _write_csv(output_dir / "biallelic_genes.csv", biallelic_ranked, gene_columns)
    _write_csv(output_dir / "incomplete_recessive.csv", [row for row in ranked if row["inheritance_compatibility"] == "incomplete_known_recessive"], gene_columns)
    _write_csv(output_dir / "sensitivity.csv", sensitivity, ["scenario", "gene", "rank", "score_rank", "base_rank", "rank_change", "score", "type", "variant_ids", "score_components"])
    _write_review(output_dir / "review.md", qc, ranked, sensitivity, biallelic_ranked)
    return {"track": 1, "strategy": "B", "output_dir": str(output_dir.resolve()), "scan_complete": reader.complete,
            "records_seen": reader.records_seen, "candidate_count": len(candidates), "ranked_gene_count": len(ranked),
            "hypothesis_count": hypothesis_qc["generated_hypothesis_count"], "top_genes": [row["gene"] for row in ranked[:20]]}


def _write_review(path: Path, qc: dict, ranked: list[dict], sensitivity: list[dict], biallelic: list[dict]):
    lines = ["# Strategy B: genome-wide genotype-first comparison", "", "Exploratory research; scores prioritize review and are not diagnoses or probabilities.", "",
             f"Scanned {qc['records_seen']:,} records (complete: {qc['scan_complete']}) against {qc['reference_gene_count']:,} selected transcripts. Retained {qc['candidate_count']:,} exon/splice annotations and ranked {qc['ranked_gene_count']:,} genes.", "",
             "No disease-gene panel determines inclusion. The earlier result was already visible when this comparison was designed.", "",
             "## Top genes under the frozen full score", "", "| Rank | Gene | Model | Genomic | Genotype | Patient HPO | Family | Clinical | Total |", "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"]
    for row in ranked[:30]:
        c = row["score_components"]
        lines.append(f"| {row['rank']} | {row['gene']} | {row['type']} | {c['genomic']:.2f} | {c['inheritance']:.2f} | {c['phenotype']:.2f} | {c['family_history']:.2f} | {c['known_clinical']:.2f} | {row['score']:.2f} |")
    lines += ["", "## Biallelic hypotheses ranked separately", "", "Single-allele scores cannot establish recessive disease. This separate list prevents a high-impact single allele from hiding a lower-mean-score pair. Known recessive-only gene annotations favor the biallelic model within that gene, but disease-specific inheritance remains unverified.", "", "| Rank | Gene | Model | Variants | Phase | Score |", "| --- | --- | --- | --- | --- | --- |"]
    for row in biallelic[:20]:
        lines.append(f"| {row['rank']} | {row['gene']} | {row['type']} | {', '.join(row['variant_ids'])} | {row['phase']} | {row['score']:.2f} |")
    lines += ["", "## Sensitivity: leading genes by scenario", "", "| Scenario | First | Second | Third |", "| --- | --- | --- | --- |"]
    scenarios = list(dict.fromkeys(row["scenario"] for row in sensitivity))
    for scenario in scenarios:
        leaders = [row["gene"] for row in sensitivity if row["scenario"] == scenario][:3]
        lines.append("| " + " | ".join([scenario, *(leaders + [""] * (3 - len(leaders)))]) + " |")
    lines += ["", "Ranks use gene-symbol order to break equal scores; score_rank in the CSV exposes ties. A changed leader under omitted evidence indicates dependence on that evidence, not proof against causality.", "", "## Scoring policy", ""]
    lines.extend(f"- **{key}:** {value}" for key, value in SCORING_POLICY.items())
    lines += ["", "## Limits and validation gates", ""]
    lines.extend("- " + item for item in qc["limitations"])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--phenotype", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--clinvar", type=Path)
    parser.add_argument("--population", type=Path)
    parser.add_argument("--max-records", type=int)
    args = parser.parse_args()
    print(json.dumps(run(args.vcf, args.reference, args.phenotype, args.output, args.clinvar, args.population, args.max_records), indent=2))
