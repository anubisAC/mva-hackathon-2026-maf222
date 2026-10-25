"""Hypothesis-led local variant triage; this is not a clinical diagnosis.

The complete VCF is audited, while consequence annotation is intentionally
limited to the supplied gene/transcript reference. Patient variants stay local.
"""
from __future__ import annotations

import csv
import itertools
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from .vcf import VCFReader, VCFRecord, normalize_chrom, number, parse_info


_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
GENETIC_CODE = {a + b + c: aa for (a, b, c), aa in zip(itertools.product(_BASES, repeat=3), _AAS)}
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
IMPACT = {
    "stop_gained": 8, "frameshift_variant": 8, "splice_donor_variant": 7,
    "splice_acceptor_variant": 7, "start_lost": 6, "stop_lost": 6,
    "inframe_insertion": 5, "inframe_deletion": 5, "missense_variant": 4,
    "coding_sequence_variant": 3, "coding_boundary_variant": 3,
    "structural_variant": 3, "splice_region_variant": 2,
    "synonymous_variant": 0, "stop_retained_variant": 0,
    "5_prime_UTR_variant": 0, "3_prime_UTR_variant": 0,
    "non_coding_transcript_exon_variant": 0, "intron_variant": 0,
}
QUALIFYING_CONSEQUENCES = {
    "stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant",
    "start_lost", "stop_lost", "inframe_insertion", "inframe_deletion", "missense_variant",
}
UNSUPPORTED_ANNOTATION_FLAGS = {
    "reference_base_mismatch_or_ambiguous_codon", "reference_sequence_mismatch",
    "cds_sequence_missing_or_length_mismatch", "transcript_model_missing",
}
LIMITATIONS = [
    "Gene selection is hypothesis-led; this is not genome-wide causal-variant discovery.",
    "Only supplied transcript models are annotated. Alternative transcripts and deep intronic effects are not evaluated.",
    "Consequence and triage scores are heuristic evidence, not pathogenicity classifications, diagnoses, or calibrated EPCR probabilities.",
    "Population allele frequency is unknown unless a source-labeled external frequency is supplied; VCF INFO/AF is never used as population AF.",
    "Same-gene heterozygous pairs remain unconfirmed phase hypotheses; supplied same-block GT/PS cis pairs are excluded and supplied trans phase still requires validation.",
    "No dedicated CNV, repeat expansion, mitochondrial, mosaicism, or structural-variant discovery is performed by this small-variant baseline.",
    "Local consequences require independent clinical annotation and validation; no complex HGVS or protein effect is invented.",
]


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


class GeneModel:
    def __init__(self, raw: dict):
        self.symbol = raw["symbol"]
        self.chrom = normalize_chrom(str(raw["chrom"]))
        self.start, self.end = int(raw["start"]), int(raw["end"])
        self.strand = int(raw["strand"])
        if self.strand not in {-1, 1} or self.start < 1 or self.end < self.start:
            raise ValueError(f"Invalid gene coordinates for {self.symbol}")
        self.transcript = raw.get("transcript_id")
        self.exons = sorted((int(exon["start"]), int(exon["end"])) for exon in raw.get("exons", []))
        if any(start < self.start or end > self.end or end < start for start, end in self.exons):
            raise ValueError(f"Invalid exon coordinates for {self.symbol}")
        self.coding_start = raw.get("translation_start")
        self.coding_end = raw.get("translation_end")
        self.cds = raw.get("cds_sequence", "").upper()
        self.segments = []
        if self.coding_start is not None and self.coding_end is not None:
            self.coding_start, self.coding_end = sorted((int(self.coding_start), int(self.coding_end)))
            self.segments = [(max(start, self.coding_start), min(end, self.coding_end))
                             for start, end in self.exons
                             if start <= self.coding_end and end >= self.coding_start]
        if self.strand == -1:
            self.segments.reverse()
        self.cds_length = sum(end - start + 1 for start, end in self.segments)
        self.sequence_valid = bool(self.cds) and len(self.cds) == self.cds_length

    def cds_offset(self, position: int) -> int | None:
        preceding = 0
        for start, end in self.segments:
            if start <= position <= end:
                return preceding + (position - start if self.strand == 1 else end - position)
            preceding += end - start + 1
        return None

    def splice_effect(self, start: int, end: int) -> str | None:
        # Only internal exon boundaries define splice junctions.
        for (_, left_end), (right_start, _) in zip(self.exons, self.exons[1:]):
            if start <= left_end + 2 and end >= left_end + 1:
                return "splice_donor_variant" if self.strand == 1 else "splice_acceptor_variant"
            if start <= right_start - 1 and end >= right_start - 2:
                return "splice_acceptor_variant" if self.strand == 1 else "splice_donor_variant"
        for (_, left_end), (right_start, _) in zip(self.exons, self.exons[1:]):
            if start <= left_end + 8 and end >= left_end - 2:
                return "splice_region_variant"
            if start <= right_start + 2 and end >= right_start - 8:
                return "splice_region_variant"
        return None

    def annotate(self, record: VCFRecord, alt: str) -> dict:
        result = {"gene": self.symbol, "transcript_id": self.transcript,
                  "consequence": "intron_variant", "additional_consequences": [],
                  "cds_position": None, "protein_change": None, "reference_codon": None,
                  "alternate_codon": None, "annotation_flags": []}
        if not self.transcript or not self.exons:
            result["annotation_flags"].append("transcript_model_missing")
        if alt.startswith("<") or alt == "*" or "[" in alt or "]" in alt:
            result["consequence"] = "structural_variant"
            result["annotation_flags"].append("symbolic_or_breakend_allele_requires_specialist_annotation")
            return result
        # Remove VCF anchor bases before classifying the affected interval.
        # Insertions span a breakpoint, so either flank can touch a splice site.
        changed_ref, changed_alt = record.ref, alt
        change_start = record.pos
        if len(record.ref) != len(alt):
            while changed_ref and changed_alt and changed_ref[0] == changed_alt[0]:
                changed_ref, changed_alt = changed_ref[1:], changed_alt[1:]
                change_start += 1
            while changed_ref and changed_alt and changed_ref[-1] == changed_alt[-1]:
                changed_ref, changed_alt = changed_ref[:-1], changed_alt[:-1]
        is_insertion = not changed_ref
        start = change_start - 1 if is_insertion else change_start
        end = change_start if is_insertion else change_start + len(changed_ref) - 1
        splice = self.splice_effect(start, end)
        offset = self.cds_offset(start)
        coding_overlap = any(start <= seg_end and end >= seg_start for seg_start, seg_end in self.segments)
        if coding_overlap:
            if self.sequence_valid:
                for ref_position, ref_base in enumerate(record.ref, record.pos):
                    ref_offset = self.cds_offset(ref_position)
                    if ref_offset is not None:
                        expected_base = ref_base if self.strand == 1 else reverse_complement(ref_base)
                        if self.cds[ref_offset] != expected_base:
                            result["consequence"] = "coding_sequence_variant"
                            result["annotation_flags"].append("reference_sequence_mismatch")
                            return result
            else:
                result["consequence"] = "coding_sequence_variant"
                result["annotation_flags"].append("cds_sequence_missing_or_length_mismatch")
                return result
            if len(record.ref) == len(alt) == 1 and offset is not None:
                result["cds_position"] = offset + 1
                if not self.sequence_valid:
                    result["consequence"] = "coding_sequence_variant"
                    result["annotation_flags"].append("cds_sequence_missing_or_length_mismatch")
                else:
                    ref_base = record.ref if self.strand == 1 else reverse_complement(record.ref)
                    alt_base = alt if self.strand == 1 else reverse_complement(alt)
                    codon_start = (offset // 3) * 3
                    codon = self.cds[codon_start:codon_start + 3]
                    mutated = codon[:offset % 3] + alt_base + codon[offset % 3 + 1:]
                    if self.cds[offset] != ref_base or codon not in GENETIC_CODE or mutated not in GENETIC_CODE:
                        result["consequence"] = "coding_sequence_variant"
                        result["annotation_flags"].append("reference_base_mismatch_or_ambiguous_codon")
                    else:
                        old_aa, new_aa = GENETIC_CODE[codon], GENETIC_CODE[mutated]
                        result["reference_codon"], result["alternate_codon"] = codon, mutated
                        result["protein_change"] = f"p.{old_aa}{offset // 3 + 1}{new_aa}"
                        if new_aa == old_aa:
                            result["consequence"] = "stop_retained_variant" if old_aa == "*" else "synonymous_variant"
                        elif new_aa == "*":
                            result["consequence"] = "stop_gained"
                        elif old_aa == "*":
                            result["consequence"] = "stop_lost"
                        elif offset < 3 and old_aa == "M":
                            result["consequence"] = "start_lost"
                        else:
                            result["consequence"] = "missense_variant"
            elif len(record.ref) != len(alt):
                fully_coding = all(self.cds_offset(pos) is not None for pos in (start, end)) and any(
                    seg_start <= start <= end <= seg_end for seg_start, seg_end in self.segments)
                if not fully_coding:
                    result["consequence"] = "coding_boundary_variant"
                elif (len(alt) - len(record.ref)) % 3:
                    result["consequence"] = "frameshift_variant"
                else:
                    result["consequence"] = "inframe_insertion" if len(alt) > len(record.ref) else "inframe_deletion"
                result["annotation_flags"].append("indel_not_left_normalized_or_independently_validated")
            else:
                result["consequence"] = "coding_sequence_variant"
                result["annotation_flags"].append("multibase_allele_requires_external_annotation")
        elif any(start <= exon_end and end >= exon_start for exon_start, exon_end in self.exons):
            if self.coding_start is None or self.coding_end is None:
                result["consequence"] = "non_coding_transcript_exon_variant"
            elif (end < self.coding_start and self.strand == 1) or (start > self.coding_end and self.strand == -1):
                result["consequence"] = "5_prime_UTR_variant"
            else:
                result["consequence"] = "3_prime_UTR_variant"
        if splice:
            if IMPACT[splice] > IMPACT[result["consequence"]]:
                if result["consequence"] != "intron_variant":
                    result["additional_consequences"].append(result["consequence"])
                result["consequence"] = splice
            elif splice != result["consequence"]:
                result["additional_consequences"].append(splice)
        return result


def _candidate(record: VCFRecord, alt_index: int, gene: GeneModel, frequencies: dict) -> dict:
    alt = record.alts[alt_index - 1]
    variant_id = f"{record.chrom}:{record.pos}:{record.ref}:{alt}"
    annotation = gene.annotate(record, alt)
    dp, gq = number(record.sample.get("DP"), True), number(record.sample.get("GQ"))
    ad = [number(value, True) for value in record.sample.get("AD", ".").split(",")]
    alt_depth = ad[alt_index] if len(ad) > alt_index else None
    ref_depth = ad[0] if ad else None
    ad_total = sum(ad) if ad and all(depth is not None and depth >= 0 for depth in ad) else None
    vaf = alt_depth / ad_total if alt_depth is not None and alt_depth >= 0 and ad_total else None
    fully_called = record.gt_valid and all(allele is not None for allele in record.gt)
    copies = record.gt.count(alt_index)
    zygosity = "unknown"
    if fully_called:
        if len(record.gt) == 1:
            zygosity = "haploid_alternate"
        elif len(record.gt) == 2:
            zygosity = "homozygous_alternate" if copies == 2 else "heterozygous"
        else:
            zygosity = "non_diploid"
    flags, fail_flags = [], []
    if zygosity == "haploid_alternate":
        flags.append("haploid_call_requires_ploidy_review")
    if record.filter_raw == ".":
        flags.append("filter_not_applied_or_missing")
    elif record.filter_raw != "PASS":
        fail_flags.append("vcf_filter_failed")
    for name, value, minimum in (("depth", dp, 10), ("genotype_quality", gq, 20), ("variant_quality", record.qual, 20)):
        if value is None:
            flags.append(name + "_missing")
        elif value < minimum:
            fail_flags.append(name + "_below_threshold")
    if not fully_called:
        fail_flags.append("genotype_partial_missing_or_invalid")
    if vaf is None:
        flags.append("allelic_depth_missing_or_invalid")
    elif zygosity == "heterozygous" and not 0.2 <= vaf <= 0.8:
        flags.append("heterozygous_allele_balance_outside_0.2_0.8")
    if alt_depth is not None and alt_depth < 3:
        fail_flags.append("alternate_read_support_below_3")
    quality_status = "fail" if fail_flags else "review" if flags else "pass"
    info = parse_info(record.info_raw)
    cohort_af = str(info.get("AF", ".")).split(",")
    population = frequencies.get(variant_id)
    population_af, population_source = None, None
    if isinstance(population, dict) and population.get("source"):
        candidate_af = number(str(population.get("af", ".")))
        if candidate_af is not None and 0 <= candidate_af <= 1:
            population_af, population_source = candidate_af, str(population["source"])
    score_parts = {"predicted_consequence": IMPACT[annotation["consequence"]],
                   "quality": 1 if quality_status == "pass" else -3 if quality_status == "fail" else 0,
                   "external_population_frequency": 0}
    if population_af is not None:
        score_parts["external_population_frequency"] = 1 if population_af < 0.001 else -3 if population_af > 0.01 else 0
    return {
        "variant_id": variant_id, "chrom": record.chrom, "pos": record.pos,
        "ref": record.ref, "alt": alt, "alt_index": alt_index, "dbsnp_id": record.identifier if record.identifier != "." else None,
        **annotation, "genotype": record.sample.get("GT", "."), "zygosity": zygosity,
        "alternate_copies": copies, "qual": record.qual, "filter": record.filter_raw,
        "depth": dp, "genotype_quality": gq, "reference_depth": ref_depth,
        "alternate_depth": alt_depth, "allele_fraction": vaf,
        "quality_status": quality_status, "quality_flags": fail_flags + flags,
        "cohort_info_af": number(cohort_af[alt_index - 1]) if len(cohort_af) >= alt_index else None,
        "population_af": population_af, "population_af_source": population_source,
        "phase_set": record.sample.get("PS"), "physical_phase_id": record.sample.get("PID"),
        "physical_genotype": record.sample.get("PGT"),
        "supplied_cis_partner_ids": [],
        "triage_score": sum(score_parts.values()), "score_components": score_parts,
        "requires_review": True,
    }


def _supplied_phase(first: dict, second: dict) -> str:
    first_block, second_block = first.get("phase_set"), second.get("phase_set")
    if first_block in {None, "", "."} or first_block != second_block or first["chrom"] != second["chrom"]:
        return "unknown"
    first_gt, second_gt = first["genotype"].split("|"), second["genotype"].split("|")
    first_alt, second_alt = str(first["alt_index"]), str(second["alt_index"])
    if len(first_gt) != 2 or len(second_gt) != 2 or first_gt.count(first_alt) != 1 or second_gt.count(second_alt) != 1:
        return "unknown"
    return "supplied_in_cis" if first_gt.index(first_alt) == second_gt.index(second_alt) else "supplied_in_trans"


def _hypotheses(candidates: list[dict]) -> list[dict]:
    hypotheses = []
    by_gene: dict[str, list[dict]] = defaultdict(list)
    seen = set()
    for candidate in candidates:
        key = (candidate["gene"], candidate["variant_id"])
        if key in seen:
            continue
        seen.add(key)
        if candidate["consequence"] not in QUALIFYING_CONSEQUENCES or candidate["quality_status"] == "fail":
            continue
        if UNSUPPORTED_ANNOTATION_FLAGS.intersection(candidate["annotation_flags"]):
            continue
        if candidate["population_af"] is not None and candidate["population_af"] > 0.01:
            continue
        if candidate["zygosity"] == "heterozygous":
            by_gene[candidate["gene"]].append(candidate)
        elif candidate["zygosity"] == "homozygous_alternate":
            hypotheses.append({"gene": candidate["gene"], "type": candidate["zygosity"],
                               "variant_ids": [candidate["variant_id"]], "triage_score": candidate["triage_score"],
                               "phase": "not_applicable", "phase_unconfirmed": False,
                               "population_frequency_status": "unknown" if candidate["population_af"] is None else "annotated",
                               "status": "unvalidated_hypothesis", "review_flags": ["disease_mechanism_and_phenotype_match_unverified"]})
    for gene, variants in by_gene.items():
        for first, second in itertools.combinations(variants, 2):
            phase = _supplied_phase(first, second)
            if phase == "supplied_in_cis":
                first["supplied_cis_partner_ids"].append(second["variant_id"])
                second["supplied_cis_partner_ids"].append(first["variant_id"])
                continue
            hypotheses.append({"gene": gene, "type": "possible_compound_heterozygous",
                               "variant_ids": [first["variant_id"], second["variant_id"]],
                               "triage_score": first["triage_score"] + second["triage_score"],
                               "phase": phase, "phase_unconfirmed": True,
                               "phase_evidence": "phased_gt_same_ps_block" if phase == "supplied_in_trans" else "unresolved",
                               "population_frequency_status": "unknown" if any(v["population_af"] is None for v in (first, second)) else "annotated",
                               "status": "unvalidated_hypothesis",
                               "review_flags": ["caller_phase_requires_independent_validation" if phase == "supplied_in_trans" else "in_trans_not_established",
                                                "disease_mechanism_and_phenotype_match_unverified"]})
    return sorted(hypotheses, key=lambda item: (-item["triage_score"], item["gene"], item["variant_ids"]))


def run(vcf_path: Path, reference_path: Path, output_dir: Path, max_records: int | None = None) -> dict:
    """Audit the stream and write local candidates, recessive hypotheses and QC."""
    vcf_path, reference_path, output_dir = Path(vcf_path), Path(reference_path), Path(output_dir)
    reference = json.loads(reference_path.read_text(encoding="utf-8-sig"))
    if reference.get("assembly") != "GRCh38":
        raise ValueError("Track 1 currently requires a GRCh38 reference")
    genes = [GeneModel(raw) for raw in reference["genes"]]
    if not genes:
        raise ValueError("The hypothesis-led reference panel must contain at least one gene")
    by_chrom: dict[str, list[GeneModel]] = defaultdict(list)
    for gene in genes:
        by_chrom[gene.chrom].append(gene)
    frequencies = reference.get("population_frequencies", {})
    reader = VCFReader(vcf_path)
    counters, chromosomes, filters, genotype_counts, types = Counter(), Counter(), Counter(), Counter(), Counter()
    quality_sums, quality_counts = Counter(), Counter()
    candidates = []
    assembly_checked = False
    for record in reader.records(max_records):
        if not assembly_checked:
            if reader.assembly in {"GRCh37", "conflicting"}:
                raise ValueError(f"VCF/reference assembly mismatch: {reader.assembly}. No candidate annotation was written.")
            assembly_checked = True
        counters["parsed_records"] += 1
        chromosomes[record.chrom] += 1
        for filter_name in record.filter_raw.split(";"):
            filters[filter_name] += 1
        genotype_counts[record.sample.get("GT", ".")] += 1
        for name, value in (("qual", record.qual), ("depth", number(record.sample.get("DP"), True)), ("genotype_quality", number(record.sample.get("GQ")))):
            if value is None:
                counters[f"missing_{name}"] += 1
            else:
                quality_sums[name] += value
                quality_counts[name] += 1
        if not record.gt_valid:
            counters["invalid_genotype_records"] += 1
        if any(allele is None for allele in record.gt):
            counters["missing_or_partial_genotype_records"] += 1
        if not record.called_alt_indices:
            counters["records_without_called_alternate"] += 1
        if len(record.alts) > 1:
            counters["multiallelic_records"] += 1
        for alt in record.alts:
            if alt.startswith("<") or alt == "*" or "[" in alt or "]" in alt:
                types["symbolic_or_breakend"] += 1
            elif len(record.ref) == len(alt) == 1:
                types["snv"] += 1
                if record.ref in "ACGT" and alt in "ACGT" and alt != record.ref:
                    counters["transitions" if {record.ref, alt} in ({"A", "G"}, {"C", "T"}) else "transversions"] += 1
            elif len(record.ref) == len(alt):
                types["multinucleotide"] += 1
            else:
                types["indel"] += 1
        overlapping = [gene for gene in by_chrom.get(record.chrom, []) if record.pos <= gene.end and record.end >= gene.start]
        if overlapping:
            counters["panel_overlapping_records"] += 1
        for alt_index in record.called_alt_indices:
            for gene in overlapping:
                candidates.append(_candidate(record, alt_index, gene, frequencies))
    if reader.assembly in {"GRCh37", "conflicting"}:
        raise ValueError(f"VCF/reference assembly mismatch: {reader.assembly}")
    candidates.sort(key=lambda c: (-c["triage_score"], c["gene"], c["chrom"], c["pos"], c["alt"]))
    hypotheses = _hypotheses(candidates)
    limitations = list(LIMITATIONS)
    if reader.assembly == "unknown":
        limitations.append("VCF assembly could not be verified from its header; coordinates are provisionally interpreted as GRCh38.")
    if not reader.complete:
        limitations.append("This run stopped at max_records; QC and candidates describe only the processed prefix.")
    qc = {
        "analysis": "hypothesis_led_gene_panel_with_streaming_vcf_audit",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_vcf": str(vcf_path.resolve()), "reference_file": str(reference_path.resolve()),
        "input_size_bytes": vcf_path.stat().st_size, "sample_count": len(reader.samples),
        "samples": reader.samples, "assembly": reader.assembly, "assembly_evidence": reader.assembly_evidence,
        "scan_complete": reader.complete, "max_records": max_records,
        "records_seen": reader.records_seen, "malformed_records": reader.malformed_records,
        "malformed_examples": reader.malformed_examples, "counts": dict(counters),
        "chromosome_counts": dict(chromosomes), "filter_counts": dict(filters),
        "genotype_counts": dict(genotype_counts), "allele_type_counts": dict(types),
        "mean_quality": {name: quality_sums[name] / count for name, count in quality_counts.items()},
        "transition_transversion_ratio": counters["transitions"] / counters["transversions"] if counters["transversions"] else None,
        "panel_genes": [gene.symbol for gene in genes], "panel_gene_count": len(genes),
        "transcript_model_warnings": [{"gene": gene.symbol, "transcript_id": gene.transcript,
                                       "expected_cds_length": gene.cds_length, "sequence_length": len(gene.cds)}
                                      for gene in genes if not gene.sequence_valid],
        "candidate_count": len(candidates), "hypothesis_count": len(hypotheses),
        "supplied_cis_pair_exclusions": sum(len(c["supplied_cis_partner_ids"]) for c in candidates) // 2,
        "consequence_counts": dict(Counter(c["consequence"] for c in candidates)),
        "candidate_quality_counts": dict(Counter(c["quality_status"] for c in candidates)),
        "limitations": limitations,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    for filename, value in (("candidates.json", candidates), ("hypotheses.json", hypotheses), ("qc.json", qc)):
        (output_dir / filename).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    columns = ["variant_id", "gene", "transcript_id", "consequence", "protein_change", "genotype", "zygosity",
               "qual", "filter", "depth", "genotype_quality", "reference_depth", "alternate_depth", "allele_fraction",
               "quality_status", "quality_flags", "population_af", "population_af_source", "cohort_info_af", "triage_score", "annotation_flags"]
    with (output_dir / "candidates.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for candidate in candidates:
            writer.writerow({key: "; ".join(value) if isinstance(value, list) else value for key, value in candidate.items()})
    _write_review(output_dir / "review.md", qc, candidates, hypotheses)
    return {"track": 1, "output_dir": str(output_dir.resolve()), "scan_complete": reader.complete,
            "records_seen": reader.records_seen, "candidate_count": len(candidates),
            "hypothesis_count": len(hypotheses), "qc": qc,
            "outputs": {name: str((output_dir / name).resolve()) for name in ("candidates.json", "candidates.csv", "hypotheses.json", "qc.json", "review.md")}}


def _write_review(path: Path, qc: dict, candidates: list[dict], hypotheses: list[dict]) -> None:
    lines = ["# Track 1 local review", "", "This is an exploratory, hypothesis-led analysis for expert review.", "",
             f"- VCF records audited: {qc['records_seen']:,}",
             f"- Complete scan: {qc['scan_complete']}", f"- Reference assembly: {qc['assembly']}",
             f"- Gene panel: {qc['panel_gene_count']} genes", f"- Gene-overlapping called alternate alleles: {len(candidates):,}",
             f"- Recessive genotype hypotheses: {len(hypotheses):,}", "", "## Highest-priority alleles", "",
             "Scores prioritize review and do not estimate pathogenicity or submission probability.", "",
             "| Gene | Variant | Consequence | Protein | Genotype | Quality | Score |", "| --- | --- | --- | --- | --- | --- | --- |"]
    for candidate in candidates[:30]:
        genotype = candidate["genotype"].replace("|", "/")
        lines.append(f"| {candidate['gene']} | {candidate['variant_id']} | {candidate['consequence']} | {candidate['protein_change'] or 'unresolved'} | {genotype} | {candidate['quality_status']} | {candidate['triage_score']} |")
    if not candidates:
        lines.append("| No called alleles intersected the panel | | | | | | |")
    lines += ["", "## Highest-priority recessive hypotheses", "",
              "Pairs require phase, population-frequency, disease-mechanism, and phenotype review. A pair score sums allele scores and is not a probability.", ""]
    for hypothesis in hypotheses[:20]:
        lines.append(f"- **{hypothesis['gene']}** — {hypothesis['type']}: {', '.join(hypothesis['variant_ids'])}; phase {hypothesis['phase']}; score {hypothesis['triage_score']}.")
    if not hypotheses:
        lines.append("No qualifying recessive hypotheses were identified in this scan and panel.")
    lines += ["", "## Next review steps", "",
              "1. Confirm candidate alleles and coverage in reads; inspect low-quality and allelic-balance flags.",
              "2. Add source-labeled population frequencies and curated disease evidence using a local reference database.",
              "3. Validate consequences with an independent transcript-aware annotator and match the phenotype and inheritance model.",
              "4. Establish phase for possible compound-heterozygous pairs and expand beyond the panel if evidence is insufficient.",
              "5. Review CNV/structural and other variant classes before ruling out a genetic diagnosis.", "", "## Limitations", ""]
    lines.extend("- " + limitation for limitation in qc["limitations"])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
